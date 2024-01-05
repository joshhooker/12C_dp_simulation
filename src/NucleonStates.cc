#include "NucleonStates.hh"

NucleonStates* NucleonStates::instance_ = NULL;

NucleonStates::NucleonStates() {
    ReadJSON();
}

NucleonStates* NucleonStates::Instance() {
    if (!instance_) {
        instance_ = new NucleonStates();
    }
    return instance_;
}

void NucleonStates::ReadJSON() {
    // Read and parse nuclear_states.json
    Json::Value config_nuclear_states;
    std::ifstream config_nuclear_states_stream("nuclear_states.json");
    ASSERT_WITH_MESSAGE(config_nuclear_states_stream.is_open(), 
        "Could not find 'nuclear_states.json'\n");
    config_nuclear_states_stream >> config_nuclear_states;
    config_nuclear_states_stream.close();

    std::vector<isotope_struct> isotopes;
    for (uint i = 0; i < config_nuclear_states["Isotopes"].size(); i++) {
        std::string name = config_nuclear_states["Isotopes"][i]["Name"].asString();
        uint charge = config_nuclear_states["Isotopes"][i]["ZA"][0].asUInt();
        uint mass = config_nuclear_states["Isotopes"][i]["ZA"][1].asUInt();

        std::vector<state_struct> states;
        for (uint j = 0; j < config_nuclear_states["Isotopes"][i]["States"].size(); j++) {
            G4double energy = config_nuclear_states["Isotopes"][i]["States"][j]["Energy"].asDouble();
            G4double width  = config_nuclear_states["Isotopes"][i]["States"][j]["Width"].asDouble();
            uint spin = 0;
            G4bool parity = true;
            G4bool use_spin_parity = false;
            if (config_nuclear_states["Isotopes"][i]["States"][j]["Spin-parity"].isString()) {
                std::string spinparity = config_nuclear_states["Isotopes"][i]["States"][j]["Spin-parity"].asString();
                char parity_char = spinparity.back();
                G4bool good_state = false;
                if (parity_char == '+') {
                    parity = true;
                    good_state = true;
                }
                else if (parity_char == '-') {
                    parity = false;
                    good_state = true;
                }

                if (!good_state) {
                    printf("Spin-parity for isotope %s with energy %f not set correctly!\n", name.c_str(), energy);
                    printf("Setting state to 0+\n");
                    parity = true;
                    good_state = false;
                }
                else {
                    std::string spinparity_copy = spinparity;
                    spinparity_copy.pop_back();
                    G4double spind = stof(spinparity_copy)*2;
                    spin = static_cast<int>(spind);
                    use_spin_parity = true;
                }
            }

            G4double probability = 1.;
            if (config_nuclear_states["Isotopes"][i]["States"][j]["Probability"].isDouble()) {
                probability = config_nuclear_states["Isotopes"][i]["States"][j]["Probability"].asDouble();
            }

            G4bool use_angle = false;
            if (config_nuclear_states["Isotopes"][i]["States"][j]["UseAngle"].isBool()) {
                use_angle = config_nuclear_states["Isotopes"][i]["States"][j]["UseAngle"].asBool();
            }

            G4bool angle_cm = false;
            if (config_nuclear_states["Isotopes"][i]["States"][j]["CMAngle"].isBool()) {
                angle_cm = config_nuclear_states["Isotopes"][i]["States"][j]["CMAngle"].asBool();
            }

            G4String angle_file = "";
            if (config_nuclear_states["Isotopes"][i]["States"][j]["AngleFile"].isString()) {
                angle_file = config_nuclear_states["Isotopes"][i]["States"][j]["AngleFile"].asString();
            }

            state_struct state;
            state.energy = energy;
            state.width = width;
            if (use_spin_parity) {
                state.spin2 = spin;
                state.parity = parity;
            }
            state.use_spin_parity = use_spin_parity;
            state.probability = probability;
            state.angle_cm = angle_cm;
            state.angle_file = angle_file;
            state.use_angle = use_angle;

            states.push_back(state);
        }

        // Get total probability of all the states
        G4double totalProbability = 0.;
        for (const auto& state: states) {
            totalProbability += state.probability;
        }

        // Iterate through individual states and calculate probabilities and obtain angular distribution
        for (size_t i = 0; i < states.size(); i++) {
            // Obtain cumulative probabilities and probabilities for rng
            states[i].cumulative_probability = states[i].probability/totalProbability;
            if (i == 0) states[i].probability_rng = states[i].cumulative_probability;
            else states[i].probability_rng = states[i - 1].probability_rng + states[i].cumulative_probability;

            // Check if angle file exists and read file if it does
            struct stat buffer;
            if (stat (states[i].angle_file.c_str(), &buffer) == 0) {
                // Read file and store angles and cross section in vectors
                std::ifstream in(states[i].angle_file.c_str());
                G4double var1, var2;
                std::vector<G4double> angle_vector;
                std::vector<G4double> cross_sec_vector;
                while (in >> var1 >> var2) {
                    angle_vector.push_back(var1);
                    cross_sec_vector.push_back(var2);
                }

                // Get spline of input angle and cross section
                CubicSpline input_spline(angle_vector, cross_sec_vector);

                // Get new vectors based on spline with better angular resolution
                std::vector<G4double> angle_vector_res;
                std::vector<G4double> cross_sec_vector_res;
                int steps = static_cast<int>((angle_vector[angle_vector.size() - 1] - angle_vector[0])/0.1);
                for (int j = 0; j <= steps; j++) {
                    G4double new_angle = angle_vector[0] + j*0.1;
                    G4double new_cross_sec = input_spline(new_angle);
                    angle_vector_res.push_back(new_angle);
                    cross_sec_vector_res.push_back(new_cross_sec);
                }

                // Calculate cumulative probability for the angular distribution
                G4double total_probability = 0.;
                for (auto cross_sec: cross_sec_vector_res) {
                    total_probability += cross_sec;
                }

                std::vector<G4double> cross_sec_cumulative_vector;
                for (size_t j = 0; j < cross_sec_vector_res.size(); j++) {
                    if (j == 0) cross_sec_cumulative_vector.push_back(cross_sec_vector_res[j]/total_probability);
                    else cross_sec_cumulative_vector.push_back(cross_sec_cumulative_vector[j - 1] + cross_sec_vector_res[j]/total_probability);
                }

                // Create spline for angle and cross_section
                CubicSpline angle_cross_sec_spline(cross_sec_cumulative_vector, angle_vector_res);
                states[i].angle_spline = angle_cross_sec_spline;

                // Obtain the low and high end of the angular distribution
                G4double angle_low = angle_vector[0];
                G4double angle_high = angle_vector[angle_vector.size() - 1];
                // Swap if low > high
                if (angle_low > angle_high) {
                    std::swap(angle_low, angle_high);
                }
                states[i].angle_low = angle_low;
                states[i].angle_high = angle_high;

                if (states[i].angle_cm) {
                    printf("CURRENTLY CM ANGLES ARE NOT SUPPORTED\n");
                    states[i].use_angle = false;
                }
            }
        }

        // Set the thresholds
        std::vector<theshold_struct> thresholds;
        for (uint j = 0; j < config_nuclear_states["Isotopes"][i]["Thresholds"].size(); j++) {
            double energy = config_nuclear_states["Isotopes"][i]["Thresholds"][j]["Energy"].asDouble();
            uint threshold_charge = config_nuclear_states["Isotopes"][i]["Thresholds"][j]["Decay"][0].asUInt();
            uint threshold_mass = config_nuclear_states["Isotopes"][i]["Thresholds"][j]["Decay"][1].asUInt();
            theshold_struct threshold = {energy, threshold_charge, threshold_mass};
            thresholds.push_back(threshold);
        }

        isotope_struct isotope = {name, charge, mass, states, thresholds};
        isotopes.push_back(isotope);
    }

    printf("Printing nuclear table to use:\n");
    for (const auto& isotope: isotopes) {
        printf("Name: %5s; Z: %3d; A: %3d\n", isotope.name.c_str(), isotope.charge, isotope.mass);
        printf("\t States:\n");
        for (const auto& state: isotope.states) {
            printf("\t\tEnergy: %7.3f; Width: %7.3f; Input Prob.: %6.3f;  Cumul. Prob: %6.3f", state.energy, state.width, state.probability, state.cumulative_probability);
            if (state.use_angle) printf("; Angular file: %s", state.angle_file.c_str());
            printf("\t\t\n");
        }
        printf("\t Thresholds:\n");
        for (auto threshold: isotope.thresholds) {
            printf("\t\tEnergy: %8.3f; Decay Product [%3d, %3d]\n", threshold.energy, threshold.decay_charge, threshold.decay_mass);
        }
        nucleons_[isotope.charge][isotope.mass] = isotope;
    }
}

std::pair<G4bool, G4double> NucleonStates::GetAngularDistribution(uint charge, uint mass, uint state_number) {
    std::vector<state_struct> states = nucleons_[charge][mass].states;

    if (states.empty()) {
        printf("No states found for nucleus Z: %d A: %d\n", charge, mass);
        return std::make_pair(false, 0.);
    }

    G4double probability = G4UniformRand();
    return std::make_pair(states[state_number].angle_cm, states[state_number].angle_spline(probability));
}

excited_state_struct NucleonStates::GetExcitedLevel(uint charge, uint mass) {
    std::vector<state_struct> states = nucleons_[charge][mass].states;

    if (states.empty()) {
        printf("No states found for nucleus Z: %d A: %d\n", charge, mass);
        excited_state_struct hit = {0, 0., 0., false};
        return hit;
    }

    G4double probability = G4UniformRand();
    G4int state_number = 0;
    G4double ex_energy = 0.;
    G4double ex_width = 0.;

    if (probability < states[0].probability_rng) {
        state_number = 0;
        ex_energy = states[0].energy;
        ex_width = states[0].width;
    }
    else {
        for (size_t i = 1; i < states.size(); i++) {
            if (probability > states[i - 1].probability_rng && probability < states[i].probability_rng) {
                state_number = i;
                ex_energy = states[i].energy;
                ex_width = states[i].width;
            }
        }
    }

    excited_state_struct hit = {state_number, ex_energy, ex_width, states[state_number].use_angle};
    return hit;
}

std::vector<theshold_struct> NucleonStates::GetThresholds(uint charge, uint mass) {
    return nucleons_[charge][mass].thresholds;
}
