#include "Calibrations.hh"

Calibrations* Calibrations::instance_ = NULL;

Calibrations::Calibrations() {
    double gamma0 = 50.; // mm
    double dist = 80.; // mm
    double w = 4.9375; // mm
    for(int i = 0; i < 16; i++) {
        double gamma_i_low = gamma0 + (static_cast<double>(i) + 1.)*w;
        double gamma_i = gamma0 + (static_cast<double>(i) + 0.5)*w;
        double gamma_i_high = gamma0 + (static_cast<double>(i))*w;
        yu_angle_low_bin_[i] = 180. - atan2(gamma_i_low, dist)*180./M_PI;
        yu_angle_[i] = 180. - atan2(gamma_i, dist)*180./M_PI;
        yu_angle_high_bin_[i] = 180. - atan2(gamma_i_high, dist)*180./M_PI;
        // G4cout << i << '\t' << gamma_i_high << '\t' << gamma_i << '\t' << gamma_i_low << '\t' << yu_angle_low_bin_[i] << '\t' << yu_angle_[i] << '\t' << yu_angle_high_bin_[i] << G4endl;
    }

    // Read parameters from config.json
    ReadJSON();
}

Calibrations* Calibrations::Instance() {
    if(!instance_) {
        instance_ = new Calibrations();
    }
    return instance_;
}

void Calibrations::ReadJSON() {
    // Read and parse config.json
    Json::Value config;
    std::ifstream config_stream("config.json");
    ASSERT_WITH_MESSAGE(config_stream.is_open(),
        "Could not find 'nuclear_states.json'\n");
    config_stream >> config;
    config_stream.close();

    output_file_ = config["outputFile"].asString();

    beam_ion_.first = config["beamIon[Z,A]"][0].asUInt();
    beam_ion_.second = config["beamIon[Z,A]"][1].asUInt();

    target_ion_.first = config["targetIon[Z,A]"][0].asUInt();
    target_ion_.second = config["targetIon[Z,A]"][1].asUInt();

    ejectile_ion_.first = config["ejectileIon[Z,A]"][0].asUInt();
    ejectile_ion_.second = config["ejectileIon[Z,A]"][1].asUInt();

    beam_energy_ = config["beamEnergy"].asDouble()*MeV;

    beam_offset_x = config["beamOffset"][0].asDouble()*mm;
    beam_offset_y = config["beamOffset"][1].asDouble()*mm;

    target_thickness_ = config["targetThickness"].asDouble()*um;

    yu_energy_fwhm_ = config["yuEnergyFWHM"].asDouble()*MeV;
    yu_threshold_ = config["yuThreshold"].asDouble();

    state_energy_fwhm_ = config["stateEnergyFWHM"].asDouble()*MeV;

    background_amount_ = config["backgroundAmount"].asDouble();

    use_energy_qvalue_ = config["energyQvalue"].asBool();

    printf("Calibration parameters:\n");
    printf("\t Output filename: %s\n", output_file_.c_str());
    printf("\t Beam ion (Z, A):     (%d, %d)\n", beam_ion_.first, beam_ion_.second);
    printf("\t Target ion (Z, A):   (%d, %d)\n", target_ion_.first, target_ion_.second);
    printf("\t Ejectile ion (Z, A): (%d, %d)\n", ejectile_ion_.first, ejectile_ion_.second);
    printf("\t Beam Energy: %7.4f; Beam offset (mm): (%f, %f);\n", beam_energy_, beam_offset_x, beam_offset_y);
    printf("\t Target thickness (um): %7.4f; Yu Energy Resolution (MeV): %7.4f; Yu Energy Threshold (MeV) : %7.4f\n", target_thickness_, yu_energy_fwhm_, yu_threshold_);
    printf("\t Background Amount (%%): %7.4f; Use Q-values instead of Excitation Energy for nuclear_states.json: %s\n", background_amount_, use_energy_qvalue_ ? "True" : "False");
}

void Calibrations::GetTargetProperties(G4Material* material) {
    G4Material* target_material = material;

    G4ParticleTable* particle_table = G4ParticleTable::GetParticleTable();
    G4String particle_name;

    G4ParticleDefinition* neutron = particle_table->FindParticle(particle_name="neutron");

    G4ParticleDefinition* beam_particle = particle_table->GetIonTable()->GetIon(beam_ion_.first, beam_ion_.second, 0.0);

    G4ParticleDefinition* ejectile_particle;
    if(ejectile_ion_.first == 0 && ejectile_ion_.second == 1) {
        ejectile_particle = neutron;
    }
    else {
        ejectile_particle = particle_table->GetIonTable()->GetIon(ejectile_ion_.first, ejectile_ion_.second, 0.0);
    }

    G4EmCalculator em_cal;

    beam_mass_pdg_ = beam_particle->GetPDGMass();
    beam_mass_dbl_ = beam_particle->GetAtomicMass();

    ejectile_mass_pdg_ = ejectile_particle->GetPDGMass();
    ejectile_mass_dbl_ = ejectile_particle->GetAtomicMass();

    G4int heavy_recoil_charge = beam_ion_.first + target_ion_.first - ejectile_ion_.first;
    G4int heavy_recoil_mass = beam_ion_.second + target_ion_.second - ejectile_ion_.second;

    G4ParticleDefinition* heavy_recoil_particle = particle_table->GetIonTable()->GetIon(heavy_recoil_charge, heavy_recoil_mass, 0.0);

    heavy_recoil_mass_pdg_ = heavy_recoil_particle->GetPDGMass();
    heavy_recoil_mass_dbl_ = heavy_recoil_particle->GetAtomicMass();

    FILE* file_beam = fopen("beam_energyloss.dat", "w");
    for(G4double energy = 0.00001; energy < 0.001; energy += 0.00001) {
        G4double energy_particle = energy*MeV;
        G4double dedx = em_cal.ComputeTotalDEDX(energy_particle, beam_particle, target_material)/(MeV/mm);
        fprintf(file_beam, "%f %f \n", energy, dedx);
    }
    for(G4double energy = 0.001; energy < 100.; energy += 0.001) {
        G4double energy_particle = energy*MeV;
        G4double dedx = em_cal.ComputeTotalDEDX(energy_particle, beam_particle, target_material)/(MeV/mm);
        fprintf(file_beam, "%f %f \n", energy, dedx);
    }
    for(G4double energy = 100.; energy < 1000.1; energy += 1.) {
        G4double energy_particle = energy*MeV;
        G4double dedx = em_cal.ComputeTotalDEDX(energy_particle, beam_particle, target_material)/(MeV/mm);
        fprintf(file_beam, "%f %f \n", energy, dedx);
    }
    fflush(file_beam);
    fclose(file_beam);

    FILE* file_ejectile = fopen("ejectile_energyloss.dat", "w");
    for(G4double energy = 0.00001; energy < 0.001; energy += 0.00001) {
        G4double energy_particle = energy*MeV;
        G4double dedx = em_cal.ComputeTotalDEDX(energy_particle, ejectile_particle, target_material)/(MeV/mm);
        fprintf(file_ejectile, "%f %f \n", energy, dedx);
    }
    for(G4double energy = 0.001; energy < 100.; energy += 0.001) {
        G4double energy_particle = energy*MeV;
        G4double dedx = em_cal.ComputeTotalDEDX(energy_particle, ejectile_particle, target_material)/(MeV/mm);
        fprintf(file_ejectile, "%f %f \n", energy, dedx);
    }
    for(G4double energy = 100.; energy < 1000.1; energy += 1.) {
        G4double energy_particle = energy*MeV;
        G4double dedx = em_cal.ComputeTotalDEDX(energy_particle, ejectile_particle, target_material)/(MeV/mm);
        fprintf(file_ejectile, "%f %f \n", energy, dedx);
    }
    fflush(file_ejectile);
    fclose(file_ejectile);
}

void Calibrations::ReaddEdxTables() {
    beam_energyloss = new EnergyLoss();
    beam_energyloss->ReadBasicdEdx("beam_energyloss.dat");

    ejectile_energyloss = new EnergyLoss();
    ejectile_energyloss->ReadBasicdEdx("ejectile_energyloss.dat");
}

