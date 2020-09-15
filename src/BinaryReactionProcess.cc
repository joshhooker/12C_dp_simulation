#include "BinaryReactionProcess.hh"

BinaryReactionProcess::BinaryReactionProcess(const G4String& process_name)
  : G4VDiscreteProcess(process_name, fHadronic), scattering_energy_(1e6) {
    SetProcessSubType(111);
}

BinaryReactionProcess::~BinaryReactionProcess() {}

G4double BinaryReactionProcess::GetMeanFreePath(const G4Track& aTrack, G4double, G4ForceCondition* condition) {

    G4double energy = aTrack.GetKineticEnergy()/MeV;

    G4int particle_mass = aTrack.GetDefinition()->GetAtomicMass();
    G4int particle_charge = aTrack.GetDefinition()->GetAtomicNumber();

    const auto* detector_construction = dynamic_cast<const DetectorConstruction*>
    (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    G4LogicalVolume* target_logical = detector_construction->GetTargetVolume();
    G4LogicalVolume* current_volume = aTrack.GetStep()->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();

    G4double position = aTrack.GetPosition().z();

    double mfp = (aTrack.GetTrackID() == 1 &&
                  current_volume == target_logical &&
                  position >= scattering_length_) ? 0. : DBL_MAX;

    // Look at excited name and see if it's in an excited state
    G4String excited_name = aTrack.GetDynamicParticle()->GetDefinition()->GetParticleName();
    size_t pos = excited_name.find('[');
    double measured_excited_energy = 0.;
    G4String beam_name_energy = "";
    if (pos > 1000) {
        measured_excited_energy = 0.;
    }
    else {
        beam_name_energy = excited_name.substr(pos + 1, std::string::npos);
        beam_name_energy.pop_back();
        measured_excited_energy = std::atof(beam_name_energy.c_str())/1000.;
    }

    // Check if above threshold and needs to decay
    NucleonStates* states = NucleonStates::Instance();
    auto thresholds = states->GetThresholds(particle_charge, particle_mass);

    if (measured_excited_energy > 0 && !thresholds.empty()) {
        if (measured_excited_energy > thresholds[0].energy) {
            mfp = 0.;
        }
    }

    if (excited_name == "Be8") {
        mfp = 0.;
    }

    if (excited_name == "Be13") {
        mfp = 0.;
    }

    *condition = NotForced;
    return mfp;
}

G4VParticleChange* BinaryReactionProcess::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep) {

    G4StepPoint* pre_step_point = aStep.GetPreStepPoint();
    G4StepPoint* post_step_point = aStep.GetPostStepPoint();

    if (pre_step_point->GetStepStatus() == fGeomBoundary || post_step_point->GetStepStatus() == fGeomBoundary) {
        return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
    }

    // Determine if particle needs to decay due to excited state
    G4String incoming_particle_name = aTrack.GetDynamicParticle()->GetDefinition()->GetParticleName();
    size_t pos = incoming_particle_name.find('[');
    double measured_excited_energy = 0.;
    G4String incoming_particle_name_energy = "";
    if (pos > 100) {
        measured_excited_energy = 0.;
    }
    else {
        incoming_particle_name_energy = incoming_particle_name.substr(pos + 1, std::string::npos);
        incoming_particle_name_energy.pop_back();
        measured_excited_energy = std::atof(incoming_particle_name_energy.c_str())/1000.;
    }

    NucleonStates* states = NucleonStates::Instance();

    // Get thresholds
    auto thresholds = states->GetThresholds(aTrack.GetParticleDefinition()->GetAtomicNumber(),
        aTrack.GetParticleDefinition()->GetAtomicMass());

    if (!thresholds.empty()) {
        auto threshold = thresholds[0];
        if (measured_excited_energy > threshold.energy) {
            uint light_charge = aTrack.GetDynamicParticle()->GetDefinition()->GetAtomicNumber() - threshold.decay_charge;
            uint light_mass = aTrack.GetDynamicParticle()->GetDefinition()->GetAtomicMass() - threshold.decay_mass;
            return Decay(aTrack, light_charge, light_mass, threshold.decay_charge, threshold.decay_mass);
        }
    }

    aParticleChange.Initialize(aTrack);

    G4ParticleTable* particle_table = G4ParticleTable::GetParticleTable();
    G4String particle_name;

    Calibrations* cal = Calibrations::Instance();
    std::pair<G4int, G4int> target_ion = cal->GetTargetIon();
    std::pair<G4int, G4int> ejectile_ion = cal->GetEjectileIon();

    G4ParticleDefinition* neutron = particle_table->FindParticle(particle_name="neutron");

    G4ParticleDefinition* target_particle;
    if (target_ion.first == 0 && target_ion.second == 1) {
        target_particle = neutron;
    }
    else {
        target_particle = particle_table->GetIonTable()->GetIon(target_ion.first, target_ion.second, 0.0);
    }

    G4ParticleDefinition* ejectile_particle;
    if (ejectile_ion.first == 0 && ejectile_ion.second == 1) {
        ejectile_particle = neutron;
    }
    else {
        ejectile_particle = particle_table->GetIonTable()->GetIon(ejectile_ion.first, ejectile_ion.second, 0.0);
    }

    // Get beam properties
    G4int beam_charge = aTrack.GetParticleDefinition()->GetAtomicNumber();
    G4int beam_mass = aTrack.GetParticleDefinition()->GetAtomicMass();
    auto beam_mass_dbl = static_cast<G4double>(aTrack.GetParticleDefinition()->GetAtomicMass());
    G4double beam_mass_pdg = aTrack.GetParticleDefinition()->GetPDGMass();

    G4int target_charge, target_mass;
    G4int light_recoil_charge, light_recoil_mass;
    G4int heavy_recoil_charge, heavy_recoil_mass;
    G4double target_mass_dbl, light_recoil_mass_dbl, heavy_recoil_mass_dbl;
    G4double target_mass_pdg, light_recoil_mass_pdg, heavy_recoil_mass_pdg;

    // (d, p) Reactions
    G4ParticleDefinition* target_def = target_particle;

    target_charge = target_def->GetAtomicNumber();
    target_mass = target_def->GetAtomicMass();
    target_mass_dbl = static_cast<G4double>(target_def->GetAtomicMass());
    target_mass_pdg = target_def->GetPDGMass();

    G4double energy = aTrack.GetKineticEnergy()/MeV;
    G4double cm_energy = energy*target_mass_dbl/(beam_mass_dbl + target_mass_dbl);

    G4ParticleDefinition* light_recoil_def = ejectile_particle;

    light_recoil_charge = light_recoil_def->GetAtomicNumber();
    light_recoil_mass = light_recoil_def->GetAtomicMass();
    light_recoil_mass_dbl = static_cast<G4double>(light_recoil_def->GetAtomicMass());
    light_recoil_mass_pdg = light_recoil_def->GetPDGMass();

    // Figure out heavy particle charge and mass
    heavy_recoil_charge = beam_charge + target_charge - light_recoil_charge;
    heavy_recoil_mass = beam_mass + target_mass - light_recoil_mass;

    G4double state_energy_fwhm = cal->GetStateEnergyFWHM();
    G4bool use_energy_qvalue = cal->GetUseEnergyQvalue();

    G4int excited_state_number = 0;
    G4double excited_energy = 0.;
    G4double excited_width = 0.;
    G4bool use_angle = false;

    ///////////////////////////
    // Get Excitation Energy //
    ///////////////////////////

    // Sample on excited states
    auto excited_level = states->GetExcitedLevel(heavy_recoil_charge, heavy_recoil_mass);
    excited_state_number = excited_level.state_number;
    excited_energy = excited_level.excited_energy;
    excited_width = excited_level.excited_width;
    use_angle = excited_level.use_angle;
    excited_energy = Lorentzian_shoot(excited_energy, excited_width);
    if (!use_energy_qvalue) {
        excited_energy = excited_energy < 0 ? 0. : excited_energy; // Make sure there are no negative excitation energies
    }

    // Always set excitation energy to 0
    // excited_energy = 0.;

    // Get Heavy Ion Particle Definition
    G4ParticleDefinition* heavy_recoil_def = particle_table->GetIonTable()->GetIon(heavy_recoil_charge, heavy_recoil_mass, excited_energy);
    heavy_recoil_mass_dbl = heavy_recoil_def->GetAtomicMass();
    heavy_recoil_mass_pdg = heavy_recoil_def->GetPDGMass();

    // If the light recoil is heavier than the heavy recoil, switch them
    if (light_recoil_mass > heavy_recoil_mass) {
        std::swap(light_recoil_charge, heavy_recoil_charge);
        std::swap(light_recoil_mass, heavy_recoil_mass);
        std::swap(light_recoil_mass_dbl, heavy_recoil_mass_dbl);
        std::swap(light_recoil_mass_pdg, heavy_recoil_mass_pdg);
        std::swap(light_recoil_def, heavy_recoil_def);
    }

    G4double qvalue = beam_mass_pdg + target_mass_pdg - (light_recoil_mass_pdg + heavy_recoil_mass_pdg);

    // Check if reaction is possible
    if (cm_energy + qvalue < 0.) {
        G4cout << "Not enough energy! Energy: " << energy << "; CM Energy: " << cm_energy << G4endl;
        G4cout << "\t Q Value: " << qvalue << "; Vertex Location: " << aTrack.GetPosition().z() << G4endl;
        return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
    }

    //////////////////////////////
    // Get Angular distribution //
    //////////////////////////////

    uniform_angle_output kinematics;

    uniform_angle_input input = {aTrack, aStep, beam_mass_dbl, target_mass_dbl,
            light_recoil_mass_dbl, heavy_recoil_mass_dbl, energy, qvalue};

    if (!use_angle) {
        kinematics = UniformCMAngle(input);
        // kinematics = UniformLabAngle(input);
    }
    else {
        auto excited_angle = states->GetAngularDistribution(heavy_recoil_charge, heavy_recoil_mass, excited_state_number);

        G4double lab_light_theta = excited_angle.second*deg;
        G4double lab_light_phi   = 2.*M_PI*G4UniformRand()*radian;

        kinematics = CalculateKinematicsLab(input, lab_light_theta, lab_light_phi);
    }

    if (!kinematics.correct_kinematics) return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);

    auto *root_file = RootFile::Instance();
    root_file->GetLightAngleCMDist()->Fill(kinematics.cm_light_theta/deg);
    root_file->GetLightAngleThetaLabDist()->Fill(kinematics.lab_light_theta/deg);
    root_file->GetLightAnglePhiLabDist()->Fill(kinematics.lab_light_phi/deg);
    root_file->GetExcitationEnergyDist()->Fill(excited_energy/MeV);
    root_file->GetQValueDist()->Fill(qvalue/MeV);

    auto *sec1 = new G4Track(new G4DynamicParticle(light_recoil_def, kinematics.lab_light_vector.unit(),
                             kinematics.lab_light_energy*MeV), aTrack.GetGlobalTime(), aTrack.GetPosition());
    sec1->SetUserInformation(new TrackingInformation(energy, cm_energy, kinematics.cm_light_theta,
                                                     kinematics.lab_light_theta, kinematics.cm_light_phi,
                                                     kinematics.lab_light_phi, kinematics.cm_heavy_theta,
                                                     kinematics.lab_heavy_theta, kinematics.lab_light_energy,
                                                     kinematics.lab_heavy_energy, aTrack.GetPosition(), qvalue, excited_energy,
                                                     light_recoil_def, heavy_recoil_def, false, excited_state_number));
    auto *sec2 = new G4Track(new G4DynamicParticle(heavy_recoil_def, kinematics.lab_heavy_vector.unit(), kinematics.lab_heavy_energy*MeV),
                             aTrack.GetGlobalTime(), aTrack.GetPosition());

    aParticleChange.SetNumberOfSecondaries(2);
    aParticleChange.AddSecondary(sec1);
    aParticleChange.AddSecondary(sec2);

    aParticleChange.ProposeEnergy(0.);
    aParticleChange.ProposeTrackStatus(fStopAndKill);

    return &aParticleChange;
}

void BinaryReactionProcess::StartTracking(G4Track* track) {
    G4VProcess::StartTracking(track);	// Apply base class actions

    // Make interaction happen anywhere in the target
    const auto* detectorConstruction = dynamic_cast<const DetectorConstruction*>
    (G4RunManager::GetRunManager()->GetUserDetectorConstruction());

    G4LogicalVolume* target_logical = detectorConstruction->GetTargetVolume();
    G4double target_half_length = (dynamic_cast<G4Tubs*>(target_logical->GetSolid()))->GetZHalfLength();

    scattering_length_ = G4RandFlat::shoot(-target_half_length, target_half_length);
}

G4VParticleChange* BinaryReactionProcess::Decay(const G4Track& aTrack, G4int light_charge, G4int light_mass,
                                                      G4int heavy_charge, G4int heavy_mass) {
    G4ParticleTable* particle_table = G4ParticleTable::GetParticleTable();

    aParticleChange.Initialize(aTrack);

    // If the light recoil is heavier than the heavy recoil, switch them
    if (light_mass > heavy_mass) {
        std::swap(light_charge, heavy_charge);
        std::swap(light_mass, heavy_mass);
    }

    // Setup Particle 1
    G4DynamicParticle* particle1 = new G4DynamicParticle;
    G4ParticleDefinition* particle1_def;
    if (light_charge == 0 && light_mass == 1) {
        G4String particle_name;
        particle1_def = particle_table->FindParticle(particle_name="neutron");
    }
    else {
        if(particle_table->GetIonTable()->FindIon(light_charge, light_mass, 0.)) {
            particle1_def = particle_table->GetIonTable()->FindIon(light_charge, light_mass, 0.);
        }
        else particle1_def = particle_table->GetIonTable()->GetIon(light_charge, light_mass, 0.);
    }
    particle1->SetDefinition(particle1_def);
    G4double particle1_mass_dbl = static_cast<G4double>(particle1_def->GetAtomicMass());
    G4double particle1_mass_pdg = particle1_def->GetPDGMass()/CLHEP::amu_c2;

    // Setup Particle 2
    G4DynamicParticle* particle2 = new G4DynamicParticle;
    G4ParticleDefinition* particle2_def = G4IonTable::GetIonTable()->GetIon(heavy_charge, heavy_mass, 0.);
    particle2->SetDefinition(particle2_def);
    G4double particle2_mass_dbl = static_cast<G4double>(particle2_def->GetAtomicMass());
    G4double particle2_mass_pdg = particle2_def->GetPDGMass()/CLHEP::amu_c2;

    G4double qvalue = aTrack.GetDynamicParticle()->GetDefinition()->GetPDGMass() -
        (particle1_def->GetPDGMass() + particle2_def->GetPDGMass());
    G4double cm_energy = qvalue;

    if (cm_energy < 0.) return &aParticleChange; // Below the threshold

    // Generate random CM Angles
    G4double cm_theta = M_PI*G4UniformRand();; // 0 to pi
    G4double cm_phi = 2.*M_PI*G4UniformRand(); // 0 to 2pi

    G4double sqrt_factor = 2.*particle1->GetMass()*cm_energy*particle2_mass_dbl/(particle1_mass_dbl + particle2_mass_dbl);
    if (sqrt_factor < 0.) return &aParticleChange;

    G4double p1 = sqrt(2.*particle1->GetMass()*cm_energy*particle2_mass_dbl/(particle1_mass_dbl + particle2_mass_dbl));
    G4double p2 = sqrt(2.*particle2->GetMass()*cm_energy*particle1_mass_dbl/(particle1_mass_dbl + particle2_mass_dbl));

    // Get momentum directions
    G4ThreeVector p_new_1 = G4ThreeVector(p1*sin(cm_theta)*sin(cm_phi), p1*sin(cm_theta)*cos(cm_phi), p1*cos(cm_theta));
    G4ThreeVector p_new_2 = -p_new_1;
    G4ThreeVector p_parent = aTrack.GetMomentum();
    p_new_1 += p_parent*(particle1_mass_dbl/(particle1_mass_dbl + particle2_mass_dbl));
    p_new_2 += p_parent*(particle2_mass_dbl/(particle1_mass_dbl + particle2_mass_dbl));
    particle1->SetMomentum(p_new_1);
    particle2->SetMomentum(p_new_2);

    G4double total_mom_1 = p_new_1.getR();
    G4double total_mom_2 = p_new_2.getR();
    particle1->SetKineticEnergy((total_mom_1*total_mom_1)/(2.*particle1->GetMass()));
    particle2->SetKineticEnergy((total_mom_2*total_mom_2)/(2.*particle2->GetMass()));

    G4Track* sec1 = new G4Track(particle1, aTrack.GetGlobalTime(), aTrack.GetPosition());
    G4Track* sec2 = new G4Track(particle2, aTrack.GetGlobalTime(), aTrack.GetPosition());

    aParticleChange.SetNumberOfSecondaries(2);
    aParticleChange.AddSecondary(sec1);
    aParticleChange.AddSecondary(sec2);

    aParticleChange.ProposeEnergy(0.);
    aParticleChange.ProposeTrackStatus(fStopAndKill);

    return &aParticleChange;
}

uniform_angle_output BinaryReactionProcess::UniformCMAngle(uniform_angle_input input) {
    G4double beam_mass         = input.beam_mass;
    G4double target_mass       = input.target_mass;
    G4double light_recoil_mass = input.light_recoil_mass;
    G4double heavy_recoil_mass = input.heavy_recoil_mass;

    G4double energy       = input.beam_energy;
    G4double q_value      = input.q_value;
    G4double total_energy = energy + q_value;

    G4double cm_energy       = energy*target_mass/(beam_mass + target_mass);
    G4double cm_total_energy = cm_energy + q_value;

    const G4Track& aTrack = input.aTrack;
    const G4Step& aStep   = input.aStep;

    G4double cm_light_theta = M_PI*G4UniformRand()*radian;
    G4double cm_light_phi   = 2.*M_PI*G4UniformRand()*radian;

    G4ThreeVector momentum_direction = aTrack.GetMomentumDirection();
    G4ThreeVector v = G4ThreeVector(0., 0., 1.).cross(momentum_direction);
    G4double rot_angle = acos(momentum_direction.z());
    G4ThreeVector direction = G4ThreeVector(sin(cm_light_theta)*cos(cm_light_phi), sin(cm_light_theta)*sin(cm_light_phi), cos(cm_light_theta));
    if (v.getR() > 0) direction.rotate(v, rot_angle);

    G4double momentum_amplitude_light = sqrt(2.*light_recoil_mass*CLHEP::amu_c2*cm_total_energy*(heavy_recoil_mass/(light_recoil_mass + heavy_recoil_mass)));
    G4double momentum_amplitude_heavy = sqrt(2.*heavy_recoil_mass*CLHEP::amu_c2*cm_total_energy*(light_recoil_mass/(light_recoil_mass + heavy_recoil_mass)));

    G4ThreeVector momentum_vector_light = momentum_amplitude_light*direction;
    G4ThreeVector momentum_vector_heavy = -momentum_vector_light;

    G4ThreeVector pN = aTrack.GetMomentum();

    momentum_vector_light += pN*(light_recoil_mass/(light_recoil_mass + heavy_recoil_mass));
    momentum_vector_heavy += pN*(heavy_recoil_mass/(light_recoil_mass + heavy_recoil_mass));

    G4double lab_light_theta = momentum_vector_light.theta();
    G4double lab_light_phi   = momentum_vector_light.phi();

    return CalculateKinematicsLab(input, lab_light_theta, lab_light_phi);
}

uniform_angle_output BinaryReactionProcess::UniformLabAngle(uniform_angle_input input) {
    G4double beam_mass         = input.beam_mass;
    G4double target_mass       = input.target_mass;
    G4double light_recoil_mass = input.light_recoil_mass;
    G4double heavy_recoil_mass = input.heavy_recoil_mass;

    G4double energy       = input.beam_energy;
    G4double q_value      = input.q_value;
    G4double total_energy = energy + q_value;

    G4double b = beam_mass*light_recoil_mass/(beam_mass + target_mass)/(light_recoil_mass + heavy_recoil_mass)*
               (energy/total_energy);
    G4double d = target_mass*heavy_recoil_mass/(beam_mass + target_mass)/(light_recoil_mass + heavy_recoil_mass)*
               (1. + beam_mass/target_mass*q_value/total_energy);

    G4double max_angle = (b < d) ? M_PI : asin(sqrt(d/b));

    G4double lab_light_theta = max_angle*G4UniformRand()*radian;
    G4double lab_light_phi   = 2.*M_PI*G4UniformRand()*radian;

    return CalculateKinematicsLab(input, lab_light_theta, lab_light_phi);
}

uniform_angle_output BinaryReactionProcess::CalculateKinematicsLab(uniform_angle_input input, G4double lab_light_theta, G4double lab_light_phi) {
    G4double beam_mass         = input.beam_mass;
    G4double target_mass       = input.target_mass;
    G4double light_recoil_mass = input.light_recoil_mass;
    G4double heavy_recoil_mass = input.heavy_recoil_mass;

    G4double energy       = input.beam_energy;
    G4double q_value      = input.q_value;
    G4double total_energy = energy + q_value;

    G4double cm_energy       = energy*target_mass/(beam_mass + target_mass);
    G4double cm_total_energy = cm_energy + q_value;

    const G4Track& aTrack = input.aTrack;
    const G4Step& aStep   = input.aStep;

    G4double mass_factor = (beam_mass + target_mass)*(light_recoil_mass + heavy_recoil_mass);

    G4double B = beam_mass*light_recoil_mass/mass_factor*(energy/total_energy);
    G4double D = target_mass*heavy_recoil_mass/mass_factor*
        (1. + beam_mass/target_mass*q_value/total_energy);

    G4ThreeVector light_vector(sin(lab_light_theta)*cos(lab_light_phi), sin(lab_light_theta)*sin(lab_light_phi),
                               cos(lab_light_theta));
    G4ThreeVector beam_vector(0., 0., 1.);

    G4double sqrt_factor = D/B - sin(lab_light_theta)*sin(lab_light_theta);
    if (sqrt_factor < 0.) {
        uniform_angle_output output = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        G4ThreeVector(0., 0., 0.), G4ThreeVector(0., 0., 0.), 0., 0., false};

        return output;
    }

    G4double lab_light_energy  = (B <= D) ? total_energy*B*pow(cos(lab_light_theta) + sqrt(D/B - sin(lab_light_theta)*sin(lab_light_theta)), 2.) :
                                total_energy*B*pow(cos(lab_light_theta) - sqrt(D/B - sin(lab_light_theta)*sin(lab_light_theta)), 2.);
    G4double lab_light_energy2 = (B <= D) ? total_energy*B*pow(cos(lab_light_theta + 0.001) + sqrt(D/B - sin(lab_light_theta + 0.001)*sin(lab_light_theta + 0.001)), 2.) :
                                 total_energy*B*pow(cos(lab_light_theta + 0.001) - sqrt(D/B - sin(lab_light_theta + 0.001)*sin(lab_light_theta + 0.001)), 2.);

    G4double lab_heavy_energy = total_energy - lab_light_energy;

    // Make sure the energies of both recoils are greater than 0
    if (lab_light_energy <= 0. || lab_heavy_energy <= 0.) {
        G4cout << "NEGATIVE ENERGIES! light energy: " << lab_light_energy << '\t' << "; heavy energy: " << lab_heavy_energy << G4endl;
        uniform_angle_output output = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        G4ThreeVector(0., 0., 0.), G4ThreeVector(0., 0., 0.), 0., 0., false};

        return output;
    }

    // Don't continue if light energy less than 0.1 MeV. Si threshold is 0.2 MeV
    if (lab_light_energy < 0.05) {
        // G4cout << "lab_light_energy too small: " << lab_light_energy << " MeV" << G4endl;
        uniform_angle_output output = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        G4ThreeVector(0., 0., 0.), G4ThreeVector(0., 0., 0.), 0., 0., false};

        return output;
    }


    if (lab_heavy_energy < 0.001) {
        // G4cout << "lab_heavy_energy too small: " << lab_heavy_energy << " MeV" << G4endl;
        uniform_angle_output output = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        G4ThreeVector(0., 0., 0.), G4ThreeVector(0., 0., 0.), 0., 0., false};

        return output;
    }

    G4double val_sqrt_factor = lab_light_energy/total_energy/D;
    if (val_sqrt_factor < 0.) {
        uniform_angle_output output = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        G4ThreeVector(0., 0., 0.), G4ThreeVector(0., 0., 0.), 0., 0., false};

        return output;
    }

    G4double val2_sqrt_factor = lab_light_energy2/total_energy/D;
    if (val2_sqrt_factor < 0.) {
        uniform_angle_output output = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        G4ThreeVector(0., 0., 0.), G4ThreeVector(0., 0., 0.), 0., 0., false};

        return output;
    }

    G4double val  = sqrt(lab_light_energy/total_energy/D)*sin(lab_light_theta);
    G4double val2 = sqrt(lab_light_energy2/total_energy/D)*sin(lab_light_theta + 0.001);

    G4double cm_light_theta = (val2 > val) ? asin(val) : M_PI - asin(val);
    G4double cm_light_phi   = lab_light_phi;

    G4double asin_factor = sqrt(light_recoil_mass/(heavy_recoil_mass)*
                                 lab_light_energy/lab_heavy_energy)*sin(lab_light_theta);
    if (fabs(asin_factor) > 1.) {
        uniform_angle_output output = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        G4ThreeVector(0., 0., 0.), G4ThreeVector(0., 0., 0.), 0., 0., false};

        return output;
    }

    G4double lab_heavy_theta = asin(sqrt(light_recoil_mass/(heavy_recoil_mass)*
                                 lab_light_energy/lab_heavy_energy)*sin(lab_light_theta));
    G4ThreeVector heavy_vector(-1.*sin(lab_heavy_theta)*cos(lab_light_phi), -1.*sin(lab_heavy_theta)*sin(lab_light_phi),
                               cos(lab_heavy_theta));
    G4double cm_heavy_theta = M_PI - cm_light_theta;

    G4ThreeVector momentum_direction = aTrack.GetMomentumDirection();
    G4ThreeVector v = beam_vector.cross(momentum_direction);
    G4double rot_angle = acos(momentum_direction.z());

    G4ThreeVector lab_light_vector(sin(lab_light_theta)*cos(lab_light_phi),
                            sin(lab_light_theta)*sin(lab_light_phi),
                            cos(lab_light_theta));
    if (v.getR() > 0) lab_light_vector = lab_light_vector.rotate(v, rot_angle);

    G4ThreeVector lab_heavy_vector(-1.*sin(lab_heavy_theta)*cos(lab_light_phi),
                            -1.*sin(lab_heavy_theta)*sin(lab_light_phi),
                            cos(lab_heavy_theta));
    if (v.getR() > 0) lab_heavy_vector = lab_heavy_vector.rotate(v, rot_angle);

    uniform_angle_output output = {energy, cm_energy, q_value, total_energy, cm_light_theta,
        cm_light_phi, cm_heavy_theta, lab_light_theta, lab_light_phi, lab_heavy_theta,
        lab_light_vector, lab_heavy_vector, lab_light_energy, lab_heavy_energy, true};

    return output;
}

G4double BinaryReactionProcess::Lorentzian_shoot(G4double mean, G4double gamma) {
    double u;
    do {
        u = G4UniformRand();
    }
    while (u == 0.5);

    return mean + gamma*tan(M_PI*u);
}