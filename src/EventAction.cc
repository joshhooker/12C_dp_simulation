#include "EventAction.hh"

auto predYy = [] (const yy_detect_struct& lhs, const yy_detect_struct& rhs) {return lhs.energy > rhs.energy;};
auto predSd = [] (const sd_detect_struct& lhs, const sd_detect_struct& rhs) {return lhs.energy > rhs.energy;};

EventAction::EventAction() :
    G4UserEventAction() {

    for(int & i : yu_hcid_) {
        i = -1;
    }
    for(G4int i = 0; i < 24; i++) {
        sd1_hcid_[i] = -1;
        sd2_hcid_[i] = -1;
    }
}

EventAction::~EventAction() = default;

void EventAction::BeginOfEventAction(const G4Event*) {
    char name[256];

    G4SDManager* sd_manager = G4SDManager::GetSDMpointer();
    for (G4int i = 0; i < 16; i++) {
        sprintf(name, "yu%d/genCollection", i + 1);
        yu_hcid_[i] = sd_manager->GetCollectionID(name);
    }

    for (G4int i = 0; i < 24; i++) {
        sprintf(name, "sd1%d/genCollection", i + 1);
        sd1_hcid_[i] = sd_manager->GetCollectionID(name);

        sprintf(name, "sd2%d/genCollection", i + 1);
        sd2_hcid_[i] = sd_manager->GetCollectionID(name);
    }

    energy_ = 0.;
    light_angle_theta_cm_ = 0.;
    light_angle_theta_lab_ = 0.;
    light_energy_ = 0.;
    vertex_z_ = 0.;
}

void EventAction::EndOfEventAction(const G4Event* event) {
    G4HCofThisEvent* hce = event->GetHCofThisEvent();
    if (!hce) {
        G4ExceptionDescription msg;
        msg << "No hits collection of this event found.\n";
            G4Exception("EventAction::EndOfEventAction()",
		    "Code001", JustWarning, msg);
        return;
    }

    GenHitsCollection* yu_hc[16];
    for (int i = 0; i < 16; i++) {
        yu_hc[i] = dynamic_cast<GenHitsCollection*>(hce->GetHC(yu_hcid_[i]));
    }

    GenHitsCollection* sd1_hc[24];
    GenHitsCollection* sd2_hc[24];

    for (int i = 0; i < 24; i++) {
        sd1_hc[i] = dynamic_cast<GenHitsCollection *>(hce->GetHC(sd1_hcid_[i]));
        sd2_hc[i] = dynamic_cast<GenHitsCollection *>(hce->GetHC(sd2_hcid_[i]));
    }

    std::vector<yy_detect_struct> yu_detect_vec;

    Calibrations* cal = Calibrations::Instance();
    double yu_energy_fwhm = cal->GetYuEnergyFWHM();

    // Loop through Yu
    for (G4int i = 0; i < 16; i++) {
        for (size_t j = 0; j < yu_hc[i]->entries(); ++j) {
            G4double energy = (*yu_hc[i])[j]->GetTotalEnergy();
            energy = G4RandGauss::shoot(energy, yu_energy_fwhm/2.355);
            if (energy < 0.0001) continue;
            yy_detect_struct hit = {(*yu_hc[i])[j]->GetID(), i, energy, (*yu_hc[i])[j]->GetPosition()};
            yu_detect_vec.push_back(hit);
        }
    }

    if (yu_detect_vec.empty()) return;

    // Loop through Sd1
    std::vector<sd_detect_struct> sd1_detect_vec;
    for (G4int i = 0; i < 24; i++) {
        for (size_t j = 0; j < sd1_hc[i]->entries(); ++j) {
            G4double energy = (*sd1_hc[i])[j]->GetTotalEnergy();
            energy = G4RandGauss::shoot(energy, yu_energy_fwhm/2.355);
            if (energy < 0.0001) continue;
            sd_detect_struct hit = {(*sd1_hc[i])[j]->GetID(), i, energy, (*sd1_hc[i])[j]->GetPosition()};
            sd1_detect_vec.push_back(hit);
        }
    }

    if (sd1_detect_vec.empty()) return;

    // Loop through Sd2
    std::vector<sd_detect_struct> sd2_detect_vec;
    for (G4int i = 0; i < 24; i++) {
        for (size_t j = 0; j < sd2_hc[i]->entries(); ++j) {
            G4double energy = (*sd2_hc[i])[j]->GetTotalEnergy();
            energy = G4RandGauss::shoot(energy, yu_energy_fwhm/2.355);
            if (energy < 0.0001) continue;
            sd_detect_struct hit = {(*sd2_hc[i])[j]->GetID(), i, energy, (*sd2_hc[i])[j]->GetPosition()};
            sd2_detect_vec.push_back(hit);
        }
    }

    if (sd2_detect_vec.empty()) return;

    std::sort(yu_detect_vec.begin(), yu_detect_vec.end(), predYy);
    std::sort(sd1_detect_vec.begin(), sd1_detect_vec.end(), predSd);
    std::sort(sd2_detect_vec.begin(), sd2_detect_vec.end(), predSd);

    yu_det_ = yu_detect_vec[0].det;
    yu_ring_ = yu_detect_vec[0].ring;
    yu_position_ = yu_detect_vec[0].position;

    yu_energy_ = 0.;
    for (const auto& yuDet: yu_detect_vec) {
        yu_energy_ += yuDet.energy;
    }

    sd1_det_ = sd1_detect_vec[0].det;
    sd1_ring_ = sd1_detect_vec[0].ring;
    sd1_energy_ = sd1_detect_vec[0].energy;

    sd2_det_ = sd2_detect_vec[0].det;
    sd2_ring_ = sd2_detect_vec[0].ring;
    sd2_energy_ = sd2_detect_vec[0].energy;

    if (yu_energy_ < cal->GetYuThreshold()) return;

    G4PrimaryVertex* primary_vertex = event->GetPrimaryVertex();
    G4PrimaryParticle* primary_particle = primary_vertex->GetPrimary();
    beam_energy_ = primary_particle->GetKineticEnergy();

    // Get position of proton
    double yu_ring_width = 4.9375; // in mm
    double yu_angle = 21. + 45.*yu_det_; // in degrees
    double yu_radius = 50. + yu_ring_width*(yu_ring_ + 0.5); // in mm
    double yuX = yu_radius*sin(yu_angle*M_PI/180.);
    double yuY = yu_radius*cos(yu_angle*M_PI/180.);
    // yuX = yu_position_.x();
    // yuY = yu_position_.y();

    auto beam_offset = cal->GetBeamOffset();

    TVector3 v_target(0., 0., 1.);
    TVector3 v_hit(yuX - beam_offset.first, yuY - beam_offset.second, 80.8);
    double angle = 180. - v_target.Angle(v_hit)*180./M_PI;

    EnergyLoss* beam_energyloss = cal->GetBeamEnergyLoss();
    EnergyLoss* ejectile_energyloss = cal->GetEjectileEnergyLoss();

    G4double beam_mass_dbl = cal->GetBeamMassDBL();

    G4double ejectile_mass_dbl = cal->GetEjectileMassDBL();

    G4double heavy_recoil_mass_dbl = cal->GetHeavyRecoilMassDBL();

    // D2 target
    // double half_target_thickness = cal->GetTargetThickness()/4.;
    double half_target_thickness = cal->GetTargetThickness()/2.;
    double distance = half_target_thickness;
    double distance_proton = distance/fabs(cos((180. - angle)*M_PI/180.));
    double proton_energy = ejectile_energyloss->AddBack(yu_energy_, distance_proton);

    // Beam energy at center of target
    G4double beam_e = cal->GetBeamEnergy();
    beam_e = beam_energyloss->CalcRemainder(beam_e, half_target_thickness*2.0);

    double e1 = beam_energyloss->CalcRemainder(cal->GetBeamEnergy(), distance);
    double e3 = proton_energy;

    double mass1 = beam_mass_dbl;
    double mass3 = ejectile_mass_dbl;
    double mass4 = heavy_recoil_mass_dbl;

    // ideal case
    // e1 = energy_/MeV;
    // e3 = light_energy_/MeV;
    // angle = light_angle_theta_lab_/deg;

    double qvalue = (e3*(mass3 + mass4)/mass4 - e1*(mass4 - mass1)/mass4 - (2.*sqrt(mass1*mass3*e1*e3)*cos(angle*M_PI/180.))/mass4)*MeV;

    auto *root_file = RootFile::Instance();

    root_file->GetLightAngleCM()->Fill(light_angle_theta_cm_/deg);
    root_file->GetLightAngleThetaLab()->Fill(light_angle_theta_lab_/deg);
    root_file->GetLightAnglePhiLab()->Fill(light_angle_phi_lab_/deg);
    root_file->GetLightEnergy()->Fill(light_energy_/MeV);
    root_file->GetEnergy()->Fill(energy_/MeV);
    root_file->GetVertexZ()->Fill(vertex_z_/um);

    root_file->GetExcitationEnergy()->Fill(excitation_energy_/MeV);
    root_file->GetQValue()->Fill(qvalue_/MeV);
    root_file->GetCalcQValue()->Fill(qvalue/MeV);
    if (background_) {
        root_file->GetCalcQValueBackground()->Fill(qvalue/MeV);
    }
    else {
        root_file->GetCalcQValueStates()->Fill(qvalue/MeV);
        root_file->GetCalcQvalueState(state_number_)->Fill(qvalue/MeV);
    }
    root_file->GetCalcQValueFull()->Fill(qvalue/MeV);
    root_file->GetQValueDiff()->Fill(qvalue/MeV - qvalue_/MeV);

    if (yu_ring_ < 8) {
        root_file->GetCalcQValue07()->Fill(qvalue/MeV);
    }
    else {
        root_file->GetCalcQValue815()->Fill(qvalue/MeV);
    }

    root_file->GetProtonEnergyError()->Fill(proton_energy - light_energy_/MeV);
    // root_file->GetProtonAngleError()->Fill();
    root_file->GetBeamEnergyError()->Fill(e1 - energy_/MeV);
    root_file->GetQvalueError()->Fill(qvalue/MeV - qvalue_/MeV);

    root_file->GetSd1Energy()->Fill(sd1_energy_/MeV);
    root_file->GetSd2Energy()->Fill(sd2_energy_/MeV);
    root_file->GetSdTotalEnergy()->Fill(sd1_energy_/MeV + sd2_energy_/MeV);

    root_file->GetQValueVertex()->Fill(qvalue_/MeV, vertex_z_/um);

    // fill ttree
    auto analysis_manager = G4RootAnalysisManager::Instance();

    analysis_manager->FillNtupleIColumn(0, yu_det_);
    analysis_manager->FillNtupleIColumn(1, yu_ring_);
    analysis_manager->FillNtupleDColumn(2, yu_energy_/MeV);

    analysis_manager->FillNtupleIColumn(3, sd1_det_);
    analysis_manager->FillNtupleIColumn(4, sd1_ring_);
    analysis_manager->FillNtupleDColumn(5, sd1_energy_/MeV);

    analysis_manager->FillNtupleIColumn(6, sd2_det_);
    analysis_manager->FillNtupleIColumn(7, sd2_ring_);
    analysis_manager->FillNtupleDColumn(8, sd2_energy_/MeV);

    analysis_manager->FillNtupleDColumn(9, energy_/MeV);
    analysis_manager->FillNtupleDColumn(10, light_energy_/MeV);
    analysis_manager->FillNtupleDColumn(11, qvalue_/MeV);
    analysis_manager->FillNtupleDColumn(12, vertex_z_/um);

    G4int background = -1;
    G4int excited_state_number = -1;
    if (background_) {
        background = 1;
    }
    else {
        excited_state_number = state_number_;
    }

    analysis_manager->FillNtupleIColumn(13, background);
    analysis_manager->FillNtupleIColumn(14, excited_state_number);
    analysis_manager->FillNtupleDColumn(15, proton_energy);
    analysis_manager->FillNtupleDColumn(16, qvalue/MeV);

    analysis_manager->AddNtupleRow();

}
