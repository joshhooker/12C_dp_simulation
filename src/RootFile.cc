#include "RootFile.hh"

RootFile* RootFile::instance_ = nullptr;

RootFile* RootFile::Instance() {
    if(!instance_) {
        instance_ = new RootFile();
    }
    return instance_;
}

void RootFile::Initialize() {
    h_light_angle_cm_ = new TH1F("lightAngleCM", "Light Recoil Angle in CM Measured; Angle (deg); Counts", 900, 0, 180);
    h_light_angle_cm_dist_ = new TH1F("lightAngleCMDist", "Light Recoil Angle in CM Distribution; Angle (deg); Counts", 900, 0, 180);
    h_light_angle_theta_lab_ = new TH1F("lightAngleThetaLab", "Light Recoil Angle Theta (Polar) Measured; Angle (deg); Counts", 900, 0, 180);
    h_light_angle_theta_lab_dist_ = new TH1F("lightAngleThetaLabDist", "Light Recoil Angle Theta (Polar) Distribution; Angle (deg); Counts", 900, 0, 180);
    h_light_angle_phi_lab_ = new TH1F("lightAnglePhiLab", "Light Recoil Angle Phi (Azimuthal) Measured; Angle (deg); Counts", 900, 0, 180);
    h_light_angle_phi_lab_dist_ = new TH1F("lightAnglePhiLabDist", "Light Recoil Angle Phi (Azimuthal) Distribution; Angle (deg); Counts", 900, 0, 180);
    h_light_energy_ = new TH1F("lightEnergy", "Light Recoil Energy Measured; Energy (MeV); Counts", 750, 0, 150);
    h_energy_ = new TH1F("energy", "Beam Energy at Reaction Measured; Energy (MeV); Counts", 500, 100, 130);
    h_vertex_z_ = new TH1F("vertexZ", "Vertex on Z axis Measured; Vertex (um); Counts", 600, -60, 60);

    h_excitation_energy_ = new TH1F("excitationEnergy", "Generated Excitation Energy Distribution Measured; Excitation Energy (MeV); Counts", 100, 0, 6);
    h_excitation_energy_dist_ = new TH1F("excitationEnergyDist", "Generated Excitation Energy Distribution; Excitation Energy (MeV); Counts", 100, 0, 6);
    h_qvalue_ = new TH1F("qValue", "Generated Q-Value Distribution Measured; Q-Value (MeV); Counts", 150, -7.0, 4.0);
    h_qvalue_dist_ = new TH1F("qValueDist", "Generated Q-Value Distribution; Q-Value (MeV); Counts", 150, -7.0, 4.0);

    h_calc_qvalue_ = new TH1F("calcQValue", "Calculated Q-Value Spectrum; Q-Value (MeV); Counts", 150, -7.0, 4.0);
    h_calc_qvalue_->GetXaxis()->CenterTitle(); h_calc_qvalue_->GetYaxis()->CenterTitle(); h_calc_qvalue_->Sumw2();
    h_calc_qvalue_background_ = new TH1F("calcQValueBackground", "Calculated Q-Value Spectrum Background", 150, -7.0, 4.0);
    h_calc_qvalue_background_->GetXaxis()->CenterTitle(); h_calc_qvalue_background_->GetYaxis()->CenterTitle(); h_calc_qvalue_background_->Sumw2();
    h_calc_qvalue_states_ = new TH1F("calcQValueStates", "Calculated Q-Value Spectrum States", 150, -7.0, 4.0);
    h_calc_qvalue_states_->GetXaxis()->CenterTitle(); h_calc_qvalue_states_->GetYaxis()->CenterTitle(); h_calc_qvalue_states_->Sumw2();
    h_calc_qvalue_full_ = new TH1F("calcQValueFull", "Calculated Q-Value Spectrum States", 200, -7, 2);
    h_qvalue_diff_ = new TH1F("qValueDiff", "Calculated Q-Value - Generated Q-Value; Q-Value Difference (MeV); Counts", 100, -3, 3);

    h_calc_qvalue_07_ = new TH1F("calcQValue07", "Calculated Q-Value Spectrum Rings 0-7; Q-Value (MeV); Counts", 150, -7.0, 4.0);
    h_calc_qvalue_815_ = new TH1F("calcQValue815", "Calculated Q-Value Spectrum Rings 8-15; Q-Value (MeV); Counts", 150, -7.0, 4.0);

    h_proton_energy_error_ = new TH1F("protonEnergyError", "Calculated Proton Energy at Reaction (Center) - Actual Proton Energy; Energy Difference (MeV); Counts", 500, -5., 5.);
    h_proton_angle_error_  = new TH1F("protonAngleError", "Calculated Proton Angle at Reaction (Center) - Actual Proton Angle; Angle Difference (deg); Counts", 500, -10., 10.);
    h_beam_energy_error_   = new TH1F("beamEnergyError", "Calculated Beam Energy at Reaction (Center) - Actual Beam Energy; Energy Difference (MeV); Counts", 500, -10., 10.);
    h_qvalue_error_        = new TH1F("qvalueError", "Calculated Qvalue at Reaction (Center) - Actual Qvalue; Qvalue Difference (MeV); Counts", 500, -5., 5.);

    h_sd1_energy_ = new TH1F("sd1Energy", "sd1Energy; Energy (MeV); Counts", 500, 0, 50);
    h_sd2_energy_ = new TH1F("sd2Energy", "sd2Energy; Energy (MeV); Counts", 500, 0, 120);
    h_sd_total_energy_ = new TH1F("sdTotalEnergy", "sdTotalEnergy; Energy (MeV); Counts", 500, 0, 120);

    h_qvalue_vertex_ = new TH2F("qvalue_vertex", "Q-value vs Vertex; Q-value (MeV); Vertex (um)", 500, -5, -2, 500, -50, 50);

    for (int i = 0; i < 10; i++) {
        h_angular_lab_state_[i] = new TH1F(Form("angular_lab_state%d", i), Form("Measured Angular Distribution in Lab Frame for State #%d; Angle (deg); Counts", i), 360, 0, 180);
        h_angular_lab_state_dist_[i] = new TH1F(Form("angular_lab_state_dist%d", i), Form("Angular Distribution in Lab Frame for State #%d Distribution; Angle (deg); Counts", i), 360, 0, 180);
	    h_calc_qvalue_state_[i] = new TH1F(Form("calcQValueState%d", i), Form("Calculated Q-Value Spectrum State #%d", i), 100, -3, -2);
    }
}

void RootFile::WriteToFile() {
    auto *cal = Calibrations::Instance();
    TString output_file = Form("%sHistograms.root", cal->GetOutputFileName().c_str());

    auto *out = new TFile(output_file, "recreate");

    h_light_angle_cm_->Write();
    h_light_angle_cm_dist_->Write();
    h_light_angle_theta_lab_->Write();
    h_light_angle_theta_lab_dist_->Write();
    h_light_angle_phi_lab_->Write();
    h_light_angle_phi_lab_dist_->Write();
    h_light_energy_->Write();
    h_energy_->Write();
    h_vertex_z_->Write();

    h_excitation_energy_->Write();
    h_excitation_energy_dist_->Write();
    h_qvalue_->Write();
    h_qvalue_dist_->Write();
    h_calc_qvalue_->Write();
    h_calc_qvalue_background_->Write();
    h_calc_qvalue_states_->Write();
    h_calc_qvalue_full_->Write();
    h_qvalue_diff_->Write();

    h_calc_qvalue_07_->Write();
    h_calc_qvalue_815_->Write();

    h_proton_energy_error_->Write();
    h_proton_angle_error_->Write();
    h_beam_energy_error_->Write();
    h_qvalue_error_->Write();

    h_sd1_energy_->Write();
    h_sd2_energy_->Write();
    h_sd_total_energy_->Write();

    h_qvalue_vertex_->Write();

    for (int i = 0; i < 10; i++) {
        h_angular_lab_state_[i]->Write();
        h_angular_lab_state_dist_[i]->Write();
	h_calc_qvalue_state_[i]->Write();
    }

    out->Close();
}
