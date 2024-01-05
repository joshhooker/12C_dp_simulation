#pragma once

#include <cassert>
#include <cstring>
#include <fstream>
#include <vector>

#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TF2.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>

#include "Calibrations.hh"

class RootFile {
public:
    RootFile() = default;
    static RootFile* Instance();

    void Initialize();
    void WriteToFile();

    TH1F* GetLightAngleCM() {return h_light_angle_cm_;}
    TH1F* GetLightAngleCMDist() {return h_light_angle_cm_dist_;}
    TH1F* GetLightAngleThetaLab() {return h_light_angle_theta_lab_;}
    TH1F* GetLightAngleThetaLabDist() {return h_light_angle_theta_lab_dist_;}
    TH1F* GetLightAnglePhiLab() {return h_light_angle_phi_lab_;}
    TH1F* GetLightAnglePhiLabDist() {return h_light_angle_phi_lab_dist_;}
    TH1F* GetLightEnergy() {return h_light_energy_;}
    TH1F* GetEnergy() {return h_energy_;}
    TH1F* GetVertexZ() {return h_vertex_z_;}

    TH1F* GetExcitationEnergy() {return h_excitation_energy_;}
    TH1F* GetExcitationEnergyDist() {return h_excitation_energy_dist_;}
    TH1F* GetQValue() {return h_qvalue_;}
    TH1F* GetQValueDist() {return h_qvalue_dist_;}
    TH1F* GetCalcQValue() {return h_calc_qvalue_;}
    TH1F* GetCalcQValueBackground() {return h_calc_qvalue_background_;}
    TH1F* GetCalcQValueStates() {return h_calc_qvalue_states_;}
    TH1F* GetCalcQValueFull() {return h_calc_qvalue_full_;}
    TH1F* GetQValueDiff() {return h_qvalue_diff_;}

    TH1F* GetCalcQValue07() {return h_calc_qvalue_07_;}
    TH1F* GetCalcQValue815() {return h_calc_qvalue_815_;}

    TH1F* GetProtonEnergyError() {return h_proton_energy_error_;}
    TH1F* GetProtonAngleError() {return h_proton_angle_error_;}
    TH1F* GetBeamEnergyError() {return h_beam_energy_error_;}
    TH1F* GetQvalueError() {return h_qvalue_error_;}

    TH1F* GetSd1Energy() {return h_sd1_energy_;}
    TH1F* GetSd2Energy() {return h_sd2_energy_;}
    TH1F* GetSdTotalEnergy() {return h_sd_total_energy_;}

    TH1F* GetAngularLabState(uint state) {return h_angular_lab_state_[state];}
    TH1F* GetAngularLabStateDist(uint state) {return h_angular_lab_state_dist_[state];}
    TH1F* GetCalcQvalueState(uint state) {return h_calc_qvalue_state_[state];}

    TH2F* GetQValueVertex() {return h_qvalue_vertex_;}

private:
    static RootFile* instance_;

    TH1F* h_light_angle_cm_ { nullptr };
    TH1F* h_light_angle_cm_dist_ { nullptr };
    TH1F* h_light_angle_theta_lab_ { nullptr };
    TH1F* h_light_angle_theta_lab_dist_ { nullptr };
    TH1F* h_light_angle_phi_lab_ { nullptr };
    TH1F* h_light_angle_phi_lab_dist_ { nullptr };
    TH1F* h_light_energy_ { nullptr };
    TH1F* h_energy_ { nullptr };
    TH1F* h_vertex_z_ { nullptr };

    TH1F* h_excitation_energy_ { nullptr };
    TH1F* h_excitation_energy_dist_ { nullptr };
    TH1F* h_qvalue_ { nullptr };
    TH1F* h_qvalue_dist_ { nullptr };
    TH1F* h_calc_qvalue_ { nullptr };
    TH1F* h_calc_qvalue_background_ { nullptr };
    TH1F* h_calc_qvalue_states_ { nullptr };
    TH1F* h_calc_qvalue_full_ { nullptr };
    TH1F* h_qvalue_diff_ { nullptr };

    TH1F* h_calc_qvalue_07_ { nullptr };
    TH1F* h_calc_qvalue_815_ { nullptr };

    TH1F* h_proton_energy_error_ { nullptr };
    TH1F* h_proton_angle_error_ { nullptr };
    TH1F* h_beam_energy_error_ { nullptr };
    TH1F* h_qvalue_error_ { nullptr };

    TH1F* h_sd1_energy_ { nullptr };
    TH1F* h_sd2_energy_ { nullptr };
    TH1F* h_sd_total_energy_ { nullptr };

    TH1F* h_angular_lab_state_[10] { nullptr };
    TH1F* h_angular_lab_state_dist_[10] { nullptr };
    TH1F* h_calc_qvalue_state_[10] { nullptr };

    TH2F* h_qvalue_vertex_ { nullptr };
};
