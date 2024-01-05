#pragma once

#include <cassert>
#include <cstdio>
#include <map>
#include <vector>

#include <G4Event.hh>
#include <G4ParticleDefinition.hh>
#include <G4RootAnalysisManager.hh>
#include <G4RunManager.hh>
#include <G4SDManager.hh>
#include <G4SystemOfUnits.hh>
#include <G4ThreeVector.hh>
#include <G4Types.hh>
#include <G4UserEventAction.hh>
#include <Randomize.hh>

#include <TF1.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TRandom3.h>
#include <TROOT.h>
#include <TVector3.h>

#include "Calibrations.hh"
#include "GenHit.hh"
#include "EnergyLoss.hh"
#include "RootFile.hh"
#include "TypeDef.hh"

class EventAction : public G4UserEventAction {

public:
    EventAction();
    ~EventAction();

    void BeginOfEventAction(const G4Event*);
    void EndOfEventAction(const G4Event*);

    std::vector<G4int>& GetYuDetect() {return yu_detect_vec_;}
    std::vector<G4int>& GetYuRing() {return yu_ring_vec_;}
    std::vector<G4double>& GetYuEnergy() {return yu_energy_vec_;}

    void SetEnergy(G4double energy) {energy_ = energy;}
    void SetExcitationEnergy(G4double energy) {excitation_energy_ = energy;}
    void SetQValue(G4double qValue) {qvalue_ = qValue;}
    void SetLightEnergy(G4double energy) {light_energy_ = energy;}
    void SetLightAngleThetaCM(G4double angle) {light_angle_theta_cm_ = angle;}
    void SetLightAngleThetaLab(G4double angle) {light_angle_theta_lab_ = angle;}
    void SetLightAnglePhiCM(G4double angle) {light_angle_phi_cm_ = angle;}
    void SetLightAnglePhiLab(G4double angle) {light_angle_phi_lab_ = angle;}
    void SetVertexX(G4double vertex) {vertex_x_ = vertex;}
    void SetVertexY(G4double vertex) {vertex_y_ = vertex;}
    void SetVertexZ(G4double vertex) {vertex_z_ = vertex;}
    void SetBackground(G4bool background) {background_ = background;}
    void SetStateNumber(G4int state_number) {state_number_ = state_number;}

private:
    G4int yu_hcid_[16] { 0 };
    G4int sd1_hcid_[24] { 0 };
    G4int sd2_hcid_[24] { 0 };

    G4double beam_energy_ { 0. };
    std::vector<G4int> yu_detect_vec_;
    std::vector<G4int> yu_ring_vec_;
    std::vector<G4double> yu_energy_vec_;

    G4int yu_det_ { 0 };
    G4int yu_ring_ { 0 };
    G4double yu_energy_ { 0. };
    G4ThreeVector yu_position_;

    G4int sd1_det_ { 0 };
    G4int sd1_ring_ { 0 };
    G4double sd1_energy_ { 0. };

    G4int sd2_det_ { 0 };
    G4int sd2_ring_ { 0 };
    G4double sd2_energy_ { 0. };

    G4double energy_ { 0. };
    G4double excitation_energy_ { 0. };
    G4double qvalue_ { 0. };
    G4double light_energy_ { 0. };
    G4double light_angle_theta_cm_ { 0. };
    G4double light_angle_theta_lab_ { 0. };
    G4double light_angle_phi_cm_ { 0. };
    G4double light_angle_phi_lab_ { 0. };
    G4double vertex_x_ { 0. };
    G4double vertex_y_ { 0. };
    G4double vertex_z_ { 0. };

    G4bool background_ { false };
    G4int state_number_ { 0 };
};
