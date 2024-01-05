#pragma once

#include <G4ParticleDefinition.hh>
#include <G4ThreeVector.hh>
#include <G4VUserTrackInformation.hh>

class TrackingInformation : public G4VUserTrackInformation {

public:
    TrackingInformation(G4double, G4double, G4double, G4double, G4double,
                        G4double, G4double, G4double, G4double, G4double,
                        G4ThreeVector, G4double, G4double, G4ParticleDefinition*,
                        G4ParticleDefinition*, G4bool, G4int);
    ~TrackingInformation();

    void Print() const;

    G4double GetEnergy() const {return energy_;}
    G4double GetCMEnergy() const {return cm_energy_;}
    G4double GetCMLightTheta() const {return cm_light_theta_;}
    G4double GetLabLightTheta() const {return lab_light_theta_;}
    G4double GetCMLightPhi() const {return cm_light_phi_;}
    G4double GetLabLightPhi() const {return lab_light_phi_;}
    G4double GetCMHeavyTheta() const {return cm_heavy_theta_;}
    G4double GetLabHeavyTheta() const {return lab_heavy_theta_;}
    G4double GetLightEnergy() const {return light_energy_;}
    G4double GetHeavyEnergy() const {return heavy_energy_;}
    G4ThreeVector GetVertex() const {return vertex_;}
    G4double GetQValue() const {return qvalue_;}
    G4double GetExcitedEnergy() const {return excited_energy_;}
    G4ParticleDefinition* GetLightRecoil() const {return light_recoil_;}
    G4ParticleDefinition* GetHeavyRecoil() const {return heavy_recoil_;}
    G4bool GetBackground() const {return background_;}
    G4int GetStateNumber() const {return state_number_;}

    void SetEnergy(G4double energy) {energy_ = energy;}
    void SetCMEnergy(G4double cmEnergy) {cm_energy_ = cmEnergy;}
    void SetCMLightTheta(G4double theta) {cm_light_theta_ = theta;}
    void SetLabLightTheta(G4double theta) {lab_light_theta_ = theta;}
    void SetCMLightPhi(G4double phi) {cm_light_phi_ = phi;}
    void SetLabLightPhi(G4double phi) {lab_light_phi_ = phi;}
    void SetCMHeavyTheta(G4double angle) {cm_heavy_theta_ = angle;}
    void SetLabHeavyTheta(G4double energy) {lab_heavy_theta_ = energy;}
    void SetLightEnergy(G4double energy) {light_energy_ = energy;}
    void SetHeavyEnergy(G4double energy) {heavy_energy_ = energy;}
    void SetVertex(G4ThreeVector v) {vertex_ = v;}
    void SetQValue(G4double qValue) {qvalue_ = qValue;}
    void SetExcitedEnergy(G4double excitedEnergy) {excited_energy_ = excitedEnergy;}
    void SetLightRecoil(G4ParticleDefinition* particle) {light_recoil_ = particle;}
    void SetHeavyRecoil(G4ParticleDefinition* particle) {heavy_recoil_ = particle;}

private:
    G4double energy_ { 0. };
    G4double cm_energy_ { 0. };
    G4double cm_light_theta_ { 0. };
    G4double lab_light_theta_ { 0. };
    G4double cm_light_phi_ { 0. };
    G4double lab_light_phi_ { 0. };
    G4double cm_heavy_theta_ { 0. };
    G4double lab_heavy_theta_ { 0. };
    G4double light_energy_ { 0. };
    G4double heavy_energy_ { 0. };
    G4ThreeVector vertex_;
    G4double qvalue_ { 0. };
    G4double excited_energy_;
    G4ParticleDefinition* light_recoil_ { nullptr };
    G4ParticleDefinition* heavy_recoil_ { nullptr };
    G4bool background_ { false };
    G4int state_number_ { 0 };
};
