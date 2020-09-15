#include "TrackingInformation.hh"

TrackingInformation::TrackingInformation(G4double energy, G4double cm_energy, G4double cm_light_theta,
                                         G4double lab_light_theta, G4double cm_light_phi, G4double lab_light_phi,
                                         G4double cm_heavy_theta, G4double lab_heavy_theta, G4double light_energy,
                                         G4double heavy_energy, G4ThreeVector vertex, G4double qvalue, G4double excited_energy,
                                         G4ParticleDefinition* light_recoil, G4ParticleDefinition* heavy_recoil,
                                         G4bool background, G4int state_number) {
    energy_ = energy;
    cm_energy_ = cm_energy;
    cm_light_theta_ = cm_light_theta;
    lab_light_theta_ = lab_light_theta;
    cm_light_phi_ = cm_light_phi;
    lab_light_phi_ = lab_light_phi;
    cm_heavy_theta_ = cm_heavy_theta;
    lab_heavy_theta_ = lab_heavy_theta;
    light_energy_ = light_energy;
    heavy_energy_ = heavy_energy;
    vertex_ = vertex;
    qvalue_ = qvalue;
    excited_energy_ = excited_energy;
    light_recoil_ = light_recoil;
    heavy_recoil_ = heavy_recoil;
    background_ = background;
    state_number_ = state_number;
}

TrackingInformation::~TrackingInformation() = default;

void TrackingInformation::Print() const {}
