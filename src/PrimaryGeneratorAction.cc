#include "PrimaryGeneratorAction.hh"

PrimaryGeneratorAction::PrimaryGeneratorAction() :
    G4VUserPrimaryGeneratorAction() , particle_(0) {

    G4int n_particle = 1;
    particle_gun_  = new G4ParticleGun(n_particle);

    DefineCommands();
}

PrimaryGeneratorAction::~PrimaryGeneratorAction() {
    delete particle_gun_;
    delete messenger_;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* g4Event) {

    G4IonTable* particle_table = G4IonTable::GetIonTable();
    G4String particle_name;
    G4double kin_e;

    Calibrations* cal = Calibrations::Instance();

    std::pair<G4int, G4int> beam_ion = cal->GetBeamIon();
    particle_ = particle_table->GetIon(beam_ion.first, beam_ion.second, 0.);
    kin_e = cal->GetBeamEnergy();

    std::pair<double, double> beam_offset = cal->GetBeamOffset();

    G4double beam_x_pos = G4RandGauss::shoot(beam_offset.first, 3./2.355);
    G4double beam_y_pos = G4RandGauss::shoot(beam_offset.second, 3./2.355);
    // G4double beam_x_pos = beam_offset.first;
    // G4double beam_y_pos = beam_offset.second;
    particle_gun_->SetParticlePosition(G4ThreeVector(beam_x_pos, beam_y_pos, -50.0*cm));
    particle_gun_->SetParticleMomentumDirection(G4ThreeVector(0, 0, 1));

    G4double sigma = kin_e*0.0195;
    kin_e = G4RandGauss::shoot(kin_e, sigma/2.355);
    particle_gun_->SetParticleEnergy(kin_e);

    particle_gun_->SetParticleDefinition(particle_);
    particle_gun_->GeneratePrimaryVertex(g4Event);
}

void PrimaryGeneratorAction::DefineCommands() {
    messenger_ = new G4GenericMessenger(this,
                                        "/generator/",
                                        "Primary generator control");

}
