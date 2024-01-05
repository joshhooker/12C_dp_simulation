#pragma once

#include <ctime>

#include <G4GenericMessenger.hh>
#include <G4IonTable.hh>
#include <G4ParticleGun.hh>
#include <G4ParticleTable.hh>
#include <G4RunManager.hh>
#include <G4SystemOfUnits.hh>
#include <G4VUserPrimaryGeneratorAction.hh>
#include <Randomize.hh>

#include "Calibrations.hh"
#include "EventAction.hh"

class G4ParticleGun;
class G4ParticleDefinition;
class G4GenericMessenger;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {

public:
    PrimaryGeneratorAction();
    virtual ~PrimaryGeneratorAction();

    virtual void GeneratePrimaries(G4Event*);

    G4ParticleGun* GetParticleGun() {return particle_gun_;};

private:
    G4ParticleGun* particle_gun_ { nullptr };
    G4ParticleDefinition* particle_ { nullptr };
    G4GenericMessenger* messenger_ { nullptr };

    void DefineCommands();
};
