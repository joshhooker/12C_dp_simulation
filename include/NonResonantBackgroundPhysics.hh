#pragma once

#include <map>
#include <string>

#include <G4GenericIon.hh>
#include <G4PhysicsListHelper.hh>
#include <G4VPhysicsConstructor.hh>
#include <globals.hh>

#include "NonResonantBackgroundProcess.hh"

class NonResonantBackgroundPhysics : public G4VPhysicsConstructor {
public:
    NonResonantBackgroundPhysics(G4int verbose = 0);
    NonResonantBackgroundPhysics(const G4String& name);
    virtual ~NonResonantBackgroundPhysics();

    virtual void ConstructParticle();
    virtual void ConstructProcess();

private:
};
