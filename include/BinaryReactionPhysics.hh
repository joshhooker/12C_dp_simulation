#pragma once

#include <map>
#include <string>

#include <G4GenericIon.hh>
#include <G4PhysicsConstructorFactory.hh>
#include <G4PhysicsListHelper.hh>
#include <G4VPhysicsConstructor.hh>
#include <globals.hh>

#include "BinaryReactionProcess.hh"

class BinaryReactionPhysics : public G4VPhysicsConstructor {
public:
 	BinaryReactionPhysics(G4int verbose = 0);
 	BinaryReactionPhysics(const G4String& name);
 	virtual ~BinaryReactionPhysics();

 	virtual void ConstructParticle();
 	virtual void ConstructProcess();

private:
};
