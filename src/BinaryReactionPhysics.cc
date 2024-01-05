#include "BinaryReactionPhysics.hh"

G4_DECLARE_PHYSCONSTR_FACTORY(BinaryReactionPhysics);

BinaryReactionPhysics::BinaryReactionPhysics(G4int)
:  G4VPhysicsConstructor("BinaryReactionPhysics") {
}

BinaryReactionPhysics::BinaryReactionPhysics(const G4String& name)
:  G4VPhysicsConstructor(name) {
}

BinaryReactionPhysics::~BinaryReactionPhysics() = default;

void BinaryReactionPhysics::ConstructParticle() {
 	G4GenericIon::GenericIon();
}

void BinaryReactionPhysics::ConstructProcess() {
 	auto* reaction_process = new BinaryReactionProcess();
 	G4PhysicsListHelper::GetPhysicsListHelper()->RegisterProcess(reaction_process, G4GenericIon::GenericIon());
}

