#include "NonResonantBackgroundPhysics.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(NonResonantBackgroundPhysics);

NonResonantBackgroundPhysics::NonResonantBackgroundPhysics(G4int)
: G4VPhysicsConstructor("NonResonantBackgroundPhysics") {}

NonResonantBackgroundPhysics::NonResonantBackgroundPhysics(const G4String& name)
: G4VPhysicsConstructor(name) {}

NonResonantBackgroundPhysics::~NonResonantBackgroundPhysics() = default;

void NonResonantBackgroundPhysics::ConstructParticle() {
    G4GenericIon::GenericIon();
}

void NonResonantBackgroundPhysics::ConstructProcess() {
    auto* reaction_process = new NonResonantBackgroundProcess();
    G4PhysicsListHelper::GetPhysicsListHelper()->RegisterProcess(reaction_process, G4GenericIon::GenericIon());
}

