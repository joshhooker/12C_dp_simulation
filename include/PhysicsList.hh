#ifndef PhysicsList_hh
#define PhysicsList_hh

#include <G4DecayPhysics.hh>
// #include <G4DeexPrecoParameters.hh>
#include <G4EmExtraPhysics.hh>
#include <G4EmParameters.hh>
#include <G4EmStandardPhysics.hh>
#include <G4HadronElasticPhysics.hh>
#include <G4HadronPhysicsFTFP_BERT.hh>
#include <G4HadronInelasticQBBC.hh>
#include <G4HadronPhysicsINCLXX.hh>
#include <G4IonConstructor.hh>
#include <G4IonElasticPhysics.hh>
#include <G4IonINCLXXPhysics.hh>
#include <G4IonPhysics.hh>
#include <G4LossTableManager.hh>
// #include <G4NuclearLevelData.hh>
#include <G4NuclideTable.hh>
#include <G4NeutronTrackingCut.hh>
#include <G4ParticleTypes.hh>
#include <G4PhysicsListHelper.hh>
// #include <G4Radioactivation.hh>
#include <G4StoppingPhysics.hh>
#include <G4SystemOfUnits.hh>
#include <G4UAtomicDeexcitation.hh>
#include <G4UnitsTable.hh>
#include <G4VModularPhysicsList.hh>
#include <G4VUserPhysicsList.hh>

// particles

#include <G4BosonConstructor.hh>
#include <G4LeptonConstructor.hh>
#include <G4MesonConstructor.hh>
#include <G4BosonConstructor.hh>
#include <G4BaryonConstructor.hh>
#include <G4IonConstructor.hh>
#include <G4ShortLivedConstructor.hh>

#include <globals.hh>

class PhysicsList: public G4VModularPhysicsList {
public:
    PhysicsList();
    ~PhysicsList();

protected:
    // Construct particle and physics
    virtual void ConstructParticle();
    virtual void SetCuts();
};

#endif