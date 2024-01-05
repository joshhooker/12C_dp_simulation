#include "PhysicsList.hh"

PhysicsList::PhysicsList() : G4VModularPhysicsList() {
    G4int verb = 0;
    SetVerboseLevel(verb);

    // Add new units for radioactive decays
    const G4double minute = 60*second;
    const G4double hour   = 60*minute;
    const G4double day    = 24*hour;
    const G4double year   = 365*day;

    new G4UnitDefinition("minute", "min", "Time", minute);
    new G4UnitDefinition("hour",   "h",   "Time", hour);
    new G4UnitDefinition("day",    "d",   "Time", day);
    new G4UnitDefinition("year",   "y",   "Time", year);

    new G4UnitDefinition( "millielectronVolt", "meV", "Energy", 1.e-3*eV);

    // Mandatory for G4NuclideTable
    // G4NuclideTable::GetInstance()->SetThresholdOfHalfLife(0.1*picosecond);
    // G4NuclideTable::GetInstance()->SetLevelTolerance(1.0*eV);

    // EM Physics
    RegisterPhysics(new G4EmStandardPhysics(verb));
    G4EmParameters* param = G4EmParameters::Instance();
    // param->SetAugerCascade(true);
    // param->SetStepFunction(1., 1*CLHEP::mm);
    // param->SetStepFunctionMuHad(1., 1*CLHEP::mm);

    // Decay
    RegisterPhysics(new G4DecayPhysics(verb));

    // Hadron Elastic scattering
    // RegisterPhysics(new G4HadronElasticPhysics(verb));

    // Hadron Inelastic physics
    // RegisterPhysics(new G4HadronPhysicsFTFP_BERT(verb));

    // Ion Elastic Scattering replaced by NonResonantBackgroundProcess
    // RegisterPhysics( new G4IonElasticPhysics(verb));

    // Ion Inelastic physics
    // RegisterPhysics(new G4IonPhysics(verb));

    // Gamma-Nuclear Physics
    auto* gnuc = new G4EmExtraPhysics(verb);
    // gnuc->ElectroNuclear(false);
    gnuc->MuonNuclear(false);
    RegisterPhysics(gnuc);
}

PhysicsList::~PhysicsList() = default;

void PhysicsList::ConstructParticle() {
    G4BosonConstructor  boson_constructor;
    boson_constructor.ConstructParticle();

    G4LeptonConstructor lepton_constructor;
    lepton_constructor.ConstructParticle();

    G4MesonConstructor meson_constructor;
    meson_constructor.ConstructParticle();

    G4BaryonConstructor baryon_constructor;
    baryon_constructor.ConstructParticle();

    G4IonConstructor ion_constructor;
    ion_constructor.ConstructParticle();

    G4ShortLivedConstructor short_lived_constructor;
    short_lived_constructor.ConstructParticle();
}

void PhysicsList::SetCuts() {
    SetCutValue(0*mm, "proton");
    SetCutValue(10*km, "e-");
    SetCutValue(10*km, "e+");
    SetCutValue(10*km, "gamma");
}
