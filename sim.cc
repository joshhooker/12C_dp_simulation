#include <ctime>

#include <G4HadronicProcessStore.hh>
#include <G4MTRunManager.hh>
#include <G4RadioactiveDecayPhysics.hh>
#include <G4RunManager.hh>
#include <G4StepLimiterPhysics.hh>
#include <G4UIExecutive.hh>
#include <G4UImanager.hh>
#include <G4VisExecutive.hh>
#include <G4VModularPhysicsList.hh>
#include <globals.hh>
#include <QGSP_BERT.hh>
#include <Randomize.hh>

#ifdef G4VIS_USE
#include <G4VisExecutive.hh>
#endif

#ifdef G4UI_USE
#include <G4UIExecutive.hh>
#endif

#include <TChain.h>
#include <TFile.h>
#include <TFileMerger.h>
#include <TH1.h>
#include <TTree.h>
#include <TKey.h>
#include <Riostream.h>

#include "ActionInitialization.hh"
#include "BinaryReactionPhysics.hh"
#include "Calibrations.hh"
#include "DetectorConstruction.hh"
#include "EnergyLoss.hh"
#include "json/json.h"
#include "NonResonantBackgroundPhysics.hh"
#include "NucleonStates.hh"
#include "PhysicsList.hh"

int main(int argc,char** argv) {
    if(argc < 2) {
        std::cout << "Usage: sim config-file" << std::endl;
        return 0;
    }

    // Get the pointer to the User Interface manager
    G4UImanager* ui_manager = G4UImanager::GetUIpointer();

    // Create and read json config
    Json::Value config;
    std::string config_filename = argv[1];
    std::ifstream config_stream(config_filename.c_str());
    config_stream >> config;
    config_stream.close();

    // Parse JSON
    G4String macro_name   = config["macroName"].asString();
    G4bool is_interactive = config["interactive"].asBool();

    // Choose the Random engine
    CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine());
    // Set random seed with system time
    G4long seed = time(nullptr);
    if(argc>2) seed += 473879*atoi(argv[2]);
    CLHEP::HepRandom::setTheSeed(seed);

#ifdef G4MULTITHREADED
    G4int n_threads = 0;
#endif
    for(G4int i = 1; i < argc; i++) {
        G4cout << argv[i] << G4endl;
#ifdef G4MULTITHREADED
        if(G4String(argv[i]) == "-t") {
            n_threads = G4UIcommand::ConvertToInt(argv[i + 1]);
        }
#endif
    }

    // Construct the default run manager
#ifdef G4MULTITHREADED
    auto* run_manager = new G4MTRunManager;
    if(n_threads > 0) {
        run_manager->SetNumberOfThreads(n_threads);
    }
#else
    auto* run_manager = new G4RunManager;
#endif

    // Load calibrations
    Calibrations* calibration = Calibrations::Instance();
    calibration->CheckLoaded();

    // Load Nucleon States
    NucleonStates* states = NucleonStates::Instance();
    states->CheckLoaded();

    // Mandatory user initialization classes
    auto* detector = new DetectorConstruction();
    run_manager->SetUserInitialization(detector);

    // G4VModularPhysicsList* physicsList = new QGSP_BERT(0);
    G4VModularPhysicsList* physics_list = new PhysicsList;
    auto* reaction_physics = new BinaryReactionPhysics();
    // NonResonantBackgroundPhysics* reactionPhysics = new NonResonantBackgroundPhysics();
    physics_list->RegisterPhysics(new G4StepLimiterPhysics());
    physics_list->RegisterPhysics(reaction_physics);
    run_manager->SetUserInitialization(physics_list);
    G4HadronicProcessStore::Instance()->SetVerbose(0);

    // User action initialization
    auto* action_init = new ActionInitialization(detector);
    run_manager->SetUserInitialization(action_init);

    // Initialize Geant4 kernel
    run_manager->Initialize();

    // Generate energy loss tables for beam and ejectile
    calibration->GetTargetProperties(detector->GetTargetMaterial());
    calibration->ReaddEdxTables();

    // Visualization manager construction
    G4VisManager* vis_manager = new G4VisExecutive("Quiet");
    vis_manager->Initialize();

    if (!is_interactive) {
        // Batch mode
        G4String command = "/control/execute ";
        ui_manager->ApplyCommand(command + macro_name);
    }
    else {
        // Interactive mode
        G4UIExecutive* ui = new G4UIExecutive(argc, argv);
        ui_manager->ApplyCommand("/control/execute init_vis.mac");
        ui->SessionStart();
        delete ui;
    }

    // Job termination
    // Free the store: user actions, physics_list and detector_description are
    // owned and deleted by the run manager, so they should not be deleted
    // in the main() program !

    delete vis_manager;
    delete run_manager;

    if (!is_interactive) {
        // Write histograms/trees
        auto *root_file = RootFile::Instance();
        root_file->WriteToFile();

        TFileMerger *fm = new TFileMerger(false);
        fm->OutputFile("sim.root");
        fm->AddFile("simHistograms.root");
        fm->AddFile("simTree.root");
        fm->Merge();

        if (remove("simHistograms.root") != 0) perror("Error deleting simHistograms.root");
        if (remove("simTree.root") != 0) perror("Error deleting simTree.root");
    }

    return 0;
}
