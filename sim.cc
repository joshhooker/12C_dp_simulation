#include <ctime>

#include <G4HadronicProcessStore.hh>
#include <G4MTRunManager.hh>
#include <G4RadioactiveDecayPhysics.hh>
#include <G4RunManager.hh>
#include <G4StepLimiterPhysics.hh>
#include <G4UImanager.hh>
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
    G4UImanager* UImanager = G4UImanager::GetUIpointer();

    // Create and read json config
    Json::Value config;
    std::string configFileName = argv[1];
    std::ifstream configStream(configFileName.c_str());
    configStream >> config;
    configStream.close();

    // Parse JSON
    G4String macroName   = config["macroName"].asString();
    G4bool is_interactive = config["interactive"].asBool();

    // Choose the Random engine
    CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine());
    // Set random seed with system time
    G4long seed = time(NULL);
    if(argc>2) seed += 473879*atoi(argv[2]);
    CLHEP::HepRandom::setTheSeed(seed);

#ifdef G4MULTITHREADED
    G4int nThreads = 0;
#endif
    for(G4int i = 1; i < argc; i++) {
        G4cout << argv[i] << G4endl;
#ifdef G4MULTITHREADED
        if(G4String(argv[i]) == "-t") {
            nThreads = G4UIcommand::ConvertToInt(argv[i + 1]);
        }
#endif
    }

    // Construct the default run manager
#ifdef G4MULTITHREADED
    G4MTRunManager* runManager = new G4MTRunManager;
    if(nThreads > 0) {
        runManager->SetNumberOfThreads(nThreads);
    }
#else
    G4RunManager* runManager = new G4RunManager;
#endif

    // Load calibrations
    Calibrations* calibration = Calibrations::Instance();
    calibration->CheckLoaded();

    // Load Nucleon States
    NucleonStates* states = NucleonStates::Instance();
    states->CheckLoaded();

    // Mandatory user initialization classes
    DetectorConstruction* detector = new DetectorConstruction();
    runManager->SetUserInitialization(detector);

    // G4VModularPhysicsList* physicsList = new QGSP_BERT(0);
    G4VModularPhysicsList* physicsList = new PhysicsList;
    BinaryReactionPhysics* reactionPhysics = new BinaryReactionPhysics();
    // NonResonantBackgroundPhysics* reactionPhysics = new NonResonantBackgroundPhysics();
    physicsList->RegisterPhysics(new G4StepLimiterPhysics());
    physicsList->RegisterPhysics(reactionPhysics);
    runManager->SetUserInitialization(physicsList);
    G4HadronicProcessStore::Instance()->SetVerbose(0);

    // User action initialization
    ActionInitialization* actionInit = new ActionInitialization(detector);
    runManager->SetUserInitialization(actionInit);

    // Initialize Geant4 kernel
    runManager->Initialize();

    // Generate energy loss tables for beam and ejectile
    calibration->GetTargetProperties(detector->GetTargetMaterial());
    calibration->ReaddEdxTables();

#ifdef G4VIS_USE
    // Visualization manager construction
    // G4VisManager* visManager = new G4VisExecutive();
    // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
    G4VisManager* visManager = new G4VisExecutive("Quiet");
    visManager->Initialize();
#endif

    if(!is_interactive) {
        // execute an argument macro file if exists
        G4String command = "/control/execute ";
        G4String filename = macroName;
        UImanager->ApplyCommand(command + filename);
    }
    else {
        // start interactive session
#ifdef G4UI_USE
        G4UIExecutive* ui = new G4UIExecutive(argc, argv);
#ifdef G4VIS_USE
        UImanager->ApplyCommand("/control/execute init_vis.mac");
#else
        UImanager->ApplyCommand("/control/execute init.mac");
#endif
        if(ui->IsGUI()) {
            UImanager->ApplyCommand("/control/execute gui.mac");
            ui->SessionStart();
            delete ui;
        }
#endif
    }

    // Job termination
    // Free the store: user actions, physics_list and detector_description are
    // owned and deleted by the run manager, so they should not be deleted
    // in the main() program !

#ifdef G4VIS_USE
    delete visManager;
#endif
    delete runManager;

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
