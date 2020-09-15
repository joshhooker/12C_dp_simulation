#include "RunAction.hh"

RunAction::RunAction(DetectorConstruction* detector, PrimaryGeneratorAction* primary) :
    G4UserRunAction(), detector_(detector), primary_(primary) {
    // set printing event number per each event
    G4RunManager::GetRunManager()->SetPrintProgress(1);

    // Create analysis manager
    auto analysisManager = G4AnalysisManager::Instance();
    G4cout << "Using " << analysisManager->GetType() << G4endl;

    analysisManager->SetVerboseLevel(1);
    analysisManager->SetNtupleMerging(true);
    analysisManager->SetFileName("simTree");
}

RunAction::~RunAction() {
    delete G4AnalysisManager::Instance();
}

G4Run* RunAction::GenerateRun() {return (new RunData);}

void RunAction::BeginOfRunAction(const G4Run* run) {
    auto analysis_manager = G4AnalysisManager::Instance();

    const EventAction* kEventAction = dynamic_cast<const EventAction*>(G4RunManager::GetRunManager()->GetUserEventAction());
    EventAction* event_action = const_cast<EventAction*>(kEventAction);

    if(G4Threading::G4GetThreadId() == -1) {
        auto *root_file = RootFile::Instance();
        root_file->Initialize();
    }

    auto analysisManager = G4AnalysisManager::Instance();
    analysisManager->CreateNtuple("simData", "simulation data");

    analysisManager->CreateNtupleIColumn("yu_det");    // columnID = 0
    analysisManager->CreateNtupleIColumn("yu_ring");   // columnID = 1
    analysisManager->CreateNtupleDColumn("yu_energy"); // columnID = 2

    analysisManager->CreateNtupleIColumn("sd1_det");      // columnID = 3
    analysisManager->CreateNtupleIColumn("sd1_ring");     // columnID = 4
    analysisManager->CreateNtupleDColumn("sd1_energy");   // columnID = 5

    analysisManager->CreateNtupleIColumn("sd2_det");      // columnID = 6
    analysisManager->CreateNtupleIColumn("sd2_ring");     // columnID = 7
    analysisManager->CreateNtupleDColumn("sd2_energy");   // columnID = 8

    analysisManager->CreateNtupleDColumn("actual_beam_energy");   // columnID = 9
    analysisManager->CreateNtupleDColumn("actual_proton_energy"); // columnID = 10
    analysisManager->CreateNtupleDColumn("actual_q_value");       // columnID = 11
    analysisManager->CreateNtupleDColumn("actual_vertex_z");      // columnID = 12

    analysisManager->CreateNtupleIColumn("background");             // columnID = 13
    analysisManager->CreateNtupleIColumn("excited_state_number");   // columnID = 14
    analysisManager->CreateNtupleDColumn("measured_proton_energy"); // columnID = 15
    analysisManager->CreateNtupleDColumn("measured_q_value");       // columnID = 16

    analysisManager->FinishNtuple();

    analysisManager->OpenFile();
}

void RunAction::EndOfRunAction(const G4Run*) {
    if(G4Threading::G4GetThreadId() == -1) {
        G4cout << "End of Run" << G4endl;
        auto *root_file = RootFile::Instance();
        G4cout << "Writing to sim.root...";
        root_file->WriteToFile();
        G4cout << " done!" << G4endl;
    }

    // save histograms & ntuple
    auto analysisManager = G4AnalysisManager::Instance();
    analysisManager->Write();
    analysisManager->CloseFile();
}
