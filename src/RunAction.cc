#include "RunAction.hh"

RunAction::RunAction(DetectorConstruction* detector, PrimaryGeneratorAction* primary) :
    G4UserRunAction(), detector_(detector), primary_(primary) {
    // set printing event number per each event
    G4RunManager::GetRunManager()->SetPrintProgress(1);

    // Create analysis manager
    auto analysis_manager = G4RootAnalysisManager::Instance();
    G4cout << "Using " << analysis_manager->GetType() << G4endl;

    analysis_manager->SetVerboseLevel(1);
    analysis_manager->SetNtupleMerging(true);
    analysis_manager->SetFileName("simTree");
}

RunAction::~RunAction() {
    delete G4RootAnalysisManager::Instance();
}

G4Run* RunAction::GenerateRun() {return (new RunData);}

void RunAction::BeginOfRunAction(const G4Run* run) {
    const EventAction* kEventAction = dynamic_cast<const EventAction*>(G4RunManager::GetRunManager()->GetUserEventAction());
    EventAction* event_action = const_cast<EventAction*>(kEventAction);

    if(G4Threading::G4GetThreadId() == -1) {
        auto *root_file = RootFile::Instance();
        root_file->Initialize();
    }

    auto analysis_manager = G4RootAnalysisManager::Instance();
    analysis_manager->CreateNtuple("simData", "simulation data");

    analysis_manager->CreateNtupleIColumn("yu_det");    // columnID = 0
    analysis_manager->CreateNtupleIColumn("yu_ring");   // columnID = 1
    analysis_manager->CreateNtupleDColumn("yu_energy"); // columnID = 2

    analysis_manager->CreateNtupleIColumn("sd1_det");      // columnID = 3
    analysis_manager->CreateNtupleIColumn("sd1_ring");     // columnID = 4
    analysis_manager->CreateNtupleDColumn("sd1_energy");   // columnID = 5

    analysis_manager->CreateNtupleIColumn("sd2_det");      // columnID = 6
    analysis_manager->CreateNtupleIColumn("sd2_ring");     // columnID = 7
    analysis_manager->CreateNtupleDColumn("sd2_energy");   // columnID = 8

    analysis_manager->CreateNtupleDColumn("actual_beam_energy");   // columnID = 9
    analysis_manager->CreateNtupleDColumn("actual_proton_energy"); // columnID = 10
    analysis_manager->CreateNtupleDColumn("actual_q_value");       // columnID = 11
    analysis_manager->CreateNtupleDColumn("actual_vertex_z");      // columnID = 12

    analysis_manager->CreateNtupleIColumn("background");             // columnID = 13
    analysis_manager->CreateNtupleIColumn("excited_state_number");   // columnID = 14
    analysis_manager->CreateNtupleDColumn("measured_proton_energy"); // columnID = 15
    analysis_manager->CreateNtupleDColumn("measured_q_value");       // columnID = 16

    analysis_manager->FinishNtuple();

    analysis_manager->OpenFile();
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
    auto analysis_manager = G4RootAnalysisManager::Instance();
    analysis_manager->Write();
    analysis_manager->CloseFile();
}
