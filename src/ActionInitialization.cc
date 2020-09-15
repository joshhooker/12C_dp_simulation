#include "ActionInitialization.hh"

ActionInitialization::ActionInitialization(DetectorConstruction* detector) :
    G4VUserActionInitialization(), detector_(detector) {

}

ActionInitialization::~ActionInitialization() {}

void ActionInitialization::BuildForMaster() const {
    SetUserAction(new RunAction(detector_, NULL));
}

void ActionInitialization::Build() const {
    PrimaryGeneratorAction* primary = new PrimaryGeneratorAction();
    SetUserAction(primary);

    SetUserAction(new RunAction(detector_, primary));

    auto event_action = new EventAction();
    SetUserAction(event_action);

    SetUserAction(new TrackingAction(event_action));
}
