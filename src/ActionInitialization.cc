#include "ActionInitialization.hh"

ActionInitialization::ActionInitialization(DetectorConstruction* detector) :
    G4VUserActionInitialization(), detector_(detector) {

}

ActionInitialization::~ActionInitialization() = default;

void ActionInitialization::BuildForMaster() const {
    SetUserAction(new RunAction(detector_, nullptr));
}

void ActionInitialization::Build() const {
    auto* primary = new PrimaryGeneratorAction();
    SetUserAction(primary);

    SetUserAction(new RunAction(detector_, primary));

    auto event_action = new EventAction();
    SetUserAction(event_action);

    SetUserAction(new TrackingAction(event_action));
}
