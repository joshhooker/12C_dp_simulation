#pragma once

#include <G4String.hh>
#include <G4UserTrackingAction.hh>
#include <G4VProcess.hh>

#include "BinaryReactionProcess.hh"
#include "EventAction.hh"
#include "NonResonantBackgroundProcess.hh"
#include "TrackingInformation.hh"

class TrackingAction : public G4UserTrackingAction {
public:
    TrackingAction(EventAction*);
    ~TrackingAction();

    void PreUserTrackingAction(const G4Track*);

private:
    G4String name_;
    EventAction* event_action_ { nullptr };
};
