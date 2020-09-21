#include "TrackingAction.hh"

TrackingAction::TrackingAction(EventAction* eventAction)
    : event_action_(eventAction) {
    name_ = "BinaryReaction";
}

TrackingAction::~TrackingAction() = default;

void TrackingAction::PreUserTrackingAction(const G4Track* track) {
    const G4VProcess* creatorProcess = track->GetCreatorProcess();
    if(!creatorProcess) return;

    if(creatorProcess->GetProcessName() != name_) return;

    if(track->GetTrackID() != 2) return;

    auto* info = (TrackingInformation*) track->GetUserInformation();

    event_action_->SetEnergy(info->GetEnergy());
    event_action_->SetExcitationEnergy(info->GetExcitedEnergy());
    event_action_->SetQValue(info->GetQValue());
    event_action_->SetLightEnergy(info->GetLightEnergy());
    event_action_->SetLightAngleThetaCM(info->GetCMLightTheta());
    event_action_->SetLightAngleThetaLab(info->GetLabLightTheta());
    event_action_->SetLightAnglePhiCM(info->GetCMLightPhi());
    event_action_->SetLightAnglePhiLab(info->GetLabLightPhi());
    event_action_->SetVertexX(info->GetVertex().x()/mm);
    event_action_->SetVertexY(info->GetVertex().y()/mm);
    event_action_->SetVertexZ(info->GetVertex().z()/mm);
    event_action_->SetBackground(info->GetBackground());
    event_action_->SetStateNumber(info->GetStateNumber());
}
