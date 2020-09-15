#ifndef ActionInitialization_h
#define ActionInitialization_h

#include <cstdio>
#include <map>

#include <G4VUserActionInitialization.hh>

#include "DetectorConstruction.hh"
#include "EventAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "TrackingAction.hh"

class DetectorConstruction;
class EventAction;

class ActionInitialization : public G4VUserActionInitialization {
public:
    ActionInitialization(DetectorConstruction*);
    virtual ~ActionInitialization();

    virtual void BuildForMaster() const;
    virtual void Build() const;

private:
    DetectorConstruction* detector_;
};

#endif
