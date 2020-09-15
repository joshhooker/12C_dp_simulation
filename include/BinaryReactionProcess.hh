#ifndef BinaryReactionProcess_h
#define BinaryReactionProcess_h

#include <cmath>

#include <G4Event.hh>
#include <G4IonTable.hh>
#include <G4ParticleTable.hh>
#include <G4PhysicalConstants.hh>
#include <G4RunManager.hh>
#include <G4SystemOfUnits.hh>
#include <G4VDiscreteProcess.hh>
#include <G4UserEventAction.hh>
#include <Randomize.hh>

#include "Analysis.hh"
#include "Calibrations.hh"
#include "DetectorConstruction.hh"
#include "NucleonStates.hh"
#include "RootFile.hh"
#include "TrackingInformation.hh"
#include "TypeDef.hh"

class BinaryReactionProcess : public G4VDiscreteProcess {
public:
 	BinaryReactionProcess(const G4String& process_name = "BinaryReaction");
 	~BinaryReactionProcess();

 	G4double GetMeanFreePath(const G4Track&, G4double, G4ForceCondition*);
 	G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);

 	G4VParticleChange* Decay(const G4Track&, G4int, G4int, G4int, G4int);

  	void StartTracking(G4Track*);

private:
    G4double scattering_energy_;
 	G4double scattering_length_;

    uniform_angle_output UniformCMAngle(uniform_angle_input input);
    uniform_angle_output UniformLabAngle(uniform_angle_input input);
    uniform_angle_output CalculateKinematicsLab(uniform_angle_input input, G4double angle_theta, G4double angle_phi);

    G4double Lorentzian_shoot(G4double mean, G4double gamma);
};

#endif
