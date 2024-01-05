#pragma once

#include <G4EmCalculator.hh>
#include <G4ParticleGun.hh>
#include <G4ProcessManager.hh>
#include <G4RootAnalysisManager.hh>
#include <G4SystemOfUnits.hh>
#include <G4Threading.hh>
#include <G4UnitsTable.hh>
#include <G4UserRunAction.hh>

#include "DetectorConstruction.hh"
#include "EventAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RootFile.hh"
#include "RunData.hh"

class DetectorConstruction;
class PrimaryGeneratorAction;

class RunAction : public G4UserRunAction {
public:
 	RunAction(DetectorConstruction*, PrimaryGeneratorAction*);
 	~RunAction();

  	virtual G4Run* GenerateRun();

 	void BeginOfRunAction(const G4Run*);
 	void EndOfRunAction(const G4Run*);

private:
 	DetectorConstruction* detector_ { nullptr };
 	PrimaryGeneratorAction* primary_ { nullptr };
};
