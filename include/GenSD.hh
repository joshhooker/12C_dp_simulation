#pragma once

#include <G4SDManager.hh>
#include <G4SystemOfUnits.hh>
#include <G4ThreeVector.hh>
#include <G4Types.hh>
#include <G4VSensitiveDetector.hh>

#include "GenHit.hh"

class GenSD : public G4VSensitiveDetector {
public:
 	GenSD(G4String name);
 	virtual ~GenSD();

 	virtual void Initialize(G4HCofThisEvent* HCE);
 	virtual G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist);

private:
 	GenHitsCollection* hits_collection_ { nullptr };
 	G4int hcid_ { 0 };
};
