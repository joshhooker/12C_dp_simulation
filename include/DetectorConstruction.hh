#ifndef DetectorConstruction_h
#define DetectorConstruction_h

#include <algorithm>
#include <iostream>
#include <string>

#include <G4Box.hh>
#include <G4ChordFinder.hh>
#include <G4ClassicalRK4.hh>
#include <G4Colour.hh>
#include <G4EqMagElectricField.hh>
#include <G4ExtrudedSolid.hh>
#include <G4FieldManager.hh>
#include <G4GenericMessenger.hh>
#include <G4IntersectionSolid.hh>
#include <G4LogicalVolume.hh>
#include <G4MagIntegratorDriver.hh>
#include <G4NistManager.hh>
#include <G4PVPlacement.hh>
#include <G4RotationMatrix.hh>
#include <G4RunManager.hh>
#include <G4SDManager.hh>
#include <G4SubtractionSolid.hh>
#include <G4SystemOfUnits.hh>
#include <G4TransportationManager.hh>
#include <G4Tubs.hh>
#include <G4TwoVector.hh>
#include <G4UniformElectricField.hh>
#include <G4UserLimits.hh>
#include <G4VisAttributes.hh>
#include <G4VPhysicalVolume.hh>
#include <G4VUserDetectorConstruction.hh>
#include <globals.hh>

#include "Calibrations.hh"
#include "GenSD.hh"

class DetectorConstruction : public G4VUserDetectorConstruction {

public:
    DetectorConstruction();
    virtual ~DetectorConstruction();

    virtual G4VPhysicalVolume* Construct();
    virtual void ConstructSDandField();

    G4LogicalVolume* GetWorldVolume() const {return world_logical_;}
    G4LogicalVolume* GetTargetVolume() const {return target_logical_;}

    G4Material* GetTargetMaterial() {return target_material_;}

private:
    void ConstructMaterials();
    void SetAttributes();

    G4LogicalVolume* world_logical_;
    G4LogicalVolume* target_logical_;
    G4LogicalVolume* ag_logical_;

    G4LogicalVolume* yu_logical_[16];
    G4LogicalVolume* yd_logical_[16];

    G4LogicalVolume* sd1_d1_al1_logical_;
    G4LogicalVolume* sd1_d1_sio2_logical_;
    G4LogicalVolume* sd1_d1_al2_logical_;
    G4LogicalVolume* sd1_d1_b_logical_;
    G4LogicalVolume* sd1_d2_p_logical_;
    G4LogicalVolume* sd1_d2_al_logical_;

    G4LogicalVolume* sd2_d2_p_logical_;
    G4LogicalVolume* sd2_d2_al_logical_;

    G4LogicalVolume* sd1_logical_[24];
    G4LogicalVolume* sd2_logical_[24];

    G4UserLimits* step_limit_;

    G4Material* target_material_;
};

#endif
