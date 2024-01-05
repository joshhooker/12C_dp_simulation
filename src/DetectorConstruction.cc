#include "DetectorConstruction.hh"

DetectorConstruction::DetectorConstruction() :
    G4VUserDetectorConstruction()
{}

DetectorConstruction::~DetectorConstruction() = default;

G4VPhysicalVolume* DetectorConstruction::Construct() {
    ConstructMaterials();

    Calibrations* cal = Calibrations::Instance();

    // Overlaps flag
    G4bool check_overlaps = false;

    G4int nel, natoms, Ze, N;
    G4double a, z, density, fraction_mass;

    // Deuteron isotope
    G4Isotope* deuteron = new G4Isotope("deuteron", Ze=1, N=2, a=2.0141018*g/mole);

    // Deuterium element
    G4Element* deuterium = new G4Element("deuterium", "deuterium", 1);
    deuterium->AddIsotope(deuteron, 1);

    // Solid deuterium
    G4Material* D2 = new G4Material("D2", 0.201*g/cm3, 1, kStateSolid, 4.0*kelvin);
    D2->AddElement(deuterium, 2);

    target_material_ = D2;

    // Silicon Detectors
    G4Material* si_material = G4Material::GetMaterial("G4_Si");

    // Carbon target
    G4Material* c_material = G4Material::GetMaterial("G4_C");

    // Aluminum
    G4Material* al_material = G4Material::GetMaterial("G4_Al");

    // Silicon Dioxide
    G4Material* sio2_material = G4Material::GetMaterial("G4_SILICON_DIOXIDE");

    // Boron
    G4Material* b_material = G4Material::GetMaterial("G4_B");

    // Phosphorus
    G4Material* p_material = G4Material::GetMaterial("G4_P");

    // Silver foil
    G4Material* ag_material = G4Material::GetMaterial("G4_Ag");

    auto* vacuum =
        new G4Material("Vacuum",      //Name as String
            1,             //Atomic Number,  in this case we use 1 for hydrogen
            1.008*g/mole,  //Mass per Mole "Atomic Weight"  1.008*g/mole for Hydoren
            1.e-25*g/cm3,  //Density of Vaccuum  *Cant be Zero, Must be small instead
            kStateGas,     //kStateGas for Gas
            2.73*kelvin,   //Temperatuer for gas
            1.e-25*g/cm3); //Pressure for Vaccum

    // Create vacuum filled world
    G4VSolid* world_solid = new G4Box("worldBox", 0.4*m, 0.4*m, 1.*m);
    world_logical_ = new G4LogicalVolume(world_solid, vacuum, "worldLogical");
    G4VPhysicalVolume* world_physical = new G4PVPlacement(nullptr, G4ThreeVector(), world_logical_, "worldPhysical", nullptr,
                                                          false, 0, check_overlaps);

    // D2 Target
    G4double target_thickness = cal->GetTargetThickness();
    G4VSolid* target_solid = new G4Tubs("targetSolid", 0., 50.*mm, target_thickness/2., 0., 360.*deg);
    target_logical_ = new G4LogicalVolume(target_solid, D2, "targetLogical");
    new G4PVPlacement(nullptr, G4ThreeVector(0., 0., 0.), target_logical_, "targetPhysical", world_logical_,
                      false, 0, check_overlaps);

    // Carbon target instead of D2
    // G4double target_thickness = 8.*um;
    // G4VSolid* target_solid = new G4Tubs("target_solid", 0., 50.*mm, target_thickness/2., 0., 360.*deg);
    // target_logical_ = new G4LogicalVolume(target_solid, c_material, "targetLogical");
    // new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), target_logical_, "targetPhysical", world_logical_,
    //                   false, 0, check_overlaps);

    G4double max_step = target_thickness/100.;
    step_limit_ = new G4UserLimits(max_step);
    target_logical_->SetUserLimits(step_limit_);

    // Ag Backing
    G4double ag_thickness = 4.64*um;
    G4VSolid* ag_solid = new G4Tubs("agSolid", 0., 50.*mm, ag_thickness/2., 0., 360.*deg);
    ag_logical_ = new G4LogicalVolume(ag_solid, ag_material, "agLogical");
    new G4PVPlacement(nullptr, G4ThreeVector(0., 0., target_thickness/2. + ag_thickness/2.), ag_logical_, "agPhysical", world_logical_,
                      false, 0, check_overlaps);

    // Yu Detector
    G4double ring_width = 4.9375*mm;
    G4double inner_radius = 50.*mm;
    G4double yu_thickness = 500.*um;
    G4double yu_distance = -80.8*mm;
    G4ThreeVector yu_vector(0., 0., yu_distance);
    char yu_name[256];
    for(G4int i = 0; i < 16; i++) {
        sprintf(yu_name, "yu%d", i + 1);
        G4VSolid* yu_solid;
        if(i < 13) {
            yu_solid = new G4Tubs(yu_name, inner_radius + i*ring_width, inner_radius + (i + 1)*ring_width, yu_thickness/2., 0., 42.*deg);
        }
        else if(i < 14) {
            yu_solid = new G4Tubs(yu_name, inner_radius + i*ring_width, inner_radius + (i + 1)*ring_width, yu_thickness/2., 3.36*deg, 35.28*deg);
        }
        else if(i < 15) {
            yu_solid = new G4Tubs(yu_name, inner_radius + i*ring_width, inner_radius + (i + 1)*ring_width, yu_thickness/2., 7.72*deg, 26.56*deg);
        }
        else {
            yu_solid = new G4Tubs(yu_name, inner_radius + i*ring_width, inner_radius + (i + 1)*ring_width, yu_thickness/2., 13.50*deg, 15.01*deg);
        }

        G4double angle_rot = 42.;

        sprintf(yu_name, "yuLogical%d", i + 1);
        yu_logical_[i] = new G4LogicalVolume(yu_solid, si_material, yu_name);
        sprintf(yu_name, "yuPhysical%d", i + 1);
        for(G4int j = 0; j < 8; j++) {
            auto rot = new G4RotationMatrix;
            rot->rotateZ((angle_rot + 45.*j - 45.*2)*deg);
            new G4PVPlacement(rot, yu_vector, yu_logical_[i], yu_name, world_logical_, false, j, check_overlaps);
        }
    }

    // Yd Detector
    G4double yd_thickness = 500.*um;
    G4double yd_distance = 86*mm;
    G4ThreeVector yd_vector(0., 0., yd_distance);
    char yd_name[256];
    for(G4int i = 0; i < 16; i++) {
        sprintf(yd_name, "yd%d", i + 1);
        G4VSolid* yd_solid;
        if(i < 13) {
            yd_solid = new G4Tubs(yd_name, inner_radius + i*ring_width, inner_radius + (i + 1)*ring_width, yd_thickness/2., 0., 42.*deg);
        }
        else if(i < 14) {
            yd_solid = new G4Tubs(yd_name, inner_radius + i*ring_width, inner_radius + (i + 1)*ring_width, yd_thickness/2., 3.36*deg, 35.28*deg);
        }
        else if(i < 15) {
            yd_solid = new G4Tubs(yd_name, inner_radius + i*ring_width, inner_radius + (i + 1)*ring_width, yd_thickness/2., 7.72*deg, 26.56*deg);
        }
        else {
            yd_solid = new G4Tubs(yd_name, inner_radius + i*ring_width, inner_radius + (i + 1)*ring_width, yd_thickness/2., 13.50*deg, 15.01*deg);
        }

        G4double angle_rot = 42.;

        sprintf(yd_name, "ydLogical%d", i + 1);
        yd_logical_[i] = new G4LogicalVolume(yd_solid, si_material, yd_name);
        sprintf(yd_name, "ydPhysical%d", i + 1);
        for(G4int j = 0; j < 8; j++) {
            auto rot = new G4RotationMatrix;
            rot->rotateZ((angle_rot + 45.*j)*deg);
            new G4PVPlacement(rot, yd_vector, yd_logical_[i], yd_name, world_logical_, false, j, check_overlaps);
        }
    }

    G4double s3_inner_radius = 11.*mm;
    G4double s3_outer_radius = 35.*mm;
    G4double s3_ring_width = 1.*mm;
    G4double s3_angle_rot = 11.25;

    // Sd1
    G4double sd1_distance = 600.*mm;
    G4double sd1_thickness = 61.*um;

    char sd1_d1_al1_name[256];
    char sd1_d1_sio2_name[256];
    char sd1_d1_al2_name[256];
    char sd1_d1_b_name[256];
    char sd1_d2_p_name[256];
    char sd1_d2_al_name[256];
    char sd_name[256];

    double sd1_d1_al1_thickness = 1.5*um;
    double sd1_d1_sio2_thickness = 3.5*um;
    double sd1_d1_al2_thickness = 0.3*um;
    double sd1_d1_b_thickness = 0.5*um;
    double sd1_d2_p_thickness = 0.5*um;
    double sd1_d2_al_thickness = 0.3*um;

    // Sd1 dead layers D1
    sprintf(sd1_d1_al1_name, "sd1_d1_al1");
    G4VSolid *sd1_d1_al1_solid = new G4Tubs(sd1_d1_al1_name, s3_inner_radius, s3_outer_radius, sd1_d1_al1_thickness/2., 0., 360.*deg);
    sprintf(sd1_d1_al1_name, "sd1_d1_al1Logical");
    sd1_d1_al1_logical_ = new G4LogicalVolume(sd1_d1_al1_solid, al_material, sd1_d1_al1_name);
    sprintf(sd1_d1_al1_name, "sd1_d1_al1Physical");
    new G4PVPlacement(nullptr, G4ThreeVector(0., 0., sd1_distance - sd1_thickness/2. - sd1_d1_b_thickness - sd1_d1_al2_thickness -
        sd1_d1_sio2_thickness - sd1_d1_al1_thickness/2.), sd1_d1_al1_logical_, sd1_d1_al1_name, world_logical_, false, 0, check_overlaps);

    sprintf(sd1_d1_sio2_name, "sd1_d1_sio2");
    G4VSolid *sd1_d1_sio2_solid = new G4Tubs(sd1_d1_sio2_name, s3_inner_radius, s3_outer_radius, sd1_d1_sio2_thickness/2., 0., 360.*deg);
    sprintf(sd1_d1_sio2_name, "sd1_d1_sio2Logical");
    sd1_d1_sio2_logical_ = new G4LogicalVolume(sd1_d1_sio2_solid, sio2_material, sd1_d1_sio2_name);
    sprintf(sd1_d1_sio2_name, "sd1_d1_sio2Physical");
    new G4PVPlacement(nullptr, G4ThreeVector(0., 0., sd1_distance - sd1_thickness/2. - sd1_d1_b_thickness - sd1_d1_al2_thickness -
        sd1_d1_sio2_thickness/2.), sd1_d1_sio2_logical_, sd1_d1_sio2_name, world_logical_, false, 0, check_overlaps);

    sprintf(sd1_d1_al2_name, "sd1_d1_al2");
    G4VSolid *sd1_d1_al2_solid = new G4Tubs(sd1_d1_al2_name, s3_inner_radius, s3_outer_radius, sd1_d1_al2_thickness/2., 0., 360.*deg);
    sprintf(sd1_d1_al2_name, "sd1_d1_al2Logical");
    sd1_d1_al2_logical_ = new G4LogicalVolume(sd1_d1_al2_solid, al_material, sd1_d1_al2_name);
    sprintf(sd1_d1_al2_name, "sd1_d1_al2Physical");
    new G4PVPlacement(nullptr, G4ThreeVector(0., 0., sd1_distance - sd1_thickness/2. - sd1_d1_b_thickness - sd1_d1_al2_thickness/2.),
        sd1_d1_al2_logical_, sd1_d1_al2_name, world_logical_, false, 0, check_overlaps);

    sprintf(sd1_d1_b_name, "sd1_d1_b");
    G4VSolid *sd1_d1_b_solid = new G4Tubs(sd1_d1_b_name, s3_inner_radius, s3_outer_radius, sd1_d1_b_thickness/2., 0., 360.*deg);
    sprintf(sd1_d1_b_name, "sd1_d1_bLogical");
    sd1_d1_b_logical_ = new G4LogicalVolume(sd1_d1_b_solid, b_material, sd1_d1_b_name);
    sprintf(sd1_d1_b_name, "sd1_d1_bPhysical");
    new G4PVPlacement(nullptr, G4ThreeVector(0., 0., sd1_distance - sd1_thickness/2. - sd1_d1_b_thickness/2.),
        sd1_d1_b_logical_, sd1_d1_b_name, world_logical_, false, 0, check_overlaps);

    // Sd1 active Si layer
    for (G4int i = 0; i < 24; i++) {
        sprintf(sd_name, "sd1%d", i + 1);
        G4VSolid *sd_solid = new G4Tubs(sd_name, s3_inner_radius + i * s3_ring_width,
                                        s3_inner_radius + (i + 1) * s3_ring_width, sd1_thickness / 2., 0.,
                                        s3_angle_rot * deg);
        sprintf(sd_name, "sd1Logical%d", i + 1);
        sd1_logical_[i] = new G4LogicalVolume(sd_solid, si_material, sd_name);
        sprintf(sd_name, "sd1Physical%d", i + 1);
        for (G4int j = 0; j < 32; j++) {
            auto rot = new G4RotationMatrix;
            rot->rotateZ((-s3_angle_rot * j - s3_angle_rot * 8) * deg);
            new G4PVPlacement(rot, G4ThreeVector(0., 0., sd1_distance), sd1_logical_[i], sd_name, world_logical_, false,
                              j, check_overlaps);
        }
    }

    // Sd1 dead layers D2
    sprintf(sd1_d2_p_name, "sd1_d2_p");
    G4VSolid *sd1_d2_p_solid = new G4Tubs(sd1_d2_p_name, s3_inner_radius, s3_outer_radius, sd1_d2_p_thickness/2., 0., 360.*deg);
    sprintf(sd1_d2_p_name, "sd1_d2_pLogical");
    sd1_d2_p_logical_ = new G4LogicalVolume(sd1_d2_p_solid, p_material, sd1_d2_p_name);
    sprintf(sd1_d2_p_name, "sd1_d2_pPhysical");
    new G4PVPlacement(nullptr, G4ThreeVector(0., 0., sd1_distance + sd1_thickness/2. + sd1_d2_p_thickness/2.),
        sd1_d2_p_logical_, sd1_d2_p_name, world_logical_, false, 0, check_overlaps);

    sprintf(sd1_d2_al_name, "sd1_d2_al");
    G4VSolid *sd1_d2_al_solid = new G4Tubs(sd1_d2_al_name, s3_inner_radius, s3_outer_radius, sd1_d2_al_thickness/2., 0., 360.*deg);
    sprintf(sd1_d2_al_name, "sd1_d2_alLogical");
    sd1_d2_al_logical_ = new G4LogicalVolume(sd1_d2_al_solid, al_material, sd1_d2_al_name);
    sprintf(sd1_d2_al_name, "sd1_d2_alPhysical");
    new G4PVPlacement(nullptr, G4ThreeVector(0., 0., sd1_distance + sd1_thickness/2. + sd1_d2_p_thickness + sd1_d2_al_thickness/2.),
        sd1_d2_al_logical_, sd1_d2_al_name, world_logical_, false, 0, check_overlaps);

    // Sd2
    G4double sd2_distance = 690.*mm;
    G4double sd2_thickness = 500.*um;

    // Sd2 dead layers D2
    char sd2_d2_p_name[256];
    char sd2_d2_al_name[256];

    double sd2_d2_p_thickness = 0.5*um;
    double sd2_d2_al_thickness = 0.3*um;

    sprintf(sd2_d2_al_name, "sd2_d2_al");
    G4VSolid *sd2_d2_al_solid = new G4Tubs(sd2_d2_al_name, s3_inner_radius, s3_outer_radius, sd2_d2_al_thickness/2., 0., 360.*deg);
    sprintf(sd2_d2_al_name, "sd2_d2_alLogical");
    sd2_d2_al_logical_ = new G4LogicalVolume(sd2_d2_al_solid, al_material, sd2_d2_al_name);
    sprintf(sd2_d2_al_name, "sd2_d2_alPhysical");
    new G4PVPlacement(nullptr, G4ThreeVector(0., 0., sd2_distance - sd2_thickness/2. - sd2_d2_p_thickness - sd2_d2_al_thickness/2.),
        sd2_d2_al_logical_, sd2_d2_al_name, world_logical_, false, 0, check_overlaps);

    sprintf(sd2_d2_p_name, "sd2_d2_p");
    G4VSolid *sd2_d2_p_solid = new G4Tubs(sd2_d2_p_name, s3_inner_radius, s3_outer_radius, sd2_d2_p_thickness/2., 0., 360.*deg);
    sprintf(sd2_d2_p_name, "sd2_d2_pLogical");
    sd2_d2_p_logical_ = new G4LogicalVolume(sd2_d2_p_solid, p_material, sd2_d2_p_name);
    sprintf(sd2_d2_p_name, "sd2_d2_pPhysical");
    new G4PVPlacement(nullptr, G4ThreeVector(0., 0., sd2_distance - sd2_thickness/2. - sd1_d2_p_thickness/2.),
        sd2_d2_p_logical_, sd2_d2_p_name, world_logical_, false, 0, check_overlaps);

    // Sd2 active Si layer
    for (G4int i = 0; i < 24; i++) {
        sprintf(sd_name, "sd2%d", i + 1);
        G4VSolid *sd_solid = new G4Tubs(sd_name, s3_inner_radius + i * s3_ring_width,
                                        s3_inner_radius + (i + 1) * s3_ring_width, sd2_thickness / 2., 0.,
                                        s3_angle_rot * deg);
        sprintf(sd_name, "sd2Logical%d", i + 1);
        sd2_logical_[i] = new G4LogicalVolume(sd_solid, si_material, sd_name);
        sprintf(sd_name, "sd2Physical%d", i + 1);
        for (G4int j = 0; j < 32; j++) {
            auto rot = new G4RotationMatrix;
            rot->rotateZ((-s3_angle_rot * j - s3_angle_rot * 8) * deg);
            new G4PVPlacement(rot, G4ThreeVector(0., 0., sd2_distance), sd2_logical_[i], sd_name, world_logical_, false,
                              j, check_overlaps);
        }
    }

    SetAttributes();

    return world_physical;
}

void DetectorConstruction::ConstructMaterials() {
    G4NistManager* man = G4NistManager::Instance();
    man->FindOrBuildMaterial("G4_C");
    man->FindOrBuildMaterial("G4_Si");
    man->FindOrBuildMaterial("G4_Ag");
    man->FindOrBuildMaterial("G4_Al");
    man->FindOrBuildMaterial("G4_P");
    man->FindOrBuildMaterial("G4_B");
    man->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
}

void DetectorConstruction::ConstructSDandField() {
    G4SDManager* sd_man = G4SDManager::GetSDMpointer();
    G4String sd_name;

    // Yu Detector
    for (G4int i = 0; i < 16; i++) {
        char name_yu[256];
        sprintf(name_yu, "yu%d", i + 1);
        G4VSensitiveDetector* yu_detector = new GenSD(sd_name = name_yu);
        sd_man->AddNewDetector(yu_detector);
        yu_logical_[i]->SetSensitiveDetector(yu_detector);
    }

    // Yd Detector
    for (G4int i = 0; i < 16; i++) {
        char name_yd[256];
        sprintf(name_yd, "yd%d", i + 1);
        G4VSensitiveDetector* yd_detector = new GenSD(sd_name = name_yd);
        sd_man->AddNewDetector(yd_detector);
        yd_logical_[i]->SetSensitiveDetector(yd_detector);
    }

    // Sd1 Detector
    char name_sd[256];
    for (G4int i = 0; i < 24; i++) {
        sprintf(name_sd, "sd1%d", i + 1);
        G4VSensitiveDetector *sd1_detector = new GenSD(sd_name = name_sd);
        sd_man->AddNewDetector(sd1_detector);
        sd1_logical_[i]->SetSensitiveDetector(sd1_detector);
    }

    // Sd2 Detector
    for (G4int i = 0; i < 24; i++) {
        sprintf(name_sd, "sd2%d", i + 1);
        G4VSensitiveDetector *sd2_detector = new GenSD(sd_name = name_sd);
        sd_man->AddNewDetector(sd2_detector);
        sd2_logical_[i]->SetSensitiveDetector(sd2_detector);
    }
}

void DetectorConstruction::SetAttributes() {
    auto* world_attr = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));
    world_attr->SetVisibility(true);
    world_attr->SetForceWireframe(true);
    world_logical_->SetVisAttributes(world_attr);

    auto* target_attr = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));
    target_attr->SetVisibility(true);
    target_attr->SetForceSolid(true);
    target_logical_->SetVisAttributes(target_attr);

    auto* ag_attr = new G4VisAttributes(G4Colour::White());
    ag_attr->SetVisibility(true);
    ag_attr->SetForceSolid(true);
    ag_logical_->SetVisAttributes(ag_attr);
}
