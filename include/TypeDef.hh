#ifndef TypeDef_h
#define TypeDef_h

#include <G4ThreeVector.hh>
#include <G4Track.hh>

#include "CubicSpline.hh"

typedef struct theshold_struct {
    G4double energy;
    uint decay_charge;
    uint decay_mass;
} threshold_struct;

typedef struct state_struct {
    G4double energy;
    G4double width;
    uint spin2;
    G4bool parity;
    G4bool use_spin_parity;
    G4double probability;
    G4double cumulative_probability;
    G4double probability_rng;
    G4bool angle_cm;
    G4String angle_file;
    CubicSpline angle_spline;
    G4bool use_angle;
    G4double angle_low;
    G4double angle_high;
} state_struct;

typedef struct isotope_struct {
    std::string name;
    uint charge;
    uint mass;
    std::vector<state_struct> states;
    std::vector<threshold_struct> thresholds;
} isotope_struct;

typedef struct yy_detect_struct {
    G4int det;
    G4int ring;
    G4double energy;
    G4ThreeVector position;
} yyDetect_struct;

typedef struct sd_detect_struct {
    G4int det;
    G4int ring;
    G4double energy;
    G4ThreeVector position;
} sd_detect_struct;

typedef struct excited_state_struct {
    G4int state_number;
    G4double excited_energy;
    G4double excited_width;
    G4bool use_angle;
} excited_state_struct;

typedef struct uniform_angle_input {
    const G4Track& aTrack;
    const G4Step& aStep;
    G4double beam_mass;
    G4double target_mass;
    G4double light_recoil_mass;
    G4double heavy_recoil_mass;
    G4double beam_energy;
    G4double q_value;
} uniform_angle_input;

typedef struct uniform_angle_output {
    G4double energy;
    G4double cm_energy;
    G4double qvalue;
    G4double total_energy;
    G4double cm_light_theta;
    G4double cm_light_phi;
    G4double cm_heavy_theta;
    G4double lab_light_theta;
    G4double lab_light_phi;
    G4double lab_heavy_theta;
    G4ThreeVector lab_light_vector;
    G4ThreeVector lab_heavy_vector;
    G4double lab_light_energy;
    G4double lab_heavy_energy;
    G4bool correct_kinematics;
} uniform_angle_output;

#endif //TypeDef_h
