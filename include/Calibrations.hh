#ifndef Calibrations_h
#define Calibrations_h

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <utility>
#include <vector>

#include <G4EmCalculator.hh>
#include <G4IonTable.hh>
#include <G4NistManager.hh>
#include <G4ParticleGun.hh>
#include <G4ParticleTable.hh>
#include <G4SystemOfUnits.hh>
#include <G4Types.hh>
#include <G4UnitsTable.hh>
#include <globals.hh>

#include "EnergyLoss.hh"
#include "json/json.h"

class Calibrations {
public:
    Calibrations();
    static Calibrations* Instance();

    void CheckLoaded() {G4cout << "Loaded Calibrations!" << G4endl;}

    double GetYuAngle(int ring) {return yu_angle_[ring];}
    std::pair<double, double> GetYuBins(int ring) {return std::make_pair(yu_angle_low_bin_[ring], yu_angle_high_bin_[ring]);}

    std::pair<G4int, G4int> GetBeamIon() {return beam_ion_;}
    std::pair<G4int, G4int> GetTargetIon() {return target_ion_;}
    std::pair<G4int, G4int> GetEjectileIon() {return ejectile_ion_;}
    double GetBeamEnergy() {return beam_energy_;}
    std::pair<double, double> GetBeamOffset() {return std::make_pair(beam_offset_x, beam_offset_y);}
    double GetTargetThickness() {return target_thickness_;}
    double GetYuEnergyFWHM() {return yu_energy_fwhm_;}
    double GetYuThreshold() {return yu_threshold_;}
    double GetStateEnergyFWHM() {return state_energy_fwhm_;}
    double GetBackgroundAmount() {return background_amount_;}
    bool GetUseEnergyQvalue() {return use_energy_qvalue_;}

    std::string GetOutputFileName() {return output_file_;}

    void GetTargetProperties(G4Material* material);
    void ReaddEdxTables();

    EnergyLoss* GetBeamEnergyLoss() {return beam_energyloss;}
    EnergyLoss* GetEjectileEnergyLoss() {return ejectile_energyloss;}

    G4double GetBeamMassPDG() {return beam_mass_pdg_;}
    G4double GetBeamMassDBL() {return beam_mass_dbl_;}
    G4double GetEjectileMassPDG() {return ejectile_mass_pdg_;}
    G4double GetEjectileMassDBL() {return ejectile_mass_dbl_;}
    G4double GetHeavyRecoilMassPDG() {return heavy_recoil_mass_pdg_;}
    G4double GetHeavyRecoilMassDBL() {return heavy_recoil_mass_dbl_;}

private:
    static Calibrations* instance_;

    EnergyLoss* beam_energyloss;
    EnergyLoss* ejectile_energyloss;

    G4double beam_mass_pdg_;
    G4double beam_mass_dbl_;
    G4double ejectile_mass_pdg_;
    G4double ejectile_mass_dbl_;
    G4double heavy_recoil_mass_pdg_;
    G4double heavy_recoil_mass_dbl_;

    double yu_angle_low_bin_[16];
    double yu_angle_high_bin_[16];
    double yu_angle_[16];

    void ReadJSON();
    std::pair<G4int, G4int> beam_ion_;
    std::pair<G4int, G4int> target_ion_;
    std::pair<G4int, G4int> ejectile_ion_;
    double beam_energy_;
    double beam_offset_x;
    double beam_offset_y;
    double target_thickness_;
    double yu_energy_fwhm_;
    double yu_threshold_;
    double state_energy_fwhm_;
    double background_amount_;

    bool use_energy_qvalue_;

    std::string output_file_;
};


#endif
