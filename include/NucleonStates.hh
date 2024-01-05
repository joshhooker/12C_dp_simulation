#pragma once

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <sys/stat.h>
#include <vector>

#include <G4Types.hh>
#include <globals.hh>
#include <Randomize.hh>

#include "CubicSpline.hh"
#include "json/json.h"
#include "TypeDef.hh"

class NucleonStates {
public:
 	NucleonStates();
 	static NucleonStates* Instance();

 	void CheckLoaded() {G4cout << "Loaded NucleonStates!" << G4endl;}

 	void ReadJSON();
    std::pair<G4bool, G4double> GetAngularDistribution(uint, uint, uint);
 	excited_state_struct GetExcitedLevel(uint, uint);
 	std::vector<theshold_struct> GetThresholds(uint, uint);

private:
 	static NucleonStates* instance_;

 	std::map<uint, std::map<uint, isotope_struct> > nucleons_;
};
