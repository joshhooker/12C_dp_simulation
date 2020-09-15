#ifndef RunData_h
#define RunData_h

#include <G4Run.hh>
#include <G4RunManager.hh>
#include <G4UnitsTable.hh>
#include <globals.hh>

#include "Analysis.hh"

class RunData : public G4Run {
public:
 	RunData();
 	virtual ~RunData();

private:
};

#endif