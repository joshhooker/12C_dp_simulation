#pragma once

#include <G4Run.hh>
#include <G4RunManager.hh>
#include <G4UnitsTable.hh>
#include <globals.hh>

class RunData : public G4Run {
public:
    RunData();
    virtual ~RunData() = default;

private:
};
