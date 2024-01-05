# 12C(d, p) Simulation of IRIS

GEANT4 Simulation for the 12C(d, p) reaction with IRIS

## Getting Started

These instructions will get you a copy of the project up and running on your machine.

### Prerequisites

Required Packages:

```
GEANT4 (Working on v11.1)
ROOT
CMake
```

### Build and Run

To build and run the simulation with the above packages installed, create and move to a build directory

```
mkdir build && cd build
```

Use CMake to then build and make the simulation executable

```
cmake .. && make
```

To run the simulation, use the created executable with a configuration file (a sample is provided)

```
./sim config.json -t N
```
where N is the number of threads to run for multi-threading.

### Configuration Files

There are three configuration files provided needed to run the simulation: config.json, nuclear_states.json and run.mac.

#### run.mac

The current example run.mac file looks like this:
```
/run/printProgress 1000
/run/beamOn 50000
```

This will run 50,000 events and print the progress every 1,000 events.

#### config.json

The current example config.json file looks like this:
```
{
  "interactive"      : false,
  "macroName"        : "run.mac",
  "outputFile"       : "sim",
  "beamIon[Z,A]"     : [6, 12],
  "targetIon[Z,A]"   : [1, 2],
  "ejectileIon[Z,A]" : [1, 1],
  "beamEnergy"       : 112.75,
  "beamOffset"       : [0.0, 0.0],
  "targetThickness"  : 43.0,
  "yuEnergyFWHM"     : 0.035,
  "yuThreshold"      : 0.2,
  "use13BeStates"    : false,
  "stateEnergyFWHM"  : 0.1,
  "backgroundAmount" : 0.725,
  "energyQvalue"     : true
}
```

An explanation of each configuration parameter is described below:

- "interactive": When this is set to true, it will load up a GUI to visualize the detector setup and geometry. This can be only set to true or false.
- "macroName": This is macro file to run automatically when "interactive" is false. It's parameters are described above.
- "outputFile": This is the name of the output file that will be generate (with a .root extension).
- "beamIon[Z,A]": This is the Z and A of the beam ion.
- "targetIon[Z,A]": This is the Z and A of the target ion.
- "ejectileIon[Z,A]": This is the Z and A of the ejectile ion (normally the light ion).
- "beamEnergy": Energy of the beam ion (in MeV). This is after the ionization chamber and before the target.
- "beamOffset": Position offset of the beam (in mm).
- "targetThickness": Target thickness of D2 target (in um).
- "yuEnergyFWHM": FWHM resolution of the YYU detector (in MeV).
- "yuThreshold": Threshold of the YYU detector (in MeV).
- "use13BeStates": Not used?
- "stateEnergyFWHM": Not used?
- "backgroundAmount": Amount of events generate that are part of the non-resonant background. Should be between 0 and 1.
- "energyQvalue": If true, energies in nuclear_states.json are in Q-value. If false, energies in nuclear_states.json are in excitation energy.

#### nuclear_states.json

The current example nuclear_states.json file is the following:
```
{
    "Isotopes": [
        {
            "Name": "13Be",
            "ZA": [4, 13],
            "States": [
                {"Energy": -2.20, "Probability": 0.00},
                {"Energy": -2.82, "Probability": 0.213},
                {"Energy": -4.34, "Probability": 0.525}
            ],
            "Thresholds": [
                {"Energy": -1.0, "Decay": [4, 12]}
            ]
        }
    ]
}
```

A general idea of the nuclear_states.json is to input various states to be populated by some probability for a given nucleus. In the 12Be(d, p) simulation, the nucleus we are populating is 13Be.
As "energyQvalue" is set to true in config.json, the energies in this configuration are the Q-value and not the exictation energy. The probabilities that are assigned for each state are automatically normalized
and can range from any value (granted it's positive).

Only the "Energy" for each state is required.

nuclear_states.json has several parameters that can be used:
- "Probability": Probability that the state will be generated in relation to the other states.
- "Spin-parity": Not working.
- "CMAngle": If using angular distribution, should be set to false. True is not working.
- "AngleFile": Points to an angular distribution file. File should have two columns: lab angle and cross section (or some kind of probability).


Another example of nuclear_states.json:
```
{
    "Isotopes": [
        {
            "Name": "13Be",
            "ZA": [4, 13],
            "States": [
                {"Energy": -2.82, "Probability": 0.213, "CMAngle": false, "AngleFile": "../s_wave.dat"},
                {"Energy": -4.34, "Probability": 0.525, "CMAngle": false, "AngleFile": "../d_wave.dat"}
            ],
            "Thresholds": [
                {"Energy": -1.0, "Decay": [4, 12]}
            ]
        }
    ]
}
```