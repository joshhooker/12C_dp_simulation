#ifndef EnergyLoss_h
#define EnergyLoss_h

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

#include "CubicSpline.hh"

class EnergyLoss {
public:
    EnergyLoss();
    EnergyLoss(const char*);
    EnergyLoss(const char*, double);
    EnergyLoss(const char*, bool);
    EnergyLoss(const char*, double , bool);

    void ReadBasicdEdx(const char*);
    void SetDebug(bool);

    void ReadInitParams();
    double CalcRemainder(double, double);
    double AddBack(double, double);
    double CalcRange(double, double);
    void CalcRemainderError(double);
    void AddBackError(double);
    void AddBackHigh(double);

    friend double GetdEdx(const EnergyLoss&, double);

    void UseGL16();
    void UseGL32();
    void UseGL64();
    void UseGL128();
    void UseGL256();
    void UseGL512();
    void UseGL1024();

    void UseCalcRangeDebug();

private:
    bool debug_ = false;
    bool calc_range_debug_ = false;

    void ReadFile(const char*);
    void ReadFileSimple(const char*, double);

    std::vector<double> energy_vec_;
    std::vector<double> dedx_vec_;
    CubicSpline energy_spline_;

    double x16[8], w16[8], x32[16], w32[16], x64[32], w64[32];
    double x128[64], w128[64], x256[128], w256[128];
    double x512[256], w512[256], x1024[512], w1024[512];

    double calc_remainder_err_;
    double add_back_err_;
    double add_back_high_point_;

    bool gl16_, gl32_, gl64_, gl128_, gl256_, gl512_, gl1024_;

    double CalcRangeGL16(double, double);
    double CalcRangeGL32(double, double);
    double CalcRangeGL64(double, double);
    double CalcRangeGL128(double, double);
    double CalcRangeGL256(double, double);
    double CalcRangeGL512(double, double);
    double CalcRangeGL1024(double, double);

    // Variables to read SRIM File
    const int kMaxCharsPerLine_ = 1024000;
    const char* const kDelimiter_ = " ";
    bool before_multiply_ = true;
    std::vector<char*> token;
    double stopping_conversion_power_;
    std::vector<double> stopping_conversion_power_vec_;
    std::vector<std::string> stopping_conversion_unit1_vec_;
    std::vector<std::string> stopping_conversion_unit2_vec_;

    std::string PrintOutput(std::string Output, std::string Color);
};

inline EnergyLoss::EnergyLoss() {
    debug_ = false;
}

inline EnergyLoss::EnergyLoss(const char* srimFile) {
    debug_ = false;
    ReadFile(srimFile);
}

inline EnergyLoss::EnergyLoss(const char* srimFile, bool debugger) {
    debug_ = debugger;
    ReadFile(srimFile);
}

inline EnergyLoss::EnergyLoss(const char* srimFile, double stoppingPower) {
    debug_ = false;
    ReadFileSimple(srimFile, stoppingPower);
}

inline EnergyLoss::EnergyLoss(const char* srimFile, double stoppingPower, bool debugger) {
    debug_ = debugger;
    ReadFileSimple(srimFile, stoppingPower);
}

inline void EnergyLoss::ReadFile(const char* srimFile) {
    std::ifstream inFile(srimFile);
    ASSERT_WITH_MESSAGE(inFile.is_open(), "Cannot find SRIM file!\n");

    if(debug_) std::cout << PrintOutput("DEBUGGING: READING IN SRIM FILE", "red") << std::endl;

    while(!inFile.eof()) {
        char buf[1024000];
        inFile.getline(buf, kMaxCharsPerLine_);
        char* segment = strtok(buf, kDelimiter_);
        token.push_back(segment);

        if(token[0]) { //check if the line is good
            if(atof(token[0]) && before_multiply_) { // check if number

                double energy = atof(token[0]);
                std::string energyUnit = strtok(0, kDelimiter_);
                double dEElec = atof(strtok(0, kDelimiter_));
                double dENuclear  = atof(strtok(0, kDelimiter_));

                if(debug_) {
                    char* buffer[256];
                    sprintf(reinterpret_cast<char *>(buffer), "Energy: %f %s dEElec: %f dENuclear: %f", energy, energyUnit.c_str(), dEElec, dENuclear);
                    std::cout << PrintOutput(reinterpret_cast<const char *>(buffer), "green") << std::endl;
                }

                if(energyUnit.compare("eV") == 0) energy_vec_.push_back(energy*1.e-6);
                else if(energyUnit.compare("keV") == 0) energy_vec_.push_back(energy*1.e-3);
                else if(energyUnit.compare("MeV") == 0) energy_vec_.push_back(energy);
                else if(energyUnit.compare("GeV") == 0) energy_vec_.push_back(energy*1.e3);
                dedx_vec_.push_back(dEElec+dENuclear);
            }

            if(!strcmp(token[0], "Multiply")) {
                before_multiply_ = false;
            }

            if(atof(token[0]) && !before_multiply_) {
                stopping_conversion_power_vec_.push_back(atof(token[0]));
                stopping_conversion_unit1_vec_.push_back(strtok(0, kDelimiter_));
                strtok(0, kDelimiter_);
                stopping_conversion_unit2_vec_.push_back(strtok(0, kDelimiter_));
            }
        }

        if(!token.empty()) {
            token.clear(); // clear token vector
        }
    }

    for(size_t i = 0; i < stopping_conversion_power_vec_.size(); i++) {
        if(stopping_conversion_unit1_vec_[i].compare("MeV") == 0 && stopping_conversion_unit2_vec_[i].compare("mm") == 0) {
            stopping_conversion_power_ = stopping_conversion_power_vec_[i];
            if(debug_) {
                char* buffer[256];
                sprintf(reinterpret_cast<char *>(buffer), "Conversion: %f %s/%s", stopping_conversion_power_, stopping_conversion_unit1_vec_[i].c_str(), stopping_conversion_unit2_vec_[i].c_str());
                std::cout << PrintOutput(reinterpret_cast<const char *>(buffer), "green") << std::endl;
            }
        }
    }
    std::transform(dedx_vec_.begin(), dedx_vec_.end(), dedx_vec_.begin(),
                   std::bind(std::multiplies<double>(), std::placeholders::_1, stopping_conversion_power_));
    energy_spline_.SetPoints(energy_vec_, dedx_vec_);

    gl16_ = gl32_ = gl64_ = gl128_ = gl256_ = gl1024_ = false;
    gl512_ = true;

    calc_remainder_err_ = 1.e-8;
    add_back_err_ = 1.e-8;
    add_back_high_point_ = 250.;

    ReadInitParams();
}

inline void EnergyLoss::ReadFileSimple(const char* srimFile, double stoppingPower) {
    std::ifstream inFile(srimFile);
    ASSERT_WITH_MESSAGE(inFile.is_open(), "Cannot find SRIM file!\n");

    if(debug_) std::cout << PrintOutput("DEBUGGING: READING IN SRIM FILE", "red") << std::endl;

    while(!inFile.eof()) {
        char buf[1024000];
        inFile.getline(buf, kMaxCharsPerLine_);
        char* segment = strtok(buf, kDelimiter_);
        token.push_back(segment);

        if(token[0]) { //check if the line is good
            if(atof(token[0]) && before_multiply_) { // check if number
                double energy = atof(token[0]);
                std::string energyUnit = strtok(0, kDelimiter_);
                double dEElec = atof(strtok(0, kDelimiter_));
                double dENuclear = atof(strtok(0, kDelimiter_));

                if(debug_) {
                    char* buffer[256];
                    sprintf(reinterpret_cast<char *>(buffer), "Energy: %f %s dEElec: %f dENuclear: %f", energy, energyUnit.c_str(), dEElec, dENuclear);
                    std::cout << PrintOutput(reinterpret_cast<const char *>(buffer), "green") << std::endl;
                }

                if(energyUnit.compare("eV") == 0) energy_vec_.push_back(energy*1.e-6);
                else if(energyUnit.compare("keV") == 0) energy_vec_.push_back(energy*1.e-3);
                else if(energyUnit.compare("MeV") == 0) energy_vec_.push_back(energy);
                else if(energyUnit.compare("GeV") == 0) energy_vec_.push_back(energy*1.e3);
                dedx_vec_.push_back(dEElec+dENuclear);
            }
        }

        if(!token.empty()) {
            token.clear(); // clear token vector
        }
    }

    std::transform(dedx_vec_.begin(), dedx_vec_.end(), dedx_vec_.begin(),
                   std::bind(std::multiplies<double>(), std::placeholders::_1, stoppingPower));
    energy_spline_.SetPoints(energy_vec_, dedx_vec_);

    gl16_ = gl32_ = gl64_ = gl128_ = gl256_ = gl1024_ = false;
    gl512_ = true;

    calc_remainder_err_ = 1.e-8;
    add_back_err_ = 1.e-8;
    add_back_high_point_ = 250.;

    ReadInitParams();
}

inline void EnergyLoss::ReadBasicdEdx(const char* inputFile) {
    std::ifstream inFile(inputFile);
    ASSERT_WITH_MESSAGE(inFile.is_open(), "Cannot find input file!\n");

    double energy, dEdx;
    while(inFile >> energy >> dEdx) {
        energy_vec_.push_back(energy);
        dedx_vec_.push_back(dEdx);
    }

    energy_spline_.SetPoints(energy_vec_, dedx_vec_);

    gl16_ = gl32_ = gl64_ = gl128_ = gl256_ = gl1024_ = false;
    gl512_ = true;

    calc_remainder_err_ = 1.e-8;
    add_back_err_ = 1.e-8;
    add_back_high_point_ = 250.;

    ReadInitParams();
}

inline void EnergyLoss::SetDebug(bool flag) {
    debug_ = flag;
}

inline double EnergyLoss::CalcRemainder(double initialEnergy, double distance) {
    if(distance == 0.) return initialEnergy;

    distance = fabs(distance);

    if(initialEnergy < calc_remainder_err_) return 0.;

    if(debug_) std::cout << PrintOutput("DEBUGGING: CalcRemainder", "red") << std::endl;

    distance = fabs(distance);

    double maxRange = CalcRange(initialEnergy, 0.);
    if(debug_) std::cout << PrintOutput("Max Range of particle in the material (mm): ", "green") << maxRange << std::endl;

    if(distance > maxRange) return 0.;

    double lowEnergy = 0.;
    double highEnergy = initialEnergy;
    double guessEnergy = (highEnergy + lowEnergy)/2.0;

    double range = CalcRange(initialEnergy, guessEnergy);
    while(fabs(range - distance) > calc_remainder_err_) {
        if(debug_) {
            char* buffer[256];
            sprintf(reinterpret_cast<char *>(buffer), "Low Energy: %f; Guess Energy: %f; High Energy: %f Range: %f; Distance: %f",
                    lowEnergy, guessEnergy, highEnergy, range, distance);
            std::cout << PrintOutput(reinterpret_cast<const char *>(buffer), "green") << std::endl;
        }
        if(range > distance) {
            lowEnergy = guessEnergy;
            guessEnergy = (lowEnergy+highEnergy)/2.0;
        } else {
            highEnergy = guessEnergy;
            guessEnergy = (lowEnergy+highEnergy)/2.0;
        }
        range = CalcRange(initialEnergy, guessEnergy);
    }
    return guessEnergy;
}

inline double EnergyLoss::AddBack(double finalEnergy, double distance) {
    long counter = 0;
    if(distance == 0.) return finalEnergy;

    distance = fabs(distance);

    if(debug_) std::cout << PrintOutput("DEBUGGING: AddBack", "red") << std::endl;

    double lowEnergy = finalEnergy;
    double highEnergy = add_back_high_point_;
    double guessEnergy = (highEnergy + lowEnergy)/2.0;

    double range = CalcRange(guessEnergy, finalEnergy);
    while(fabs(range - distance) > add_back_err_) {
        if(debug_) {
            char* buffer[256];
            sprintf(reinterpret_cast<char *>(buffer), "Low Energy: %f; Guess Energy: %f; High Energy: %f Range: %f; Distance: %f",
                    lowEnergy, guessEnergy, highEnergy, range, distance);
            std::cout << PrintOutput(reinterpret_cast<const char *>(buffer), "green") << std::endl;
        }
        if(range < distance) {
            lowEnergy = guessEnergy;
            guessEnergy = (lowEnergy + highEnergy)/2.0;
        } else {
            highEnergy = guessEnergy;
            guessEnergy = (lowEnergy + highEnergy)/2.0;
        }
        range = CalcRange(guessEnergy, finalEnergy);

        if(guessEnergy > (add_back_high_point_-0.05)) {
            printf("Error: EnergyLoss::AddBack above starting high guess!\n");
            printf("Use the AddBackHigh to set the starting high guess higher!\n");
            return guessEnergy;
        }

        if(counter > 10000) {
            printf("There was a problem: Counter for AddBack > 10,000!\n");
            printf("Final Energy: %f Distance: %f\n", finalEnergy, distance);
            return guessEnergy;
        }

        counter++;
    }
    return guessEnergy;
}

inline double EnergyLoss::CalcRange(double initialEnergy, double remainder) {
    if(initialEnergy < 0.) return 0.;
    if(remainder < 0.) return 0.;

    if(remainder > initialEnergy) return 0;

    if(initialEnergy == remainder) return 0;
    double distance;

    if(gl16_) distance = CalcRangeGL16(remainder, initialEnergy);
    else if(gl32_) distance = CalcRangeGL32(remainder, initialEnergy);
    else if(gl64_) distance = CalcRangeGL64(remainder, initialEnergy);
    else if(gl128_) distance = CalcRangeGL128(remainder, initialEnergy);
    else if(gl256_) distance = CalcRangeGL256(remainder, initialEnergy);
    else if(gl512_) distance = CalcRangeGL512(remainder, initialEnergy);
    else if(gl1024_) distance = CalcRangeGL1024(remainder, initialEnergy);
    else distance = CalcRangeGL512(remainder, initialEnergy);

    return distance;
}

inline void EnergyLoss::CalcRemainderError(double err) {
    calc_remainder_err_ = err;
}

inline void EnergyLoss::AddBackError(double err) {
    add_back_err_ = err;
}

inline void EnergyLoss::AddBackHigh(double high) {
    add_back_high_point_ = high;
}

inline void EnergyLoss::UseGL16() {
    gl32_ = gl64_ = gl128_ = gl256_ = gl512_ = gl1024_ = false;
    gl16_ = true;
}

inline void EnergyLoss::UseGL32() {
    gl16_ = gl64_ = gl128_ = gl256_ = gl512_ = gl1024_ = false;
    gl32_ = true;
}

inline void EnergyLoss::UseGL64() {
    gl16_ = gl32_ = gl128_ = gl256_ = gl512_ = gl1024_ = false;
    gl64_ = true;
}

inline void EnergyLoss::UseGL128() {
    gl16_ = gl32_ = gl64_ = gl256_ = gl512_ = gl1024_ = false;
    gl128_ = true;
}

inline void EnergyLoss::UseGL256() {
    gl16_ = gl32_ = gl64_ = gl128_ = gl512_ = gl1024_ = false;
    gl256_ = true;
}

inline void EnergyLoss::UseGL512() {
    gl16_ = gl32_ = gl64_ = gl128_ = gl256_ = gl1024_ = false;
    gl512_ = true;
}

inline void EnergyLoss::UseGL1024() {
    gl16_ = gl32_ = gl64_ = gl128_ = gl256_ = gl512_ = false;
    gl1024_ = true;
}

inline double EnergyLoss::CalcRangeGL16(double a, double b) {
    double t0 = (a + b)/2.0;
    double dt = (b - a)/2.0;
    double result = 0.;
    if(calc_range_debug_) std::cout << PrintOutput("DEBUGGING: CalcRange", "red") << std::endl;
    for(int i = 7; i > -1; i--) {
        if(calc_range_debug_) {
            char* buffer[256];
            sprintf(reinterpret_cast<char *>(buffer), "i: %d; w16[i]: %f; -x16[i]: %f energySpline(t0 + dt*(-x16[i])): %f",
                    i, w16[i], -x16[i], energy_spline_(t0 + dt*(-x16[i])));
            std::cout << PrintOutput(reinterpret_cast<const char *>(buffer), "green") << std::endl;
        }
        result += w16[i]*(1./energy_spline_(t0 + dt*(-x16[i])));
    }
    for(int i = 0; i < 8; i++) {
        if(calc_range_debug_) {
            char* buffer[256];
            sprintf(reinterpret_cast<char *>(buffer), "i: %d; w16[i]: %f; x16[i]: %f energySpline(t0 + dt*(x16[i])): %f",
                    i, w16[i], x16[i], energy_spline_(t0 + dt*(x16[i])));
            std::cout << PrintOutput(reinterpret_cast<const char *>(buffer), "green") << std::endl;
        }
        result += w16[i]*(1./energy_spline_(t0 + dt*x16[i]));
    }
    return result*dt;
}

inline double EnergyLoss::CalcRangeGL32(double a, double b) {
    double t0 = (a + b)/2.0;
    double dt = (b - a)/2.0;
    double result = 0.;
    for(int i = 15; i > -1; i--) {
        if(calc_range_debug_) {
            char* buffer[256];
            sprintf(reinterpret_cast<char *>(buffer), "i: %d; w32[i]: %f; -x32[i]: %f energySpline(t0 + dt*(-x32[i])): %f",
                    i, w32[i], -x32[i], energy_spline_(t0 + dt*(-x32[i])));
            std::cout << PrintOutput(reinterpret_cast<const char *>(buffer), "green") << std::endl;
        }
        result += w32[i]*(1./energy_spline_(t0 + dt*(-x32[i])));
    }
    for(int i = 0; i < 16; i++) {
        if(calc_range_debug_) {
            char* buffer[256];
            sprintf(reinterpret_cast<char *>(buffer), "i: %d; w32[i]: %f; x32[i]: %f energySpline(t0 + dt*(x32[i])): %f",
                    i, w32[i], x32[i], energy_spline_(t0 + dt*(x32[i])));
            std::cout << PrintOutput(reinterpret_cast<const char *>(buffer), "green") << std::endl;
        }
        result += w32[i]*(1./energy_spline_(t0 + dt*x32[i]));
    }
    return result*dt;
}

inline double EnergyLoss::CalcRangeGL64(double a, double b) {
    double t0 = (a + b)/2.0;
    double dt = (b - a)/2.0;
    double result = 0.;
    for(int i = 31; i > -1; i--) {
        if(calc_range_debug_) {
            char* buffer[256];
            sprintf(reinterpret_cast<char *>(buffer), "i: %d; w64[i]: %f; -x64[i]: %f energySpline(t0 + dt*(-x64[i])): %f",
                    i, w64[i], -x64[i], energy_spline_(t0 + dt*(-x64[i])));
            std::cout << PrintOutput(reinterpret_cast<const char *>(buffer), "green") << std::endl;
        }
        result += w64[i]*(1./energy_spline_(t0 + dt*(-x64[i])));
    }
    for(int i = 0; i < 32; i++) {
        if(calc_range_debug_) {
            char* buffer[256];
            sprintf(reinterpret_cast<char *>(buffer), "i: %d; w64[i]: %f; x64[i]: %f energySpline(t0 + dt*(x64[i])): %f",
                    i, w64[i], x64[i], energy_spline_(t0 + dt*(x64[i])));
            std::cout << PrintOutput(reinterpret_cast<const char *>(buffer), "green") << std::endl;
        }
        result += w64[i]*(1./energy_spline_(t0 + dt*x64[i]));
    }
    return result*dt;
}

inline double EnergyLoss::CalcRangeGL128(double a, double b) {
    double t0 = (a + b)/2.0;
    double dt = (b - a)/2.0;
    double result = 0.;
    for(int i = 63; i > -1; i--) {
        if(calc_range_debug_) {
            char* buffer[256];
            sprintf(reinterpret_cast<char *>(buffer), "i: %d; w128[i]: %f; -x128[i]: %f energySpline(t0 + dt*(-x128[i])): %f",
                    i, w128[i], -x128[i], energy_spline_(t0 + dt*(-x128[i])));
            std::cout << PrintOutput(reinterpret_cast<const char *>(buffer), "green") << std::endl;
        }
        result += w128[i]*(1./energy_spline_(t0 + dt*(-x128[i])));
    }
    for(int i = 0; i < 64; i++) {
        if(calc_range_debug_) {
            char* buffer[256];
            sprintf(reinterpret_cast<char *>(buffer), "i: %d; w128[i]: %f; x128[i]: %f energySpline(t0 + dt*(x128[i])): %f",
                    i, w128[i], x128[i], energy_spline_(t0 + dt*(x128[i])));
            std::cout << PrintOutput(reinterpret_cast<const char *>(buffer), "green") << std::endl;
        }
        result += w128[i]*(1./energy_spline_(t0 + dt*x128[i]));
    }
    return result*dt;
}

inline double EnergyLoss::CalcRangeGL256(double a, double b) {
    double t0 = (a + b)/2.0;
    double dt = (b - a)/2.0;
    double result = 0.;
    for(int i = 127; i > -1; i--) {
        if(calc_range_debug_) {
            char* buffer[256];
            sprintf(reinterpret_cast<char *>(buffer), "i: %d; w256[i]: %f; -x256[i]: %f energySpline(t0 + dt*(-x256[i])): %f",
                    i, w256[i], -x256[i], energy_spline_(t0 + dt*(-x256[i])));
            std::cout << PrintOutput(reinterpret_cast<const char *>(buffer), "green") << std::endl;
        }
        result += w256[i]*(1./energy_spline_(t0 + dt*(-x256[i])));
    }
    for (int i = 0; i < 128; i++) {
        if(calc_range_debug_) {
            char* buffer[256];
            sprintf(reinterpret_cast<char *>(buffer), "i: %d; w256[i]: %f; x256[i]: %f energySpline(t0 + dt*(x256[i])): %f",
                    i, w256[i], x256[i], energy_spline_(t0 + dt*(x256[i])));
            std::cout << PrintOutput(reinterpret_cast<const char *>(buffer), "green") << std::endl;
        }
        result += w256[i]*(1./energy_spline_(t0 + dt*x256[i]));
    }
    return result*dt;
}

inline double EnergyLoss::CalcRangeGL512(double a, double b) {
    double t0 = (a + b)/2.0;
    double dt = (b - a)/2.0;
    double result = 0.;
    for(int i = 255; i > -1; i--) {
        if(calc_range_debug_) {
            char* buffer[256];
            sprintf(reinterpret_cast<char *>(buffer), "i: %d; w512[i]: %f; -x512[i]: %f energySpline(t0 + dt*(-x512[i])): %f",
                    i, w512[i], -x512[i], energy_spline_(t0 + dt*(-x512[i])));
            std::cout << PrintOutput(reinterpret_cast<const char *>(buffer), "green") << std::endl;
        }
        result += w512[i]*(1./energy_spline_(t0 + dt*(-x512[i])));
    }
    for(int i = 0; i < 256; i++) {
        if(calc_range_debug_) {
            char* buffer[256];
            sprintf(reinterpret_cast<char *>(buffer), "i: %d; w512[i]: %f; x512[i]: %f energySpline(t0 + dt*(x512[i])): %f",
                    i, w512[i], x512[i], energy_spline_(t0 + dt*(x512[i])));
            std::cout << PrintOutput(reinterpret_cast<const char *>(buffer), "green") << std::endl;
        }
        result += w512[i]*(1./energy_spline_(t0 + dt*x512[i]));
    }
    return result*dt;
}

inline double EnergyLoss::CalcRangeGL1024(double a, double b) {
    double t0 = (a + b)/2.0;
    double dt = (b - a)/2.0;
    double result = 0.;
    for(int i = 511; i > -1; i--) {
        if(calc_range_debug_) {
            char* buffer[256];
            sprintf(reinterpret_cast<char *>(buffer), "i: %d; w1024[i]: %f; -x1024[i]: %f energySpline(t0 + dt*(-x1024[i])): %f",
                    i, w1024[i], -x1024[i], energy_spline_(t0 + dt*(-x1024[i])));
            std::cout << PrintOutput(reinterpret_cast<const char *>(buffer), "green") << std::endl;
        }
        result += w1024[i]*(1./energy_spline_(t0 + dt*(-x1024[i])));
    }
    for(int i = 0; i < 512; i++) {
        if(calc_range_debug_) {
            char* buffer[256];
            sprintf(reinterpret_cast<char *>(buffer), "i: %d; w1024[i]: %f; x1024[i]: %f energySpline(t0 + dt*(x1024[i])): %f",
                    i, w1024[i], x1024[i], energy_spline_(t0 + dt*(x1024[i])));
            std::cout << PrintOutput(reinterpret_cast<const char *>(buffer), "green") << std::endl;
        }
        result += w1024[i]*(1./energy_spline_(t0 + dt*x1024[i]));
    }
    return result * dt;
}

inline void EnergyLoss::ReadInitParams() {
    std::ifstream inFile("EnergyLoss.dat", std::ios::in | std::ifstream::binary);
    inFile.read(reinterpret_cast<char*> (x16), sizeof(x16));
    inFile.read(reinterpret_cast<char*> (w16), sizeof(w16));
    inFile.read(reinterpret_cast<char*> (x32), sizeof(x32));
    inFile.read(reinterpret_cast<char*> (w32), sizeof(w32));
    inFile.read(reinterpret_cast<char*> (x64), sizeof(x64));
    inFile.read(reinterpret_cast<char*> (w64), sizeof(w64));
    inFile.read(reinterpret_cast<char*> (x128), sizeof(x128));
    inFile.read(reinterpret_cast<char*> (w128), sizeof(w128));
    inFile.read(reinterpret_cast<char*> (x256), sizeof(x256));
    inFile.read(reinterpret_cast<char*> (w256), sizeof(w256));
    inFile.read(reinterpret_cast<char*> (x512), sizeof(x512));
    inFile.read(reinterpret_cast<char*> (w512), sizeof(w512));
    inFile.read(reinterpret_cast<char*> (x1024), sizeof(x1024));
    inFile.read(reinterpret_cast<char*> (w1024), sizeof(w1024));
    inFile.close();
}

inline std::string EnergyLoss::PrintOutput(std::string Output, std::string Color) {
    int ColorCode = 0;
    if(Color.compare("red") == 0){
        ColorCode = 31;
    } else if(Color.compare("green") == 0){
        ColorCode = 32;
    } else if(Color.compare("yellow") == 0){
        ColorCode = 33;
    } else if(Color.compare("blue") == 0) {
        ColorCode = 34;
    } else if(Color.compare("magenta") == 0) {
        ColorCode = 35;
    } else if(Color.compare("cyan") == 0) {
        ColorCode = 36;
    } else {
        return Output;
    }
    char buffer[256];
    sprintf(buffer, "\033[1;%dm%s\033[0m", ColorCode, Output.c_str());
    return buffer;
}

inline void EnergyLoss::UseCalcRangeDebug() {
    calc_range_debug_ = true;
}

inline double GetdEdx(const EnergyLoss& elClass, double energy) {
    return elClass.energy_spline_(energy);
}

#endif