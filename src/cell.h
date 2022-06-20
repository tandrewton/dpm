#ifndef CELL_H
#define CELL_H

/*

        HEADER FILE FOR CELL CLASS

                -- Inherits from DPM class
                -- For cells with variable force routines and variable simulation boundaries
*/

#include <math.h>
#include <algorithm>
#include <stdexcept>
#include "dpm.h"

// namespace
using namespace std;

class cell;
typedef void (cell::*cellMemFn)(void);

class cell : public dpm {
 protected:
  // specific output objects
  std::ofstream enout;      // output stream for energy data
  std::ofstream stressout;  // output stream for stress data

  // simulation-scale stored quantities
  double simclock;         // accumulator for simulation time
  std::vector<double> VL;  // wall velocity if simulating with wall forces

  // energy/force modifiers
  std::vector<double> cellID;  // 0 = tissue, 1 = cell (can set other numbers for other identities,
                               // which can be used to control the interaction types)
                               //

 public:
  cell(int n, double att1, double att2, int seed)
      : dpm(n, seed) {
    // set values here
    vector<double> zerosNCELLS(NCELLS, 0);
    l1 = att1;
    l2 = att2;
    simclock = 0.0;
    VL = zerosNCELLS;
    cout << "Initializing epi2D object, l1 = " << l1 << ", l2 = " << l2 << '\n';
  }

  // epi cell interactions
  void repulsiveForceUpdateWithWalls();
  // void repulsiveForceUpdateWithoutTopWall();
  void vertexAttractiveForces2D_2();
  void attractiveForceUpdate();
  void repulsiveWithPolarityForceUpdate();
  void attractiveWithPolarityForceUpdate();

  void wallForces(bool top, bool bottom, bool left, bool right, double& forceTop, double& forceBottom, double& forceLeft, double& forceRight, double appliedUniaxialPressure = 0.0);
  void cellPolarityForces(int ci, double k_polarity, std::string direction = "y");

  // File openers
  void openEnergyObject(std::string& str) {
    enout.open(str.c_str());
    if (!enout.is_open()) {
      std::cout << "	ERROR: file could not open " << str << "..."
                << std::endl;
      exit(1);
    } else
      std::cout << "** Opening file " << str << " ..." << std::endl;
  }

  void openStressObject(std::string& str) {
    stressout.open(str.c_str());
    if (!stressout.is_open()) {
      std::cout << "	ERROR: file could not open " << str << "..."
                << std::endl;
      exit(1);
    } else
      std::cout << "** Opening file " << str << " ..." << std::endl;
  }

  // routines
  void vertexCompress2Target2D(dpmMemFn forceCall, double Ftol, double dt0, double phi0Target, double dphi0);
  void simulateDampedWithWalls(dpmMemFn forceCall, double B, double dt0, double duration, double printInterval, bool wallsOn, bool topOn, bool bottomOn, bool leftOn, bool rightOn, double trueStrainRateX = 0.0, double trueStrainRateY = 0.0, double appliedUniaxialPressure = 0.0);

  // printouts
  void printConfiguration2D();
};

#endif