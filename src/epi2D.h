#ifndef EPI2D_H
#define EPI2D_H

/*

	HEADER FILE FOR epi CLASS

		-- Inherits from DPM class 
		-- For collections of epithelial cells
		-- Incorporates ?
		-- ONLY FOR 2D

*/

#include <algorithm>
#include <stdexcept>
#include "dpm.h"

// namespace
using namespace std;

class epi2D;
typedef void (epi2D::*epi2DMemFn)(void);

class epi2D : public dpm {
 protected:
  // bending energy per vertex
  // NOTE: will need to add different Hessian computation
  std::vector<double> kbi;

  // vertex-vertex contact network
  std::vector<bool> gij;

  // vv contacts per cell
  std::vector<int> z;

  // adhesion strength
  //double att;

  // directors (direction vectors) for activity
  std::vector<double> psi;

  // rotational diffusion constant
  double Dr0;

  // specific output objects
  std::ofstream enout;
  std::ofstream stressout;

  // simulation time keeper (accumulates elapsed simulation time during MD routines)
  double simclock;

 public:
  // constructor and destructor
  epi2D(int n, double att1, double att2, double Dr, int seed)
      : dpm(n, seed) {
    z.resize(n);
    //att = attraction;
    l1 = att1;
    l2 = att2;
    Dr0 = Dr;
    vector<double> temp(NCELLS, 0.0);
    psi = temp;
    for (int ci = 0; ci < NCELLS; ci++) {
      //psi.at(ci) = 2 * PI * (drand48() - 0.5);  // [-pi,pi) coordinates
      //psi.at(ci) = 2 * PI * drand48();  // [0, 2pi) coordinates
      psi.at(ci) = PI / 2 * (ci % 2) - PI / 2 * ((ci + 1) % 2);  //should be : up if odd, down if even
    }
    cout << "Initializing epi2D object, l1 = " << l1 << ", l2 = " << l2 << ", Dr0 = " << Dr0 << '\n';
    simclock = 0.0;
  };

  // File openers
  void openEnergyObject(std::string& str) {
    enout.open(str.c_str());
    if (!enout.is_open()) {
      std::cout << "	ERROR: file could not open " << str << "..." << std::endl;
      exit(1);
    } else
      std::cout << "** Opening file " << str << " ..." << std::endl;
  }

  void openStressObject(std::string& str) {
    stressout.open(str.c_str());
    if (!stressout.is_open()) {
      std::cout << "	ERROR: file could not open " << str << "..." << std::endl;
      exit(1);
    } else
      std::cout << "** Opening file " << str << " ..." << std::endl;
  }

  // setters
  void setkbi(double val) { fill(kbi.begin(), kbi.end(), val); };
  void setRandPsi() {
    for (int ci = 0; ci < NCELLS; ci++)
      psi[ci] = (2.0 * drand48() - 1.0) * PI;
  }
  void seta0() {
  }

  // getters

  // main double
  double getPsi(int ci) { return psi.at(ci); };

  // editors & updates
  double meanl0();
  double meancalA0();
  double meankb();

  // epi cell interactions
  void vertexAttractiveForces2D_2();
  void attractiveForceUpdate_2();
  void activeAttractiveForceUpdate();

  // protocols
  void vertexCompress2Target2D(dpmMemFn forceCall, double Ftol, double dt0, double phi0Target, double dphi0);
  void tensileLoading(double scaleFactorX, double scaleFactorY);
  void zeroMomentum();
  void scaleBoxSize(double boxLengthScale, double scaleFactorX, double scaleFactorY);
  void dampedNVE2D(dpmMemFn forceCall, double B, double dt0, int NT, int NPRINTSKIP);

  int getIndexOfCellLocatedHere(double xLoc, double yLoc);
  void deleteCell(double sizeRatio, int nsmall, double xLoc, double yLoc);
  void laserAblate(int numCellsAblated, double sizeRatio, int nsmall, double xLoc, double yLoc);

  void notchTest(int numCellsToDelete, double boxLengthScale, double sizeRatio, int nsmall, dpmMemFn forceCall, double B, double dt0, int NT, int NPRINTSKIP, int maxit = 10, std::string loadingType = "uniaxial");
  void orientDirector(int ci, double xLoc, double yLoc);
  void deflectOverlappingDirectors();

  // polymorphism: write configuration information to file
  void printConfiguration2D();
};

#endif