#ifndef EPI2D_H
#define EPI2D_H

/*

	HEADER FILE FOR epi CLASS

		-- Inherits from DPM class 
		-- For collections of epithelial cells
		-- Incorporates ?
		-- ONLY FOR 2D

*/

#include <math.h>
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

  // active vertex velocity scale
  double v0;

  // directors (direction vectors) for activity
  std::vector<double> psi;

  // rotational diffusion constant
  double Dr0;

  // whether cell ci has planted a flag, to which it will anchor one of its vertices to the flag with a spring
  std::vector<bool> flag;

  // position of flags
  std::vector<std::vector<double>> flagPos;

  // motility scale factor for contact inhibition repolarization (active 'kick' to help with repulsion)
  std::vector<double> activePropulsionFactor;

  // counter for number of deflected polarizations
  int polarizationCounter;

  // boolean for whether cells execute CIL or not
  bool boolCIL;

  // stores x length of box, used for tracking box size when appling tensile strain
  // set it after equilibrating before tensile strain, and never otherwise.
  double initialLx;

  std::vector<double> VL;

  // stores nearest neighbors indices of each vertex according to adjacency (intracell) and adhesion (intercell)
  std::vector<std::vector<int>> vnn;
  std::vector<int> vnn_label;

  // specific output objects
  std::ofstream enout;
  std::ofstream stressout;
  std::ofstream bout;

  // simulation time keeper (accumulates elapsed simulation time during MD routines)
  double simclock;

  // simclock's last value before substrate adhesion springs update
  double previousUpdateSimclock;

 public:
  // constructor and destructor
  epi2D(int n, double att1, double att2, double Dr, int seed)
      : dpm(n, seed) {
    z.resize(n);
    //att = attraction;
    l1 = att1;
    l2 = att2;
    Dr0 = Dr;
    boolCIL = false;
    vector<double> temp(NCELLS, 0.0);
    vector<double> temp2(NCELLS, 1.0);
    //psi = temp;
    //activePropulsionFactor = temp2;
    psi.resize(NCELLS);
    activePropulsionFactor.resize(NCELLS);
    flag.resize(NCELLS);
    flagPos.resize(NCELLS);
    polarizationCounter = 0;
    for (int ci = 0; ci < NCELLS; ci++) {
      psi.at(ci) = PI / 2 * (ci % 2) - PI / 2 * ((ci + 1) % 2);  //should be : up if odd, down if even
      flagPos[ci].resize(NDIM);
      flag[ci] = false;
    }
    vector<double> temp3(NDIM * 2, 0.0);
    VL = temp3;
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

  // File openers
  void openBoundaryObject(std::string& str) {
    bout.open(str.c_str());
    if (!bout.is_open()) {
      std::cerr << "	ERROR: file could not open " << str << "..." << std::endl;
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
  void setv0(double val) { v0 = val; };
  void scalea0(double scaleFactor) {
    for (int ci = 0; ci < NCELLS; ci++)
      a0[ci] *= scaleFactor;
  }
  void setboolCIL(bool val) { boolCIL = val; };
  void setInitialLx(double val) { initialLx = val; };

  // getters
  double getInitialLx() { return initialLx; };

  // main double
  double getPsi(int ci) { return psi.at(ci); };

  // initialization
  void monodisperse2D(double calA0, int n);

  // editors & updates
  double meanl0();
  double meancalA0();
  double meankb();
  double distanceLineAndPoint(double x1, double y1, double x2, double y2, double x0, double y0);
  void directorDiffusion();

  // epi cell interactions
  void repulsiveForceUpdateWithWalls();
  void vertexAttractiveForces2D_2();
  void attractiveForceUpdate_2();
  void activeAttractiveForceUpdate();
  void substrateadhesionAttractiveForceUpdate();

  // protocols
  void vertexCompress2Target2D(dpmMemFn forceCall, double Ftol, double dt0, double phi0Target, double dphi0);
  void expandBoxAndCenterParticles(double boxLengthScaleFactor, double boxLengthScale);
  void ageCellAreas(double areaScaleFactor);
  void tensileLoading(double scaleFactorX, double scaleFactorY);
  void updateSubstrateSprings(double refreshInterval);
  void zeroMomentum();
  void scaleBoxSize(double boxLengthScale, double scaleFactorX, double scaleFactorY);
  void dampedNVE2D(dpmMemFn forceCall, double B, double dt0, double duration, double printInterval);
  void dampedNP0(dpmMemFn forceCall, double B, double dt0, double duration, double printInterval, bool wallsOn);
  void wallForces(bool top, bool bottom, bool left, bool right, double& forceTop, double& forceBottom, double& forceLeft, double& forceRight);

  int getIndexOfCellLocatedHere(double xLoc, double yLoc);
  // note: whenever adding member-level data structures that depend on NVTOT/NCELLS, need to make sure to modify the size in deleteCell appropriately
  void deleteCell(double sizeRatio, int nsmall, double xLoc, double yLoc);
  void laserAblate(int numCellsAblated, double sizeRatio, int nsmall, double xLoc, double yLoc);

  // void detection algorithms (Newman-Ziff)
  void initializevnn();
  void boundaries();
  std::vector<int> refineBoundaries();
  void printBoundaries();
  int findRoot(int i, std::vector<int>& ptr);

  void notchTest(int numCellsToDelete, double strain, double strainRate, double boxLengthScale, double sizeRatio, int nsmall, dpmMemFn forceCall, double B, double dt0, double printInterval, std::string loadingType);
  void orientDirector(int ci, double xLoc, double yLoc);
  void deflectOverlappingDirectors();

  // polymorphism: write configuration information to file
  void printConfiguration2D();
};

#endif