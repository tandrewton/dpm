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

  // perimeter energy
  double kL;

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

  // whether cell ci has planted a flag, to which it will anchor one of its
  // vertices to the flag with a spring
  std::vector<bool> flag;

  // position of flags
  std::vector<std::vector<double>> flagPos;

  // motility scale factor for contact inhibition repolarization (active 'kick'
  // to help with repulsion)
  std::vector<double> activePropulsionFactor;

  // counter for number of deflected polarizations
  int polarizationCounter;

  // boolean for whether cells execute CIL or not
  bool boolCIL;

  // stores x length of box, used for tracking box size when appling tensile
  // strain set it after equilibrating before tensile strain, and never
  // otherwise.
  double initialLx;

  std::vector<double> VL;

  // stores nearest neighbors indices of each vertex according to adjacency
  // (intracell) and adhesion (intercell)
  std::vector<std::vector<int>> vnn;
  std::vector<int> vnn_label;
  std::vector<int> order;

  // specific output objects
  std::ofstream enout;
  std::ofstream stressout;
  std::ofstream bout;
  std::ofstream edgeout;
  std::ofstream purseout;

  // simulation time keeper (accumulates elapsed simulation time during MD
  // routines)
  double simclock;

  // simclock's last value before substrate adhesion springs update
  double previousUpdateSimclock;

  bool boolStopDeleting = false;

  // initial vertex radii before purseString shrinkage
  std::vector<double> initialRadius;
  std::vector<double> initiall0;
  double initialPreferredPerimeter;

  // flag for vertex repulsion (if a cell has only 1 wound vertex, then turn off repulsion so that it gets sucked into the bulk)
  std::vector<int> listTurnOffRepulsion;
  std::vector<int> sortedWoundIndices;
  bool woundIsClosedPolygon;
  bool regridChecker;
  std::vector<int> cellsLeavingPurseString;

  // virtual purseString variables
  std::vector<double> x_ps;
  std::vector<double> v_ps;
  std::vector<double> F_ps;
  std::vector<double> l0_ps;
  std::vector<int> psContacts;

 public:
  // constructor and destructor
  epi2D(int n, double att1, double att2, double Dr, int seed)
      : dpm(n, seed) {
    z.resize(n);
    // att = attraction;
    l1 = att1;
    l2 = att2;
    Dr0 = Dr;
    boolCIL = false;
    vector<double> temp(NCELLS, 0.0);
    vector<double> temp2(NCELLS, 1.0);
    // psi = temp;
    // activePropulsionFactor = temp2;
    psi.resize(NCELLS);
    activePropulsionFactor.resize(NCELLS);
    flag.resize(NCELLS);
    flagPos.resize(NCELLS);
    polarizationCounter = 0;
    for (int ci = 0; ci < NCELLS; ci++) {
      psi.at(ci) =
          PI / 2 * (ci % 2) -
          PI / 2 * ((ci + 1) % 2);  // should be : up if odd, down if even
      flagPos[ci].resize(NDIM);
      flag[ci] = false;
    }
    vector<double> temp3(NDIM * 2, 0.0);
    VL = temp3;
    cout << "Initializing epi2D object, l1 = " << l1 << ", l2 = " << l2
         << ", Dr0 = " << Dr0 << '\n';
    simclock = 0.0;
  };

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

  // File openers
  void openBoundaryObject(std::string& str) {
    bout.open(str.c_str());
    if (!bout.is_open()) {
      std::cerr << "	ERROR: file could not open " << str << "..."
                << std::endl;
      exit(1);
    } else
      std::cout << "** Opening file " << str << " ..." << std::endl;
  }

  // File openers
  void openEdgeObject(std::string& str) {
    edgeout.open(str.c_str());
    if (!edgeout.is_open()) {
      std::cerr << "	ERROR: file could not open " << str << "..."
                << std::endl;
      exit(1);
    } else
      std::cout << "** Opening file " << str << " ..." << std::endl;
  }

  // File openers
  void openPurseStringObject(std::string& str) {
    purseout.open(str.c_str());
    if (!purseout.is_open()) {
      std::cerr << "	ERROR: file could not open " << str << "..."
                << std::endl;
      exit(1);
    } else
      std::cout << "** Opening file " << str << " ..." << std::endl;
  }

  // setters
  void setkbi(double val) { fill(kbi.begin(), kbi.end(), val); };
  void setkL(double val) { kL = val; }
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
  double vertDistNoPBC(int gi, int gj);

  // initialization
  void monodisperse2D(double calA0, int n);

  // editors & updates
  double meanl0();
  double meancalA0();
  double meankb();
  double getPreferredPerimeter(int ci);
  double distanceLineAndPoint(double x1, double y1, double x2, double y2, double x0, double y0);
  void directorDiffusion();
  std::vector<int> regridSegment(int ci, double vrad);

  // epi cell interactions
  void repulsiveForceUpdateWithWalls();
  void vertexAttractiveForces2D_2();
  void attractiveForceUpdate_2();
  void epi_shapeForces2D();
  void activeAttractiveForceUpdate();
  void substrateadhesionAttractiveForceUpdate();
  void repulsiveForceWithCircularApertureWall();

  // protocols
  void vertexCompress2Target2D(dpmMemFn forceCall, double Ftol, double dt0, double phi0Target, double dphi0);
  void expandBoxAndCenterParticles(double boxLengthScaleFactor,
                                   double boxLengthScale);
  void ageCellAreas(double areaScaleFactor);
  void tensileLoading(double scaleFactorX, double scaleFactorY);
  void updateSubstrateSprings(double refreshInterval);
  void zeroMomentum();
  void scaleBoxSize(double boxLengthScale, double scaleFactorX, double scaleFactorY);
  void dampedNVE2D(dpmMemFn forceCall, double B, double dt0, double duration, double printInterval);
  void dampedNP0(dpmMemFn forceCall, double B, double dt0, double duration, double printInterval, bool wallsOn);
  void wallForces(bool top, bool bottom, bool left, bool right, double& forceTop, double& forceBottom, double& forceLeft, double& forceRight);
  void circularApertureForces(double radius);

  int getIndexOfCellLocatedHere(double xLoc, double yLoc);
  // note: whenever adding member-level data structures that depend on
  // NVTOT/NCELLS, need to make sure to modify the size in deleteCell
  // appropriately
  void deleteCell(double sizeRatio, int nsmall, double xLoc, double yLoc);
  void deleteVertex(std::vector<int>& deleteList);
  void laserAblate(int numCellsAblated, double sizeRatio, int nsmall, double xLoc, double yLoc);

  // void detection algorithms (Newman-Ziff)
  void initializevnn();
  void boundaries();
  std::vector<int> refineBoundaries();
  void NewmanZiff(std::vector<int>& ptr, int empty, int& mode, int& big);
  void printBoundaries(int nthLargestCluster = 1);
  void getWoundVertices(int nthLargestCluster = 1);
  bool checkWoundClosedPolygon(std::vector<int>& listOfIndices);
  int findRoot(int i, std::vector<int>& ptr);

  void notchTest(int numCellsToDelete, double strain, double strainRate, double boxLengthScale, double sizeRatio, int nsmall, dpmMemFn forceCall, double B, double dt0, double printInterval, std::string loadingType);
  void orientDirector(int ci, double xLoc, double yLoc);
  void deflectOverlappingDirectors();
  void evaluateGhostDPForces(std::vector<int>& giList, double trate);
  double rotateAndCalculateArcLength(int ci, std::vector<int>& woundIndicesBelongingToCi);
  void purseStringContraction(double trate);
  void purseStringContraction2(double trate, double deltaSq, double k_wp, double B);
  void initializePurseStringVariables();
  void updatePurseStringContacts();
  void evaluatePurseStringForces(double deltasq, double k_wp, double B);
  void integratePurseString(double deltaSq, double k_wp, double B);
  // polymorphism: write configuration information to file
  void printConfiguration2D(double deltaSq = 1);
};

#endif