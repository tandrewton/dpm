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

  std::vector<double> VL; // velocity of box lengths L (there are 4 box lengths)

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
  std::ofstream vout;
  std::ofstream bulkout;
  std::ofstream innerout;
  std::ofstream woundPropertiesout;
  std::ofstream cellIDout;

  // simulation time keeper (accumulates elapsed simulation time during MD
  // routines)
  double simclock;

  // simclock's last value before substrate adhesion springs update
  double previousUpdateSimclock;

  // initial vertex info before purseString shrinkage
  std::vector<double> initialRadius;
  std::vector<double> initiall0;
  double initialPreferredPerimeter;
  std::vector<int> initialWoundCellIndices;
  double woundCenterX;
  double woundCenterY;
  double woundArea = 1e10;

  // flag for vertex repulsion (if a cell has only 1 wound vertex, then turn off repulsion so that it gets sucked into the bulk)
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
  std::vector<bool> isSpringBroken;
  bool isPurseStringDoneShrinking;

  double strainRate_ps, k_ps, k_LP, tau_LP, deltaSq, maxProtrusionLength;
  // purse string (ps) strain rate (constant true strain rate)
  // ps spring constant between wound vertex and ps vertex
  // lamellipodia (lp) spring constant between protruded spring anchor (flag) and a colinear vertex
  // lp timescale (10% chance of decaying per tau_LP)
  // squared max distance for wound-ps vertex spring (units of vdiam)
  // max length of lp (units of vdiam)
  std::vector<double> timeElapsedSinceFlagPlanted;
  std::vector<double> restLengthLPx;  // rest length is x_vertex - x_flag, can be negative
  std::vector<double> restLengthLPy;
  std::vector<int> giConnectedToFlag;  // will have to modify this when regridding, haven't done that yet.

  double shapeRelaxationRate;
  double U_crawling;
  double U_ps;

 public:
  // constructor and destructor
  epi2D(int n, double att1, double att2, double Dr, double omega_ps, double kps, double kLP, double tauLP, double deltaSquared, double maxCrawlLength, int seed)
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
    timeElapsedSinceFlagPlanted.resize(NCELLS);
    restLengthLPx.resize(NCELLS);
    restLengthLPy.resize(NCELLS);
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
    strainRate_ps = omega_ps;
    k_ps = kps;
    k_LP = kLP;
    tau_LP = tauLP;
    deltaSq = deltaSquared;
    maxProtrusionLength = maxCrawlLength;
    isPurseStringDoneShrinking = false;
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

  void openBoundaryObject(std::string& str) {
    bout.open(str.c_str());
    if (!bout.is_open()) {
      std::cerr << "	ERROR: file could not open " << str << "..."
                << std::endl;
      exit(1);
    } else
      std::cout << "** Opening file " << str << " ..." << std::endl;
  }

  void openVoidAreaObject(std::string& str) {
    vout.open(str.c_str());
    if (!vout.is_open()) {
      std::cerr << "	ERROR: file could not open " << str << "..."
                << std::endl;
      exit(1);
    } else
      std::cout << "** Opening file " << str << " ..." << std::endl;
  }

  void openEdgeObject(std::string& str) {
    edgeout.open(str.c_str());
    if (!edgeout.is_open()) {
      std::cerr << "	ERROR: file could not open " << str << "..."
                << std::endl;
      exit(1);
    } else
      std::cout << "** Opening file " << str << " ..." << std::endl;
  }

  void openPurseStringObject(std::string& str) {
    purseout.open(str.c_str());
    if (!purseout.is_open()) {
      std::cerr << "	ERROR: file could not open " << str << "..."
                << std::endl;
      exit(1);
    } else
      std::cout << "** Opening file " << str << " ..." << std::endl;
  }

  void openBulkCellShapeObject(std::string& str) {
    bulkout.open(str.c_str());
    if (!bulkout.is_open()) {
      std::cerr << "	ERROR: file could not open " << str << "..."
                << std::endl;
      exit(1);
    } else
      std::cout << "** Opening file " << str << " ..." << std::endl;
  }

  void openInnerCellShapeObject(std::string& str) {
    innerout.open(str.c_str());
    if (!innerout.is_open()) {
      std::cerr << "	ERROR: file could not open " << str << "..."
                << std::endl;
      exit(1);
    } else
      std::cout << "** Opening file " << str << " ..." << std::endl;
  }

  void openWoundPropertiesObject(std::string& str) {
    woundPropertiesout.open(str.c_str());
    if (!woundPropertiesout.is_open()) {
      std::cerr << "	ERROR: file could not open " << str << "..."
                << std::endl;
      exit(1);
    } else
      std::cout << "** Opening file " << str << " ..." << std::endl;
  }

  void openCellIDObject(std::string& str) {
    cellIDout.open(str.c_str());
    if (!woundPropertiesout.is_open()) {
      std::cerr << "	ERROR: file could not open " << str << "..."
                << std::endl;
      exit(1);
    } else
      std::cout << "** Opening file " << str << " ..." << std::endl;
  }


  // setters
  void setkbi(double val) { fill(kbi.begin(), kbi.end(), val); };
  void setShapeRelaxationRate(double val) { shapeRelaxationRate = val; }
  void setRandPsi() {
    for (int ci = 0; ci < NCELLS; ci++)
      psi[ci] = (2.0 * drand48() - 1.0) * PI;
  }
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

  // editors & updates
  double meanl0();
  double meancalA0();
  double meankb();
  double getPreferredPerimeter(int ci);
  void directorDiffusion();
  std::vector<int> regridSegment(int ci, double vrad);
  void resetActiveEnergy();

  // epi cell interactions
  void repulsiveForceUpdateWithWalls();
  void vertexAttractiveForces2D_2();
  void circuloLineAttractiveForces();
  void attractiveForceUpdate_circulo();
  void attractiveForceUpdate_2();
  void activeAttractiveForceUpdate();
  void substrateadhesionAttractiveForceUpdate();
  void repulsiveForceWithCircularApertureWall();
  void repulsiveForceUpdateWithPolyWall();
  void attractiveForceUpdateWithPolyWall();
  void crawlingWithPurseStringAndCircularWalls();

  // protocols
  void expandBoxAndCenterParticles(double boxLengthScaleFactor,
                                   double boxLengthScale);
  void ageCellAreas(double areaScaleFactor);
  void ageCellPerimeters(double shapeRelaxationRate, double dt);
  void tensileLoading(double scaleFactorX, double scaleFactorY);
  void updateSubstrateSprings();
  void zeroMomentum();
  void scaleBoxSize(double boxLengthScale, double scaleFactorX, double scaleFactorY);
  void dampedNVETest(dpmMemFn forceCall, double T, double dt0, int NT, int NPRINTSKIP);
  void dampedNVE2D(dpmMemFn forceCall, double B, double dt0, double duration, double printInterval);
  void dampedNP0(dpmMemFn forceCall, double B, double dt0, double duration, double printInterval, int purseStringOn = 0);
  void wallForces(bool left, bool bottom, bool right, bool top, double& forceLeft, double& forceBottom, double& forceRight, double& forceTop, int forceOption = 0);
  void circularApertureForces(double radius);
  void vertexCompress2Target2D_polygon(dpmMemFn forceCall, double Ftol, double dt0, double phi0Target, double dphi0);

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
  void NewmanZiff(std::vector<int>& ptr, int empty, int& mode, int& big, std::vector<int>& order);
  void printBoundaries(int nthLargestCluster = 1);
  void getWoundVertices(int nthLargestCluster = 1);
  bool checkWoundClosedPolygon(std::vector<int>& listOfIndices);
  double computeWoundVerticesUsingRays(double& woundCenterX, double& woundCenterY, int numRays);
  int findRoot(int i, std::vector<int>& ptr);
  double calculateWoundArea(double& woundPointX, double& woundPointY);
  bool isPointInPolygons(double xloc, double yloc);
  int pnpoly(int nvert, std::vector<double> vertx, std::vector<double> verty, double testx, double testy);
  double calculateArea(std::vector<double>& vertx, std::vector<double>& verty);
  double calculateAreaFlattened(std::vector<double>& vertPosFlattened);

  void notchTest(int numCellsToDelete, double strain, double strainRate, double boxLengthScale, double sizeRatio, int nsmall, dpmMemFn forceCall, double B, double dt0, double printInterval, std::string loadingType);
  void orientDirector(int ci, double xLoc, double yLoc);
  void deflectOverlappingDirectors();
  double getDistanceToVertexAtAnglePsi(int ci, double psi_ci, double cx, double cy, int& gi);
  double getDistanceToRandomUnadheredVertex(int ci, double cx, double cy, int&gi);
  double rotateAndCalculateArcLength(int ci, std::vector<int>& woundIndicesBelongingToCi);
  void purseStringContraction(double B);
  void initializePurseStringVariables();
  void updatePurseStringContacts();
  void evaluatePurseStringForces(double B);
  void integratePurseString(double B);
  // polymorphism: write configuration information to file
  void printConfiguration2D();
};

#endif