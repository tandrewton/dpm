#ifndef CELL_H
#define CELL_H

#include <algorithm>
#include <random>
#include <stdexcept>
#include "dpm.h"

// namespace
using namespace std;

class cell;
typedef void (cell::*cellMemFn)(void);

class cell : public dpm {
 protected:
  // specific output objects
  std::ofstream enout;         // output stream for energy data
  std::ofstream stressout;     // output stream for stress data
  std::ofstream tissueout;     // output stream for tissue measurements over time
  std::ofstream catchBondOut;  // output stream for catch bond positions

  double simclock;             // accumulator for simulation time
  std::vector<double> VL, XL;  // wall velocity and positions if simulating with wall forces.
  // position of box boundaries (there are 4 box boundaries) : xrange (XL[0] XL[2]), yrange (XL[1], XL[3])
  // XL modifies the box size by mapping 0 -> 0 + XL[0], L[0] -> L[0] + XL[2], 0 -> 0 + XL[1], L[1] -> L[1] + XL[3]

  // energy/force modifiers
  std::vector<int> cellID;                          // 0/1/2/... = NT/PSM1/NC/PSM2/cell0/cell1/cell2/cell3
  std::vector<std::vector<double>> cellTypeIntMat;  // size(cellID) x size(cellID) matrix where (i,j)th entry is the force modifier for interactions between cell types i and j
  double originalVertexRadius;

  // activity data
  std::vector<bool> cellTouchesWallsLeft;
  std::vector<bool> cellTouchesWallsRight;

  // active brownian variables and parameters
  std::vector<double> psi;  // directors
  double v0_ABP, tau_ABP;   // velocity parameter, rotation time constant

  // catch bonds for substrate adhesion parameters
  std::vector<bool> isGiCatchBonded;
  std::vector<std::vector<double>> catchBondPosition;
  double k_off, k_ecm;

 public:
  cell(int n, double att1, double att2, int seed, int numCellTypes)
      : dpm(n, seed) {
    // set values here
    l1 = att1;
    l2 = att2;
    originalVertexRadius = 0.0;
    simclock = 0.0;
    VL = vector<double>(NDIM * 2, 0.0);
    XL = VL;
    cellID = vector<int>(n, 0);  // cellID determines how to evaluate cell-cell forces using cellTypeIntMat
    // make sure to add and subtract from cellID when adding or removing cells
    cellTypeIntMat = vector<vector<double>>(numCellTypes, vector<double>(numCellTypes, 1.0));  // all interactions start with multiplicative modifier 1.0
    cellTouchesWallsLeft = vector<bool>(n, false);
    cellTouchesWallsRight = vector<bool>(n, false);
    v0_ABP = 0.0;
    tau_ABP = 1.0;
    k_off = 1.0;
    k_ecm = 1.0;
    psi = vector<double>(n, 0.0);
    cout << "Initializing epi2D object, l1 = " << l1 << ", l2 = " << l2 << '\n';
  }

  // test routines for force calculation
  void moveVertex(int gi, double xpos, double ypos);

  // cell-cell interaction matrix routines
  void setCellTypeAttractionModifiers(int i, int j, double val) {
    cellTypeIntMat[i][j] = val;
    cellTypeIntMat[j][i] = val;
  }
  void removeCellIDFromInteractionMatrix(int cellID);
  void addCellIDToInteractionMatrix(int cellID);
  void printInteractionMatrix();
  // consider writing a function that uses cell IDs

  void monodisperseSmooth(double calA0, int n);
  void initializeVertexShapeParametersSmooth(double calA0, int nref);
  // cell interactions
  void maxwellRelaxationRestLengths(std::vector<double>& l, std::vector<int> cellTypes);
  void shapeForces2D();
  void repulsiveForceUpdateWithWalls();
  // void repulsiveForceUpdateWithoutTopWall();
  void attractiveForceUpdate();
  void attractiveForceUpdateWithCrawling();
  void attractiveSmoothForceUpdateWithCrawling();
  void attractiveSmoothActiveCatchBonds();
  void attractiveSmoothForceUpdate();
  void attractiveSmoothForceUpdateWithPolyWall();
  void repulsiveWithPolarityForceUpdate();
  void attractiveWithPolarityForceUpdate();
  void attractiveWithPolarityForceAndWallCrawlingUpdate();
  void attractiveWallCrawlingForceUpdate();
  void repulsiveForceUpdateWithPolyWall();
  void attractiveForceUpdateWithPolyWall();
  void vertexAttractiveForces2D_2();
  void vertexAttractiveForces2D_test(double& energy);
  void circuloLineAttractiveForces();
  void calculateSmoothInteraction(double rx, double ry, double sij, double shellij, double cutij, double kint, double kc, int gi, int gj, double& projection, int& ci, int& cj);
  void wallForces(bool left, bool bottom, bool right, bool top, double& forceLeft, double& forceBottom, double& forceRight, double& forceTop, double appliedUniaxialPressure = 0.0);
  void wallCrawlingForces();
  void cellPolarityForces(int ci, double k_polarity, std::string direction = "y");
  void setActiveBrownianParameters(double v0, double tau) {
    v0_ABP = v0;
    tau_ABP = tau;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 2 * PI);
    std::generate(psi.begin(), psi.end(), [&]() { return dis(gen); });
  }
  void setkecm(double val) { k_ecm = val; }
  void setkoff(double val) { k_off = val; }
  void brownianCrawlingUpdate();
  void directorDiffusion();
  void catchBondsUpdate();
  void resizeCatchBonds();

  // File openers
  void openEnergyObject(const std::string& filename) { openFile(enout, filename); }
  void openStressObject(const std::string& filename) { openFile(stressout, filename); }
  void openTissueObject(const std::string& filename) { openFile(tissueout, filename); }
  void openCatchBondObject(const std::string& filename) { openFile(catchBondOut, filename); }

  // boundary routines
  void replacePolyWallWithDP(int numCellTypes);
  void addDP(int numVerts, const vector<double>& dp_x, const vector<double>& dp_y, int cellTypeIndex, int numCellTypes);

  // routines
  void initializeFourTransverseTissues(double phi0, double Ftol);
  void initializeTransverseTissue(double phi0, double Ftol, int polyShapeID = 0);
  void vertexCompress2Target2D(dpmMemFn forceCall, double Ftol, double dt0, double phi0Target, double dphi0, bool isFIRE = true);
  void vertexCompress2Target2D_polygon(dpmMemFn forceCall, double Ftol, double dt0, double phi0Target, double dphi0, bool isFIRE = true);
  void shrinkCellVertices(dpmMemFn forceCall, double dt0, double shrinkRatio);
  void simulateDampedWithWalls(dpmMemFn forceCall, double dt0, double duration, double printInterval, double pressureRate, double adhesionRate, bool wallsOn, bool leftOpen, bool bottomOpen, bool rightOpen, bool topOpen, double trueStrainRateX = 0.0, double trueStrainRateY = 0.0, double appliedUniaxialPressure = 0.0);
  void vertexNVE(dpmMemFn forceCall, double T, double dt0, double duration, double printInterval);
  void dampedVertexNVE(dpmMemFn forceCall, double dt0, double duration, double printInterval);
  void takeTissueMeasurements(int cellBoundaryType);

  // printouts
  void printConfiguration2D();
};

#endif