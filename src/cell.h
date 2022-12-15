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
  std::ofstream tissueout;  // output stream for tissue measurements over time

  // simulation-scale stored quantities
  double simclock;         // accumulator for simulation time
  std::vector<double> VL;  // wall velocity if simulating with wall forces
  std::vector<double> XL;  // position of box boundaries (there are 4 box boundaries) : xrange (XL[0] XL[2]), yrange (XL[1], XL[3])
  // XL modifies the box size by mapping 0 -> 0 + XL[0], L[0] -> L[0] + XL[2], 0 -> 0 + XL[1], L[1] -> L[1] + XL[3]

  // energy/force modifiers
  std::vector<int> cellID;                          // 0/1/2/... = NT/PSM1/NC/PSM2/cell0/cell1/cell2/cell3
  std::vector<std::vector<double>> cellTypeIntMat;  // size(cellID) x size(cellID) matrix where (i,j)th entry is the force modifier for interactions between cell types i and j

  // activity data
  std::vector<bool> cellTouchesWallsLeft;
  std::vector<bool> cellTouchesWallsRight;

  // active brownian variables and parameters
  std::vector<double> psi;  // directors
  double v0_ABP, tau_ABP;   // velocity parameter, rotation time constant

 public:
  cell(int n, double att1, double att2, int seed, int numCellTypes)
      : dpm(n, seed) {
    // set values here
    vector<double> zerosNCELLS(NCELLS, 0.0);
    vector<int> zerosNCELLS_int(NCELLS, 0);
    vector<int> zerosNumCellTypes(numCellTypes, 0);
    vector<bool> falseNCELLS(NCELLS, false);
    vector<double> temp3(NDIM * 2, 0.0);
    l1 = att1;
    l2 = att2;
    simclock = 0.0;
    VL = temp3;
    XL = VL;
    cellID = zerosNCELLS_int;  // cellID is a vector of size NCELLS that determines how to evaluate cell-cell forces using cellTypeIntMat
    // make sure to add and subtract from cellID when adding or removing cells
    vector<vector<double>> onesCellIDxCellID(numCellTypes, vector<double>(vector<double>(numCellTypes, 1.0)));
    cellTypeIntMat = onesCellIDxCellID;  // all interactions start with multiplicative modifier 1.0
    cellTouchesWallsLeft = falseNCELLS;
    cellTouchesWallsRight = falseNCELLS;
    v0_ABP = 0.0;
    tau_ABP = 1.0;
    psi = zerosNCELLS;
    cout << "Initializing epi2D object, l1 = " << l1 << ", l2 = " << l2 << '\n';
  }

  /*// destructor
  ~cell() {
    // clear all private vectors
    // should update this soon
    cout << "destructor\n";
    L.clear();
    pbc.clear();
    a0.clear();
    l0.clear();
    t0.clear();
    nv.clear();
    szList.clear();
    im1.clear();
    ip1.clear();
    r.clear();
    x.clear();
    v.clear();
    F.clear();
    stress.clear();
    sb.clear();
    lb.clear();
    cout << "nn size = " << nn.size() << '\n';
    cout << "NBX = " << NBX << '\n';
    for (int i = 0; i < NBX; i++)
      nn.at(i).clear();
    nn.clear();
    head.clear();
    last.clear();
    list.clear();

    if (posout.is_open())
      posout.close();
    if (enout.is_open())
      enout.close();
    if (stressout.is_open())
      stressout.close();

    cellID.clear();
    cellTypeIntMat.clear();
    cellTouchesWallsLeft.clear();
    cellTouchesWallsRight.clear();
  }*/

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

  // cell interactions
  void shapeForces2D();
  void repulsiveForceUpdateWithWalls();
  // void repulsiveForceUpdateWithoutTopWall();
  void attractiveForceUpdate();
  void attractiveForceUpdateWithCrawling();
  void attractiveSmoothForceUpdate();
  void repulsiveWithPolarityForceUpdate();
  void attractiveWithPolarityForceUpdate();
  void attractiveWithPolarityForceAndWallCrawlingUpdate();
  void attractiveWallCrawlingForceUpdate();
  void repulsiveForceUpdateWithPolyWall();
  void attractiveForceUpdateWithPolyWall();
  void attractiveForceUpdatePrint(double& forceX, double& forceY, double& energy);
  void attractiveForceUpdateSmoothPrint(double& forceX, double& forceY, double& energy);
  void vertexAttractiveForces2D_2();
  void vertexAttractiveForces2D_test(double& energy);
  void circuloLineAttractiveForces();
  void calculateSmoothInteraction(double& rx, double& ry, double& sij, double& shellij, double& cutij, double& kint, double& kc, int& gi, int& gj, double& projection, int& ci, int& cj);
  void wallForces(bool left, bool bottom, bool right, bool top, double& forceLeft, double& forceBottom, double& forceRight, double& forceTop, double appliedUniaxialPressure = 0.0);
  void wallCrawlingForces();
  void cellPolarityForces(int ci, double k_polarity, std::string direction = "y");
  void setActiveBrownianParameters(double v0, double tau) {
    v0_ABP = v0;
    tau_ABP = tau;
  }
  void brownianCrawlingUpdate();
  void directorDiffusion();

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

  void openTissueObject(std::string& str) {
    tissueout.open(str.c_str());
    if (!tissueout.is_open()) {
      std::cout << "	ERROR: file could not open " << str << "..."
                << std::endl;
      exit(1);
    } else
      std::cout << "** Opening file " << str << " ..." << std::endl;
  }

  // boundary routines
  void replacePolyWallWithDP(int numCellTypes);
  void addDP(int numVerts, vector<double>& dp_x, vector<double>& dp_y, int cellTypeIndex, int numCellTypes);

  // routines
  void initializeFourTransverseTissues(double phi0, double Ftol);
  void initializeTransverseTissue(double phi0, double Ftol);
  void vertexCompress2Target2D(dpmMemFn forceCall, double Ftol, double dt0, double phi0Target, double dphi0);
  void vertexCompress2Target2D_polygon(dpmMemFn forceCall, double Ftol, double dt0, double phi0Target, double dphi0);
  void simulateDampedWithWalls(dpmMemFn forceCall, double dt0, double duration, double printInterval, double pressureRate, double adhesionRate, bool wallsOn, bool leftOpen, bool bottomOpen, bool rightOpen, bool topOpen, double trueStrainRateX = 0.0, double trueStrainRateY = 0.0, double appliedUniaxialPressure = 0.0);
  void vertexNVE(dpmMemFn forceCall, double T, double dt0, double duration, double printInterval);
  void dampedVertexNVE(dpmMemFn forceCall, double dt0, double duration, double printInterval);
  void takeTissueMeasurements(int cellBoundaryType);

  // printouts
  void printConfiguration2D();
};

#endif