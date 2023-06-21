// simulate transverse section, focusing on tissue structure
// 2D justification: cells don't move much in the Anterior-Posterior plane of the PSM
//   the measurements are taken in the 2D transverse section, so a 2D simulation would help
//   our understanding of tissue structure in those same 2D sections.
// Scientific question: Why do cells gain extracellular space in cadherin AND/OR fibronectin/fibrillin mutants?
//  cadherin mutants lose intercellular adhesion. fN/fbn mutants lose cell-ECM adhesion.
// My understanding - cell-ECM adhesion acts as a compression. The ECM proteins link up with the actin skeleton
//   of cells on the boundary of the PSM. This causes action like a purse-string at the boundary of the tissue
//   which might explain why PSM is rounded up like a cylinder.

// Compilation command:
// g++ -O3 --std=c++11 -g -I src main/cell/psm2D.cpp src/dpm.cpp src/cell.cpp -o main/cell/psm2D.o
// run command:

/*
./main/cell/psm2D.o   12   25 1.05 0.85 0.01  25.0   0.01  1.0   1.0   1.0   1     400    test1
./main/cell/psm2D.o   12   25 1.05 0.85 0.05  25.0   0.01  1.0   1.0   1.0   1    400    test2
./main/cell/psm2D.o   12   25 1.05 0.85 0.1   25.0   0.01  1.0   1.0   1.0   1    400    test3
./main/cell/psm2D.o   12   25 1.05 0.85 0.2   25.0   0.01  1.0   1.0   1.0   1    400    test4

./main/cell/psm2D.o   12   16 1.05 0.75 0.0   0.0   0.05   50.0  1.0     1.0   1    50    test5
./main/cell/psm2D.o   12   16 1.05 0.75 0.01  0.0   0.05   50.0  1.0     1.0   1    50    test6
./main/cell/psm2D.o   12   16 1.05 0.75 0.1   0.0   0.05   50.0  1.0     1.0   1    50    test7
./main/cell/psm2D.o   8   16 1.05 0.75 0.2    0.0   0.05   50.0  1.0    0.1   1    50    test8

./main/cell/psm2D.o   40   16 1.05 0.9 0.1    0.0   0.1   50.0   1.0     1.0   1    100    test9
./main/cell/psm2D.o   40   16 1.05 0.9 0.1    0.0   0.05   50.0  1.0     1.0   1    1000    test10
./main/cell/psm2D.o   40   16 1.05 0.9 0.1    0.0   0.01   50.0  1.0     1.0   1    1000    test11

./main/cell/psm2D.o   40   16 1.05 0.9 0.05   0.0   0.1   50.0   1.0     1.0   1    1000    test12
./main/cell/psm2D.o   40   16 1.05 0.9 0.05   0.0   0.05   50.0  1.0     1.0   1    1000    test13
./main/cell/psm2D.o   40   16 1.05 0.9 0.05   0.0   0.01   50.0  1.0     1.0   1    1000    test14

./main/cell/psm2D.o   40   16 1.05 0.9 0.01   0.0   0.1   50.0   1.0     1.0   1    1000    test15
./main/cell/psm2D.o   40   16 1.05 0.9 0.01   0.0   0.05   50.0  1.0     1.0   1    1000    test16
./main/cell/psm2D.o   40   16 1.05 0.9 0.01   0.0   0.01   50.0  1.0     1.0   1    1000    test17

./main/cell/psm2D.o   40   16 1.05 0.9   0    0.0   0.1   50.0    1.0     1.0   1    1000    test18
./main/cell/psm2D.o   40   16 1.05 0.9   0    0.0   0.05   50.0   1.0     1.0   1    1000    test19
./main/cell/psm2D.o   40   16 1.05 0.9   0    0.0   0.01   50.0   1.0     1.0   1    1000    test20

*/
//                  NCELLS NV  A0  phi att t_maxwell v0  tau_abp k_ecm k_off seed duration outFileStem

#include <sstream>
#include "cell.h"

#define NDIM 2
using namespace std;

// helper functions to parse arguments
template <typename T>
T parseArg(const std::string& arg) {
  T value;
  std::stringstream(arg) >> value;
  return value;
}

// global constants
const bool plotCompression = 0;     // whether or not to plot configuration during compression protocol (0 saves memory)
const double dphi0 = 0.01;          // packing fraction increment
const double ka = 1.0;              // area force spring constant (should be unit)
const double kc = 1.0;              // interaction force spring constant (should be unit)
const double kb = 0.01;             // bending energy spring constant (should be zero)
const double kl = 1.0;              // segment length interaction force (should be unit)
const double boxLengthScale = 2.5;  // neighbor list box size in units of initial l0
const double phi0 = 0.91;           // initial preferred packing fraction
const double dt0 = 0.01;            // initial magnitude of time step in units of MD time
const double Ptol = 1e-5;
const double Ftol = 1e-6;
const double att_range = 0.3;

int main(int argc, char const* argv[]) {
  // local variables to be read in
  double B = 1.0;
  // Read command-line arguments into corresponding variables
  int NCELLS = parseArg<int>(argv[1]);
  int nv = parseArg<int>(argv[2]);
  double calA0 = parseArg<double>(argv[3]);
  double phi = parseArg<double>(argv[4]);
  double att = parseArg<double>(argv[5]);
  double t_stress = parseArg<double>(argv[6]);
  double v0_abp = parseArg<double>(argv[7]);
  double tau_abp = parseArg<double>(argv[8]);
  double k_ecm = parseArg<double>(argv[9]);
  double k_off = parseArg<double>(argv[10]);
  int seed = parseArg<int>(argv[11]);
  double runTime = parseArg<double>(argv[12]);
  std::string outFileStem = argv[13];

  string positionFile = outFileStem + ".pos";
  string tissueFile = outFileStem + ".tissue";
  string catchBondFile = outFileStem + ".catchBond";

  int numCellTypes = 2;  // 0 = interior cell type (PSM) and 1 = exterior cell type (boundary)
  cell cell2D(NCELLS, 0.0, 0.0, seed, numCellTypes);
  cell2D.openPosObject(positionFile);
  cell2D.openTissueObject(tissueFile);
  cell2D.openCatchBondObject(catchBondFile);

  // set spring constants
  cell2D.setka(ka);
  cell2D.setkl(kl);
  cell2D.setkb(kb);
  cell2D.setkc(kc);
  cell2D.setkecm(k_ecm);
  cell2D.setkoff(k_off);
  cell2D.setB(B);
  if (t_stress > 0.0)
    cell2D.setMaxwellRelaxationTime(t_stress);  // t_stress is infinity unless this is uncommented
  cell2D.setpbc(0, false);                      //  specify non-periodic boundaries
  cell2D.setpbc(1, false);
  cell2D.setl1(att);  // set adhesion scales
  cell2D.setl2(att_range);
  assert(att < att_range);  // required to have a differentiable, finite adhesive potential

  dpmMemFn repulsiveForceUpdate = &dpm::repulsiveForceUpdate;
  dpmMemFn attractiveForceUpdate = static_cast<void (dpm::*)()>(&cell::attractiveForceUpdate);
  dpmMemFn attractionWithActiveBrownianUpdate = static_cast<void (dpm::*)()>(&cell::attractiveForceUpdateWithCrawling);
  dpmMemFn attractionSmoothWithActiveBrownianUpdate = static_cast<void (dpm::*)()>(&cell::attractiveSmoothForceUpdateWithCrawling);
  dpmMemFn attractionSmoothActiveBrownianCatchBondsUpdate = static_cast<void (dpm::*)()>(&cell::attractiveSmoothActiveCatchBonds);
  dpmMemFn attractiveSmoothForceUpdate = static_cast<void (dpm::*)()>(&cell::attractiveSmoothForceUpdate);
  dpmMemFn attractiveSmoothWithPolyWalls = static_cast<void (dpm::*)()>(&cell::attractiveSmoothForceUpdateWithPolyWall);
  dpmMemFn attractivePolarityForceUpdate = static_cast<void (dpm::*)()>(&cell::attractiveWithPolarityForceUpdate);
  dpmMemFn repulsiveForceUpdateWithPolyWalls = static_cast<void (dpm::*)()>(&cell::repulsiveForceUpdateWithPolyWall);
  dpmMemFn attractiveForceUpdateWithPolyWalls = static_cast<void (dpm::*)()>(&cell::attractiveForceUpdateWithPolyWall);

  // cellTypeIntMat = ones to begin with
  for (int cellIDi = 0; cellIDi < numCellTypes; cellIDi++) {
    for (int cellIDj = 0; cellIDj < numCellTypes; cellIDj++) {
      // other than boundaries, off-diagonals are zero (only cells of same type interact)
      // if either cellID is a boundary, set it to 0 (no attraction with boundaries)
      // otherwise, one of the cellIDs must not be a boundary. If the cells are the same, they should attract (set to 1).
      double boundaryID = numCellTypes - 1;
      double attractionModifier =
          (cellIDi == boundaryID || cellIDj == boundaryID) ? 0.0 : (cellIDi == cellIDj) ? 1.0
                                                                                        : 0.0;
      cell2D.setCellTypeAttractionModifiers(cellIDi, cellIDj, attractionModifier);
    }
  }
  cell2D.printInteractionMatrix();

  // initialize particles with the same number of vertices and the same preferred shape parameter calA0
  cell2D.monodisperse2D(calA0, nv);
  int circleID = 0, rectangleID = 1;
  cell2D.initializeTransverseTissue(phi0, Ftol, circleID);  // initialize within a ring boundary
  cell2D.initializeNeighborLinkedList2D(boxLengthScale);

  // switch ring boundary for rectangular boundary
  // cell2D.replaceCircularBoundary(rectangleID, 2.0);

  // compress to desired density
  bool isFIRE = true;  // use damped NVE to quench
  cell2D.resizeNeighborLinkedList2D();
  cell2D.vertexCompress2Target2D_polygon(attractiveSmoothWithPolyWalls, Ftol, dt0, phi, 2 * dphi0, isFIRE);

  cell2D.replacePolyWallWithDP(numCellTypes);
  cell2D.resizeCatchBonds();
  cell2D.resizeNeighborLinkedList2D();

  double relaxTime = 50.0;
  cell2D.dampedVertexNVE(attractiveSmoothForceUpdate, dt0, relaxTime, 0);
  cell2D.setl00();  // set l00 to be l0 before setting maxwell relaxation time
  cell2D.setActiveBrownianParameters(v0_abp, tau_abp);

  cell2D.dampedVertexNVE(attractionSmoothActiveBrownianCatchBondsUpdate, dt0, relaxTime, 0);
  cell2D.dampedVertexNVE(attractionSmoothActiveBrownianCatchBondsUpdate, dt0, runTime, runTime / 20.0);
  cout << "\n** Finished psm.cpp (2D transverse section of pre-somitic mesoderm), ending. " << endl;

  return 0;
}