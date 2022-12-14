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

// ./main/cell/psm2D.o   24   25 1.05 0.05 0.1  1.0    0    1    test1
// ./main/cell/psm2D.o   24   25 1.05 0.1  0.1  1.0    0    1    test2
// ./main/cell/psm2D.o   24   25 1.05 0.15 0.1  1.0    0    1    test3
// ./main/cell/psm2D.o   24   25 1.05 0.2  0.1  1.0    0    1    test4

// ./main/cell/psm2D.o   24   25 1.05 0.05 0.5  1.0    0    1    test5
// ./main/cell/psm2D.o   24   25 1.05 0.1  0.5  1.0    0    1    test6
// ./main/cell/psm2D.o   24   25 1.05 0.15 0.5  1.0    0    1    test7
// ./main/cell/psm2D.o   24   25 1.05 0.2  0.5  1.0    0    1    test8

// ./main/cell/psm2D.o   24   25 1.05 0.05 1.0  1.0    0    1    test9
// ./main/cell/psm2D.o   24   25 1.05 0.1  1.0  1.0    0    1    test10
// ./main/cell/psm2D.o   24   25 1.05 0.15 1.0  1.0    0    1    test11
// ./main/cell/psm2D.o   24   25 1.05 0.2  1.0  1.0    0    1    test12
//                     NCELLS NV  A0  att   v0  tau_abp sm  seed outFileStem

#include <sstream>
#include "cell.h"

#define NDIM 2
using namespace std;

// global constants
const bool plotCompression = 0;     // whether or not to plot configuration during compression protocol (0 saves memory)
const double dphi0 = 0.005;         // packing fraction increment
const double ka = 1.0;              // area force spring constant (should be unit)
const double kc = 1.0;              // interaction force spring constant (should be unit)
const double kb = 0.01;             // bending energy spring constant (should be zero)
const double kl = 1.0;              // segment length interaction force (should be unit)
const double boxLengthScale = 2.5;  // neighbor list box size in units of initial l0
const double phi0 = 0.7;            // initial packing fraction
const double phiMax = 0.6;
const double smallfrac = 1.0;  // fraction of small particles
const double sizeratio = 1.0;  // size ratio between small and large particles
const double dt0 = 0.05;       // initial magnitude of time step in units of MD time
const double Ptol = 1e-8;
const double Ftol = 1e-12;
const double att_range = 0.3;

int main(int argc, char const* argv[]) {
  // local variables to be read in
  int NCELLS, nv, seed, sm;
  double calA0, att;
  double v0_abp, tau_abp;

  // read in parameters from command line input
  string NCELLS_str = argv[1];
  string nv_str = argv[2];
  string calA0_str = argv[3];
  string att_str = argv[4];
  string v0_str = argv[5];
  string tau_abp_str = argv[6];
  string sm_str = argv[7];
  string seed_str = argv[8];
  string outFileStem = argv[9];

  string positionFile = outFileStem + ".pos";
  string tissueFile = outFileStem + ".tissue";

  // using sstreams to get parameters
  stringstream NCELLSss(NCELLS_str);
  stringstream nvss(nv_str);
  stringstream calA0ss(calA0_str);
  stringstream attss(att_str);
  stringstream v0ss(v0_str);
  stringstream tau_abpss(tau_abp_str);
  stringstream smss(sm_str);
  stringstream seedss(seed_str);

  // read into data
  NCELLSss >> NCELLS;
  nvss >> nv;
  calA0ss >> calA0;
  attss >> att;
  v0ss >> v0_abp;
  tau_abpss >> tau_abp;
  smss >> sm;
  seedss >> seed;

  int numCellTypes = 2;  // 0 = interior cell type (PSM) and 1 = exterior cell type (boundary)
  cell cell2D(NCELLS, 0.0, 0.0, seed, numCellTypes);
  cell2D.openPosObject(positionFile);
  cell2D.openTissueObject(tissueFile);

  // set spring constants
  cell2D.setka(ka);
  cell2D.setkl(kl);
  cell2D.setkb(kb);
  cell2D.setkc(kc);
  cout << "ka, kl, kb, kc = " << ka << '\t' << kl << '\t' << kb << '\t' << kc << '\n';

  // specify non-periodic boundaries
  cell2D.setpbc(0, false);
  cell2D.setpbc(1, false);

  // set adhesion scale
  cell2D.setl1(att);
  cell2D.setl2(att_range);
  if (att > att_range) {
    cout << "attraction stronger than attraction range; discontinuous adhesive potential; error, exiting!\n";
    return 1;
  }
  if (fabs(att - att_range) < 1e-5) {
    cout << "or, adhesion lengthscales are the same, so we will divide by zero; error' exiting!\n";
    return 1;
  }

  dpmMemFn repulsiveForceUpdate = &dpm::repulsiveForceUpdate;
  dpmMemFn repulsiveForceUpdateWithWalls = static_cast<void (dpm::*)()>(&cell::repulsiveForceUpdateWithWalls);
  dpmMemFn attractiveForceUpdate = static_cast<void (dpm::*)()>(&cell::attractiveForceUpdate);
  dpmMemFn attractionWithActiveBrownianUpdate = static_cast<void (dpm::*)()>(&cell::attractiveForceUpdateWithCrawling);
  // dpmMemFn attractionSmoothWithActiveBrownianUpdate = static_cast<void (dpm::*)()>(&cell::attractiveSmoothForceUpdateWithCrawling);
  dpmMemFn attractiveSmoothForceUpdate = static_cast<void (dpm::*)()>(&cell::attractiveSmoothForceUpdate);
  dpmMemFn repulsivePolarityForceUpdate = static_cast<void (dpm::*)()>(&cell::repulsiveWithPolarityForceUpdate);
  dpmMemFn attractivePolarityForceUpdate = static_cast<void (dpm::*)()>(&cell::attractiveWithPolarityForceUpdate);
  dpmMemFn repulsiveForceUpdateWithPolyWalls = static_cast<void (dpm::*)()>(&cell::repulsiveForceUpdateWithPolyWall);
  dpmMemFn attractiveForceUpdateWithPolyWalls = static_cast<void (dpm::*)()>(&cell::attractiveForceUpdateWithPolyWall);

  // cellTypeIntMat = ones to begin with
  for (int i = 0; i < numCellTypes; i++) {
    for (int j = 0; j < numCellTypes; j++) {
      if (i == numCellTypes - 1 || j == numCellTypes - 1) {  // boundaries interact with everyone
        cell2D.setCellTypeAttractionModifiers(i, j, 1.0);
      } else if (i != j) {  // other than boundaries, off diagonals are zero (only cells of same type interact)
                            // else // no cells interact
        cell2D.setCellTypeAttractionModifiers(i, j, 1.0);
      }
    }
  }
  cell2D.printInteractionMatrix();

  // initialize particles with the same number of vertices and the same preferred shape parameter calA0
  cell2D.monodisperse2D(calA0, nv);
  // initialize particle positions
  cell2D.initializeTransverseTissue(phi0, Ftol);
  cell2D.printConfiguration2D();

  cell2D.initializeNeighborLinkedList2D(boxLengthScale);
  cell2D.printConfiguration2D();

  // compress to target packing fraction
  cell2D.vertexCompress2Target2D_polygon(attractiveForceUpdateWithPolyWalls, Ftol, dt0, phiMax, dphi0);
  cout << "done compressing to target packing fraction\n";
  cell2D.printConfiguration2D();

  bool wallsBool = true;
  double relaxTimeShort = 5.0;
  double relaxTime = 100.0;
  double runTime = 400.0;
  double B = 1.0;
  std::vector<double> savedPositions;

  // dpmMemFn customForceUpdate = repulsivePolarityForceUpdate;
  // dpmMemFn customForceUpdate = attractivePolarityForceUpdate;
  dpmMemFn customForceUpdate;
  if (sm) {
    if (v0_abp <= 0.0)
      customForceUpdate = attractiveSmoothForceUpdate;
    else
      assert(false);  // don't have this yet
                      // customForceUpdate = attractionSmoothWithActiveBrownianUpdate;
  } else {
    // bumpy
    if (v0_abp <= 0.0)
      customForceUpdate = attractiveForceUpdate;
    else {
      customForceUpdate = attractionWithActiveBrownianUpdate;
      cell2D.setActiveBrownianParameters(v0_abp, tau_abp);
    }
  }
  cell2D.dampedVertexNVE(attractiveForceUpdateWithPolyWalls, B, dt0, relaxTimeShort, relaxTimeShort / 2);
  cell2D.replacePolyWallWithDP(numCellTypes);
  cout << "after replacePolyWallWithDP\n";
  // cell2D.dampedVertexNVE(customForceUpdate, B, dt0, relaxTime, relaxTime / 15);
  if (v0_abp <= 0.0) {
    // thermal simulation, no activity, no damping
    cell2D.vertexNVE(customForceUpdate, 1e-2, dt0, runTime, runTime / 50.0);
  } else {
    // active simulation, damping
    cell2D.dampedVertexNVE(customForceUpdate, B, dt0, runTime, runTime / 50.0);
  }
  // cell2D.saveConfiguration(savedPositions);
  // cell2D.loadConfiguration(savedPositions);

  /*cell2D.dampedVertexNVE(customForceUpdate, B, dt0, relaxTimeShort, relaxTimeShort/5);
  cell2D.loadConfiguration(savedPositions);

  cell2D.setCellTypeAttractionModifiers(0,3,10.0);

  cell2D.dampedVertexNVE(customForceUpdate, B, dt0, relaxTimeShort, relaxTimeShort/5);
  cell2D.loadConfiguration(savedPositions);
  */
  cout
      << "\n** Finished transverseIntertissue.cpp, ending. " << endl;

  return 0;
}