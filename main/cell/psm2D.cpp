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
./main/cell/psm2D.o   12   25 1.05 0.01  25.0   0.01  1.0    1   1      400    test1
./main/cell/psm2D.o   12   25 1.05 0.05  25.0   0.01  1.0    1   1      400    test2
./main/cell/psm2D.o   12   25 1.05 0.1   25.0   0.01  1.0    1   1      400    test3
./main/cell/psm2D.o   12   25 1.05 0.2   25.0   0.01  1.0    1   1      400    test4

./main/cell/psm2D.o   12   25 1.05 0.01  25.0   0.05  1.0    1   1      400    test5
./main/cell/psm2D.o   12   25 1.05 0.05  25.0   0.05  1.0    1   1      400    test6
./main/cell/psm2D.o   12   25 1.05 0.1   25.0   0.05  1.0    1   1      400    test7
./main/cell/psm2D.o   12   25 1.05 0.2   25.0   0.05  1.0    1   1      400    test8

./main/cell/psm2D.o   12   25 1.05 0.01  25.0   0.1   1.0    1   1      400    test9
./main/cell/psm2D.o   12   25 1.05 0.05  25.0   0.1   1.0    1   1      400    test10
./main/cell/psm2D.o   12   25 1.05 0.1   25.0   0.1   1.0    1   1      400    test11
./main/cell/psm2D.o   12   25 1.05 0.2   0.0   0.0   1.0    1   1      100    test12

./main/cell/psm2D.o   12   25 1.05 0.0  25.0   0.05   1.0    1   1      400    test1
./main/cell/psm2D.o   12   25 1.05 0.01 25.0   0.05   1.0    1   1      400    test2
./main/cell/psm2D.o   12   25 1.05 0.1  25.0   0.05   1.0    1   1      400    test3
*/
//                  NCELLS NV  A0  att t_maxwell v0  tau_abp sm seed duration outFileStem

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
const double phi0 = 0.91;           // initial packing fraction
const double phiMax = 0.8;
const double smallfrac = 1.0;  // fraction of small particles
const double sizeratio = 1.0;  // size ratio between small and large particles
const double dt0 = 0.01;       // initial magnitude of time step in units of MD time
const double Ptol = 1e-8;
const double Ftol = 1e-12;
const double att_range = 0.3;
const double maxwellRelaxationTime = 10.0;

int main(int argc, char const* argv[]) {
  // local variables to be read in
  int NCELLS, nv, seed, sm;
  double calA0, att, B = 1.0;
  double t_stress, runTime;
  double v0_abp, tau_abp;

  // read in parameters from command line input
  string NCELLS_str = argv[1];
  string nv_str = argv[2];
  string calA0_str = argv[3];
  string att_str = argv[4];
  string t_stress_str = argv[5];
  string v0_str = argv[6];
  string tau_abp_str = argv[7];
  string sm_str = argv[8];
  string seed_str = argv[9];
  string duration_str = argv[10];
  string outFileStem = argv[11];

  string positionFile = outFileStem + ".pos";
  string tissueFile = outFileStem + ".tissue";

  // using sstreams to get parameters
  stringstream NCELLSss(NCELLS_str);
  stringstream nvss(nv_str);
  stringstream calA0ss(calA0_str);
  stringstream attss(att_str);
  stringstream t_stressss(t_stress_str);
  stringstream v0ss(v0_str);
  stringstream tau_abpss(tau_abp_str);
  stringstream smss(sm_str);
  stringstream durationss(duration_str);
  stringstream seedss(seed_str);

  // read into data
  NCELLSss >> NCELLS;
  nvss >> nv;
  calA0ss >> calA0;
  attss >> att;
  t_stressss >> t_stress;
  v0ss >> v0_abp;
  tau_abpss >> tau_abp;
  smss >> sm;
  durationss >> runTime;
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

  cell2D.setB(B);
  if (t_stress > 0.0)
    cell2D.setMaxwellRelaxationTime(maxwellRelaxationTime);  // t_stress is infinity unless this is uncommented
  //  specify non-periodic boundaries
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
  dpmMemFn attractionSmoothWithActiveBrownianUpdate = static_cast<void (dpm::*)()>(&cell::attractiveSmoothForceUpdateWithCrawling);
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
  // cell2D.monodisperseSmooth(calA0, nv);
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
  std::vector<double> savedPositions;

  // dpmMemFn customForceUpdate = repulsivePolarityForceUpdate;
  // dpmMemFn customForceUpdate = attractivePolarityForceUpdate;
  dpmMemFn customForceUpdate;
  if (sm) {
    if (v0_abp <= 0.0)
      customForceUpdate = attractiveSmoothForceUpdate;
    else {
      // assert(false);  // don't have this yet
      customForceUpdate = attractionSmoothWithActiveBrownianUpdate;
      cell2D.setActiveBrownianParameters(v0_abp, tau_abp);
    }
  } else {
    // bumpy
    if (v0_abp <= 0.0)
      customForceUpdate = attractiveForceUpdate;
    else {
      customForceUpdate = attractionWithActiveBrownianUpdate;
      cell2D.setActiveBrownianParameters(v0_abp, tau_abp);
    }
  }
  // cell2D.dampedVertexNVE(attractiveForceUpdateWithPolyWalls, dt0, relaxTimeShort, relaxTimeShort / 2);
  // cell2D.replacePolyWallWithDP(numCellTypes);
  cout << "after replacePolyWallWithDP, about to run NVE for duration " << runTime << "\n";
  cell2D.resizeNeighborLinkedList2D();
  // cell2D.dampedVertexNVE(customForceUpdate, dt0, relaxTime, relaxTime / 15);
  if (v0_abp <= 0.0) {
    // thermal simulation, no activity, no damping
    cell2D.vertexNVE(customForceUpdate, 1e-2, dt0, runTime, runTime / 20.0);
  } else {
    // active simulation, damping
    cell2D.dampedVertexNVE(customForceUpdate, dt0, runTime, runTime / 100.0);
  }
  // cell2D.saveConfiguration(savedPositions);
  // cell2D.loadConfiguration(savedPositions);

  cout
      << "\n** Finished psm.cpp (2D transverse section of pre-somitic mesoderm), ending. " << endl;

  return 0;
}