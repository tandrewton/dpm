// simulate neural tube aspect of convergent extension
// WALL WALL WALL WALL WALL WALL
// WALL ||  DPM DPM DPM  || WALL
// WALL ||  DPM DPM DPM  || WALL
// WALL ||  DPM DPM DPM  || WALL
// WALL ||  DPM DPM DPM  || WALL
// WALL ||  DPM DPM DPM  || WALL
// WALL ||  DPM DPM DPM  || WALL
// .............................
// .............................
// .............................

// Compilation command:
// g++ -O3 --std=c++11 -g -I src main/cell/neuralTube.cpp src/dpm.cpp src/cell.cpp -o main/cell/NT.o
// run command:
// ./main/cell/NT.o    10   20 1.0 0.01      -0.01          0.002         0.002        1    30     test
//                  NCELLS NV  A0  att initial_pressure pressure_rate adhesion_rate  seed runtime outFileStem

#include <sstream>
#include "cell.h"

#define NDIM 2
using namespace std;

// global constants
const bool plotCompression = 0;     // whether or not to plot configuration during compression protocol (0 saves memory)
const double dphi0 = 0.005;         // packing fraction increment
const double ka = 1.0;              // area force spring constant (should be unit)
const double kc = 1.0;              // interaction force spring constant (should be unit)
const double kb = 0.0;              // bending energy spring constant (should be zero)
const double kl = 1.0;              // segment length interaction force (should be unit)
const double boxLengthScale = 2.5;  // neighbor list box size in units of initial l0
const double phi0 = 0.81;           // initial packing fraction
const double phiMax = 0.815;
const double smallfrac = 1.0;  // fraction of small particles
const double sizeratio = 1.0;  // size ratio between small and large particles
const double dt0 = 1e-2;       // initial magnitude of time step in units of MD time
const double Ptol = 1e-8;
const double Ftol = 1e-12;
const double att_range = 0.5;

int main(int argc, char const* argv[]) {
  // local variables to be read in
  int NCELLS, nv, seed;
  double calA0, att;
  double pressureRate, adhesionRate, appliedUniaxialPressure, runTime;

  // read in parameters from command line input
  string NCELLS_str = argv[1];
  string nv_str = argv[2];
  string calA0_str = argv[3];
  string att_str = argv[4];
  string init_pressure_str = argv[5];
  string pressure_rate_str = argv[6];
  string adhesion_rate_str = argv[7];
  string seed_str = argv[8];
  string runtime_str = argv[9];
  string outFileStem = argv[10];

  string positionFile = outFileStem + ".pos";

  // using sstreams to get parameters
  stringstream NCELLSss(NCELLS_str);
  stringstream nvss(nv_str);
  stringstream calA0ss(calA0_str);
  stringstream attss(att_str);
  stringstream init_pressuress(init_pressure_str);
  stringstream pressure_ratess(pressure_rate_str);
  stringstream adh_ratess(adhesion_rate_str);
  stringstream seedss(seed_str);
  stringstream runtimess(runtime_str);

  // read into data
  NCELLSss >> NCELLS;
  nvss >> nv;
  calA0ss >> calA0;
  attss >> att;
  init_pressuress >> appliedUniaxialPressure;
  pressure_ratess >> pressureRate;
  adh_ratess >> adhesionRate;
  seedss >> seed;
  runtimess >> runTime;

  cell cell2D(NCELLS, 0.0, 0.0, seed);
  cell2D.openPosObject(positionFile);

  // set spring constants
  cell2D.setka(ka);
  cell2D.setkl(kl);
  cell2D.setkb(kb);
  cell2D.setkc(kc);

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

  // initialize particles with the same number of vertices and the same preferred shape parameter calA0
  cell2D.monodisperse2D(calA0, nv);

  // initialize particle positions
  double aspectRatioBox = 2.0;
  cell2D.initializePositions2D(phi0, Ftol, false, aspectRatioBox);


  // initialize neighbor linked list
  cell2D.initializeNeighborLinkedList2D(boxLengthScale);

  dpmMemFn repulsiveForceUpdate = &dpm::repulsiveForceUpdate;
  dpmMemFn repulsiveForceUpdateWithWalls = static_cast<void (dpm::*)()>(&cell::repulsiveForceUpdateWithWalls);
  dpmMemFn attractiveForceUpdate = static_cast<void (dpm::*)()>(&cell::attractiveForceUpdate);
  dpmMemFn repulsivePolarityForceUpdate = static_cast<void (dpm::*)()>(&cell::repulsiveWithPolarityForceUpdate);
  dpmMemFn attractivePolarityForceUpdate = static_cast<void (dpm::*)()>(&cell::attractiveWithPolarityForceUpdate);
  dpmMemFn attractiveWallCrawlingPolarityForceUpdate = static_cast<void (dpm::*)()>(&cell::attractiveWithPolarityForceAndWallCrawlingUpdate);
  dpmMemFn attractiveWallCrawlingForceUpdate = static_cast<void (dpm::*)()>(&cell::attractiveWallCrawlingForceUpdate);

  // compress to target packing fraction
  cell2D.vertexCompress2Target2D(repulsiveForceUpdateWithWalls, Ftol, dt0, phiMax, dphi0);
  cout << "done compressing to target packing fraction\n";

  bool wallsBool = true;
  double relaxTime = 200.0;
  // double runTime = 100.0;
  double B = 1.0;

  // dpmMemFn customForceUpdate = repulsivePolarityForceUpdate;
  dpmMemFn customForceUpdate = attractivePolarityForceUpdate;
  dpmMemFn customActiveForceUpdate = attractiveWallCrawlingPolarityForceUpdate;
  dpmMemFn customActiveForceUpdateWithoutPolarity = attractiveWallCrawlingForceUpdate;
  // dpmMemFn customForceUpdate = attractiveForceUpdate;

  // release top wall and let system equilibrate
  // note: repulsiveForceUpdateWithWalls is not used here because we want to control which walls are on,
  //       as well as compute forces on the walls in case we want to integrate wall dynamics too
  double printIntervalNone = 0.0;
  double adhesionRateNone = 0.0;
  double pressureRateNone = 0.0;
  cell2D.simulateDampedWithWalls(customForceUpdate, B, dt0, relaxTime, printIntervalNone, pressureRateNone, adhesionRateNone, wallsBool, true, true, true, true);
  cout << "done equilibrating under 4 dynamic walls\n";

  // equilibrate under 3 fixed walls and 1 dynamic wall
  cell2D.simulateDampedWithWalls(customForceUpdate, B, dt0, relaxTime, printIntervalNone, pressureRateNone, adhesionRateNone, wallsBool, false, false, true, false);
  cout << "done equilibrating under 1 dynamic wall\n";

  // begin lateral-medial compression
  // might want to save the initial box length
  // slow is 0.0001, runTime = 1000
  // cell2D.simulateDampedWithWalls(customForceUpdate, B, dt0, runTime, runTime / 40.0, wallsBool, true, false, false, false, 0.001, 0.0, 0.0);

  // begin uniaxial pressure simulation (P can be positive or negative, applied force to boundary is P * L)
  // top is fixed wall, sides are dynamic with applied pressure, bottom is dynamic with no pressure
  // appliedUniaxialPressure = 1 is too large, +/- 1e-2 is appropriate with timescale 200
  //double appliedUniaxialPressure = -1e-2; // note that the particle-wall interactions might overpower the pressure
  //double pressureRate = 2e-3; // exponential growth rates of pressure and adhesion
  //double adhesionRate = 2e-3;

  // note: this is the active part of the simulation, where the activeForceUpdate has crawling, and there could be a command-line pressure and adhesion rate.
  cell2D.simulateDampedWithWalls(customActiveForceUpdateWithoutPolarity, B, dt0, runTime, runTime / 40.0, pressureRate, adhesionRate, wallsBool, true, true, true, false, 0.0, 0.0, appliedUniaxialPressure);

  cout
      << "\n** Finished neuralTube.cpp, ending. " << endl;
  return 0;
}