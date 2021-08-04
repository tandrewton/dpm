// File to compress DPM into confluent epithelial layer, then laser ablate
// ** Features to add: purse-string contraction, substrate interaction
//
// Cells are bidisperse, with bending energy (broken) and vertex-vertex attraction options.
//
//
// Create bidisperse DPM particles, set constants, place particle centers,
// relax shapes + positions, compress to target phi, delete cells, conduct
// damped MD dynamics.
//
// Compilation command:
// g++ -O3 --std=c++11 -I src main/epi2D/laserAblation.cpp src/dpm.cpp src/epi2D.cpp -o main/epi2D/laserAblation.o
// ./main/epi2D/laserAblation.o 48 24 10 1.08 0.2 0.85 1.0 0.3 0.5 1.0 0.5 1 1 1000 pos.test energy.test stress.test
// ./main/epi2D/laserAblation.o 24 24 5 1.08 0.2 0.85 1.0 0.1 0 1.0 0.5 1 1 1000 pos.test energy.test stress.test
// ./main/epi2D/laserAblation.o 24 24 5 1.08 0.2 0.85 1.0 0 0 1.0 0.1 0 1 100 pos.test energy.test stress.test
// ./main/epi2D/laserAblation.o 24 24 5 1.08 0.2 0.85 1.0 0.3 0 1.0 0.1 0 1 100 pos.test energy.test stress.test
// ./main/epi2D/laserAblation.o 24 24 5 1.08 0.2 0.85 1.0 0.25 0.5 1.0 0.5 0 1 100 pos.test energy.test stress.test
//
// Parameter input list
// 1. NCELLS: 			number of particles
// 2. nsmall: 			number of vertices on small particles (larger
// particles set by 1.4:1.0 bidispersity)
// 3. ndelete:      number of cells to delete
// 4. calA0: 			  preferred shape parameter for all particles
// 5. phiMin: 			p ?
// 6. phiMax: 			p
// 7. kl: 				  perimeter spring constant
// 8. att:          attraction strength parameter
// 9. v0:           active vertex velocity scale
// 10. B:           (over)damping coefficient gamma
// 11. Dr0:         rotational diffusion constant for protrusion activity
// 12. boolCIL:     bool for whether cells conduct contact inhibition of locomotion
// 13. seed: 			  seed for random number generator
// 14. time:        amount of time (tau) to simulate
// 15. positionFile: 	string of path to output file with position/configuration data
// 16. energyFile:  string of path to output file with energy data
// 17. stressFile:  string of path to output file with stress data

// header files
#include <sstream>
#include "epi2D.h"

// preprocessor macros
#define NDIM 2

// namspace
using namespace std;

const bool plotCompression = 0;     // whether or not to plot configuration during compression protocol
const double dphi0 = 0.005;         // packing fraction increment
const double ka = 1.0;              // area force spring constant (should be unit)
const double kc = 1.0;              // interaction force spring constant (should be unit)
const double boxLengthScale = 2.5;  // neighbor list box size in units of initial l0
//const double phi0 = 0.5;            // initial packing fraction
const double smallfrac = 0.5;  // fraction of small particles
const double sizeratio = 1.4;  // size ratio between small and large particles
const double dt0 = 0.01;       // initial magnitude of time step in units of MD time
const double Ptol = 1e-8;
const double Ftol = 1e-12;
const double att_range = 0.3;

int main(int argc, char const* argv[]) {
  // local variables to be read in
  int NCELLS, nsmall, seed, gi, ndelete;
  double calA0, kl, kb = 0.0, phiMin, phiMax, att, v0, B, Dr0, time_dbl;
  bool boolCIL;

  // read in parameters from command line input
  string NCELLS_str = argv[1];
  string nsmall_str = argv[2];
  string ndelete_str = argv[3];
  string calA0_str = argv[4];
  string phiMin_str = argv[5];
  string phiMax_str = argv[6];
  string kl_str = argv[7];
  string att_str = argv[8];
  string v0_str = argv[9];
  string B_str = argv[10];
  string Dr0_str = argv[11];
  string boolCIL_str = argv[12];
  string seed_str = argv[13];
  string time_str = argv[14];
  string positionFile = argv[15];
  string energyFile = argv[16];
  string stressFile = argv[17];

  // using sstreams to get parameters
  stringstream NCELLSss(NCELLS_str);
  stringstream nsmallss(nsmall_str);
  stringstream ndeletess(ndelete_str);
  stringstream calA0ss(calA0_str);
  stringstream phiMinss(phiMin_str);
  stringstream phiMaxss(phiMax_str);
  stringstream klss(kl_str);
  stringstream attss(att_str);
  stringstream v0ss(v0_str);
  stringstream Bss(B_str);
  stringstream Dr0ss(Dr0_str);
  stringstream boolCILss(boolCIL_str);
  stringstream seedss(seed_str);
  stringstream timess(time_str);

  // read into data
  NCELLSss >> NCELLS;
  nsmallss >> nsmall;
  ndeletess >> ndelete;
  calA0ss >> calA0;
  phiMinss >> phiMin;
  phiMaxss >> phiMax;
  klss >> kl;
  attss >> att;
  v0ss >> v0;
  Bss >> B;
  Dr0ss >> Dr0;
  boolCILss >> boolCIL;
  seedss >> seed;
  timess >> time_dbl;

  // number of time steps
  int NT = int(time_dbl / dt0);

  const double phi0 = phiMin;

  // instantiate object
  epi2D epithelial(NCELLS, 0.0, 0.0, Dr0, seed);

  epithelial.openPosObject(positionFile);
  epithelial.openEnergyObject(energyFile);
  epithelial.openStressObject(stressFile);

  // set spring constants
  epithelial.setka(ka);
  epithelial.setkl(kl);
  epithelial.setkb(kb);
  epithelial.setkc(kc);

  //set activity scale, and CIL option
  epithelial.setv0(v0);
  epithelial.setboolCIL(boolCIL);
  epithelial.setpbc(0, false);
  epithelial.setpbc(1, false);

  //set adhesion scale
  epithelial.setl1(att);
  epithelial.setl2(att_range);
  if (att > att_range) {
    cout << "attraction stronger than attraction range; discontinuous adhesive potential; error, exiting!\n";
    return 1;
  }

  epithelial.monodisperse2D(calA0, nsmall);

  epithelial.initializePositions2D(phi0, Ftol, true);
  epithelial.printConfiguration2D();

  epithelial.initializeNeighborLinkedList2D(boxLengthScale);

  // set base dpm force, upcast derived epi2D forces
  dpmMemFn repulsiveForceUpdate = &dpm::repulsiveForceUpdate;
  dpmMemFn repulsiveForceUpdateWithWalls = static_cast<void (dpm::*)()>(&epi2D::repulsiveForceUpdateWithWalls);
  dpmMemFn attractiveForceUpdate = static_cast<void (dpm::*)()>(&epi2D::attractiveForceUpdate_2);
  dpmMemFn activeForceUpdate = static_cast<void (dpm::*)()>(&epi2D::activeAttractiveForceUpdate);

  epithelial.vertexCompress2Target2D(repulsiveForceUpdateWithWalls, Ftol, dt0, phiMax, dphi0);
  epithelial.printConfiguration2D();

  //after compress, turn on damped NVE
  double T = 1e-4;
  epithelial.drawVelocities2D(T);
  epithelial.dampedNP0(attractiveForceUpdate, B, dt0, 10, 0);

  // LASER ABLATION SCHEME
  double xLoc = 0.0, yLoc = 0.0;
  int numCellsToAblate = ndelete;
  epithelial.laserAblate(numCellsToAblate, sizeratio, nsmall, xLoc, yLoc);
  epithelial.setRandPsi();
  epithelial.zeroMomentum();

  epithelial.dampedNP0(activeForceUpdate, B, dt0, time_dbl, time_dbl / 20);

  cout << "\n** Finished laserAblation.cpp, ending. " << endl;
  return 0;
}