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
// g++ -O3 --std=c++11 -g -I src main/epi2D/laserAblation.cpp src/dpm.cpp src/epi2D.cpp -o main/epi2D/laserAblation.o

// ./main/epi2D/laserAblation.o 48 20 10 1.08 0.85 0.855 1.0 0.3 0.5 1.0 0.5 0 1 100 pos.test energy.test stress.test void.test

//./main/epi2D/laserAblation.o 20 20 4 1.08 0.92 0.925 1.0 0.3 0.01 2 1.0 2.0 1.0 3 1.0 0.5 0 1 300 pos.test energy.test stress.test void.test

//
// Parameter input list
// 1. NCELLS: 			number of particles
// 2. nsmall: 			number of vertices on small particles (larger
// particles set by 1.4:1.0 bidispersity)
// 3. ndelete:      number of cells to delete
// 4. calA0: 			  preferred shape parameter for all particles
// 5. phiMin: 			p
// 6. phiMax: 			p
// 7. kl: 				  perimeter spring constant
// 8. att:          attraction strength parameter
// 9. omega:        true strain rate for shrinking pursestring segments
// 10. deltaSq:     yield length squared for purse-string springs (in units of vertex diameter squared)
// 11. k_ps:        spring constant for purse string virtual particle to wound vertex
// 12. k_lp:        spring constant for flag to nearest vertex on wound edge for crawling
// 13. tau_lp:      protrusion time constant (controls stochastic lifetime of a protrusion)
// 14. d_flag:      protrusion distance from cell edge in units of vertex diameter
// 15. B:           (over)damping coefficient gamma
// 16. Dr0:         rotational diffusion constant for protrusion activity
// 17. boolCIL:     bool for whether cells conduct contact inhibition of locomotion
// 18. seed: 			  seed for random number generator
// 19. time:        amount of time (tau) to simulate
// 20. positionFile: 	string of path to output file with position/configuration data
// 21. energyFile:  string of path to output file with energy data
// 22. stressFile:  string of path to output file with stress data
// 23. voidFile:    string of path to output file with void location data

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
// const double phi0 = 0.5;            // initial packing fraction
const double smallfrac = 0.5;  // fraction of small particles
const double sizeratio = 1.4;  // size ratio between small and large particles
const double dt0 = 0.01;       // initial magnitude of time step in units of MD time
const double Ptol = 1e-8;
const double Ftol = 1e-12;
const double att_range = 0.3;
bool wallsOn = true;
bool wallsOff = false;

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
  string voidFile = argv[18];
  string edgeFile = "edge.test";
  string purseStringFile = "purseString.test";

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

  double strainRate_ps = 0.01, k_ps = 1.0, k_LP = 2.0, tau_LP = 1.0, deltaSq = 2.0, maxProtrusionLength = 3.0;
  // purse string (ps) strain rate (constant true strain rate)
  // ps spring constant between wound vertex and ps vertex
  // lamellipodia (lp) spring constant between protruded spring anchor (flag) and a colinear vertex
  // lp timescale (10% chance of decaying per tau_LP)
  // squared max distance for wound-ps vertex spring (units of vdiam)
  // max length of lp (units of vdiam)
  // instantiate object

  epi2D epithelial(NCELLS, 0.0, 0.0, Dr0, strainRate_ps, k_ps, k_LP, tau_LP, deltaSq, maxProtrusionLength, seed);

  epithelial.openPosObject(positionFile);
  epithelial.openEnergyObject(energyFile);
  epithelial.openStressObject(stressFile);
  epithelial.openBoundaryObject(voidFile);
  epithelial.openEdgeObject(edgeFile);
  epithelial.openPurseStringObject(purseStringFile);

  // set spring constants
  epithelial.setka(ka);
  epithelial.setkl(kl);
  epithelial.setkb(kb);
  epithelial.setkc(kc);
  epithelial.setkL(double(kl / nsmall));  // calculate the energy kL such that the average force on a vertex due to kL xis equal to the average force due to kl

  // set activity scale, and CIL option
  epithelial.setv0(v0);
  epithelial.setboolCIL(boolCIL);
  epithelial.setpbc(0, false);
  epithelial.setpbc(1, false);
  epithelial.setRandPsi();

  // set adhesion scale
  epithelial.setl1(att);
  epithelial.setl2(att_range);
  if (att > att_range) {
    cout << "attraction stronger than attraction range; discontinuous adhesive potential; error, exiting!\n";
    return 1;
  }

  epithelial.monodisperse2D(calA0, nsmall);
  epithelial.initializevnn();

  epithelial.initializePositions2D(phi0, Ftol, false);
  // epithelial.printConfiguration2D();

  epithelial.initializeNeighborLinkedList2D(boxLengthScale);

  // set base dpm force, upcast derived epi2D forces
  dpmMemFn repulsiveForceUpdate = &dpm::repulsiveForceUpdate;
  dpmMemFn repulsiveForceUpdateWithWalls = static_cast<void (dpm::*)()>(&epi2D::repulsiveForceUpdateWithWalls);
  dpmMemFn attractiveForceUpdate = static_cast<void (dpm::*)()>(&epi2D::attractiveForceUpdate_2);
  dpmMemFn activeForceUpdate = static_cast<void (dpm::*)()>(&epi2D::activeAttractiveForceUpdate);
  dpmMemFn substrateAdhesionForceUpdate = static_cast<void (dpm::*)()>(&epi2D::substrateadhesionAttractiveForceUpdate);
  dpmMemFn repulsiveForceUpdateWithCircularAperture = static_cast<void (dpm::*)()>(&epi2D::repulsiveForceWithCircularApertureWall);

  epithelial.vertexCompress2Target2D(repulsiveForceUpdateWithWalls, Ftol, dt0, phiMax, dphi0);
  // epithelial.vertexCompress2Target2D(repulsiveForceUpdateWithCircularAperture, Ftol, dt0, phiMax, dphi0);
  // epithelial.printConfiguration2D();

  // after compress, turn on damped NVE
  double T = 1e-4;
  double relaxTime = 10.0;
  epithelial.drawVelocities2D(T);
  epithelial.dampedNP0(attractiveForceUpdate, B, dt0, relaxTime, 0, wallsOff);
  // double boxSizeMultiplier = 1.2;
  // epithelial.expandBoxAndCenterParticles(boxSizeMultiplier, boxLengthScale);

  // LASER ABLATION SCHEME
  double xLoc = 0.0, yLoc = 0.0;
  int numCellsToAblate = ndelete;
  epithelial.laserAblate(numCellsToAblate, sizeratio, nsmall, xLoc, yLoc);
  epithelial.zeroMomentum();

  epithelial.dampedNVE2D(attractiveForceUpdate, B, dt0, relaxTime, 0);

  //  dampedNP0 already takes care of purse-string. might want to separate, or just change spring constant
  epithelial.dampedNP0(substrateAdhesionForceUpdate, B, dt0, time_dbl, time_dbl / 100.0, wallsOff);
  // epithelial.dampedNP0(attractiveForceUpdate, B, dt0, time_dbl, time_dbl / 100.0, wallsOff);

  cout << "\n** Finished laserAblation.cpp, ending. " << endl;
  return 0;
}