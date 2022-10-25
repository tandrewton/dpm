// File to compress DPM into confluent epithelial layer, then laser ablate
//
//
// Create bidisperse DPM particles, set constants, place particle centers,
// relax shapes + positions, compress to target phi, delete cells, conduct
// damped MD dynamics.
//
// Compilation command:
// g++ -O3 --std=c++11 -g -I src main/epi2D/laserAblation.cpp src/dpm.cpp src/epi2D.cpp -o main/epi2D/laserAblation.o
//
//  note: if using circular boundaries via polyWall, try pmin = 0.9 and pmax = 0.85 because pmin is soft disk density and pmax is DP preferred density
//
// below: no purse-string, no crawling (inactive simulation)
//./main/epi2D/laserAblation.o 20 20 0 1.10 0.92 0.925 1.0 0.5 0.0 0.01  0.0  4.0  4.0 1.0  0.0  1.0 0.5  0  0   0 1  100  test
// ........................... N  NV Nd A0  pMin  pMax  kl ka att  om   dsq  kps  klp tau dflag  B  Dr0 CIL bound sm sd time file
// below: no purse-string, only crawling
//./main/epi2D/laserAblation.o 40 16 6 0.9 0.94 0.85 1.0 1.0 0.2 0.001  2.0  4.0  4.0 1.0  3.0  1.0 0.5  0  0   0 1 200  test
// ........................... N  NV Nd A0  pMin  pMax  kl ka att  om   dsq  kps  klp tau dflag  B  Dr0 CIL bound sm sd time file
// below: purse-string, no crawling
//./main/epi2D/laserAblation.o 40 20 4 1.10 0.94 0.85 1.0 1.0 0.2 0.001  2.0  1.0  4.0 1.0  0.0  1.0 0.5  0  0   0 1  200  test
// ........................... N  NV Nd A0  pMin  pMax  kl ka att  om   dsq  kps  klp tau dflag  B  Dr0 CIL bound sm sd time file
// below: purse-string, and crawling
//./main/epi2D/laserAblation.o 20 20 4 1.10 0.92 0.865 1.0 1.0 0.1 0.01  2.0  4.0  4.0 1.0  3.0  1.0 0.5  0  1   1 1  110  test
// ........................... N  NV Nd A0  pMin  pMax  kl ka att  om   dsq  kps  klp tau dflag  B  Dr0 CIL bound sm sd time file

// ./main/epi2D/laserAblation.o 40 20 4 1.0 0.92 0.925 1.0 0.5 0.2 0.013  2.0  4.0  4.0 1.0  3.0  1.0 0.5  0  0   0 1  200  test

// bash bash/epi2D/submit_laserAblation.sh 40 20 6 1.10 0.92 0.925 1.0 1.0 0.2 0.01 0.0 4.0 4.0 1.0 3.0 1.0 0.5 0 0 400 pi_ohern,day,scavenge 0-24:00:00 1 1

//
// Parameter input list
// 1. NCELLS: 			number of particles
// 2. nsmall: 			number of vertices on small particles
// 3. ndelete:      number of cells to delete
// 4. calA0: 			  preferred shape parameter for all particles
// 5. phiMin: 			p
// 6. phiMax: 			p
// 7. kl: 				  perimeter spring constant
// 8. ka:           area spring constant
// 9. att:          attraction strength parameter
// 10. omega:        true strain rate for shrinking pursestring segments
// 11. deltaSq:     yield length squared for purse-string springs (in units of vertex diameter squared). auto set to 0 if omega is also 0
// 12. k_ps:        spring constant for purse string virtual particle to wound vertex
// 13. k_lp:        spring constant for flag to nearest vertex on wound edge for crawling
// 14. tau_lp:      protrusion time constant (controls stochastic lifetime of a protrusion)
// 15. d_flag:      protrusion distance from cell edge in units of vertex diameter
// 16. B:           (over)damping coefficient gamma
// 17. Dr0:         rotational diffusion constant for protrusion activity
// 18. boolCIL:     bool for whether cells conduct contact inhibition of locomotion
// 19. boundary     choice of boundary condition (0 = free, 1 = sticky circular boundaries)
// 20. smooth       choice of smooth or bumpy forces (0 = bumpy, 1 = smooth)
// 21. seed: 			  seed for random number generator
// 22. time:        amount of time (tau) to simulate
// 23. outFileStem  stem of output file names, i.e. for "test", energy.test, position.test, etc.

// header files
#include <sstream>
#include "epi2D.h"

// preprocessor macros
#define NDIM 2

using namespace std;

const bool plotCompression = 0;  // whether or not to plot configuration during compression protocol
const double dphi0 = 0.005;      // packing fraction increment
// const double ka = 1.0;              // area force spring constant (should be unit)
double kc = 1.0;                    // interaction force spring constant (should be unit)
const double boxLengthScale = 2.5;  // neighbor list box size in units of initial l0
// const double phi0 = 0.5;            // initial packing fraction
const double smallfrac = 0.5;  // fraction of small particles
const double sizeratio = 1.4;  // size ratio between small and large particles
const double dt0 = 0.04;       // initial magnitude of time step in units of MD time
const double Ptol = 1e-8;
const double Ftol = 1e-12;
const double att_range = 0.3;
bool isPbcOn = false;

int main(int argc, char const* argv[]) {
  // local variables to be read in
  int NCELLS, nsmall, seed, gi, ndelete;
  double calA0, kl, ka = 1.0, kb = 0.0, phiMin, phiMax, att, B, Dr0, time_dbl;
  bool boolCIL, boolBound, boolSmooth, purseStringOn = true;
  double strainRate_ps, k_ps, k_LP, tau_LP, deltaSq, maxProtrusionLength, shapeRelaxationRate;

  // read in parameters from command line input
  string NCELLS_str = argv[1];
  string nsmall_str = argv[2];
  string ndelete_str = argv[3];
  string calA0_str = argv[4];
  string phiMin_str = argv[5];
  string phiMax_str = argv[6];
  string kl_str = argv[7];
  string ka_str = argv[8];
  string att_str = argv[9];
  string omega_str = argv[10];
  string deltaSq_str = argv[11];
  string k_ps_str = argv[12];
  string k_lp_str = argv[13];
  string tau_lp_str = argv[14];
  string d_flag_str = argv[15];
  string B_str = argv[16];
  string Dr0_str = argv[17];
  string boolCIL_str = argv[18];
  string bound_str = argv[19];
  string smooth_str = argv[20];
  string seed_str = argv[21];
  string time_str = argv[22];
  string outFileStem = argv[23];

  string positionFile = outFileStem + ".pos";
  string energyFile = outFileStem + ".energy";
  string stressFile = outFileStem + ".stress";
  string voidFile = outFileStem + ".void";
  string edgeFile = outFileStem + ".edge";
  string purseStringFile = outFileStem + ".purseString";
  string voidAreaFile = outFileStem + ".voidArea";
  string bulkCellShapeFile = outFileStem + ".bulkCellShape";
  string innerCellShapeFile = outFileStem + ".innerCellShape";
  string woundPropertiesFile = outFileStem + ".woundProperties";
  string cellIDFile = outFileStem + ".cellID";

  // using sstreams to get parameters
  stringstream NCELLSss(NCELLS_str);
  stringstream nsmallss(nsmall_str);
  stringstream ndeletess(ndelete_str);
  stringstream calA0ss(calA0_str);
  stringstream phiMinss(phiMin_str);
  stringstream phiMaxss(phiMax_str);
  stringstream klss(kl_str);
  stringstream kass(ka_str);
  stringstream attss(att_str);
  stringstream omegass(omega_str);
  stringstream deltaSqss(deltaSq_str);
  stringstream k_psss(k_ps_str);
  stringstream k_lpss(k_lp_str);
  stringstream tau_lpss(tau_lp_str);
  stringstream d_flagss(d_flag_str);
  stringstream Bss(B_str);
  stringstream Dr0ss(Dr0_str);
  stringstream boolCILss(boolCIL_str);
  stringstream boundss(bound_str);
  stringstream smoothss(smooth_str);
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
  kass >> ka;
  attss >> att;
  omegass >> strainRate_ps;
  deltaSqss >> deltaSq;
  k_psss >> k_ps;
  k_lpss >> k_LP;
  tau_lpss >> tau_LP;
  d_flagss >> maxProtrusionLength;
  Bss >> B;
  Dr0ss >> Dr0;
  boolCILss >> boolCIL;
  boundss >> boolBound;
  smoothss >> boolSmooth;
  seedss >> seed;
  timess >> time_dbl;

  // number of time steps
  int NT = int(time_dbl / dt0);

  const double phi0 = phiMin;

  // double strainRate_ps = 0.01, k_ps = 1.0, k_LP = 2.0, tau_LP = 1.0, deltaSq = 2.0, maxProtrusionLength = 3.0;
  // purse string (ps) strain rate (constant true strain rate)
  // ps spring constant between wound vertex and ps vertex
  // lamellipodia (lp) spring constant between protruded spring anchor (flag) and a colinear vertex
  // lp timescale (10% chance of decaying per tau_LP)
  // squared max distance for wound-ps vertex spring (units of vdiam)
  // max length of lp (units of vdiam)

  if (strainRate_ps < 1e-20) {  // auto set deltasq to 0 if pursestring is disabled, so that all ps springs break off
    cout << "strain rate for pursestring is zero, setting deltaSq yield length to also be zero\n";
    deltaSq = 0.0;
  }
  epi2D epithelial(NCELLS, 0.0, 0.0, Dr0, strainRate_ps, k_ps, k_LP, tau_LP, deltaSq, maxProtrusionLength, seed);

  epithelial.openPosObject(positionFile);
  epithelial.openEnergyObject(energyFile);
  epithelial.openStressObject(stressFile);
  epithelial.openBoundaryObject(voidFile);
  epithelial.openEdgeObject(edgeFile);
  epithelial.openPurseStringObject(purseStringFile);
  epithelial.openVoidAreaObject(voidAreaFile);
  epithelial.openBulkCellShapeObject(bulkCellShapeFile);
  epithelial.openInnerCellShapeObject(innerCellShapeFile);
  epithelial.openWoundPropertiesObject(woundPropertiesFile);
  epithelial.openCellIDObject(cellIDFile);

  // set spring constants
  epithelial.setka(ka);
  epithelial.setkl(kl);
  epithelial.setkb(kb);
  epithelial.setkc(kc);
  epithelial.setShapeRelaxationRate(shapeRelaxationRate);

  // set CIL option
  epithelial.setboolCIL(boolCIL);

  epithelial.setpbc(0, isPbcOn);
  epithelial.setpbc(1, isPbcOn);

  epithelial.setRandPsi();

  // set adhesion scale
  epithelial.setl1(att);
  epithelial.setl2(att_range);
  if (att > att_range) {
    cout << "attraction stronger than attraction range; discontinuous adhesive potential; error, exiting!\n";
    return 1;
  }
  if (fabs(att - att_range) < 1e-5) {
    cout << "or, adhesion lengthscales are the same, so we will divide by zero; error' exiting!\n";
    return 1;
  }

  // set base dpm force, upcast derived epi2D forces
  dpmMemFn repulsiveForceUpdate = &dpm::repulsiveForceUpdate;
  dpmMemFn repulsiveForceUpdateWithWalls = static_cast<void (dpm::*)()>(&epi2D::repulsiveForceUpdateWithWalls);
  dpmMemFn attractiveForceUpdate = static_cast<void (dpm::*)()>(&epi2D::attractiveForceUpdate_2);
  dpmMemFn crawlingWithPSForceUpdate = static_cast<void (dpm::*)()>(&epi2D::crawlingWithPurseString);
  dpmMemFn repulsiveForceUpdateWithCircularAperture = static_cast<void (dpm::*)()>(&epi2D::repulsiveForceWithCircularApertureWall);
  dpmMemFn repulsiveForceUpdateWithCircularWalls = static_cast<void (dpm::*)()>(&epi2D::repulsiveForceUpdateWithPolyWall);
  dpmMemFn attractiveForceUpdateWithCircularWalls = static_cast<void (dpm::*)()>(&epi2D::attractiveForceUpdateWithPolyWall);
  dpmMemFn crawlingWithPSForceUpdateWithCircularWalls = static_cast<void (dpm::*)()>(&epi2D::crawlingWithPurseStringAndCircularWalls);
  dpmMemFn circuloLineAttraction = static_cast<void (dpm::*)()>(&epi2D::attractiveForceUpdate_circulo);
  dpmMemFn crawlingWithPSSmooth = static_cast<void (dpm::*)()>(&epi2D::crawlingWithPurseStringCirculo);
  dpmMemFn crawlingWithPSSmooth_Walls = static_cast<void (dpm::*)()>(&epi2D::crawlingWithPurseStringCirculoWalls);
  dpmMemFn circuloLineAttractionWithCircularWalls = static_cast<void (dpm::*)()>(&epi2D::circuloLineAttractionWithCircularWalls);

  epithelial.monodisperse2D(calA0, nsmall);
  epithelial.initializevnn();

  double boxAspectRatio = 1.0;
  bool setUpCircularBoundary = true;
  // initialize positions and setup polygonal boundary condition if setUpCircularBoundary is enabled
  epithelial.initializePositions2D(phi0, Ftol, false, boxAspectRatio, setUpCircularBoundary);
  epithelial.printConfiguration2D();

  epithelial.initializeNeighborLinkedList2D(boxLengthScale);

  // epithelial.vertexCompress2Target2D(repulsiveForceUpdateWithWalls, Ftol, dt0, phiMax, dphi0);
  //  epithelial.vertexCompress2Target2D(repulsiveForceUpdateWithCircularAperture, Ftol, dt0, phiMax, dphi0);
  if (isPbcOn)
    epithelial.vertexCompress2Target2D_polygon(repulsiveForceUpdate, Ftol, dt0, phiMax, dphi0);
  else
    epithelial.vertexCompress2Target2D_polygon(repulsiveForceUpdateWithCircularWalls, Ftol, dt0, phiMax, dphi0);
  epithelial.printConfiguration2D();

  // after compress, turn on damped NVE
  double T = 1e-10;
  double relaxTime = 10.0;
  assert(relaxTime <= 10.0);  // if relaxTime > 10.0, then dampedNP0 will run over the hardcoded limit. Could pass a parameter to dampedNP0 to tell it how long to wait, or just hardcode for now.
  double printInterval = relaxTime / 2.0;
  double runTime = 25.0;
  epithelial.drawVelocities2D(T);

  dpmMemFn customForceUpdate_inactive;
  dpmMemFn customForceUpdate_active;
  dpmMemFn customForceUpdate_inactive_with_circular_walls;
  dpmMemFn customForceUpdate_active_with_circular_walls;
  if (!boolSmooth) {
    cout << "bumpy forces\n";
    customForceUpdate_inactive = attractiveForceUpdate;
    customForceUpdate_active = crawlingWithPSForceUpdate;
    customForceUpdate_inactive_with_circular_walls = attractiveForceUpdateWithCircularWalls;
    customForceUpdate_active_with_circular_walls = crawlingWithPSForceUpdateWithCircularWalls;
  } else {
    cout << "smooth forces\n";
    customForceUpdate_inactive = circuloLineAttraction;
    customForceUpdate_active = crawlingWithPSSmooth;
    customForceUpdate_inactive_with_circular_walls = circuloLineAttractionWithCircularWalls;
    customForceUpdate_active_with_circular_walls = crawlingWithPSSmooth_Walls;
  }

  if (isPbcOn)
    epithelial.dampedNP0(customForceUpdate_inactive, B, dt0, relaxTime, printInterval);
  else
    epithelial.dampedNP0(customForceUpdate_inactive_with_circular_walls, B, dt0, relaxTime, printInterval);

  // epithelial.dampedNP0(customForceUpdate_inactive_with_circular_walls, B, dt0, runTime, runTime/10.0);

  // LASER ABLATION SCHEME
  double xLoc = 0.0, yLoc = 0.0;
  int numCellsToAblate = ndelete;
  epithelial.laserAblate(numCellsToAblate, sizeratio, nsmall, xLoc, yLoc);
  epithelial.zeroMomentum();

  //  dampedNP0 runs simulation with purse-string integration and crawling

  cout << "boolBound = " << boolBound << '\n';
  // boolBound is used here to do wound healing simulations with walls, simulating a static bulk medium
  if (boolBound) {
    for (int i = 0; i < 2; i++)
      epithelial.dampedNP0(customForceUpdate_inactive_with_circular_walls, B, dt0, relaxTime, printInterval);

    epithelial.dampedNP0(customForceUpdate_active_with_circular_walls, B, dt0, time_dbl, printInterval, purseStringOn);
  } else {
    for (int i = 0; i < 2; i++)
      epithelial.dampedNP0(customForceUpdate_inactive, B, dt0, relaxTime, printInterval);

    epithelial.dampedNP0(customForceUpdate_active, B, dt0, time_dbl, printInterval, purseStringOn);
  }
  cout << "\n** Finished laserAblation.cpp, ending. " << endl;
  return 0;
}