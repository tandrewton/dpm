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
//./main/epi2D/laserAblation.o 20 20 0 1.10 0.92 0.925 1.0 4.0 0.01 0.0 1.0  0.0  4.0  4.0 1.0  0.0   10.0     0    0 1  100  test
// ........................... N  NV Nd A0  pMin  pMax  kl  ka  kb  att  om   dsq  kps  klp tau dflag t_stress tauRatio sm sd time file
// below: no purse-string, only crawling
//./main/epi2D/laserAblation.o 50 30 3 1.05 0.94 0.85 1.0 1.0 0.01 0.2 1.0  0.0  1.0  4.0 1.0  3.0     25.0     0   0 1  500  test
// ........................... N  NV Nd A0  pMin  pMax  kl ka  kb  att  om   dsq  kps  klp tau dflag  t_stress tauRatio sm sd time file
// below: purse-string, no crawling
//./main/epi2D/laserAblation.o 50 30  5 1.20 0.94 0.99 1.0 16.0 0.01 0.1 1.0  4.0  4.0  4.0 1.0  0.0   76.8  0  1  1 500  test
// ........................... N  NV Nd A0  pMin  pMax  kl  ka   kb  att  om  dsq  kps  klp tau dflag t_stress tauRatio sm sd time file
// below: purse-string, and crawling
//./main/epi2D/laserAblation.o 20 20 4 1.10 0.92 0.865 1.0 4.0 0.01 0.1 1.0   2.0  4.0  4.0 1.0  3.0     10.0    0   0 1 1  110  test
// ........................... N  NV Nd A0  pMin  pMax  kl ka  kb  att  om   dsq  kps  klp tau dflag  t_stress tauRatio sm sd time file

// bash bash/epi2D/submit_laserAblation.sh 40 20 6 1.10 0.92 0.925 1.0 1.0 0.2 1.0  0.0 4.0 4.0 1.0 3.0 1.0 0.5 0 0 400 pi_ohern,day,scavenge 0-24:00:00 1 1
/*
seeds=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20)
for sd in ${seeds[@]}; do
  echo "./main/epi2D/laserAblation.o 30 20 3 1.20 0.94 0.99 1.0 1.0 0.01 0.1 1.0 4.0 4.0 4.0 1.0 0.0 76.8 0 1 $sd 40 test > out_$sd.txt"
done


*/

//
// Parameter input list
// 1. NCELLS: 			number of DPs
// 2. nsmall: 			minimum number of vertices (small DPs)
// 3. ndelete:      number of DPs to delete
// 4. calA0: 			  preferred shape parameter for all DPs
// 5. phiMin: 			packing fraction initial
// 6. phiMax: 			packing fraction final after compression
// 7. kl: 				  perimeter spring constant
// 8. ka:           area spring constant
// 9. kb:           bending spring constant
// 10. att:         attraction strength parameter
// 11. omega:       true strain rate for shrinking pursestring segments
// 12. deltaSq:     yield length squared for purse-string springs (units=vertex diameter squared). auto 0 if omega is 0
// 13. k_ps:        spring constant for purse string virtual particle to wound vertex
// 14. k_lp:        spring constant for flag to nearest vertex on wound edge for crawling
// 15. tau_lp:      protrusion time constant (controls stochastic lifetime of a protrusion)
// 16. d_flag:      protrusion distance from cell edge in units of vertex diameter
// 17. tau_maxwell: stress relaxation timescale for preferred lengths
// 18. tauRatio     ad-hoc parameter to achieve relaxation of cell shapes in the wing disc case.
// 19. smooth       choice of smooth or bumpy forces (0 = bumpy, 1 = smooth)
// 20. seed: 			  seed for random number generator
// 21. time:        amount of time (tau) to simulate
// 22. outFileStem  stem of output file names, i.e. for "test", energy.test, position.test, etc.

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
const double sizeratio = 1.0;  // size ratio between small and large particles
const double dt0 = 0.05;       // initial magnitude of time step in units of MD time
const double Ptol = 1e-8;
const double Ftol = 1e-12;
const double att_range = 0.3;

int main(int argc, char const* argv[]) {
  // local variables to be read in
  int NCELLS, nsmall, seed, ndelete;
  double calA0, kl, ka = 1.0, kb, phiMin, phiMax, att, B = 1.0, t_stress, time_dbl, tauRatio;
  double ka_for_equilibration = 2.0, kb_for_equilibration = 0.0;
  bool boolSmooth, purseStringOn = true;
  double strainRate_ps, k_ps, k_LP, tau_LP, deltaSq, maxProtrusionLength;

  // read in parameters from command line input
  string NCELLS_str = argv[1];
  string nsmall_str = argv[2];
  string ndelete_str = argv[3];
  string calA0_str = argv[4];
  string phiMin_str = argv[5];
  string phiMax_str = argv[6];
  string kl_str = argv[7];
  string ka_str = argv[8];
  string kb_str = argv[9];
  string att_str = argv[10];
  string omega_str = argv[11];
  string deltaSq_str = argv[12];
  string k_ps_str = argv[13];
  string k_lp_str = argv[14];
  string tau_lp_str = argv[15];
  string d_flag_str = argv[16];
  string t_stress_str = argv[17];
  string tauRatio_str = argv[18];
  string smooth_str = argv[19];
  string seed_str = argv[20];
  string time_str = argv[21];
  string outFileStem = argv[22];

  // using sstreams to get parameters
  stringstream NCELLSss(NCELLS_str);
  stringstream nsmallss(nsmall_str);
  stringstream ndeletess(ndelete_str);
  stringstream calA0ss(calA0_str);
  stringstream phiMinss(phiMin_str);
  stringstream phiMaxss(phiMax_str);
  stringstream klss(kl_str);
  stringstream kass(ka_str);
  stringstream kbss(kb_str);
  stringstream attss(att_str);
  stringstream omegass(omega_str);
  stringstream deltaSqss(deltaSq_str);
  stringstream k_psss(k_ps_str);
  stringstream k_lpss(k_lp_str);
  stringstream tau_lpss(tau_lp_str);
  stringstream d_flagss(d_flag_str);
  stringstream t_stressss(t_stress_str);
  stringstream tauRatioss(tauRatio_str);
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
  kbss >> kb;
  attss >> att;
  omegass >> strainRate_ps;
  deltaSqss >> deltaSq;
  k_psss >> k_ps;
  k_lpss >> k_LP;
  tau_lpss >> tau_LP;
  d_flagss >> maxProtrusionLength;
  t_stressss >> t_stress;
  tauRatioss >> tauRatio;
  smoothss >> boolSmooth;
  seedss >> seed;
  timess >> time_dbl;

  const double phi0 = phiMin;
  // auto set deltasq to 0 if pursestring is disabled, so that all ps springs break off
  if (strainRate_ps < 1e-20)
    deltaSq = 0.0;

  epi2D epithelial(NCELLS, 0.0, 0.0, strainRate_ps, k_ps, k_LP, tau_LP, deltaSq, maxProtrusionLength, seed);

  epithelial.openFileStreams(outFileStem);

  // set spring constants
  epithelial.setka(ka_for_equilibration);
  epithelial.setkl(kl);
  epithelial.setkb(kb_for_equilibration);
  epithelial.setkc(kc);
  epithelial.setB(B);
  epithelial.setpbc(0, false);
  epithelial.setpbc(1, false);
  epithelial.setRandPsi();

  // set adhesion scale
  epithelial.setl1(att);
  epithelial.setl2(att_range);
  if (att > att_range || fabs(att - att_range) < 1e-5) {
    cout << "attraction longer than attraction range; discontinuous adhesive potential; error, exiting\n";
    cout << ", or,  adhesion lengthscales are the same, so we will divide by zero; error, exiting \n";
    return 1;
  }

  // set base dpm force, upcast derived epi2D forces
  dpmMemFn repulsiveForceUpdate = &dpm::repulsiveForceUpdate;
  dpmMemFn attractiveForceUpdate = static_cast<void (dpm::*)()>(&epi2D::attractiveForceUpdate_2);
  dpmMemFn crawlingWithPSForceUpdate = static_cast<void (dpm::*)()>(&epi2D::crawlingWithPurseString);
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
  epithelial.setka(ka);
  epithelial.setkb(kb);

  // epithelial.vertexCompress2Target2D_polygon(repulsiveForceUpdateWithCircularWalls, Ftol, dt0, phiMax, dphi0);
  epithelial.shrinkPolyWall(attractiveForceUpdateWithCircularWalls, Ftol, dt0, phiMax, dphi0);
  epithelial.moveSimulationToPositiveCoordinates();  // positive coordinates make the neighbor list storage work better
  epithelial.printConfiguration2D();

  if (t_stress > 9e4) {  //     cout << "t_stress is large, setting to infinity\n";
    t_stress = INFINITY;
  }
  if (tauRatio < 1e-5) {
    tauRatio = INFINITY;
  }

  // after compress, turn on damped NVE
  double T = 1e-2;
  double relaxTime = 50.0;
  double printInterval = relaxTime / 5.0;
  epithelial.drawVelocities2D(T);
  epithelial.setl00();  // set l00 to be l0 before setting relaxation time

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

  // rigidify cell areas before running NVE dynamics
  epithelial.vertexNVE(customForceUpdate_inactive, dt0, 1.0 * relaxTime / dt0, 0);
  // energy minimize without walls to get proper cell shapes
  epithelial.dampedNP0(customForceUpdate_inactive, dt0, 1.0 * relaxTime, 0);

  // save image right before wounding
  epithelial.printConfiguration2D();

  // LASER ABLATION SCHEME
  double xLoc = 0.0, yLoc = 0.0;
  int numCellsToAblate = ndelete;
  epithelial.laserAblate(numCellsToAblate, sizeratio, nsmall, xLoc, yLoc);
  epithelial.zeroMomentum();

  //  dampedNP0 runs simulation with purse-string integration and crawling
  for (int i = 0; i < 2; i++)
    epithelial.dampedNP0(customForceUpdate_inactive, dt0, relaxTime, 0);

  epithelial.setMaxwellRelaxationTime(t_stress);
  epithelial.setTaus(tauRatio * t_stress);
  epithelial.dampedNP0(customForceUpdate_active, dt0, time_dbl, printInterval, purseStringOn);

  cout << "ending simulation: printing one last configuration\n";
  epithelial.printConfiguration2D();
  cout << "\n** Finished laserAblation.cpp, ending. " << endl;
  return 0;
}