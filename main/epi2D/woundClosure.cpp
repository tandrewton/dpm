// File to simulate wound healing starting from a set of initial positions
//
// Compilation command:
// g++ -O3 --std=c++11 -g -I src main/epi2D/woundClosure.cpp src/dpm.cpp src/epi2D.cpp -o main/epi2D/woundClosure.o

// below: no purse-string, only crawling
//./main/epi2D/woundClosure.o 100 20 4 1.10 0.01 0.925 1.0 0.2 0.01  0.0  1.0  2.0 1.0  3.0  1.0 0.5  0  0.00   1  200  test
// ........................... N  NV Nd A0  pMin pMax kl att  om   dsq  kps  klp tau dflag  B  Dr0 CIL prate  sd time file
// below: purse-string, no crawling
//./main/epi2D/woundClosure.o 100 20 4 1.10 0.01 0.925 1.0 0.2 0.01  2.0  1.0  2.0 1.0  0.0  1.0 0.5  0  0.00   1  200  test
// ........................... N  NV Nd A0  pMin pMax kl att  om   dsq  kps  klp tau dflag  B  Dr0 CIL prate  sd time file
// below: purse-string, and crawling
// ./main/epi2D/woundClosure.o 100 20 4 1.10 0.01 0.925 1.0 0.2 0.01  2.0  1.0  2.0 1.0  3.0  1.0 0.5  0  0.00   1  400  test
// ........................... N  NV Nd A0  pMin pMax  kl att  om   dsq  kps  klp tau dflag  B  Dr0 CIL prate  sd time file
// below: no purse-string, and no crawling
//./main/epi2D/woundClosure.o 100 20 4 1.10 0.01 0.925 1.0 0.2 0.01  0.0  1.0  2.0 1.0  0.0  1.0 0.5  0  0.00   1  200  test
// ........................... N  NV Nd A0  pMin pMax  kl att  om   dsq  kps  klp tau dflag  B  Dr0 CIL prate  sd time file

// make sure the corresponding bash file writes to a different output and pipeline folder, instead of ablate it'll be
//  woundClosure or something

//
// Parameter input list
// 1. NCELLS: 			number of particles
// 2. nv: 			number of vertices on particles
// 3. ndelete:      number of cells to delete
// 4. calA0: 			  preferred shape parameter for all particles
// 5. phiMin: 			phi (total particle area/box area)
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
// 18. shapeRelaxRate rate for how quickly cells relax their perimeters
// 19. seed: 			  seed for random number generator
// 20. time:        amount of time (tau) to simulate
// 21. outFileStem  stem of output file names, i.e. for "test", energy.test, position.test, etc.

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
  int NCELLS, nv, seed, gi, ndelete;
  double calA0, kl, kb = 0.0, phiMin, phiMax, att, B, Dr0, time_dbl;
  bool boolCIL;
  double strainRate_ps, k_ps, k_LP, tau_LP, deltaSq, maxProtrusionLength, shapeRelaxationRate;

  // read in parameters from command line input
  string NCELLS_str = argv[1];
  string nv_str = argv[2];
  string ndelete_str = argv[3];
  string calA0_str = argv[4];
  string phiMin_str = argv[5];
  string phiMax_str = argv[6];
  string kl_str = argv[7];
  string att_str = argv[8];
  string omega_str = argv[9];
  string deltaSq_str = argv[10];
  string k_ps_str = argv[11];
  string k_lp_str = argv[12];
  string tau_lp_str = argv[13];
  string d_flag_str = argv[14];
  string B_str = argv[15];
  string Dr0_str = argv[16];
  string boolCIL_str = argv[17];
  string shapeRelax_str = argv[18];
  string seed_str = argv[19];
  string time_str = argv[20];
  string outFileStem = argv[21];

  string positionFile = outFileStem + ".pos";
  string energyFile = outFileStem + ".energy";
  string stressFile = outFileStem + ".stress";
  string voidFile = outFileStem + ".void";
  string edgeFile = outFileStem + ".edge";
  string purseStringFile = outFileStem + ".purseString";
  string voidAreaFile = outFileStem + ".voidArea";

  // using sstreams to get parameters
  stringstream NCELLSss(NCELLS_str);
  stringstream nvss(nv_str);
  stringstream ndeletess(ndelete_str);
  stringstream calA0ss(calA0_str);
  stringstream phiMinss(phiMin_str);
  stringstream phiMaxss(phiMax_str);
  stringstream klss(kl_str);
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
  stringstream shapeRelaxss(shapeRelax_str);
  stringstream seedss(seed_str);
  stringstream timess(time_str);

  // read into data
  NCELLSss >> NCELLS;
  nvss >> nv;
  ndeletess >> ndelete;
  calA0ss >> calA0;
  phiMinss >> phiMin;
  phiMaxss >> phiMax;
  klss >> kl;
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
  shapeRelaxss >> shapeRelaxationRate;
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
  // instantiate object

  // read in file
  // cin int double double
  // reset NCELLS

  std::string inputConfigFile = "testDPMConfigurationFile_Verboon_wounded.txt";

  // need to read in file and then construct epi2D object so I have initialized areas, or just use more cells so that I have enough memory allocated, and hopefully all the resizing operations appropriately manage that (is area resized?)

  epi2D epithelial(NCELLS, 0.0, 0.0, Dr0, strainRate_ps, k_ps, k_LP, tau_LP, deltaSq, maxProtrusionLength, seed);

  // read in configuration file (cells, vertices, and coordinates), establish shape parameters
  epithelial.initializeFromConfigurationFile(inputConfigFile, phi0);

  cout << "before initializevnn\n";
  epithelial.initializevnn();
  cout << "after initializevnn\n";
  // epithelial.initializePositions2D(phi0, Ftol, false);
  epithelial.initializeFieldStress();

  epithelial.openPosObject(positionFile);
  epithelial.openEnergyObject(energyFile);
  epithelial.openStressObject(stressFile);
  epithelial.openBoundaryObject(voidFile);
  epithelial.openEdgeObject(edgeFile);
  epithelial.openPurseStringObject(purseStringFile);
  epithelial.openVoidAreaObject(voidAreaFile);

  // set spring constants
  epithelial.setka(ka);
  epithelial.setkl(kl);
  epithelial.setkb(kb);
  epithelial.setkc(kc);
  epithelial.setShapeRelaxationRate(shapeRelaxationRate);

  // set CIL option
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
  if (fabs(att - att_range) < 1e-5) {
    cout << "or, adhesion lengthscales are the same, so we will divide by zero; error' exiting!\n";
    return 1;
  }

  epithelial.printConfiguration2D();

  epithelial.initializeNeighborLinkedList2D(boxLengthScale);

  // set base dpm force, upcast derived epi2D forces
  dpmMemFn repulsiveForceUpdate = &dpm::repulsiveForceUpdate;
  dpmMemFn repulsiveForceUpdateWithWalls = static_cast<void (dpm::*)()>(&epi2D::repulsiveForceUpdateWithWalls);
  dpmMemFn attractiveForceUpdate = static_cast<void (dpm::*)()>(&epi2D::attractiveForceUpdate_2);
  dpmMemFn substrateAdhesionForceUpdate = static_cast<void (dpm::*)()>(&epi2D::substrateadhesionAttractiveForceUpdate);
  dpmMemFn repulsiveForceUpdateWithCircularAperture = static_cast<void (dpm::*)()>(&epi2D::repulsiveForceWithCircularApertureWall);

  // after initialization, turn on damped NVE
  double T = 1e-4;
  double relaxTime = 10.0;
  epithelial.drawVelocities2D(T);
  // epithelial.dampedNP0(attractiveForceUpdate, B, dt0, relaxTime, 0, wallsOff);

  cout << "after damped NP0\n";
  epithelial.zeroMomentum();

  epithelial.dampedNVE2D(attractiveForceUpdate, B, dt0, relaxTime, 0);

  epithelial.dampedNVE2D(attractiveForceUpdate, B, dt0, time_dbl, 0);
  cout << "after dampedNVE2D\n";

  //  dampedNP0 already takes care of purse-string. might want to separate, or just change spring constant
  epithelial.dampedNP0(substrateAdhesionForceUpdate, B, dt0, time_dbl, time_dbl / 100.0, wallsOff);
  // epithelial.dampedNP0(attractiveForceUpdate, B, dt0, time_dbl, time_dbl / 100.0, wallsOff);

  cout << "\n** Finished woundClosure.cpp, ending. " << endl;
  return 0;
}