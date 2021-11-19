// File to compress DPM to jamming, then conduct a notch test (introduce defect and pull on it with tension)
//
// Cells are bidisperse, with bending energy (broken) and vertex-vertex attraction options.
//
//
// Create bidisperse DPM particles, set constants, place particle centers,
// relax shapes + positions, compress to target phi, delete cells, conduct
// damped MD dynamics.
//
// Compilation command:
// g++ -O3 --std=c++11 -I src main/epi2D/notchTest.cpp src/dpm.cpp src/epi2D.cpp -o main/epi2D/notchTest.o
// ./main/epi2D/notchTest.o 12 24 1.08 0.8 0.805 1.0 0.3 1.0 0 0.5 0.001 1 100 pos.test energy.test stress.test void.test uniaxial
// ./main/epi2D/notchTest.o 24 24 1.08 0.8 0.805 1.0 0.3 1.0 0 0.5 0.01 1 100 pos.test energy.test stress.test void.test uniaxial
// ./main/epi2D/notchTest.o 48 24 1.08 0.8 0.8 1.0 0.3 1.0 0 0.5 0.01 2 10 pos.test energy.test stress.test void.test uniaxial
// ./main/epi2D/notchTest.o 96 24 1.08 0.8 0.805 1.0 0.3 1.0 0 0.5 0.01 1 100 pos.test energy.test stress.test void.test uniaxial
// ./main/epi2D/notchTest.o 192 24 1.08 0.8 0.805 1.0 0.3 1.0 0 0.5 0.001 1 100 pos.test energy.test stress.test void.test uniaxial
// ./main/epi2D/notchTest.o 960 24 1.08 0.8 0.805 1.0 0.3 1.0 0 0.5 0.001 1 100 pos.test energy.test stress.test void.test uniaxial
//
//
// Parameter input list
// 1. NCELLS: 			number of particles
// 2. nsmall: 			number of vertices on small particles (larger
// particles set by 1.4:1.0 bidispersity)
// 3. calA0: 			  preferred shape parameter for all particles
// 4. phiMin: 			p ?
// 5. phiMax: 			p
// 6. kl: 				  perimeter spring constant
// 7. att:          attraction range and strength parameter
// 8. B:            (over)damping coefficient gamma
// 9. Dr0:          rotational diffusion constant for protrusion activity (unused for notchTest code)
// 10. strain: (L-L0)/L0, corresponds to % deformation (0.5 = 50%) due to applied strain
// 11. strainRate:  constant true strain rate = constant fractional deformation per unit time
// 12. seed: 			  seed for random number generator
// 13. time:        amount of time (tau) to simulate
// 14. positionFile:string of path to output file with position/configuration data'
// 15. energyFile:  string of path to output file with energy data
// 16. stressFile:  string of path to output file with stress data
// 17. voidFile:    string of path to output file with void location data
// 18. loadingType: string representing whether to conduct a uniaxial loading or isotropic loading experiment

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
  int NCELLS, nsmall, seed;
  double calA0, kl, kb = 0.0, phiMin, phiMax, att, B, Dr0, time_dbl, strain, strainRate;
  int loadingType;  // 0 = uniaxial (tensile), 1 = isotropic (tensile)

  // read in parameters from command line input
  string NCELLS_str = argv[1];
  string nsmall_str = argv[2];
  string calA0_str = argv[3];
  string phiMin_str = argv[4];
  string phiMax_str = argv[5];
  string kl_str = argv[6];
  string att_str = argv[7];
  string B_str = argv[8];
  string Dr0_str = argv[9];
  string strain_str = argv[10];
  string strainRate_str = argv[11];
  string seed_str = argv[12];
  string time_str = argv[13];
  string positionFile = argv[14];
  string energyFile = argv[15];
  string stressFile = argv[16];
  string voidFile = argv[17];
  string loadingType_str = argv[18];

  // using sstreams to get parameters
  stringstream NCELLSss(NCELLS_str);
  stringstream nsmallss(nsmall_str);
  stringstream calA0ss(calA0_str);
  stringstream phiMinss(phiMin_str);
  stringstream phiMaxss(phiMax_str);
  stringstream klss(kl_str);
  stringstream attss(att_str);
  stringstream Bss(B_str);
  stringstream Dr0ss(Dr0_str);
  stringstream strainss(strain_str);
  stringstream strainRatess(strainRate_str);
  stringstream seedss(seed_str);
  stringstream timess(time_str);

  // read into data
  NCELLSss >> NCELLS;
  nsmallss >> nsmall;
  calA0ss >> calA0;
  phiMinss >> phiMin;
  phiMaxss >> phiMax;
  klss >> kl;
  Dr0ss >> Dr0;
  attss >> att;
  Bss >> B;
  strainss >> strain;
  strainRatess >> strainRate;
  seedss >> seed;
  timess >> time_dbl;
  if (loadingType_str == "uniaxial") {
    loadingType = 0;
  } else if (loadingType_str == "isotropic") {
    loadingType = 1;
  } else {
    cout << "ERROR: loadingType not valid, exiting:\n";
    return 1;
  }

  // number of time steps
  int NT = int(time_dbl / dt0);

  const double phi0 = phiMin;

  // instantiate object
  epi2D epithelial(NCELLS, 0.0, 0.0, Dr0, seed);
  epithelial.setl1(att);
  epithelial.setl2(att_range);
  if (att > att_range) {
    cout << "att, att_range = " << att << ',' << att_range << '\n';
    cout << "attraction stronger than attraction range; discontinuous adhesive potential; error, exiting!\n";
    return 1;
  }

  epithelial.openPosObject(positionFile);
  epithelial.openEnergyObject(energyFile);
  epithelial.openStressObject(stressFile);
  epithelial.openBoundaryObject(voidFile);

  // set spring constants
  epithelial.setka(ka);
  epithelial.setkl(kl);
  epithelial.setkb(kb);
  epithelial.setkc(kc);

  epithelial.monodisperse2D(calA0, nsmall);
  epithelial.initializevnn();

  epithelial.initializePositions2D(phi0, Ftol);

  epithelial.initializeNeighborLinkedList2D(boxLengthScale);

  // set base dpm force, upcast derived epi2D forces
  dpmMemFn repulsiveForceUpdate = &dpm::repulsiveForceUpdate;
  dpmMemFn attractiveForceUpdate = static_cast<void (dpm::*)()>(&epi2D::attractiveForceUpdate_2);
  dpmMemFn activeForceUpdate = static_cast<void (dpm::*)()>(&epi2D::activeAttractiveForceUpdate);

  epithelial.vertexCompress2Target2D(repulsiveForceUpdate, Ftol, dt0, phiMax, dphi0);

  //after compress, turn on damped NVE
  double T = 1e-4;
  epithelial.drawVelocities2D(T);
  epithelial.dampedNVE2D(attractiveForceUpdate, B, dt0, time_dbl / 10, 0);
  epithelial.setInitialLx(epithelial.getL(0));

  //ELASTIC SHEET FRACTURE SCHEME (introduce defect + tensile loading)
  int numCellsToDelete = 0;
  if (NCELLS < 50)
    numCellsToDelete = 2;
  else if (NCELLS < 200)
    numCellsToDelete = 5;
  else
    numCellsToDelete = 10;
  double printInterval = 1.0;  // do not set larger than 10 tau (relax time set in notchTest, both this and that are hardcoded)!

  if (loadingType == 0) {
    epithelial.notchTest(numCellsToDelete, strain, strainRate, boxLengthScale, sizeratio,
                         nsmall, attractiveForceUpdate, B, dt0, printInterval, "uniaxial");
  } else {
    cout << "loadingType not found. Closing.\n";
    return 1;
  }

  cout << "\n** Finished notchTest.cpp, ending. " << endl;

  return 0;
}