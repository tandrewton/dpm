// File to of just one DPM for figure
//
//
// Compilation command:
// g++ -O3 --std=c++11 -g -I src main/epi2D/oneDP.cpp src/dpm.cpp src/epi2D.cpp -o main/epi2D/oneDP.o
// ./main/epi2D/oneDP.o 1.0 10.0 oneDP.pos oneDP.energy oneDP.stress oneDP.void
//   ^ run this with dampedCompression to get the illustration of compression and decompression under plastic/elastic
//
//                      ka  tau_maxwell
//
// ./main/epi2D/oneDP.o 1.0 10000 oneDP.pos oneDP.energy oneDP.stress oneDP.void
/*
ka_arr=(0.01 0.05 0.1 0.2 0.4 1.6 12.8)
tau_arr=(10.0 20.0 40.0 80.0 320.0)
for ka in ${ka_arr[@]}; do
    for t_stress in ${tau_arr[@]}; do
        echo ./main/epi2D/oneDP.o $ka $t_stress oneDP.pos oneDP.energy oneDP.stress oneDP.void
    done
done
*/
//
// positionFile: string of path to output file with position/configuration data

// header files
#include <sstream>
#include "epi2D.h"

// preprocessor macros
#define NDIM 2

// namspace
using namespace std;

const bool plotCompression = 0;     // whether or not to plot configuration during compression protocol
const double dphi0 = 0.005;         // packing fraction increment
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
  int NCELLS = 1, nsmall = 24, seed = 1, gi;
  double calA0 = 1.0, ka, kl = 1.0, kb = 0.01, phiMin = 0.6, phiMax = 0.4, att = 0.1;
  double t_plastic, B = 1.0;
  bool boolCIL;

  // read in parameters from command line input
  string ka_str = argv[1];
  string tau_str = argv[2];
  string positionFile = argv[3];
  string energyFile = argv[4];
  string stressFile = argv[5];
  string voidFile = argv[6];
  string edgeFile = "edge.test";

  stringstream kass(ka_str);
  stringstream tauss(tau_str);

  kass >> ka;
  tauss >> t_plastic;
  if (t_plastic > 1e6)
    t_plastic = INFINITY;

  const double phi0 = phiMin;

  // instantiate object
  epi2D epithelial(NCELLS, 0.0, 0.0, 0.01, 4.0, 4.0, 1.0, 0.0, 0.0, seed);

  epithelial.openPosObject(positionFile);
  epithelial.openEnergyObject(energyFile);
  epithelial.openStressObject(stressFile);
  epithelial.openBoundaryObject(voidFile);
  epithelial.openEdgeObject(edgeFile);

  // set spring constants
  epithelial.setka(ka);
  epithelial.setkl(kl);
  epithelial.setkb(kb);
  epithelial.setkc(kc);
  epithelial.setB(B);

  // set activity scale, and CIL option
  epithelial.setpbc(0, false);
  epithelial.setpbc(1, false);

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

  epithelial.initializeNeighborLinkedList2D(boxLengthScale);

  // set base dpm force, upcast derived epi2D forces
  dpmMemFn repulsiveForceUpdate = &dpm::repulsiveForceUpdate;
  dpmMemFn repulsiveForceUpdateWithWalls = static_cast<void (dpm::*)()>(&epi2D::repulsiveForceUpdateWithWalls);
  dpmMemFn attractiveForceUpdate = static_cast<void (dpm::*)()>(&epi2D::attractiveForceUpdate_2);

  dpmMemFn circuloLineAttraction = static_cast<void (dpm::*)()>(&epi2D::attractiveForceUpdate_circulo);
  dpmMemFn customForceUpdate_inactive = circuloLineAttraction;

  epithelial.vertexCompress2Target2D(repulsiveForceUpdateWithWalls, Ftol, dt0, phiMax, dphi0);

  epithelial.setMaxwellRelaxationTime(t_plastic);

  // epithelial.dampedNVE2D(repulsiveForceUpdate, dt0, 500, 10);
  for (int i = 0; i < 3; i++) {
    epithelial.displaceVertex(0, i, 0.005, 0.0);
  }
  epithelial.printConfiguration2D();

  epithelial.dampedCompression(repulsiveForceUpdate, dt0, 500, 25);  // used for my compression experiments for wound healing explanation of plasticity

  std::string forceDipoleFilename = "forceDipole/forceDipole_ka_" + ka_str + "_tau_" + tau_str + ".txt";
  double forceMoment = 0.01;
  // epithelial.dampedForceDipoleExperiment(repulsiveForceUpdate, 0.01, dt0, 500, 25, forceDipoleFilename);  // used for my force dipole experiments to explain how hard plasticity deforms more than elastic or soft plastic

  // epithelial.drawVelocities2D(1e-3);
  //   epithelial.vertexNVE(repulsiveForceUpdate, dt0, 20000, 250);
  //  epithelial.vertexNVE(customForceUpdate_inactive, dt0, 100000, 2000);
  // epithelial.dampedNVE2D(repulsiveForceUpdate, dt0, 500, 100);

  epithelial.printConfiguration2D();
  cout << "\n** Finished oneDP.cpp, ending. " << endl;
  return 0;
}