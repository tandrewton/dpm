// testing smooth DPM adhesion force
// Compilation command:
// g++ -O3 --std=c++11 -g -I src main/cell/testAdhesionForces.cpp src/dpm.cpp src/cell.cpp -o main/cell/testAdhesion.o
// run command:
// ./main/cell/testAdhesion.o   2   10 1.0 0.3  1    test
//                          NCELLS NV  A0  att seed outFileStem

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
const double phi0 = 0.85;           // initial packing fraction
const double phiMax = 0.75;
const double smallfrac = 1.0;  // fraction of small particles
const double sizeratio = 1.0;  // size ratio between small and large particles
const double dt0 = 1e-2;       // initial magnitude of time step in units of MD time
const double Ptol = 1e-8;
const double Ftol = 1e-12;
const double att_range = 0.5;  // att_range > 0 might cause problems with current scheme for smooth forces

int main(int argc, char const* argv[]) {
  // local variables to be read in
  int NCELLS, nv, seed;
  double calA0, att;

  // read in parameters from command line input
  string NCELLS_str = argv[1];
  string nv_str = argv[2];
  string calA0_str = argv[3];
  string att_str = argv[4];
  string seed_str = argv[5];
  string outFileStem = argv[6];

  string positionFile = outFileStem + ".pos";

  // using sstreams to get parameters
  stringstream NCELLSss(NCELLS_str);
  stringstream nvss(nv_str);
  stringstream calA0ss(calA0_str);
  stringstream attss(att_str);
  stringstream seedss(seed_str);

  // read into data
  NCELLSss >> NCELLS;
  nvss >> nv;
  calA0ss >> calA0;
  attss >> att;
  seedss >> seed;

  int numCellTypes = 1;
  cell cell2D(NCELLS, 0.0, 0.0, seed, numCellTypes);
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

  dpmMemFn repulsiveForceUpdate = &dpm::repulsiveForceUpdate;
  dpmMemFn repulsiveForceUpdateWithWalls = static_cast<void (dpm::*)()>(&cell::repulsiveForceUpdateWithWalls);
  dpmMemFn attractiveForceUpdate = static_cast<void (dpm::*)()>(&cell::attractiveForceUpdate);
  dpmMemFn repulsivePolarityForceUpdate = static_cast<void (dpm::*)()>(&cell::repulsiveWithPolarityForceUpdate);
  dpmMemFn attractivePolarityForceUpdate = static_cast<void (dpm::*)()>(&cell::attractiveWithPolarityForceUpdate);
  dpmMemFn repulsiveForceUpdateWithPolyWalls = static_cast<void (dpm::*)()>(&cell::repulsiveForceUpdateWithPolyWall);
  dpmMemFn attractiveForceUpdateWithPolyWalls = static_cast<void (dpm::*)()>(&cell::attractiveForceUpdateWithPolyWall);
  
  // initialize particles with the same number of vertices and the same preferred shape parameter calA0
  cell2D.monodisperse2D(calA0, nv);
  // initialize particle positions
  cell2D.initializePositions2D(phi0, Ftol);
  // cell2D.printConfiguration2D();

  cell2D.initializeNeighborLinkedList2D(boxLengthScale);

  // compress to target packing fraction
  //cell2D.vertexCompress2Target2D(repulsiveForceUpdate, Ftol, dt0, phiMax, dphi0);
  //cout << "done compressing to target packing fraction\n";
  //cell2D.printConfiguration2D();

  bool wallsBool = true;
  double relaxTimeShort = 10.0;
  double relaxTime = 100.0;
  double runTime = 3000.0;
  double B = 1.0;
  std::vector<double> savedPositions;

  // dpmMemFn customForceUpdate = repulsivePolarityForceUpdate;
  // dpmMemFn customForceUpdate = attractivePolarityForceUpdate;
  dpmMemFn customForceUpdate = attractiveForceUpdate;
  //cell2D.dampedVertexNVE(customForceUpdate, B, dt0, relaxTimeShort, relaxTimeShort/5);

  // expecting exactly 2 cells with 10 vertices each here
  double diameter = 2*cell2D.getr(0);
  double offset = - diameter * 0.25;

  for (int gi = 0; gi < 20; gi++){  
    if (gi < 10)
      cell2D.moveVertex(gi, 0, 0); // move most cell 0 vertices out of the way
  }

  // dpm is written counterclockwise for vertex numbering
  bool testConcave = false;
  bool test90Degrees = false;
  bool testAttraction = true;
  double psi, theta, dx, dy, dx2, dy2;
  int numDirections = 4; // numDirections sets how many directions we're scanning over when testing overlaps
  if (!testConcave){
    // move vertices in cell 1 into a rectangle for testing forces on vertex 0 cell 0 as a function of distance
    cell2D.moveVertex(19, 1+0*diameter, 1+0*diameter);
    cell2D.moveVertex(18, 1+0*diameter, 1+1*diameter);
    cell2D.moveVertex(17, 1+0*diameter, 1+2*diameter);
    cell2D.moveVertex(16, 1+0*diameter, 1+3*diameter);
    cell2D.moveVertex(15, 1+1*diameter, 1+3*diameter);
    cell2D.moveVertex(14, 1+2*diameter, 1+3*diameter);
    cell2D.moveVertex(13, 1+2*diameter, 1+2*diameter);
    cell2D.moveVertex(12, 1+2*diameter, 1+1*diameter);
    cell2D.moveVertex(11, 1+2*diameter, 1+0*diameter);
    cell2D.moveVertex(10, 1+1*diameter, 1+0*diameter);
  } else if (testConcave && test90Degrees) {
    // move vertices in cell 1 into a concave shape for concave testing
    cell2D.moveVertex(19, 1+1*diameter, 1+1*diameter);
    cell2D.moveVertex(18, 1+1*diameter, 1+2*diameter);
    cell2D.moveVertex(17, 1+0*diameter, 1+2*diameter);
    cell2D.moveVertex(16, 1+0*diameter, 1+3*diameter);
    cell2D.moveVertex(15, 1+1*diameter, 1+3*diameter);
    cell2D.moveVertex(14, 1+2*diameter, 1+3*diameter);
    cell2D.moveVertex(13, 1+2*diameter, 1+2*diameter);
    cell2D.moveVertex(12, 1+2*diameter, 1+1*diameter);
    cell2D.moveVertex(11, 1+2*diameter, 1+0*diameter);
    cell2D.moveVertex(10, 1+1*diameter, 1+0*diameter);
  } else if (testConcave && !test90Degrees) {
    numDirections = 2;
    psi = 3*PI/4;
    theta = PI/4;
    dx = diameter*cos(theta), dy = diameter*sin(theta);
    dx2 = diameter*cos(5*PI/4 - psi), dy2 = diameter*sin(5*PI/4 - psi);
    // move vertices in cell 1 into a concave shape for concave testing
    cell2D.moveVertex(19, 1+0*diameter, 1+0*diameter);
    cell2D.moveVertex(18, 1+0*diameter, 1+1*diameter);
    cell2D.moveVertex(17, 1+0*diameter, 1+2*diameter);
    cell2D.moveVertex(16, 1+0*diameter+dx, 1+2*diameter+dy);
    cell2D.moveVertex(15, 1+0*diameter+dx+dx2, 1+2*diameter+dy+dy2);
    cell2D.moveVertex(14, 1+2*diameter, 1+3*diameter);
    cell2D.moveVertex(13, 1+2*diameter, 1+2*diameter);
    cell2D.moveVertex(12, 1+2*diameter, 1+1*diameter);
    cell2D.moveVertex(11, 1+2*diameter, 1+0*diameter);
    cell2D.moveVertex(10, 1+1*diameter, 1+0*diameter);
  }

  /*// test case : move 17 16 15 out of the way to see if the presence of 17 affects forces
  cell2D.moveVertex(17, 1+1*diameter, 1+3*diameter);
  cell2D.moveVertex(16, 1+1*diameter, 1+4*diameter);
  cell2D.moveVertex(15, 1+2*diameter, 1+4*diameter);*/

  if (testAttraction){
   offset = -offset/2;
  }
  double origin = 1-diameter-offset;

  cell2D.moveVertex(0, origin, 1); // place vertex 0 from cell 0 next to vertex 0-2 in cell 1.

  //cell2D.printConfiguration2D();
  cell2D.saveConfiguration(savedPositions);

  std::vector<double> forcex, forcey, energy;

  for (int i = 0; i < numDirections; i++){
    cout << "att, att_range = " << att << '\t' << att_range << '\n';
    // revert to original position, then change attraction
    cell2D.loadConfiguration(savedPositions);
    cell2D.setl1(att);
    cell2D.setl2(att_range);
    
    int epsilonInv = 50;
    if (!testConcave){
      for (int j = 0; j < epsilonInv + 1; j++){
        double fx, fy, u;
        // move vertex 0 in x or y direction along the rectangle made by cell 1
        if (i == 0){
          cell2D.moveVertex(0, origin, origin + j*(diameter*5.0 + 2*offset)/epsilonInv);
        }
        if (i == 1){
          cell2D.moveVertex(0, origin + j*(diameter*4.0 + 2*offset) /epsilonInv, origin + diameter*5.0 + 2*offset);
        }        
        if (i == 2){
          cell2D.moveVertex(0, origin + diameter*4.0 + 2*offset, origin + diameter*5.0 + 2*offset - j*(diameter*5.0 + 2*offset)/epsilonInv);
        }
        if (i == 3){
          cell2D.moveVertex(0, origin + diameter*4.0 + 2*offset - j*(diameter*4.0 + 2*offset)/epsilonInv, origin);
        }
        cout <<  "in frame = " << i*epsilonInv + j + 1 << "\n\n\n";
        // print configuration and calculate its forces
        cell2D.printConfiguration2D();
        //cell2D.attractiveForceUpdatePrint(fx, fy, u);
        cell2D.attractiveForceUpdateSmoothPrint(fx, fy, u);
        if (fy == 0)
          cout << "frame " << i*epsilonInv + j + 1 << " has fy = 0!\n\n";
        forcex.push_back(fx);
        forcey.push_back(fy);
        energy.push_back(u);
      }
    } else if (testConcave && test90Degrees){
      for (int j = 0; j < epsilonInv + 1; j++){
        double fx, fy, u;
        // move vertex 0 in x or y direction along the rectangle made by cell 1
        if (i == 0){
          cell2D.moveVertex(0, origin+ 1*diameter, origin + j*(diameter*2.0)/epsilonInv);
        }
        if (i == 1){
          cell2D.moveVertex(0, origin+ 1*diameter - j*(diameter) /epsilonInv, origin + diameter*2.0);
        }
        if (i == 2){
          cell2D.moveVertex(0, origin, origin + diameter*2.0 + j*(diameter*3.0 + 2*offset)/epsilonInv);
        }
        if (i == 3){
          cell2D.moveVertex(0, origin + j*(diameter*4.0 + 2*offset)/epsilonInv, origin + 2*diameter + 3*diameter + 2*offset);
        }
        cout <<  "in frame = " << i*epsilonInv + j + 1 << "\n\n\n";
        // print configuration and calculate its forces
        cell2D.printConfiguration2D();
        //cell2D.attractiveForceUpdatePrint(fx, fy, u);
        cell2D.attractiveForceUpdateSmoothPrint(fx, fy, u);
        if (fy == 0)
          cout << "frame " << i*epsilonInv + j + 1 << " has fy = 0!\n\n";
        forcex.push_back(fx);
        forcey.push_back(fy);
        energy.push_back(u);
      }
    } else if (testConcave) {
      for (int j = 0; j < epsilonInv + 1; j++){
        double fx, fy, u;
        // move vertex 0 in x or y direction along the rectangle made by cell 1
        if (i == 0){
          cell2D.moveVertex(0, 1 + 4.0*offset + j*1.0/epsilonInv*dx, 1 + 2*diameter + 0*offset + j*1.0/epsilonInv*dy);
        }
        if (i == 1){
          cell2D.moveVertex(0, 1 + 4.0*offset + dx + j*1.0/epsilonInv*dx2, 1 + 2*diameter + 0*offset + dy + j*1.0/epsilonInv*dy2);
        }
        if (i == 2){
          cell2D.moveVertex(0, origin + diameter*4.0 + 2*offset, origin + diameter*5.0 + 2*offset - j*(diameter*5.0 + 2*offset)/epsilonInv);
        }
        if (i == 3){
          cell2D.moveVertex(0, origin + diameter*4.0 + 2*offset - j*(diameter*4.0 + 2*offset)/epsilonInv, origin);
        }
        cout <<  "in frame = " << i*epsilonInv + j + 1 << "\n\n\n";
        // print configuration and calculate its forces
        cell2D.printConfiguration2D();
        //cell2D.attractiveForceUpdatePrint(fx, fy, u);
        cell2D.attractiveForceUpdateSmoothPrint(fx, fy, u);
        if (fy == 0)
          cout << "frame " << i*epsilonInv + j + 1 << " has fy = 0!\n\n";
        forcex.push_back(fx);
        forcey.push_back(fy);
        energy.push_back(u);
      }
    }
  }

  std::ofstream energyAndForce("energyAndForceTest.txt");
  for (int i = 0; i < forcex.size(); i++){
    energyAndForce << forcex[i] << '\t' << forcey[i] << '\t' << energy[i] << '\n';
  }

  cout
      << "\n** Finished testAdhesion.cpp, ending. " << endl;
      
  return 0;
}

