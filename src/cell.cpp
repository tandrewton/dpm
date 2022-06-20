#include "cell.h"

// namespace
using namespace Eigen;
using namespace std;

void cell::repulsiveForceUpdateWithWalls() {
  double a, b, c, d;
  resetForcesAndEnergy();
  shapeForces2D();
  vertexRepulsiveForces2D();
  wallForces(false, false, false, false, a, b, c, d);
}

/*void cell::repulsiveForceUpdateWithoutTopWall() {
  double a, b, c, d;
  resetForcesAndEnergy();
  shapeForces2D();
  vertexRepulsiveForces2D();
  wallForces(true, false, false, false, a, b, c, d);
}*/

void cell::repulsiveWithPolarityForceUpdate() {
  resetForcesAndEnergy();
  shapeForces2D();
  // give each cell a y-penalized polarity
  for (int ci = 0; ci < NCELLS; ci++) {
    cellPolarityForces(ci, 0.1, "y");
  }
  vertexRepulsiveForces2D();
}

void cell::attractiveWithPolarityForceUpdate() {
  resetForcesAndEnergy();
  shapeForces2D();
  // give each cell a y-penalized polarity
  for (int ci = 0; ci < NCELLS; ci++) {
    cellPolarityForces(ci, 0.1, "y");
  }
  vertexAttractiveForces2D_2();
}

void cell::vertexAttractiveForces2D_2() {
  // altered from dpm attractive force code, because it works with larger l2
  // values. (warning: probably won't work with bending.) local variables
  int ci, cj, gi, gj, vi, vj, bi, bj, pi, pj;
  double sij, rij, dx, dy, rho0;
  double ftmp, fx, fy;

  // attraction shell parameters
  double shellij, cutij, xij, kint = (kc * l1) / (l2 - l1);
  // cout << "kc / kint = " << kc / kint << '\t' << kc << '\t' << kint << '\n';

  // sort particles
  sortNeighborLinkedList2D();

  // get fundamental length
  rho0 = sqrt(a0[0]);

  // reset contact network
  fill(cij.begin(), cij.end(), 0);

  // loop over boxes in neighbor linked list
  for (bi = 0; bi < NBX; bi++) {
    // get start of list of vertices
    pi = head[bi];

    // loop over linked list
    while (pi > 0) {
      // real particle index
      gi = pi - 1;

      // cell index of gi
      cindices(ci, vi, gi);

      // next particle in list
      pj = list[pi];

      // loop down neighbors of pi in same cell
      while (pj > 0) {
        // real index of pj
        gj = pj - 1;

        // cell index of j
        cindices(cj, vj, gj);

        if (gj == ip1[gi] || gj == im1[gi]) {
          pj = list[pj];
          continue;
        }

        // contact distance
        sij = r[gi] + r[gj];

        // attraction distances
        shellij = (1.0 + l2) * sij;
        cutij = (1.0 + l1) * sij;

        // particle distance
        dx = x[NDIM * gj] - x[NDIM * gi];
        if (pbc[0])
          dx -= L[0] * round(dx / L[0]);
        if (dx < shellij) {
          dy = x[NDIM * gj + 1] - x[NDIM * gi + 1];
          if (pbc[1])
            dy -= L[1] * round(dy / L[1]);
          if (dy < shellij) {
            rij = sqrt(dx * dx + dy * dy);
            if (rij < shellij) {
              // scaled distance
              xij = rij / sij;

              // pick force based on vertex-vertex distance
              if (ci == cj) {
                // if vertices (not neighbors) are in same cell, compute
                // repulsions
                if (rij < sij) {
                  ftmp = kc * (1 - (rij / sij)) * (rho0 / sij);
                  cellU[ci] += 0.5 * kc * pow((1 - (rij / sij)), 2.0);
                } else
                  ftmp = 0;
              } else if (rij > cutij) {
                // force scale
                ftmp = kint * (xij - 1.0 - l2) / sij;

                // increase potential energy
                U += -0.5 * kint * pow(1.0 + l2 - xij, 2.0);
                cellU[ci] += -0.5 * kint * pow(1.0 + l2 - xij, 2.0) / 2.0;
                cellU[cj] += -0.5 * kint * pow(1.0 + l2 - xij, 2.0) / 2.0;
              } else {
                // force scale
                ftmp = kc * (1 - xij) / sij;

                // increase potential energy
                U += 0.5 * kc * (pow(1.0 - xij, 2.0) - l1 * l2);
                cellU[ci] += 0.5 * kc * (pow(1.0 - xij, 2.0) - l1 * l2) / 2.0;
                cellU[cj] += 0.5 * kc * (pow(1.0 - xij, 2.0) - l1 * l2) / 2.0;
              }

              // force elements
              fx = ftmp * (dx / rij);
              fy = ftmp * (dy / rij);

              // add to forces
              F[NDIM * gi] -= fx;
              F[NDIM * gi + 1] -= fy;

              F[NDIM * gj] += fx;
              F[NDIM * gj + 1] += fy;

              // add to virial stress
              // note: 4/7/22 I'm using -dx/2 instead of dx and same for dy for stress calculation, since
              //  I want to calculate force times separation from geometric center of interaction
              stress[0] += -dx * fx;
              stress[1] += -dy * fy;
              stress[2] += -0.5 * (dx * fy + dy * fx);

              fieldStress[gi][0] += -dx / 2 * fx;
              fieldStress[gi][1] += -dy / 2 * fy;
              fieldStress[gi][2] += -0.5 * (dx / 2 * fy + dy / 2 * fx);

              // stress on gj should be the same as on gi, since it's the opposite separation and also opposite force
              fieldStress[gj][0] += -dx / 2 * fx;
              fieldStress[gj][1] += -dy / 2 * fy;
              fieldStress[gj][2] += -0.5 * (dx / 2 * fy + dy / 2 * fx);

              if (ci > cj)
                cij[NCELLS * cj + ci - (cj + 1) * (cj + 2) / 2]++;
              else if (ci < cj)
                cij[NCELLS * ci + cj - (ci + 1) * (ci + 2) / 2]++;
            }
          }
        }

        // update pj
        pj = list[pj];
      }
      // test overlaps with forward neighboring cells
      for (bj = 0; bj < NNN; bj++) {
        // only check if boundaries permit
        if (nn[bi][bj] == -1)
          continue;

        // get first particle in neighboring cell
        pj = head[nn[bi][bj]];

        // loop down neighbors of pi in same cell
        while (pj > 0) {
          // real index of pj
          gj = pj - 1;

          // cell index of j
          cindices(cj, vj, gj);

          if (gj == ip1[gi] || gj == im1[gi]) {
            pj = list[pj];
            continue;
          }

          // contact distance
          sij = r[gi] + r[gj];

          // attraction distances
          shellij = (1.0 + l2) * sij;
          cutij = (1.0 + l1) * sij;

          // particle distance
          dx = x[NDIM * gj] - x[NDIM * gi];
          if (pbc[0])
            dx -= L[0] * round(dx / L[0]);
          if (dx < shellij) {
            dy = x[NDIM * gj + 1] - x[NDIM * gi + 1];
            if (pbc[1])
              dy -= L[1] * round(dy / L[1]);
            if (dy < shellij) {
              rij = sqrt(dx * dx + dy * dy);
              if (rij < shellij) {
                // scaled distance
                xij = rij / sij;

                // pick force based on vertex-vertex distance
                if (ci == cj) {
                  // if vertices (not neighbors) are in same cell, compute
                  // repulsions
                  if (rij < sij) {
                    ftmp = kc * (1 - (rij / sij)) * (rho0 / sij);
                    cellU[ci] += 0.5 * kc * pow((1 - (rij / sij)), 2.0);
                  } else
                    ftmp = 0;
                } else if (rij > cutij) {
                  // force scale
                  ftmp = kint * (xij - 1.0 - l2) / sij;

                  // increase potential energy
                  U += -0.5 * kint * pow(1.0 + l2 - xij, 2.0);
                  cellU[ci] += -0.5 * kint * pow(1.0 + l2 - xij, 2.0) / 2.0;
                  cellU[cj] += -0.5 * kint * pow(1.0 + l2 - xij, 2.0) / 2.0;
                } else {
                  // force scale
                  ftmp = kc * (1 - xij) / sij;

                  // increase potential energy
                  U += 0.5 * kc * (pow(1.0 - xij, 2.0) - l1 * l2);
                  cellU[ci] += 0.5 * kc * (pow(1.0 - xij, 2.0) - l1 * l2) / 2.0;
                  cellU[cj] += 0.5 * kc * (pow(1.0 - xij, 2.0) - l1 * l2) / 2.0;
                }

                // force elements
                fx = ftmp * (dx / rij);
                fy = ftmp * (dy / rij);

                // add to forces
                F[NDIM * gi] -= fx;
                F[NDIM * gi + 1] -= fy;

                F[NDIM * gj] += fx;
                F[NDIM * gj + 1] += fy;

                // add to virial stress
                stress[0] += -dx * fx;
                stress[1] += -dy * fy;
                stress[2] += -0.5 * (dx * fy + dy * fx);

                fieldStress[gi][0] += -dx / 2 * fx;
                fieldStress[gi][1] += -dy / 2 * fy;
                fieldStress[gi][2] += -0.5 * (dx / 2 * fy + dy / 2 * fx);

                // stress on gj should be the same as on gi, since it's the opposite separation and also opposite force
                fieldStress[gj][0] += -dx / 2 * fx;
                fieldStress[gj][1] += -dy / 2 * fy;
                fieldStress[gj][2] += -0.5 * (dx / 2 * fy + dy / 2 * fx);

                if (ci != cj) {
                  if (ci > cj)
                    cij[NCELLS * cj + ci - (cj + 1) * (cj + 2) / 2]++;
                  else if (ci < cj)
                    cij[NCELLS * ci + cj - (ci + 1) * (ci + 2) / 2]++;
                }
              }
            }
          }

          // update pj
          pj = list[pj];
        }
      }

      // update pi index to be next
      pi = list[pi];
    }
  }
  // normalize stress by box area, make dimensionless
  stress[0] *= (rho0 / (L[0] * L[1]));
  stress[1] *= (rho0 / (L[0] * L[1]));
  stress[2] *= (rho0 / (L[0] * L[1]));

  // normalize per-cell stress by preferred cell area
  for (int ci = 0; ci < NCELLS; ci++) {
    for (int vi = 0; vi < nv[ci]; vi++) {
      int gi = gindex(ci, vi);
      fieldStressCells[ci][0] += fieldStress[gi][0];
      fieldStressCells[ci][1] += fieldStress[gi][1];
      fieldStressCells[ci][2] += fieldStress[gi][2];
    }
    fieldStressCells[ci][0] *= rho0 / a0[ci];
    fieldStressCells[ci][1] *= rho0 / a0[ci];
    fieldStressCells[ci][2] *= rho0 / a0[ci];
  }
}

void cell::attractiveForceUpdate() {
  resetForcesAndEnergy();
  shapeForces2D();
  vertexAttractiveForces2D_2();
}

void cell::wallForces(bool left, bool bottom, bool right, bool top, double& forceLeft, double& forceBottom, double& forceRight, double& forceTop, double appliedUniaxialPressure) {
  // compute particle-wall forces and wall-particle forces. Only the latter
  // exists, unless bool is
  //  set to true, in which case the wall can be pushed on. bool set to true =
  //  barostat e.g. for 3 fixed walls and 1 moving wall, set true false false
  //  false forceTop, etc., are forces on the walls due to newton's 3rd law +
  //  collisions.

  // TODO: compute wall-cell attraction energy per cell (pain to go from vi to
  // ci, so not doing this)
  bool collideTopOrRight, collideBottomOrLeft;
  int vi = 0, gi = 0;
  double cx = 0, cy = 0;
  double boxL, K = 1, fmag = 0;
  double shell, cut, s, distLower, distUpper, scaledDist, ftmp, f;
  double kint = (kc * l1) / (l2 - l1);
  forceTop = 0;
  forceBottom = 0;
  forceLeft = 0;
  forceRight = 0;

  // if any cells have their centers outside of the box, force them towards the
  // center
  for (int ci = 0; ci < NCELLS; ci++) {
    com2D(ci, cx, cy);
    if (cx < 0 + XL[0]|| cx > L[0]) { // 0 + XL[0] is the left boundary, L[0] = L[0] + XL[3] is the right boundary. XL lets the boundaries be movable
      for (int vi = 0; vi < nv[ci]; vi++) {
        gi = gindex(ci, vi);
        F[gi * NDIM] += -0.5 * (x[gi * NDIM] - L[0] / 2);
      }
    }
    if (cy < 0 + XL[1]|| cy > L[1]) {
      for (int vi = 0; vi < nv[ci]; vi++) {
        gi = gindex(ci, vi);
        F[gi * NDIM + 1] += -0.5 * (x[gi * NDIM + 1] - L[1] / 2);
      }
    }
  }

  for (int i = 0; i < vertDOF; i++) {
    // XL is positive if original left/bottom boundaries have moved right/up, negative if they've moved left/down, respectively
    vi = i / 2;
    boxL = L[i % NDIM]; 

    s = 2 * r[vi];
    cut = (1.0 + l1) * s;
    shell = (1.0 + l2) * s;
    distLower = fabs(x[i] - XL[i % NDIM]);  // XL[0] is new position of left x boundary, XL[1] is new position of bottom y boundary
    distUpper = fabs(boxL - x[i]);          // boxL is new position of right x or top y boundaries

    // check if vertex is within attracting distance
    if (distLower < shell) {
      // scaled distance
      scaledDist = distLower / s;
      // check within attraction ring
      if (distLower > cut) {
        // force scale
        ftmp = kint * (scaledDist - 1.0 - l2) / s;

        U += -0.5 * 1 * pow(1.0 + l2 - scaledDist, 2.0);
      } else {
        // force scale
        ftmp = kc * (1 - scaledDist) / s;

        U += 0.5 * kc * (pow(1.0 - scaledDist, 2.0) - l1 * l2);
      }
      f = ftmp * s;
      F[i] += f;
      if (i % NDIM == 0 && left)
        forceLeft -= f;
      if (i % NDIM == 1 && bottom)
        forceBottom -= f;
    }
    // check if vertex is within attracting distance
    if (distUpper < shell) {
      // scaled distance
      scaledDist = distUpper / s;
      // check within attraction ring
      if (distUpper > cut) {
        // force scale
        ftmp = kint * (scaledDist - 1.0 - l2) / s;

        U += -0.5 * 1 * pow(1.0 + l2 - scaledDist, 2.0);
      } else {
        // force scale
        ftmp = kc * (1 - scaledDist) / s;

        U += 0.5 * kc * (pow(1.0 - scaledDist, 2.0) - l1 * l2);
      }
      f = -ftmp * s;
      F[i] += f;
      if (i % NDIM == 0 && right)
        forceRight -= f;
      if (i % NDIM == 1 && top)
        forceTop -= f;
    }
  }

  if (L[0] < r[0]){
    cout << "forceLeft = " << forceLeft << ", added force = " << appliedUniaxialPressure * L[1] << '\n';
    cout << "L[0] = " << L[0] << " < r[0] = " << r[0] << ", there is no room left to compress or simclock < 410\n";
    cout << "XL[0] = " << XL[0] << '\n';
    //assert(false);
  }
    
  forceLeft += appliedUniaxialPressure * L[1];
  forceRight += -appliedUniaxialPressure * L[1];
}

void cell::cellPolarityForces(int ci, double k_polarity, std::string direction) {
  // compute energy and force representing cell polarity (due to cytoskeleton)
  // in its current state, this function penalizes any deviation between the y-values of the cell center and the cell's vertices
  // parameter ci is included in case I only want a subset of cells to have a polarity

  double rho0, cx, cy, rx, ry;
  double dDim, cDim, rDim;  // dimension-agnostic (x or y) coordinate of center-vertex separation, center, and vertex
  int gi, dimensionOffset;  // dimensionOffset adds to gi * NDIM in order to access x (0) or y (1) coordinates

  rho0 = sqrt(a0.at(0));  // define lengthscale off of the same particle every time
  com2D(ci, cx, cy);

  if (direction == "x") {
    dimensionOffset = 0;
    cDim = cx;
  } else if (direction == "y") {
    dimensionOffset = 1;
    cDim = cy;
  } else {
    cerr << "invalid input into cellPolarityForces : neither 'x' nor 'y' \n";
    assert(false);
  }

  for (int vi = 0; vi < nv.at(ci); vi++) {
    gi = gindex(ci, vi);
    rDim = x[NDIM * gi + dimensionOffset];
    dDim = rDim - cDim;

    if (pbc[dimensionOffset])
      dDim -= L[dimensionOffset] * round(dDim / L[dimensionOffset]);

    F[gi * NDIM + dimensionOffset] += -k_polarity * dDim;
    cellU[ci] += 0.5 * k_polarity * pow(dDim, 2);

    // does not factor into stress calculation
    // does not factor into field stress calculation (may do this later)
    U += 0.5 * k_polarity * pow(dDim, 2);
  }
}

void cell::vertexCompress2Target2D(dpmMemFn forceCall, double Ftol, double dt0, double phi0Target, double dphi0) {
  // same as for dpm, but with less printing at the end
  // local variables
  int it = 0, itmax = 1e4;
  double phi0 = vertexPreferredPackingFraction2D();
  double scaleFactor = 1.0, P, Sxy;

  // loop while phi0 < phi0Target
  while (phi0 < phi0Target && it < itmax) {
    // scale particle sizes
    scaleParticleSizes2D(scaleFactor);

    // update phi0
    phi0 = vertexPreferredPackingFraction2D();
    // relax configuration (pass member function force update)
    vertexFIRE2D(forceCall, Ftol, dt0);

    // get scale factor
    scaleFactor = sqrt((phi0 + dphi0) / phi0);

    // get updated pressure
    P = 0.5 * (stress[0] + stress[1]);
    Sxy = stress[2];

    // print to console
    if (it % 5 == 0) {
      cout << endl
           << endl;
      cout << "===============================" << endl;
      cout << "								"
           << endl;
      cout << " 	C O M P R E S S I O N 		" << endl;
      cout << "								"
           << endl;
      cout << "	P R O T O C O L 	  		" << endl;
      cout << "								"
           << endl;
      cout << "===============================" << endl;
      cout << endl;
      cout << "	** it 			= " << it << endl;
      cout << "	** phi0 curr	= " << phi0 << endl;
      if (phi0 + dphi0 < phi0Target)
        cout << "	** phi0 next 	= " << phi0 + dphi0 << endl;
      cout << "	** P 			= " << P << endl;
      cout << "	** Sxy 			= " << Sxy << endl;
      cout << "	** U 			= " << U << endl;
      // printConfiguration2D();
      cout << endl
           << endl;
    }

    // update iterate
    it++;
  }
}

void cell::simulateDampedWithWalls(dpmMemFn forceCall, double B, double dt0, double duration, double printInterval, bool wallsOn, bool leftOpen, bool bottomOpen, bool rightOpen, bool topOpen, double trueStrainRateX, double trueStrainRateY, double appliedUniaxialPressure) {
  // make sure velocities exist or are already initialized before calling this
  // assuming zero temperature - ignore thermostat (not implemented)
  // allow box lengths to move as a dynamical variable - rudimentary barostat,
  // doesn't matter for non-equilibrium anyway. 
  // NOTE: this part of the code begins to modify all 4 boundary lengths of the simulation through the variables XL. boundaries will no longer be based on the old convention bottom left = (0,0)

  // open bools indicate whether a wall is static (closed) or dynamic (open)
  // trueStrainRate (not engineering strain) multiplies corresponding box dimension by some value

  int i;
  double K, t0 = simclock;
  double temp_simclock = simclock;
  double FT, FB, FL, FR;
  double oldFT, oldFB, oldFL, oldFR;
  std::vector<double> boundaryDisplacement(NDIM, 0.0);


  // set time step magnitude
  setdt(dt0);
  int NPRINTSKIP = printInterval / dt;

  // loop over time, print energy
  while (simclock - t0 < duration) {
    XL[2] = L[0]; // set the dynamic boundary lengths to be equal to the current static boundary lengths. Both will be dynamic from here on out.
    XL[3] = L[1];
    // VV POSITION UPDATE

    for (i = 0; i < vertDOF; i++) {
      // update position
      x[i] += dt * v[i] + 0.5 * dt * dt * F[i];
      // recenter in box
      if (x[i] > L[i % NDIM] && pbc[i % NDIM])
        x[i] -= L[i % NDIM];
      else if (x[i] < 0 && pbc[i % NDIM])
        x[i] += L[i % NDIM];
    }

    // FORCE UPDATE
    std::vector<double> F_old = F;

    CALL_MEMBER_FN(*this, forceCall)
    ();  // calls main force routine

    // VV VELOCITY UPDATE #2
    for (i = 0; i < vertDOF; i++) {
      F[i] -= (B * v[i] + B * F_old[i] * dt / 2);
      F[i] /= (1 + B * dt / 2);
      v[i] += 0.5 * (F[i] + F_old[i]) * dt;
    }

    if (wallsOn == true) {
      // apply strain to walls
      boundaryDisplacement[0] = XL[2] * (1 - exp(-trueStrainRateX * dt));
      boundaryDisplacement[1] = XL[3] * (1 - exp(-trueStrainRateY * dt));

      // VV position update (4 walls)
      XL[0] += dt * VL[0] + 0.5 * dt * dt * FL - boundaryDisplacement[0]/2;
      XL[1] += dt * VL[1] + 0.5 * dt * dt * FB - boundaryDisplacement[1]/2;
      XL[2] += dt * VL[2] + 0.5 * dt * dt * FR + boundaryDisplacement[0]/2;
      XL[3] += dt * VL[3] + 0.5 * dt * dt * FT + boundaryDisplacement[1]/2;

      L[0] = XL[2]; 
      L[1] = XL[3];

      /*cout << "(velocity\tforce): left vs right\n";
      cout << VL[0] << '\t' << FL << '\n';
      cout << VL[2] << '\t' << FR << '\n';*/

      // VV force update (walls)
      oldFT = FT;
      oldFB = FB;
      oldFL = FL;
      oldFR = FR;
      // true false false false means 3 walls and an open top
      wallForces(leftOpen, bottomOpen, rightOpen, topOpen, FL, FB, FR, FT, appliedUniaxialPressure);
      // if bool is false in wallForces, then corresponding force is zero

      FL -= (B * VL[0] + B * oldFL * dt / 2);
      FL /= (1 + B * dt / 2);
      FB -= (B * VL[1] + B * oldFB * dt / 2);
      FB /= (1 + B * dt / 2);
      FR -= (B * VL[2] + B * oldFR * dt / 2);
      FR /= (1 + B * dt / 2);
      FT -= (B * VL[3] + B * oldFT * dt / 2);
      FT /= (1 + B * dt / 2);

      VL[0] += 0.5 * (FL + oldFL) * dt;
      VL[1] += 0.5 * (FB + oldFB) * dt;
      VL[2] += 0.5 * (FR + oldFR) * dt;
      VL[3] += 0.5 * (FT + oldFT) * dt;
    }

    // update sim clock
    simclock += dt;
    // print to console and file
    if (int(printInterval) != 0) {
      if (int((simclock - t0) / dt) % NPRINTSKIP == 0) {
        temp_simclock = simclock;
        // compute kinetic energy
        K = vertexKineticEnergy();

        // print to console
        cout << endl
             << endl;
        cout << "===============================" << endl;
        cout << "	D P M  						"
             << endl;
        cout << " 			 				"
                "	"
             << endl;
        cout << "		N P 0 (DAMPED) 				"
                "	"
             << endl;
        cout << "===============================" << endl;
        cout << endl;
        cout << "	** simclock - t0 / duration	= " << simclock - t0
             << " / " << duration << endl;
        cout << " ** simclock = " << simclock << endl;
        cout << "	** U 		= " << setprecision(12) << U << endl;
        cout << "	** K 		= " << setprecision(12) << K << endl;
        cout << "	** E 		= " << setprecision(12) << U + K
             << endl;

        if (enout.is_open()) {
          // print to energy file
          cout << "** printing energy" << endl;
          enout << setw(wnum) << left << simclock;
          enout << setw(wnum) << setprecision(12) << U;
          enout << setw(wnum) << setprecision(12) << K;
          enout << endl;
        }

        if (stressout.is_open()) {
          double shapeStressXX = 0.0, shapeStressYY = 0.0, shapeStressXY = 0.0;
          for (int ci = 0; ci < NCELLS; ci++) {
            shapeStressXX += fieldShapeStressCells[ci][0];
            shapeStressYY += fieldShapeStressCells[ci][1];
            shapeStressXY += fieldShapeStressCells[ci][2];
          }
          // print to stress file
          cout << "** printing stress" << endl;
          cout << "field shape stresses: " << -shapeStressXX << '\t'
               << -shapeStressYY << '\t' << -shapeStressXY << '\n';
          stressout << setw(wnum) << left << simclock;
          stressout << setw(wnum) << -stress[0];
          stressout << setw(wnum) << -stress[1];
          stressout << setw(wnum) << -stress[2];
          stressout << setw(wnum) << -shapeStressXX;
          stressout << setw(wnum) << -shapeStressYY;
          stressout << setw(wnum) << -shapeStressXY;
          stressout << endl;
        }

        // print to configuration only if position file is open
        if (posout.is_open()) {
          printConfiguration2D();
          cerr << "done printing in simulateDampedConstantPressure\n";
        }
      }
    }
  }
}

void cell::printConfiguration2D() {
  // overloaded to print out specific quantities of interest
  // local variables
  int ci, cj, vi, gi, ctmp, zc, zv;
  double xi, yi, dx, dy;

  // check if pos object is open
  if (!posout.is_open()) {
    cerr << "** ERROR: in printConfiguration2D, posout is not open, but "
            "function call will try to use. Ending here."
         << endl;
    exit(1);
  } else
    cout << "** In printConfiguration2D, printing particle positions to file..."
         << endl;

  // save box sizes

  double L_left = XL[0];
  double L_right = L[0];
  double L_bottom = XL[1];
  double L_top = L[1];
  //cout << L_left << '\t' << L_right << '\t' << L_bottom << '\t' << L_top << '\n';

  // print information starting information
  posout << setw(w) << left << "NEWFR"
         << " " << endl;
  posout << setw(w) << left << "NUMCL" << setw(w) << left << NCELLS << endl;
  posout << setw(w) << left << "PACKF" << setw(wnum) << setprecision(pnum)
         << left << vertexPackingFraction2D() << endl;

  // print box sizes
  posout << setw(w) << left << "BOXSZ";
  posout << setw(wnum) << setprecision(pnum) << left << L_left;
  posout << setw(wnum) << setprecision(pnum) << left << L_bottom;
  posout << setw(wnum) << setprecision(pnum) << left << L_right;
  posout << setw(wnum) << setprecision(pnum) << left << L_top;
  posout << endl;

  // print stress info
  posout << setw(w) << left << "STRSS";
  posout << setw(wnum) << setprecision(pnum) << left << -stress[0];
  posout << setw(wnum) << setprecision(pnum) << left << -stress[1];
  posout << setw(wnum) << setprecision(pnum) << left << -stress[2];
  posout << endl;

  // print simulation time
  posout << setw(w) << left << "TIME";
  posout << setw(wnum) << setprecision(pnum) << left << simclock;
  posout << endl;

  // print coordinate for rest of the cells
  for (ci = 0; ci < NCELLS; ci++) {
    // get cell contact data
    zc = 0;
    zv = 0;
    for (cj = 0; cj < NCELLS; cj++) {
      if (ci != cj) {
        // contact info from entry ci, cj
        if (ci < cj)
          ctmp = cij[NCELLS * ci + cj - (ci + 1) * (ci + 2) / 2];
        else
          ctmp = cij[NCELLS * cj + ci - (cj + 1) * (cj + 2) / 2];

        // add to contact information
        zv += ctmp;
        if (ctmp > 0)
          zc++;
      }
    }

    // cell information
    posout << setw(w) << left << "CINFO";
    posout << setw(w) << left << nv[ci];
    posout << setw(w) << left << zc;
    posout << setw(w) << left << zv;
    posout << setw(wnum) << left << a0[ci];
    posout << setw(wnum) << left << area(ci);
    posout << setw(wnum) << left << perimeter(ci);
    posout << setw(wnum) << left
           << -fieldStressCells[ci][0];  // virial contribution to stress tensor
                                         // comes with a minus sign
    posout << setw(wnum) << left << -fieldStressCells[ci][1];
    posout << setw(wnum) << left << -fieldStressCells[ci][2];
    posout << setw(wnum) << left
           << -(fieldStressCells[ci][0] + fieldStressCells[ci][1]);
    posout << setw(wnum) << left << -fieldShapeStressCells[ci][0];
    posout << setw(wnum) << left << -fieldShapeStressCells[ci][1];
    posout << setw(wnum) << left << -fieldShapeStressCells[ci][2];
    posout << setw(wnum) << left << cellU[ci];
    posout << endl;

    // get initial vertex positions
    gi = gindex(ci, 0);
    xi = x[NDIM * gi];
    yi = x[NDIM * gi + 1];

    // place back in box center
    if (pbc[0])
      xi = fmod(xi, L_right);
    if (pbc[1])
      yi = fmod(yi, L_top);

    posout << setw(w) << left << "VINFO";
    posout << setw(w) << left << ci;
    posout << setw(w) << left << 0;

    // output initial vertex information
    posout << setw(wnum) << setprecision(pnum) << right << xi;
    posout << setw(wnum) << setprecision(pnum) << right << yi;
    posout << setw(wnum) << setprecision(pnum) << right << gi;
    posout << setw(wnum) << setprecision(pnum) << right << r[gi];
    posout << setw(wnum) << setprecision(pnum) << right << l0[gi];
    posout << setw(wnum) << setprecision(pnum) << right << t0[gi];
    posout << endl;

    // vertex information for next vertices
    for (vi = 1; vi < nv[ci]; vi++) {
      // get global vertex index for next vertex
      gi++;

      // get next vertex positions
      dx = x[NDIM * gi] - xi;
      if (pbc[0])
        dx -= L_right * round(dx / L_right);
      xi += dx;

      dy = x[NDIM * gi + 1] - yi;
      if (pbc[1])
        dy -= L_top * round(dy / L_top);
      yi += dy;

      // Print indexing information
      posout << setw(w) << left << "VINFO";
      posout << setw(w) << left << ci;
      posout << setw(w) << left << vi;

      // output vertex information
      posout << setw(wnum) << setprecision(pnum) << right << xi;
      posout << setw(wnum) << setprecision(pnum) << right << yi;
      posout << setw(wnum) << setprecision(pnum) << right << gi;
      posout << setw(wnum) << setprecision(pnum) << right << r[gi];
      posout << setw(wnum) << setprecision(pnum) << right << l0[gi];
      posout << setw(wnum) << setprecision(pnum) << right << t0[gi];
      posout << setw(wnum) << setprecision(pnum) << right << fieldStress[gi][0];
      posout << setw(wnum) << setprecision(pnum) << right << fieldStress[gi][1];
      posout << setw(wnum) << setprecision(pnum) << right << fieldStress[gi][2];
      posout << setw(wnum) << setprecision(pnum) << right << fieldShapeStress[gi][0];
      posout << setw(wnum) << setprecision(pnum) << right << fieldShapeStress[gi][1];
      posout << setw(wnum) << setprecision(pnum) << right << fieldShapeStress[gi][2];
      posout << endl;
    }
  }

  // print end frame
  posout << setw(w) << left << "ENDFR"
         << " " << endl;

  cout << "leaving printPosCell\n";
}