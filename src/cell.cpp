#include "cell.h"

// namespace
using namespace Eigen;
using namespace std;

void cell::moveVertex(int gi, double xpos, double ypos) {
  x[gi * NDIM] = xpos;
  x[gi * NDIM + 1] = ypos;
}

void cell::removeCellIDFromInteractionMatrix(int cellID) {
  cellTypeIntMat.erase(cellTypeIntMat.begin() + cellID);  // remove ci row
  // remove ci column
  for (int i = 0; i < cellTypeIntMat.size(); i++) {
    if (cellTypeIntMat[i].size() > cellID) {
      cellTypeIntMat[i].erase(cellTypeIntMat[i].begin() + cellID);
    }
  }
}

void cell::addCellIDToInteractionMatrix(int cellID) {
  assert(cellTypeIntMat.size() >= cellID);
  assert(cellTypeIntMat[0].size() >= cellID);
  // add column
  cout << "adding column of zeros to Interaction Matrix\n";
  for (auto& row : cellTypeIntMat) {
    row.insert(row.begin() + cellID, 0);
  }
  // add row
  vector<double> tempVec(cellTypeIntMat[0].size(), 0.0);
  cellTypeIntMat.insert(cellTypeIntMat.begin() + cellID, tempVec);
  cout << "adding row of zeros to InteractionMatrix\n please use setCellCellIntMat to add interactions to the newly added cell";
}

void cell::printInteractionMatrix() {
  cout << "printing interaction matrix\n";
  for (auto i : cellTypeIntMat) {
    for (auto j : i)
      cout << j << '\t';
    cout << '\n';
  }
}

void cell::shapeForces2D() {
  // uses cellID to determine whether to compute certain parts of the shape forces
  // for example, the ECM cellID shouldn't have any bending energy.

  // local variables
  int ci, ci_real, gi, vi, nvtmp;
  double fa, fli, flim1, fb, cx, cy, xi, yi;
  double rho0, l0im1, l0i, a0tmp, atmp;
  double dx, dy, da, dli, dlim1, dtim1, dti, dtip1;
  double lim2x, lim2y, lim1x, lim1y, lix, liy, lip1x, lip1y, li, lim1;
  double rim2x, rim2y, rim1x, rim1y, rix, riy, rip1x, rip1y, rip2x, rip2y;
  double nim1x, nim1y, nix, niy, sinim1, sini, sinip1, cosim1, cosi, cosip1;
  double ddtim1, ddti;
  double forceX, forceY;
  double unwrappedX, unwrappedY;

  // loop over vertices, add to force
  rho0 = sqrt(a0.at(0));
  ci = 0;
  for (gi = 0; gi < NVTOT; gi++) {
    // -- Area force (and get cell index ci)
    if (ci < NCELLS) {
      if (gi == szList[ci]) {
        //  shape information
        nvtmp = nv[ci];
        a0tmp = a0[ci];

        // preferred segment length of last segment
        l0im1 = l0[im1[gi]];

        // compute area deviation
        atmp = area(ci);
        da = (atmp / a0tmp) - 1.0;

        // update potential energy
        U += 0.5 * ka * (da * da);
        cellU[ci] += 0.5 * ka * (da * da);

        // shape force parameters
        fa = ka * da * (rho0 / a0tmp);
        fb = kb * rho0;

        // compute cell center of mass
        xi = x[NDIM * gi];
        yi = x[NDIM * gi + 1];
        cx = xi;
        cy = yi;
        for (vi = 1; vi < nvtmp; vi++) {
          // get distances between vim1 and vi
          dx = x[NDIM * (gi + vi)] - xi;
          dy = x[NDIM * (gi + vi) + 1] - yi;
          if (pbc[0])
            dx -= L[0] * round(dx / L[0]);
          if (pbc[1])
            dy -= L[1] * round(dy / L[1]);

          // add to centers
          xi += dx;
          yi += dy;

          cx += xi;
          cy += yi;
        }
        cx /= nvtmp;
        cy /= nvtmp;

        // get coordinates relative to center of mass
        rix = x[NDIM * gi] - cx;
        riy = x[NDIM * gi + 1] - cy;

        // get prior adjacent vertices
        rim2x = x[NDIM * im1[im1[gi]]] - cx;
        rim2y = x[NDIM * im1[im1[gi]] + 1] - cy;
        if (pbc[0])
          rim2x -= L[0] * round(rim2x / L[0]);
        if (pbc[1])
          rim2y -= L[1] * round(rim2y / L[1]);

        rim1x = x[NDIM * im1[gi]] - cx;
        rim1y = x[NDIM * im1[gi] + 1] - cy;
        if (pbc[0])
          rim1x -= L[0] * round(rim1x / L[0]);
        if (pbc[1])
          rim1y -= L[1] * round(rim1y / L[1]);

        // get prior segment vectors
        lim2x = rim1x - rim2x;
        lim2y = rim1y - rim2y;

        lim1x = rix - rim1x;
        lim1y = riy - rim1y;

        // increment cell index
        ci++;
      }
    }

    ci_real = ci - 1;  // ci is one too large from above code, but more efficient than calling cindices repeatedly

    // unwrapped vertex coordinate
    unwrappedX = cx + rix;
    unwrappedY = cy + riy;

    // preferred segment length
    l0i = l0[gi];

    // get next adjacent vertices
    rip1x = x[NDIM * ip1[gi]] - cx;
    rip1y = x[NDIM * ip1[gi] + 1] - cy;
    if (pbc[0])
      rip1x -= L[0] * round(rip1x / L[0]);
    if (pbc[1])
      rip1y -= L[1] * round(rip1y / L[1]);

    // -- Area force (comes from a cross product)
    forceX = 0.5 * fa * (rim1y - rip1y);
    forceY = 0.5 * fa * (rip1x - rim1x);

    // cellIndex is separate from ci because ci is used to get adjacent vertices
    // int cellIndex, vertexIndex;
    // cindices(cellIndex, vertexIndex, gi);
    F[NDIM * gi] += forceX;
    F[NDIM * gi + 1] += forceY;

    fieldShapeStress[gi][0] += unwrappedX * forceX;
    fieldShapeStress[gi][1] += unwrappedY * forceY;
    fieldShapeStress[gi][2] += unwrappedX * forceY;

    // -- Perimeter force
    lix = rip1x - rix;
    liy = rip1y - riy;

    // segment lengths
    lim1 = sqrt(lim1x * lim1x + lim1y * lim1y);
    li = sqrt(lix * lix + liy * liy);

    // segment deviations (note: m is prior vertex, p is next vertex i.e. gi - 1, gi + 1 mod the right number of vertices)
    dlim1 = (lim1 / l0im1) - 1.0;
    dli = (li / l0i) - 1.0;

    // segment forces
    flim1 = kl * (rho0 / l0im1);
    fli = kl * (rho0 / l0i);

    // add to forces
    forceX = (fli * dli * lix / li) - (flim1 * dlim1 * lim1x / lim1);
    forceY = (fli * dli * liy / li) - (flim1 * dlim1 * lim1y / lim1);
    F[NDIM * gi] += forceX;
    F[NDIM * gi + 1] += forceY;

    // note - Andrew here, confirmed that the shape stress matrix is diagonal as written
    fieldShapeStress[gi][0] += unwrappedX * forceX;
    fieldShapeStress[gi][1] += unwrappedY * forceY;
    fieldShapeStress[gi][2] += unwrappedX * forceY;

    // update potential energy
    U += 0.5 * kl * (dli * dli);
    cellU[ci_real] += 0.5 * kl * (dli * dli);

    // -- Bending force
    if (kb > 0 && cellID[ci_real] == 0) {
      // if (kb > 0) {
      //  get ip2 for third angle
      rip2x = x[NDIM * ip1[ip1[gi]]] - cx;
      rip2y = x[NDIM * ip1[ip1[gi]] + 1] - cy;
      if (pbc[0])
        rip2x -= L[0] * round(rip2x / L[0]);
      if (pbc[1])
        rip2y -= L[1] * round(rip2y / L[1]);

      // get last segment length
      lip1x = rip2x - rip1x;
      lip1y = rip2y - rip1y;

      // get angles
      sinim1 = lim1x * lim2y - lim1y * lim2x;
      cosim1 = lim1x * lim2x + lim1y * lim2y;

      sini = lix * lim1y - liy * lim1x;
      cosi = lix * lim1x + liy * lim1y;

      sinip1 = lip1x * liy - lip1y * lix;
      cosip1 = lip1x * lix + lip1y * liy;

      // get normal vectors
      nim1x = lim1y;
      nim1y = -lim1x;

      nix = liy;
      niy = -lix;

      // get change in angles
      dtim1 = atan2(sinim1, cosim1) - t0[im1[gi]];
      dti = atan2(sini, cosi) - t0[gi];
      dtip1 = atan2(sinip1, cosip1) - t0[ip1[gi]];

      // get delta delta theta's
      ddtim1 = (dti - dtim1) / (lim1 * lim1);
      ddti = (dti - dtip1) / (li * li);

      // add to force
      F[NDIM * gi] += fb * (ddtim1 * nim1x + ddti * nix);
      F[NDIM * gi + 1] += fb * (ddtim1 * nim1y + ddti * niy);

      // update potential energy
      U += 0.5 * kb * (dti * dti);
      cellU[ci] += 0.5 * kb * (dti * dti);
    }

    // update old coordinates
    rim2x = rim1x;
    rim1x = rix;
    rix = rip1x;

    rim2y = rim1y;
    rim1y = riy;
    riy = rip1y;

    // update old segment vectors
    lim2x = lim1x;
    lim2y = lim1y;

    lim1x = lix;
    lim1y = liy;

    // update old preferred segment length
    l0im1 = l0i;
  }

  // normalize per-cell stress by preferred cell area
  for (int ci = 0; ci < NCELLS; ci++) {
    for (int vi = 0; vi < nv[ci]; vi++) {
      int gi = gindex(ci, vi);
      fieldShapeStressCells[ci][0] += fieldShapeStress[gi][0];
      fieldShapeStressCells[ci][1] += fieldShapeStress[gi][1];
      fieldShapeStressCells[ci][2] += fieldShapeStress[gi][2];
    }
    // nondimensionalize the stress
    fieldShapeStressCells[ci][0] *= rho0 / a0[ci];
    fieldShapeStressCells[ci][1] *= rho0 / a0[ci];
    fieldShapeStressCells[ci][2] *= rho0 / a0[ci];
  }
}

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

void cell::attractiveWithPolarityForceAndWallCrawlingUpdate() {
  resetForcesAndEnergy();
  shapeForces2D();
  // give each cell a y-penalized polarity
  for (int ci = 0; ci < NCELLS; ci++) {
    cellPolarityForces(ci, 0.1, "y");
  }
  vertexAttractiveForces2D_2();
  wallCrawlingForces();
}

void cell::attractiveWallCrawlingForceUpdate() {
  resetForcesAndEnergy();
  shapeForces2D();
  vertexAttractiveForces2D_2();
  wallCrawlingForces();
}

void cell::repulsiveForceUpdateWithPolyWall() {
  resetForcesAndEnergy();
  shapeForces2D();
  vertexRepulsiveForces2D();
  for (int i = 0; i < poly_bd_x.size(); i++) {
    evaluatePolygonalWallForces(poly_bd_x[i], poly_bd_y[i]);
  }
}

void cell::attractiveForceUpdateWithPolyWall() {
  resetForcesAndEnergy();
  shapeForces2D();
  vertexAttractiveForces2D_2();
  for (int i = 0; i < poly_bd_x.size(); i++) {
    evaluatePolygonalWallForces(poly_bd_x[i], poly_bd_y[i]);
  }
}

void cell::attractiveForceUpdate() {
  resetForcesAndEnergy();
  shapeForces2D();
  vertexAttractiveForces2D_2();
}

void cell::attractiveForceUpdateWithCrawling() {
  resetForcesAndEnergy();
  shapeForces2D();
  vertexAttractiveForces2D_2();
  brownianCrawlingUpdate();
}

void cell::brownianCrawlingUpdate() {
  int gi;
  // propel at constant speed v_0
  // v_0 should have a force scale comparable to the shape energy? or other energy scale. check that and make v0_ABP scale with one of the spring constants
  // printf("v0_ABP = %f, kc * rho0 / 2*r = %f \n", v0_ABP, kc * sqrt(a0[0]) / (2 * r[0]));
  for (int ci = 0; ci < NCELLS; ci++) {
    double director = psi[ci];
    for (int vi = 0; vi < nv[ci]; vi++) {
      gi = gindex(ci, vi);
      F[gi * NDIM] += v0_ABP * cos(director);
      F[gi * NDIM + 1] += v0_ABP * sin(director);
    }
  }

  // shake up the polarities
  directorDiffusion();
}

void cell::directorDiffusion() {
  double r1, r2, grv;
  double Dr0 = 1 / tau_ABP;
  for (int ci = 0; ci < NCELLS; ci++) {
    // propagate diffusion of directors psi
    r1 = drand48();
    r2 = drand48();
    grv = sqrt(-2.0 * log(r1)) * cos(2.0 * PI * r2);
    // if flag is present, do not diffuse.
    psi[ci] += sqrt(dt * 2.0 * Dr0) * grv;
    psi[ci] -= 2 * PI * round(psi[ci] / (2 * PI));
  }
}

void cell::attractiveForceUpdatePrint(double& forceX, double& forceY, double& energy) {
  resetForcesAndEnergy();
  // shapeForces2D();
  energy = 0.0;
  vertexAttractiveForces2D_test(energy);
  forceX = F[0];
  forceY = F[1];
}

void cell::attractiveSmoothForceUpdate() {
  resetForcesAndEnergy();
  shapeForces2D();
  circuloLineAttractiveForces();
}

/*void cell::attractiveForceUpdateSmoothPrint(double& forceX, double& forceY, double& energy) {
  resetForcesAndEnergy();
  // shapeForces2D();
  energy = 0.0;
  smoothAttractiveForces2D_test(energy);
  forceX = F[0];
  forceY = F[1];
}*/

void cell::circuloLineAttractiveForces() {
  // altered from vertexAttractiveForces2D_2, here we use vertex-vertex and vertex-line segment distances to make a smooth interaction
  // models sliding adhesion and repulsion.
  // different from epi2D smooth forces because I comment out vnn
  int ci, cj, gi, gj, vi, vj, bi, bj, pi, pj;
  double sij, rij, rho0;
  double d, rx, ry, dx, dy;         // distance and components of separation between projection of point p onto line segment v-w
  double d_arg, y10, x10, d21;      // for calculating 3-body forces for projection 1 (vertex-line-segment)
  double projection;                // parameterized projection value. if between [0,1] then it's circulo-line, if < 0 or > 1 then it is either nothing or end-end.
  double endCapAngle, endEndAngle;  // endCapAngle is PI minus interior angle of vertices, endEndAngle is between interaction centers and circulo-line endpoints
  double ftmp, fx, fy, energytmp;
  bool isConcaveInteraction, isConvexInteraction, isSelfInteraction = false;
  int sign = 1;  // used to flip the sign of force and energy in the concave interaction case for negative endCap vertex-vertex interactions

  // attraction shell parameters
  double shellij, cutij, xij, kint = (kc * l1) / (l2 - l1);

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

        isSelfInteraction = false;
        // calculate self-penetration: if self penetrating, compute self-repulsion and move on
        if (ci == cj) {
          isSelfInteraction = true;
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
                if (rij < sij) {
                  ftmp = kc * (1 - (rij / sij)) * (1 / sij);
                  cellU[ci] += 0.5 * kc * pow((1 - (rij / sij)), 2.0);
                  U += 0.5 * kc * pow((1 - (rij / sij)), 2.0);
                  // force elements
                  fx = ftmp * (dx / rij);
                  fy = ftmp * (dy / rij);

                  // add to forces
                  F[NDIM * gi] -= fx;
                  F[NDIM * gi + 1] -= fy;

                  F[NDIM * gj] += fx;
                  F[NDIM * gj + 1] += fy;
                }
              }
            }
          }
        }
        // calculate d, rx, ry, and projection.
        // d = distance from point [gi] to line segment (next[gj] to gj)
        // rx, ry = x,y components of d
        // projection = parametrization value of the projection of gi onto the line segment.

        for (int swapii = 0; swapii < 2; swapii++) {
          d = linePointDistancesAndProjection(x[NDIM * im1[gj]], x[NDIM * im1[gj] + 1], x[NDIM * gj], x[NDIM * gj + 1], x[NDIM * gi], x[NDIM * gi + 1], rx, ry, projection, x10, y10);
          if (!isSelfInteraction) {
            if (projection < 1 || d < shellij) {
              // check that the projection falls within the interacting portion of vertex i
              // checking vertex j means checking a circulo-line between j (P=1) and j-1 (P=0)
              // projection <= 0 means that p0 projection onto p1-p2 is behind p1, which is a potential v-v interaction
              // 0 < projection < 1 means that p0 projection onto p1-p2 falls between p1 and p2, so it's a vertex-line-segment contact
              // projection > 1 means p0 projection falls ahead of p2, so ignore
              // unless d < shellij, in which case we need to check if i-(j-1) is close enough to need an inverse v-v interaction to patch discontinuities
              calculateSmoothInteraction(rx, ry, sij, shellij, cutij, kint, kc, gi, gj, projection, ci, cj);
            }
          }

          int gi_temp = gi;
          gi = gj;
          gj = gi_temp;
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

          isSelfInteraction = false;
          // calculate self-penetration: if self penetrating, compute self-repulsion and move on
          if (ci == cj) {
            isSelfInteraction = true;
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
                  if (rij < sij) {
                    ftmp = kc * (1 - (rij / sij)) * (1 / sij);
                    cellU[ci] += 0.5 * kc * pow((1 - (rij / sij)), 2.0);
                    U += 0.5 * kc * pow((1 - (rij / sij)), 2.0);
                    // force elements
                    fx = ftmp * (dx / rij);
                    fy = ftmp * (dy / rij);

                    // add to forces
                    F[NDIM * gi] -= fx;
                    F[NDIM * gi + 1] -= fy;

                    F[NDIM * gj] += fx;
                    F[NDIM * gj + 1] += fy;
                  }
                }
              }
            }
          }

          // calculate d, rx, ry, and projection.
          // d = distance from point [gi] to line segment (next[gj] to gj)
          // rx, ry = x,y components of d
          // projection = parametrization value of the projection of gi onto the line segment.

          for (int swapii = 0; swapii < 2; swapii++) {
            d = linePointDistancesAndProjection(x[NDIM * im1[gj]], x[NDIM * im1[gj] + 1], x[NDIM * gj], x[NDIM * gj + 1], x[NDIM * gi], x[NDIM * gi + 1], rx, ry, projection, x10, y10);
            if (!isSelfInteraction) {
              if (projection < 1 || d < shellij) {
                // check that the projection falls within the interacting portion of vertex i
                // checking vertex j means checking a circulo-line between j (P=1) and j-1 (P=0)
                // projection <= 0 means that p0 projection onto p1-p2 is behind p1, which is a potential v-v interaction
                // 0 < projection < 1 means that p0 projection onto p1-p2 falls between p1 and p2, so it's a vertex-line-segment contact
                // projection > 1 means p0 projection falls ahead of p2, so ignore
                // unless d < shellij, in which case we need to check if i-(j-1) is close enough to need an inverse v-v interaction to patch discontinuities
                calculateSmoothInteraction(rx, ry, sij, shellij, cutij, kint, kc, gi, gj, projection, ci, cj);
              }
            }
            int gi_temp = gi;
            gi = gj;
            gj = gi_temp;
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

void cell::calculateSmoothInteraction(double& rx, double& ry, double& sij, double& shellij, double& cutij, double& kint, double& kc, int& gi, int& gj, double& projection, int& ci, int& cj) {
  double rij, xij, ftmp, energytmp;
  double dx, dy, fx, fy;
  double endEndAngle, endCapAngle, endEndAngle2;
  bool isConvexInteractionWithCap, isConcaveInteractionWithCap, isCapInteraction;
  bool isConcaveInteraction, isInverseInteraction;
  double d_arg, y21, x21, y20, x20, y10, x10, norm_P12, prefix, prefix2;  // for calculating 3-body forces for projection 1 (vertex-line-segment)
  int sign = 1;

  dx = -rx;
  if (pbc[0])
    dx -= L[0] * round(dx / L[0]);
  if (dx < shellij) {  // check that d is within the interaction shell
    dy = -ry;
    if (pbc[1])
      dy -= L[1] * round(dy / L[1]);
    if (dy < shellij) {
      rij = sqrt(dx * dx + dy * dy);
      if (rij < shellij) {
        // scaled distance
        xij = rij / sij;
        // pick force based on vertex-vertex distance
        if (rij > cutij) {
          // force scale for inner layer of interaction shell
          ftmp = kint * (xij - 1.0 - l2) / sij;
          energytmp = -0.5 * kint * pow(1.0 + l2 - xij, 2.0);
        } else {
          // force scale for outer layer of interaction shell
          ftmp = kc * (1 - xij) / sij;
          energytmp = 0.5 * kc * (pow(1.0 - xij, 2.0) - l1 * l2);
        }
        // endEndAngle is the angle between the separation between interacting vertices and the endcap edge closest to the circulo-line.
        // endCapAngle is the angle between the far edge of the endcap and the near edge of the endcap
        int left = gj;             // i
        int middle = im1[gj];      // i-1
        int right = im1[im1[gj]];  // i-2

        double drx = x[left * NDIM] - x[middle * NDIM];
        double dry = x[left * NDIM + 1] - x[middle * NDIM + 1];
        double drx_prev = x[right * NDIM] - x[middle * NDIM];
        double dry_prev = x[right * NDIM + 1] - x[middle * NDIM + 1];
        double vv_rx = x[NDIM * gi] - x[NDIM * middle];
        double vv_ry = x[NDIM * gi + 1] - x[NDIM * middle + 1];

        if (pbc[0]) {
          drx -= L[0] * round(drx / L[0]);
          drx_prev -= L[0] * round(drx_prev / L[0]);
          vv_rx -= L[0] * round(vv_rx / L[0]);
        }

        if (pbc[1]) {
          dry -= L[1] * round(dry / L[1]);
          dry_prev -= L[1] * round(dry_prev / L[1]);
          vv_ry -= L[1] * round(vv_ry / L[1]);
        }

        // note that using r dot dr_prev only gives the correct endEndangle for the convex case
        //  because the vertex-line distance in the convex case will be measured from the end of the projection cutoff which is P = 0
        // When generalizing to the concave case, the measurement needs to explicitly be from vertex-vertex (vv) at P = 0

        // r is the separation from vertex to vertex
        // A is the vector from middle vertex to the outer (when convex) corner of middle-left's circulo line segment
        // B is the vector from middle vertex to the outer (when convex) corner of middle-right's circulo line segment

        endEndAngle = atan2(vv_rx * dry - drx * vv_ry, vv_rx * drx + vv_ry * dry);
        if (endEndAngle < 0) {
          endEndAngle += 2 * PI;
        }
        endEndAngle = endEndAngle - PI / 2;  // theta', measures angle from r ccw to A

        endCapAngle = atan2(drx_prev * dry - drx * dry_prev, drx_prev * drx + dry_prev * dry);
        if (endCapAngle < 0)
          endCapAngle += 2 * PI;
        endCapAngle = endCapAngle - PI;  // phi, measures angle from B to A

        endEndAngle2 = PI - endEndAngle + fabs(endCapAngle);  // theta'', measures angle from -B to r

        endEndAngle2 -= 2 * PI * (endEndAngle2 > 2 * PI);
        endEndAngle2 += 2 * PI * (endEndAngle2 < 0);  // roll over [0,2pi]

        // if convex: theta' is within the vertex end-cap region, so compute a bumpy interaction
        // if concave: theta' is within the interior vertex end-cap region, so compute a bumpy interaction (only happens with attraction and severe concavity)
        isConvexInteractionWithCap = (endEndAngle >= 0 && endEndAngle <= endCapAngle);
        isConcaveInteractionWithCap = (endCapAngle < 0 && endEndAngle > PI - fabs(endCapAngle) && endEndAngle < PI);
        // isConcaveInteractionWithCap = false;
        isCapInteraction = (isConvexInteractionWithCap || isConcaveInteractionWithCap);

        // end-cap region is on the inside of the particle, and theta' indicates a vertex might be in the concave overlap region
        double endEndAngleNewRange = (endEndAngle - 2 * PI * (endEndAngle > PI));  // roll over to [-pi, pi]
        isConcaveInteraction = (endCapAngle < 0 && endEndAngleNewRange < 0 && endEndAngleNewRange >= endCapAngle);

        // unstable concave overlap is present, or unstable convex overlap is present, so compute a correction energy for numerical stability
        isInverseInteraction = (isConcaveInteraction || (endCapAngle > 0 && (endEndAngle2 < endCapAngle)));

        if (projection > 0 && projection < 1) {  // projection less than 1 and greater than 0, so projection is on the main line segment
          // Force on particle 0,1,2 is determined by F = - dU/dr = (partials) dU/dr * <dr/dxi , dr/dyi>
          // 3-body contact, 6 forces (3 pairs of forces)
          // y21, x21, y20, x20, y10, x10, norm_P12, d_arg
          int g2 = im1[gj];
          int g2_ind = NDIM * g2;
          int g1_ind = NDIM * gj;
          x21 = x[g2_ind] - x[g1_ind];
          y21 = x[g2_ind + 1] - x[g1_ind + 1];
          x20 = x[g2_ind] - x[NDIM * gi];
          y20 = x[g2_ind + 1] - x[NDIM * gi + 1];
          x10 = x[g1_ind] - x[NDIM * gi];
          y10 = x[g1_ind + 1] - x[NDIM * gi + 1];
          if (pbc[0]) {
            x21 -= L[0] * round(x21 / L[0]);
            x20 -= L[0] * round(x20 / L[0]);
            x10 -= L[0] * round(x10 / L[0]);
          }
          if (pbc[1]) {
            y21 -= L[1] * round(y21 / L[1]);
            y20 -= L[1] * round(y20 / L[1]);
            y10 -= L[1] * round(y10 / L[1]);
          }
          d_arg = x21 * y10 - x10 * y21;
          norm_P12 = sqrt(pow(x21, 2) + pow(y21, 2));
          prefix = d_arg / fabs(d_arg) / norm_P12;
          prefix2 = fabs(d_arg) / pow(norm_P12, 3);

          F[NDIM * gi] += ftmp * prefix * y21;
          F[NDIM * gi + 1] += ftmp * prefix * -x21;

          F[NDIM * gj] += ftmp * (prefix * -y20 + x21 * prefix2);
          F[NDIM * gj + 1] += ftmp * (prefix * x20 + y21 * prefix2);

          F[NDIM * g2] += ftmp * (prefix * y10 - x21 * prefix2);
          F[NDIM * g2 + 1] += ftmp * (prefix * -x10 - y21 * prefix2);

          cellU[ci] += energytmp / 2;
          cellU[cj] += energytmp / 2;
          U += energytmp;

          // add to virial stress - not including this code now because I haven't worked out the stress of a 3-body interaction
        }
        // projection is either on the endpoint or outside the endpoint, i.e. not on the line segment
        if ((projection <= 0 && isCapInteraction) || (projection > 0 && isInverseInteraction)) {
          //|| (!isConcaveInteraction && isInverseInteraction)) {
          // pure 2-body contact determined by angles and distances between contact points or by self interaction
          if (isInverseInteraction) {
            //  if concave, compute interaction between vertex and inverse vertex. sign = -1 to explicitly demonstrate that the only difference between vertex and inverse vertex is a minus sign
            //  have to reevaluate the distances because previous code uses vertex-line distance, whereas we need vertex-vertex distance
            //   for the special case of concave interactions
            sign = 0;
            // if (not close enough to interact) sign = 0;
            dx = -vv_rx;
            if (pbc[0])
              dx -= L[0] * round(dx / L[0]);
            if (dx < shellij) {
              dy = -vv_ry;
              if (pbc[1])
                dy -= L[1] * round(dy / L[1]);
              if (dy < shellij) {
                rij = sqrt(dx * dx + dy * dy);
                if (rij < shellij) {
                  sign = -1;  // confirmed contact with negative potential vertex, so flip sign
                  xij = rij / sij;
                  if (rij > cutij) {
                    ftmp = kint * (xij - 1.0 - l2) / sij;
                    energytmp = -0.5 * kint * pow(1.0 + l2 - xij, 2.0);
                  } else {
                    ftmp = kc * (1 - xij) / sij;
                    energytmp = 0.5 * kc * (pow(1.0 - xij, 2.0) - l1 * l2);
                  }
                }
              }
            }
          } else {  // if (isCapInteraction) {
            sign = 1;
          }
          // above, if concave, altered dx, dy, rij, xij, ftmp, energytmp.

          // force elements
          fx = sign * ftmp * (dx / rij);  // dx/rij comes from the chain rule (dU/dx1 = dU/dr * dr/dx1)
          fy = sign * ftmp * (dy / rij);
          F[NDIM * gi] -= fx;
          F[NDIM * gi + 1] -= fy;

          F[NDIM * middle] += fx;  // the relevant degrees of freedom are gi and middle, not gi and gj. This is due to the choice of counterclockwise coordinates, and projection < 0 condition.
          F[NDIM * middle + 1] += fy;

          cellU[ci] += sign * energytmp / 2;
          cellU[cj] += sign * energytmp / 2;
          U += sign * energytmp;

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
          fieldStress[middle][0] += -dx / 2 * fx;
          fieldStress[middle][1] += -dy / 2 * fy;
          fieldStress[middle][2] += -0.5 * (dx / 2 * fy + dy / 2 * fx);
        }

        if (ci > cj)
          cij[NCELLS * cj + ci - (cj + 1) * (cj + 2) / 2]++;
        else if (ci < cj)
          cij[NCELLS * ci + cj - (ci + 1) * (ci + 2) / 2]++;
      }
    }
  }
}

void cell::vertexAttractiveForces2D_test(double& energy) {
  // attractive forces calculation, focuses on forces and energy involving vertex gi = 0
  energy = 0.0;

  // altered from dpm attractive force code, because it works with larger l2
  // in cell class, also includes cell type specific interactions through cellTypeIntMat
  // values. (warning: probably won't work with bending. Francesco says it should be fine though, haven't tested.) local variables
  int ci, cj, gi, gj, vi, vj, bi, bj, pi, pj;
  double sij, rij, dx, dy, rho0;
  double ftmp, fx, fy;
  double cellTypeIntModifier = 1.0, cellTypeIntModifier_repulsion = 1.0;

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
        cellTypeIntModifier = cellTypeIntMat[cellID[ci]][cellID[cj]];
        cellTypeIntModifier_repulsion = 1.0;

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
                  ftmp = cellTypeIntModifier_repulsion * kc * (1 - (rij / sij)) * (rho0 / sij);
                  cellU[ci] += 0.5 * cellTypeIntModifier_repulsion * kc * pow((1 - (rij / sij)), 2.0);
                } else
                  ftmp = 0;
              } else if (rij > cutij) {
                // force scale
                ftmp = kint * cellTypeIntModifier * (xij - 1.0 - l2) / sij;

                // increase potential energy
                U += -0.5 * cellTypeIntModifier * kint * pow(1.0 + l2 - xij, 2.0);
                cellU[ci] += -0.5 * cellTypeIntModifier * kint * pow(1.0 + l2 - xij, 2.0) / 2.0;
                cellU[cj] += -0.5 * cellTypeIntModifier * kint * pow(1.0 + l2 - xij, 2.0) / 2.0;
                if (gi == 0 || gj == 0) {
                  cout << "rij > cut interaction between " << gi << '\t' << gj << '\n';
                  energy += -0.5 * cellTypeIntModifier * kint * pow(1.0 + l2 - xij, 2.0);
                }
              } else {
                // force scale
                ftmp = cellTypeIntModifier * kc * (1 - xij) / sij;

                // increase potential energy
                U += 0.5 * cellTypeIntModifier * kc * (pow(1.0 - xij, 2.0) - l1 * l2);
                cellU[ci] += 0.5 * cellTypeIntModifier * kc * (pow(1.0 - xij, 2.0) - l1 * l2) / 2.0;
                cellU[cj] += 0.5 * cellTypeIntModifier * kc * (pow(1.0 - xij, 2.0) - l1 * l2) / 2.0;
                if (gi == 0 || gj == 0) {
                  cout << "rij < cut interaction between " << gi << '\t' << gj << '\n';
                  energy += 0.5 * cellTypeIntModifier * kc * (pow(1.0 - xij, 2.0) - l1 * l2);
                }
              }

              // force elements
              fx = ftmp * (dx / rij);
              fy = ftmp * (dy / rij);

              if (gi == 0 || gj == 0)
                cout << "fx, fy = " << fx << '\t' << fy << '\n';

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
          cellTypeIntModifier = cellTypeIntMat[cellID[ci]][cellID[cj]];
          cellTypeIntModifier_repulsion = 1.0;

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
                    ftmp = cellTypeIntModifier_repulsion * kc * (1 - (rij / sij)) * (rho0 / sij);
                    cellU[ci] += 0.5 * cellTypeIntModifier_repulsion * kc * pow((1 - (rij / sij)), 2.0);
                  } else
                    ftmp = 0;
                } else if (rij > cutij) {
                  // force scale
                  ftmp = kint * cellTypeIntModifier * (xij - 1.0 - l2) / sij;

                  // increase potential energy
                  U += -0.5 * cellTypeIntModifier * kint * pow(1.0 + l2 - xij, 2.0);
                  cellU[ci] += -0.5 * cellTypeIntModifier * kint * pow(1.0 + l2 - xij, 2.0) / 2.0;
                  cellU[cj] += -0.5 * cellTypeIntModifier * kint * pow(1.0 + l2 - xij, 2.0) / 2.0;
                  if (gi == 0 || gj == 0) {
                    cout << "rij > cut interaction between " << gi << '\t' << gj << '\n';
                    energy += -0.5 * cellTypeIntModifier * kint * pow(1.0 + l2 - xij, 2.0);
                  }
                } else {
                  // force scale
                  ftmp = cellTypeIntModifier * kc * (1 - xij) / sij;

                  // increase potential energy
                  U += 0.5 * cellTypeIntModifier * kc * (pow(1.0 - xij, 2.0) - l1 * l2);
                  cellU[ci] += 0.5 * cellTypeIntModifier * kc * (pow(1.0 - xij, 2.0) - l1 * l2) / 2.0;
                  cellU[cj] += 0.5 * cellTypeIntModifier * kc * (pow(1.0 - xij, 2.0) - l1 * l2) / 2.0;
                  if (gi == 0 || gj == 0) {
                    cout << "rij < cut interaction between " << gi << '\t' << gj << '\n';
                    energy += 0.5 * cellTypeIntModifier * kc * (pow(1.0 - xij, 2.0) - l1 * l2);
                  }
                }

                // force elements
                fx = ftmp * (dx / rij);
                fy = ftmp * (dy / rij);

                if (gi == 0 || gj == 0) {
                  cout << "fx, fy = " << fx << '\t' << fy << '\n';
                }

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

void cell::vertexAttractiveForces2D_2() {
  // altered from dpm attractive force code, because it works with larger l2
  // in cell class, also includes cell type specific interactions through cellTypeIntMat
  // values. (warning: probably won't work with bending. Francesco says it should be fine though, haven't tested.) local variables
  int ci, cj, gi, gj, vi, vj, bi, bj, pi, pj;
  double sij, rij, dx, dy, rho0;
  double ftmp, utmp, fx, fy;
  double cellTypeIntModifier = 1.0, cellTypeIntModifier_repulsion = 1.0;

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
        // cout << "ci, cj, cellIDci, cellIDcj = " << ci << '\t' << cj << '\t' << cellID[ci] << '\t' << cellID[cj] << '\n';
        // cout << "cellTypeIntMat size = " << cellTypeIntMat.size() << '\t' << cellTypeIntMat[0].size() << '\n';
        cellTypeIntModifier = cellTypeIntMat[cellID[ci]][cellID[cj]];
        cellTypeIntModifier_repulsion = 1.0;

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
                  ftmp = cellTypeIntModifier_repulsion * kc * (1 - (rij / sij)) * (rho0 / sij);
                  utmp = 0.5 * cellTypeIntModifier_repulsion * kc * pow((1 - (rij / sij)), 2.0);
                  U += utmp;
                  cellU[ci] += utmp;
                } else
                  ftmp = 0;
              } else if (rij > cutij) {
                // force scale
                ftmp = kint * cellTypeIntModifier * (xij - 1.0 - l2) / sij;
                utmp = -0.5 * cellTypeIntModifier * kint * pow(1.0 + l2 - xij, 2.0);
                U += utmp;
                cellU[ci] += utmp / 2.0;
                cellU[cj] += utmp / 2.0;
              } else {
                // force scale
                ftmp = cellTypeIntModifier * kc * (1 - xij) / sij;
                utmp = 0.5 * cellTypeIntModifier * kc * (pow(1.0 - xij, 2.0) - l1 * l2);
                U += utmp;
                cellU[ci] += utmp / 2.0;
                cellU[cj] += utmp / 2.0;
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
          cellTypeIntModifier = cellTypeIntMat[cellID[ci]][cellID[cj]];
          cellTypeIntModifier_repulsion = 1.0;

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
                    ftmp = cellTypeIntModifier_repulsion * kc * (1 - (rij / sij)) * (rho0 / sij);
                    utmp = 0.5 * cellTypeIntModifier_repulsion * kc * pow((1 - (rij / sij)), 2.0);
                    U += utmp;
                    cellU[ci] += utmp;
                  } else
                    ftmp = 0;
                } else if (rij > cutij) {
                  // force scale
                  ftmp = kint * cellTypeIntModifier * (xij - 1.0 - l2) / sij;
                  utmp = -0.5 * cellTypeIntModifier * kint * pow(1.0 + l2 - xij, 2.0);
                  U += utmp;
                  cellU[ci] += utmp / 2.0;
                  cellU[cj] += utmp / 2.0;
                } else {
                  // force scale
                  ftmp = cellTypeIntModifier * kc * (1 - xij) / sij;
                  utmp = 0.5 * cellTypeIntModifier * kc * (pow(1.0 - xij, 2.0) - l1 * l2);
                  U += utmp;
                  cellU[ci] += utmp / 2.0;
                  cellU[cj] += utmp / 2.0;
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

void cell::wallForces(bool left, bool bottom, bool right, bool top, double& forceLeft, double& forceBottom, double& forceRight, double& forceTop, double appliedUniaxialPressure) {
  // compute particle-wall forces and wall-particle forces. Only the latter
  // exists, unless bool is
  //  set to true, in which case the wall can be pushed on. bool set to true =
  //  barostat e.g. for 3 fixed walls and 1 moving wall, set true false false
  //  false forceTop, etc., are forces on the walls due to newton's 3rd law +
  //  collisions.
  // also fills cellTouchesWallsLeft and Right, which tells whether a cell is touching either wall

  // TODO: compute wall-cell attraction energy per cell (pain to go from vi to
  // ci, so not doing this)
  bool collideTopOrRight, collideBottomOrLeft;
  int vi = 0, gi = 0, ci_temp, vi_temp;
  double cx = 0, cy = 0;
  double boxL, force_multiplier = 10, fmag = 0;
  if (fabs(appliedUniaxialPressure) < 1e-10)
    force_multiplier = 1;  // if no applied pressure, then don't modify wall forces
  double kc_temp = kc * force_multiplier;
  double shell, cut, s, distLower, distUpper, scaledDist, ftmp, f;
  double kint = (kc_temp * l1) / (l2 - l1);
  forceTop = 0;
  forceBottom = 0;
  forceLeft = 0;
  forceRight = 0;

  std::fill(cellTouchesWallsLeft.begin(), cellTouchesWallsLeft.end(), false);  // reset cellTouchesWallsLeft and Right
  cellTouchesWallsRight = cellTouchesWallsLeft;

  // if any cells have their centers outside of the box, force them towards the
  // center
  for (int ci = 0; ci < NCELLS; ci++) {
    com2D(ci, cx, cy);
    if (cx < 0 + XL[0] || cx > L[0]) {  // 0 + XL[0] is the left boundary, L[0] = L[0] + XL[3] is the right boundary. XL lets the boundaries be movable
      for (int vi = 0; vi < nv[ci]; vi++) {
        gi = gindex(ci, vi);
        F[gi * NDIM] += -0.5 * (x[gi * NDIM] - L[0] / 2);
      }
    }
    if (cy < 0 + XL[1] || cy > L[1]) {
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
      if (i % 2 == 0) {  // only care about touching left or right walls at this stage
        cindices(ci_temp, vi_temp, vi);
        cellTouchesWallsLeft[ci_temp] = true;
      }
      // scaled distance
      scaledDist = distLower / s;
      // check within attraction ring
      if (distLower > cut) {
        // force scale
        ftmp = kint * (scaledDist - 1.0 - l2) / s;

        U += -0.5 * kint * pow(1.0 + l2 - scaledDist, 2.0);
      } else {
        // force scale
        ftmp = kc_temp * (1 - scaledDist) / s;

        U += 0.5 * kc_temp * (pow(1.0 - scaledDist, 2.0) - l1 * l2);
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
      if (i % 2 == 0) {
        cindices(ci_temp, vi_temp, vi);
        cellTouchesWallsRight[ci_temp] = true;
      }
      // scaled distance
      scaledDist = distUpper / s;
      // check within attraction ring
      if (distUpper > cut) {
        // force scale
        ftmp = kint * (scaledDist - 1.0 - l2) / s;

        U += -0.5 * 1 * pow(1.0 + l2 - scaledDist, 2.0);
      } else {
        // force scale
        ftmp = kc_temp * (1 - scaledDist) / s;

        U += 0.5 * kc_temp * (pow(1.0 - scaledDist, 2.0) - l1 * l2);
      }
      f = -ftmp * s;
      F[i] += f;
      if (i % NDIM == 0 && right)
        forceRight -= f;
      if (i % NDIM == 1 && top)
        forceTop -= f;
    }
  }

  forceLeft += appliedUniaxialPressure * L[1];
  forceRight += -appliedUniaxialPressure * L[1];

  if (L[0] < r[0]) {
    cout << "forceLeft = " << forceLeft << ", added force = " << appliedUniaxialPressure * L[1] << '\n';
    cout << "L[0] = " << L[0] << " < r[0] = " << r[0] << ", there is no room left to compress or simclock < 410\n";
    cout << "XL[0] = " << XL[0] << '\n';
    // assert(false);
  }
}

void cell::wallCrawlingForces() {
  // compute crawling forces for cells that are touching either the left or right walls, on their polarized end away from the wall.
  double kint = 10 * (kc * l1) / (l2 - l1);
  for (int ci = 0; ci < NCELLS; ci++) {
    double nearestXValueToCenter = 0.0;
    int nearestGiToCenter = -1;
    if (cellTouchesWallsLeft[ci]) {
      nearestXValueToCenter = XL[0];
      for (int vi = 0; vi < nv[ci]; vi++) {
        int gi = gindex(ci, vi);
        if (x[gi * NDIM] > nearestXValueToCenter) {
          nearestXValueToCenter = x[gi * NDIM];
          nearestGiToCenter = gi;
        }
      }
    } else if (cellTouchesWallsRight[ci]) {
      nearestXValueToCenter = L[0];
      for (int vi = 0; vi < nv[ci]; vi++) {
        int gi = gindex(ci, vi);
        if (x[gi * NDIM] < nearestXValueToCenter) {
          nearestXValueToCenter = x[gi * NDIM];
          nearestGiToCenter = gi;
        }
      }
    }
    bool nearestVertexToCenterFound = (nearestGiToCenter != -1);
    if (nearestVertexToCenterFound) {
      // cout << "adding force to vertex gi = " << nearestGiToCenter << ", at XValue " << nearestXValueToCenter << ", simclock = " << simclock << '\n';
      // cout << "force is equal to " << kint * ((XL[0] + L[0])/2.0 - nearestXValueToCenter) << '\n';
      F[nearestGiToCenter * NDIM] += kint * ((XL[0] + L[0]) / 2.0 - nearestXValueToCenter);
      // cout << "total force is equal to " << F[nearestGiToCenter*NDIM] << '\n';
    }
  }
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

// boundary routines
void cell::replacePolyWallWithDP(int numCellTypes) {
  // take a polyWall in poly_bd_x and poly_bd_y, create data to replace the polywall with a DP of the same size and shape, then delete the polywall.
  // must account for and modify all vectors with size dependent on NCELLS and NVTOT, such as szList
  double polyPerimeter = 0.0;
  vector<double> dp_x, dp_y;

  for (int tissueIt = 0; tissueIt < poly_bd_x.size(); tissueIt++) {
    vector<double> poly_x = poly_bd_x[tissueIt];
    vector<double> poly_y = poly_bd_y[tissueIt];
    for (int i = 0; i < poly_x.size(); i++) {
      polyPerimeter += sqrt(pow(poly_x[i] - poly_x[(i + 1) % poly_x.size()], 2) + pow(poly_y[i] - poly_y[(i + 1) % poly_y.size()], 2));
    }
    int numVerts = polyPerimeter / l0[0];
    int cellTypeIndex = numCellTypes - 1;
    // interpolate poly_bd_x,poly_bd_y into dp_x,dp_y
    vector<double> dp_coords = resample_polygon(poly_x, poly_y, polyPerimeter, numVerts);
    for (int i = 0; i < dp_coords.size(); i += 2) {
      dp_x.push_back(dp_coords[i]);
      dp_y.push_back(dp_coords[i + 1]);
    }
    addDP(numVerts, dp_x, dp_y, cellTypeIndex, numCellTypes);

    polyPerimeter = 0.0;
    dp_x = {};
    dp_y = {};
  }
}

void cell::addDP(int numVerts, vector<double>& dp_x, vector<double>& dp_y, int cellTypeIndex, int numCellTypes) {
  // add a DP with cell type = cellID, with numVerts vertices, each of which are located according to (dp_x,dp_y)
  // must account for and modify all vectors with size dependent on NCELLS and NVTOT, such as szList
  int gi, vim1, vip1;
  double dist, polyArea = 0.0;

  NCELLS += 1;
  NVTOT += numVerts;
  vertDOF += NDIM * numVerts;
  szList.push_back(szList.back() + nv.back());
  nv.push_back(numVerts);
  // next: a0, l0, t0, r; im1, ip1, x, v, f in initializeVertexShapeParameters and initializeVertexIndexing
  gi = szList[NCELLS - 1];
  cout << "gi = " << gi << '\t' << ", which should = szList[NCELLS] = " << szList[NCELLS - 1] << "\t, with NCELLS = " << NCELLS << '\n';

  int j = dp_x.size() - 1;
  for (int vi = 0; vi < dp_x.size(); vi++) {
    dist = sqrt(pow(dp_x[vi] - dp_x[(vi + 1) % dp_x.size()], 2) + pow(dp_y[vi] - dp_y[(vi + 1) % dp_y.size()], 2));
    l0.push_back(dist);
    t0.push_back(0.0);
    r.push_back(0.5 * dist);

    // assign vertex degrees of freedom
    x.push_back(dp_x[vi]);
    x.push_back(dp_y[vi]);
    for (int d = 0; d < NDIM; d++) {
      v.push_back(0.0);
      F.push_back(0.0);
    }

    // calculate area of DP for a0 via shoelace method
    polyArea += (dp_x[j] + dp_x[vi]) * (dp_y[j] - dp_y[vi]);
    j = vi;
  }
  polyArea = abs(polyArea / 2.0);
  a0.push_back(polyArea);

  // save list of adjacent vertices
  im1.resize(NVTOT);
  ip1.resize(NVTOT);
  for (int ci = 0; ci < NCELLS; ci++) {
    // vertex indexing
    for (int vi = 0; vi < nv.at(ci); vi++) {
      // wrap local indices
      vim1 = (vi - 1 + nv.at(ci)) % nv.at(ci);
      vip1 = (vi + 1) % nv.at(ci);

      // get global wrapped indices
      gi = gindex(ci, vi);
      im1.at(gi) = gindex(ci, vim1);
      ip1.at(gi) = gindex(ci, vip1);
    }
  }

  // linked list variables
  //  list should have size NVTOT + 1
  list.resize(NVTOT + 1);

  // other variables
  // cellU, stresses, cij in main DPM class need to be resized
  cellU.resize(NCELLS);
  // resize instead of push_back leads to exception, not sure why.
  /*fieldStress.resize(NVTOT, vector<double>(3));
  fieldStressCells.resize(NCELLS, vector<double>(3));
  fieldShapeStress.resize(NVTOT, vector<double>(3));
  fieldShapeStressCells.resize(NCELLS, vector<double>(3));*/
  for (int i = 0; i < numVerts; i++) {
    fieldStress.push_back(vector<double>(3, 0.0));
    fieldShapeStress.push_back(vector<double>(3, 0.0));
  }
  fieldStressCells.push_back(vector<double>(3, 0.0));
  fieldShapeStressCells.push_back(vector<double>(3, 0.0));
  cij.resize(NCELLS * (NCELLS - 1) / 2);

  // cellID, cellTouchesLeft, cellTouchesRight in Cell class need to be resized, then I should be good
  cellID.push_back(cellTypeIndex);
  if (cellTypeIndex >= numCellTypes) {
    cout << "error: cellTypeIndex specified in addDP is larger than any existing cellTypes. This is unexpected, either add a cellType or choose a smaller cellType.\n";
    assert(false);
  }
  cellTouchesWallsLeft.push_back(false);
  cellTouchesWallsRight.push_back(false);
}

// routines

// eventually, call this 4x in order to replace initializeFourTransverseTissues
void cell::initializeTransverseTissue(double phi0, double Ftol) {
  // initialize starting positions and all related data for a large DP (ECM) and numCellsInside smaller DPs
  // isFixedBoundary is an optional bool argument that tells cells to stay away from boundary during initialization
  // aspect ratio is L[0]/L[1]
  int i, d, ci, cj, vi, vj, gi, cellDOF = NDIM * NCELLS, cumNumCells = 0;
  int numEdges = 10;  // number of edges in the polygonal walls to approximate a circle
  double areaSum, xtra = 1.1;

  // local disk vectors
  vector<double> drad(NCELLS, 0.0);
  vector<double> dpos(cellDOF, 0.0);
  vector<double> dv(cellDOF, 0.0);
  vector<double> dF(cellDOF, 0.0);

  // print to console
  cout << "** initializing particle positions using 2D SP model and FIRE relaxation ..." << endl;

  // initialize stress field
  initializeFieldStress();

  // initial transverse geometry vectors in units of some lengthscale
  int numTissues = 1, totalNumCellsCheck = 0;
  double totalArea = 0.0;
  vector<double> cx = {1 / 2.0}, cy = {1 / 2.0};
  vector<double> tissueRadii = {1 / 2.0}, cellFractionPerTissue, numCellsInTissue;
  for (int i = 0; i < tissueRadii.size(); i++) {
    totalArea += tissueRadii[i] * tissueRadii[i];
    cout << totalArea << '\t' << tissueRadii[i] << '\n';
  }

  for (int i = 0; i < tissueRadii.size(); i++) {
    cellFractionPerTissue.push_back(tissueRadii[i] * tissueRadii[i] / totalArea);
    numCellsInTissue.push_back(round(cellFractionPerTissue[i] * NCELLS));
    cout << "cellFraction = " << tissueRadii[i] * tissueRadii[i] / totalArea << '\n';
    cout << "initializing " << NCELLS << " cells, with tissue " << i << " having cell fraction = " << cellFractionPerTissue[i] << '\n';
    cout << "NCELLS * cellFraction = " << NCELLS * cellFractionPerTissue[i] << ", which is " << round(NCELLS * cellFractionPerTissue[i]) << " when rounded\n";
    totalNumCellsCheck += round(NCELLS * cellFractionPerTissue[i]);
  }

  assert(totalNumCellsCheck == NCELLS);

  // initialize box size based on packing fraction
  areaSum = 0.0;
  for (ci = 0; ci < NCELLS; ci++)
    areaSum += a0.at(ci) + 0.25 * PI * pow(l0.at(ci), 2.0) * (0.5 * nv.at(ci) - 1);

  // set box size : phi_0 = areaSum / A => A = areaSum/phi_0 which gives us the following formulas for L
  for (d = 0; d < NDIM; d++) {
    L.at(d) = pow(4 / PI * areaSum * cellFractionPerTissue[0] / phi0, 1.0 / NDIM);
  }

  // multiply lengths by lengthscale
  for (int i = 0; i < numTissues; i++) {
    double tissueLengthscale = L[0];
    cx[i] *= tissueLengthscale;
    cy[i] *= tissueLengthscale;
    tissueRadii[i] *= tissueLengthscale;
  }

  ofstream boundaryStream("polyBoundary.txt");  // clear boundary text file, for visualizing or for debugging
  ofstream ipStream("initPosSP.txt");           // clear initialPos text file, for visualizing or for debugging
  ofstream ip2Stream("initPosSP2.txt");         // clear initialPos text file, for visualizing or for debugging
  for (int n = 0; n < numTissues; n++) {
    cout << "initializing cell centers randomly but rejecting if further than R from the center for tissue " << n << "out of " << numTissues - 1 << "\n";
    double scale_radius = 1.1;                   // make the polygon radius slightly larger so that it encompasses the circle that points are initialized in
    poly_bd_x.push_back(std::vector<double>());  // make new data for generateCircularBoundary to write a polygon
    poly_bd_y.push_back(std::vector<double>());
    generateCircularBoundary(numEdges, scale_radius * tissueRadii[n], cx[n], cy[n], poly_bd_x[n], poly_bd_y[n]);

    for (i = cumNumCells; i < cumNumCells + numCellsInTissue[n]; i++) {
      cout << "i = " << i << "\t, dpos.size() = " << dpos.size() << ", cumNumCells = " << cumNumCells << "\t, numCellsInTissue[n] = " << numCellsInTissue[n] << '\n';
      double dpos_x = tissueRadii[n] * (2 * drand48() - 1) + cx[n], dpos_y = tissueRadii[n] * (2 * drand48() - 1) + cy[n];
      while (pow(dpos_x - cx[n], 2) + pow(dpos_y - cy[n], 2) > pow(tissueRadii[n] / (scale_radius), 2)) {
        dpos_x = tissueRadii[n] * (2 * drand48() - 1) + cx[n];
        dpos_y = tissueRadii[n] * (2 * drand48() - 1) + cy[n];
      }
      dpos.at(i * NDIM) = dpos_x;
      dpos.at(i * NDIM + 1) = dpos_y;
      ipStream << dpos[i * NDIM] << '\t' << dpos[i * NDIM + 1] << '\n';
      cellID[i] = n;
      cout << "cell " << i << " is in tissue number " << n << "with cell type = " << cellID[i] << '\n';
      cout << " with position " << dpos_x << '\t' << dpos_y << '\n';
    }
    cout << "before adding to cumNumCells\n";
    cumNumCells += numCellsInTissue[n];
    cout << "just after adding to cumNumCells\n";
  }

  cout << "setting radii of SP disks\n";
  // set radii of SP disks
  for (ci = 0; ci < NCELLS; ci++) {
    xtra = 1.1;  // disks should have radius similar to the final particle radius, or could modify vrad[i] condition in wall calculation later
    drad.at(ci) = xtra * sqrt((2.0 * a0.at(ci)) / (nv.at(ci) * sin(2.0 * PI / nv.at(ci))));
    cout << "drad = " << drad[ci] << '\n';
  }

  // FIRE VARIABLES
  double P = 0;
  double fnorm = 0;
  double vnorm = 0;
  double alpha = alpha0;

  double dt0 = 1e-2;
  double dtmax = 10 * dt0;
  double dtmin = 1e-8 * dt0;

  int npPos = 0;
  int npNeg = 0;

  int fireit = 0;
  double fcheck = 10 * Ftol;

  // interaction variables
  double rij, sij, dtmp, ftmp, vftmp;
  double dr[NDIM];

  // initial step size
  dt = dt0;

  // loop until force relaxes
  while ((fcheck > Ftol) && fireit < itmax) {
    // FIRE step 1. Compute P
    P = 0.0;
    for (i = 0; i < cellDOF; i++)
      P += dv[i] * dF[i];

    // FIRE step 2. adjust simulation based on net motion of degrees of freedom
    if (P > 0) {
      // increase positive counter
      npPos++;

      // reset negative counter
      npNeg = 0;

      // alter simulation if enough positive steps have been taken
      if (npPos > NMIN) {
        // change time step
        if (dt * finc < dtmax)
          dt *= finc;

        // decrease alpha
        alpha *= falpha;
      }
    } else {
      // reset positive counter
      npPos = 0;

      // increase negative counter
      npNeg++;

      // check if simulation is stuck
      if (npNeg > NNEGMAX) {
        cerr << "	** ERROR: During initial FIRE minimization, P < 0 for too long, so ending." << endl;
        exit(1);
      }

      // take half step backwards, reset velocities
      for (i = 0; i < cellDOF; i++) {
        // take half step backwards
        dpos[i] -= 0.5 * dt * dv[i];

        // reset velocities
        dv[i] = 0.0;
      }

      // decrease time step if past initial delay
      if (fireit > NDELAY) {
        // decrease time step
        if (dt * fdec > dtmin)
          dt *= fdec;

        // reset alpha
        alpha = alpha0;
      }
    }

    // FIRE step 3. First VV update
    for (i = 0; i < cellDOF; i++)
      dv[i] += 0.5 * dt * dF[i];

    // FIRE step 4. adjust velocity magnitude
    fnorm = 0.0;
    vnorm = 0.0;
    for (i = 0; i < cellDOF; i++) {
      fnorm += dF[i] * dF[i];
      vnorm += dv[i] * dv[i];
    }
    fnorm = sqrt(fnorm);
    vnorm = sqrt(vnorm);
    if (fnorm > 0) {
      for (i = 0; i < cellDOF; i++)
        dv[i] = (1 - alpha) * dv[i] + alpha * (vnorm / fnorm) * dF[i];
    }

    // FIRE step 4. Second VV update
    for (i = 0; i < cellDOF; i++) {
      dpos[i] += dt * dv[i];
      dF[i] = 0.0;
    }

    // FIRE step 5. Update forces
    for (ci = 0; ci < NCELLS; ci++) {
      for (cj = ci + 1; cj < NCELLS; cj++) {
        // contact distance
        sij = drad[ci] + drad[cj];

        // true distance
        rij = 0.0;
        for (d = 0; d < NDIM; d++) {
          // get distance element
          dtmp = dpos[NDIM * cj + d] - dpos[NDIM * ci + d];
          if (pbc[d])
            dtmp -= L[d] * round(dtmp / L[d]);

          // add to true distance
          rij += dtmp * dtmp;

          // save in distance array
          dr[d] = dtmp;
        }
        rij = sqrt(rij);

        // check distances
        if (rij < sij) {
          // force magnitude
          ftmp = kc * (1.0 - (rij / sij)) / sij;

          // add to vectorial force
          for (d = 0; d < NDIM; d++) {
            vftmp = ftmp * (dr[d] / rij);
            dF[NDIM * ci + d] -= vftmp;
            dF[NDIM * cj + d] += vftmp;
          }
        }
      }
    }
    // FIRE step 4.1 Compute wall force
    for (int k = 0; k < poly_bd_x.size(); k++) {
      std::vector<double> poly_x = poly_bd_x[k];
      std::vector<double> poly_y = poly_bd_y[k];
      int n = poly_x.size();
      double distanceParticleWall, Rx, Ry, dw, K = 1;
      double bound_x1, bound_x2, bound_y1, bound_y2;
      // loop over boundary bars
      // loop over particles
      //  compute particle-boundary bar overlaps
      //  if overlap, Fx += K * dw * Rx/R, where K is a constant, dw = diameter/2 - R, Rx = x - px, R = sqrt(Rx^2 + Ry^2)
      for (int bound_i = 0; bound_i < n; bound_i++) {
        // use distanceLineAndPoint to get R, Rx, and Ry
        bound_x1 = poly_x[bound_i];
        bound_x2 = poly_x[(bound_i + 1) % n];
        bound_y1 = poly_y[bound_i];
        bound_y2 = poly_y[(bound_i + 1) % n];
        for (i = 0; i < cellDOF / NDIM; i++) {
          distanceParticleWall = distanceLinePointComponents(bound_x1, bound_y1, bound_x2, bound_y2, dpos[i * NDIM], dpos[i * NDIM + 1], Rx, Ry);
          dw = drad[i] - distanceParticleWall;
          if (distanceParticleWall <= drad[i]) {
            dF[i * NDIM] += K * dw * Rx / distanceParticleWall;
            dF[i * NDIM + 1] += K * dw * Ry / distanceParticleWall;
          }
        }
      }
    }
    // FIRE step 5. Final VV update
    for (i = 0; i < cellDOF; i++)
      dv[i] += 0.5 * dt * dF[i];

    // update forces to check
    fcheck = 0.0;
    for (i = 0; i < cellDOF; i++)
      fcheck += dF[i] * dF[i];
    fcheck = sqrt(fcheck / NCELLS);

    // print to console
    if (fireit % NSKIP == 0) {
      cout << endl
           << endl;
      cout << "===========================================" << endl;
      cout << "		I N I T I A L  S P 			" << endl;
      cout << " 	F I R E 						" << endl;
      cout << "		M I N I M I Z A T I O N 	" << endl;
      cout << "===========================================" << endl;
      cout << endl;
      cout << "	** fireit = " << fireit << endl;
      cout << "	** fcheck = " << fcheck << endl;
      cout << "	** fnorm = " << fnorm << endl;
      cout << "	** vnorm = " << vnorm << endl;
      cout << "	** dt = " << dt << endl;
      cout << "	** P = " << P << endl;
      cout << "	** Pdir = " << P / (fnorm * vnorm) << endl;
      cout << "	** alpha = " << alpha << endl;
    }

    // update iterate
    fireit++;
  }
  // check if FIRE converged
  if (fireit == itmax) {
    cout << "	** FIRE minimization did not converge, fireit = " << fireit << ", itmax = " << itmax << "; ending." << endl;
    exit(1);
  } else {
    cout << endl
         << endl;
    cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl;
    cout << "===========================================" << endl;
    cout << " 	F I R E 						" << endl;
    cout << "		M I N I M I Z A T I O N 	" << endl;
    cout << "	C O N V E R G E D! 				" << endl
         << endl;

    cout << "	(for initial disk minimization) " << endl;
    cout << "===========================================" << endl;
    cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl;
    cout << endl;
    cout << "	** fireit = " << fireit << endl;
    cout << "	** fcheck = " << fcheck << endl;
    cout << "	** vnorm = " << vnorm << endl;
    cout << "	** dt = " << dt << endl;
    cout << "	** P = " << P << endl;
    cout << "	** alpha = " << alpha << endl;
  }

  // initialize vertex positions based on cell centers
  for (ci = 0; ci < NCELLS; ci++) {
    ip2Stream << dpos[ci * NDIM] << '\t' << dpos[ci * NDIM + 1] << '\n';
    for (vi = 0; vi < nv.at(ci); vi++) {
      // get global vertex index
      gi = gindex(ci, vi);

      // length from center to vertex minus 1/2 the vertex radius to prevent overlaps
      dtmp = sqrt((2.0 * a0.at(ci)) / (nv.at(ci) * sin((2.0 * PI) / nv.at(ci)))) - r[gi] / 2.0;

      // set positions
      x.at(NDIM * gi) = dtmp * cos((2.0 * PI * vi) / nv.at(ci)) + dpos.at(NDIM * ci) + 1e-2 * l0[gi] * drand48();
      x.at(NDIM * gi + 1) = dtmp * sin((2.0 * PI * vi) / nv.at(ci)) + dpos.at(NDIM * ci + 1) + 1e-2 * l0[gi] * drand48();
    }
  }
}

void cell::initializeFourTransverseTissues(double phi0, double Ftol) {
  cout << "entered initializeTransverseTissue\n";
  // initialize the starting positions and all related data for a large deformable particle and numCellsInside smaller deformable particles
  // isFixedBoundary is an optional bool argument that tells cells to stay away from the boundary during initialization
  // aspectRatio is the ratio L[0] / L[1]
  int i, d, ci, cj, vi, vj, gi, cellDOF = NDIM * NCELLS, cumNumCells = 0;
  int numEdges = 10;  // number of edges in the polygonal walls to approximate a circle
  double areaSum, xtra = 1.1;

  // local disk vectors
  vector<double> drad(NCELLS, 0.0);
  vector<double> dpos(cellDOF, 0.0);
  vector<double> dv(cellDOF, 0.0);
  vector<double> dF(cellDOF, 0.0);

  // print to console
  cout << "** initializing particle positions using 2D SP model and FIRE relaxation ..." << endl;

  // initialize stress field
  initializeFieldStress();

  // initial transverse geometry vectors in units of some lengthscale
  int numTissues = 4, totalNumCellsCheck = 0;
  double totalArea = 0.0;
  // double offset = 1/8.0/2.0
  double offset = -1 / 8.0;
  double yoffset = offset * 2;
  vector<double> cx = {1 / 2.0, 5 / 4.0 + offset, 2.0 + 2 * offset, 5 / 4.0 + offset}, cy = {1 / 2.0, 3 / 8.0, 1 / 2.0, 5 / 4.0 + yoffset};
  vector<double> tissueRadii = {1 / 2.0, 1 / 4.0, 1 / 2.0, 1 / 2.0}, cellFractionPerTissue, numCellsInTissue;
  for (int i = 0; i < tissueRadii.size(); i++) {
    totalArea += tissueRadii[i] * tissueRadii[i];
    cout << totalArea << '\t' << tissueRadii[i] << '\n';
  }

  for (int i = 0; i < tissueRadii.size(); i++) {
    cellFractionPerTissue.push_back(tissueRadii[i] * tissueRadii[i] / totalArea);
    numCellsInTissue.push_back(round(cellFractionPerTissue[i] * NCELLS));
    cout << "cellFraction = " << tissueRadii[i] * tissueRadii[i] / totalArea << '\n';
    cout << "initializing " << NCELLS << " cells, with tissue " << i << " having cell fraction = " << cellFractionPerTissue[i] << '\n';
    cout << "NCELLS * cellFraction = " << NCELLS * cellFractionPerTissue[i] << ", which is " << round(NCELLS * cellFractionPerTissue[i]) << " when rounded\n";
    totalNumCellsCheck += round(NCELLS * cellFractionPerTissue[i]);
  }

  assert(totalNumCellsCheck == NCELLS);

  // initialize box size based on packing fraction
  areaSum = 0.0;
  for (ci = 0; ci < NCELLS; ci++)
    areaSum += a0.at(ci) + 0.25 * PI * pow(l0.at(ci), 2.0) * (0.5 * nv.at(ci) - 1);

  // set box size : phi_0 = areaSum / A => A = areaSum/phi_0 which gives us the following formulas for L
  for (d = 0; d < NDIM; d++) {
    L.at(d) = pow(4 / PI * areaSum * cellFractionPerTissue[0] / phi0, 1.0 / NDIM);
  }

  // multiply lengths by lengthscale
  for (int i = 0; i < numTissues; i++) {
    double tissueLengthscale = L[0];
    cx[i] *= tissueLengthscale;
    cy[i] *= tissueLengthscale;
    tissueRadii[i] *= tissueLengthscale;
  }

  ofstream boundaryStream("polyBoundary.txt");  // clear boundary text file, for visualizing or for debugging
  ofstream ipStream("initPosSP.txt");           // clear initialPos text file, for visualizing or for debugging
  ofstream ip2Stream("initPosSP2.txt");         // clear initialPos text file, for visualizing or for debugging
  for (int n = 0; n < numTissues; n++) {
    cout << "initializing cell centers randomly but rejecting if further than R from the center for tissue " << n << "out of " << numTissues - 1 << "\n";
    double scale_radius = 1.1;                   // make the polygon radius slightly larger so that it encompasses the circle that points are initialized in
    poly_bd_x.push_back(std::vector<double>());  // make new data for generateCircularBoundary to write a polygon
    poly_bd_y.push_back(std::vector<double>());
    generateCircularBoundary(numEdges, scale_radius * tissueRadii[n], cx[n], cy[n], poly_bd_x[n], poly_bd_y[n]);

    for (i = cumNumCells; i < cumNumCells + numCellsInTissue[n]; i++) {
      cout << "i = " << i << "\t, dpos.size() = " << dpos.size() << ", cumNumCells = " << cumNumCells << "\t, numCellsInTissue[n] = " << numCellsInTissue[n] << '\n';
      double dpos_x = tissueRadii[n] * (2 * drand48() - 1) + cx[n], dpos_y = tissueRadii[n] * (2 * drand48() - 1) + cy[n];
      while (pow(dpos_x - cx[n], 2) + pow(dpos_y - cy[n], 2) > pow(tissueRadii[n], 2)) {
        dpos_x = tissueRadii[n] * (2 * drand48() - 1) + cx[n];
        dpos_y = tissueRadii[n] * (2 * drand48() - 1) + cy[n];
      }
      dpos.at(i * NDIM) = dpos_x;
      dpos.at(i * NDIM + 1) = dpos_y;
      ipStream << dpos[i * NDIM] << '\t' << dpos[i * NDIM + 1] << '\n';
      cellID[i] = n;
      cout << "cell " << i << " is in tissue number " << n << "with cell type = " << cellID[i] << '\n';
    }
    cout << "before adding to cumNumCells\n";
    cumNumCells += numCellsInTissue[n];
    cout << "just after adding to cumNumCells\n";
  }

  cout << "setting radii of SP disks\n";
  // set radii of SP disks
  for (ci = 0; ci < NCELLS; ci++) {
    xtra = 1.1;  // disks should have radius similar to the final particle radius, or could modify vrad[i] condition in wall calculation later
    drad.at(ci) = xtra * sqrt((2.0 * a0.at(ci)) / (nv.at(ci) * sin(2.0 * PI / nv.at(ci))));
    cout << "drad = " << drad[ci] << '\n';
  }

  // FIRE VARIABLES
  double P = 0;
  double fnorm = 0;
  double vnorm = 0;
  double alpha = alpha0;

  double dt0 = 1e-2;
  double dtmax = 10 * dt0;
  double dtmin = 1e-8 * dt0;

  int npPos = 0;
  int npNeg = 0;

  int fireit = 0;
  double fcheck = 10 * Ftol;

  // interaction variables
  double rij, sij, dtmp, ftmp, vftmp;
  double dr[NDIM];

  // initial step size
  dt = dt0;

  // loop until force relaxes
  while ((fcheck > Ftol) && fireit < itmax) {
    // FIRE step 1. Compute P
    P = 0.0;
    for (i = 0; i < cellDOF; i++)
      P += dv[i] * dF[i];

    // FIRE step 2. adjust simulation based on net motion of degrees of freedom
    if (P > 0) {
      // increase positive counter
      npPos++;

      // reset negative counter
      npNeg = 0;

      // alter simulation if enough positive steps have been taken
      if (npPos > NMIN) {
        // change time step
        if (dt * finc < dtmax)
          dt *= finc;

        // decrease alpha
        alpha *= falpha;
      }
    } else {
      // reset positive counter
      npPos = 0;

      // increase negative counter
      npNeg++;

      // check if simulation is stuck
      if (npNeg > NNEGMAX) {
        cerr << "	** ERROR: During initial FIRE minimization, P < 0 for too long, so ending." << endl;
        exit(1);
      }

      // take half step backwards, reset velocities
      for (i = 0; i < cellDOF; i++) {
        // take half step backwards
        dpos[i] -= 0.5 * dt * dv[i];

        // reset velocities
        dv[i] = 0.0;
      }

      // decrease time step if past initial delay
      if (fireit > NDELAY) {
        // decrease time step
        if (dt * fdec > dtmin)
          dt *= fdec;

        // reset alpha
        alpha = alpha0;
      }
    }

    // FIRE step 3. First VV update
    for (i = 0; i < cellDOF; i++)
      dv[i] += 0.5 * dt * dF[i];

    // FIRE step 4. adjust velocity magnitude
    fnorm = 0.0;
    vnorm = 0.0;
    for (i = 0; i < cellDOF; i++) {
      fnorm += dF[i] * dF[i];
      vnorm += dv[i] * dv[i];
    }
    fnorm = sqrt(fnorm);
    vnorm = sqrt(vnorm);
    if (fnorm > 0) {
      for (i = 0; i < cellDOF; i++)
        dv[i] = (1 - alpha) * dv[i] + alpha * (vnorm / fnorm) * dF[i];
    }

    // FIRE step 4. Second VV update
    for (i = 0; i < cellDOF; i++) {
      dpos[i] += dt * dv[i];
      dF[i] = 0.0;
    }

    // FIRE step 5. Update forces
    for (ci = 0; ci < NCELLS; ci++) {
      for (cj = ci + 1; cj < NCELLS; cj++) {
        // contact distance
        sij = drad[ci] + drad[cj];

        // true distance
        rij = 0.0;
        for (d = 0; d < NDIM; d++) {
          // get distance element
          dtmp = dpos[NDIM * cj + d] - dpos[NDIM * ci + d];
          if (pbc[d])
            dtmp -= L[d] * round(dtmp / L[d]);

          // add to true distance
          rij += dtmp * dtmp;

          // save in distance array
          dr[d] = dtmp;
        }
        rij = sqrt(rij);

        // check distances
        if (rij < sij) {
          // force magnitude
          ftmp = kc * (1.0 - (rij / sij)) / sij;

          // add to vectorial force
          for (d = 0; d < NDIM; d++) {
            vftmp = ftmp * (dr[d] / rij);
            dF[NDIM * ci + d] -= vftmp;
            dF[NDIM * cj + d] += vftmp;
          }
        }
      }
    }
    // FIRE step 4.1 Compute wall force
    for (int k = 0; k < poly_bd_x.size(); k++) {
      std::vector<double> poly_x = poly_bd_x[k];
      std::vector<double> poly_y = poly_bd_y[k];
      int n = poly_x.size();
      double distanceParticleWall, Rx, Ry, dw, K = 1;
      double bound_x1, bound_x2, bound_y1, bound_y2;
      // loop over boundary bars
      // loop over particles
      //  compute particle-boundary bar overlaps
      //  if overlap, Fx += K * dw * Rx/R, where K is a constant, dw = diameter/2 - R, Rx = x - px, R = sqrt(Rx^2 + Ry^2)
      for (int bound_i = 0; bound_i < n; bound_i++) {
        // use distanceLineAndPoint to get R, Rx, and Ry
        bound_x1 = poly_x[bound_i];
        bound_x2 = poly_x[(bound_i + 1) % n];
        bound_y1 = poly_y[bound_i];
        bound_y2 = poly_y[(bound_i + 1) % n];
        for (i = 0; i < cellDOF / NDIM; i++) {
          distanceParticleWall = distanceLinePointComponents(bound_x1, bound_y1, bound_x2, bound_y2, dpos[i * NDIM], dpos[i * NDIM + 1], Rx, Ry);
          dw = drad[i] - distanceParticleWall;
          if (distanceParticleWall <= drad[i]) {
            dF[i * NDIM] += K * dw * Rx / distanceParticleWall;
            dF[i * NDIM + 1] += K * dw * Ry / distanceParticleWall;
          }
        }
      }
    }
    // FIRE step 5. Final VV update
    for (i = 0; i < cellDOF; i++)
      dv[i] += 0.5 * dt * dF[i];

    // update forces to check
    fcheck = 0.0;
    for (i = 0; i < cellDOF; i++)
      fcheck += dF[i] * dF[i];
    fcheck = sqrt(fcheck / NCELLS);

    // print to console
    if (fireit % NSKIP == 0) {
      cout << endl
           << endl;
      cout << "===========================================" << endl;
      cout << "		I N I T I A L  S P 			" << endl;
      cout << " 	F I R E 						" << endl;
      cout << "		M I N I M I Z A T I O N 	" << endl;
      cout << "===========================================" << endl;
      cout << endl;
      cout << "	** fireit = " << fireit << endl;
      cout << "	** fcheck = " << fcheck << endl;
      cout << "	** fnorm = " << fnorm << endl;
      cout << "	** vnorm = " << vnorm << endl;
      cout << "	** dt = " << dt << endl;
      cout << "	** P = " << P << endl;
      cout << "	** Pdir = " << P / (fnorm * vnorm) << endl;
      cout << "	** alpha = " << alpha << endl;
    }

    // update iterate
    fireit++;
  }
  // check if FIRE converged
  if (fireit == itmax) {
    cout << "	** FIRE minimization did not converge, fireit = " << fireit << ", itmax = " << itmax << "; ending." << endl;
    exit(1);
  } else {
    cout << endl
         << endl;
    cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl;
    cout << "===========================================" << endl;
    cout << " 	F I R E 						" << endl;
    cout << "		M I N I M I Z A T I O N 	" << endl;
    cout << "	C O N V E R G E D! 				" << endl
         << endl;

    cout << "	(for initial disk minimization) " << endl;
    cout << "===========================================" << endl;
    cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl;
    cout << endl;
    cout << "	** fireit = " << fireit << endl;
    cout << "	** fcheck = " << fcheck << endl;
    cout << "	** vnorm = " << vnorm << endl;
    cout << "	** dt = " << dt << endl;
    cout << "	** P = " << P << endl;
    cout << "	** alpha = " << alpha << endl;
  }

  // initialize vertex positions based on cell centers
  for (ci = 0; ci < NCELLS; ci++) {
    ip2Stream << dpos[ci * NDIM] << '\t' << dpos[ci * NDIM + 1] << '\n';
    for (vi = 0; vi < nv.at(ci); vi++) {
      // get global vertex index
      gi = gindex(ci, vi);

      // length from center to vertex minus 1/2 the vertex radius to prevent overlaps
      dtmp = sqrt((2.0 * a0.at(ci)) / (nv.at(ci) * sin((2.0 * PI) / nv.at(ci)))) - r[gi] / 2.0;

      // set positions
      x.at(NDIM * gi) = dtmp * cos((2.0 * PI * vi) / nv.at(ci)) + dpos.at(NDIM * ci) + 1e-2 * l0[gi] * drand48();
      x.at(NDIM * gi + 1) = dtmp * sin((2.0 * PI * vi) / nv.at(ci)) + dpos.at(NDIM * ci + 1) + 1e-2 * l0[gi] * drand48();
    }
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
    cout << "compressing! it = " << it << '\n';
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

void cell::vertexCompress2Target2D_polygon(dpmMemFn forceCall, double Ftol, double dt0, double phi0Target, double dphi0) {
  // TEMPORARY JUST TO USE OVERLOADED PRINTCONFIGURATION2D
  // same as vertexCompress2Target2D, but with polygonal boundaries (affects packing fraction calculation, and expects forceCall to
  //  account for polygonal boundary forces
  // local variables
  int it = 0, itmax = 1e4;
  double phi0 = vertexPreferredPackingFraction2D_polygon();
  double scaleFactor = 1.0, P, Sxy;

  // loop while phi0 < phi0Target
  while (phi0 < phi0Target && it < itmax) {
    if (it >= itmax)
      cout << "inside vertexCompress2Target2D_polygon, reached maxit. exiting compression steps\n";
    // scale particle sizes
    scaleParticleSizes2D(scaleFactor);

    // update phi0
    phi0 = vertexPreferredPackingFraction2D_polygon();
    // relax configuration (pass member function force update)
    // make sure that forceCall is a force routine that includes a call to evaluatePolygonalWallForces
    vertexFIRE2D(forceCall, Ftol, dt0);

    // get scale factor
    scaleFactor = sqrt((phi0 + dphi0) / phi0);

    // get updated pressure
    P = 0.5 * (stress[0] + stress[1]);
    Sxy = stress[2];

    // print to console
    if (it % NSKIP == 0) {
      cout << " 	C O M P R E S S I O N 		" << endl;
      cout << "	** it 			= " << it << endl;
      cout << "	** phi0 curr	= " << phi0 << endl;
      if (phi0 + dphi0 < phi0Target)
        cout << "	** phi0 next 	= " << phi0 + dphi0 << endl;
      cout << "	** P 			= " << P << endl;
      cout << "	** Sxy 			= " << Sxy << endl;
      cout << "	** U 			= " << U << endl;
      printConfiguration2D();
      cout << endl
           << endl;
    }

    // update iterate
    it++;
  }
}

// used for running MD on neuralTube simulation
void cell::simulateDampedWithWalls(dpmMemFn forceCall, double B, double dt0, double duration, double printInterval, double pressureRate, double adhesionRate, bool wallsOn, bool leftOpen, bool bottomOpen, bool rightOpen, bool topOpen, double trueStrainRateX, double trueStrainRateY, double appliedUniaxialPressure) {
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
  double FT = 0, FB = 0, FL = 0, FR = 0;
  double oldFT = FT, oldFB = FB, oldFL = FL, oldFR = FR;
  std::vector<double> boundaryDisplacement(NDIM, 0.0);

  // set time step magnitude
  setdt(dt0);
  int NPRINTSKIP = printInterval / dt;

  // loop over time, print energy
  while (simclock - t0 < duration) {
    // set the dynamic boundary lengths to be equal to the current static boundary lengths. Both will be dynamic from here on out.
    XL[2] = L[0];
    XL[3] = L[1];

    appliedUniaxialPressure *= exp(pressureRate * dt);
    l1 *= exp(adhesionRate * dt);

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
      XL[0] += dt * VL[0] + 0.5 * dt * dt * FL - boundaryDisplacement[0] / 2;
      XL[1] += dt * VL[1] + 0.5 * dt * dt * FB - boundaryDisplacement[1] / 2;
      XL[2] += dt * VL[2] + 0.5 * dt * dt * FR + boundaryDisplacement[0] / 2;
      XL[3] += dt * VL[3] + 0.5 * dt * dt * FT + boundaryDisplacement[1] / 2;

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

// for testing numerical stability
void cell::vertexNVE(dpmMemFn forceCall, double T, double dt0, double duration, double printInterval) {
  // local variables
  int t, i;
  double K, t0 = simclock;
  double temp_simclock = simclock;

  // set time step magnitude
  setdt(dt0);
  int NPRINTSKIP = printInterval / dt;

  // initialize time keeper
  simclock = 0.0;

  // initialize velocities
  drawVelocities2D(T);

  // loop over time, print energy
  while (simclock - t0 < duration) {
    // VV VELOCITY UPDATE #1
    for (i = 0; i < vertDOF; i++)
      v[i] += 0.5 * dt * F[i];

    // VV POSITION UPDATE
    for (i = 0; i < vertDOF; i++) {
      // update position
      x[i] += dt * v[i];

      // recenter in box
      if (x[i] > L[i % NDIM] && pbc[i % NDIM])
        x[i] -= L[i % NDIM];
      else if (x[i] < 0 && pbc[i % NDIM])
        x[i] += L[i % NDIM];
    }

    // FORCE UPDATE
    CALL_MEMBER_FN(*this, forceCall)
    ();

    // VV VELOCITY UPDATE #2
    for (i = 0; i < vertDOF; i++)
      v[i] += 0.5 * F[i] * dt;

    // update sim clock
    simclock += dt;

    // print to console and file
    if (printInterval > dt) {
      if (int((simclock - t0) / dt) % NPRINTSKIP == 0 &&
          (simclock - temp_simclock) > printInterval / 2.0) {
        temp_simclock = simclock;
        // compute kinetic energy
        K = vertexKineticEnergy();

        // print to console
        cout << endl
             << endl;
        cout << "===============================" << endl;
        cout << "	D P M  						" << endl;
        cout << " 			 					" << endl;
        cout << "		N V E 					" << endl;
        cout << "===============================" << endl;
        cout << endl;
        cout << "	** simclock - t0 / duration	= " << simclock - t0
             << " / " << duration << endl;
        cout << " **  dt   = " << setprecision(12) << dt << endl;
        cout << "	** U 		= " << setprecision(12) << U << endl;
        cout << "	** K 		= " << setprecision(12) << K << endl;
        cout << "	** E 		= " << setprecision(12) << U + K
             << endl;

        // print to energy file
        cout << "** printing energy" << endl;
        enout << setw(w) << left << t;
        enout << setw(wnum) << left << simclock;
        enout << setw(wnum) << setprecision(12) << U;
        enout << setw(wnum) << setprecision(12) << K;
        enout << setw(wnum) << setprecision(12) << U + K;
        enout << endl;

        // print to configuration only if position file is open
        if (posout.is_open())
          printConfiguration2D();
      }
    }
  }
}

void cell::dampedVertexNVE(dpmMemFn forceCall, double B, double dt0, double duration, double printInterval) {
  // make sure velocities exist or are already initialized before calling this
  int i;
  double K, t0 = simclock;
  double temp_simclock = simclock;
  cout << "inside dampedVertexNVE in Cell class\n";
  cout << "vertDOF = " << vertDOF << '\n';
  cout << "x.size() = " << x.size() << '\n';
  cout << "l0.size() = " << l0.size() << '\n';
  cout << "a0.size() = " << a0.size() << '\n';
  cout << "cellID.size() = " << cellID.size() << '\n';

  // set time step magnitude
  setdt(dt0);
  int NPRINTSKIP = printInterval / dt;

  // loop over time, print energy
  while (simclock - t0 < duration) {
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
    ();

    // VV VELOCITY UPDATE #2
    for (i = 0; i < vertDOF; i++) {
      F[i] -= (B * v[i] + B * F_old[i] * dt / 2);
      F[i] /= (1 + B * dt / 2);
      v[i] += 0.5 * (F[i] + F_old[i]) * dt;
    }

    // update sim clock
    simclock += dt;

    // print to console and file
    if (int(printInterval) != 0) {
      if (int((simclock - t0) / dt) % NPRINTSKIP == 0 &&
          (simclock - temp_simclock) > printInterval / 2) {
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
        cout << "		N V E (DAMPED) 				"
                "	"
             << endl;
        cout << "===============================" << endl;
        cout << endl;
        cout << "	** simclock - t0 / duration	= " << simclock - t0
             << " / " << duration << endl;
        cout << " **  dt   = " << setprecision(12) << dt << endl;
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
          cout << "done printing in NVE\n";
        }

        if (tissueout.is_open()) {
          int maxCellID = *std::max_element(cellID.begin(), cellID.end());
          if (maxCellID > 0) {
            cout << "printing maxCellID boundary ID tissueMeasurements to file\n";
            // hardcoding boundary ID for now.
            takeTissueMeasurements(maxCellID);
          }
        }
      }
    }
  }
}

void cell::takeTissueMeasurements(int cellBoundaryType) {
  // compute boundary area - total cell area = extracellular space
  // write to file
  double boundaryArea = 0.0, totalCellArea = 0.0;
  int cellType;
  for (int ci = 0; ci < NCELLS; ci++) {
    cellType = cellID[ci];
    if (cellType == cellBoundaryType)
      boundaryArea = area(ci);
    else
      totalCellArea += area(ci);
  }
  if (boundaryArea == 0.0 || totalCellArea == 0.0) {
    cout << "error, could not find either boundary or cells in takeTissueMeasurements\n";
  }
  tissueout << simclock << '\t' << boundaryArea << '\t' << totalCellArea << '\n';
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
  // cout << L_left << '\t' << L_right << '\t' << L_bottom << '\t' << L_top << '\n';

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
    if (cellID.size() > 0)
      posout << setw(wnum) << left << cellID[ci];
    else
      posout << setw(wnum) << left << 1;
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