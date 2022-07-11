/*

        FUNCTION DEFINITIONS for epi (epithelial) class
  EPITHELIA FEATURE FLUCTUATING NUMBERS OF PARTICLES, ESPECIALLY SUBJECT
    TO EXPERIMENTAL PROCEDURES LIKE LASER ABLATION
  EPI2D CLASS INHERITS FROM DPM, WITH ALTERNATIVE ATTRACTION FORCES THAT AREN'T
  COMPATIBLE WITH THE CURRENT BENDING ENERGY FORCES EPI2D CLASS MUST TOLERATE
  PARTICLE DELETION AND ADDITION PROTOCOLS, SO SHOULD HAVE EASY ACCESS TO
  RESIZING ITS SIZE-DEPENDENT VECTORS

Warning to the new reader and myself - I use a few global-style variables defined inside of epi2D.h that are repeatedly modified throughout the code. This is bad engineering, and it would be good to instead pass needed information through function arguments and only through function arguments. This will be a focus of a major refactoring effort once I've managed to get a working wound simulation.
*/

#include "epi2D.h"

// namespace
using namespace Eigen;
using namespace std;

/******************************

        I N I T I A L -

        I Z A T I O N

*******************************/

/******************************

        E P I 2 D

        E D I T O R S  &

        U P D A T E S

*******************************/

double epi2D::meanl0() {
  int gi;
  double val;

  val = 0.0;
  for (gi = 0; gi < NVTOT; gi++)
    val += l0[gi];
  val /= NVTOT;

  return val;
}

double epi2D::meancalA0() {
  int gi, ci, vi;
  double calA0tmp, val;

  val = 0.0;
  for (ci = 0; ci < NCELLS; ci++) {
    calA0tmp = 0.0;
    for (vi = 0; vi < nv[ci]; vi++)
      calA0tmp += l0[szList[ci] + vi];
    calA0tmp *= calA0tmp;
    calA0tmp /= 4.0 * PI * a0[ci];
    val += calA0tmp;
  }
  val /= NCELLS;

  return val;
}

double epi2D::meankb() {
  int gi;
  double val;

  val = 0.0;
  for (gi = 0; gi < NVTOT; gi++)
    val += kbi[gi];
  val /= NVTOT;

  return val;
}

double epi2D::getPreferredPerimeter(int ci) {
  int gi0 = gindex(ci, 0);
  double prefPerimeter = 0.0;
  for (int gi = gi0; gi < gi0 + nv[ci]; gi++) {
    prefPerimeter += l0[gi];
  }
  return prefPerimeter;
}

double epi2D::vertDistNoPBC(int gi, int gj) {
  // not valid for PBC
  double dist = 0;
  dist = pow(x[gi * NDIM] - x[gj * NDIM], 2) +
         pow(x[gi * NDIM + 1] - x[gj * NDIM + 1], 2);
  return sqrt(dist);
}

void epi2D::directorDiffusion() {
  double r1, r2, grv;
  for (int ci = 0; ci < NCELLS; ci++) {
    // propagate diffusion of directors psi
    r1 = drand48();
    r2 = drand48();
    grv = sqrt(-2.0 * log(r1)) * cos(2.0 * PI * r2);
    // if flag is present, do not diffuse.
    psi[ci] += sqrt(dt * 2.0 * Dr0 * !flag[ci]) * grv;
    psi[ci] -= 2 * PI * round(psi[ci] / (2 * PI));
  }
}

std::vector<int> epi2D::regridSegment(int wVertIndex, double vrad) {
  // need to subtract indices from giConnectedToFlag if using this function in the future
  // return value is the indices of the deleted vertices for proper gi counting later
  /*if (vnn_label[wVertIndex] != 0) {
    if (simclock > 320 && simclock < 360) {
      cout << "vnn_label[wVertIndex] = " << vnn_label[wVertIndex] << ", wVertIndex = " << wVertIndex << '\n';
    }
    return {};
  }*/

  // regridSegment will adjust a cell's wound-adjacent # vertices and vertex positions to not overlap, deleting vertices and shifting as needed.
  //   to maintain constant energy before and after regridding, need to calculate a number c and use it
  //   to adjust preferred total perimeter after regridding?

  // wVertIndex = index of one of the wound cells. from here we'll calculate the
  // rest, take in a list of vertices that form a chain of line segments, and
  // return a linearly interpolated equivalent list of vertices
  int ci, vi;
  double d_T, newx, newy, T = 0, lerp, vdiam = vrad * 2, arc_length = 0;
  double cx, cy, rip1x, rip1y, rix, riy, lix, liy, li;
  cindices(ci, vi, wVertIndex);
  if (std::find(cellsLeavingPurseString.begin(), cellsLeavingPurseString.end(), ci) != cellsLeavingPurseString.end()) {
    return {};
  }
  if (simclock > 320 && simclock < 360) {
    cout << "inside regrid, we are working on cell ci = " << ci << '\n';
  }
  double L0tmp = getPreferredPerimeter(ci), L0new;
  double LCurrent = perimeter(ci);

  std::vector<int> deleteList;
  std::vector<double> result;
  std::vector<int> unorderedWoundList;

  // use sortedWoundIndices to find ci's wound-adjacent vertices
  int cj, vj, gi, gj;
  for (int j = 0; j < sortedWoundIndices.size(); j++) {
    cindices(cj, vj, sortedWoundIndices[j]);
    if (cj == ci) {  // collect all indices for cell ci that are wound-adjacent
      unorderedWoundList.push_back(sortedWoundIndices[j]);
    }
  }

  std::vector<int> orderedWVerts = unorderedWoundList;

  arc_length = rotateAndCalculateArcLength(ci, orderedWVerts);
  if (orderedWVerts.size() == 0) {
    cout << "why the heck does orderedWVerts.size() = 0?\n";
    cout << "because we've been asked to regrid a cell that has no vertices in the woundList\n";
    cout << "how is it that regridSegment has been called on a cell with no wound vertices?\n";
    cout << "sounds like i have cells ci in woundCells with no vertices in the wound list.\n";
    cout << "how to handle this?\n";
    cout << "also, at simclock = 259ish, wound cell 14 keeps attempting a regrid but never succeeds.";
    // to handle the last bit, let me first try printing out the number of vertices there.
    return {};
  }

  // nv_new = total length of segments divided by diameter of particle
  int nv_new = ceil(arc_length / vdiam);
  arc_length /= nv_new;  // arc_length could be infinity (nv_new could be zero)

  int numVertsToDelete = orderedWVerts.size() - nv_new;
  if (simclock > 320 && simclock < 360) {
    cout << "arc length during regrid = " << arc_length << " for wound cell " << ci << " simclock = " << simclock << '\n';
    cout << "position of first wound vertex is = " << x[orderedWVerts[0] * NDIM] << '\t' << x[orderedWVerts[0] * NDIM + 1] << '\n';
    cout << " orderedWVerts.size() = " << orderedWVerts.size() << '\n';
    cout << " numVertsToDelete = " << numVertsToDelete << '\n';
    cout << " nv_new = " << nv_new << '\n';
  }

  if (numVertsToDelete <= 0)
    return {};

  if (nv_new == 1 && orderedWVerts.size() > 0) {
    double avgx = 0.0, avgy = 0.0;
    // special case, just remove all but one vertex, placed at the center of all the wounded vertices
    for (int i = 0; i < orderedWVerts.size(); i++) {
      gi = orderedWVerts[i];
      deleteList.push_back(gi);
      avgx += x[gi * NDIM];
      avgy += x[gi * NDIM + 1];
      cout << "in nv_new=1 special case. x[gi*NDIM], x[gi*NDIM + 1] = " << x[gi * NDIM] << '\t' << x[gi * NDIM + 1] << '\n';
    }
    cout << "special case nv_new = 1, deleting all vertices but one\n";
    for (auto i : szList)
      cout << "szList : " << i << '\n';
    deleteList.pop_back();  // save one vertex
    avgx /= orderedWVerts.size();
    avgy /= orderedWVerts.size();
    x[gi * NDIM] = avgx;
    x[gi * NDIM + 1] = avgy;
    cout << "modifying position of vertex " << gi << ", to " << x[gi * NDIM] << '\t' << x[gi * NDIM + 1] << ", simclock = " << simclock << '\n';
    // since I've deleted all of its wound vertices, I need to force this cell to be excommunicated
    cellsLeavingPurseString.push_back(ci);
  } else {
    for (int i = 0; i < int(orderedWVerts.size()) - 1; i++) {
      gi = orderedWVerts[i];
      gj = orderedWVerts[i + 1];
      d_T = vertDistNoPBC(gi, gj) / arc_length;
      while (T + d_T >= result.size() / 2 && result.size() / 2 < nv_new) {
        lerp = (result.size() / 2 - T) / d_T;
        newx = (1 - lerp) * x[gi * NDIM] + lerp * x[gj * NDIM];
        newy = (1 - lerp) * x[gi * NDIM + 1] + lerp * x[gj * NDIM + 1];
        result.push_back(newx);
        result.push_back(newy);
      }
      T += d_T;
    }
    if (simclock > 320 && simclock < 360) {
      if (orderedWVerts.size() < 3) {
        cout << "orderedWVerts.size = " << orderedWVerts.size() << " for cell " << ci << '\n';
        cout << "result.size() = " << result.size() << '\n';
      }
    }
    if (numVertsToDelete >= 1) {
      // int(orderedWVerts.size()) - 1 : don't delete the last vertex, it's not accounted for
      int maxIt = int(orderedWVerts.size()) - 1;
      for (int i = 1; i < int(orderedWVerts.size()) - 1; i++) {
        gi = orderedWVerts[i];
        if (i < result.size() / 2 && maxIt - 1 >= result.size() / 2) {  // if i is small enough, and if its max value ensures deleteList will be populated
          x[gi * NDIM] = result[2 * i];
          x[gi * NDIM + 1] = result[2 * i + 1];
          cout << "modifying position of vertex " << gi << ", to " << x[gi * NDIM] << '\t' << x[gi * NDIM + 1] << ", simclock = " << simclock << '\n';
        } else if (i >= result.size() / 2)
          deleteList.push_back(gi);
      }
      if (!deleteList.empty()) {
        cout << "numVertsToDelete =  " << numVertsToDelete << "!\n";
        cout << "rotated orderedWVerts looks like (before deletion): \n";
        for (auto i : orderedWVerts) {
          cout << "gi = " << i << "\t, ip1[gi] = " << ip1[i] << "\t, im1[gi] = " << im1[i] << '\n';
        }
      }
    }
  }

  if (!deleteList.empty()) {
    for (auto i : deleteList)
      cout << "vertex to delete: " << i << '\t' << x[NDIM * i] << '\t' << x[NDIM * i + 1] << '\n';
    cout << "wVertIndex = " << wVertIndex << '\n';
    if (wVertIndex == 338) {
      for (auto j : szList)
        cout << "szList : " << j << '\n';
      for (auto j : unorderedWoundList)
        cout << "unorderedWoundList : " << j << '\n';
      cout << "ci = " << ci << '\n';
    }
    for (auto i : orderedWVerts)
      cout << "orderedWVerts : " << i << '\t' << x[NDIM * i] << '\t' << x[NDIM * i + 1] << '\n';

    deleteVertex(deleteList);

    for (auto i : deleteList) {
      orderedWVerts.erase(std::remove(orderedWVerts.begin(), orderedWVerts.end(), i), orderedWVerts.end());
      for (int j = 0; j < orderedWVerts.size(); j++)
        if (orderedWVerts[j] > i)
          orderedWVerts[j]--;
    }

    // after regridding and deleting, have all wound vertices set their preferred lengths to their current length. gives rigidity.

    for (int i = 0; i < orderedWVerts.size(); i++) {
      gi = orderedWVerts[i];
      cout << "gi, l0(before), l0(after) = " << gi << '\t' << l0[gi] << "\t" << vertDistNoPBC(gi, ip1[gi]) << '\n';
      l0[gi] = vertDistNoPBC(gi, ip1[gi]);
      cout << "orderedWVerts after regrid: " << gi << ", ip1[gi] = " << ip1[gi] << '\t' << x[NDIM * gi] << '\t' << x[NDIM * gi + 1] << '\n';
    }
  }
  return deleteList;
}

void epi2D::resetActiveEnergy() {
  U_ps = 0.0;
  U_crawling = 0.0;
}

void epi2D::repulsiveForceUpdateWithWalls() {
  double a, b, c, d;
  resetForcesAndEnergy();
  shapeForces2D();
  vertexRepulsiveForces2D();
  wallForces(false, false, false, false, a, b, c, d);
}

void epi2D::repulsiveForceWithCircularApertureWall() {
  double radius = 2.0;
  resetForcesAndEnergy();
  shapeForces2D();
  vertexRepulsiveForces2D();
  circularApertureForces(radius);
}

void epi2D::repulsiveForceUpdateWithPolyWall() {
  resetForcesAndEnergy();
  shapeForces2D();
  vertexRepulsiveForces2D();
  for (int i = 0; i < poly_bd_x.size(); i++){
    evaluatePolygonalWallForces(poly_bd_x[i], poly_bd_y[i]);
  }
}

void epi2D::attractiveForceUpdateWithPolyWall() {
  resetForcesAndEnergy();
  shapeForces2D();
  vertexAttractiveForces2D_2();
  for (int i = 0; i < poly_bd_x.size(); i++){
    evaluatePolygonalWallForces(poly_bd_x[i], poly_bd_y[i]);
  }
}

void epi2D::vertexAttractiveForces2D_2() {
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

              // add to contacts
              for (int i = 0; i < vnn[gi].size(); i++) {
                if (ci == cj)
                  break;

                if (vnn[gi][i] < 0) {
                  vnn[gi][i] = gj;  // set the first unused array element to gj, in gi's neighbor list

                  for (int j = 0; j < vnn[gj].size(); j++) {
                    if (vnn[gj][j] < 0) {
                      vnn[gj][j] = gi;  // set the first unused array element to gi, in gj's neighbor list
                      break;
                    }
                  }

                  break;
                }
              }
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
                  for (int i = 0; i < vnn[gi].size(); i++) {
                    if (vnn[gi][i] < 0) {
                      vnn[gi][i] = gj;  // set the first unused array element to gj, in gi's neighbor list

                      for (int j = 0; j < vnn[gj].size(); j++) {
                        if (vnn[gj][j] < 0) {
                          vnn[gj][j] = gi;  // set the first unused array element to gi, in gj's neighbor list
                          break;
                        }
                      }

                      break;
                    }
                  }
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

void epi2D::attractiveForceUpdate_2() {
  resetForcesAndEnergy();
  shapeForces2D();
  vertexAttractiveForces2D_2();
}

void epi2D::substrateadhesionAttractiveForceUpdate() {
  // compute forces for shape, attractive, and substrate adhesion contributions
  int gi = 0, argmin, flagcount = 0;
  int numVerticesAttracted = 1;
  double dx, dy;
  double fx, fy;

  // reset forces, then get shape and attractive forces.
  attractiveForceUpdate_2();
  resetActiveEnergy();
  // directorDiffusion();
  updateSubstrateSprings();

  for (auto ci : initialWoundCellIndices) {
    std::vector<double> distSqVertToFlag(nv[0], 0.0);
    // check for protrusions (unimpeded flag toss, and flag location inside box)
    if (flag[ci]) {
      // find nearest vertices
      /*for (int vi = 0; vi < nv[ci]; vi++) {
        gi = gindex(ci, vi);
        distSqVertToFlag[vi] = pow(x[gi * NDIM] - flagPos[ci][0], 2) +
                               pow(x[gi * NDIM + 1] - flagPos[ci][1], 2);
      }
      argmin = std::distance(
          distSqVertToFlag.begin(),
          std::min_element(distSqVertToFlag.begin(), distSqVertToFlag.end()));
      gi = gindex(ci, argmin);*/
      gi = giConnectedToFlag[ci];

      // evaluate force for spring-vertex interaction between nearest vertex and
      // flag position, i.e. force due to crawling
      dx = x[gi * NDIM] - flagPos[ci][0] - restLengthLPx[ci];
      dy = x[gi * NDIM + 1] - flagPos[ci][1] - restLengthLPy[ci];
      fx = -k_LP * dx;
      fy = -k_LP * dy;
      F[gi * NDIM] += fx;
      F[gi * NDIM + 1] += fy;

      fieldStress[gi][0] += (x[gi * NDIM] - flagPos[ci][0]) / 2 * fx;
      fieldStress[gi][1] += (x[gi * NDIM + 1] - flagPos[ci][1]) / 2 * fy;
      fieldStress[gi][2] += 0.5 * ((x[gi * NDIM] - flagPos[ci][0]) / 2 * fy + (x[gi * NDIM + 1] - flagPos[ci][1]) / 2 * fx);

      flagcount++;
      // strain the rest lengths towards zero by exp(-t/tau_LP)
      restLengthLPx[ci] *= exp(-dt / tau_LP);
      restLengthLPy[ci] *= exp(-dt / tau_LP);
      timeElapsedSinceFlagPlanted[ci] += dt;  // unused for now
      // update energy due to crawling
      U += 0.5 * k_LP * (dx * dx + dy * dy);
      U_crawling += 0.5 * k_LP * (dx * dx + dy * dy);
    }
  }
}

/******************************

        EPI

                P R O T O C O L S

*******************************/

void epi2D::ageCellAreas(double areaScaleFactor) {
  int gi, ci, vi;
  double lengthScaleFactor = sqrt(areaScaleFactor);

  // loop over cells, scale
  for (ci = 0; ci < NCELLS; ci++) {
    // scale preferred area
    a0[ci] *= areaScaleFactor;

    // first global index for ci
    gi = szList[ci];

    for (vi = 0; vi < nv[ci]; vi++) {
      // scale vertex radii
      // r[gi + vi] *= lengthScaleFactor;
      // l0[gi + vi] *= lengthScaleFactor;
    }
  }
}

void epi2D::ageCellPerimeters(double shapeRelaxationRate, double dt) {
  int gi;
  double current_perimeter, current_preferred_perimeter, new_preferred_perimeter;
  double shapeRelaxFactor = exp(-shapeRelaxationRate * dt);
  double scaleFactor = 1;
  for (int ci = 0; ci < NCELLS; ci++) {
    current_perimeter = perimeter(ci);
    current_preferred_perimeter = 0.0;
    for (int vi = 0; vi < nv[ci]; vi++) {
      gi = gindex(ci, vi);
      current_preferred_perimeter += l0[gi];
    }
    new_preferred_perimeter = current_perimeter * (1 - shapeRelaxFactor) + current_preferred_perimeter * shapeRelaxFactor;
    scaleFactor = new_preferred_perimeter / current_preferred_perimeter;
    for (int vi = 0; vi < nv[ci]; vi++) {
      gi = gindex(ci, vi);
      l0[gi] *= scaleFactor;
    }
  }
}

void epi2D::expandBoxAndCenterParticles(double boxLengthScaleFactor,
                                        double boxLengthScale) {
  if (boxLengthScaleFactor >= 1) {
    for (int gi = 0; gi < NVTOT; gi++) {
      x[gi * NDIM] += L[0] * (boxLengthScaleFactor - 1) / 2;
      x[gi * NDIM + 1] += L[1] * (boxLengthScaleFactor - 1) / 2;
    }
    L[0] *= boxLengthScaleFactor;
    L[1] *= boxLengthScaleFactor;

  } else {
    cout << "canceling expandBoxAndCenterParticles because of invalid scale "
            "factor\n";
  }
  initializeNeighborLinkedList2D(boxLengthScale);
}

// scale particle sizes by specified ratio in x and y directions
// compression in X and Y corresponds to isotropic tensile loading
void epi2D::tensileLoading(double scaleFactorX, double scaleFactorY) {
  // local variables
  int gi, ci, vi, xind, yind;
  double xi, yi, cx, cy, dx, dy;

  double totalScaleFactor = scaleFactorX * scaleFactorY;

  // loop over cells, scale
  for (ci = 0; ci < NCELLS; ci++) {
    // scale preferred area
    a0[ci] *= totalScaleFactor;

    // first global index for ci
    gi = szList[ci];

    com2D(ci, cx, cy);

    for (vi = 0; vi < nv[ci]; vi++) {
      // x and y inds
      xind = NDIM * (gi + vi);
      yind = xind + 1;

      // closest relative position
      dx = x[xind] - cx;
      if (pbc[0])
        dx -= L[0] * round(dx / L[0]);

      dy = x[yind] - cy;
      if (pbc[1])
        dy -= L[1] * round(dy / L[1]);

      // update vertex positions
      x[xind] += (scaleFactorX - 1.0) * dx;
      x[yind] += (scaleFactorY - 1.0) * dy;

      // scale vertex radii
      r[gi + vi] *= sqrt(totalScaleFactor);
      l0[gi + vi] *= sqrt(totalScaleFactor);
    }
  }
}

void epi2D::updateSubstrateSprings() {
  // refresh spring contacts based on protrusion lifetime tau_LP and protrusion length maxProtrusionLength
  // a cell will protrude as far out as maxProtrusionLength if unobstructed, otherwise it will choose a shorter length
  // check to see if enough time has passed for us to update springs again
  double refreshInterval = 1.0;
  double woundAreaCutOff = 30 * PI * r[0] * r[0];  // 30 particle areas is roughly the size of a wound where each vertex is fully adhered

  // refreshInterval = tau_LP;
  if (simclock - previousUpdateSimclock > refreshInterval) {
    double cx, cy, gj, minFlagDistance, flagDistance;

    if (giConnectedToFlag.size() != NCELLS)
      giConnectedToFlag.resize(NCELLS);

    // flagDistance = 1.5 * sqrt(a0[0] / PI);
    bool cancelFlagToss;
    std::vector<std::vector<double>> center(NCELLS,
                                            std::vector<double>(2, 0.0));
    for (int ci = 0; ci < NCELLS; ci++) {
      // find cell centers
      com2D(ci, cx, cy);
      center[ci][0] = cx;
      center[ci][1] = cy;
    }
    // sweep through cells. If no flag, try to plant one. Else if flag, try to
    // dissociate.
    for (auto ci : initialWoundCellIndices) {
      cancelFlagToss = false;
      if (!flag[ci]) {
        // pick a direction, throw a flag in that direction
        int gi;  // nearest vertex to the flag
        // minFlagDistance = getDistanceToVertexAtAnglePsi(ci, psi[ci], center[ci][0], center[ci][1], gi);

        // randomly chooses an unadhered vertex, updates psi[ci] using that vertex, and gives the distance to that vertex.
        minFlagDistance = getDistanceToRandomUnadheredVertex(ci, center[ci][0], center[ci][1], gi);

        if (minFlagDistance < 0 && woundArea > woundAreaCutOff)  // failed to find a valid vertex to throw a flag from
          continue;
        else if (minFlagDistance < 0 && woundArea <= woundAreaCutOff)  // wound is small, let every cell crawl towards its last known polarity to close the wound
          minFlagDistance = getDistanceToVertexAtAnglePsi(ci, psi[ci], center[ci][0], center[ci][1], gi);

        // flagDistance += 3 * 2 * r[gi];
        double fractionOfDiameter = 4.0;  // try protruding in increments of diameter / fractionOfDiameter
        for (int i = 0; i < fractionOfDiameter * floor(maxProtrusionLength) - 1 && !flag[ci]; i++) {
          cancelFlagToss = false;
          flagDistance = minFlagDistance + (maxProtrusionLength - i / fractionOfDiameter) * 2.0 * r[gi];
          // start with larger flag distances. if those are blocked, then loop tries smaller values
          flagPos[ci][0] = center[ci][0] + flagDistance * cos(psi[ci]);
          flagPos[ci][1] = center[ci][1] + flagDistance * sin(psi[ci]);
          // loop over cells near enough to block flag
          for (int cj = 0; cj < NCELLS; cj++) {
            if (cancelFlagToss == true)
              continue;
            if (ci != cj) {
              // check if centers are near enough to possibly block flag
              if (pow(center[ci][0] - center[cj][0], 2) +
                      pow(center[ci][1] - center[cj][1], 2) <
                  3 * pow(flagDistance, 2)) {
                // check if vertices actually block flag (flag is allowed to penetrate the vertex but not pass it, because basal protrusions have more freedom than 3D membranes)
                for (int vj = 0; vj < nv[cj]; vj++) {
                  gj = gindex(cj, vj);
                  if (distanceLineAndPoint(center[ci][0], center[ci][1],
                                           flagPos[ci][0], flagPos[ci][1],
                                           x[gj * NDIM],
                                           x[gj * NDIM + 1]) < r[gj]) {
                    // yes, flag has been blocked at this length. Try new flag length.
                    cancelFlagToss = true;
                    break;
                  }
                }
              }
            }
          }
          if (cancelFlagToss == false) {  // flag planted successfully, move on to next cell's flag toss attempt
            flag[ci] = true;
            timeElapsedSinceFlagPlanted[ci] = 0.0;
            // restLength can be negative, indicating relative direction to x-x_flag
            restLengthLPx[ci] = x[gi * NDIM] - flagPos[ci][0];
            restLengthLPy[ci] = x[gi * NDIM + 1] - flagPos[ci][1];
            giConnectedToFlag[ci] = gi;
            // restLengthLPx[ci] = 0;
            // restLengthLPy[ci] = 0;
            // cout << "flagDistance = " << flagDistance << "\t, 2*r[gi] = " << 2*r[gi] << '\n';
            // cout << "i = " << i << "\t, minFlagDistance = " << minFlagDistance << "\t, psi[ci] = " << psi[ci] << '\n';
          }
        }
      } else if (flag[ci]) {
        // dissociation rate determines if existing flag is destroyed this step
        // 10% chance of dissociating per tau_LP, but also automatically dissociate if spring has shrunk below some threshhold ~ vdiam/10
        if (drand48() < 0.0 || fabs(restLengthLPx[ci]) + fabs(restLengthLPx[ci]) < r[0] / 10.0) {
          flag[ci] = false;
        }
      } else
        cerr << "Error: no flag bool value found\n";
    }
    // update time tracker
    previousUpdateSimclock = simclock;
  }
}

void epi2D::dampedNVE2D(dpmMemFn forceCall, double B, double dt0, double duration, double printInterval) {
  // make sure velocities exist or are already initialized before calling this
  int i;
  double K, t0 = simclock;
  double temp_simclock = simclock;

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
    boundaries();  // build vnn for void detection
    std::vector<double> F_old = F;
    CALL_MEMBER_FN(*this, forceCall)
    ();

    // VV VELOCITY UPDATE #2
    for (i = 0; i < vertDOF; i++) {
      F[i] -= (B * v[i] + B * F_old[i] * dt / 2);
      F[i] /= (1 + B * dt / 2);
      v[i] += 0.5 * (F[i] + F_old[i]) * dt;
      // cout << "force at simclock = " << simclock << " = ," << F[i] << '\n';
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
          enout << setw(wnum) << left << L[0] / initialLx - 1;
          enout << setw(wnum) << setprecision(12) << U;
          enout << setw(wnum) << setprecision(12) << K;
          enout << setw(wnum) << setprecision(12) << U_ps;
          enout << setw(wnum) << setprecision(12) << U_crawling;
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
          stressout << setw(wnum) << left << L[0] / initialLx - 1;
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
          printBoundaries();
          cout << "done printing in NVE\n";
        }

        cerr << "Number of polarization deflections: " << polarizationCounter
             << '\n';
      }
    }
  }
}

// currently, dampedNP0 is my way of coding a simulation without boundaries, i.e. 0 pressure simulation
void epi2D::dampedNP0(dpmMemFn forceCall, double B, double dt0, double duration, double printInterval, int wallsOn) {
  // make sure velocities exist or are already initialized before calling this
  // assuming zero temperature - ignore thermostat (not implemented)
  // allow box lengths to move as a dynamical variable - rudimentary barostat,
  // doesn't matter for non-equilibrium anyway

  // need to erase the lines that recenter stuff? also needs to be set as an
  // option in e.g. the force calls? in dpm? and in epi2D?
  int i;
  double K, t0 = simclock;
  double temp_simclock = simclock;
  double FT, FB, FL, FR;
  double oldFT, oldFB, oldFL, oldFR;

  // set time step magnitude
  setdt(dt0);
  int NPRINTSKIP = printInterval / dt;
  int nthLargestCluster = 2;

  // set purse-string spring breaking distance
  double vertDiameterSq = pow(2 * r[0], 2);

  initialRadius = r;
  initiall0 = l0;
  initialPreferredPerimeter = 0;
  for (int i = 0; i < nv[0]; i++) {
    initialPreferredPerimeter += l0[i];
  }

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
    boundaries();  // build vnn for void detection, only intracellularly. The rest (due to adhesion) is done in forceCall.
    std::vector<double> F_old = F;

    CALL_MEMBER_FN(*this, forceCall)
    ();  // calls main force routine

    if (simclock - t0 > 10) {
      //.clear();
      initialPreferredPerimeter = 0;
      for (int i = 0; i < nv[0]; i++) {
        initialPreferredPerimeter += l0[i];
      }
      if (psContacts.size() == 0 && simclock - t0 < 100) {
        getWoundVertices(nthLargestCluster);
        woundCenterX = 0;
        woundCenterY = 0;
        for (auto gi : sortedWoundIndices) {
          woundCenterX += x[gi * NDIM];
          woundCenterY += x[gi * NDIM + 1];
        }
        woundCenterX /= sortedWoundIndices.size();
        woundCenterY /= sortedWoundIndices.size();
        initializePurseStringVariables();
        cout << "initialized purseString variables!\n";
      }

      // get max and min of x coords of purse-string; if max-min is near zero, then purse-string should be dissolved
      double max_ps = x_ps[0], min_ps = x_ps[0];
      double min_allowed_ps_length = r[0];
      // double min_allowed_ps_length = 0.0;
      for (int psi = 0; psi < x_ps.size() / 2; psi += 2) {
        if (x_ps[psi] < min_ps)
          min_ps = x_ps[psi];
        if (x_ps[psi] > max_ps)
          max_ps = x_ps[psi];
      }

      if (max_ps - min_ps > min_allowed_ps_length) {
        purseStringContraction(B);
      } else {  // purse-string can't shrink anymore
        // cout << "purse-string is too small, disassemble it\n";
        // k_ps = 0.0;
        // deltaSq = 0.0;
        if (isPurseStringDoneShrinking == false) {  // if purse-string can't shrink anymore, set this flag to true and stop straining the pursestring
          cout << "disassembling purse-string, lengths are too short to shrink further \n";
          isPurseStringDoneShrinking = true;
          strainRate_ps = 0.0;
        }
        purseStringContraction(B);
      }

      // vout << simclock << '\t' << woundArea << '\n';
      ageCellPerimeters(shapeRelaxationRate, dt);
      if (int(simclock / dt) % 500 == 0) {
        cout << "woundCenterX, Y before calculating area = " << woundCenterX << '\t' << woundCenterY << '\n';
        woundArea = calculateWoundArea(woundCenterX, woundCenterY);
        vout << simclock << '\t' << woundArea << '\n';
      }
    }

    // VV VELOCITY UPDATE #2
    for (i = 0; i < vertDOF; i++) {
      F[i] -= (B * v[i] + B * F_old[i] * dt / 2);
      F[i] /= (1 + B * dt / 2);
      v[i] += 0.5 * (F[i] + F_old[i]) * dt;
    }

    if (wallsOn > 0) {
      // VV position update (walls)
      L[0] += dt * VL[0] + 0.5 * dt * dt * FT;
      L[1] += dt * VL[1] + 0.5 * dt * dt * FR;

      // VV force update (walls)
      oldFT = FT;
      oldFB = FB;
      oldFL = FL;
      oldFR = FR;
      if (wallsOn == 1)
        wallForces(false, false, true, true, FL, FB, FR, FT);
      if (wallsOn == 2)
        wallForces(false, false, false, false, FL, FB, FR, FT, wallsOn);  // all walls fixed, and force_multiplier option makes them super sticky
      FT -= (B * VL[1] + B * oldFT * dt / 2);
      FT /= (1 + B * dt / 2);
      FB -= (B * VL[1] + B * oldFB * dt / 2);
      FB /= (1 + B * dt / 2);
      FL -= (B * VL[0] + B * oldFL * dt / 2);
      FL /= (1 + B * dt / 2);
      FR -= (B * VL[0] + B * oldFR * dt / 2);
      FR /= (1 + B * dt / 2);
      VL[0] += 0.5 * (FR + oldFR) * dt;
      VL[1] += 0.5 * (FT + oldFT) * dt;
    }

    // update sim clock
    simclock += dt;
    // print to console and file
    if (int(printInterval) != 0) {
      /*if (int((simclock - t0) / dt) % NPRINTSKIP == 0 &&
          (simclock - temp_simclock) > printInterval / 2) {*/
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
          enout << setw(wnum) << left << L[0] / initialLx - 1;
          enout << setw(wnum) << setprecision(12) << U;
          enout << setw(wnum) << setprecision(12) << K;
          enout << setw(wnum) << setprecision(12) << U_ps;
          enout << setw(wnum) << setprecision(12) << U_crawling;
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
          stressout << setw(wnum) << left << L[0] / initialLx - 1;
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
          int nthLargestCluster = 2;
          printConfiguration2D();
          printBoundaries(nthLargestCluster);
          cerr << "done printing in NP0\n";
        }

        cerr << "Number of polarization deflections: " << polarizationCounter
             << '\n';
      }
    }
  }
}

void epi2D::circularApertureForces(double radius) {
  // if any cells have their centers inside the aperture, force them away from
  // the aperture
  double cx, cy, gi, s, cut, shell, rx, ry, disti, overlap, theta;
  for (int ci = 0; ci < NCELLS; ci++) {
    com2D(ci, cx, cy);
    cx -= L[0] / 2;
    cy -= L[1] / 2;
    if (cx * cx + cy * cy < radius * radius) {
      for (int vi = 0; vi < nv[ci]; vi++) {
        gi = gindex(ci, vi);
        F[gi * NDIM] += -0.5 * (L[0] / 2 - cx);
        F[gi * NDIM + 1] += -0.5 * (L[1] / 2 - cy);
      }
    }
  }

  // aperture has volume exclusion with any vertices.
  for (int gi = 0; gi < NVTOT; gi++) {
    s = 2 * r[gi];
    cut = (1.0 + l1) * s;
    shell = (1.0 + l2) * s;
    rx = x[gi * NDIM];
    ry = x[gi * NDIM + 1];
    rx -= L[0] / 2;
    ry -= L[1] / 2;
    disti = sqrt(rx * rx + ry * ry);
    overlap = disti - fabs(radius - r[gi]);
    if (disti < 0) {
      theta = atan2(rx, ry);
      F[gi * NDIM] += 0.5 * overlap / s * cos(theta);
      F[gi * NDIM + 1] += 0.5 * overlap / s * sin(theta);
    }
  }
}

void epi2D::vertexCompress2Target2D_polygon(dpmMemFn forceCall, double Ftol, double dt0, double phi0Target, double dphi0) {
  // TEMPORARY JUST TO USE PRINTCONFIGURATION2D in epi2D overload
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

void epi2D::wallForces(bool left, bool bottom, bool right, bool top, double& forceLeft, double& forceBottom, double& forceRight, double& forceTop, int forceOption) {
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
  double force_multiplier = 1;
  if (forceOption == 2)
    force_multiplier = 10.0;
  double kc_temp = kc * force_multiplier;
  double shell, cut, s, distLower, distUpper, scaledDist, ftmp, f;
  double kint = (kc_temp * l1) / (l2 - l1);
  forceTop = 0;
  forceBottom = 0;
  forceLeft = 0;
  forceRight = 0;

  // if any cells have their centers outside of the box, force them towards the
  // center
  for (int ci = 0; ci < NCELLS; ci++) {
    com2D(ci, cx, cy);
    if (cx < 0 || cx > L[0]) {
      for (int vi = 0; vi < nv[ci]; vi++) {
        gi = gindex(ci, vi);
        F[gi * NDIM] += -0.5 * (x[gi * NDIM] - L[0] / 2);
      }
    }
    if (cy < 0 || cy > L[1]) {
      for (int vi = 0; vi < nv[ci]; vi++) {
        gi = gindex(ci, vi);
        F[gi * NDIM + 1] += -0.5 * (x[gi * NDIM + 1] - L[1] / 2);
      }
    }
  }

  for (int i = 0; i < vertDOF; i++) {
    vi = i / 2;
    boxL = L[i % NDIM];

    s = 2 * r[vi];
    cut = (1.0 + l1) * s;
    shell = (1.0 + l2) * s;
    distLower = fabs(x[i]);
    distUpper = fabs(boxL - x[i]);

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
}

void epi2D::zeroMomentum() {
  // subtract off any linear momentum by reducing momentum of each particle by
  // momentum/#vertices
  double v_cm_x = 0.0, v_cm_y = 0.0;

  // accumulate global vertex velocities
  for (int gi = 1; gi < NVTOT; gi++) {
    v_cm_x += v[NDIM * gi];
    v_cm_y += v[NDIM * gi + 1];
  }

  // subtract off center of mass drift
  for (int gi = 0; gi < NVTOT; gi++) {
    v[NDIM * gi] -= v_cm_x / NVTOT;
    v[NDIM * gi + 1] -= v_cm_y / NVTOT;
  }
}

void epi2D::scaleBoxSize(double boxLengthScale, double scaleFactorX, double scaleFactorY) {
  // scale box size while accounting for periodic boundary conditions
  // takes in the original coordinates, alters the box length, and updates the
  // coordinates according to the new boundary conditions local variables
  int gi, ci, vi, xind, yind;
  double xi, yi, cx, cy, dx, dy;
  double newLx = scaleFactorX * L[0];
  double newLy = scaleFactorY * L[1];
  // loop over cells, scale
  for (ci = 0; ci < NCELLS; ci++) {
    // first global index for ci
    gi = szList[ci];

    com2D(ci, cx, cy);

    for (vi = 0; vi < nv[ci]; vi++) {
      // x and y inds
      xind = NDIM * (gi + vi);
      yind = xind + 1;

      // closest relative position
      dx = x[xind] - cx;
      if (pbc[0])
        dx -= L[0] * round(dx / L[0]);

      dy = x[yind] - cy;
      if (pbc[1])
        dy -= L[1] * round(dy / L[1]);

      // compute vertex positions relative to center of box
      x[xind] = (cx - L[0] / 2 + dx);
      x[yind] = (cy - L[1] / 2 + dy);

      // wrap coordinates relative to center of box
      if (pbc[0])
        x[xind] -= L[0] * round(x[xind] / L[0]);
      if (pbc[1])
        x[yind] -= L[1] * round(x[yind] / L[1]);

      // scale coordinates relative to center of box
      x[xind] *= scaleFactorX;
      x[yind] *= scaleFactorY;

      // put coordinates back to bottom left origin convention
      x[xind] += newLx / 2;
      x[yind] += newLy / 2;

      // wrap coordinates relative to bottom left corner, with new scaled box
      // lengths
      x[xind] -= newLx * round(x[xind] / newLx);
      x[yind] -= newLy * round(x[yind] / newLy);
    }
  }
  L[0] = newLx;
  L[1] = newLy;
  initializeNeighborLinkedList2D(boxLengthScale);
}

int epi2D::getIndexOfCellLocatedHere(double xLoc, double yLoc) {
  // return cell index of cell nearest specified (xLoc,yLoc) coordinate.
  // loop through global indices, compute all centers of mass (stored in
  // simulation as unbounded coordinates) distance (squared) of centers of mass
  // to origin
  vector<double> distanceSq = vector<double>(NCELLS, 1e7);
  double dx = 0.0, dy = 0.0, cx, cy;
  for (int ci = 0; ci < NCELLS; ci++) {
    // first global index for ci
    int gi = szList[ci];

    com2D(ci, cx, cy);  // cx, cy use coordinates such that bottom left corner of
                        // sim box = origin

    distanceSq[ci] =
        pow(cx - (L[0] / 2 + xLoc), 2) + pow(cy - (L[1] / 2 + yLoc), 2);
  }
  // compute argmin
  int argmin =
      std::distance(distanceSq.begin(),
                    std::min_element(distanceSq.begin(), distanceSq.end()));

  com2D(argmin, cx, cy);

  std::cout << "\ndistanceSq[argmin] = " << distanceSq[argmin] << ", armgin = " << argmin << '\n';
  std::cout << "cx[argmin], cy[argmin] = " << cx << '\t' << cy << '\n';
  // std::cout << "\nexiting getCellIndexHere\n";
  return argmin;
}

void epi2D::deleteCell(double sizeRatio, int nsmall, double xLoc, double yLoc) {
  /*need to touch smallN, largeN, NCELLS, NVTOT, cellDOF, vertDOF, szList, nv,
  list, vvel, vpos, vF, vrad, im1, ip1, vim1, vip1, a0, l0, NCTCS, cij, calA0,
  psi, DrList (this list of items is from jamFracture.cpp) and possibly others

  deleteCell effectively deletes a cell by erasing 1 element from vectors who
  have size = NCELLS and by erasing largeNV or smallNV elements from vectors who
  have size = NVTOT
  */
  int vim1, vip1;
  int smallNV = nsmall;
  int largeNV = round(sizeRatio * smallNV);

  // cell index of center cell, to be deleted
  int deleteIndex = getIndexOfCellLocatedHere(xLoc, yLoc);
  double cx, cy;
  com2D(deleteIndex, cx, cy);
  cout << "deleting cell " << deleteIndex << " at position " << cx << '\t' << cy << '\n';

  // isDeleteLarge is true if deleting a large particle, false if small.
  bool isDeleteLarge = (nv[deleteIndex] == largeNV);
  if (nv[deleteIndex] != largeNV && nv[deleteIndex] != smallNV)
    throw std::invalid_argument("nv does not correspond to large or small\n");
  NCELLS -= 1;

  // re-count contact matrix after decrementing NCELLS. Delete cell => delete
  // row => delete NCELLS-1 entries location of deletion doesn't matter because
  // contact matrix is reset each integration step.
  std::cout << "NCELLS = " << NCELLS << '\n';

  cij.erase(cij.begin(), cij.begin() + (NCELLS - 1));

  // number of vertices to delete
  int numVertsDeleted = largeNV * isDeleteLarge + smallNV * !isDeleteLarge;

  // total number of vertices
  NVTOT -= numVertsDeleted;

  // degree of freedom counts
  vertDOF = NDIM * NVTOT;

  // sum up number of vertices of each cell until reaching the cell to delete
  int sumVertsUntilGlobalIndex = szList[deleteIndex];

  // adjust szList and nv, which keep track of global vertex indices
  // szList stores gi of each cell. To account for a deleted particle, delete
  // one index, then subtract numVertsDeleted from successive indices

  for (auto i = szList.begin() + deleteIndex; i != szList.end(); i++) {
    *i -= numVertsDeleted;
  }

  szList.erase(szList.begin() + deleteIndex);

  // nv,a0,calA0 have dimension (NCELLS), so need to remove the correct cell
  // entry
  nv.erase(nv.begin() + deleteIndex);
  a0.erase(a0.begin() + deleteIndex);
  psi.erase(psi.begin() + deleteIndex);
  flag.erase(flag.begin() + deleteIndex);
  flagPos.erase(flagPos.begin() + deleteIndex);
  activePropulsionFactor.erase(activePropulsionFactor.begin() + deleteIndex);
  cellU.erase(cellU.begin() + deleteIndex);
  fieldStressCells.erase(fieldStressCells.begin() + deleteIndex);
  fieldShapeStressCells.erase(fieldShapeStressCells.begin() + deleteIndex);
  timeElapsedSinceFlagPlanted.erase(timeElapsedSinceFlagPlanted.begin() + deleteIndex);
  restLengthLPx.erase(restLengthLPx.begin() + deleteIndex);
  restLengthLPy.erase(restLengthLPy.begin() + deleteIndex);

  int deleteIndexGlobal = gindex(deleteIndex, 0);
  cout << "deleteIndexGlobal = " << deleteIndexGlobal << '\n';

  // remove one largeNV or smallNV worth of indices from vectors with dimension
  // (NVTOT)
  list.erase(list.begin() + deleteIndexGlobal,
             list.begin() + deleteIndexGlobal + numVertsDeleted);
  r.erase(r.begin() + deleteIndexGlobal,
          r.begin() + deleteIndexGlobal + numVertsDeleted);
  l0.erase(l0.begin() + deleteIndexGlobal,
           l0.begin() + deleteIndexGlobal + numVertsDeleted);
  fieldStress.erase(fieldStress.begin() + deleteIndexGlobal,
                    fieldStress.begin() + deleteIndexGlobal + numVertsDeleted);
  fieldShapeStress.erase(fieldShapeStress.begin() + deleteIndexGlobal,
                         fieldShapeStress.begin() + deleteIndexGlobal +
                             numVertsDeleted);
  vnn.erase(vnn.begin() + deleteIndexGlobal,
            vnn.begin() + deleteIndexGlobal + numVertsDeleted);
  vnn_label.erase(vnn_label.begin() + deleteIndexGlobal,
                  vnn_label.begin() + deleteIndexGlobal + numVertsDeleted);

  for (int i = 0; i < nv[deleteIndex]; i++) {
    cout << x[NDIM * sumVertsUntilGlobalIndex + 2 * i] << '\t' << x[NDIM * sumVertsUntilGlobalIndex + 2 * i + 1] << '\n';
  }

  // remove an entire cell of indices (NDIM*numVertsDeleted), for vectors of
  // dimension (vertDOF)
  v.erase(v.begin() + NDIM * sumVertsUntilGlobalIndex,
          v.begin() + NDIM * (sumVertsUntilGlobalIndex + numVertsDeleted));
  x.erase(x.begin() + NDIM * sumVertsUntilGlobalIndex,
          x.begin() + NDIM * (sumVertsUntilGlobalIndex + numVertsDeleted));
  F.erase(F.begin() + NDIM * sumVertsUntilGlobalIndex,
          F.begin() + NDIM * (sumVertsUntilGlobalIndex + numVertsDeleted));

  /*for (int i = 0; i < nv[deleteIndex]; i++){
    cout << x[NDIM*sumVertsUntilGlobalIndex+2*i] << '\t' << x[NDIM*sumVertsUntilGlobalIndex + 2*i + 1] << '\n';
  }*/

  cout << "done printing copy \n";

  // save list of adjacent vertices
  im1 = vector<int>(NVTOT, 0);
  ip1 = vector<int>(NVTOT, 0);

  for (int ci = 0; ci < NCELLS; ci++) {
    for (int vi = 0; vi < nv[ci]; vi++) {
      // wrap local indices
      vim1 = (vi - 1 + nv[ci]) % nv[ci];
      vip1 = (vi + 1) % nv[ci];

      // get global wrapped indices
      int gi = gindex(ci, vi);
      im1[gi] = gindex(ci, vim1);
      ip1[gi] = gindex(ci, vip1);
    }
  }
}

void epi2D::deleteVertex(std::vector<int>& deleteList) {
  // delete the vertices with indices in deleteList. option to interpolate
  // (attempt to distribute remaining vertices into space left behind by deleted
  // old vertex)
  cout << "\ndeleting vertices!!!\n\n simclock = " << simclock
       << '\n';
  cout << "NVTOT = " << NVTOT << '\n';
  cout << "r.size() = " << r.size() << '\n';
  cout << "deleteList size = " << deleteList.size() << '\n';
  // sort descending so deletion doesn't interfere with itself
  std::sort(deleteList.begin(), deleteList.end(), std::greater<int>());
  for (auto i : deleteList)
    cout << " element of deleteList = " << i << " with vnn_label "
         << vnn_label[i] << '\n';

  // need to delete an index from anything that depends on NVTOT, nv, szList
  // deleteList holds gi of indices to delete
  int vim1, vip1, ci, vi;
  for (auto i : deleteList) {
    int ci, vi;
    cindices(ci, vi, i);
    nv[ci] -= 1;
    NVTOT -= 1;
    vertDOF -= NDIM;
    // cij deletion is wrong atm
    cij.erase(cij.begin());
    list.erase(list.begin() + i);
    r.erase(r.begin() + i);
    l0.erase(l0.begin() + i);
    vnn.erase(vnn.begin() + i);
    vnn_label.erase(vnn_label.begin() + i);
    // initialRadius, initiall0 are handled in dampedNP0
    sortedWoundIndices.erase(std::remove(sortedWoundIndices.begin(), sortedWoundIndices.end(), i), sortedWoundIndices.end());
    // sortedWoundIndices relies on gi, so need to decrement after deletion
    for (auto& j : sortedWoundIndices)
      if (j > i)
        j--;

    fieldStress.pop_back();
    fieldShapeStress.pop_back();

    // vector.erase erases [start,stop) so this deletes 2 elements
    v.erase(v.begin() + NDIM * i, v.begin() + NDIM * (i + 1));
    x.erase(x.begin() + NDIM * i, x.begin() + NDIM * (i + 1));
    F.erase(F.begin() + NDIM * i, F.begin() + NDIM * (i + 1));

    for (int j = 1; j < NCELLS; j++) {
      szList[j] = szList[j - 1] + nv[j - 1];
    }
  }

  // save list of adjacent vertices
  im1.resize(NVTOT);
  ip1.resize(NVTOT);

  for (int ci = 0; ci < NCELLS; ci++) {
    for (int vi = 0; vi < nv[ci]; vi++) {
      // wrap local indices
      vim1 = (vi - 1 + nv[ci]) % nv[ci];
      vip1 = (vi + 1) % nv[ci];

      // get global wrapped indices
      int gi = gindex(ci, vi);
      im1[gi] = gindex(ci, vim1);
      ip1[gi] = gindex(ci, vip1);
    }
  }
  cout << "exiting deleteVertex\n";
}

void epi2D::laserAblate(int numCellsAblated, double sizeRatio, int nsmall, double xLoc, double yLoc) {
  for (int i = 0; i < numCellsAblated; i++) {
    deleteCell(sizeRatio, nsmall, xLoc, yLoc);
  }
  cout << "deleted " << numCellsAblated << "cells!\n";
  zeroMomentum();
}

void epi2D::initializevnn() {
  vnn.resize(NVTOT);
  for (int gi = 0; gi < NVTOT; gi++) {
    vnn[gi].resize(7);
  }
  vnn_label.resize(NVTOT);
}

void epi2D::boundaries() {
  // construct vertex neighbor list for void detection. This routine only
  // handles vertices in the same cell. Neighbors via adhesion are handled
  // externally, in the force routine
  int gi = 0, gi0 = 0;
  for (int gi = 0; gi < vnn.size(); gi++) {
    fill(vnn[gi].begin(), vnn[gi].end(), -1);
  }

  for (int ci = 0; ci < NCELLS; ci++) {
    // get first index of cell ci
    gi0 = gindex(ci, 0);

    // loop over vertices in cell vi, put forward and backward neighbors into
    // nn[0], nn[1]
    for (int vi = 0; vi < nv[ci]; vi++) {
      gi = gindex(ci, vi);
      vnn[gi][0] = gi0 + (vi + 1 + nv[ci]) % nv[ci];
      vnn[gi][1] = gi0 + (vi - 1 + nv[ci]) % nv[ci];
    }
  }
}

std::vector<int> epi2D::refineBoundaries() {
  // refine the boundaries, return site occupations for Newman-Ziff.
  std::vector<int> voidFacingVertexIndices;
  int c1, c2, c3, v1, v2, v3;  // stores dangling_end_labels for neighbor decision later

  int counter, k;
  fill(vnn_label.begin(), vnn_label.end(), -1);
  // first refinement: vertices with exactly 2 positive labels are non-adhering,
  // therefore void-adjacent (label 0)
  for (int gi = 0; gi < NVTOT; gi++) {
    // tally up how many neighbors each cell has
    counter = 0;
    for (int i = 0; i < vnn[gi].size(); i++) {
      counter += (vnn[gi][i] >= 0);
    }
    if (counter < 2)
      cout << "warning: definitely should not have less than 2 neighbors for "
              "this vertex!!\n";
    else if (counter == 2) {
      vnn_label[gi] = 0;
    }
  }

  // second refinement: a vertex gets a 1 label if it is unlabeled and is next
  // to a void-adjacent vertex. It is an edge.
  for (int gi = 0; gi < NVTOT; gi++) {
    // check if gi is a non-adhering vertex
    if (vnn_label[gi] == 0) {
      for (int j = 0; j < vnn[gi].size(); j++) {
        k = vnn[gi][j];
        // give 1 labels to all non-zero labeled valid neighbors
        if (k >= 0 && vnn_label[k] != 0)
          vnn_label[k] = 1;
      }
    }
  }

  // third refinement: a vertex gets a 2 label if it is 1-labeled and next to a
  // 1 or a 2. It is a corner (2 or more edges).
  for (int gi = 0; gi < NVTOT; gi++) {
    if (vnn_label[gi] == 1) {
      for (auto i : vnn[gi]) {
        if (vnn_label[i] == 1 || vnn_label[i] == 2) {
          vnn_label[gi] = 2;
          vnn_label[i] = 2;
        }
      }
    }
  }

  // NEW: fourth refinement: a dangling end is defined as a 1-labeled vertex at this stage.
  // To rescue two corner cases: cells with 1 or 2 vertices on the boundary,
  // which aren't caught by the previous refinements.
  // Proceed by looking at the unlabeled neighbors of dangling ends. If an
  // unlabeled neighbor attached to dangling
  //    end gi has another unlabeled neighbor attached to dangling end gj != gi
  //    Then we have found 3 or 4 corners, so we mark both dangling ends and 1
  //    or 2 vertices in-between as corners.

  std::vector<int> dangling_end_label(NVTOT, -1);  // stores which vertex the dangling neighbor is a neighbor of
  for (int gi = 0; gi < NVTOT; gi++) {
    if (vnn_label[gi] == 1) {
      for (auto i : vnn[gi]) {
        // assign unlabeled neighbors of dangling ends a -2 label.
        if (vnn_label[i] == -1 || vnn_label[i] == -2) {
          vnn_label[i] = -2;

          // if i (a dangling neighbor) already belongs to another vertex
          //  let i belong to whichever vertex is in the same cell.
          // Here, gi is the new potential owner, and dangling_end_label[i] is the old owner
          if (dangling_end_label[i] != -1) {
            cindices(c1, v1, i);
            cindices(c2, v2, dangling_end_label[i]);
            cindices(c3, v3, gi);
            if (c1 == c3 && c1 != c2) {  // i is in same cell as gi, but not its old owner
              dangling_end_label[i] = gi;
              dangling_end_label[gi] = gi;
            }
            // otherwise, do nothing (i.e. keep the old owner)
          }
          // otherwise, assign labels as normal.
          else {
            dangling_end_label[i] = gi;
            dangling_end_label[gi] = gi;
          }
        }
      }
    }
  }

  for (int gi = 0; gi < NVTOT; gi++) {
    // for all dangling ends, if it has more than 3 dangling end
    //  neighbors then it is not a corner and should be ignored (set to -3)
    // 11/29/21 : not sure if this is a good criterion..
    if (vnn_label[gi] == 1) {
      int danglingCounter = 0;
      for (auto i : vnn[gi]) {
        if (vnn_label[i] == -2 || vnn_label[i] == -3)
          danglingCounter++;
      }
      if (danglingCounter > 3) {
        vnn_label[gi] = -3;
      }
    }
  }

  for (int gi = 0; gi < NVTOT; gi++) {
    // if vertex is an unlabeled neighbor of a dangling end
    if (vnn_label[gi] == -2) {
      // check that vertex's neighbors
      for (auto j : vnn[gi]) {
        // if its neighbor is also an unlabeled neighbor of a dangling end, and
        // of a different dangling end, all are corners
        if (vnn_label[j] == -2 &&
            dangling_end_label[gi] != dangling_end_label[j] &&
            dangling_end_label[j] > 0) {
          vnn_label[gi] = 2;
          vnn_label[j] = 2;
          vnn_label[dangling_end_label[j]] = 2;
          vnn_label[dangling_end_label[gi]] = 2;
        }
      }
    }
  }
  // Done with refinements. Store void adjacent vertices and corner vertices.
  for (int gi = 0; gi < NVTOT; gi++) {
    // get an occupation order for void-facing indices. 0 is edge, 2 is corner
    if (vnn_label[gi] == 0 || vnn_label[gi] == 2)
      voidFacingVertexIndices.push_back(gi);
  }

  return voidFacingVertexIndices;
}

void epi2D::NewmanZiff(std::vector<int>& ptr, int empty, int& mode, int& big, std::vector<int>& order) {
  // union/find algorithm that seeks contiguous clusters of points.
  //  the root (ptr) of each node is the first vertex in its cluster, and the
  //  root of that node is the negative size of the cluster
  // order is the occupation order for void-facing indices. It is a set of
  // void-facing vi.
  //  mode is the cluster with largest size, big is the largest size
  int s1, s2, r1, r2;
  for (int i = 0; i < order.size(); i++) {
    r1 = s1 = order[i];
    ptr[s1] = -1;
    for (int j = 0; j < vnn[s1].size(); j++) {
      // key adjustment from Newman-Ziff: skip if neighbor list has -1 label
      // (means no neighbor is there)
      if (vnn[s1][j] == -1)
        continue;
      s2 = vnn[s1][j];
      if (ptr[s2] != empty) {
        r2 = findRoot(s2, ptr);
        if (r2 != r1) {
          if (ptr[r1] > ptr[r2]) {
            ptr[r2] += ptr[r1];
            ptr[r1] = r2;
            r1 = r2;
          } else {
            ptr[r1] += ptr[r2];
            ptr[r2] = r1;
          }
          if (-ptr[r1] > big) {
            big = -ptr[r1];
            mode = r1;
          }
        }
      }
    }
  }
}

// printBoundaries will determine the (nth) largest cluster size, save it to the epi2D object, and print it.
//  affects where purse-string is
//  also saves sortedWoundIndices
void epi2D::printBoundaries(int nthLargestCluster) {
  cerr << "entering printBoundaries\n";
  //  empty is the maximum negative cluster size
  int empty = -NVTOT - 1;
  // mode is the statistical mode of the cluster labels
  int mode;
  // order is the occupation/population order for our neighbor topology
  order = refineBoundaries();
  cout << "order size in printBoundaries = " << order.size() << '\n';
  // ptr of a root node gives the negative cluster size, of a non-root gives the
  // next level in the tree, of an empty gives empty
  std::vector<int> ptr(NVTOT, empty);
  int big = 0;
  bool choseCorner;
  // fill ptr with a clustered network of nodes
  NewmanZiff(ptr, empty, mode, big, order);

  cout << "mode = " << mode << '\n';
  cout << "with size = " << big << '\n';
  // if I want something other than the largest cluster, change mode to identify
  // one of the smaller sized clusters
  if (nthLargestCluster > 1) {
    cout << "requested " << nthLargestCluster
         << "-th most populated cluster, so recalculating " << nthLargestCluster
         << "-mode value\n";
    std::vector<int> clusterSize;
    std::vector<int> clusterRoot;
    for (int gi = 0; gi < NVTOT; gi++) {
      if (ptr[gi] < 0 && ptr[gi] != empty) {
        clusterSize.push_back(-ptr[gi]);
        clusterRoot.push_back(gi);
      }
    }
    sort(clusterSize.begin(), clusterSize.end(), greater<int>());
    while (nthLargestCluster - 1 >= clusterSize.size()) {
      // requested out of range access, decrement the index until in range
      // has the effect of not giving the nth largest cluster, but rather the largest cluster
      nthLargestCluster--;
    }
    int nthClusterSize = clusterSize[nthLargestCluster - 1];
    for (auto i : clusterRoot) {
      if (-ptr[i] == nthClusterSize)
        mode = i;
    }
    cout << "new mode = " << mode << '\n';
    big = nthClusterSize;
    cout << "new size (big) = " << big << '\n';
  }

  // in no particular order, print out wrapped locations of all main cluster vertices
  for (int gi = 0; gi < NVTOT; gi++) {
    if (vnn_label[gi] == 0 || vnn_label[gi] == 2 || vnn_label[gi] == 1 || vnn_label[gi] == -1 || vnn_label[gi] == -2 || vnn_label[gi] == -3) {
      edgeout << x[gi * NDIM] << '\t' << x[gi * NDIM + 1] << '\t' << vnn_label[gi] + (findRoot(gi, ptr) == mode) * 10 << '\n';
    }
  }

  // in order, print out the unwrapped locations of all main cluster vertices,
  // starting with a vertex in the big cluster
  cerr << "before ordering vertices\n";
  int it = 0, middleit = 1000, maxit = 2000;
  int current_vertex = -2;
  bool stopSignal = false;
  double currentx, currenty, nextx, nexty;
  std::vector<int> previous_vertex;
  std::vector<double> unwrapped_x;
  std::vector<int> unwrapped_x_gi;  // indices of vertices in unwrapped_x
  for (int gi = 0; gi < NVTOT; gi++) {
    if (findRoot(gi, ptr) == mode) {
      currentx = x[gi * NDIM];
      currenty = x[gi * NDIM + 1];
      unwrapped_x.push_back(currentx);
      unwrapped_x.push_back(currenty);
      unwrapped_x_gi.push_back(gi);
      current_vertex = gi;
      break;
    }
  }
  while (previous_vertex.size() != big && stopSignal == false) {
    // move current vertex to list of old vertices
    if (previous_vertex.size() < 2 ||
        previous_vertex.end()[-1] != current_vertex) {
      it = 0;
      previous_vertex.push_back(current_vertex);
    }
    // for debugging:
    int badVertexIndex;
    // end debugging

    for (auto i : vnn[current_vertex]) {
      badVertexIndex = i;
      // if i == -1, then we've found a dead end. backtrack (go back by one current_vertex) and try again
      if (i == -1) {
        current_vertex = previous_vertex.end()[-2];
        // unwrapped_x.pop_back();
        // unwrapped_x.pop_back();
        break;
      }

      // check i valid neighbor, points to the biggest cluster, and is not
      // already in the previous vertex list
      if (i >= 0 && findRoot(i, ptr) == mode &&
          std::find(previous_vertex.begin(), previous_vertex.end(), i) ==
              previous_vertex.end()) {
        // switch focus to new vertex
        current_vertex = i;
        // (nextx, nexty) is the nearest image location of the next vertex
        nextx = x[current_vertex * NDIM] - currentx;
        nexty = x[current_vertex * NDIM + 1] - currenty;

        nextx -= L[0] * round(nextx / L[0]);
        nexty -= L[1] * round(nexty / L[1]);

        nextx += currentx;
        nexty += currenty;

        unwrapped_x.push_back(nextx);
        unwrapped_x.push_back(nexty);
        unwrapped_x_gi.push_back(current_vertex);

        currentx = nextx;
        currenty = nexty;
        break;

      } else if (it > middleit) {  // we have formed a closed loop, but have not
                                   // traversed the entire cluster.
                                   // if the closed loop is large enough, we'll accept it as the main
                                   // closed loop and stop adding vertices
                                   // if (previous_vertex.size() >= big / 2) {
        if (previous_vertex.size() >= big / 2) {
          stopSignal = true;
          cout << "truncating; found a closed loop that's large enough!\n";
          break;
        }
        // If cluster is too small to be the main loop, seek an alternate route.
        for (auto j : previous_vertex) {
          // go through previously added vertices, look for their neighbors that
          // are not already added.
          for (auto k : vnn[j]) {
            if (k >= 0 && findRoot(k, ptr) == mode &&
                std::find(previous_vertex.begin(), previous_vertex.end(), k) ==
                    previous_vertex.end()) {
              it = 0;
              current_vertex = k;
              nextx = x[current_vertex * NDIM] - currentx;
              nexty = x[current_vertex * NDIM + 1] - currenty;

              nextx -= L[0] * round(nextx / L[0]);
              nexty -= L[1] * round(nexty / L[1]);

              nextx += currentx;
              nexty += currenty;

              // forget the closed subloop, move onto a new one. Keep
              // previous_vertex, which forbids subloop from joining new loop
              unwrapped_x.clear();
              unwrapped_x.push_back(nextx);
              unwrapped_x.push_back(nexty);
              unwrapped_x_gi.push_back(current_vertex);

              currentx = nextx;
              currenty = nexty;
              cout << "cluster is too small, seeking alternate route\n";
              break;
            }
          }
        }
      }
    }
    it++;
    if (it >= maxit) {
      for (auto i : vnn[current_vertex])
        cout << "broken vertex is: " << current_vertex
             << "\nneighbors of broken vertex are: " << i << ", with root "
             << findRoot(i, ptr)
             << ", and bool for whether it's in current list: "
             << (std::find(previous_vertex.begin(), previous_vertex.end(), i) !=
                 previous_vertex.end())
             << '\n'
             << "at position " << x[i * NDIM] << '\t' << x[i * NDIM + 1] << '\n'
             << "size of vertex list = " << previous_vertex.size() << '\n';
      break;
    }
  }
  if (checkWoundClosedPolygon(unwrapped_x_gi)) {
    for (int i = 0; i < unwrapped_x.size(); i += 2)
      bout << unwrapped_x[i] << '\t' << unwrapped_x[i + 1] << '\n';
  } else {  // unwrapped_x_gi is not closed, hence does not represent what I'm using for the void segmentation
    for (int i = 0; i < sortedWoundIndices.size() / 2; i += 2) {
      bout << x[NDIM * sortedWoundIndices[i]] << '\t' << x[NDIM * sortedWoundIndices[i] + 1] << '\n';
    }
  }
  // indicate end of block, i.e. end of boundary matrix for this frame
  bout << "*EOB\n";
  edgeout << "*EOB\n";
  cout << "leaving printBoundaries\n";
}

void epi2D::getWoundVertices(int nthLargestCluster) {
  std::vector<int> temporaryWoundIndices;
  // empty is the maximum negative cluster size
  int empty = -NVTOT - 1;
  // mode is the statistical mode of the cluster labels
  int mode;
  // order is the occupation/population order for our neighbor topology
  order = refineBoundaries();
  // ptr of a root node gives the negative cluster size, of a non-root gives the
  // next level in the tree, of an empty gives empty
  std::vector<int> ptr(NVTOT, empty);
  int big = 0;
  // fill ptr with a clustered network of nodes
  NewmanZiff(ptr, empty, mode, big, order);

  if (nthLargestCluster > 1) {
    std::vector<int> clusterSize;
    std::vector<int> clusterRoot;
    for (int gi = 0; gi < NVTOT; gi++) {
      if (ptr[gi] < 0 && ptr[gi] != empty) {
        clusterSize.push_back(-ptr[gi]);
        clusterRoot.push_back(gi);
      }
    }
    sort(clusterSize.begin(), clusterSize.end(), greater<int>());
    int nthClusterSize = clusterSize[nthLargestCluster - 1];
    for (auto i : clusterRoot) {
      if (-ptr[i] == nthClusterSize)
        mode = i;
    }
    // cout << "new mode in getWoundVertices() = " << mode << '\n';
    big = nthClusterSize;
    // cout << "new size (big) in getWoundVertices() = " << big << '\n';
  }

  // in order, print out the unwrapped locations of all main cluster vertices,
  // starting with a vertex in the big cluster
  // cerr << "before ordering vertices\n";
  int it = 0, middleit = 1000, maxit = 2000;
  int current_vertex = -2;
  bool stopSignal = false, tempWoundIsClosedPolygon = false;
  double currentx, currenty, nextx, nexty;
  std::vector<int> previous_vertex;
  // std::vector<double> unwrapped_x;
  for (int gi = 0; gi < NVTOT; gi++) {  // the lowest gi is the first entry of temporaryWoundIndices, the rest are sorted in based on connectivity
    if (findRoot(gi, ptr) == mode) {
      temporaryWoundIndices.push_back(gi);
      current_vertex = gi;
      break;
    }
  }
  while (previous_vertex.size() != big && stopSignal == false) {
    // move current vertex to list of old vertices
    if (previous_vertex.size() < 2 ||
        previous_vertex.end()[-1] != current_vertex) {
      it = 0;
      previous_vertex.push_back(current_vertex);
    }
    // for debugging:
    int badVertexIndex;
    // end debugging

    for (auto i : vnn[current_vertex]) {
      badVertexIndex = i;
      // if i == -1, then we've found a dead end. backtrack (go back by one current_vertex) and try again
      // another possibility is if i==-1, then we've closed the wound. Check for this
      if (i == -1) {
        current_vertex = previous_vertex.end()[-2];
        if (checkWoundClosedPolygon(temporaryWoundIndices)) {
          tempWoundIsClosedPolygon = true;
        }
        break;
      }

      // check i valid neighbor, points to the biggest cluster, and is not
      // already in the previous vertex list
      if (i >= 0 && findRoot(i, ptr) == mode &&
          std::find(previous_vertex.begin(), previous_vertex.end(), i) ==
              previous_vertex.end()) {
        // switch focus to new vertex
        current_vertex = i;

        temporaryWoundIndices.push_back(current_vertex);
        break;

      } else if (it > middleit) {  // we have formed a closed loop, but have not
                                   // traversed the entire cluster.
                                   // if the closed loop is large enough, we'll accept it as the main
                                   // closed loop and stop adding vertices
                                   // if (previous_vertex.size() >= big / 2) {
        if (previous_vertex.size() >= big / 2) {
          stopSignal = true;
          cout << "truncating; found a closed loop that's large enough!\n";
          break;
        }
        // If cluster is too small to be the main loop, seek an alternate route.
        for (auto j : previous_vertex) {
          // go through previously added vertices, look for their neighbors that
          // are not already added.
          for (auto k : vnn[j]) {
            if (k >= 0 && findRoot(k, ptr) == mode &&
                std::find(previous_vertex.begin(), previous_vertex.end(), k) ==
                    previous_vertex.end()) {
              it = 0;
              current_vertex = k;

              // forget the closed subloop, move onto a new one. Keep
              // previous_vertex, which forbids subloop from joining new loop
              temporaryWoundIndices.clear();
              temporaryWoundIndices.push_back(current_vertex);

              cout << "cluster is too small, seeking alternate route\n";
              break;
            }
          }
        }
      }
    }
    // if wound indices form a closed shape, we can move on. otherwise, move to escape clause
    if (tempWoundIsClosedPolygon)
      break;

    it++;
    if (it >= maxit) {
      for (auto i : vnn[current_vertex])
        cout << "broken vertex is: " << current_vertex
             << "\nneighbors of broken vertex are: " << i << ", with root "
             << findRoot(i, ptr)
             << ", and bool for whether it's in current list: "
             << (std::find(previous_vertex.begin(), previous_vertex.end(), i) !=
                 previous_vertex.end())
             << '\n'
             << "at position " << x[i * NDIM] << '\t' << x[i * NDIM + 1] << '\n'
             << "size of vertex list = " << previous_vertex.size() << '\n';
      break;
    }
  }
  if (checkWoundClosedPolygon(temporaryWoundIndices)) {
    sortedWoundIndices = temporaryWoundIndices;
  } else {
    // cout << "temporaryWoundIndices was not a closed shape, use old shape for sortedWoundIndices.\n";
  }
}

bool epi2D::checkWoundClosedPolygon(std::vector<int>& listOfIndices) {
  // true if listOfIndices (wound) is a closed polygon, false if not
  if (listOfIndices.size() == 0)
    return false;
  int firstVertex = listOfIndices[0];
  int lastVertex = listOfIndices.back();
  return (std::find(vnn[firstVertex].begin(), vnn[firstVertex].end(), lastVertex) !=
          vnn[firstVertex].end());
}

double epi2D::computeWoundVerticesUsingRays(double& woundCenterX, double& woundCenterY, int numRays) {
  // draw rays from the wound center in order to determine the nearest vertex at each angle, which approximates the wound
  // returns the area of the wound
  double theta, dist;
  std::vector<int> newWoundList;
  for (int i = 0; i < numRays; i++) {
    theta = i * 2 * PI / numRays;
    for (auto gi : sortedWoundIndices) {
      // this restriction only considers vertices in the original or current? segmentation
      dist = distanceLineAndPoint(woundCenterX, woundCenterY, 1000 * r[0] * cos(theta), -1000 * r[0] * sin(theta), x[gi * NDIM], x[gi * NDIM + 1]);  // could try to not use sqrt in function, later
      if (dist < r[gi]) {
        if (newWoundList.size() == 0 || newWoundList.back() != gi)
          newWoundList.push_back(gi);
        break;
      }
    }
  }
  sortedWoundIndices = newWoundList;

  // calculate area

  double xi = x[NDIM * sortedWoundIndices[0]];
  double yi = x[NDIM * sortedWoundIndices[0] + 1];
  double xip1, yip1, areaVal;
  for (int i = 1; i < sortedWoundIndices.size(); i++) {
    double dx = x[NDIM * sortedWoundIndices[i]] - xi;
    double dy = x[NDIM * sortedWoundIndices[i] + 1] - yi;
    xip1 = xi + dx;
    yip1 = yi + dy;
    areaVal += xi * yip1 - xip1 * yi;
    xi = xip1;
    yi = yip1;
  }
  areaVal *= 0.5;
  return abs(areaVal);
  cout << "wound area = " << areaVal << ", simclock = " << simclock << '\n';
}

int epi2D::findRoot(int i, std::vector<int>& ptr) {
  if (ptr[i] < 0)
    return i;
  return ptr[i] = findRoot(ptr[i], ptr);
}

double epi2D::calculateWoundArea(double& woundPointX, double& woundPointY) {
  // input: a point previously inside the wound.
  //  we will nucleate the area calculation around this point
  // this algorithm gives the area of a wound by dividing up the simulation box into a grid
  //  each grid point is checked for being within any cell and for proximity to cell vertices
  //  if either is true, give that point a value of 1 (true)
  // then take woundPoint and calculate the largest continuous blob of 0 values in occupancyMatrix.
  //  multiply by grid point area to get the total wound area.

  // note: running this every 100 or 1000 timesteps should be fine.
  std::vector<double> posX(NVTOT / 2), posY(NVTOT / 2);
  for (int i = 0; i < NVTOT; i++) {
    if (i % 2 == 0)
      posX[i / 2] = x[NDIM * i];
    else
      posY[(i - 1) / 2] = x[NDIM * (i - 1) + 1];
  }
  // resolution is the unit box length that we'll use to map occupancyMatrix to real coordinates
  double resolution = r[0] / 2.0 / 2.0;
  if (resolution <= 1e-10)
    cerr << "bug: resolution in calculateWoundArea is zero\n";

  double xLow = *std::min_element(posX.begin(), posX.end());
  double xHigh = *std::max_element(posX.begin(), posX.end());
  double yLow = *std::min_element(posY.begin(), posY.end());
  double yHigh = *std::max_element(posY.begin(), posY.end());
  if (fabs(woundPointX) > xHigh) {
    cout << "woundPoint does not lie within the bounds of the simulation box, probably failed to find the center of a wound\n returning 0 area, skipping.\n";
    return 0.0;
  }

  int xResolution = (xHigh - xLow) / resolution;
  int yResolution = (yHigh - yLow) / resolution;
  std::vector<std::vector<bool>> occupancyMatrix(xResolution, std::vector<bool>(yResolution, 0));
  for (int i = 0; i < xResolution; i++) {
    for (int j = 0; j < yResolution; j++) {
      // first pass: occupancy is 1 if inside a cell, 0 if not inside a cell (pnpoly point inclusion in polygon test algorithm)
      occupancyMatrix[i][j] = !isPointInPolygons(i * resolution, j * resolution);

      // second pass: if occupancy is 0, set occupancy back to 1 if within vrad (previously resolution) of a vertex
      if (occupancyMatrix[i][j] == 0) {
        for (int k = 0; k < NVTOT; k++) {
          if (fabs(i * resolution - x[NDIM * k]) < r[0] && fabs(j * resolution - x[NDIM * k + 1]) < r[0]) {
            occupancyMatrix[i][j] = 1;
            break;
          }
        }
      }
    }
  }
  /*if (simclock > 100)
    for (int i = 0; i < xResolution; i++) {
      cout << "[ ";
      for (int j = 0; j < yResolution; j++) {
        cout << occupancyMatrix[i][j];
      }
      cout << "]" << '\n';
    }
  */

  // now occupancy matrix is filled. find the largest cluster of open space, given a point within the wound.
  int woundPointXIndex = woundPointX / resolution;
  int woundPointYIndex = woundPointY / resolution;

  // check if the given point is not within the wound
  if (occupancyMatrix[woundPointXIndex][woundPointYIndex] == 1) {
    // given point is not in the wound, so we need to look for a nearby point that's hopefully in the wound (searching diagonally for now)
    int searchRange = 5, offset = 0;
    int newXIndex = woundPointXIndex, newYIndex = woundPointYIndex;
    while (occupancyMatrix[newXIndex][newYIndex] == 1 && offset <= searchRange) {
      offset++;
      if (occupancyMatrix[woundPointXIndex + offset][woundPointYIndex + offset] == 0) {
        newXIndex = woundPointXIndex + offset;
        newYIndex = woundPointYIndex + offset;
      }
      if (occupancyMatrix[woundPointXIndex - offset][woundPointYIndex - offset] == 0) {
        newXIndex = woundPointXIndex - offset;
        newYIndex = woundPointYIndex - offset;
      }
      if (occupancyMatrix[woundPointXIndex - offset][woundPointYIndex + offset] == 0) {
        newXIndex = woundPointXIndex - offset;
        newYIndex = woundPointYIndex + offset;
      }
      if (occupancyMatrix[woundPointXIndex + offset][woundPointYIndex - offset] == 0) {
        newXIndex = woundPointXIndex + offset;
        newYIndex = woundPointYIndex - offset;
      }
    }
    if (offset > searchRange) {
      cout << "failed to find a nearby point identifiable as a wound, returning 0.0 for area\n";
      return 0.0;
    }
    // set wound point to the newly identified point
    woundPointX = newXIndex * resolution;
    woundPointY = newYIndex * resolution;
    woundPointXIndex = woundPointX / resolution;
    woundPointYIndex = woundPointY / resolution;
  }

  std::vector<int> emptyGridIndices;
  bool done = false;

  // algorithm: connected component labeling
  int current_label = 1;
  std::vector<std::vector<int>> labels(xResolution, std::vector<int>(yResolution, 0));

  // choose initial pixel
  int i = woundPointXIndex, j = woundPointYIndex, nni, nnj;
  // if pixel has value 0 and is unlabeled, give it current label and add it to queue.
  if (occupancyMatrix[i][j] == 0 && labels[i][j] == 0) {
    emptyGridIndices.push_back(i * yResolution + j);
    labels[i][j] = current_label;

    // pop out an element from queue, look at its neighbors. If neighbor has pixel value 0 and is unlabeled, give it current label and add to queue. repeat until no more elements are in the queue.
    while (emptyGridIndices.size() > 0) {
      int current_element = emptyGridIndices.back();
      emptyGridIndices.pop_back();
      i = floor(current_element / yResolution);
      j = current_element % yResolution;
      std::vector<int> nnx = {i, i, i - 1, i + 1};
      std::vector<int> nny = {j - 1, j + 1, j, j};
      for (int k = 0; k < 4; k++) {
        nni = nnx[k];
        nnj = nny[k];
        if (nni < 0 || nni >= occupancyMatrix.size())
          continue;
        if (nnj < 0 || nnj >= occupancyMatrix[nni].size())
          continue;
        if (occupancyMatrix[nni][nnj] == 0 && labels[nni][nnj] == 0) {
          emptyGridIndices.push_back(nni * yResolution + nnj);
          labels[nni][nnj] = current_label;
        }
      }
    }
  }

  /*if (simclock > 100)
    for (int i = 0; i < xResolution; i++) {
      cout << "[ ";
      for (int j = 0; j < yResolution; j++) {
        if (i == woundPointXIndex && j == woundPointYIndex)
          cout << "x" << '\n';
        else
          cout << labels[i][j];
      }
      cout << "]" << '\n';
    }
  */

  double sum = 0.0;
  for (int i = 0; i < labels.size(); i++) {
    for (int j = 0; j < labels[i].size(); j++) {
      sum += labels[i][j];
    }
  }

  //   area should be the number of boxes times the box area
  //   alternatively, I could get the exact area by locating vertices on the edge of my newly segmented void area, but that's more computational work
  return sum * pow(resolution, 2);
}

bool epi2D::isPointInPolygons(double xloc, double yloc) {
  // check whether (xloc, yloc) is within any dp polygonal areas
  // pnpoly must be false every time in order for isInside to register as true.
  int gi;
  bool isInside = true;
  for (int ci = 0; ci < NCELLS; ci++) {
    std::vector<double> vertx(nv[ci]), verty(nv[ci]);
    for (int vi = 0; vi < nv[ci]; vi++) {
      gi = gindex(ci, vi);
      vertx[vi] = x[NDIM * gi];
      verty[vi] = x[NDIM * gi + 1];
    }
    isInside *= !pnpoly(nv[ci], vertx, verty, xloc, yloc);
  }
  return isInside;
}

int epi2D::pnpoly(int nvert, std::vector<double> vertx, std::vector<double> verty, double testx, double testy) {
  // nvert = number of vertices in polygon of interest, vert = vertex position arrays, test is the test point.
  /*
  I run a semi-infinite ray horizontally (increasing x, fixed y) out from the test point, and count how many edges it crosses. At each crossing, the ray switches between inside and outside. This is called the Jordan curve theorem.
  */
  int i, j, c = 0;
  for (i = 0, j = nvert - 1; i < nvert; j = i++) {
    if (((verty[i] > testy) != (verty[j] > testy)) &&
        (testx < (vertx[j] - vertx[i]) * (testy - verty[i]) / (verty[j] - verty[i]) + vertx[i]))
      c = !c;
  }
  return c;
}

double epi2D::calculateArea(std::vector<double>& vertx, std::vector<double>& verty) {
  // calculate area of a list of points given by (vertx,verty)
  double xi = vertx[0], yi = verty[0];
  double dx, dy, xip1, yip1, areaVal;
  for (int i = 0; i < vertx.size(); i++) {  // doesn't account for periodic boundaries yet
    dx = vertx[(i + 1) % vertx.size()] - xi;
    if (pbc[0])
      dx -= L[0] * round(dx / L[0]);
    xip1 = xi + dx;

    dy = verty[(i + 1) % verty.size()] - yi;
    if (pbc[1])
      dy -= L[1] * round(dy / L[1]);
    yip1 = yi + dy;

    // increment area
    areaVal += xi * yip1 - xip1 * yi;

    // set next coordinates
    xi = xip1;
    yi = yip1;
  }
  areaVal *= 0.5;
  return abs(areaVal);
}

double epi2D::calculateAreaFlattened(std::vector<double>& vertPosFlattened) {
  // calculate area of a list of points given by (vert[i],vert[i+1])
  double xi = vertPosFlattened[0], yi = vertPosFlattened[1];
  double dx, dy, xip1, yip1, areaVal;
  int numberOfVerts = vertPosFlattened.size() / 2;
  for (int i = 0; i < numberOfVerts; i++) {  // doesn't account for periodic boundaries yet
    dx = vertPosFlattened[NDIM * ((i + 1) % numberOfVerts)] - xi;
    if (pbc[0])
      dx -= L[0] * round(dx / L[0]);
    xip1 = xi + dx;

    dy = vertPosFlattened[NDIM * ((i + 1) % numberOfVerts) + 1] - yi;
    if (pbc[1])
      dy -= L[1] * round(dy / L[1]);
    yip1 = yi + dy;

    // increment area
    areaVal += xi * yip1 - xip1 * yi;

    // set next coordinates
    xi = xip1;
    yi = yip1;
  }
  areaVal *= 0.5;
  return abs(areaVal);
}

void epi2D::notchTest(int numCellsToDelete, double strain, double strainRate, double boxLengthScale, double sizeRatio, int nsmall, dpmMemFn forceCall, double B, double dt0, double printInterval, std::string loadingType) {
  // select xloc, yloc. delete the nearest cell. scale particle sizes, loop.
  // our proxy for isotropic stretching is to scale the particle sizes down.
  // Inter-vertex distances change,
  double xLoc, yLoc;
  double tauRelax = 1.0;
  std::cout << "inside notch test!\n";

  for (int i = 0; i < 1; i++) {
    xLoc = 0.0;
    yLoc = 0.0;
    laserAblate(numCellsToDelete, sizeRatio, nsmall, xLoc, yLoc);

    if (loadingType == "uniaxial") {
      while (L[0] / initialLx - 1 < strain) {
        cout << "current strain = " << L[0] / initialLx - 1 << '\n';
        scaleBoxSize(boxLengthScale, 1 + strainRate * tauRelax, 1);
        cout << "scaling box\n";
        dampedNVE2D(forceCall, B, dt0, tauRelax, printInterval);
        cout << "finished relaxing\n";
      }
    } else
      std::cout << "Issue: loadingType not understood\n";
  }
}

void epi2D::orientDirector(int ci, double xLoc, double yLoc) {
  // point the director of cell ci towards (xLoc, yLoc)
  double dx, dy, cx, cy;
  // compute cell center of mass
  com2D(ci, cx, cy);
  // compute angle needed for psi to point towards (xLoc,yLoc) - for now, just
  // towards origin
  double theta = atan2(cy - 0.5 * L[1], cx - 0.5 * L[0]) + PI;
  theta -= 2 * PI * round(theta / (2 * PI));
  psi[ci] = theta;
}

void epi2D::deflectOverlappingDirectors() {
  // if any two cells are overlapping (according to cij), they are candidates to
  // have their directors swapped in the same timestep directors will be swapped
  // if they are approximately 180 degrees out of phase AND pointed towards each
  // other this doesn't really catch a lot of edge cases. try a different thing.
  double dot_product = 0.0, dist1_sq, dist2_sq, psi_temp;
  int gi;
  double rix, riy, rjx, rjy, nix, niy, njx, njy;

  // tolerance for valid CIL head to head collision angles
  double angle_cutoff = cos(7 * PI / 8);

  // count number of head-to-head collisions per timestep
  polarizationCounter = 0;

  // compute directors for all cells
  vector<vector<double>> director(NCELLS, vector<double>(2, 0.0));
  for (int ci = 0; ci < NCELLS; ci++) {
    director[ci][0] = cos(psi[ci]);
    director[ci][1] = sin(psi[ci]);
  }

  for (int ci = 0; ci < NCELLS; ci++) {
    int counter = 0;

    // reset active propulsion factor
    activePropulsionFactor[ci] = 1.0;
    for (int cj = ci + 1; cj < NCELLS; cj++) {
      if (cij[NCELLS * ci + cj - (ci + 1) * (ci + 2) / 2] > 0) {
        // check if cells are overlapping
        nix = director[ci][0];
        niy = director[ci][1];
        njx = director[cj][0];
        njy = director[cj][1];
        dot_product = nix * njx + niy * njy;
        if (dot_product <= angle_cutoff && dot_product >= -1) {
          // compute center of mass of cells i and j
          com2D(ci, rix, riy);
          com2D(cj, rjx, rjy);

          // add director
          dist1_sq =
              pow(rix + nix - rjx - njx, 2) + pow(riy + niy - rjy - njy, 2);
          dist2_sq =
              pow(rix - nix - rjx + njx, 2) + pow(riy - niy - rjy + njy, 2);

          if (dist1_sq < dist2_sq) {
            // then current directors point towards each other, so deflect them
            /*psi_temp = psi[ci];
            psi[ci] = psi[cj];
            psi[cj] = psi_temp;*/
            psi[ci] += PI + (2 * drand48() - 1) * PI / 4;
            psi[ci] -= 2 * PI * round(psi[ci] / (2 * PI));
            psi[cj] += PI + (2 * drand48() - 1) * PI / 4;
            psi[cj] -= 2 * PI * round(psi[cj] / (2 * PI));

            activePropulsionFactor[ci] = 1.0;
            activePropulsionFactor[cj] = 1.0;
            counter++;  // indicate that a swap has occurred for cell ci
            polarizationCounter++;
          }
        }
      }
    }
    if (counter > 1) {
      std::cout << "deflected cell # " << ci << " " << counter
                << " times in the same timestep!\n";
    }
  }
}

double epi2D::getDistanceToVertexAtAnglePsi(int ci, double psi_ci, double cx, double cy, int& gi) {
  // given a cell ci and its director angle psi_ci, its com cx cy, return the distance to the nearest vertex at psi_ci and its index gi
  double dx, dy;
  std::vector<double> psi_i;
  double angular_displacement;
  for (int gj = gindex(ci, 0); gj < gindex(ci, 0) + nv[ci]; gj++) {
    dx = x[gj * NDIM] - cx;
    dy = x[gj * NDIM + 1] - cy;
    angular_displacement = atan2(dy, dx) - psi_ci;
    angular_displacement -= 2 * PI * round(angular_displacement / (2 * PI));
    psi_i.push_back(fabs(angular_displacement));
  }
  int vi = std::distance(psi_i.begin(), std::min_element(psi_i.begin(), psi_i.end()));  // argmin
  gi = vi + gindex(ci, 0);
  return sqrt(pow(x[gi * NDIM] - cx, 2) + pow(x[gi * NDIM + 1] - cy, 2));
}

double epi2D::getDistanceToRandomUnadheredVertex(int ci, double cx, double cy, int& gi) {
  // given cell index ci, and center cx cy, modify gi and psi[ci], and return the distance to a random unadhered vertex
  // unadhered is determined by vnn, which is modified by the attractive force routine
  int gj;
  double dx, dy;
  std::vector<int> listOfVoidVertices;
  for (int vi = 0; vi < nv[ci]; vi++) {
    gj = vi + gindex(ci, 0);
    if (vnn[gj][2] < 0) {
      // if vnn[gj][2] is negative, then vnn has no external neighbors, so it is a void-adjacent vertex
      listOfVoidVertices.push_back(gj);
    }
  }
  // select a vertex from this vector randomly
  if (listOfVoidVertices.size() > 0) {
    gi = listOfVoidVertices[rand() % listOfVoidVertices.size()];
    dx = x[gi * NDIM] - cx;
    dy = x[gi * NDIM + 1] - cy;
    psi[ci] = atan2(dy, dx);
    return sqrt(pow(dx, 2) + pow(dy, 2));
  } else
    return -1;
}

double epi2D::rotateAndCalculateArcLength(int ci, std::vector<int>& woundIndicesBelongingToCi) {
  // compute arc length of ci's contribution to the wound, and rotate woundIndicesBelongingToCi into sequential order according to ip1 or im1
  double arc_length = 0.0;
  for (int i = 0; i < woundIndicesBelongingToCi.size(); i++) {
    // check if ip1[woundIndicesBelongingToCi[j]] == woundIndicesBelongingToCi[j+1] for all j
    // if not, then do a cyclic permutation

    bool isOrderedRight = true;
    bool isOrderedLeft = true;
    for (int j = 0; j < int(woundIndicesBelongingToCi.size()) - 1; j++) {
      isOrderedRight *= (ip1[woundIndicesBelongingToCi[j]] == woundIndicesBelongingToCi[j + 1]);
      isOrderedLeft *= (im1[woundIndicesBelongingToCi[(j + woundIndicesBelongingToCi.size()) % woundIndicesBelongingToCi.size()]] == woundIndicesBelongingToCi[(j - 1 + woundIndicesBelongingToCi.size()) % woundIndicesBelongingToCi.size()]);
    }
    if (isOrderedRight || isOrderedLeft) {
      break;
    } else {
      std::rotate(woundIndicesBelongingToCi.begin(), woundIndicesBelongingToCi.begin() + 1, woundIndicesBelongingToCi.end());
    }
  }

  for (int i = 0; i < int(woundIndicesBelongingToCi.size()) - 1; i++) {
    arc_length += vertDistNoPBC(woundIndicesBelongingToCi[i], woundIndicesBelongingToCi[i + 1]);
  }
  return arc_length;
}

void epi2D::purseStringContraction(double B) {
  // updatePurseStringContacts();
  integratePurseString(B);  // evaluate forces on and due to purse-string, and integrate its position
  for (int psi = 0; psi < psContacts.size(); psi++) {
    l0_ps[psi] *= exp(-strainRate_ps * dt);
  }
}

void epi2D::initializePurseStringVariables() {
  // using sortedWoundIndices, establish a list of virtual purse-string particles
  int ci, vi;
  std::vector<double> x_ps_temp;
  std::vector<double> v_ps_temp;
  std::vector<double> F_ps_temp;
  std::vector<double> l0_ps_temp;
  std::vector<int> psContacts_temp;
  std::vector<bool> isSpringBroken_temp;
  int gi, gnext;
  for (int i = 0; i < sortedWoundIndices.size(); i++) {
    gi = sortedWoundIndices[i];
    gnext = sortedWoundIndices[(i + 1 + sortedWoundIndices.size()) % sortedWoundIndices.size()];
    x_ps_temp.push_back(x[gi * NDIM]);
    x_ps_temp.push_back(x[gi * NDIM + 1]);
    v_ps_temp.push_back(v[gi * NDIM]);
    v_ps_temp.push_back(v[gi * NDIM + 1]);
    F_ps_temp.push_back(F[gi * NDIM]);
    F_ps_temp.push_back(F[gi * NDIM + 1]);
    psContacts_temp.push_back(gi);
    l0_ps_temp.push_back(vertDistNoPBC(gi, gnext));
    isSpringBroken_temp.push_back(false);
    cindices(ci, vi, gi);
    if (std::find(initialWoundCellIndices.begin(), initialWoundCellIndices.end(), ci) == initialWoundCellIndices.end()) {
      initialWoundCellIndices.push_back(ci);  // push back unique ci's
      cout << "ci in initial wound indices = " << ci << '\n';
    }
  }
  x_ps = x_ps_temp;
  v_ps = v_ps_temp;
  F_ps = F_ps_temp;
  l0_ps = l0_ps_temp;
  psContacts = psContacts_temp;
  isSpringBroken = isSpringBroken_temp;
  if (x_ps.size() <= 1)
    cout << "WARNING: purse-string vector has size 1 or smaller, in initializePurseStringVariables.\n";
}

void epi2D::evaluatePurseStringForces(double B) {
  // using psContacts and the preferred length interaction, calculate forces on purse-string virtual particles
  // deltaSq is the breaking point for the virtual particle-wound particle interaction
  // k_ps is the wound-purseString spring force constant
  // B is the damping constant
  int gi, ipi, imi;
  double xp, yp, xw, yw;  // purse-string and wound coordinates
  bool cutoff;            // cutoff distance
  double dx, dy, fx, fy, l0i, l0im1;

  double fli, flim1, cx, cy, xi, yi, gip1, xip1, yip1;
  double rho0 = sqrt(a0.at(0));  // shouldn't have shape parameter for the wound. what's my alternative?
  double lim1x, lim1y, lix, liy, lip1x, lip1y, li, lim1, dli, dlim1;
  double rim1x, rim1y, rix, riy, rip1x, rip1y;

  // zero forces here, since no other function contributes to their values each time step
  std::fill(F_ps.begin(), F_ps.end(), 0);

  // compute center of mass for wound segment interaction
  xi = x_ps[0];
  yi = x_ps[1];
  cx = xi;
  cy = yi;
  for (int psi = 0; psi < psContacts.size() - 1; psi++) {
    gip1 = (psi + 1 + psContacts.size()) % psContacts.size();
    dx = x_ps[NDIM * gip1] - xi;
    dy = x_ps[NDIM * gip1 + 1] - yi;
    if (pbc[0])
      dx -= L[0] * round(dx / L[0]);
    if (pbc[1])
      dy -= L[1] * round(dy / L[1]);
    xip1 = xi + dx;
    yip1 = yi + dy;
    cx += xip1;
    cy += yip1;
    xi = xip1;
    yi = yip1;
  }
  cx /= psContacts.size();
  cy /= psContacts.size();

  for (int psi = 0; psi < psContacts.size(); psi++) {
    // compute forces due to the spring between wound and virtual particles
    gi = psContacts[psi];
    xp = x_ps[psi * NDIM];
    yp = x_ps[psi * NDIM + 1];
    xw = x[gi * NDIM];
    yw = x[gi * NDIM + 1];
    dx = xp - xw;
    dy = yp - yw;
    cutoff = ((dx * dx + dy * dy) < deltaSq * pow(r[0] * 2, 2));
    isSpringBroken[psi] = !cutoff;  // record broken springs for printouts
    fx = k_ps * dx * cutoff;
    fy = k_ps * dy * cutoff;
    F[gi * NDIM] += fx;
    F[gi * NDIM + 1] += fy;
    F_ps[psi * NDIM] -= fx;
    F_ps[psi * NDIM + 1] -= fy;

    // stress on gj should be the same as on gi, since it's the opposite separation and also opposite force
    fieldStress[gi][0] += -dx / 2 * fx;
    fieldStress[gi][1] += -dy / 2 * fy;
    fieldStress[gi][2] += -0.5 * (dx / 2 * fy + dy / 2 * fx);

    // update energy from purse-string
    U += 0.5 * fx * dx + 0.5 * fy * dy;
    U_ps += 0.5 * fx * dx + 0.5 * fy * dy;

    // compute forces due to preferred segment lengths between virtual particles
    // ipi, im1 are the next and previous virtual particle indices, respectively
    ipi = (psi + 1 + psContacts.size()) % psContacts.size();
    imi = (psi - 1 + psContacts.size()) % psContacts.size();
    l0i = l0_ps[psi];
    l0im1 = l0_ps[imi];

    rix = x_ps[NDIM * psi] - cx;
    riy = x_ps[NDIM * psi + 1] - cy;
    rip1x = x_ps[NDIM * ipi] - cx;
    rip1y = x_ps[NDIM * ipi + 1] - cy;
    if (pbc[0])
      rip1x -= L[0] * round(rip1x / L[0]);
    if (pbc[1])
      rip1y -= L[1] * round(rip1y / L[1]);

    rim1x = x_ps[NDIM * imi] - cx;
    rim1y = x_ps[NDIM * imi + 1] - cy;
    if (pbc[0])
      rim1x -= L[0] * round(rim1x / L[0]);
    if (pbc[1])
      rim1y -= L[1] * round(rim1y / L[1]);

    // Perimeter force
    lix = rip1x - rix;
    liy = rip1y - riy;

    lim1x = rix - rim1x;
    lim1y = riy - rim1y;

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
    fx = (fli * dli * lix / li) - (flim1 * dlim1 * lim1x / lim1);
    fy = (fli * dli * liy / li) - (flim1 * dlim1 * lim1y / lim1);
    F_ps[NDIM * psi] += fx;
    F_ps[NDIM * psi + 1] += fy;

    // choosing not to update potential energy of the purse-string itself
    // U += 0.5 * kl * (dli * dli);
  }
}

void epi2D::integratePurseString(double B) {
  // velocity verlet force update with damping
  // VV position update
  int virtualDOF = psContacts.size() * NDIM;
  for (int i = 0; i < virtualDOF; i++) {
    x_ps[i] += dt * v_ps[i] + 0.5 * dt * dt * F_ps[i];

    // recenter in box
    if (x_ps[i] > L[i % NDIM] && pbc[i % NDIM])
      x_ps[i] -= L[i % NDIM];
    else if (x_ps[i] < 0 && pbc[i % NDIM])
      x_ps[i] += L[i % NDIM];
  }

  // force update
  std::vector<double> F_ps_old = F_ps;
  evaluatePurseStringForces(B);

  // VV VELOCITY UPDATE #2
  for (int i = 0; i < virtualDOF; i++) {
    F_ps[i] -= (B * v_ps[i] + B * F_ps_old[i] * dt / 2);
    F_ps[i] /= (1 + B * dt / 2);
    v_ps[i] += 0.5 * (F_ps[i] + F_ps_old[i]) * dt;
  }
}

void epi2D::printConfiguration2D() {
  // overloaded to print out psi and other very specific quantities of interest
  // local variables
  int ci, cj, vi, gi, ctmp, zc, zv;
  double xi, yi, dx, dy, Lx, Ly;

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
  Lx = L[0];
  Ly = L[1];

  // print information starting information
  posout << setw(w) << left << "NEWFR"
         << " " << endl;
  posout << setw(w) << left << "NUMCL" << setw(w) << left << NCELLS << endl;
  posout << setw(w) << left << "PACKF" << setw(wnum) << setprecision(pnum)
         << left << vertexPackingFraction2D() << endl;

  // print box sizes
  posout << setw(w) << left << "BOXSZ";
  posout << setw(wnum) << setprecision(pnum) << left << Lx;
  posout << setw(wnum) << setprecision(pnum) << left << Ly;
  posout << endl;

  // print stress info
  posout << setw(w) << left << "STRSS";
  posout << setw(wnum) << setprecision(pnum) << left << -stress[0];
  posout << setw(wnum) << setprecision(pnum) << left << -stress[1];
  posout << setw(wnum) << setprecision(pnum) << left << -stress[2];
  posout << setw(wnum) << setprecision(pnum) << left << L[0] / initialLx - 1;
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
    posout << setw(wnum) << left << psi[ci];
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
    posout << setw(wnum) << left << flag[ci];
    posout << setw(wnum) << left << flagPos[ci][0];
    posout << setw(wnum) << left << flagPos[ci][1];
    posout << endl;

    // get initial vertex positions
    gi = gindex(ci, 0);
    xi = x[NDIM * gi];
    yi = x[NDIM * gi + 1];

    // place back in box center
    if (pbc[0])
      xi = fmod(xi, Lx);
    if (pbc[1])
      yi = fmod(yi, Ly);

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
        dx -= Lx * round(dx / Lx);
      xi += dx;

      dy = x[NDIM * gi + 1] - yi;
      if (pbc[1])
        dy -= Ly * round(dy / Ly);
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

  cout << "in printPos, writing to purseout!\n";

  if (psContacts.size() != 0) {
    for (int j = 0; j < psContacts.size(); j++) {
      // print position of purse-string virtual vertex and its bonded real wound vertex
      double xp = x_ps[j * NDIM];
      double yp = x_ps[j * NDIM + 1];
      double xw = x[psContacts[j] * NDIM];
      double yw = x[psContacts[j] * NDIM + 1];

      // if bond is unbroken and ps is still active, print a line from virtual ps vertex to wound vertex.
      if (!isSpringBroken[j] && deltaSq + k_ps > 1e-10) {
        purseout << xp << '\t' << yp << '\n';
        purseout << xw << '\t' << yw << '\n';
      }
    }
  }
  purseout << "*EOB\n";

  cout << "leaving printPos\n";
}