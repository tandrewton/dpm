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

double epi2D::vertDistSqNoPBC(int gi, int gj) {
  // not valid for PBC
  double distsq = 0;
  distsq = pow(x[gi * NDIM] - x[gj * NDIM], 2) +
           pow(x[gi * NDIM + 1] - x[gj * NDIM + 1], 2);
  return distsq;
}

std::vector<int> epi2D::regridSegment(int wVertIndex, double vrad) {
  // need to subtract indices from giConnectedToFlag if using this function in the future
  // return value is the indices of the deleted vertices for proper gi counting later

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
    cout << "why does orderedWVerts.size() = 0?\n";
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
  for (int i = 0; i < poly_bd_x.size(); i++) {
    evaluatePolygonalWallForces(poly_bd_x[i], poly_bd_y[i]);
  }
}

void epi2D::attractiveForceUpdateWithPolyWall() {
  resetForcesAndEnergy();
  shapeForces2D();
  vertexAttractiveForces2D_2();
  for (int i = 0; i < poly_bd_x.size(); i++) {
    evaluatePolygonalWallForces(poly_bd_x[i], poly_bd_y[i]);
  }
}

void epi2D::vertexAttractiveForces2D_2() {
  // altered from dpm attractive force code, because it works with larger l2
  // values. (warning: probably won't work with bending.) local variables
  int ci, cj, gi, gj, vi, vj, bi, bj, pi, pj;
  double sij, rij, dx, dy, rho0;
  double ftmp, utmp, fx, fy;

  // attraction shell parameters
  double shellij, cutij, xij, kint = (kc * l1) / (l2 - l1);

  // sort particles
  sortNeighborLinkedList2D();
  // cout << "r[0] = " << r[0] << '\n';

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
                  ftmp = kc * (1 - (rij / sij)) * (1 / sij);
                  utmp = 0.5 * kc * pow((1 - (rij / sij)), 2.0);
                  cellU[ci] += utmp;
                  U += utmp;
                } else
                  ftmp = 0;
              } else if (rij > cutij) {
                // force scale
                ftmp = kint * (xij - 1.0 - l2) / sij;
                utmp = -0.5 * kint * pow(1.0 + l2 - xij, 2.0);
                U += utmp;
                cellU[ci] += utmp / 2.0;
                cellU[cj] += utmp / 2.0;
              } else {
                // force scale
                ftmp = kc * (1 - xij) / sij;
                utmp = 0.5 * kc * (pow(1.0 - xij, 2.0) - l1 * l2);
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
                    ftmp = kc * (1 - (rij / sij)) * (1 / sij);
                    utmp = 0.5 * kc * pow((1 - (rij / sij)), 2.0);
                    cellU[ci] += utmp;
                    U += utmp;
                  } else
                    ftmp = 0;
                } else if (rij > cutij) {
                  // force scale
                  ftmp = kint * (xij - 1.0 - l2) / sij;
                  utmp = -0.5 * kint * pow(1.0 + l2 - xij, 2.0);
                  U += utmp;
                  cellU[ci] += utmp / 2.0;
                  cellU[cj] += utmp / 2.0;
                } else {
                  // force scale
                  ftmp = kc * (1 - xij) / sij;
                  utmp = 0.5 * kc * (pow(1.0 - xij, 2.0) - l1 * l2);
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

void epi2D::circuloLineAttractiveForces() {
  // altered from vertexAttractiveForces2D_2, here we use vertex-vertex and vertex-line segment distances to make a smooth interaction
  // models sliding adhesion and repulsion.
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

void epi2D::calculateSmoothInteraction(double& rx, double& ry, double& sij, double& shellij, double& cutij, double& kint, double& kc, int& gi, int& gj, double& projection, int& ci, int& cj) {
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

        // add to contacts
        for (int i = 0; i < vnn[gi].size(); i++) {
          if (ci == cj)
            break;

          // gj should probably be "middle" here
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

void epi2D::attractiveForceUpdate_2() {
  resetForcesAndEnergy();
  shapeForces2D();
  vertexAttractiveForces2D_2();
}

void epi2D::attractiveForceUpdate_circulo() {
  resetForcesAndEnergy();
  shapeForces2D();
  // std::vector<double> shapeForce(NDIM*NVTOT,0.0);
  // shapeForce = F;

  circuloLineAttractiveForces();
}

void epi2D::substrateadhesionAttractiveForceUpdate(bool isCirculoLine) {
  // compute forces for shape, attractive, and substrate adhesion contributions
  int gi = 0, argmin, flagcount = 0;
  int numVerticesAttracted = 1;
  double dx, dy;
  double fx, fy;

  // reset forces, then get shape and attractive forces.
  if (isCirculoLine)
    attractiveForceUpdate_circulo();
  else
    attractiveForceUpdate_2();

  resetActiveEnergy();
  // directorDiffusion();
  updateSubstrateSprings();

  for (auto ci : initialWoundCellIndices) {
    std::vector<double> distSqVertToFlag(nv[0], 0.0);
    // check for protrusions (unimpeded flag toss, and flag location inside box)
    if (flag[ci]) {
      // find nearest vertices
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

void epi2D::crawlingWithPurseString() {
  bool isCirculoLine = false;
  substrateadhesionAttractiveForceUpdate(isCirculoLine);
}

void epi2D::crawlingWithPurseStringAndCircularWalls() {
  bool attractionOn = true;
  substrateadhesionAttractiveForceUpdate();
  for (int i = 0; i < poly_bd_x.size(); i++) {
    evaluatePolygonalWallForces(poly_bd_x[i], poly_bd_y[i], attractionOn);
  }
}

void epi2D::crawlingWithPurseStringCirculo() {
  bool isCirculoLine = true;
  substrateadhesionAttractiveForceUpdate(isCirculoLine);
}

void epi2D::crawlingWithPurseStringCirculoWalls() {
  bool isCirculoLine = true;
  bool attractionOn = true;
  substrateadhesionAttractiveForceUpdate(isCirculoLine);
  for (int i = 0; i < poly_bd_x.size(); i++) {
    evaluatePolygonalWallForces(poly_bd_x[i], poly_bd_y[i], attractionOn);
  }
}

void epi2D::circuloLineAttractionWithCircularWalls() {
  attractiveForceUpdate_circulo();
  bool attractionOn = true;
  for (int i = 0; i < poly_bd_x.size(); i++) {
    evaluatePolygonalWallForces(poly_bd_x[i], poly_bd_y[i], attractionOn);
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

void epi2D::dampedNVETest(dpmMemFn forceCall, double T, double dt0, int NT, int NPRINTSKIP) {
  // local variables
  int t, i;
  double K, simclock, B = 1.0;

  // set time step magnitude
  setdt(dt0);

  // loop over time, print energy
  for (t = 0; t < NT; t++) {
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

      if (i / 2 >= szList[0] && i / 2 < szList[1] && i % 2 == 0 && t == 0) {
        v[i] += 1.0;
        cout << 0.5 * (F[i] + F_old[i]) * dt << '\n';
      }
    }

    // update sim clock
    simclock += dt;

    // print to console and file
    if (t % NPRINTSKIP == 0) {
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
      cout << "	** t / NT	= " << t << " / " << NT << endl;
      cout << "	** U 		= " << setprecision(12) << U << endl;
      cout << "	** K 		= " << setprecision(12) << K << endl;
      cout << "	** E 		= " << setprecision(12) << U + K << endl;

      // print to configuration only if position file is open
      if (posout.is_open())
        printConfiguration2D();
    }
  }
}

// for testing numerical stability
void epi2D::vertexNVE(dpmMemFn forceCall, double dt0, int NT, int NPRINTSKIP) {
  // local variables
  int t, i;
  double K;

  // set time step magnitude
  setdt(dt0);

  // initialize time keeper
  simclock = 0.0;

  // loop over time, print energy
  for (t = 0; t < NT; t++) {
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

    if (NPRINTSKIP != 0) {
      if (t % NPRINTSKIP == 0) {
        //                         compute kinetic energy
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
        cout << "	** t / NT	= " << t << " / " << NT << endl;
        cout << " ** simclock = " << setprecision(12) << simclock << endl;
        cout << "	** U 		= " << setprecision(12) << U << endl;
        cout << "	** K 		= " << setprecision(12) << K << endl;
        cout << "	** E 		= " << setprecision(12) << U + K << endl;

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

void epi2D::dampedNVE2D(dpmMemFn forceCall, double dt0, double duration, double printInterval) {
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
    }

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

void epi2D::dampedCompression(dpmMemFn forceCall, double dt0, double duration, double printInterval) {
  // make sure velocities exist or are already initialized before calling this
  int i;
  double K, t0 = simclock;
  double temp_simclock = simclock;

  // set time step magnitude
  setdt(dt0);
  int NPRINTSKIP = printInterval / dt;

  // initial coordinate of walls
  double lowerWallPos = 0, upperWallPos = L[1], leftWallPos = -1e10, rightWallPos = 1e10;
  double comx, comy;
  com2D(0, comx, comy);
  displaceCell(0, L[0] / 2 - comx, L[1] / 2 - comy);
  // center the particle in the middle of the walls so the compression is symmetric

  std::ofstream wallout("wallPositions.txt");
  // print wall positions for easy debugging of wall locations

  // loop over time, print energy
  while (simclock - t0 < duration) {
    if (int((simclock - t0) / dt) % 100 == 0) {
      wallout << simclock - t0 << '\t' << lowerWallPos << '\t' << upperWallPos << '\t' << leftWallPos << '\t' << rightWallPos << '\n';
      double increment = r[0] * 0.025;
      if (simclock - t0 < 0.1 * duration) {
        // bring the upper and lower walls toward each other
        upperWallPos -= increment;
        lowerWallPos += increment;
      } else if (simclock - t0 < 0.6 * duration) {
        // hold the wall at max displacement for a little bit
      } else if (simclock - t0 >= 0.6 * duration) {
        // release the wall
        upperWallPos += r[0];
        lowerWallPos -= r[0];
      }
    }

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

    // compute additional wall forces
    computeWallForce(lowerWallPos, upperWallPos, leftWallPos, rightWallPos);

    // VV VELOCITY UPDATE #2
    for (i = 0; i < vertDOF; i++) {
      F[i] -= (B * v[i] + B * F_old[i] * dt / 2);
      F[i] /= (1 + B * dt / 2);
      v[i] += 0.5 * (F[i] + F_old[i]) * dt;
    }

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

// simulation without boundaries, i.e. 0 pressure simulation
void epi2D::dampedNP0(dpmMemFn forceCall, double dt0, double duration, double printInterval, int purseStringOn, double relaxTime) {
  // make sure velocities exist or are already initialized before calling this
  // assuming zero temperature - ignore thermostat (not implemented)
  int i;
  double K, t0 = simclock, shape_ci;
  double temp_simclock = simclock;
  double initialWoundArea = 1e10, healingTime = NAN;
  bool alreadyRecordedFinalCells = false;

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

    if (simclock - t0 > relaxTime && purseStringOn == 1) {  // purseStringOn refers to whether it's been initialized, not its parameters. so dsq = 0 has a nonfunctional pursestring, but still has purseStringOn = 1
      if (psContacts.size() == 0 && std::isnan(woundArea)) {
        cout << "inside psContacts.size() == 0 and woundArea == NAN case, which should only occur once!\n";
        getWoundVertices(nthLargestCluster);

        if (sortedWoundIndices.size() < nv[0]) {
          cout << "getWoundVertices has returned a cluster smaller than a single cell. print 0th cluster, see if it works\n";
          getWoundVertices(0);
          cout << sortedWoundIndices.size() << '\n';
        }

        woundCenterX = 0;
        woundCenterY = 0;
        for (auto gi : sortedWoundIndices) {
          woundCenterX += x[gi * NDIM];
          woundCenterY += x[gi * NDIM + 1];
          cout << "x,y of wound vertex gi = " << gi << " is " << x[gi * NDIM] << '\t' << x[gi * NDIM + 1] << '\n';
        }
        woundCenterX /= sortedWoundIndices.size();
        woundCenterY /= sortedWoundIndices.size();

        cout << "wound center before calculateWoundArea (initial wound detection) = " << woundCenterX << '\t' << woundCenterY << '\n';
        woundArea = calculateWoundArea(woundCenterX, woundCenterY);
        vout << simclock - t0 << '\t' << woundArea << '\t' << purseStringTension << '\t' << purseStringTransmittedTension << '\n';
        // cout << simclock - t0 << '\t' << woundArea << '\n';
        cout << "wound center after calculateWoundArea (initial wound detection) = " << woundCenterX << '\t' << woundCenterY << '\n';
        initialWoundArea = woundArea;
        assert(!std::isnan(initialWoundArea));
        initializePurseStringVariables();
      }

      // get max and min of x coords of purse-string; if max-min is near zero, then purse-string should be dissolved
      double max_ps = x_ps[0], min_ps = x_ps[0];
      double min_allowed_ps_length = 0.01 * r[0];
      // double min_allowed_ps_length = 0.0;
      for (int psi = 0; psi < x_ps.size(); psi += 2) {
        if (x_ps[psi] < min_ps)
          min_ps = x_ps[psi];
        if (x_ps[psi] > max_ps)
          max_ps = x_ps[psi];
      }

      if (deltaSq == 0.0) {
        // do nothing - ignore pursestring
      } else if (max_ps - min_ps > min_allowed_ps_length) {
        purseStringContraction();
      } else {                                      // purse-string can't shrink anymore
        if (isPurseStringDoneShrinking == false) {  // if purse-string can't shrink anymore, set this flag to true and stop straining the pursestring
          cout << "disassembling purse-string, lengths are too short to shrink further \n";
          isPurseStringDoneShrinking = true;
          strainRate_ps = 0.0;
          cout << "max_ps = " << max_ps << '\t' << ", min_ps = " << min_ps << '\n';
          cout << "simclock = " << simclock << '\n';
          assert(false);
        }
        purseStringContraction();
      }

      if (int(simclock / dt) % 50 == 0) {
        // choice of woundArea cutoff is 4 vertex areas. % of initial wound area is unreliable, 10 vertex areas is too many, 1 vertex area will lead to issues with void segmentation
        woundAreaCutoffEndSimulation = 4 * 0.5 * PI * r[0] * r[0];
        woundArea = calculateWoundArea(woundCenterX, woundCenterY);
        bool oldWoundSeedsTooSmall = false;
        double newAreaFraction = fabs(woundArea - previousWoundArea) / previousWoundArea;
        while (std::isnan(woundArea) || woundArea < woundAreaCutoffEndSimulation || (newAreaFraction > 0.5 && previousWoundArea > woundAreaCutoffEndSimulation)) {
          cout << "in while loop, detected wound area is nan or very small, or has a discontinuous-looking drop in area above the cutoff criterion (which is relevant for determining healing time). repeating wound area calculation with a stored wound seed from previous configurations. \n oldWoundLocations.size() = " << oldWoundLocations.size() << "\n";
          cout << "current (proposed) wound area = " << woundArea << ", previous wound area = " << previousWoundArea << '\n';
          // every iteration, remove an element from oldWoundLocations. terminate when woundArea is valid, or when oldWoundLocations is empty
          if (oldWoundLocations.size() == 0) {
            cout << "reached end of oldWoundLocations.. wound area is still nan or below simulation end threshold.\n";
            break;
          }
          woundCenterX = oldWoundLocations[0][0];
          woundCenterY = oldWoundLocations[0][1];
          cout << "potential wound location: " << woundCenterX << '\t' << woundCenterY << '\n';
          // do not record old wound points, since we're looping over them
          woundArea = calculateWoundArea(woundCenterX, woundCenterY, false);
          newAreaFraction = fabs(woundArea - previousWoundArea) / previousWoundArea;

          if (oldWoundLocations.size() < 50)
            oldWoundLocations.erase(oldWoundLocations.begin());
          else
            oldWoundLocations.erase(oldWoundLocations.begin(), oldWoundLocations.begin() + int(0.05 * oldWoundLocations.size()));

          if (!std::isnan(woundArea) && woundArea < woundAreaCutoffEndSimulation)
            oldWoundSeedsTooSmall = true;  // non-NAN values found for area, but were too small and were discarded
        }
        if (previousWoundArea < woundAreaCutoffEndSimulation || (std::isnan(woundArea) && oldWoundSeedsTooSmall)) {
          // if previous area is small and current area is small, then end simulation. if current area is NaN and alternative areas were small, then end simulation.
          woundArea = 0;
        } else if (std::isnan(woundArea) && !oldWoundSeedsTooSmall) {
          cout << "woundArea is nan! exiting\n";
          assert(!std::isnan(woundArea));
        }

        vout << simclock - t0 << '\t' << woundArea << '\t' << purseStringTension << '\t' << purseStringTransmittedTension << '\n';
        previousWoundArea = woundArea;

        // write shape information to files
        innerout << simclock - t0 << '\t';
        bulkout << simclock - t0 << '\t';
        for (int ci = 0; ci < NCELLS; ci++) {
          shape_ci = pow(perimeter(ci), 2) / (4 * PI * area(ci));
          // if ci is an initial wound-edge cell
          /*if (std::find(initialWoundCellIndices.begin(), initialWoundCellIndices.end(), ci) != initialWoundCellIndices.end()) {
            innerout << shape_ci << '\t';
          } else {
            bulkout << shape_ci << '\t';
          }*/
          innerout << shape_ci << '\t';
          bulkout << shape_ci << '\t';
        }
        innerout << '\n';
        bulkout << '\n';

        if (woundArea < 0.05 * initialWoundArea && std::isnan(healingTime)) {
          cout << "wound area is less than 5 percent of initial!\n";
          cout << simclock - t0 << '\t' << woundArea << '\n';
          healingTime = simclock - t0;
        }
      }
    }

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
          int nthLargestCluster = 1;
          printConfiguration2D();
          printBoundaries(nthLargestCluster);
          cerr << "done printing in NP0\n";
        }

        cerr << "Number of polarization deflections: " << polarizationCounter
             << '\n';
      }
    }
    // do things after all the simulation steps have been taken, or if the wound is closed (area=0 or area= 4 vertex areas)
    if (duration > 100 && (simclock - t0 > duration - dt || woundArea == 0.0 || woundArea < woundAreaCutoffEndSimulation) && !std::isnan(healingTime) && !alreadyRecordedFinalCells) {
      // exit conditions
      cout << "woundAreaCutoff to end simulation = " << woundAreaCutoffEndSimulation << '\n';
      cout << "final calculated wound center = " << woundCenterX << '\t' << woundCenterY << '\n';
      cout << "woundArea, simclock - t0, duration-dt = " << woundArea << " " << simclock - t0 << " " << duration - dt << '\n';
      break;
    }
  }

  if (duration > 100) {  // production run, not equilibration run
    // to get here, simulation must have (time > duration) or (wound is detected to be closed)
    cout << "breaking from wound healing integration, recording final quantities.\n";

    int ci, vi;
    double fivePercentWoundArea_radius_sq = 0.05 * initialWoundArea / PI;
    double distance;
    std::vector<int> finalCellsInCenterOfWound;

    for (int gi = 0; gi < NVTOT; gi++) {
      distance = pow(x[gi * NDIM] - woundCenterX, 2) + pow(x[gi * NDIM + 1] - woundCenterY, 2);
      if (distance < fivePercentWoundArea_radius_sq) {  // if vertex is within radius set by 5% of wound area, record the cell
        cindices(ci, vi, gi);
        finalCellsInCenterOfWound.push_back(ci);
      }
    }
    // store the final cells in the center of the wound
    std::sort(finalCellsInCenterOfWound.begin(), finalCellsInCenterOfWound.end());
    finalCellsInCenterOfWound.erase(std::unique(finalCellsInCenterOfWound.begin(), finalCellsInCenterOfWound.end()), finalCellsInCenterOfWound.end());

    for (auto i : finalCellsInCenterOfWound)
      cout << "final cell : " << i << '\n';

    // if ci is an initial wound cell, record it as in the first row. if not, record it as in the second row.
    std::vector<int> firstRow, secondRow;
    for (int ci = 0; ci < NCELLS; ci++) {
      if (std::find(initialWoundCellIndices.begin(), initialWoundCellIndices.end(), ci) != initialWoundCellIndices.end())
        firstRow.push_back(ci);
      else
        secondRow.push_back(ci);
    }

    // record each cellID, whether it's initially wound adjacent (firstRow) or not, and whether its finally wound adjacent (finalCellsInCenterOfWound) or not.
    for (int ci = 0; ci < NCELLS; ci++) {
      bool isInInitialWound = 0, isInFinalWound = 0;

      isInInitialWound = std::find(initialWoundCellIndices.begin(), initialWoundCellIndices.end(), ci) != initialWoundCellIndices.end();

      isInFinalWound = std::find(finalCellsInCenterOfWound.begin(), finalCellsInCenterOfWound.end(), ci) != finalCellsInCenterOfWound.end();

      cellIDout << ci << '\t' << isInInitialWound << '\t' << isInFinalWound << '\n';
    }

    cout << "number of cells in center = " << finalCellsInCenterOfWound.size() << '\n';
    cout << "for 5% wound area of radius " << sqrt(fivePercentWoundArea_radius_sq) << '\n';
    if (!std::isnan(healingTime)) {
      cout << "healingTime = " << healingTime << '\n';
      cout << "wound center = " << woundCenterX << '\t' << woundCenterY << '\n';
    } else {
      cout << "wound did not close.\n";
      cout << "psContacts : ";
      healingTime = duration;
    }
    woundPropertiesout << healingTime << '\t' << finalCellsInCenterOfWound.size() << '\n';
    alreadyRecordedFinalCells = true;
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

void epi2D::computeWallForce(double lowerWall, double upperWall, double leftWall, double rightWall) {
  // simple routine to compute the forces between all vertices and the coordinates making up a hollow rectangular wall
  double sij, rij, ftmp, rho0 = sqrt(a0[0]);
  double dx, dy, fx, fy, wall;
  double forceMultiplier = 10;  // make walls more impenetrable than ordinary particles
  for (int i = 0; i < vertDOF; i++) {
    sij = 2 * r[floor(i / 2)];
    // calculate separations between particle degree of freedom and each wall.
    for (int j = 0; j < 2; j++) {  // for each degree of freedom, check the two walls (behind and in front of) it
      if (i % 2 == 0) {
        if (j == 0)
          wall = leftWall;
        else
          wall = rightWall;
      } else {
        if (j == 0)
          wall = lowerWall;
        else
          wall = upperWall;
      }
      dx = x[i] - wall;
      rij = fabs(dx);

      if (rij < sij) {
        ftmp = kc * forceMultiplier * (1 - (rij / sij)) * (rho0 / sij);
        fx = ftmp * (dx / rij);
        F[i] -= fx;
      }
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

    // we measure distances from the center of the box L/2, L/2, so args=0,0 is really xloc,yloc=L/2,L/2
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
  psi, and possibly others

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
  vl0.erase(vl0.begin() + deleteIndexGlobal,
            vl0.begin() + deleteIndexGlobal + numVertsDeleted);
  Fl0.erase(Fl0.begin() + deleteIndexGlobal,
            Fl0.begin() + deleteIndexGlobal + numVertsDeleted);
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
  cout << "\ndeleting vertices!\n\n simclock = " << simclock
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
    vl0.erase(vl0.begin() + i);
    Fl0.erase(Fl0.begin() + i);
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

  // If any vertices have both of their same-cell neighbors in voidFacingVertexIndices, let them also be in the list to rescue certain edge cases
  for (int gi = 0; gi < NVTOT; gi++) {
    // if gi is not in the voidFacingVertexIndices list
    if (std::find(voidFacingVertexIndices.begin(), voidFacingVertexIndices.end(), gi) == voidFacingVertexIndices.end()) {
      // if same-cell neighbors are in voidFacingVertexIndices
      int left = im1[gi];
      int right = ip1[gi];
      auto indOfLeft = std::find(voidFacingVertexIndices.begin(), voidFacingVertexIndices.end(), left);
      if (indOfLeft != voidFacingVertexIndices.end()) {
        auto indOfRight = std::find(voidFacingVertexIndices.begin(), voidFacingVertexIndices.end(), right);
        if (indOfRight != voidFacingVertexIndices.end()) {
          // then gi is probably mistakenly excluded from voidFacingVertexIndices, so add it.
          int offset;
          if (gi < right) {
            // then right is gi+1
            // so we want to insert gi into the location preceding right
            // i.e. ip1 = 19, gi = 0, ip1 = 1
            offset = 0;
            voidFacingVertexIndices.insert(indOfRight, gi);
          } else {
            // right is gi + 1 - nv[ci]
            // so we want to insert gi into the location after left
            // i.e. ip1 = 18, gi = 19, im1 = 0
            voidFacingVertexIndices.insert(indOfLeft + 1, gi);
          }
        }
      }
    }
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
  cout << "woundArea = " << woundArea << ", vertex area = " << 0.5 * PI * r[0] * r[0] << '\n';
  cout << "woundCenterX, woundCenterY = " << woundCenterX << '\t' << woundCenterY << '\n';
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

    for (auto i : vnn[current_vertex]) {
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
              cout << "before seeking alternate route, unwrapped_x_gi = ";
              for (auto unwrappedIt : unwrapped_x_gi) {
                cout << unwrappedIt << '\t';
              }
              cout << '\n';

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
  // used to get the initial wound vertices for initialization of the purse-string

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

  for (int gi = 0; gi < NVTOT; gi++) {
    // for debugging purposes
    if (findRoot(gi, ptr) == mode) {
      cout << "for debugging: indices in big cluster = " << gi << '\n';
    }
  }

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

    for (auto i : vnn[current_vertex]) {
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
  cout << "temporaryWoundIndices size = " << temporaryWoundIndices.size() << '\n';
  cout << "previous_vertex size = " << previous_vertex.size() << '\n';
  if (checkWoundClosedPolygon(temporaryWoundIndices)) {
    sortedWoundIndices = temporaryWoundIndices;
  } else {
    cout << "temporaryWoundIndices size = " << temporaryWoundIndices.size() << '\n';
    for (auto i : temporaryWoundIndices)
      cout << i << '\t';
    cout << '\n';
    cout << "wound was not a closed polygon, not returning any values for sortedWoundIndices, program should terminate soon\n";
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

double epi2D::calculateWoundArea(double& woundPointX, double& woundPointY, bool recordOldWoundPoints) {
  // input: a point previously inside the wound.
  //  we will nucleate the area calculation around this point
  // this algorithm gives the area of a wound by dividing up the simulation box into a grid
  //  each grid point is checked for being within any cell and for proximity to cell vertices
  //  if either is true, give that point a value of 1 (true)
  // then take woundPoint and calculate the largest continuous blob of 0 values in occupancyMatrix.
  //  multiply by grid point area to get the total wound area.

  // note: NAN will be ignored in plots, so I use it here for plotting purposes. Make sure that woundArea does not collide with important things (i.e. get multiplied into meaningful simulation quantities)
  // cout << "calculating woundArea, simclock = " << simclock << '\n';

  // note: running this every 100 or 1000 timesteps should be fine.
  std::vector<double> posX(NVTOT / 2), posY(NVTOT / 2);
  double iResolution, jResolution;  // stores scaled coordinates to pass as reference to isPointInPolygons
  for (int i = 0; i < NVTOT; i++) {
    if (i % 2 == 0)
      posX[i / 2] = x[NDIM * i];
    else
      posY[(i - 1) / 2] = x[NDIM * (i - 1) + 1];
  }
  // resolution is the unit box length that we'll use to map occupancyMatrix to real coordinates
  // double resolution = r[0] / 2.0 / 2.0;
  double resolution = r[0] / 2;
  if (resolution <= 1e-10)
    cerr << "bug: resolution in calculateWoundArea is zero\n";

  double xLow = *std::min_element(posX.begin(), posX.end());
  double xHigh = *std::max_element(posX.begin(), posX.end());
  double yLow = *std::min_element(posY.begin(), posY.end());
  double yHigh = *std::max_element(posY.begin(), posY.end());
  if (fabs(woundPointX) > xHigh) {
    cout << "woundPoint does not lie within the bounds of the simulation box, probably failed to find the center of a wound\n returning NAN area, skipping.\n";
    return NAN;
  }

  int xResolution = (xHigh - xLow) / resolution;
  int yResolution = (yHigh - yLow) / resolution;
  std::vector<std::vector<bool>> occupancyMatrix(xResolution, std::vector<bool>(yResolution, 0));
  for (int i = 0; i < xResolution; i++) {
    iResolution = i * resolution;
    for (int j = 0; j < yResolution; j++) {
      jResolution = j * resolution;
      // first pass: occupancy is 1 if inside a cell, 0 if not inside a cell (pnpoly point inclusion in polygon test algorithm)
      // note: area calculation occurs around xResolution^2 times = ~ 1-10 million times per call to calculateWoundArea
      occupancyMatrix[i][j] = !isPointInPolygons(iResolution, jResolution);

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

  //  check if the given point is not within the wound
  if (occupancyMatrix[woundPointXIndex][woundPointYIndex] == 1) {  // 1 is not in wound, 0 is in wound
    // given point is not in the wound, so we need to look for a nearby point that's hopefully in the wound
    int searchRange = 5, offset = 0;
    // look up to searchRange boxes away.
    int newXIndex = woundPointXIndex, newYIndex = woundPointYIndex;
    while (occupancyMatrix[newXIndex][newYIndex] == 1 && offset <= searchRange) {
      offset++;
      int rightIndex = woundPointXIndex + offset;
      int leftIndex = woundPointXIndex - offset;
      int aboveIndex = woundPointYIndex + offset;
      int belowIndex = woundPointYIndex - offset;
      std::vector<int> xIndex = {leftIndex, woundPointXIndex, rightIndex};
      std::vector<int> yIndex = {belowIndex, woundPointYIndex, aboveIndex};
      // using all combinations of xIndex and yIndex, e.g. could search diagonally (4 directions) or in all 8 directions
      for (int xii = 0; xii < 3; xii++) {
        for (int yii = 0; yii < 3; yii++) {
          if (occupancyMatrix[xIndex[xii]][yIndex[yii]] == 0) {
            newXIndex = xIndex[xii];
            newYIndex = yIndex[yii];
          }
        }
      }
    }

    if (offset > searchRange) {
      cout << "failed to find a nearby point identifiable as a wound, returning NAN for area\n";
      return NAN;
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
        // make sure row and column value is valid
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

  // now we've assigned 0 to all non-main-cluster boxes, so just sum over the number of 1s to get the number of boxes in main cluster
  double sum = 0.0;
  double xLocs = 0.0, yLocs = 0.0;
  for (int i = 0; i < labels.size(); i++) {
    for (int j = 0; j < labels[i].size(); j++) {
      sum += labels[i][j];
      xLocs += i * resolution * labels[i][j];
      yLocs += j * resolution * labels[i][j];
    }
  }

  // coarseness controls how many oldWoundLocations we store. currently using sqrt(area)/2 which gives a small fraction of the total number of boxes in each dimension
  int coarseness = sqrt(sum) / 2.0;
  if (coarseness < 1)
    coarseness = 1;
  if (recordOldWoundPoints && woundArea > woundAreaCutoffEndSimulation && sum > 100 && sum * pow(resolution, 2) > woundAreaCutoffEndSimulation) {
    // only clear rescue conditions if wound area is larger than needed. if it's less than needed, I should be keeping the old rescue condition, because sum will be small, so I won't have any stored points to remember
    oldWoundLocations.clear();
  }

  if (recordOldWoundPoints) {
    // coarseness must be > 0 or leads to floating point exception for modulo operation
    for (int i = 0; i < labels.size(); i++) {
      for (int j = 0; j < labels[i].size(); j++) {
        if (labels[i][j])
          if (i % coarseness == 0 && j % coarseness == 0)
            oldWoundLocations.push_back({i * resolution, j * resolution});
        // if this is a point in the wound, consider saving it
      }
    }
  }

  // cout << "inside calculateWoundArea, recording=" << recordOldWoundPoints << ", sum = " << sum << ", oldWoundLocations stored " << oldWoundLocations.size() << " locations\n";

  currentWoundIndices.clear();
  bool iAndjAreInDebug = false;
  int woundSearchRange = r[0] / resolution + 2;  // r[0] / resolution gives # boxes representing a vertex radius. + 1 allows it to search for immediate contacts to the vertex radius. + 2 gives buffer in edge case (vertex just barely passes edge of box). not totally sure why +1 isn't enough, but +2 seems to be the right amount to find the void
  int woundSearchRange_debug = woundSearchRange + 2;
  for (int gi = 0; gi < NVTOT; gi++) {
    double dist;
    int ci, vi;
    int i, j, sum_labels = 0;
    cindices(ci, vi, gi);
    if (std::find(initialWoundCellIndices.begin(), initialWoundCellIndices.end(), ci) != initialWoundCellIndices.end()) {
      // cell is indeed one of the initially wound-adjacent cells
      i = x[NDIM * gi] / resolution;
      j = x[NDIM * gi + 1] / resolution;
      if (i < 0)
        i = 0;
      if (j < 0)
        j = 0;
      if (i >= labels.size())
        i = labels.size() - 1;
      if (j >= labels[0].size())
        j = labels[0].size() - 1;

      if (labels[i][j] == 1) {
        cout << "location of bad label = " << i * resolution << '\t' << j * resolution << '\n';
        assert(labels[i][j] == 0);
      }
      // count up all neighbors of i,j. they are 1-valued if in wound, 0 otherwise
      for (int xInd = -woundSearchRange; xInd < woundSearchRange; xInd++) {
        for (int yInd = -woundSearchRange; yInd < woundSearchRange; yInd++) {
          if (i + xInd >= 0 && j + yInd >= 0 && i + xInd < labels.size() && j + yInd < labels[0].size())  // make sure it's in bounds
            sum_labels += labels[i + xInd][j + yInd];
        }
      }
      if (sum_labels > 0) {
        currentWoundIndices.push_back(gi);
      }
    }
  }

  //   area should be the number of boxes times the box area
  //   alternatively, I could get the exact area by locating vertices on the edge of my newly segmented void area, but that's more computational work
  return sum * pow(resolution, 2);
}

bool epi2D::isPointInPolygons(double& xloc, double& yloc) {
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

int epi2D::pnpoly(int& nvert, std::vector<double>& vertx, std::vector<double>& verty, double& testx, double& testy) {
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

void epi2D::updatePurseStringContacts() {
  // check if a vertex is wound adjacent but does not have a PS bond.
  // give that vertex a PS bond if it does not already have one.

  if (psContacts.size() < 2) {
    return;  // invalid size for sorting
  }

  /*cout << "psContacts : ";
  for (auto i : psContacts)
    cout << i << '\t';
  cout << '\n';*/

  int ci, cj, vi, c_first, c_second;
  double first = 1e10, second = 1e10;  // used to find two closest elements (minimum distances) in one traversal of psContacts
  int firstInd = INT_MAX, secondInd = INT_MAX, insert_index;
  int indexOfPsContacts_first = 0, indexOfPsContacts_second = 0, indexOfPsContacts_third = 0;
  std::vector<std::pair<double, int>> distanceToPsContacts;
  double l0_insert;
  bool isBetweenNeighbors = false;   // for debugging
  bool isTopologicalInsert = false;  // true when doing insertion based on topology, false when based on geometry

  for (auto gi : currentWoundIndices) {
    if (std::find(psContacts.begin(), psContacts.end(), gi) == psContacts.end()) {
      // could not find gi in psContacts, so we have a current wound adjacent vertex with no PS bond - figure out where to insert gi into psContacts
      cindices(ci, vi, gi);
      auto im1_ind = std::find(psContacts.begin(), psContacts.end(), im1[gi]);
      auto ip1_ind = std::find(psContacts.begin(), psContacts.end(), ip1[gi]);
      isBetweenNeighbors = false;
      isTopologicalInsert = false;

      if (im1_ind != psContacts.end() && ip1_ind != psContacts.end()) {
        // case 1 (easy): both of gi's same-cell neighbors are in psContacts, just insert gi there.
        isBetweenNeighbors = true;
        isTopologicalInsert = true;
        indexOfPsContacts_first = im1_ind - psContacts.begin();
        indexOfPsContacts_second = ip1_ind - psContacts.begin();
        assert(im1[gi] == psContacts[indexOfPsContacts_first]);
      } else if (im1_ind != psContacts.end() || ip1_ind != psContacts.end()) {
        // case 2 (harder): if one of gi's same-cell neighbors are in psContacts, check if it can be inserted next to that neighbor
        // if this works, set isTopologicalInsert to true
        int neighborIndex, neighbor;
        int leftNeighbor, rightNeighbor, leftIndex, rightIndex, c_left, c_right;
        if (im1_ind != psContacts.end())
          neighborIndex = im1_ind - psContacts.begin();
        else
          neighborIndex = ip1_ind - psContacts.begin();
        neighbor = psContacts[neighborIndex];
        // indices of left and right of neighbor
        leftIndex = (neighborIndex - 1 + psContacts.size()) % psContacts.size();
        rightIndex = (neighborIndex + 1 + psContacts.size()) % psContacts.size();
        leftNeighbor = psContacts[leftIndex];
        rightNeighbor = psContacts[rightIndex];
        cindices(c_left, vi, leftNeighbor);
        cindices(c_right, vi, rightNeighbor);
        if ((c_left != ci && c_right == ci) || (c_left == ci && c_right != ci)) {
          // case 1 : one is in the same-cell and one is in a diff-cell
          // if gi fits between neighborIndex and same-cell, insert it there.
          // otherwise, leave it to a distance calculation.
          if (c_left == ci) {
            if (isFitBetween(gi, neighbor, leftNeighbor, ci)) {
              indexOfPsContacts_first = neighborIndex;
              indexOfPsContacts_second = leftIndex;
              isTopologicalInsert = true;
            }
          } else {
            if (isFitBetween(gi, neighbor, rightNeighbor, ci)) {
              indexOfPsContacts_first = neighborIndex;
              indexOfPsContacts_second = rightIndex;
              isTopologicalInsert = true;
            }
          }
        } else if ((c_left == ci && c_right == ci)) {
          // case 2 : both are same-cell
          // if gi fits between neighbor and either left or right, insert it there.
          // otherwise, leave it to the distance calculation
          if (isFitBetween(gi, neighbor, leftNeighbor, ci)) {
            indexOfPsContacts_first = neighborIndex;
            indexOfPsContacts_second = leftIndex;
            isTopologicalInsert = true;
          } else if (isFitBetween(gi, neighbor, rightNeighbor, ci)) {
            indexOfPsContacts_first = neighborIndex;
            indexOfPsContacts_second = rightIndex;
            isTopologicalInsert = true;
          }
        }
        // case 3 : both are diff-cell, give up on topological and move to distance calculation
      }

      if (!isTopologicalInsert) {  // if neither of the topological cases are satisfied above, try a distance calculation to determine neighbors
        distanceToPsContacts.clear();
        // distanceToPsContacts = <distance, index> pair, used to sort distances
        for (int i = 0; i < psContacts.size(); i++) {
          cindices(cj, vi, psContacts[i]);
          distanceToPsContacts.push_back(std::pair<double, int>(vertDistSqNoPBC(gi, psContacts[i]), psContacts[i]));
        }
        // save 0,1,2 for potential insertion locations
        sort(distanceToPsContacts.begin(), distanceToPsContacts.end());
        indexOfPsContacts_first = std::find(psContacts.begin(), psContacts.end(), distanceToPsContacts[0].second) - psContacts.begin();
        indexOfPsContacts_second = std::find(psContacts.begin(), psContacts.end(), distanceToPsContacts[1].second) - psContacts.begin();
        indexOfPsContacts_third = std::find(psContacts.begin(), psContacts.end(), distanceToPsContacts[2].second) - psContacts.begin();

        int diffOfIndices, first_gi, second_gi, realDiffOfIndices;
        bool isIndexFirstAndSecondSameCellNeighbors, isGiNeighborOfFirst, isGiNeighborOfSecond, isGiNeighborOfFirstOrSecond, isFirstSecondInSameCell;

        // attempts to choose the right indices: between 0-1, 0-2, and 1-2 among the shortest distances in the sorted list distanceToPsContacts
        diffOfIndices = abs(indexOfPsContacts_first - indexOfPsContacts_second);
        int diffOfIndicesOneThree = abs(indexOfPsContacts_first - indexOfPsContacts_third);
        int diffOfIndicesTwoThree = abs(indexOfPsContacts_second - indexOfPsContacts_third);

        first_gi = psContacts[indexOfPsContacts_first];
        second_gi = psContacts[indexOfPsContacts_second];
        cindices(c_first, vi, first_gi);
        cindices(c_second, vi, second_gi);
        realDiffOfIndices = abs(first_gi - second_gi);
        isIndexFirstAndSecondSameCellNeighbors = (c_first == c_second && (realDiffOfIndices == 1 || realDiffOfIndices == nv[ci] - 1));

        // cout << "trying to insert " << gi << " between " << first_gi << '\t' << second_gi << ", simclock = " << simclock << '\n';

        // just making sure that the selected indices are next to each other in the psContacts list
        if ((diffOfIndices != 1 && diffOfIndices != psContacts.size() - 1) || isIndexFirstAndSecondSameCellNeighbors) {
          // diffOfIndices was not adjacent, or indices are adjacent (in psContacts) and neighbors (topologically)
          // cout << "trying to switch out first and second for other pairs\n";
          // cout << "diffOfIndicesOneThree, TwoThree = " << diffOfIndicesOneThree << '\t' << diffOfIndicesTwoThree << '\n';
          // switch out first and second for other pairs if those pairs are adjacent in psContacts
          // set new first and second indices
          if (diffOfIndicesOneThree == 1 || diffOfIndicesOneThree == psContacts.size() - 1) {
            // cout << "trying 1-3\n";
            indexOfPsContacts_second = indexOfPsContacts_third;
          } else if (diffOfIndicesTwoThree == 1 || diffOfIndicesTwoThree == psContacts.size() - 1) {
            // cout << "trying 2-3\n";
            indexOfPsContacts_first = indexOfPsContacts_second;
            indexOfPsContacts_second = indexOfPsContacts_third;
          } else {
            // cout << "none worked, refusing to insert gi. continue with next gi.\n";
            //  this happens when a lot of vertices are very close to gi, and it can't differentiate which is the correct one
            //  in this case, check if gi has a neighbor in psContacts.
            if (ip1_ind != psContacts.end()) {
              int neighbor_index = ip1_ind - psContacts.begin();
              int neighbor_gi = psContacts[neighbor_index];
              // cout << "gi of right neighbor = " << neighbor_gi << ", differences of right neighbor with its neighbors = " << neighbor_gi - psContacts[(neighbor_index + 1 + psContacts.size()) % psContacts.size()] << ", " << neighbor_gi - psContacts[(neighbor_index - 1 + psContacts.size()) % psContacts.size()] << '\n';
            }
            continue;
            // assert(false);
          }

          // reset variables to match new first and second indices
          diffOfIndices = abs(indexOfPsContacts_first - indexOfPsContacts_second);
          first_gi = psContacts[indexOfPsContacts_first];
          second_gi = psContacts[indexOfPsContacts_second];
          cindices(c_first, vi, first_gi);
          cindices(c_second, vi, second_gi);
          realDiffOfIndices = abs(first_gi - second_gi);
          isIndexFirstAndSecondSameCellNeighbors = (c_first == c_second && (realDiffOfIndices == 1 || realDiffOfIndices == nv[ci] - 1));

          // if first and second are not adjacent in psContacts, then inserting between them will be a large discontinuous break in the shape of the PS cable. so don't allow insertion between non-adjacent elements
          if (diffOfIndices != 1 && diffOfIndices != psContacts.size() - 1) {
            cout << "simclock = " << simclock << ", skipping insertion of " << gi << " between " << first_gi << ", and " << second_gi;
            cout << " because diff(indices) = " << diffOfIndices << ", psContacts.size - 1 = " << psContacts.size() - 1 << '\n';
            assert(!isBetweenNeighbors);
            continue;
          }
        }

        insert_index = 0;
        isGiNeighborOfFirst = (gi == ip1[first_gi] || gi == im1[first_gi]);
        isGiNeighborOfSecond = (gi == ip1[second_gi] || gi == im1[second_gi]);
        isGiNeighborOfFirstOrSecond = (isGiNeighborOfFirst || isGiNeighborOfSecond);
        isFirstSecondInSameCell = (c_first == c_second);

        // having made sure that first and second are adjacent in psContacts, consider whether they are the right indices to insert between.
        if (isFirstSecondInSameCell) {
          // first and second in same cell, make sure gi is inserted between them.
          if (isIndexFirstAndSecondSameCellNeighbors) {
            // if first and second end up being same cell neighbors, then we can't find an insert location and we should move on to the next iteration
            continue;
            // assert(!isIndexFirstAndSecondSameCellNeighbors);
          }
          // if we forbid this case, we simplify our insertion analysis
          // since first and second cannot be neighbors by this assumption, we only need to check if gi is 2 indices apart from first and second, in which case we are good for an insertion
          int diffGiFirst = abs(gi - first_gi);
          int diffGiSecond = abs(gi - second_gi);
          int sizeDiffFirst = nv[ci] - diffGiFirst;
          int sizeDiffSecond = nv[ci] - diffGiSecond;
          int numberOfSatisfiedConditions = ((diffGiFirst <= 2) + (diffGiSecond <= 2) + (sizeDiffFirst <= 2) + (sizeDiffSecond <= 2));
          // we are good for an insertion if any two of these conditions are <= 2.
          if (numberOfSatisfiedConditions >= 2) {
            // cout << "first and second are in the same cell, and we've detected that gi is within 2 indices of first and second in " << numberOfSatisfiedConditions << " instances\n";
          } else {
            cout << "first and second are in the same cell, and we haven't detected a good insertion point for gi. calling continue.\n";
            cout << "gi, first_gi, second_gi = " << gi << '\t' << first_gi << '\t' << second_gi << '\n';
            continue;
          }
        } else {
          // first and second are not in same cell, but are adjacent in psContacts. this should be ok
          // cout << "first and second are adjacent in psContact, but not in same cell. I think that's ok?\n";
        }
      }

      // insert gi into psContacts in the middle of these adjacent elements
      if ((indexOfPsContacts_first == 0 && indexOfPsContacts_second == psContacts.size() - 1) || (indexOfPsContacts_second == 0 && indexOfPsContacts_first == psContacts.size() - 1)) {
        // edge case - insert gi at the beginning of the entire vector
        // gi, first, ... , second
        insert_index = 0;
      } else if (indexOfPsContacts_first < indexOfPsContacts_second) {
        // ... first, gi, second, ...
        insert_index = indexOfPsContacts_first + 1;
      } else {
        // ... second, gi, first
        insert_index = indexOfPsContacts_second + 1;
      }
      // cout << "inserting " << gi << " at index " << insert_index << " between " << psContacts[(insert_index - 1 + psContacts.size()) % psContacts.size()] << " and " << psContacts[insert_index] << '\n';

      psContacts.insert(psContacts.begin() + insert_index, gi);

      // next and previous represent the index of the next and previous neighbor of the newly inserted index
      int next = (insert_index + 1 + psContacts.size()) % psContacts.size();
      int prev = (insert_index - 1 + psContacts.size()) % psContacts.size();

      // x_ps, v_ps, F_ps, l0, isSpringBroken
      x_ps.insert(x_ps.begin() + NDIM * insert_index, x.begin() + NDIM * gi, x.begin() + NDIM * gi + 2);
      v_ps.insert(v_ps.begin() + NDIM * insert_index, v.begin() + NDIM * gi, v.begin() + NDIM * gi + 2);
      F_ps.insert(F_ps.begin() + NDIM * insert_index, F.begin() + NDIM * gi, F.begin() + NDIM * gi + 2);
      isSpringBroken.insert(isSpringBroken.begin() + insert_index, false);

      // new l0 value is dist(x_ps(i + 1 + size % size), x[gi])
      // compute l0 last because it depends on new additions to x_ps
      l0_insert = sqrt(pow(x_ps[NDIM * prev] - x[NDIM * gi], 2) + pow(x_ps[NDIM * prev + 1] - x[NDIM * gi + 1], 2));
      l0_ps[insert_index] = sqrt(pow(x_ps[NDIM * next] - x[NDIM * gi], 2) + pow(x_ps[NDIM * next + 1] - x[NDIM * gi + 1], 2));
      l0_ps.insert(l0_ps.begin() + insert_index, l0_insert);
    }
  }
}

void epi2D::purseStringContraction() {
  purseStringTension = 0.0;
  purseStringTransmittedTension = 0.0;
  updatePurseStringContacts();
  integratePurseString();  // evaluate forces on and due to purse-string, and integrate its position
  for (int psi = 0; psi < psContacts.size(); psi++) {
    if (l0_ps[psi] <= 0.01 * r[0]) {
      l0_ps[psi] = 0.01 * r[0];
      // cout << "l0_ps belonging to " << psContacts[psi] << " is less than the threshold, setting it to minimum!\n";
    } else {
      // l0_ps[psi] *= exp(-strainRate_ps * dt); // constant strain rate
      l0_ps[psi] -= strainRate_ps * dt;  // constantly increasing tension until length < r[0], the smallest physical lengthscale
    }
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

void epi2D::evaluatePurseStringForces() {
  // using psContacts and the preferred length interaction, calculate forces on purse-string virtual particles
  // deltaSq is the breaking point for the virtual particle-wound particle interaction
  // k_ps is the wound-purseString spring force constant
  // B is the damping constant
  int gi, ipi, imi;
  double xp, yp, xw, yw;  // purse-string and wound coordinates
  bool isCutoff;          // cutoff distance
  double yieldLengthSquared = deltaSq * pow(r[0] * 2, 2);
  double M = sqrt(yieldLengthSquared) / 2.0;  // M is max distance at which the spring force saturates, sqrt(M) < yieldLengthSquared
  double dx, dy, fx, fy, l0i, l0im1;

  double fli, flim1, cx, cy, xi, yi, gip1, xip1, yip1;
  double rho0 = sqrt(a0.at(0));
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
    isCutoff = ((dx * dx + dy * dy) >= yieldLengthSquared);
    isSpringBroken[psi] = isCutoff;  // record broken springs for printouts
    if (!isCutoff) {
      // spring force saturates at a max distance of M, M < cutoff distance
      if (dx > M)
        dx = M;
      if (dy > M)
        dy = M;
      fx = k_ps * dx;
      fy = k_ps * dy;

    } else {  // fx, dx could be nan. in which case this will protect from nans rolling over into my actual vertices
      fx = 0;
      fy = 0;
      dx = 0;
      dy = 0;
    }
    if (std::isnan(fx)) {
      cout << "from purse-string interaction, fx from psi = " << psi << " is NaN! this affects vertex " << gi << '\n';
    }
    F[gi * NDIM] += fx;
    F[gi * NDIM + 1] += fy;
    F_ps[psi * NDIM] -= fx;
    F_ps[psi * NDIM + 1] -= fy;
    // cout << "force on pursestring due to virtual-real bonds = " << -fx << '\t' << -fy << '\n';
    purseStringTransmittedTension += sqrt(pow(fx, 2) + pow(fy, 2));

    if (std::isnan(F_ps[NDIM * psi])) {
      cout << "force components responsible for F_ps nan: " << '\n';
      cout << "fx = " << fx << ", fy = " << fy << '\n';
      cout << "isCutoff = " << isCutoff << ", isSpringBroken = " << isSpringBroken[psi] << ", dx dy = " << dx << '\t' << dy << '\n';
    }

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
    // cout << "kl inside ps = " << kl << '\n';

    // add to forces
    fx = (fli * dli * lix / li) - (flim1 * dlim1 * lim1x / lim1);
    fy = (fli * dli * liy / li) - (flim1 * dlim1 * lim1y / lim1);

    F_ps[NDIM * psi] += fx;
    F_ps[NDIM * psi + 1] += fy;
    purseStringTension += sqrt(pow(fx, 2) + pow(fy, 2));
    // cout << "F_ps_x due to segment length = " << fx << '\t' << fy << '\n';
    if (l0_ps.size() != psContacts.size()) {
      cout << "contradiction: l0_ps.size != psContacts.size() : " << l0_ps.size() << '\t' << psContacts.size() << '\n';
      assert(l0_ps.size() == psContacts.size());
    }
    if (std::isnan(F_ps[NDIM * psi])) {
      cout << "force components responsible for F_ps nan: " << '\n';
      cout << "fx = " << fx << ", fy = " << fy << '\n';
      cout << "fli, dli, lix, liy, li = " << fli << '\t' << dli << '\t' << lix << '\t' << liy << '\t' << li << '\n';
      cout << "flim1, dlim1, lim1x, lim1y, lim1 = " << flim1 << '\t' << dlim1 << '\t' << lim1x << '\t' << lim1y << '\t' << lim1 << '\n';
      cout << "l0i, l0im1 = " << l0i << '\t' << l0im1 << '\n';
      cout << "l0_ps.size() = " << l0_ps.size() << '\n';
      cout << "l0_ps = ";
      for (auto i : l0_ps) {
        cout << i << '\t';
      }
      cout << "\npsi, imi, psContacts.size = " << psi << '\t' << imi << '\t' << psContacts.size() << '\n';
    }

    // choosing not to update potential energy of the purse-string itself
    // U += 0.5 * kl * (dli * dli);
  }
}

void epi2D::integratePurseString() {
  // velocity verlet force update with damping

  // debugging test case : ./main/epi2D/laserAblation.o 36 36 3 1.10 0.94 0.85 1.0 1.0 0.05 0.005  2.0  4.0  4.0 1.0  0.0  1.0 0.5  0  0   0 1  1000  test

  // first step: delete virtual vertices if the virtual-real bond has yielded.
  for (int i = 0; i < psContacts.size(); i++) {
    if (isSpringBroken[i]) {
      cout << "spring broken on particle " << psContacts[i] << " with psContacts index " << i << '\n';
      int prev = (i - 1 + psContacts.size()) % psContacts.size();
      int next = (i + 1 + psContacts.size()) % psContacts.size();
      // cout << "marking a spring on gi = psContact[i] = " << psContacts[i] << " for deletion!\n";
      //  mark all doubles associated with yielded virtual vertex for deletion
      x_ps[NDIM * i] = NAN;
      x_ps[NDIM * i + 1] = NAN;
      v_ps[NDIM * i] = NAN;
      v_ps[NDIM * i + 1] = NAN;
      F_ps[NDIM * i] = NAN;
      F_ps[NDIM * i + 1] = NAN;
      // l0_ps[(i - 1 + psContacts.size()) % psContacts.size()] += l0_ps[i];
      l0_ps[prev] = vertDistNoPBC(psContacts[prev], psContacts[next]);  // adjust previous indexed l0 to match new configuration with current l0 deleted
    }
  }
  // set l0_ps = NAN and psContacts = 9999 (deletion criteria) after the rest to not interfere with l0_ps adjustment
  for (int i = 0; i < psContacts.size(); i++) {
    if (isSpringBroken[i]) {
      psContacts[i] = 9999;
      l0_ps[i] = NAN;
    }
  }

  // delete all NANs, the mark for deletion
  x_ps.erase(remove_if(x_ps.begin(), x_ps.end(), [](const double& value) { return std::isnan(value); }), x_ps.end());
  v_ps.erase(remove_if(v_ps.begin(), v_ps.end(), [](const double& value) { return std::isnan(value); }), v_ps.end());
  F_ps.erase(remove_if(F_ps.begin(), F_ps.end(), [](const double& value) { return std::isnan(value); }), F_ps.end());
  l0_ps.erase(remove_if(l0_ps.begin(), l0_ps.end(), [](const double& value) { return std::isnan(value); }), l0_ps.end());
  psContacts.erase(remove(psContacts.begin(), psContacts.end(), 9999), psContacts.end());

  // delete all bools associated with yielded virtual vertices for deletion
  isSpringBroken.erase(remove(isSpringBroken.begin(), isSpringBroken.end(), true), isSpringBroken.end());

  // VV position update
  int virtualDOF = psContacts.size() * NDIM;
  for (int i = 0; i < virtualDOF; i++) {
    x_ps[i] += dt * v_ps[i] + 0.5 * dt * dt * F_ps[i];
    if (std::isnan(x_ps[i])) {
      cout << "x_ps[i] = nan\n";
      cout << "v_ps, F_ps = " << v_ps[i] << '\t' << F_ps[i] << '\n';
    }
    if (std::isnan(F_ps[i])) {
      cout << "VV pos update: F_ps[" << i << "] = nan\n";
    }

    // recenter in box
    if (x_ps[i] > L[i % NDIM] && pbc[i % NDIM])
      x_ps[i] -= L[i % NDIM];
    else if (x_ps[i] < 0 && pbc[i % NDIM])
      x_ps[i] += L[i % NDIM];
  }

  // force update
  std::vector<double> F_ps_old = F_ps;
  evaluatePurseStringForces();

  // VV VELOCITY UPDATE #2
  for (int i = 0; i < virtualDOF; i++) {
    F_ps[i] -= (B * v_ps[i] + B * F_ps_old[i] * dt / 2);
    F_ps[i] /= (1 + B * dt / 2);
    v_ps[i] += 0.5 * (F_ps[i] + F_ps_old[i]) * dt;
    if (std::isnan(F_ps[i])) {
      cout << "VV vel update: F_ps[" << i << "] = nan\n";
    }
  }
}

void epi2D::sortPurseStringVariables() {
  // call this periodically because inserting vertices has some unwanted edge cases
  // list of data to handle: x_ps, v_ps, F_ps, l0_ps, isSpringBroken, psContacts
  // psContacts is an out of order vector, such as
  //    psContacts = 117 220 221 222 224 223 390 391 225 226 207 392
  //  it should be sorted to 117 220 221 222 223 224 225 226 207 390 391 392 ...
  //    or something equivalent to that.
  // also, psContacts[i] corresponds to l0_ps[i] and x_ps[NDIM*i] and x_ps[NDIM*i + 1].
  //  when sorting psContacts, make sure to also sort l0_ps, x_ps, v_ps, F_ps, isSpringBroken accordingly
  //
}

bool epi2D::isFitBetween(int gi, int gl, int gr, int ci) {
  // check if index gi fits between gl and gr. hardcoded the range = 4. Can change this, but definitely do not let it approach nv/2.
  // NOTE: gi gl and gr MUST be in the same cell for this to work.
  int numVerts = nv[ci];
  int range = numVerts / 3.0, rightNeighbor, leftNeighbor;
  int currentIndex = gi;
  bool hasLeft = false, hasRight = false;
  // want to check validity of gl gi gr
  for (int i = 0; i < range; i++) {
    currentIndex = ip1[currentIndex];
    if (gr == currentIndex || gl == currentIndex) {
      hasRight = true;
    }
  }
  currentIndex = gi;
  for (int i = 0; i < range; i++) {
    currentIndex = im1[currentIndex];
    if (gr == currentIndex || gl == currentIndex) {
      hasLeft = true;
    }
  }
  return (hasLeft && hasRight);
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