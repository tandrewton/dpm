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

// initialize monodisperse cell system, single calA0
void epi2D::monodisperse2D(double calA0, int n) {
  // local variables
  double calA0tmp, calAntmp, rtmp, areaSum;
  int vim1, vip1, gi, ci, vi, nlarge, smallN, largeN, NVSMALL;

  // print to console
  cout << "** initializing monodisperse DPM particles in 2D ..." << endl;

  // total number of vertices
  NVTOT = n * NCELLS;
  vertDOF = NDIM * NVTOT;

  // szList and nv (keep track of global vertex indices)
  nv.resize(NCELLS);
  szList.resize(NCELLS);

  nv.at(0) = n;
  for (ci = 1; ci < NCELLS; ci++) {
    nv.at(ci) = n;
    szList.at(ci) = szList.at(ci - 1) + nv.at(ci - 1);
  }

  // initialize vertex shape parameters
  initializeVertexShapeParameters(calA0, n);

  // initialize vertex indexing
  initializeVertexIndexing2D();
}

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

double epi2D::distanceLineAndPoint(double x1, double y1, double x2, double y2, double x0, double y0) {
  // get the distance from a line segment going through (x1,y1), (x2,y2) and a
  // point located at (x0,y0)
  double l2 = pow(x2 - x1, 2) + pow(y2 - y1, 2);  // |(pt2 - pt1)|^2
  if (l2 == 0.0)                                  // pt2 == pt1 case
    return sqrt(pow(x0 - x1, 2) + pow(y0 - y1, 2));

  double dot = (x0 - x1) * (x2 - x1) +
               (y0 - y1) * (y2 - y1);  // (pt0 - pt1) dot (pt2 - pt1)
  const double t = max(0.0, min(1.0, dot / l2));
  const double projectionx = x1 + t * (x2 - x1);
  const double projectiony = y1 + t * (y2 - y1);
  const double distance =
      sqrt(pow(x0 - projectionx, 2) + pow(y0 - projectiony, 2));
  return distance;
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

std::vector<int> epi2D::regridSegment(int wVertIndex, double vrad) {  // return value is the indices of the deleted vertices for proper gi counting later
  if (vnn_label[wVertIndex] != 0)
    return {};

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
  double L0tmp = getPreferredPerimeter(ci), L0new;
  double LCurrent = perimeter(ci);
  double l0new;

  // compute cell center of mass - don't necessarily need it, but it makes it easier to handle pbc if they are present.
  com2D(ci, cx, cy);
  // need ip1 to get rip1x, rip1y. need ? to get rix, riy. then rip1x, rix give lix. then li = lix^2 + liy^2.

  double c_regrid = 0;
  double deltaESegment = 0;
  for (int gi = gindex(ci, 0); gi < gindex(ci, 0) + nv[ci]; gi++) {
    // get coordinates relative to center of mass
    rix = x[NDIM * gi];
    riy = x[NDIM * gi + 1];

    // get next adjacent vertices
    rip1x = x[NDIM * ip1[gi]];
    rip1y = x[NDIM * ip1[gi] + 1];

    lix = rip1x - rix;
    liy = rip1y - riy;

    li = sqrt(lix * lix + liy * liy);

    deltaESegment += pow(li / l0[gi] - 1, 2);
  }

  std::vector<int> deleteList;
  std::vector<double> result;
  std::vector<int> unorderedWoundList;

  // use unsortedWoundIndices to find ci's wound-adjacent vertices
  int cj, vj, gi, gj;
  for (int j = 0; j < unsortedWoundIndices.size(); j++) {
    cindices(cj, vj, unsortedWoundIndices[j]);
    if (cj == ci) {  // collect all indices for cell ci that are wound-adjacent
      unorderedWoundList.push_back(unsortedWoundIndices[j]);
    }
  }

  std::vector<int> orderedWVerts = unorderedWoundList;

  cout << "NOT rotated orderedWVerts looks like (before deletion): \n";
  for (auto i : orderedWVerts) {
    cout << "gi = " << i << "\t, ip1[gi] = " << ip1[i] << "\t, im1[gi] = " << im1[i] << '\n';
  }

  for (int i = 0; i < orderedWVerts.size(); i++) {
    // check if ip1[orderedWVerts[j]] == orderedWVerts[j+1] for all j
    // if not, then do a cyclic permutation

    bool isOrderedRight = true;
    bool isOrderedLeft = true;
    for (int j = 0; j < orderedWVerts.size() - 1; j++) {
      isOrderedRight *= (ip1[orderedWVerts[j]] == orderedWVerts[j + 1]);
      isOrderedLeft *= (im1[orderedWVerts[(j + orderedWVerts.size()) % orderedWVerts.size()]] == orderedWVerts[(j - 1 + orderedWVerts.size()) % orderedWVerts.size()]);
    }
    if (isOrderedRight || isOrderedLeft) {
      break;
    } else {
      std::rotate(orderedWVerts.begin(), orderedWVerts.begin() + 1, orderedWVerts.end());
    }
  }

  cout << "rotated orderedWVerts looks like (before deletion): \n";
  for (auto i : orderedWVerts) {
    cout << "gi = " << i << "\t, ip1[gi] = " << ip1[i] << "\t, im1[gi] = " << im1[i] << '\n';
  }

  // orderedWVerts.push_back(ip1[gj]); // try extending the vertices to interpolate by one on each end
  // orderedWVerts.insert(orderedWVerts.begin(), im1[gi]);

  for (int i = 0; i < orderedWVerts.size() - 1; i++) {
    arc_length += vertDistNoPBC(orderedWVerts[i], orderedWVerts[i + 1]);
  }
  // nv_new = total length of segments divided by diameter of particle
  int nv_new = ceil(arc_length / vdiam);
  /*cout << "perimeter = " << arc_length
       << ", arc_length = " << arc_length / nv_new
       << ", actual arc_length should be = " << vdiam << '\n';*/
  arc_length /= nv_new;

  int numVertsToDelete = orderedWVerts.size() - nv_new;
  if (numVertsToDelete <= 0)
    return {};

  for (int i = 0; i < orderedWVerts.size() - 1; i++) {
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

  if (numVertsToDelete > 1) {
    cout << "deleting " << numVertsToDelete << " vertices in regridSegment!\n";
    // printConfiguration2D();
    //   printBoundaries();
  }

  // orderedWVerts.size() - 1 : don't delete the last vertex, it's not accounted for
  for (int i = 1; i < orderedWVerts.size() - 1; i++) {
    gi = orderedWVerts[i];
    if (i < result.size() / 2) {
      x[gi * NDIM] = result[2 * i];
      x[gi * NDIM + 1] = result[2 * i + 1];
      // cout << "modifying position of vertex " << gi << '\n';
    } else
      deleteList.push_back(gi);
  }

  if (!deleteList.empty()) {
    cout << "before delete\t" << nv[ci] << '\n';

    cout << "deleteList has " << deleteList.size()
         << " vertices in it, and we are asked to delete " << numVertsToDelete
         << " vertices\n";
    for (auto i : deleteList)
      cout << "vertex to delete: " << i << '\t' << x[NDIM * i] << '\t' << x[NDIM * i + 1] << '\n';
    cout << "wVertIndex = " << wVertIndex << '\n';
    for (auto i : orderedWVerts)
      cout << "orderedWVerts : " << i << '\t' << x[NDIM * i] << '\t' << x[NDIM * i + 1] << '\n';

    deleteVertex(deleteList);

    for (auto i : deleteList) {
      orderedWVerts.erase(std::remove(orderedWVerts.begin(), orderedWVerts.end(), i), orderedWVerts.end());
      for (int j = 0; j < orderedWVerts.size(); j++) {
        if (orderedWVerts[j] > i)
          orderedWVerts[j]--;
      }
    }

    // after regridding and deleting, have all wound vertices set their preferred lengths to their current length. gives rigidity.

    for (auto i : szList) {
      cout << "szList = " << i << '\n';
    }

    for (int i = 0; i < orderedWVerts.size(); i++) {
      gi = orderedWVerts[i];
      if (i == 0)
        l0new = vertDistNoPBC(gi, ip1[gi]);
      cout << "gi, l0(before), l0(after) = " << gi << '\t' << l0[gi] << "\t" << vertDistNoPBC(gi, ip1[gi]) << '\n';
      l0[gi] = vertDistNoPBC(gi, ip1[gi]);
      /*cout << "gi, l0(before), l0(after) = " << gi << '\t' << l0[gi] << "\t" << l0new << '\n';
      l0[gi] = l0new;*/
      cout << "orderedWVerts after regrid: " << gi << ", ip1[gi] = " << ip1[gi] << '\t' << x[NDIM * gi] << '\t' << x[NDIM * gi + 1] << '\n';
    }
  }
  return deleteList;
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

void epi2D::epi_shapeForces2D() {
  // difference from dpm::shapeForces2D is that this has an actual perimeter (L) force
  //   as opposed to a segment length force
  // local variables
  int ci, gi, vi, nvtmp;
  double fa, fL, fli, flim1, fb, cx, cy, xi, yi;
  double rho0, l0im1, l0i, L0tmp, Ltmp, a0tmp, atmp;
  double dx, dy, da, dL, dli, dlim1, dtim1, dti, dtip1;
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
        // shape information
        nvtmp = nv[ci];
        a0tmp = a0[ci];
        L0tmp = getPreferredPerimeter(ci);

        // preferred segment length of last segment
        l0im1 = l0[im1[gi]];

        // compute area deviation
        atmp = area(ci);
        da = (atmp / a0tmp) - 1.0;

        // compute perimeter deviation
        Ltmp = perimeter(ci);
        dL = (Ltmp / L0tmp) - 1.0;

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

    // perimeter force
    fL = kL * (rho0 / L0tmp);

    // add to forces
    forceX = (fli * dli * lix / li) - (flim1 * dlim1 * lim1x / lim1);
    forceY = (fli * dli * liy / li) - (flim1 * dlim1 * lim1y / lim1);
    F[NDIM * gi] += forceX;
    F[NDIM * gi + 1] += forceY;

    // fL, dL - Andrew 10/28/21 - total perimeter force for cts. energy while regridding
    //   (distinct from original segment force)
    forceX = fL * dL * (lix - lim1x) / L0tmp;
    forceY = fL * dL * (liy - lim1y) / L0tmp;
    F[NDIM * gi] += forceX;
    F[NDIM * gi + 1] += forceY;

    // note - Andrew here, confirmed that the shape stress matrix is diagonal as written
    fieldShapeStress[gi][0] += unwrappedX * forceX;
    fieldShapeStress[gi][1] += unwrappedY * forceY;
    fieldShapeStress[gi][2] += unwrappedX * forceY;

    // update potential energy
    U += 0.5 * kl * (dli * dli);
    cellU[ci] += 0.5 * kl * (dli * dli);

    // -- Bending force
    if (kb > 0) {
      // get ip2 for third angle
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

void epi2D::vertexAttractiveForces2D_2() {
  // altered from dpm attractive force code, because it works with larger l2
  // values. (warning: probably won't work with bending.) local variables
  int ci, cj, gi, gj, vi, vj, bi, bj, pi, pj;
  double sij, rij, dx, dy, rho0;
  double ftmp, fx, fy;

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
              } /*else if (std::find(listTurnOffRepulsion.begin(), listTurnOffRepulsion.end(), gi) != listTurnOffRepulsion.end()) {
                // if purseStringContraction determines that a vertex should not repel others intercellularly, turn it off here.
                ftmp = 0;
              } */
              else if (rij > cutij) {
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
              stress[0] += dx * fx;
              stress[1] += dy * fy;
              stress[2] += 0.5 * (dx * fy + dy * fx);

              fieldStress[gi][0] += dx * fx;
              fieldStress[gi][1] += dy * fy;
              fieldStress[gi][2] += 0.5 * (dx * fy + dy * fx);

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
                } /*else if (std::find(listTurnOffRepulsion.begin(), listTurnOffRepulsion.end(), gi) != listTurnOffRepulsion.end()) {
                  // if purseStringContraction determines that a vertex should not repel others intercellularly, turn it off here.
                  ftmp = 0;
                } */
                else if (rij > cutij) {
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
                stress[0] += dx * fx;
                stress[1] += dy * fy;
                stress[2] += 0.5 * (dx * fy + dy * fx);

                fieldStress[gi][0] += dx * fx;
                fieldStress[gi][1] += dy * fy;
                fieldStress[gi][2] += 0.5 * (dx * fy + dy * fx);

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
  epi_shapeForces2D();
  vertexAttractiveForces2D_2();
}

void epi2D::activeAttractiveForceUpdate() {
  // compute forces for shape, attractive, and active contributions

  // reset forces, then get shape and attractive forces.
  attractiveForceUpdate_2();

  // compute active forces
  int gi = 0, ci = 0, vi = 0;
  double dpsi = 0.0, psitmp = 0.0;
  double xi, yi, nvtmp, dx, dy, cx, cy, r1, r2, grv, v0tmp, vmin, Ds, rnorm, ux,
      uy, rix, riy;

  vmin = 1e-1 * v0;  // min velocity
  Ds = 0.2;          // active velocity spread parameter

  directorDiffusion();

  for (int gi = 0; gi < NVTOT; gi++) {
    // get coordinates relative to center of mass
    rix = x[NDIM * gi] - cx;
    riy = x[NDIM * gi + 1] - cy;
    if (pbc[0] == false && round(rix / L[0]) != 0) {
      cout << "error! rix is outside the box even though pbc are turned off!";
    }
    if (pbc[1] == false && round(riy / L[1]) != 0) {
      cout << "error! riy is outside the box even though pbc are turned off!";
    }
    if (pbc[0])
      rix -= L[0] * round(rix / L[0]);
    if (pbc[1])
      riy -= L[1] * round(riy / L[1]);

    // get angular distance from psi
    psitmp = atan2(riy, rix);
    dpsi = psitmp - psi[ci - 1];
    dpsi -= 2 * PI * round(dpsi / (2 * PI));

    // get velocity scale
    v0tmp = vmin + (v0 - vmin) * exp(-pow(dpsi, 2.0) / (2.0 * Ds * Ds));

    v0tmp *= activePropulsionFactor[ci];

    // get unit vectors
    rnorm = sqrt(rix * rix + riy * riy);
    ux = rix / rnorm;
    uy = riy / rnorm;

    // add to forces
    F[NDIM * gi] += v0tmp * ux;
    F[NDIM * gi + 1] += v0tmp * uy;
  }
  if (boolCIL == true)
    deflectOverlappingDirectors();
}

void epi2D::substrateadhesionAttractiveForceUpdate() {
  // compute forces for shape, attractive, and substrate adhesion contributions
  int gi = 0, argmin, flagcount = 0;
  int numVerticesAttracted = 1;
  double k = 2;
  double refreshInterval = 1;

  // reset forces, then get shape and attractive forces.
  attractiveForceUpdate_2();
  directorDiffusion();
  updateSubstrateSprings(refreshInterval);

  for (int ci = 0; ci < NCELLS; ci++) {
    std::vector<double> distSqVertToFlag(nv[0], 0.0);
    // check for protrusions (unimpeded flag toss, and flag location inside box)
    if (flag[ci]) {
      // find nearest vertices
      for (int vi = 0; vi < nv[ci]; vi++) {
        gi = gindex(ci, vi);
        distSqVertToFlag[vi] = pow(x[gi * NDIM] - flagPos[ci][0], 2) +
                               pow(x[gi * NDIM + 1] - flagPos[ci][1], 2);
      }
      argmin = std::distance(
          distSqVertToFlag.begin(),
          std::min_element(distSqVertToFlag.begin(), distSqVertToFlag.end()));
      gi = gindex(ci, argmin);

      // evaluate force for spring-vertex interaction between nearest vertex and
      // flag position
      F[gi * NDIM] += -k * (x[gi * NDIM] - flagPos[ci][0]);
      F[gi * NDIM + 1] += -k * (x[gi * NDIM + 1] - flagPos[ci][1]);
      flagcount++;
    }
  }
  // if (boolCIL == true)
  //   deflectOverlappingDirectors();
  /*if (flagcount >= 1) {
    cout << "simclock = " << simclock << '\t' << ", number of flags = " <<
  flagcount << '\n'; for (int ci = 0; ci < NCELLS; ci++) { if (flag[ci]) cout <<
  "cell " << ci << " has a flag!\n";
    }
  }*/
}

/******************************

        EPI

                P R O T O C O L S

*******************************/

void epi2D::vertexCompress2Target2D(dpmMemFn forceCall, double Ftol, double dt0, double phi0Target, double dphi0) {
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
    if (it % 50 == 0) {
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

void epi2D::updateSubstrateSprings(double refreshInterval) {
  // check to see if enough time has passed for us to update springs again
  if (simclock - previousUpdateSimclock > refreshInterval) {
    double cx, cy, gj;
    double flagDistance = 1.5 * sqrt(a0[0] / PI);
    // issue: if cell becomes too elongated, then flagDistance isn't large
    // enough to anchor outside the cell idea: since updateSubstrate is only
    // called in the force update, have the force update pass a vector of
    // lengths instead of using flagDistance here
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
    for (int ci = 0; ci < NCELLS; ci++) {
      cancelFlagToss = false;
      if (!flag[ci]) {
        // probably not a good idea to vary flag distance at this stage
        // double randFlagDistance = flagDistance * (drand48() / 2 + 1);
        double randFlagDistance = flagDistance;
        // pick a direction, throw a flag in that direction
        flagPos[ci][0] = center[ci][0] + randFlagDistance * cos(psi[ci]);
        flagPos[ci][1] = center[ci][1] + randFlagDistance * sin(psi[ci]);
        // loop over cells near enough to block flag
        for (int cj = 0; cj < NCELLS; cj++) {
          if (cancelFlagToss == true)
            continue;
          if (ci != cj) {
            // check if centers are near enough to possibly block flag
            if (pow(center[ci][0] - center[cj][0], 2) +
                    pow(center[ci][1] - center[cj][1], 2) <
                3 * pow(randFlagDistance, 2)) {
              // check if vertices actually block flag
              for (int vj = 0; vj < nv[cj]; vj++) {
                gj = gindex(cj, vj);
                if (distanceLineAndPoint(center[ci][0], center[ci][1],
                                         flagPos[ci][0], flagPos[ci][1],
                                         x[gj * NDIM],
                                         x[gj * NDIM + 1]) < 2 * r[gj]) {
                  // yes, flag has been blocked. Move on to next cell's flag
                  // toss attempt.
                  cancelFlagToss = true;
                  // inhibited, so director goes in opposite direction.
                  // psi[cj] += PI;
                  // psi[cj] -= 2 * PI * round(psi[cj] / (2 * PI));
                  // maybe this should be a break instead of a continue
                  break;
                }
              }
            }
          }
        }
        // if flag throw is successful, set flag[ci] = true and record flagPos
        if (cancelFlagToss == false) {
          flag[ci] = true;
        }
      } else if (flag[ci]) {
        // dissociation rate determines if existing flag is destroyed this step
        if (drand48() < 0.1) {
          flag[ci] = false;
        }
      } else
        cerr << "Error: no flag bool value found\n";
    }
    //
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
          stressout << setw(wnum) << -(stress[0] + shapeStressXX);
          stressout << setw(wnum) << -(stress[1] + shapeStressYY);
          stressout << setw(wnum) << -(stress[2] + shapeStressXY);
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

void epi2D::dampedNP0(dpmMemFn forceCall, double B, double dt0, double duration, double printInterval, bool wallsOn) {
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
      // purseStringContraction(0.01) is fast (~100 tau), 0.001 is slow (~1000 tau).
      purseStringContraction(0.01);  // must be called after forceCall, beacuse forceCall reestablishes the boundaries() properly.
      // evaluateGhostDPForces(woundVertexList, 0.001);
      if (regridChecker) {
        printConfiguration2D();
        printBoundaries();
        regridChecker = false;
      }
    }

    // VV VELOCITY UPDATE #2
    for (i = 0; i < vertDOF; i++) {
      F[i] -= (B * v[i] + B * F_old[i] * dt / 2);
      F[i] /= (1 + B * dt / 2);
      v[i] += 0.5 * (F[i] + F_old[i]) * dt;
    }

    if (wallsOn == true) {
      // VV position update (walls)
      L[0] += dt * VL[0] + 0.5 * dt * dt * FT;
      L[1] += dt * VL[1] + 0.5 * dt * dt * FR;

      // VV force update (walls)
      oldFT = FT;
      oldFB = FB;
      oldFL = FL;
      oldFR = FR;
      wallForces(true, false, false, true, FT, FB, FL, FR);
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
          stressout << setw(wnum) << -(stress[0] + shapeStressXX);
          stressout << setw(wnum) << -(stress[1] + shapeStressYY);
          stressout << setw(wnum) << -(stress[2] + shapeStressXY);
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

void epi2D::wallForces(bool top, bool bottom, bool left, bool right, double& forceTop, double& forceBottom, double& forceLeft, double& forceRight) {
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

  // need to delete an index from anything that depends on NVTOT, nv, szList?
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
    unsortedWoundIndices.erase(std::remove(unsortedWoundIndices.begin(), unsortedWoundIndices.end(), i), unsortedWoundIndices.end());

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

  if (simclock < 104.994 && simclock > 104.993) {
    cout << "debug report: simclock = " << simclock << '\n';
    for (int gi = 0; gi < NVTOT; gi++) {
      if (x[gi * NDIM] < 4.5 && x[gi * NDIM] > 4.3) {
        if (x[gi * NDIM + 1] < 3.4 && x[gi * NDIM + 1] > 3.2) {
          cout << "gi = " << gi << ", position of gi = " << x[gi * NDIM] << '\t' << x[gi * NDIM + 1] << "\n";
          for (auto i : vnn[gi]) {
            cout << "neighbor of gi = " << vnn[gi][i] << "\t, label[neighbor] = " << vnn_label[i] << '\n';
            cout << "position = " << x[i * NDIM] << '\t' << x[i * NDIM + 1] << '\n';
          }
        }
      }
    }
  }

  return voidFacingVertexIndices;
}

void epi2D::NewmanZiff(std::vector<int>& ptr, int empty, int& mode, int& big) {
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
//  also saves unsortedWoundIndices
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
  NewmanZiff(ptr, empty, mode, big);

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
    int nthClusterSize = clusterSize[nthLargestCluster - 1];
    for (auto i : clusterRoot) {
      if (-ptr[i] == nthClusterSize)
        mode = i;
    }
    cout << "new mode = " << mode << '\n';
    big = nthClusterSize;
    cout << "new size (big) = " << big << '\n';
  }

  // in no particular order, print out the wrapped locations of all main
  //  cluster vertices
  for (int gi = 0; gi < NVTOT; gi++) {
    // if (findRoot(gi, ptr) == mode) {
    // if (vnn_label[gi] == 0 || vnn_label[gi] == 2 || vnn_label[gi] == 1 || vnn_label[gi] == -1 || vnn_label[gi] == -2) {
    bool isDisabled = std::find(listTurnOffRepulsion.begin(), listTurnOffRepulsion.end(), gi) !=
                      listTurnOffRepulsion.end();

    if (vnn_label[gi] == 0 || vnn_label[gi] == 2 || vnn_label[gi] == 1 || vnn_label[gi] == -1 || vnn_label[gi] == -2 || vnn_label[gi] == -3) {
      edgeout << x[gi * NDIM] << '\t' << x[gi * NDIM + 1] << '\t' << vnn_label[gi] + (findRoot(gi, ptr) == mode) * 10 + isDisabled * 100 << '\n';
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
  if (checkWoundClosed(unwrapped_x_gi)) {
    for (int i = 0; i < unwrapped_x.size(); i += 2)
      bout << unwrapped_x[i] << '\t' << unwrapped_x[i + 1] << '\n';
  } else {  // unwrapped_x_gi is not closed, hence does not represent what I'm using for the void segmentation
    for (int i = 0; i < unsortedWoundIndices.size() / 2; i += 2) {
      bout << x[NDIM * unsortedWoundIndices[i]] << '\t' << x[NDIM * unsortedWoundIndices[i] + 1] << '\n';
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
  NewmanZiff(ptr, empty, mode, big);

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
  bool stopSignal = false, tempWoundIsClosed = false;
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
        if (checkWoundClosed(temporaryWoundIndices)) {
          tempWoundIsClosed = true;
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
    if (tempWoundIsClosed)
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
  if (checkWoundClosed(temporaryWoundIndices)) {
    unsortedWoundIndices = temporaryWoundIndices;
  } else {
    // cout << "temporaryWoundIndices was not a closed shape, use old shape for unsortedWoundIndices.\n";
  }
}

bool epi2D::checkWoundClosed(std::vector<int>& listOfIndices) {  // true if listOfIndices (wound) is closed, false if not
  if (listOfIndices.size() == 0)
    return false;
  int firstVertex = listOfIndices[0];
  int lastVertex = listOfIndices.back();
  return (std::find(vnn[firstVertex].begin(), vnn[firstVertex].end(), lastVertex) !=
          vnn[firstVertex].end());
}

int epi2D::findRoot(int i, std::vector<int>& ptr) {
  if (ptr[i] < 0)
    return i;
  return ptr[i] = findRoot(ptr[i], ptr);
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

void epi2D::purseStringContraction(double trate) {
  woundIsClosed = checkWoundClosed(unsortedWoundIndices);
  // cout << "wound is closed? = " << woundIsClosed << '\n';
  //    in epithelia, actomyosin accumulates in cells at wound edges, forms a ring
  //    that shrinks - sewing the wound shut model by shrinking vertices and length
  //    between vertices on the wound edge, pulls cells into wound by cortical
  //    tension
  //     no PBC with wound healing, so PBC aren't accounted for here
  int empty = -NVTOT - 1;
  int mode, big = 0;
  int ci, cj, vi, vj, gi;
  std::vector<int> ptr(NVTOT, empty);
  std::vector<int> deleteList;
  // (squared) distance between vertices i and j, and 1/2 the contact distance
  // between vertices
  double dij, rij;
  order = refineBoundaries();
  NewmanZiff(ptr, empty, mode, big);

  // sort clusters, pick the 2nd largest to get the wound segment
  std::vector<int> clusterSize;
  std::vector<int> clusterRoot;
  for (int gi = 0; gi < NVTOT; gi++) {
    if (ptr[gi] < 0 && ptr[gi] != empty) {
      clusterSize.push_back(-ptr[gi]);
      clusterRoot.push_back(gi);
    }
  }

  sort(clusterSize.begin(), clusterSize.end(), greater<int>());

  int nthLargestCluster = 2;
  int nthClusterSize = clusterSize[nthLargestCluster - 1];
  for (auto i : clusterRoot) {
    if (-ptr[i] == nthClusterSize)
      mode = i;
  }

  getWoundVertices(nthLargestCluster);

  std::vector<std::vector<int>> woundVerts(NCELLS);  //, vector<int> (int (NVTOT/NCELLS), -1));

  // get a list of all indices for cells with wound-facing vertices
  std::vector<int> woundCells;
  for (int gi = 0; gi < NVTOT; gi++) {
    if (std::find(unsortedWoundIndices.begin(), unsortedWoundIndices.end(), gi) != unsortedWoundIndices.end()) {
      cindices(ci, vi, gi);
      woundCells.push_back(ci);
    }
  }

  std::sort(woundCells.begin(), woundCells.end());
  auto last = std::unique(woundCells.begin(), woundCells.end());
  woundCells.erase(last, woundCells.end());

  // get a list of all indices for wound-facing vertices in a given cell ci
  for (auto ci : woundCells) {
    for (int vi = 0; vi < nv[ci]; vi++) {
      int gi = gindex(ci, vi);
      if (std::find(unsortedWoundIndices.begin(), unsortedWoundIndices.end(), gi) != unsortedWoundIndices.end()) {
        woundVerts[ci].push_back(gi);
      }
    }
  }

  // tol represents how much the vertices in the same cell are allowed to
  // overlap (tol * contact distance) before we start deleting them
  double tol = 0.7;
  std::vector<int> deleteV(1, 0);
  // shrunkV keeps track of shrunk vertices, so I can redistribute perimeter to
  // unshrunk vertices
  std::vector<int> shrunkV;

  // sort unsortedWoundIndices
  std::vector<int> woundVertexList;
  std::vector<int> woundVertexCiCounter(NCELLS, 0);  // counts # of occurrences of each ci in the wound

  woundVertexList.push_back(unsortedWoundIndices[0]);
  cindices(ci, vi, unsortedWoundIndices[0]);
  woundVertexCiCounter[ci]++;

  for (int i = 1; i < unsortedWoundIndices.size(); i++) {
    for (auto j : vnn[unsortedWoundIndices[i - 1]]) {  // check neighbors of latest entry to woundVertexList
      if (std::find(unsortedWoundIndices.begin(), unsortedWoundIndices.end(), j) != unsortedWoundIndices.end() && std::find(woundVertexList.begin(), woundVertexList.end(), j) == woundVertexList.end()) {
        // if a neighbor j is not in woundVertexList but IS in unsortedWoundIndices,
        //    then it's our next vertex (probably). Add to woundVertexList.
        woundVertexList.push_back(j);
        cindices(cj, vj, j);
        woundVertexCiCounter[cj]++;
        break;
      }
    }
  }
  if (simclock < 142.0 && simclock > 140.0) {
    cout << "woundVertexList.size() = " << woundVertexList.size() << "\t, unsortedWoundIndices.size() = " << unsortedWoundIndices.size() << '\n';

    for (auto i : woundVertexList)
      cout << "woundVertex " << i << '\n';
    for (auto i : unsortedWoundIndices)
      cout << "unsortedWoundIndex " << i << '\n';
  }
  /*
  // remove (prune) any gi in woundVertexList if it is the only such vertex of any given cell ci
  for (int i = 0; i < woundVertexList.size(); i++) {
    cindices(ci, vi, woundVertexList[i]);
    if (woundVertexCiCounter[ci] == 1 || woundVertexCiCounter[ci] == 0) {  // 0 or 1 (cell is not/hardly on wound) = remove.  2 or larger (cell is significantly on wound) = stay.
      // 0 case arises when I reject the void segmentation and use an old wound configuration.
      // excommunicate vertex i from the wound
      woundVertexList.erase(woundVertexList.begin() + i);
      woundVertexCiCounter[ci]--;
    }
  }*/

  // cout << "before evaluateGhostDPForces, woundVertexList.size() = " << woundVertexList.size() << '\n';
  evaluateGhostDPForces(woundVertexList, trate);

  //  delete vertices if significant overlap, only if >1 vertex is wound-facing
  for (auto ci : woundCells) {
    int sz = woundVerts[ci].size();
    if (sz > 1) {
      // set up initial perimeter to calculate how much needs to be distributed
      double remainderPerimeter = initialPreferredPerimeter;

      for (int i = 0; i < woundVerts[ci].size(); i++) {
        // calculate new l0 for wound vertices, and keep track of loss of
        // perimeter
        int gi = woundVerts[ci][i];
        // only shrink segments that point towards a wounded vertex
        if (vnn_label[ip1[gi]] != 0 && vnn_label[ip1[gi]] != 2)
          continue;

        // remainderPerimeter -= l0[gi];
        // l0[gi] *= exp(-trate * dt);
        // remainderPerimeter += l0[gi];

        // shrunkV.push_back(gi);

        int gj;  // j is a vertex next to i
        for (int k = -1; k < 2; k += 2) {
          if (szList[ci] == 0)
            gj = (gi + k) % nv[ci];
          else
            gj = (((gi + k) % szList[ci] % nv[ci]) + szList[ci]);
          cindices(cj, vj, gj);

          dij = pow(x[gi * NDIM] - x[gj * NDIM], 2) +
                pow(x[gi * NDIM + 1] - x[gj * NDIM + 1], 2);
          rij = pow((r[gi] + r[gj]) * tol, 2);
          if (dij < rij) {
            std::vector<int> deletedVertices = regridSegment(gi, r[gj]);  // regrid a chain of line segments around vertex gi
            // adjust woundVerts[ci] values to account for newly deleted vertices
            for (auto i : deletedVertices) {
              regridChecker = true;
              for (int j = 0; j < woundVerts.size(); j++) {
                for (auto& k : woundVerts[j]) {
                  if (k > i) {
                    k--;
                  }
                }
              }
            }
            break;
          }
        }
      }
      for (int i = 0; i < nv[ci]; i++) {
        // gi = gindex(ci, i);
        //  double remainderLength = (initialPreferredPerimeter - remainderPerimeter) / nv[ci];
        //   check that gi is not in shrunkV
        /*if (std::find(shrunkV.begin(), shrunkV.end(), gi) == shrunkV.end()) {
          // then give gi remainderPerimeter / (nv[ci] - shrunkV)
          double remainderLength = (initialPreferredPerimeter - remainderPerimeter) / (nv[ci] - shrunkV.size());
          // l0[gi] += remainderLength;
          // r[gi] += remainderLength/2;
        }*/

        // l0[gi] += remainderLength;
        // cout << "simclock = " << simclock << ", remainderLength = " << remainderLength << ", current l0[gi] = " << l0[gi] << '\n';
      }
      // remove contents of shrunkV after finishing with the regridding and reapportionment
      // shrunkV.clear();
    }  // else if (sz == 1 || sz == 2) {
       //  turn off intercellular attraction and repulsion for the vertex in question
       // listTurnOffRepulsion.push_back(woundVerts[ci][0]);
    //}
  }
}

void epi2D::evaluateGhostDPForces(vector<int>& giList, double trate) {
  // create a DP with no preferred area, no repulsions, and no attraction. Its component
  //  vertices are existing vertices within the simulation, that belong to other non-ghost DPs.
  // The ghost DP's purpose is to apply forces (i.e. purse-string forces) between vertices of
  //  different cells.
  // Ghost DP has no preferred shape, no preferred area.
  //
  // input: vector of global indices for vertices of the ghost DP. These global indices already exist
  //    and are pre-sorted in coordinate space such that giList[0] and giList[giList.size()-1] are neighbors

  std::vector<int> ghostVList = giList;
  // cout << "ghostVList.size() = " << ghostVList.size() << '\n';
  /*for (auto i : ghostVList) {
    cout << "vertex gi = " << i << '\n';
    cout << "\t at position = " << x[NDIM * i] << '\t' << x[NDIM * i + 1] << '\n';
  }*/
  int nvGhost = giList.size();
  if (nvGhost == 0)
    cout << "giList is empty, will mod by 0 leading to nan\n";
  std::vector<bool> isCorner;  // vector of ghostVList indices (i.e. 0,1,2,3, etc) that are corners (interact with ghost purse-string)
  std::vector<double> l0Ghost;
  // establish purse-string between each consecutive pair of vertices, treating segments spanning 2 cells specially
  double distBetweeniandj = 0;
  // loop over vertices in ghost wound
  // if corner, set preferred length l0Ghost[i] to be exp(-trate*dt)* current contact length
  //  else it is normal, set preferred length exp(-trate*dt) * l0[gi]
  // then compute shape forces
  for (int i = 0; i < nvGhost; i++) {
    // find cell index of vertex i and j
    int j = (i + 1) % nvGhost;
    int ci, vi, cj, vj;
    int gi = ghostVList[i];
    int gj = ghostVList[j];
    cindices(ci, vi, gi);
    cindices(cj, vj, gj);
    if (ci != cj) {  // then i and j are neighbors and wound corners belonging to different cells
      distBetweeniandj = vertDistNoPBC(gi, gj);
      // isCorner.push_back(i);
      l0Ghost.push_back(exp(-trate * dt) * distBetweeniandj);
    } else {
      // use default preferred segment length if segment is within one cell
      l0Ghost.push_back(exp(-trate * dt) * l0[gi]);
    }
  }

  // already set preferred lengths, so now just compute the shape forces for preferred lengths multiplied by the bool vector
  // check shapeForces2D(), left a draft on Atom

  double fli, flim1, cx, cy, xi, yi, gip1, xip1, yip1;
  double rho0, l0im1, l0i;
  double dx, dy, dli, dlim1;
  double lim1x, lim1y, lix, liy, lip1x, lip1y, li, lim1;
  double rim1x, rim1y, rix, riy, rip1x, rip1y;
  double forceX, forceY;

  rho0 = sqrt(a0.at(0));
  // initialize center of mass coordinates
  xi = x[NDIM * ghostVList[0]];
  yi = x[NDIM * ghostVList[0] + 1];
  cx = xi;
  cy = yi;

  // cerr << rho0 << '\t' << xi << '\n';

  // loop over vertices of cell ci, get perimeter
  for (int vi = 0; vi < nvGhost - 1; vi++) {
    // next vertex
    gip1 = ghostVList[(vi + 1) % nvGhost];

    // get positions (check minimum images)
    dx = x[NDIM * gip1] - xi;
    if (pbc[0])
      dx -= L[0] * round(dx / L[0]);
    xip1 = xi + dx;

    dy = x[NDIM * gip1 + 1] - yi;
    if (pbc[1])
      dy -= L[1] * round(dy / L[1]);
    yip1 = yi + dy;

    // add to center of mass
    cx += xip1;
    cy += yip1;

    // update coordinates
    xi = xip1;
    yi = yip1;
  }

  // cerr << xi << '\n';

  // take average to get com
  cx /= nvGhost;
  cy /= nvGhost;

  // ----

  // loop over vertices, add to force
  for (int i = 0; i < nvGhost; i++) {
    // preferred segment length
    int gi = ghostVList[i];
    int ipi = ghostVList[(i + 1 + nvGhost) % nvGhost];
    int imi = ghostVList[(i - 1 + nvGhost) % nvGhost];
    l0i = l0Ghost[i];
    l0im1 = l0Ghost[(i - 1 + nvGhost) % nvGhost];
    // cout << "l0i = " << l0i << "\t, l0[gi] = " << l0[gi] << '\n';
    l0[gi] = l0i;  // update the l0 list

    // get coordinates relative to center of mass
    rix = x[NDIM * gi] - cx;
    riy = x[NDIM * gi + 1] - cy;

    // get next adjacent vertices
    rip1x = x[NDIM * ipi] - cx;
    rip1y = x[NDIM * ipi + 1] - cy;
    if (pbc[0])
      rip1x -= L[0] * round(rip1x / L[0]);
    if (pbc[1])
      rip1y -= L[1] * round(rip1y / L[1]);

    // get next adjacent vertices
    rim1x = x[NDIM * imi] - cx;
    rim1y = x[NDIM * imi + 1] - cy;
    if (pbc[0])
      rim1x -= L[0] * round(rim1x / L[0]);
    if (pbc[1])
      rim1y -= L[1] * round(rim1y / L[1]);

    // -- Perimeter force
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

    // cerr << "kl, l0im1 = " << kl << '\t' << l0im1 << '\n';

    // cerr << fli << '\t' << dli << '\t' << lix << '\t' << li << '\n';

    // cerr << "forceX = " << forceX << '\n';

    // add to forces
    forceX = (fli * dli * lix / li) - (flim1 * dlim1 * lim1x / lim1);
    forceY = (fli * dli * liy / li) - (flim1 * dlim1 * lim1y / lim1);
    F[NDIM * gi] += forceX;
    F[NDIM * gi + 1] += forceY;

    // update potential energy
    U += 0.5 * kl * (dli * dli);
    // per cell energy will no longer be proportional to total energy due to wound energy
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
      posout << setw(wnum) << setprecision(pnum) << right << r[gi];
      posout << setw(wnum) << setprecision(pnum) << right << l0[gi];
      posout << setw(wnum) << setprecision(pnum) << right << t0[gi];
      posout << endl;
    }
  }

  // print end frame
  posout << setw(w) << left << "ENDFR"
         << " " << endl;

  cout << "leaving printPos\n";
}