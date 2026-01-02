//
//  grid.cpp
//  n-body
//
//  Created by Luis Paulo Santos on 01/01/2026.
//

#include "grid.hpp"
#include <stdlib.h>
#include <math.h>

/*
  Given a point (x,y,z), return the 8 vertex indices of the
  containing cell.

  Returns 0 on success, -1 if point is outside the grid.
*/
bool cell_vertices_from_point(float x, float y, float z, int v[]) {
    /* Compute cell indices */
    int i = (int)floor((x - xmin) / dx);
    int j = (int)floor((y - ymin) / dy);
    int k = (int)floor((z - zmin) / dz);

    /* Check bounds (cells go from 0 to nx-2, etc.) */
    if (i < 0 || i >= nx - 1 ||
        j < 0 || j >= ny - 1 ||
        k < 0 || k >= nz - 1)
    {
        return false; /* outside grid */
    }

    /* Vertex indices (VTK hex ordering) */
    v[0] = idx(i,   j,   k);
    v[1] = idx(i+1, j,   k);
    v[2] = idx(i+1, j+1, k);
    v[3] = idx(i,   j+1, k);
    v[4] = idx(i,   j,   k+1);
    v[5] = idx(i+1, j,   k+1);
    v[6] = idx(i+1, j+1, k+1);
    v[7] = idx(i,   j+1, k+1);

    return true;
}

bool init_grid (int const _nx, int const _ny, int const _nz,
                float const _xmin, float const _ymin, float const _zmin,
                float const _xmax, float const _ymax, float const _zmax) {
    nx = _nx;
    ny = _ny;
    nz = _nz;
    xmin = _xmin;    ymin = _ymin;    zmin = _zmin;
    xmax = _xmax;    ymax = _ymax;    zmax = _zmax;
    
    /* -----------------------------
           Allocate field: g = (gx,gy,gz)
           ----------------------------- */
    nVertices = nx * ny * nz;
    gx = (float *) malloc(nVertices * sizeof(float));
    gy = (float *) malloc(nVertices * sizeof(float));
    gz = (float *) malloc(nVertices * sizeof(float));

    if (!gx || !gy || !gz) {
        return false;
    }

    dx = (xmax - xmin) / (nx - 1);
    dy = (ymax - ymin) / (ny - 1);
    dz = (zmax - zmin) / (nz - 1);

    return true;
}

