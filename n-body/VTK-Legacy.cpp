//
//  VTK-Legacy.cpp
//  n-body
//
//  Created by Luis Paulo Santos on 03/12/2025.
//

#include "grid.hpp"
#include "VTK-Legacy.hpp"

#include <string.h>

#define N_CHARS 128

static char fn_rootP[N_CHARS], fn_rootG[N_CHARS];

bool VTK_Legacy_init (const char *filenameP, const char *filenameG) {
    
    strncpy(fn_rootP, filenameP, N_CHARS-15);
    strncpy(fn_rootG, filenameG, N_CHARS-15);

    return true;
}

bool VTK_Legacy_write  (int const t_stamp, PARTICLE const p) {
    
    char fn[N_CHARS];
    FILE *f;
    
    snprintf(fn, N_CHARS-1, "%s_%.8d.vtk", fn_rootP, t_stamp);
    
    f = fopen(fn, "wt");
    if (f==NULL) return false;
    
    // header
    fprintf (f, "# vtk DataFile Version 2.0\n");
    fprintf (f, "Particle sim t=%.8d\n", t_stamp);
    fprintf (f, "ASCII\n");
    
    // particles data set
    fprintf (f, "DATASET POLYDATA\n");
    fprintf (f, "POINTS %d float\n", p.N);
    for (int i = 0 ; i < p.N ; i++) {
        fprintf (f, "%f %f %f\n", p.Px[i], p.Py[i], p.Pz[i]);
    }
    // particles attributtes
    fprintf (f, "POINT_DATA %d\n", p.N);
    // mass
    fprintf (f, "SCALARS mass float 1\n");
    fprintf (f, "LOOKUP_TABLE default\n");
    for (int i = 0 ; i < p.N ; i++) {
        fprintf (f, "%f\n", p.mass[i]);
    }
    // velocity
    fprintf (f, "VECTORS V float\n");
    for (int i = 0 ; i < p.N ; i++) {
        fprintf (f, "%f %f %f\n", p.Vx[i], p.Vy[i], p.Vz[i]);
    }
    fclose(f);
    
    snprintf(fn, N_CHARS-1, "%s_%.8d.vtk", fn_rootG, t_stamp);
    f = fopen(fn, "wt");
    if (f==NULL) return false;

    // header
    fprintf (f, "# vtk DataFile Version 2.0\n");
    fprintf (f, "G field sim t=%.8d\n", t_stamp);
    fprintf (f, "ASCII\n");

    fprintf(f, "DATASET STRUCTURED_POINTS\n");
    fprintf(f, "DIMENSIONS %d %d %d\n", nx, ny, nz);
    fprintf(f, "ORIGIN %f %f %f\n", xmin, ymin, zmin);
    fprintf(f, "SPACING %f %f %f\n", dx, dy, dz);
    fprintf(f, "POINT_DATA %d\n", nVertices);
    fprintf(f, "VECTORS g float\n");

    for (int id = 0; id < nVertices; ++id) {
        fprintf(f, "%e %e %e\n", gx[id], gy[id], gz[id]);
    }

    fclose(f);

    return true;
}
