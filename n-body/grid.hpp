//
//  grid.hpp
//  n-body
//
//  Created by Luis Paulo Santos on 01/01/2026.
//

#ifndef grid_hpp
#define grid_hpp

extern int nx, ny, nz;
extern float xmin , ymin , zmin , xmax , ymax , zmax;
extern float dx , dy , dz;

extern int nVertices;
extern float *gx, *gy, *gz;

/* Indexing helper: (i,j,k) -> linear index */
inline int idx(int const i, int const j, int const k)
{
    return i + nx * (j + ny * k);
}

bool cell_vertices_from_point(float x, float y, float z, int v[]);

bool init_grid (int const _nx, int const _ny, int const _nz,
                float const _xmin, float const _ymin, float const _zmin,
                float const _xmax, float const _ymax, float const _zmax);


#endif /* grid_hpp */
