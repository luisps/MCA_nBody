//
//  main.cpp
//  n-body
//
//  Created by Luis Paulo Santos on 03/12/2025.
//

#include <iostream>

#include "n-body.hpp"
#include "grid.hpp"
#include "VTK-Legacy.hpp"
#include "constants.h"

PARTICLE dataSet;

int nx = 0, ny = 0, nz = 0;
float xmin , ymin , zmin , xmax , ymax , zmax;
float dx , dy , dz;

int nVertices;
float *gx, *gy, *gz;

int main(int argc, const char * argv[]) {
    float BBox[6], Gdx, Gdy, Gdz;
    
    // initial conditions
    if (!nBody_init (dataSet, BBox)) return 0;
    Gdx = BBox[3]-BBox[0];
    Gdy = BBox[4]-BBox[1];
    Gdz = BBox[5]-BBox[2];

    if (!init_grid (GRID_RES, GRID_RES, GRID_RES, 
                   BBox[0] - Gdx/10.f, BBox[1] - Gdy/10.f, BBox[2] - Gdz/10.f,
                    BBox[3] + Gdx/10.f, BBox[4] + Gdy/10.f, BBox[5] + Gdz/10.f)) return 0;

    G_field_tStep (dataSet);
    // save initial conditions
    if (!VTK_Legacy_init("particles", "Gfield")) return 0;
    if (!VTK_Legacy_write(0, dataSet)) return 0;
    
    // simulation loop
    for (int tStamp = 1 ; tStamp <= T_STEPS && dataSet.N>1 ; tStamp++) {
        fprintf (stdout, "\r%d", tStamp);
        fflush (stdout);
        
        // simulate time instant
        nBody_tStep(dataSet);
        
        // save data
        if (!VTK_Legacy_write(tStamp, dataSet)) return 0;
    }
    
    nBody_close (dataSet);
    
    printf ("\n\nThat's all, folks!\n");
    
    return 0;
}
