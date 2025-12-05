//
//  main.cpp
//  n-body
//
//  Created by Luis Paulo Santos on 03/12/2025.
//

#include <iostream>

#include "n-body.hpp"
#include "VTK-Legacy.hpp"
#include "constants.h"

PARTICLE dataSet;

int main(int argc, const char * argv[]) {

    // initial conditions
    if (!nBody_init (dataSet)) return 0;
    // save initial conditions
    if (!VTK_Legacy_init("particles")) return 0;
    if (!VTK_Legacy_write(0, dataSet)) return 0;
    
    // simulation loop
    for (int tStamp = 1 ; tStamp <= T_STEPS ; tStamp++) {
        fprintf (stdout, "\r%d", tStamp);
        fflush (stdout);
        
        // simulate time instant
        if (!nBody_tStep(dataSet)) return 0;
        
        // save data
        if (!VTK_Legacy_write(tStamp, dataSet)) return 0;
    }
    
    nBody_close (dataSet);
    
    printf ("\n\nThat's all, folks!\n");
    
    return 0;
}
