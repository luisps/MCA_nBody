//
//  constants.h
//  n-body
//
//  Created by Luis Paulo Santos on 04/12/2025.
//

#ifndef constants_h
#define constants_h

// gravitational constant
float const G = 6.674e-11 ;

// number of particles
int const N_initial = 2000;

// Sun mass
float const sunMass = 2000.f;

// Other bodies mean mass
float const meanMass = 1.f;
float const varMass = .25f;

// these constants are multiplied by the sun radius
// to distribute the particles distance to the SUn's center
float const Mean_Distance2Sun = 7.f;
float const Var_Distance2Sun = 1.f;

// Delta t per iteration
float const DELTA_T = 3000.f;

// number of simulation iterations
int const T_STEPS = 1000;

// note that the number of cells is GRID_RES^3
int const GRID_RES = 20;


#endif /* constants_h */
