//
//  particles.h
//  n-body
//
//  Created by Luis Paulo Santos on 03/12/2025.
//

#ifndef particles_h
#define particles_h

typedef struct {
    int N;                  // number of particles
    float *Px, *Py, *Pz;    // arrays of particles' positions
    float *mass, *radius;   // array  of particles' masses
    float *Vx, *Vy, *Vz;    // arrays of particles' velocities
} PARTICLE;

#endif /* particles_h */
