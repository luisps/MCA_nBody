//
//  particles.h
//  n-body
//
//  Created by Luis Paulo Santos on 03/12/2025.
//

#ifndef particles_h
#define particles_h

typedef struct {
    int N;
    float *Px, *Py, *Pz;
    float *mass;
    float *Vx, *Vy, *Vz;
} PARTICLE;

#endif /* particles_h */
