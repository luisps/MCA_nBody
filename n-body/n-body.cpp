//
//  n-body.cpp
//  n-body
//
//  Created by Luis Paulo Santos on 03/12/2025.
//

#include "n-body.hpp"
#include "constants.h"
#include <stdlib.h>
#include <random>
#include <cmath>

inline static void polar2cartesian (float const D, float const O, float const E, float &Px, float &Py, float &Pz) {
    Px = D * sin(E) * cos(O);
    Py = D * sin(E) * sin(O);
    Pz = D * cos(E);
}

inline static float mass2radius(float mass) {
    return powf (3.f * mass /(4.f * M_PI), 1.f/3.f);
}

bool nBody_init (PARTICLE &dataSet) {
    
    dataSet.N = N;
    
    dataSet.Px = (float *) malloc (N*sizeof(float));
    if (dataSet.Px==NULL) {
        return false;
    }
    dataSet.Py = (float *) malloc (N*sizeof(float));
    if (dataSet.Py==NULL) {
        return false;
    }
    dataSet.Pz = (float *) malloc (N*sizeof(float));
    if (dataSet.Pz==NULL) {
        return false;
    }
    dataSet.mass = (float *) malloc (N*sizeof(float));
    if (dataSet.mass==NULL) {
        return false;
    }
    dataSet.Vx = (float *) malloc (N*sizeof(float));
    if (dataSet.Vx==NULL) {
        return false;
    }
    dataSet.Vy = (float *) malloc (N*sizeof(float));
    if (dataSet.Vy==NULL) {
        return false;
    }
    dataSet.Vz = (float *) malloc (N*sizeof(float));
    if (dataSet.Vz==NULL) {
        return false;
    }

    // const particles
    dataSet.Px[0] = dataSet.Py[0] = dataSet.Pz[0] = 0.f;
    dataSet.mass[0] = sunMass;
    dataSet.Vx[0] = dataSet.Vy[0] = dataSet.Vz[0] = 0.f;
    
    float const SunRadius = mass2radius(sunMass);
    
    // generate random particles
    // Seed with a real random value, if available
    std::random_device r;
    // Get a random engine and seed it
    std::default_random_engine e1(r());

    // Distribution for mass
    float const Mean_mass = 3.0f;
    float const Var_mass = 1.0f;
    std::normal_distribution<> normal_dist_mass(Mean_mass, Var_mass);

    // Distribution for distance to the center
    float const Mean_DistanceC = 3.0f * SunRadius;
    float const Var_DistanceC = SunRadius;
    std::normal_distribution<> normal_dist_DistanceC(Mean_DistanceC, Var_DistanceC);

    // Distribution for orientation angle
    float const Min_Orientation = 0.0f;
    float const Max_Orientation = 2.0f * M_PI;
    std::uniform_real_distribution<float> uniform_dist_Orientation(Min_Orientation, Max_Orientation);

    // Distribution for elevation angle
    float const Mean_Elevation = M_PI / 2.f ;
    float const Var_Elevation = M_PI / 20.f ;
    std::normal_distribution<float> normal_dist_Elevation(Mean_Elevation, Var_Elevation);
    
    for (int i=1 ; i<N ; i++) {
        // center of the particle in polar coordinates
        float const D =
            fmax(1.5f*SunRadius,  normal_dist_DistanceC(e1));
        // orientation of the particle in polar coordinates
        float const O = uniform_dist_Orientation(e1);
        // elevation of the particle in polar coordinates
        float const E = normal_dist_Elevation(e1);
        float Px, Py, Pz;

        polar2cartesian (D, O, E, Px, Py, Pz);
        
        float const mass = fmax(normal_dist_mass(e1),0.5f);
        
        dataSet.Px[i] = Px;
        dataSet.Py[i] = Py;
        dataSet.Pz[i] = Pz;
        
        dataSet.mass[i] = mass;
        
        // velocity vector
        // tangent to the sphere at (Px,Py,Pz)
        // contained in XY
        // orbital speed (velocity magnitude)
        // V = sqrtf (GM/r)
        // G is 6.674e-11 , but here we allow it to
        // be set by the user
        float const V_magnitude = sqrtf(G*sunMass/D);

        // normalized magnitude (V=1)
        float const norm = sqrtf(Px*Px + Py*Py);
        float const norm_factor = V_magnitude / norm;
        dataSet.Vx[i] = -Py * norm_factor;
        dataSet.Vy[i] = Px * norm_factor;
        dataSet.Vz[i] = 0.f;
        

    }
    return true;
}

void nBody_close (PARTICLE dataSet) {
    if (dataSet.Px) free(dataSet.Px);
    if (dataSet.Py) free(dataSet.Py);
    if (dataSet.Pz) free(dataSet.Pz);
    if (dataSet.mass) free(dataSet.mass);
    if (dataSet.Vx) free(dataSet.Vx);
    if (dataSet.Vy) free(dataSet.Vy);
    if (dataSet.Vz) free(dataSet.Vz);
}

bool nBody_tStep (PARTICLE &dataSet) {
    
    
    // for each particle
    for (int i=0 ; i < N ; i++) {
        // compute acceleration due to all other particles
        float aX=0.f, aY=0.f, aZ=0.f;
        for (int j=0 ; j < N ; j++) {
            // compute the direction vector
            float dirX = dataSet.Px[j] - dataSet.Px[i];
            float dirY = dataSet.Py[j] - dataSet.Px[i];
            float dirZ = dataSet.Pz[j] - dataSet.Px[i];
            
            float const dist_squared =
                (i==j ? 1.f : dirX*dirX + dirY*dirY + dirZ*dirZ) ;
            
            float const dist = sqrtf(dist_squared);
            
            // normalize direction vector
            dirX /= dist;
            dirY /= dist;
            dirZ /= dist;
            
            float const acc = (i==j ? 0.f : G ) * dataSet.mass[j] / dist_squared;
            
            // add the acceleration due to j onto i
            aX += acc * dirX;
            aY += acc * dirY;
            aZ += acc * dirZ;
        }

        // update particle i position
        dataSet.Px[i] = dataSet.Px[i] + DELTA_T * dataSet.Vx[i] + 0.5f * aX * DELTA_T * DELTA_T;
        dataSet.Py[i] = dataSet.Py[i] + DELTA_T * dataSet.Vy[i] + 0.5f * aY * DELTA_T * DELTA_T;
        dataSet.Pz[i] = dataSet.Pz[i] + DELTA_T * dataSet.Vz[i] + 0.5f * aZ * DELTA_T * DELTA_T;
        
        // update particle i velocity
        dataSet.Vx[i] += aX * DELTA_T;
        dataSet.Vy[i] += aY * DELTA_T;
        dataSet.Vz[i] += aZ * DELTA_T;
    }
    
    return true;
}
