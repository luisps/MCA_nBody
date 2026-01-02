//
//  n-body.cpp
//  n-body
//
//  Created by Luis Paulo Santos on 03/12/2025.
//

#include "n-body.hpp"
#include "constants.h"
#include "grid.hpp"
#include <stdlib.h>
#include <random>
#include <cmath>

#define EPS 1e-2
#define EPS_SQ (EPS*EPS)

inline static void polar2cartesian (float const D, float const O, float const E, float &Px, float &Py, float &Pz) {
    Px = D * sin(E) * cos(O);
    Py = D * sin(E) * sin(O);
    Pz = D * cos(E);
}

inline static float mass2radius(float mass) {
    return powf (3.f * mass /(4.f * M_PI), 1.f/3.f);
}

static void nBody_merge (PARTICLE &dataSet) {
    bool merge = true;
    int const N = dataSet.N;
    int new_N = N;
    
    while (merge) {
        merge = false;
        // for each particle
        for (int i=0 ; i < new_N ; i++) {
            if (dataSet.mass[i] <= 0.f) continue;
            
            float const Radius_i = dataSet.radius[i];
            for (int j=i+1 ; j < new_N ; j++) {
                if (dataSet.mass[j] <= 0.f) continue;
                // compute the direction vector
                float dirX = dataSet.Px[j] - dataSet.Px[i];
                float dirY = dataSet.Py[j] - dataSet.Py[i];
                float dirZ = dataSet.Pz[j] - dataSet.Pz[i];
                
                float dist_squared =
                    (i==j ? 1.f : dirX*dirX + dirY*dirY + dirZ*dirZ);
                
                float const sumRadius = Radius_i + dataSet.radius[j];
                
                if (dist_squared < (sumRadius*sumRadius)) { // merge
                    merge = true;
                    
                    dataSet.mass[i] += dataSet.mass[j];
                    float const weight_j = dataSet.mass[j] / dataSet.mass[i];
                    dataSet.mass[j] = 0.f;
                    dataSet.radius[i] = mass2radius(dataSet.mass[i]);
                    
                    dataSet.Px[i] = (1-weight_j)*dataSet.Px[i] + weight_j*dataSet.Px[j];
                    dataSet.Py[i] = (1-weight_j)*dataSet.Py[i] + weight_j*dataSet.Py[j];
                    dataSet.Pz[i] = (1-weight_j)*dataSet.Pz[i] + weight_j*dataSet.Pz[j];
                    dataSet.Vx[i] = (1-weight_j)*dataSet.Vx[i] + weight_j*dataSet.Vx[j];
                    dataSet.Vy[i] = (1-weight_j)*dataSet.Vy[i] + weight_j*dataSet.Vy[j];
                    dataSet.Vz[i] = (1-weight_j)*dataSet.Vz[i] + weight_j*dataSet.Vz[j];
                    
                    // in place: copy last particle to [j]
                    if (j < (new_N-1)) {
                        dataSet.mass[j] = dataSet.mass[new_N-1];
                        dataSet.mass[new_N-1] = 0.f;
                        dataSet.radius[j] = dataSet.radius[new_N-1];
                        dataSet.Px[j] = dataSet.Px[new_N-1];
                        dataSet.Py[j] = dataSet.Py[new_N-1];
                        dataSet.Pz[j] = dataSet.Pz[new_N-1];
                        dataSet.Vx[j] = dataSet.Vx[new_N-1];
                        dataSet.Vy[j] = dataSet.Vy[new_N-1];
                        dataSet.Vz[j] = dataSet.Vz[new_N-1];
                        j--;
                    }
                    new_N--;
                }
            }
        }
    }
    dataSet.N = new_N;
}

bool nBody_init (PARTICLE &dataSet, float *BBox) {
    
    int const N = dataSet.N = N_initial;
    
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
    dataSet.radius = (float *) malloc (N*sizeof(float));
    if (dataSet.radius==NULL) {
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
    dataSet.radius[0] = SunRadius;

    // generate random particles
    // Seed with a real random value, if available
    std::random_device r;
    // Get a random engine and seed it
    std::default_random_engine e1(r());

    // Distribution for mass
    std::normal_distribution<> normal_dist_mass(meanMass, varMass);

    float const Mean_DistanceC = SunRadius * Mean_Distance2Sun;
    float const Var_DistanceC = SunRadius * Var_Distance2Sun;
    std::normal_distribution<> normal_dist_DistanceC(Mean_DistanceC, Var_DistanceC);

    // Distribution for orientation angle
    float const Min_Orientation = 0.0f;
    float const Max_Orientation = 2.0f * M_PI;
    std::uniform_real_distribution<float> uniform_dist_Orientation(Min_Orientation, Max_Orientation);

    // Distribution for elevation angle
    float const Mean_Elevation = M_PI / 2.f ;
    float const Var_Elevation = M_PI / 20.f ;
    std::normal_distribution<float> normal_dist_Elevation(Mean_Elevation, Var_Elevation);
    
    // compute the particles BBox
    BBox[0] = BBox[1] = BBox[2] = 0.f;
    BBox[3] = BBox[4] = BBox[5] = 0.f;
    
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
        
        // BBox minimum vertex
        BBox[0] = (Px< BBox[0] ? Px : BBox[0]);
        BBox[1] = (Py< BBox[1] ? Py : BBox[1]);
        BBox[2] = (Pz< BBox[2] ? Pz : BBox[2]);
        // BBox maximum vertex
        BBox[3] = (Px> BBox[3] ? Px : BBox[3]);
        BBox[4] = (Py> BBox[4] ? Py : BBox[4]);
        BBox[5] = (Pz> BBox[5] ? Pz : BBox[5]);

        dataSet.mass[i] = mass;
        dataSet.radius[i] = mass2radius(mass);
        
        // velocity vector
        // tangent to the sphere at (Px,Py,Pz)
        // contained in XY
        // orbital speed (velocity magnitude)
        // V = sqrtf (GM/r)
        // G is 6.674e-11 , but here we allow it to
        // be set by the user
        float const V_magnitude = 1.2f * sqrtf(G*sunMass/D);

        // normalized magnitude (V=1)
        float const norm = sqrtf(Px*Px + Py*Py);
        float const norm_factor = V_magnitude / norm;
        dataSet.Vx[i] = -Py * norm_factor;
        dataSet.Vy[i] = Px * norm_factor;
        dataSet.Vz[i] = 0.f;
        

    }
    
    // merge overlapping particles
    nBody_merge(dataSet);
    
    return true;
}

void nBody_close (PARTICLE dataSet) {
    if (dataSet.Px) free(dataSet.Px);
    if (dataSet.Py) free(dataSet.Py);
    if (dataSet.Pz) free(dataSet.Pz);
    if (dataSet.mass) free(dataSet.mass);
    if (dataSet.radius) free(dataSet.radius);
    if (dataSet.Vx) free(dataSet.Vx);
    if (dataSet.Vy) free(dataSet.Vy);
    if (dataSet.Vz) free(dataSet.Vz);
}

/* -----------------------------
       Compute gravitational field
       ----------------------------- */
void G_field_tStep (PARTICLE &dataSet) {
    int const N = dataSet.N;
    
    for (int k = 0; k < nz; ++k) {
        float const z = zmin + k * dz;
        for (int j = 0; j < ny; ++j) {
            float const y = ymin + j * dy;
            for (int i = 0; i < nx; ++i) {
                float const x = xmin + i * dx;
                int id = idx(i, j, k);
                
                double gxi = 0.0, gyi = 0.0, gzi = 0.0;
                
                for (int n = 0; n < N; ++n) {
                    float rx = x - dataSet.Px[n];
                    float ry = y - dataSet.Py[n];
                    float rz = z - dataSet.Pz[n];
                    
                    float r2 = rx*rx + ry*ry + rz*rz + EPS_SQ;
                    float r  = sqrt(r2);
                    float inv_r3 = 1.0 / (r2 * r);
                    
                    float coeff = -G * dataSet.mass[n] * inv_r3;
                    
                    gxi += coeff * rx;
                    gyi += coeff * ry;
                    gzi += coeff * rz;
                }
                
                gx[id] = gxi;
                gy[id] = gyi;
                gz[id] = gzi;
            }
        }
    }
}

void nBody_tStep (PARTICLE &dataSet) {
    int g_idxs[8];
    int const N = dataSet.N;
    
    G_field_tStep (dataSet);
    
    // update the particles given the grid
    // for each particle
    for (int i=0 ; i < N ; i++) {
        float const Radius_i = dataSet.radius[i];
        // compute the 8 neighbour points in the grid
        if (!cell_vertices_from_point(dataSet.Px[i], dataSet.Py[i], dataSet.Pz[i], g_idxs)) {
            continue;   // particle out of grid
        }
        // compute field in the particle by averaging over neighbours
        float gX=0.f, gY=0.f, gZ=0.f;
        for (int j=0 ; j < 8 ; j++) {
            
            // add the field due to grid vertex j
            gX += gx[g_idxs[j]];
            gY += gy[g_idxs[j]];
            gZ += gz[g_idxs[j]];
        }
        
        // normalize
        gX /= 8;    gY /= 8;    gZ /= 8;
        
        // update particle i velocity
        dataSet.Vx[i] += gX * DELTA_T;
        dataSet.Vy[i] += gY * DELTA_T;
        dataSet.Vz[i] += gZ * DELTA_T;
        // update particle i position
        dataSet.Px[i] += DELTA_T * dataSet.Vx[i];
        dataSet.Py[i] += DELTA_T * dataSet.Vy[i];
        dataSet.Pz[i] += DELTA_T * dataSet.Vz[i];
        
    }
    
    // merge overlapping particles
    nBody_merge(dataSet);
}

