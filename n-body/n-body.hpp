//
//  n-body.hpp
//  n-body
//
//  Created by Luis Paulo Santos on 03/12/2025.
//

#ifndef n_body_hpp
#define n_body_hpp

#include <stdio.h>

#include "particles.h"

bool nBody_init (PARTICLE &dataSet);
bool nBody_tStep (PARTICLE &dataSet);
void nBody_close (PARTICLE dataSet);

#endif /* n_body_hpp */
