//
//  VTK-Legacy.hpp
//  n-body
//
//  Created by Luis Paulo Santos on 03/12/2025.
//

#ifndef VTK_Legacy_hpp
#define VTK_Legacy_hpp

#include <stdio.h>
#include "particles.h"

bool VTK_Legacy_init (const char *filenameP, const char *filenameG);
bool VTK_Legacy_write  (int const t_stamp, PARTICLE const p);

#endif /* VTK_Legacy_hpp */
