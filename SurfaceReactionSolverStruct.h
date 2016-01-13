#ifndef _SURFACE_REACTION_SOLVER_STRUCT_H_
#define _SURFACE_REACTION_SOLVER_STRUCT_H_

#include "DataStructure.h"

typedef struct {
   TS    ts;   /* time-stepping solver context */
   SNES  snes; /* nonlinear solver context */
   KSP   ksp;  /* linear solver context */
   PC    pc;   /* preconditioner context */
   Vec   x, r;    /* single grid surface species concentration vector and residual vector */
                  /* x will be created with VecCreateSeqWithArray() without assigning any array, 
                   * but rather point to corresponding array memory of c_surface_global in MY_DATA for each surface grid
                   */
   Vec   c_v;    /* single grid volume species concentration */
                  /* x will be created with VecCreateSeqWithArray() without assigning any array, 
                   * but rather point to corresponding array memory of c_volume_global in MY_DATA for each surface grid
                   */
   Mat   J;    /* jacobian matrix */
} SURFACE_REACTION_SOLVER;

#endif
