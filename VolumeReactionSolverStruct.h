#ifndef _VOLUME_REACTION_SOLVER_STRUCT_H_
#define _VOLUME_REACTION_SOLVER_STRUCT_H_

#include "DataStructure.h"

typedef struct {
   TS    ts;   /* time-stepping solver context */
   SNES  snes; /* nonlinear solver context */
   KSP   ksp;  /* linear solver context */
   PC    pc;   /* preconditioner context */
   Vec   x, r;    /* single grid concentration vector and residual vector */
   Vec   x_window;/* x_window will be created with VecCreateSeqWithArray() without assigning any array, 
                   * but rather point to corresponding array memory of c_volume_global in MY_DATA for each grid
                   */
   Mat   J;    /* jacobian matrix */
} VOLUME_REACTION_SOLVER;

#endif
