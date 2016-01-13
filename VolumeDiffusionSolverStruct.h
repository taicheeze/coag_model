#ifndef _VOLUME_DIFFUSION_SOLVER_STRUCT_H_
#define _VOLUME_DIFFUSION_SOLVER_STRUCT_H_

#include "DataStructure.h"
#include "BoundaryFluxFactorComputation.h"

typedef struct {
   TS    ts;   /* time-stepping solver context */
   SNES  snes; /* nonlinear solver context */
   KSP   ksp;  /* linear solver context */
   PC    pc;   /* preconditioner context */
   Vec   x_global, r_global;    /* solution vector and residual vector for a single species*/
   Vec   x_global_backup;    /* solution vector and residual vector for a single species*/
//   Vec   x_local, r_local;    /* solution vector and residual vector for a single species, local vector for ghost cell access*/
   Mat   J;    /* jacobian matrix for a single species */
   DM    da_diffusion_operator; /* DMDA for diffusion operator */
   Vec   diffusion_operator;  /* diffusion operator vector generated from da_diffusion_operator */
   BOUNDARY_FLUX_FACTOR_STRUCT bffs; /* boundary flux factors solver */
   PetscInt currentSpecies;  /* used to pass which species we are solving now*/
} VOLUME_DIFFUSION_SOLVER;

// please refer to VolumeDiffusionOperator.h
//             _____.____
//             |         |
//             *  iwup1  *             dy_1
//    ____.____|____.____|____.____
//   |         |         |         |
//   * iwleft1 *   iw0   * iwright1*   dy_0
//   |____.____|____.____|____.____|
//             |         |
//             * iwdown1 *             dy_-1
//             |____.____|
//             |         |
//             * iwdown2 *             dy_-2
//             |____.____|
//       dx        dx         dx
//
#define VOLUME_DIFFUSION_STENCIL_SIZE 6

typedef struct {
   PetscReal iwup1, iw0, iwdown1, iwdown2, iwleft1, iwright1;
} IMPLICIT_WEIGHTS;

#endif
