#ifndef _SURFACE_REACTION_SOLVER_H_
#define _SURFACE_REACTION_SOLVER_H_

#include "DataStructure.h"
#include "SurfaceReactionSolverStruct.h"
#include "SurfaceReactionSingleGridRHS.h"

PetscErrorCode InitializeSurfaceReactionSolver(SURFACE_REACTION_SOLVER *solver_ptr, MY_DATA *user);
PetscErrorCode DestroySurfaceReactionSolver(SURFACE_REACTION_SOLVER *solver_ptr);

PetscErrorCode SurfaceReactionSolverStep(SURFACE_REACTION_SOLVER *solver_ptr, PetscReal t, PetscReal dt, MY_DATA *user);

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 *     Form and set surface reaction TS solver
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
PetscErrorCode InitializeSurfaceReactionSolver(SURFACE_REACTION_SOLVER *solver_ptr, MY_DATA *user)
{
    PetscErrorCode ierr;
    PetscInt       number_of_volume_species, number_of_surface_species;

    PetscFunctionBegin;

    number_of_volume_species = NUMBER_OF_SPECIES_IN_VOLUME;
    number_of_surface_species = NUMBER_OF_SPECIES_ON_SURFACE;

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 *     Create TS to manage hierarchical solvers
 *     and get the SNES in TS, KSP in SNES, PC in KSP
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    ierr = TSCreate(PETSC_COMM_SELF,&solver_ptr->ts);CHKERRQ(ierr);
    ierr = TSGetSNES(solver_ptr->ts,&(solver_ptr->snes));CHKERRQ(ierr);
    ierr = SNESGetKSP(solver_ptr->snes,&(solver_ptr->ksp));CHKERRQ(ierr);
    ierr = KSPGetPC(solver_ptr->ksp,&(solver_ptr->pc));CHKERRQ(ierr);

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 *     Set TS properties
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    ierr = TSSetProblemType(solver_ptr->ts,TS_NONLINEAR);CHKERRQ(ierr); // non-linear problem
    ierr = TSSetType(solver_ptr->ts,TSBEULER);CHKERRQ(ierr); // backward euler

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 *     Create x (solution vector) without first assigning data array
 *     and create r (residual vector) 
 *     both should be of size number_of_surface_species 
 *     to hold only species concentrations in one grid on each cpu
 *     (species concentration in volume)
 *     Set solution vector to x 
 *     (remmeber x is just a "single grid window" with no internal data)
 *     store right hand side in r
 *
 *     create c_v (surface species vector) in the same fashion as r
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    ierr = VecCreateSeqWithArray(MPI_COMM_SELF,1,number_of_surface_species,PETSC_NULL,&(solver_ptr->x));CHKERRQ(ierr);

    ierr = TSSetSolution(solver_ptr->ts, solver_ptr->x); CHKERRQ(ierr);

    ierr = VecCreate(MPI_COMM_SELF, &(solver_ptr->r));CHKERRQ(ierr);
    ierr = VecSetType(solver_ptr->r, VECSEQ);CHKERRQ(ierr);
    ierr = VecSetSizes(solver_ptr->r, number_of_surface_species, number_of_surface_species);CHKERRQ(ierr);

    ierr = VecCreate(MPI_COMM_SELF, &(solver_ptr->c_v));CHKERRQ(ierr);
    ierr = VecSetType(solver_ptr->c_v, VECSEQ);CHKERRQ(ierr);
    ierr = VecSetSizes(solver_ptr->c_v, number_of_volume_species, number_of_volume_species);CHKERRQ(ierr);

    ierr = TSSetRHSFunction(solver_ptr->ts,solver_ptr->r,SurfaceReactionSingleGridRHS,solver_ptr);CHKERRQ(ierr); // right-hand-side function

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 *   Create matrix data structure; set Jacobian evaluation routine
 *
 *    Set Jacobian matrix data structure and default Jacobian evaluation
 *    routine. User can override with:
 *     -snes_mf : matrix-free Newton-Krylov method with no preconditioning
 *               (unless user explicitly sets preconditioner)
 *     -snes_mf_operator : form preconditioning matrix as set by the user,
 *                         but use matrix-free approx for Jacobian-vector
 *                         products within Newton-Krylov method
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    ierr = MatCreate(MPI_COMM_SELF, &(solver_ptr->J));CHKERRQ(ierr);
    ierr = MatSetType(solver_ptr->J, MATSEQDENSE);CHKERRQ(ierr);
    ierr = MatSetSizes(solver_ptr->J, number_of_surface_species, number_of_surface_species, number_of_surface_species, number_of_surface_species);CHKERRQ(ierr);
    ierr = MatSetUp(solver_ptr->J);
    TSSetRHSJacobian(solver_ptr->ts, solver_ptr->J, solver_ptr->J, SurfaceReactionSingleGridRHSJacobian, solver_ptr); CHKERRQ(ierr);

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 *    Set up TS
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    ierr = TSSetUp(solver_ptr->ts);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 *     Destroy surface reaction TS solver
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
PetscErrorCode DestroySurfaceReactionSolver(SURFACE_REACTION_SOLVER *solver_ptr)
{
    PetscErrorCode ierr;

    PetscFunctionBegin;

    ierr = VecDestroy(&(solver_ptr->x)); CHKERRQ(ierr);
    ierr = VecDestroy(&(solver_ptr->r)); CHKERRQ(ierr);
    ierr = VecDestroy(&(solver_ptr->c_v)); CHKERRQ(ierr);
    ierr = MatDestroy(&(solver_ptr->J)); CHKERRQ(ierr);
    ierr = TSDestroy(&(solver_ptr->ts)); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 *     Execute one surface reaction step at each surface grid
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
PetscErrorCode SurfaceReactionSolverStep(SURFACE_REACTION_SOLVER *solver_ptr, PetscReal t, PetscReal dt, MY_DATA *user)
{
    PetscErrorCode   ierr;
    PetscInt         xs,xm,xe;
    PetscInt         i,k,number_of_surface_species,number_of_volume_species;
    Species_Volume   **array_c_volume;
    Species_Surface  *array_c_surface;
    PetscReal        **array_dy_volume, dy0, dy1, dy_ratio;
    PetscReal        *array_single_grid_c_v;

    PetscFunctionBegin;

    number_of_volume_species = NUMBER_OF_SPECIES_IN_VOLUME;
    number_of_surface_species = NUMBER_OF_SPECIES_ON_SURFACE;

/*
 *   Get a pointer to global vector data.
 */
    ierr = DMDAVecGetArray(user->da_c_volume, user->c_volume_global, &array_c_volume);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->da_c_surface, user->c_surface_global, &array_c_surface);CHKERRQ(ierr);

    ierr = DMDAVecGetArray(user->da_volume_dof_1, user->dy_volume_global, &array_dy_volume);CHKERRQ(ierr);

/*
 *  Get local grid boundaries (for 1-dimensional DMDA):
 *   xs   - starting grid indices (no ghost points)
 *   xm   - widths of local grid (no ghost points)
 */
    ierr = DMDAGetCorners(user->da_c_surface, &xs, PETSC_NULL, PETSC_NULL, &xm, PETSC_NULL, PETSC_NULL);CHKERRQ(ierr);

    xe = xs + xm;

/* 
 * Solve surface reaction time-stepping over the locally owned part of the grid
 * 
 */
    for (i=xs; i<xe; i++) {
        // place x to reflect surface species in single grid
        ierr = VecPlaceArray(solver_ptr->x,array_c_surface[i].c); CHKERRQ(ierr);

        // calculate c_v with linear extrapolation
        dy0 = array_dy_volume[0][i]; // y=0, boundary layer
        dy1 = array_dy_volume[1][i]; // y=1, one layer up
        dy_ratio = dy0 / (dy0 + dy1); 
        ierr = VecGetArray(solver_ptr->c_v,&array_single_grid_c_v); CHKERRQ(ierr);
        // c_boundary = ( 1 + dy0/(dy0+dy1) ) * c_0 - dy0/(dy0+dy1) * c_1
        for (k=0; k<number_of_volume_species; ++k) {
            array_single_grid_c_v[k] = ((PetscReal)1.0+dy_ratio) * array_c_volume[0][i].c[k] - dy_ratio * array_c_volume[1][i].c[k];
        }
        ierr = VecRestoreArray(solver_ptr->c_v,&array_single_grid_c_v); CHKERRQ(ierr);

        ierr = TSSetTime(solver_ptr->ts,t);CHKERRQ(ierr);
        ierr = TSSetTimeStep(solver_ptr->ts,dt);CHKERRQ(ierr);

        ierr = TSStep(solver_ptr->ts);CHKERRQ(ierr);
        ierr = VecResetArray(solver_ptr->x); CHKERRQ(ierr);
    }

/*
 *  *  Restore array to vector
 *   */
    ierr = DMDAVecRestoreArray(user->da_c_volume, user->c_volume_global, &array_c_volume);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->da_c_surface, user->c_surface_global, &array_c_surface);CHKERRQ(ierr);

    ierr = DMDAVecRestoreArray(user->da_volume_dof_1, user->dy_volume_global, &array_dy_volume);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

#endif
