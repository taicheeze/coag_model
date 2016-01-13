#ifndef _VOLUME_REACTION_SOLVER_H_
#define _VOLUME_REACTION_SOLVER_H_

#include "DataStructure.h"
#include "VolumeReactionSolverStruct.h"
#include "VolumeReactionSingleGridRHS.h"

PetscErrorCode InitializeVolumeReactionSolver(VOLUME_REACTION_SOLVER *solver_ptr, MY_DATA *user);
PetscErrorCode DestroyVolumeReactionSolver(VOLUME_REACTION_SOLVER *solver_ptr);

PetscErrorCode VolumeReactionSolverStep(VOLUME_REACTION_SOLVER *solver_ptr, PetscReal t, PetscReal dt, MY_DATA *user);

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 *     Form and set volume reaction TS solver
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
PetscErrorCode InitializeVolumeReactionSolver(VOLUME_REACTION_SOLVER *solver_ptr, MY_DATA *user)
{
    PetscErrorCode ierr;
    PetscInt       number_of_volume_species;

    PetscFunctionBegin;

    number_of_volume_species = NUMBER_OF_SPECIES_IN_VOLUME;

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

//    ierr = KSPSetTolerances(solver_ptr->ksp,PETSC_DEFAULT,0.0,PETSC_DEFAULT,PETSC_DEFAULT);

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 *     Create x_window (window to whole vector) without first assigning data array
 *     and create x (solution vector) and r (residual vector) 
 *     both should be of size number_of_volume_species 
 *     to hold only species concentrations in one grid on each cpu
 *     (species concentration in volume)
 *     Set solution vector to x 
 *     store right hand side in r
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    ierr = VecCreateSeqWithArray(MPI_COMM_SELF,1,number_of_volume_species,PETSC_NULL,&(solver_ptr->x_window));CHKERRQ(ierr);

    ierr = VecCreate(MPI_COMM_SELF, &(solver_ptr->x));CHKERRQ(ierr);
    ierr = VecSetType(solver_ptr->x, VECSEQ);CHKERRQ(ierr);
    ierr = VecSetSizes(solver_ptr->x, number_of_volume_species, number_of_volume_species);CHKERRQ(ierr);

    ierr = TSSetSolution(solver_ptr->ts, solver_ptr->x); CHKERRQ(ierr);

    ierr = VecDuplicate(solver_ptr->x, &solver_ptr->r); CHKERRQ(ierr);

    ierr = TSSetRHSFunction(solver_ptr->ts,solver_ptr->r,VolumeReactionSingleGridRHS,user);CHKERRQ(ierr); // right-hand-side function

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
    ierr = MatSetSizes(solver_ptr->J, number_of_volume_species, number_of_volume_species, number_of_volume_species, number_of_volume_species);CHKERRQ(ierr);
    ierr = MatSetUp(solver_ptr->J);
    TSSetRHSJacobian(solver_ptr->ts, solver_ptr->J, solver_ptr->J, VolumeReactionSingleGridRHSJacobian, user); CHKERRQ(ierr);

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 *    Set up TS
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    ierr = TSSetUp(solver_ptr->ts);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 *     Destroy volume reaction TS solver
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
PetscErrorCode DestroyVolumeReactionSolver(VOLUME_REACTION_SOLVER *solver_ptr)
{
    PetscErrorCode ierr;

    PetscFunctionBegin;

    ierr = VecDestroy(&(solver_ptr->x)); CHKERRQ(ierr);
    ierr = VecDestroy(&(solver_ptr->r)); CHKERRQ(ierr);
    ierr = MatDestroy(&(solver_ptr->J)); CHKERRQ(ierr);
    ierr = TSDestroy(&(solver_ptr->ts)); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 *     Execute one volume reaction step at each grid
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
PetscErrorCode VolumeReactionSolverStep(VOLUME_REACTION_SOLVER *solver_ptr, PetscReal t, PetscReal dt, MY_DATA *user)
{
    PetscErrorCode   ierr;
    PetscInt         ys,ym,ye,xs,xm,xe;
    PetscInt         i,j;
//    PetscReal        minimum, current_dt;  // for negative populatino control
//    PetscInt         step_num, step_count; // for negative population control
    Species_Volume   **array_c;

    PetscFunctionBegin;

/*
 *   Get a pointer to global vector data.
 */
    ierr = DMDAVecGetArray(user->da_c_volume, user->c_volume_global, &array_c);CHKERRQ(ierr);

/*
 *  Get local grid boundaries (for 2-dimensional DMDA):
 *   xs, ys   - starting grid indices (no ghost points)
 *   xm, ym   - widths of local grid (no ghost points)
 */
    ierr = DMDAGetCorners(user->da_c_volume, &xs, &ys, PETSC_NULL, &xm, &ym, PETSC_NULL);CHKERRQ(ierr);

    ye = ys + ym;
    xe = xs + xm;

/* 
 * Solve volume reaction time-stepping over the locally owned part of the grid
 * 
 */
    for (j=ys; j<ye; j++) {
        for (i=xs; i<xe; i++) {
//            current_dt = dt;
//            step_num = 1;
            ierr = VecPlaceArray(solver_ptr->x_window,array_c[j][i].c); CHKERRQ(ierr);
            ierr = VecCopy(solver_ptr->x_window, solver_ptr->x); CHKERRQ(ierr);
            ierr = TSSetTime(solver_ptr->ts,t);CHKERRQ(ierr);
            ierr = TSSetTimeStep(solver_ptr->ts,dt);CHKERRQ(ierr);
            ierr = TSStep(solver_ptr->ts);CHKERRQ(ierr);
            // check if this step produce negative concentrations, if so, roll back and reduce the stepsize by 1/2
/*            ierr = VecMin(solver_ptr->x,PETSC_NULL,&minimum);
            while(minimum < -CALCULATION_ERROR_TOLERANCE  && step_num < 17){
                 current_dt /= 2.0;
                 step_num *= 2;
                 minimum = 0.0;
                 ierr = VecCopy(solver_ptr->x_window, solver_ptr->x); CHKERRQ(ierr);
                 ierr = TSSetTime(solver_ptr->ts,t);CHKERRQ(ierr);
                 for(step_count=0;step_count<step_num;step_count++){
                     if(minimum > -CALCULATION_ERROR_TOLERANCE){
                        ierr = TSSetTimeStep(solver_ptr->ts,current_dt);CHKERRQ(ierr);
            		ierr = TSStep(solver_ptr->ts);CHKERRQ(ierr);
                        ierr = VecMin(solver_ptr->x,PETSC_NULL,&minimum);
                     }                     
                 }
            }
*/
	    ierr = VecCopy(solver_ptr->x, solver_ptr->x_window); CHKERRQ(ierr);
            ierr = VecResetArray(solver_ptr->x_window); CHKERRQ(ierr);
/*
            if(step_num > 16){
               SETERRQ2(PETSC_COMM_SELF,1,"VolumeReactionSolverStep(): Still produce negative concentration at grid (%d,%d) after reducing step size by 16 times.\n", i, j);
            }
*/
        }
    }

/*
 *  *  Restore array to vector
 *   */
    ierr = DMDAVecRestoreArray(user->da_c_volume, user->c_volume_global, &array_c);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

#endif
