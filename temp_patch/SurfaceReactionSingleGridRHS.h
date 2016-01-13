#ifndef _SURFACE_REACTION_SINGLE_GRID_RHS_H_
#define _SURFACE_REACTION_SINGLE_GRID_RHS_H_

#include "DataStructure.h"
#include "SurfaceReactionSolverStruct.h"
#include "SurfaceReactions.h"

/*
 * End of auto-generated surface reactions stoichiometries
 */

PetscErrorCode SurfaceReactionSingleGridRHS(TS ts,PetscReal t,Vec x,Vec f,void *ptr)
{
    SURFACE_REACTION_SOLVER  *solver_ptr = (SURFACE_REACTION_SOLVER *)ptr;  /* user-defined application context */
    PetscReal                reaction_rate[NUMBER_OF_REACTIONS_ON_SURFACE];
    PetscInt                 i, j, number_volume_species, number_surface_species, number_surface_reactions;
    PetscReal                *x_ptr, *cv_ptr, *f_ptr;
    PetscErrorCode           ierr;

    PetscFunctionBegin;

    number_volume_species   = NUMBER_OF_SPECIES_IN_VOLUME;
    number_surface_species   = NUMBER_OF_SPECIES_ON_SURFACE;
    number_surface_reactions = NUMBER_OF_REACTIONS_ON_SURFACE;

/*
 * Get the pointers pointing to data array of x, f, and c_v respectively
 *  x: surface species concentrations
 *  c_v: volume species concentrations
 */
    ierr = VecGetArray(x, &x_ptr);CHKERRQ(ierr);
    ierr = VecGetArray(f, &f_ptr);CHKERRQ(ierr);
    ierr = VecGetArray(solver_ptr->c_v, &cv_ptr);CHKERRQ(ierr);

/*
 *  Calling user-generated function from SurfaceReactions.h
 */
    ierr = SurfaceReactionRateCalculation(t, x_ptr,cv_ptr, reaction_rate); CHKERRQ(ierr);

/*
 *  Compose right-hand-side with reaction_rate from above line and stoichiometry matrix from SurfaceReactions.h
 */
    for(i=0;i<number_surface_species;++i){
        f_ptr[i] = 0.0;
        for(j=0;j<number_surface_reactions;++j){
            if(surface_stoichiometry_surface_species[i][j] != 0){
                 f_ptr[i] += (PetscReal)surface_stoichiometry_surface_species[i][j]*reaction_rate[j];
            }
        }
    }

/*
 * restore vectors x, f, and c_v from data arrays
 */
    ierr = VecRestoreArray(x, &x_ptr);CHKERRQ(ierr);
    ierr = VecRestoreArray(f, &f_ptr);CHKERRQ(ierr);
    ierr = VecRestoreArray(solver_ptr->c_v, &cv_ptr);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

PetscErrorCode SurfaceReactionSingleGridRHSJacobian(TS ts,PetscReal t,Vec x,Mat *J_ptr,Mat *Prec_ptr,MatStructure *flag,void *ptr)
{
    SURFACE_REACTION_SOLVER  *solver_ptr = (SURFACE_REACTION_SOLVER *)ptr;  /* user-defined application context */
    PetscReal                J[NUMBER_OF_SPECIES_ON_SURFACE][NUMBER_OF_SPECIES_ON_SURFACE];
    PetscInt                 i, j, number_surface_species, number_volume_species, number_surface_reactions;
    PetscInt                 rowcol[NUMBER_OF_SPECIES_ON_SURFACE];
    PetscReal                *x_ptr, *cv_ptr;
    PetscErrorCode           ierr;

    PetscFunctionBegin;

    number_volume_species   = NUMBER_OF_SPECIES_IN_VOLUME;
    number_surface_species   = NUMBER_OF_SPECIES_ON_SURFACE;
    number_surface_reactions = NUMBER_OF_REACTIONS_ON_SURFACE;

/*
 * Get the pointers pointing to data array of x and c_v
 */
    ierr = VecGetArray(x, &x_ptr);CHKERRQ(ierr);
    ierr = VecGetArray(solver_ptr->c_v, &cv_ptr);CHKERRQ(ierr);

/*
 * Initialize J to zero matrix
 */
    for(i=0;i<number_surface_species;i++)
       for(j=0;j<number_surface_species;j++)
          J[i][j] = (PetscReal)0.0;

/*
 *  Calling user-generated function from SurfaceReactions.h
 */
    ierr = SurfaceReactionJacobianCalculation(t, x_ptr, cv_ptr, J); CHKERRQ(ierr);

/*
 * insert matrix values from J
 */
    for(i=0;i<number_surface_species;i++)
       rowcol[i] = i;
    ierr = MatSetValues(*Prec_ptr,number_surface_species,rowcol,number_surface_species,rowcol,&J[0][0],INSERT_VALUES);

/*
 * restore vector x and c_v from data arrays
 */
    ierr = VecRestoreArray(x, &x_ptr);CHKERRQ(ierr);
    ierr = VecRestoreArray(solver_ptr->c_v, &cv_ptr);CHKERRQ(ierr);

/*
 * Assembly matrix
 */
    ierr = MatAssemblyBegin(*Prec_ptr,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(*Prec_ptr,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    if (*Prec_ptr != *J_ptr) {
        ierr = MatAssemblyBegin(*J_ptr,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
        ierr = MatAssemblyEnd(*J_ptr,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    }

//    *flag = SAME_PRECONDITIONER;
    *flag = SAME_NONZERO_PATTERN;

    PetscFunctionReturn(0);
}

#endif
