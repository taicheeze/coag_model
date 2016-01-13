#ifndef _VOLUME_REACTION_SINGLE_GRID_RHS_H_
#define _VOLUME_REACTION_SINGLE_GRID_RHS_H_

#include "DataStructure.h"
#include "VolumeReactions.h"

PetscErrorCode VolumeReactionSingleGridRHS(TS ts,PetscReal t,Vec x,Vec f,void *ptr)
{
//    MY_DATA        *user = (MY_DATA *)ptr;  /* user-defined application context */
    PetscReal      reaction_rate[NUMBER_OF_REACTIONS_IN_VOLUME];
    PetscInt       i, j, number_volume_species, number_volume_reactions;
    PetscReal      *x_ptr, *f_ptr;
    PetscErrorCode ierr;

    PetscFunctionBegin;

    number_volume_species   = NUMBER_OF_SPECIES_IN_VOLUME;
    number_volume_reactions = NUMBER_OF_REACTIONS_IN_VOLUME;

/*
 * Get the pointers pointing to data array of x and f respectively
 */
    ierr = VecGetArray(x, &x_ptr);CHKERRQ(ierr);
    ierr = VecGetArray(f, &f_ptr);CHKERRQ(ierr);

/*
 * Calling user-generated function from VolumeReactions.h
 */
    ierr = VolumeReactionRateCalculation(t,x_ptr,reaction_rate); CHKERRQ(ierr);

/*
 * Compose right-hand-side with reaction_rate from above line and stoichiometry matrix from VolumeReactions.h
 */
    for(i=0;i<number_volume_species;++i){
        f_ptr[i] = 0.0;
        for(j=0;j<number_volume_reactions;++j){
            if(volume_stoichiometry[i][j] != 0){
                 f_ptr[i] += (PetscReal)volume_stoichiometry[i][j]*reaction_rate[j];
            }
        }
    }

/*
 * restore vectors x and f from data arrays
 */
    ierr = VecRestoreArray(x, &x_ptr);CHKERRQ(ierr);
    ierr = VecRestoreArray(f, &f_ptr);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

PetscErrorCode VolumeReactionSingleGridRHSJacobian(TS ts,PetscReal t,Vec x,Mat *J_ptr,Mat *Prec_ptr,MatStructure *flag,void *ptr)
{
//    MY_DATA        *user = (MY_DATA *)ptr;  /* user-defined application context */
    PetscReal      J[NUMBER_OF_SPECIES_IN_VOLUME][NUMBER_OF_SPECIES_IN_VOLUME];
    PetscInt       i, j, number_volume_species, number_volume_reactions;
    PetscInt       rowcol[NUMBER_OF_SPECIES_IN_VOLUME];
    PetscReal      *x_ptr;
    PetscErrorCode ierr;

    PetscFunctionBegin;

    number_volume_species   = NUMBER_OF_SPECIES_IN_VOLUME;
    number_volume_reactions = NUMBER_OF_REACTIONS_IN_VOLUME;

/*
 * Get the pointers pointing to data array of x 
 */
    ierr = VecGetArray(x, &x_ptr);CHKERRQ(ierr);

/*
 * Initialize J to zero matrix
 */
    for(i=0;i<number_volume_species;i++)
       for(j=0;j<number_volume_species;j++)
          J[i][j] = (PetscReal)0.0;

/*
 * Calling user-generated function from VolumeReactions.h
 */
    ierr = VolumeReactionJacobianCalculation(t,x_ptr,J); CHKERRQ(ierr);

/*
 * insert matrix values from J
 */
    for(i=0;i<number_volume_species;i++)
       rowcol[i] = i;
    ierr = MatSetValues(*Prec_ptr,number_volume_species,rowcol,number_volume_species,rowcol,&J[0][0],INSERT_VALUES);

/*
 * restore vector x from data arrays
 */
    ierr = VecRestoreArray(x, &x_ptr);CHKERRQ(ierr);

/*
 * Assembly matrix
 */
    ierr = MatAssemblyBegin(*Prec_ptr,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(*Prec_ptr,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    if (*Prec_ptr != *J_ptr) {
        ierr = MatAssemblyBegin(*J_ptr,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
        ierr = MatAssemblyEnd(*J_ptr,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    }

    *flag = SAME_NONZERO_PATTERN;

    PetscFunctionReturn(0);
}

#endif
