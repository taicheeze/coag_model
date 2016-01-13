#ifndef _MAX_CHANGE_TIME_STEP_CONTROLLER_H_
#define _MAX_CHANGE_TIME_STEP_CONTROLLER_H_

#include "DataStructure.h"

typedef struct {
    PetscReal      old_max_c_volume[NUMBER_OF_SPECIES_IN_VOLUME], current_max_c_volume[NUMBER_OF_SPECIES_IN_VOLUME];
    PetscReal      old_max_c_surface[NUMBER_OF_SPECIES_ON_SURFACE], current_max_c_surface[NUMBER_OF_SPECIES_ON_SURFACE];
    PetscReal      tmp_t, tmp_dt;
    PetscReal      new_t, new_dt;
    PetscInt       flg;
} TIME_STEP_CONTROLLER; 

PetscErrorCode InitializeTimeStepController(PetscReal t, PetscReal dt, TIME_STEP_CONTROLLER *controller, MY_DATA *user);
PetscErrorCode SingleStepSizeControl(PetscReal t, PetscReal dt, TIME_STEP_CONTROLLER *controller, MY_DATA *user);
PetscErrorCode DestroyTimeStepController(TIME_STEP_CONTROLLER *controller);


PetscErrorCode InitializeTimeStepController(PetscReal t, PetscReal dt, TIME_STEP_CONTROLLER *controller, MY_DATA *user)
{
    PetscErrorCode ierr;
    PetscInt       k;
    PetscInt       number_of_volume_species = NUMBER_OF_SPECIES_IN_VOLUME;
    PetscInt       number_of_surface_species = NUMBER_OF_SPECIES_ON_SURFACE;

    PetscFunctionBegin;

    controller->tmp_t  = t;
    controller->tmp_dt = dt;
    controller->new_t  = t;
    controller->new_dt = dt;

    for(k=0;k<number_of_volume_species;++k){
        ierr = getSingleSpeciesVec(user->c_volume_global, user->c_volume_single, user->volume_species_scatter[k]); CHKERRQ(ierr);
        ierr = VecNorm(user->c_volume_single,NORM_MAX,&controller->current_max_c_volume[k]); CHKERRQ(ierr);
        controller->old_max_c_volume[k] = controller->current_max_c_volume[k];
    }
    for(k=0;k<number_of_surface_species;++k){
        ierr = getSingleSpeciesVec(user->c_surface_global, user->c_surface_single, user->surface_species_scatter[k]); CHKERRQ(ierr);
        ierr = VecNorm(user->c_surface_single,NORM_MAX,&controller->current_max_c_surface[k]); CHKERRQ(ierr);
        controller->old_max_c_surface[k] = controller->current_max_c_surface[k];
    }

    controller->flg = -1;

    PetscFunctionReturn(0);

}

/*
 * increase or decrease dt according to maximum relative change of all species
 * no rejection of any steps
 *
 */
PetscErrorCode SingleStepSizeControl(PetscReal t, PetscReal dt, TIME_STEP_CONTROLLER *controller, MY_DATA *user)
{
    PetscErrorCode ierr;
    PetscInt       k;
    PetscInt       number_of_volume_species = NUMBER_OF_SPECIES_IN_VOLUME;
    PetscInt       number_of_surface_species = NUMBER_OF_SPECIES_ON_SURFACE;
    PetscReal      dt_increase_ratio, dt_decrease_ratio;
    PetscReal      c_relative_change_min, c_relative_change_max;
    PetscReal      c_change, c_change_max;
    PetscReal      dt_min, dt_max;

    PetscFunctionBegin;

    dt_min  = DT_MIN;
    dt_max  = DT_MAX;
    dt_increase_ratio = TIMESTEP_INCREASE_RATIO;
    dt_decrease_ratio = TIMESTEP_DECREASE_RATIO;
    c_relative_change_min = C_RELATIVE_TOL_MIN;
    c_relative_change_max = C_RELATIVE_TOL_MAX;

    controller->new_t  = t + dt;

    controller->new_dt = dt;

    c_change_max = 0.0;
    for(k=0;k<number_of_volume_species;++k){
        ierr = getSingleSpeciesVec(user->c_volume_global, user->c_volume_single, user->volume_species_scatter[k]); CHKERRQ(ierr);
        ierr = VecNorm(user->c_volume_single,NORM_MAX,&controller->current_max_c_volume[k]); CHKERRQ(ierr);
        c_change = abs(controller->current_max_c_volume[k] - controller->old_max_c_volume[k]) / controller->current_max_c_volume[k];
        if(c_change > c_change_max){
            c_change_max = c_change;
        }
        controller->old_max_c_volume[k] = controller->current_max_c_volume[k];      
    }
    for(k=0;k<number_of_surface_species;++k){
        ierr = getSingleSpeciesVec(user->c_surface_global, user->c_surface_single, user->surface_species_scatter[k]); CHKERRQ(ierr);
        ierr = VecNorm(user->c_surface_single,NORM_MAX,&controller->current_max_c_surface[k]); CHKERRQ(ierr);
        c_change = abs(controller->current_max_c_surface[k] - controller->old_max_c_surface[k]) / controller->current_max_c_surface[k];
        if(c_change > c_change_max){
            c_change_max = c_change;
        }
        controller->old_max_c_surface[k] = controller->current_max_c_surface[k];      
    }
    if(c_change_max<c_relative_change_min && dt<dt_max){
        controller->new_dt = dt * dt_increase_ratio;
        if(controller->new_dt>dt_max)
            controller->new_dt = dt_max;
    } else if(c_change_max>c_relative_change_max && dt>dt_min){
        controller->new_dt = dt * dt_decrease_ratio;
        if(controller->new_dt<dt_min)
            controller->new_dt = dt_min;
    }

    controller->tmp_t  = controller->new_t;
    controller->tmp_dt = controller->new_dt;

    controller->flg = 0;

    PetscFunctionReturn(0);
}

PetscErrorCode DestroyTimeStepController(TIME_STEP_CONTROLLER *controller)
{
    PetscFunctionBegin;

    PetscFunctionReturn(0);
}
#endif
