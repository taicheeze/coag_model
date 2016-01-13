#ifndef _DESTROY_MY_DATA_H_
#define _DESTROY_MY_DATA_H_

#include "DataStructure.h"
#include "SingleSpeciesScatter.h"

// velocity is partially staggered with species, therefore 
// for vx: number of x grids is NUMBER_OF_X_GRIDS+1, number of y grids is NUMBER_OF_Y_GRIDS
// for vy: number of x grids is NUMBER_OF_X_GRIDS, number of y grids is NUMBER_OF_Y_GRIDS+1
//
//   ____._________.____
//  |         |         |
//  *    x    *    x    * dy_(j+1)
//  |____.____|____.____|
//  |         |         |
//  *    x    *    x    * dy_(j)
//  |____.____|____.____|
//      dx_i    dx_(i+1)
//    
//  x : concentration's of species  
//  * : vx's
//  . : vy's
//
//
// For one cell, we store vx_left, vx_right, vy_up, vy_down
//            vy_up
//          ____.____
//         |         |           |
// vx_left *    c    * vx_right  dy
//         |____.____|           |
//           vy_down
//          ____dx____
//
// Thus, we have relationship
//   vx_left(i,j) = vx_right(i-1,j)
//   vy_down(i,j) = vy_up(i,j-1)
//

PetscErrorCode DestroyMY_DATA(MY_DATA *user);

/*
 * destroy MY_DATA
 * 1. destroy the Vecs in MY_DATA
 * 2. destroy the DMs in MY_DATA
 *
 */
PetscErrorCode DestroyMY_DATA(MY_DATA *user)
{
  PetscErrorCode ierr;
  PetscInt       k;
  PetscInt       number_of_volume_species, number_of_surface_species;

  PetscFunctionBegin;

  number_of_volume_species = NUMBER_OF_SPECIES_IN_VOLUME;
  number_of_surface_species = NUMBER_OF_SPECIES_ON_SURFACE;

// destroy c_volume_global, global vector of species concentrations in volume
// and c_volume_local, local vector of species concentrations in volume
// and their corresponding residual vectors
  ierr = VecDestroy(&user->c_volume_global);CHKERRQ(ierr);
  ierr = VecDestroy(&user->c_volume_local);CHKERRQ(ierr);
  ierr = VecDestroy(&user->r_volume_global);CHKERRQ(ierr);
  ierr = VecDestroy(&user->r_volume_local);CHKERRQ(ierr);
// destroy x_volume_global, y_volume_global, dx_volume_global, dy_volume_global, global vectors of coordinates in volume
  ierr = VecDestroy(&user->x_volume_global);CHKERRQ(ierr);
  ierr = VecDestroy(&user->y_volume_global);CHKERRQ(ierr);
  ierr = VecDestroy(&user->dx_volume_global);CHKERRQ(ierr);
  ierr = VecDestroy(&user->dy_volume_global);CHKERRQ(ierr);
// destroy vx_left_global, vx_right_global, vy_down_global, vy_up_global, global vectors of velocities in volume
  ierr = VecDestroy(&user->vx_left_global);CHKERRQ(ierr);
  ierr = VecDestroy(&user->vx_right_global);CHKERRQ(ierr);
  ierr = VecDestroy(&user->vy_down_global);CHKERRQ(ierr);
  ierr = VecDestroy(&user->vy_up_global);CHKERRQ(ierr);
// destroy c_volume_single
  ierr = VecDestroy(&user->c_volume_single);CHKERRQ(ierr);

// destory volume species scatters
  for(k=0;k<number_of_volume_species;k++){
     ierr = destroySingleSpeciesScatter(&(user->volume_species_scatter[k])); CHKERRQ(ierr);
  }

// destroy c_surface_global, global vector of species concentrations on surface
// destroy c_surface_local, local vector of species concentrations on surface
// and their corresponding residual vectors
  ierr = VecDestroy(&user->c_surface_global);CHKERRQ(ierr);
  ierr = VecDestroy(&user->c_surface_local);CHKERRQ(ierr);
  ierr = VecDestroy(&user->r_surface_global);CHKERRQ(ierr);
  ierr = VecDestroy(&user->r_surface_local);CHKERRQ(ierr);
// destroy x_surface_global, dx_surface_global, global vectors of coordinates in surface
  ierr = VecDestroy(&user->x_surface_global);CHKERRQ(ierr);
  ierr = VecDestroy(&user->dx_surface_global);CHKERRQ(ierr);
// destroy c_surface_single
  ierr = VecDestroy(&user->c_surface_single);CHKERRQ(ierr);

// destroy x_vx, y_vx, x_vy, y_vy
  ierr = VecDestroy(&user->x_vx);CHKERRQ(ierr);
  ierr = VecDestroy(&user->y_vx);CHKERRQ(ierr);
  ierr = VecDestroy(&user->x_vy);CHKERRQ(ierr);
  ierr = VecDestroy(&user->y_vy);CHKERRQ(ierr);

// destroy da_c_volume, DM for species in volume
  ierr = DMDestroy(&user->da_c_volume);CHKERRQ(ierr);
// destroy da_volume_dof_1, DM for variables in volume whose dof=1
  ierr = DMDestroy(&user->da_volume_dof_1);CHKERRQ(ierr);
// destroy da_c_surface, DM for species in surface
  ierr = DMDestroy(&user->da_c_surface);CHKERRQ(ierr);
// destroy da_surface_dof_1, DM for variables in surface whose dof=1
  ierr = DMDestroy(&user->da_surface_dof_1);CHKERRQ(ierr);

// destory surface species scatters
  for(k=0;k<number_of_surface_species;k++){
     ierr = destroySingleSpeciesScatter(&(user->surface_species_scatter[k])); CHKERRQ(ierr);
  }

  user->initialized = false;

  PetscFunctionReturn(0);
} 

#endif
