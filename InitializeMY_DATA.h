#ifndef _INITIALIZE_MY_DATA_H_
#define _INITIALIZE_MY_DATA_H_

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

PetscErrorCode InitializeMY_DATA(MY_DATA *user);

/*
 * initialize MY_DATA
 * 1. set comm[0] and comm[1]
 * 2. form the DMs in MY_DATA
 * 3. form the Vecs in MY_DATA
 * form the DMs and Vectors in MY_DATA for both volume and surface
 */
PetscErrorCode InitializeMY_DATA(MY_DATA *user)
{
  PetscErrorCode ierr;
  PetscInt      number_of_volume_species, number_of_surface_species;
  PetscInt      number_of_x_grids, number_of_y_grids, number_of_volume_grids, number_of_surface_grids;
  PetscInt      k;

  PetscFunctionBegin;

  number_of_volume_species  = NUMBER_OF_SPECIES_IN_VOLUME;
  number_of_surface_species = NUMBER_OF_SPECIES_ON_SURFACE;
  number_of_x_grids         = NUMBER_OF_X_GRIDS;
  number_of_y_grids         = NUMBER_OF_Y_GRIDS;
  number_of_volume_grids    = number_of_x_grids * number_of_y_grids;
  number_of_surface_grids   = number_of_x_grids;

  // set both communicator to PETSC_COMM_WORLD
  user->comm[0] = PETSC_COMM_WORLD;
  user->comm[1] = PETSC_COMM_WORLD;

// create da_c_volume, DM for species in volume
// only parallel in the x direction to miminize communication to the boundary layer
#ifdef LEFT_RIGHT_PERIODIC
  ierr = DMDACreate2d(user->comm[0], DMDA_BOUNDARY_PERIODIC, DMDA_BOUNDARY_GHOSTED,DMDA_STENCIL_STAR,-number_of_x_grids,-number_of_y_grids,PETSC_DECIDE,1,number_of_volume_species,2,PETSC_NULL,PETSC_NULL,&user->da_c_volume);CHKERRQ(ierr);
#else
  ierr = DMDACreate2d(user->comm[0], DMDA_BOUNDARY_GHOSTED, DMDA_BOUNDARY_GHOSTED,DMDA_STENCIL_STAR,-number_of_x_grids,-number_of_y_grids,PETSC_DECIDE,1,number_of_volume_species,2,PETSC_NULL,PETSC_NULL,&user->da_c_volume);CHKERRQ(ierr);
#endif
// create da_volume_dof_1, DM for variables in volume whose dof=1
#ifdef LEFT_RIGHT_PERIODIC
  ierr = DMDACreate2d(user->comm[0], DMDA_BOUNDARY_PERIODIC, DMDA_BOUNDARY_GHOSTED,DMDA_STENCIL_STAR,-number_of_x_grids,-number_of_y_grids,PETSC_DECIDE,1,1,2,PETSC_NULL,PETSC_NULL,&user->da_volume_dof_1);CHKERRQ(ierr);
#else
  ierr = DMDACreate2d(user->comm[0], DMDA_BOUNDARY_GHOSTED, DMDA_BOUNDARY_GHOSTED,DMDA_STENCIL_STAR,-number_of_x_grids,-number_of_y_grids,PETSC_DECIDE,1,1,2,PETSC_NULL,PETSC_NULL,&user->da_volume_dof_1);CHKERRQ(ierr);
#endif

// create da_c_surface, DM for species on surface
  ierr = DMDACreate1d(user->comm[1], DMDA_BOUNDARY_NONE, -number_of_x_grids, number_of_surface_species,1,PETSC_NULL,&user->da_c_surface);CHKERRQ(ierr);
// create da_surface_dof_1, DM for variables in surface whose dof=1
  ierr = DMDACreate1d(user->comm[1], DMDA_BOUNDARY_NONE, -number_of_x_grids, 1,1,PETSC_NULL,&user->da_surface_dof_1);CHKERRQ(ierr);


// create c_volume_global, global vector of species concentrations in volume
// and c_volume_local, local vector of species concentrations in volume
// and their corresponding residual vectors
  ierr = DMCreateGlobalVector(user->da_c_volume,&user->c_volume_global);CHKERRQ(ierr);
  ierr = DMCreateLocalVector(user->da_c_volume,&user->c_volume_local);CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(user->da_c_volume,&user->r_volume_global);CHKERRQ(ierr);
  ierr = DMCreateLocalVector(user->da_c_volume,&user->r_volume_local);CHKERRQ(ierr);
// create x_volume_global, y_volume_global, dx_volume_global, dy_volume_global, global vectors of coordinates in volume
  ierr = DMCreateGlobalVector(user->da_volume_dof_1,&user->x_volume_global);CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(user->da_volume_dof_1,&user->y_volume_global);CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(user->da_volume_dof_1,&user->dx_volume_global);CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(user->da_volume_dof_1,&user->dy_volume_global);CHKERRQ(ierr);
// create vx_left_global, vx_right_global, vy_down_global, vy_up_global, global vectors of velocities in volume
  ierr = DMCreateGlobalVector(user->da_volume_dof_1,&user->vx_left_global);CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(user->da_volume_dof_1,&user->vx_right_global);CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(user->da_volume_dof_1,&user->vy_down_global);CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(user->da_volume_dof_1,&user->vy_up_global);CHKERRQ(ierr);
// create c_volume_single, global vector of a single species in volume, mainly for output, template, and error/change controls
  ierr = DMCreateGlobalVector(user->da_volume_dof_1,&user->c_volume_single);CHKERRQ(ierr);

// create a list of scatter from c_volume_global to c_volume_single, one for each species
  for(k=0;k<number_of_volume_species;k++){
     ierr = setSingleSpeciesScatter(user->c_volume_global, user->c_volume_single, number_of_volume_species, number_of_volume_grids, k, &user->volume_species_scatter[k]); CHKERRQ(ierr);
  }

// create c_surface_global, global vector of species concentrations on surface
// and c_surface_local, local vector of species concentrations on surface
// and their corresponding residual vectors
  ierr = DMCreateGlobalVector(user->da_c_surface,&user->c_surface_global);CHKERRQ(ierr);
  ierr = DMCreateLocalVector(user->da_c_surface,&user->c_surface_local);CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(user->da_c_surface,&user->r_surface_global);CHKERRQ(ierr);
  ierr = DMCreateLocalVector(user->da_c_surface,&user->r_surface_local);CHKERRQ(ierr);
// create x_surface_global, dx_surface_global, global vectors of coordinates in surface
  ierr = DMCreateGlobalVector(user->da_surface_dof_1,&user->x_surface_global);CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(user->da_surface_dof_1,&user->dx_surface_global);CHKERRQ(ierr);
// create c_suface_single, global vector of a single species on surface, mainly for output, template, and error/change controls
  ierr = DMCreateGlobalVector(user->da_surface_dof_1,&user->c_surface_single);CHKERRQ(ierr);

// create a list of scatter from c_volume_global to c_volume_single, one for each species
  for(k=0;k<number_of_surface_species;k++){
     ierr = setSingleSpeciesScatter(user->c_surface_global, user->c_surface_single, number_of_surface_species, number_of_surface_grids, k, &user->surface_species_scatter[k]); CHKERRQ(ierr);
  }

// create c_surface_global, global vector of species concentrations on surface
  user->D_volume = VOLUME_DIFFUSION_CONSTANT;


// create x_vx, y_vx, x_vy, and y_vy for interface with flow solver
/*
  ierr = VecCreate(user->comm[0],&user->x_vx); CHKERRQ(ierr);
  ierr = VecCreate(user->comm[0],&user->y_vx); CHKERRQ(ierr);
  ierr = VecCreate(user->comm[0],&user->x_vy); CHKERRQ(ierr);
  ierr = VecCreate(user->comm[0],&user->y_vy); CHKERRQ(ierr);
*/
/*
  ierr = VecCreate(PETSC_COMM_SELF,&user->x_vx); CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_SELF,&user->y_vx); CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_SELF,&user->x_vy); CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_SELF,&user->y_vy); CHKERRQ(ierr);

// set the size of x_vx, y_vx, x_vy, and y_vy for interface with flow solver
  ierr = VecSetSizes(user->x_vx,PETSC_DECIDE,number_of_x_grids+1); CHKERRQ(ierr); 
  ierr = VecSetSizes(user->y_vx,PETSC_DECIDE,number_of_y_grids);   CHKERRQ(ierr); 
  ierr = VecSetSizes(user->x_vy,PETSC_DECIDE,number_of_x_grids);   CHKERRQ(ierr); 
  ierr = VecSetSizes(user->y_vy,PETSC_DECIDE,number_of_y_grids+1); CHKERRQ(ierr); 

  ierr = VecSetFromOptions(user->x_vx); CHKERRQ(ierr);
  ierr = VecSetFromOptions(user->y_vx); CHKERRQ(ierr);
  ierr = VecSetFromOptions(user->x_vy); CHKERRQ(ierr);
  ierr = VecSetFromOptions(user->y_vy); CHKERRQ(ierr);
*/
// register CLASSID's
#ifdef PETSC_USE_LOG
  ierr = PetscClassIdRegister("Surface reaction solver",&user->SURFACE_REACTION_CLASSID); CHKERRQ(ierr);
  ierr = PetscClassIdRegister("Volume reaction solver",&user->VOLUME_REACTION_CLASSID); CHKERRQ(ierr);
  ierr = PetscClassIdRegister("Volume convection solver",&user->VOLUME_CONVECTION_CLASSID); CHKERRQ(ierr);
  ierr = PetscClassIdRegister("Volume diffusion solver",&user->VOLUME_DIFFUSION_CLASSID); CHKERRQ(ierr);
  ierr = PetscClassIdRegister("Step control",&user->STEP_CONTROL_CLASSID); CHKERRQ(ierr);
#endif


  user->initialized = true;

  PetscFunctionReturn(0);
} 

#endif
