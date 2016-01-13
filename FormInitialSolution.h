#ifndef _FORM_INITIAL_SOLUTION_H_
#define _FORM_INITIAL_SOLUTION_H_

#include "DataStructure.h"
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

/*
 * both FormVelocityField and FormInitialSolution has to happen after FormGrid
 */
PetscErrorCode FormInitialSolution(MY_DATA *user);
// x: grid center, dx: grid length, Lx: whole x length
// patch_left: left boundary of patch, patch_right: right boundary of patch
PetscReal findRatioOfGridInPatch(PetscReal x, PetscReal dx, PetscReal Lx, PetscReal patch_left, PetscReal patch_right);

/*
 * form the initial species concentration for both volume and surface
 */
PetscErrorCode FormInitialSolution(MY_DATA *user)
{
  PetscInt       i,j,k, patch_index;
  PetscInt       number_volume_species, number_surface_species;
  PetscInt       Nx_volume,Ny_volume,xs_volume,ys_volume,xm_volume,ym_volume;
  PetscInt       xs_ghost_volume,ys_ghost_volume,xm_ghost_volume,ym_ghost_volume;
  PetscInt       Nx_surface,xs_surface,xm_surface, Nx_surface_compare;
  PetscErrorCode ierr;
  PetscReal      Lx, Ly;
  PetscReal      x, dx, ratio_of_grid_in_patch[NUMBER_OF_PATCHES];
  Species_Volume    **array_c_volume;
  Species_Surface   *array_c_surface;
  PetscReal      *array_x_surface, *array_dx_surface;
//  PetscMPIInt    rank, size;


  PetscFunctionBegin;

  number_volume_species = NUMBER_OF_SPECIES_IN_VOLUME; number_surface_species = NUMBER_OF_SPECIES_ON_SURFACE;
  Lx = X_LENGTH; Ly = Y_LENGTH;

  ierr = DMDAGetInfo(user->da_c_volume, PETSC_IGNORE, &Nx_volume, &Ny_volume, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
                   PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE);CHKERRQ(ierr);
  ierr = DMDAGetInfo(user->da_c_surface, PETSC_IGNORE, &Nx_surface, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
                   PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE);CHKERRQ(ierr);

  ierr = DMDAGetInfo(user->da_surface_dof_1, PETSC_IGNORE, &Nx_surface_compare, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
                   PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE);CHKERRQ(ierr);

  if(Nx_volume != Nx_surface || Nx_volume != Nx_surface_compare){
      SETERRQ(PETSC_COMM_SELF,1,"FormInitialSolution(): Volume grids and surface grids do not agree, or surface concentration grids and coordinate grids do no agree, make sure they agree.\n");
  }

/*
 *  First set all the concentrations to 0.0
 */
  VecSet(user->c_volume_global,0.0);
  VecSet(user->c_surface_global,0.0);

/*
    Scatter ghost points to local vector, using the 2-step process
        DMGlobalToLocalBegin(), DMGlobalToLocalEnd().
    By placing code between these two statements, computations can be
    done while messages are in transition.
*/
    ierr = DMGlobalToLocalBegin(user->da_c_volume,user->c_volume_global,INSERT_VALUES,user->c_volume_local);CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(user->da_c_volume,user->c_volume_global,INSERT_VALUES,user->c_volume_local);CHKERRQ(ierr);

/*
 *      Get a pointer to vector data.
 *      - For default PETSc vectors, VecGetArray() returns a pointer to
 *        the data array.  Otherwise, the routine is implementation dependent.
 *      - You MUST call VecRestoreArray() when you no longer need access to
 *        the array.
 */
  ierr = DMDAVecGetArray(user->da_c_volume, user->c_volume_local, &array_c_volume);CHKERRQ(ierr);

  ierr = DMDAVecGetArray(user->da_c_surface, user->c_surface_global, &array_c_surface);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da_surface_dof_1, user->x_surface_global, &array_x_surface);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da_surface_dof_1, user->dx_surface_global, &array_dx_surface);CHKERRQ(ierr);

/*
 *      Get local grid boundaries (for 2-dimensional DMDA):
 *      xs, ys   - starting grid indices (no ghost points)
 *      xm, ym   - widths of local grid (no ghost points)
 *
 */
  ierr = DMDAGetCorners(user->da_c_volume, &xs_volume, &ys_volume, PETSC_NULL, &xm_volume, &ym_volume, PETSC_NULL);CHKERRQ(ierr);
  ierr = DMDAGetCorners(user->da_c_surface, &xs_surface, PETSC_NULL, PETSC_NULL, &xm_surface, PETSC_NULL, PETSC_NULL);CHKERRQ(ierr);

/*
 *      Get local grid boundaries (for 2-dimensional DMDA):
 *      xs_ghost, ys_ghost   - starting grid indices (with ghost points)
 *      xm_ghost, ym_ghost   - widths of local grid (with ghost points)
 *
 */
  ierr = DMDAGetGhostCorners(user->da_c_volume, &xs_ghost_volume, &ys_ghost_volume, PETSC_NULL, &xm_ghost_volume, &ym_ghost_volume, PETSC_NULL);CHKERRQ(ierr);

/*
 *      Form the initial concentration of species over the locally owned part of the grid
 *      the concentrations at everywhere else in the volume are 0
 *      (Convection bring species in afterwards)
 */
/*
 * originally we set the leftmost grids have initialconcentrations
 * later decide to set it to 0
  MPI_Comm_rank(user->comm[0],&rank);
  if(rank==0){
    for (j=ys_volume; j<ys_volume+ym_volume; j++) {
      for (k=0; k<number_volume_species; k++) {
         array_c_volume[j][0].c[k] = VOLUME_INITIAL_CONCENTRATIONS[k];
      }    
    }
  }
*/

  // fill the entire channel with initial concentrations
  for (i=xs_volume; i<xs_volume+xm_volume; i++) {
    for (j=ys_volume; j<ys_volume+ym_volume; j++) {
      for (k=0; k<number_volume_species; k++) {
         array_c_volume[j][i].c[k] = VOLUME_INITIAL_CONCENTRATIONS[k];
      }    
    }
  }

  // fill the inlet ghost cells with initial concentrations
  // fill the outlet ghost cells with OUTLET_GHOST_CELL_VALUE
  // fill the wall ghost cells with WALL_GHOST_CELL_VALUE
  for (i=xs_ghost_volume; i<xs_ghost_volume+xm_ghost_volume; i++) {
      // upper boundary
      if( ys_ghost_volume + ym_ghost_volume > Ny_volume ) {
          // first set every ghost cell value to wall
          for (k=0; k<number_volume_species; k++) {
             array_c_volume[Ny_volume][i].c[k] = WALL_GHOST_CELL_VALUE;
          }
          // inlet is on upper boundary
          if( INLET_Y_BOUNDARIES[0] == Ny_volume-1 && INLET_Y_BOUNDARIES[1] == Ny_volume-1 ){
              if( i >= INLET_X_BOUNDARIES[0] && i <= INLET_X_BOUNDARIES[1] ){
                 for (k=0; k<number_volume_species; k++) {
                     array_c_volume[Ny_volume][i].c[k] = VOLUME_INITIAL_CONCENTRATIONS[k];
                 }    
              }
          }
          // outlet is on upper boundary
          if( OUTLET_Y_BOUNDARIES[0] == Ny_volume-1 && OUTLET_Y_BOUNDARIES[1] == Ny_volume-1 ){
              if( i >= OUTLET_X_BOUNDARIES[0] && i <= OUTLET_X_BOUNDARIES[1] ){
                 for (k=0; k<number_volume_species; k++) {
                     array_c_volume[Ny_volume][i].c[k] = OUTLET_GHOST_CELL_VALUE;
                 }    
              }
          }
      }
      
      // lower boundary should all be walls
      if( ys_ghost_volume < 0 ) {
          // set every ghost cell value to wall
          for (k=0; k<number_volume_species; k++) {
             array_c_volume[-1][i].c[k] = WALL_GHOST_CELL_VALUE;
          }
      }
  }

#ifndef LEFT_RIGHT_PERIODIC
  for (j=ys_ghost_volume; j<ys_ghost_volume+ym_ghost_volume; j++) {
      // left boundary
      if (xs_ghost_volume < 0) {
          // first set every ghost cell value to wall
          for (k=0; k<number_volume_species; k++) {
             array_c_volume[j][-1].c[k] = WALL_GHOST_CELL_VALUE;
          }
          // inlet is on left boundary
          if( INLET_X_BOUNDARIES[0] == 0 && INLET_X_BOUNDARIES[1] == 0 ){
              if( j >= INLET_Y_BOUNDARIES[0] && j <= INLET_Y_BOUNDARIES[1] ){
                 for (k=0; k<number_volume_species; k++) {
                     array_c_volume[j][-1].c[k] = VOLUME_INITIAL_CONCENTRATIONS[k];
                 }    
              }
          }
          // outlet is on left boundary
          if( OUTLET_X_BOUNDARIES[0] == 0 && OUTLET_X_BOUNDARIES[1] == 0 ){
              if( j >= OUTLET_Y_BOUNDARIES[0] && j <= OUTLET_Y_BOUNDARIES[1] ){
                 for (k=0; k<number_volume_species; k++) {
                     array_c_volume[j][-1].c[k] = OUTLET_GHOST_CELL_VALUE;
                 }     
              }
          }  
      }
      
      // right boundary
      if( xs_ghost_volume + xm_ghost_volume > Nx_volume ) {
          // first set every ghost cell value to wall
          for (k=0; k<number_volume_species; k++) {
             array_c_volume[j][Nx_volume].c[k] = WALL_GHOST_CELL_VALUE;
          }
          // inlet is on left boundary
          if( INLET_X_BOUNDARIES[0] == Nx_volume-1 && INLET_X_BOUNDARIES[1] == Nx_volume-1 ){
              if( j >= INLET_Y_BOUNDARIES[0] && j <= INLET_Y_BOUNDARIES[1] ){
                 for (k=0; k<number_volume_species; k++) {
                     array_c_volume[j][Nx_volume].c[k] = VOLUME_INITIAL_CONCENTRATIONS[k];
                 }    
              }
          }
          // outlet is on left boundary
          if( OUTLET_X_BOUNDARIES[0] == Nx_volume-1 && OUTLET_X_BOUNDARIES[1] == Nx_volume-1 ){
              if( j >= OUTLET_Y_BOUNDARIES[0] && j <= OUTLET_Y_BOUNDARIES[1] ){
                 for (k=0; k<number_volume_species; k++) {
                     array_c_volume[j][Nx_volume].c[k] = OUTLET_GHOST_CELL_VALUE;
                 }    
              }
          }
      }
      
  }
#endif

  for (i=xs_surface; i<xs_surface+xm_surface; i++) {
    x =  array_x_surface[i];
    dx =  array_dx_surface[i];
    for (k=0; k<number_surface_species; k++) {
       array_c_surface[i].c[k] = 0.0;
    }    
    for ( patch_index=0; patch_index<NUMBER_OF_PATCHES; ++patch_index){
       ratio_of_grid_in_patch[patch_index] = findRatioOfGridInPatch(x,dx,Lx,PATCH_LEFT_BOUNDARIES[patch_index],PATCH_RIGHT_BOUNDARIES[patch_index]);
       if (fabs(ratio_of_grid_in_patch[patch_index] - 0.0) > 1.0e-6
             && fabs(ratio_of_grid_in_patch[patch_index] - 1.0) > 1.0e-6) {
           ierr = PetscPrintf(PETSC_COMM_WORLD,"WARNING: The boundaries of patch #%d does not coincide with one of the grid's boundary, be aware of unwanted surface diffusion in the grids containing the boundaries of patch #%d\n",patch_index+1,patch_index+1);CHKERRQ(ierr);
       }
       for (k=0; k<number_surface_species; k++) {
          array_c_surface[i].c[k] += ratio_of_grid_in_patch[patch_index] * SURFACE_INITIAL_CONCENTRATIONS_PATCHES[patch_index][k] ;
       }    
    }
  }

/*
 *      Restore vector
 */
  ierr = DMDAVecRestoreArray(user->da_c_volume, user->c_volume_local, &array_c_volume);CHKERRQ(ierr);

  ierr = DMDAVecRestoreArray(user->da_c_surface, user->c_surface_global, &array_c_surface);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da_surface_dof_1, user->x_surface_global, &array_x_surface);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da_surface_dof_1, user->dx_surface_global, &array_dx_surface);CHKERRQ(ierr);

/*
     Insert values from the local OUTPUT vector into the global 
     output vector
*/
    ierr = DMLocalToGlobalBegin(user->da_c_volume,user->c_volume_local,INSERT_VALUES,user->c_volume_global);CHKERRQ(ierr);
    ierr = DMLocalToGlobalEnd(user->da_c_volume,user->c_volume_local,INSERT_VALUES,user->c_volume_global);CHKERRQ(ierr);

  PetscFunctionReturn(0);
} 

// x: grid center, dx: grid length, Lx: whole x length
// patch_left: left boundary of patch, patch_right: right boundary of patch
PetscReal findRatioOfGridInPatch(PetscReal x, PetscReal dx, PetscReal Lx, PetscReal patch_left, PetscReal patch_right)
{
  PetscReal ratio, grid_left, grid_right, over_left, over_right;

  PetscFunctionBegin;

  grid_left  = x - dx/2.0;
  grid_right = x + dx/2.0;

  if(grid_left < patch_right && grid_right > patch_left){
    // calculate the boudaries of overlapping length
     over_left = grid_left>patch_left?grid_left:patch_left;
     over_right = grid_right<patch_right?grid_right:patch_right;

     ratio = (over_right - over_left)/dx;
  } else {
     ratio = 0.0;
  }

  PetscFunctionReturn(ratio);
}
#endif
