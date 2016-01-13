#ifndef _FORM_GRID_H_
#define _FORM_GRID_H_

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

PetscErrorCode FormGrid(MY_DATA *user);
PetscReal findGridRatio(PetscReal dy_bottom, PetscReal Ly, PetscInt Ny, PetscInt Ny_bottom);
PetscReal ratioGeometricSum(PetscReal r, PetscInt Ny);

/*
 * calculate the grid structure and label coordinate fields for both volume and surface
 */
PetscErrorCode FormGrid(MY_DATA *user)
{
  PetscInt       i,j;
  PetscInt       Nx_volume,Ny_volume,xs_volume,ys_volume,xm_volume,ym_volume, Ny_volume_bottom;
  PetscInt       Nx_surface,xs_surface,xm_surface;
  PetscErrorCode ierr;
  PetscReal      Lx, Ly;
  PetscReal      dx_volume, dy_volume_bottom, dy_volume_ratio, dy_temp_ratio;
  PetscReal      dx_surface;
  PetscReal      **array_x_volume, **array_y_volume, **array_dx_volume, **array_dy_volume;
  PetscReal      *array_x_surface, *array_dx_surface;


  PetscFunctionBegin;

  Ny_volume_bottom = NUMBER_OF_Y_GRIDS_BOTTOM;
  Lx = X_LENGTH; Ly = Y_LENGTH;
  dy_volume_bottom = BOTTOM_LAYER_DY;

  // set dy0 in MY_DATA
  user->dy0 = dy_volume_bottom;

  ierr = DMDAGetInfo(user->da_volume_dof_1, PETSC_IGNORE, &Nx_volume, &Ny_volume, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
                   PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE);CHKERRQ(ierr);
  ierr = DMDAGetInfo(user->da_surface_dof_1, PETSC_IGNORE, &Nx_surface, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
                   PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE);CHKERRQ(ierr);

  if(Nx_volume != Nx_surface){
      SETERRQ(PETSC_COMM_SELF,1,"FormGrid(): Volume grids and surface grids do not agree, make sure they agree.\n");
  }

/*
 *      Get a pointer to vector data.
 *      - For default PETSc vectors, VecGetArray() returns a pointer to
 *        the data array.  Otherwise, the routine is implementation dependent.
 *      - You MUST call VecRestoreArray() when you no longer need access to
 *        the array.
 */
 // get array for volume coordinate vectors
  ierr = DMDAVecGetArray(user->da_volume_dof_1, user->x_volume_global, &array_x_volume);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da_volume_dof_1, user->y_volume_global, &array_y_volume);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da_volume_dof_1, user->dx_volume_global, &array_dx_volume);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da_volume_dof_1, user->dy_volume_global, &array_dy_volume);CHKERRQ(ierr);

 // get array for surface coordinate vectors
  ierr = DMDAVecGetArray(user->da_surface_dof_1, user->x_surface_global, &array_x_surface);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da_surface_dof_1, user->dx_surface_global, &array_dx_surface);CHKERRQ(ierr);

/*
 *      Get local grid boundaries (for 2-dimensional DMDA):
 *      xs, ys   - starting grid indices (no ghost points)
 *      xm, ym   - widths of local grid (no ghost points)
 *
 */
  ierr = DMDAGetCorners(user->da_volume_dof_1, &xs_volume, &ys_volume, PETSC_NULL, &xm_volume, &ym_volume, PETSC_NULL);CHKERRQ(ierr);
  ierr = DMDAGetCorners(user->da_surface_dof_1, &xs_surface, PETSC_NULL, PETSC_NULL, &xm_surface, PETSC_NULL, PETSC_NULL);CHKERRQ(ierr);

/*
 *      Form coordinate system over the locally owned part of the grid
 */
  dx_volume = Lx/(PetscReal)Nx_volume;
  dx_surface = Lx/(PetscReal)Nx_surface;
  if( dy_volume_bottom > Ly/(PetscReal)Ny_volume ){
      SETERRQ(PETSC_COMM_SELF,1,"FormGrid(): dy0 > Ly/Ny, make sure dy0 < Ly/Ny.\n");
  }
  dy_volume_ratio = findGridRatio(dy_volume_bottom, Ly, Ny_volume, Ny_volume_bottom);
  if( dy_volume_ratio > 2.0 ){
      SETERRQ(PETSC_COMM_SELF,1,"FormGrid(): dy_ratio > 2, this might cause instability in diffusion solver.\n");
  }
  // set dy_ratio in MY_DATA
  user->dy_ratio = dy_volume_ratio;
  
  for (j=ys_volume; j<ys_volume+ym_volume; j++) {
    for (i=xs_volume; i<xs_volume+xm_volume; i++) {
      array_x_volume[j][i]  = dx_volume * ((PetscReal)i + 0.5);
      array_dx_volume[j][i] = dx_volume;

      if (j<Ny_volume_bottom) {
         array_y_volume[j][i]  = dy_volume_bottom * ((PetscReal)j + 0.5);
         array_dy_volume[j][i] = dy_volume_bottom;
      } else {
         dy_temp_ratio = pow(dy_volume_ratio,j-Ny_volume_bottom);
         array_y_volume[j][i]  = dy_volume_bottom * Ny_volume_bottom 
                                              + dy_volume_bottom * dy_volume_ratio * (dy_temp_ratio - 1.0) / (dy_volume_ratio - 1.0)
                                               + dy_volume_bottom * dy_temp_ratio * dy_volume_ratio * 0.5;
         array_dy_volume[j][i] = dy_temp_ratio * dy_volume_ratio * dy_volume_bottom;
      }
    }
  }

  for (i=xs_surface; i<xs_surface+xm_surface; i++) {
      array_x_surface[i]  = dx_surface * ((PetscReal)i + 0.5);
      array_dx_surface[i] = dx_surface;
  }

/*
 *      Restore vector
 */
  ierr = DMDAVecRestoreArray(user->da_volume_dof_1, user->x_volume_global, &array_x_volume);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da_volume_dof_1, user->y_volume_global, &array_y_volume);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da_volume_dof_1, user->dx_volume_global, &array_dx_volume);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da_volume_dof_1, user->dy_volume_global, &array_dy_volume);CHKERRQ(ierr);

  ierr = DMDAVecRestoreArray(user->da_surface_dof_1, user->x_surface_global, &array_x_surface);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da_surface_dof_1, user->dx_surface_global, &array_dx_surface);CHKERRQ(ierr);

  PetscFunctionReturn(0);
} 

/*
 *  _                                                   
 *   |                               
 *   |  dy_bottom * r^(Ny-Ny_bottom)
 *   |  ...                 
 *  _| 
 *   |  dy_bottom * r
 *  _|                
 *  _|  dy_bottom   
 *  _|  dy_bottom
 *
 * find the ratio r that satisfy equation
 *    dy_bottom * (r + r^2 + ... + r^(Ny-Ny_bottom))
 *       = Ly - dy_bottom * Ny_bottom
 *
 */
PetscReal findGridRatio(PetscReal dy_bottom, PetscReal Ly, PetscInt Ny, PetscInt Ny_bottom)
{
    PetscReal rMin, rMax, r;
    PetscReal Ly_up;
    PetscInt  Ny_up;
    PetscReal ratioSum;
  
    PetscFunctionBegin; 

    rMin = 0.0; rMax = 2.0; r = 1.0;
    Ly_up = Ly - dy_bottom * (PetscReal)Ny_bottom;
    Ny_up = Ny-Ny_bottom;
    ratioSum = Ly_up/dy_bottom;
     
    while(ratioGeometricSum(rMax, Ny_up) < ratioSum){
        rMax *= 2.0;
    }

    // now rMin and rMax bound exact r
    r = 0.5*(rMin + rMax);
    while(rMax-rMin > CALCULATION_ERROR_TOLERANCE){
        if(ratioGeometricSum(r, Ny_up) > ratioSum){
            rMax = r;
        } else{
            rMin = r;
        }
        r = 0.5*(rMin + rMax);
    }

    PetscFunctionReturn(r);
}

/* 
 * calculate r*(r^Ny-1)/(r-1) = r+r^2+...+r^Ny
 */
PetscReal ratioGeometricSum(PetscReal r, PetscInt Ny)
{
    PetscReal ratioSum;

    PetscFunctionBegin;

    if(r == 1.0){
        ratioSum = (PetscReal)Ny;
    } else {
        ratioSum = r*(pow(r,Ny)-1.0)/(r-1.0);
    }

    PetscFunctionReturn(ratioSum);
}

#endif
