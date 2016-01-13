#ifndef _UPDATE_VELOCITY_FIELD_
#define _UPDATE_VELOCITY_FIELD_

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
 */
PetscErrorCode UpdateVelocityField(PetscReal t, MY_DATA *user);
PetscReal updateVx(PetscReal x, PetscReal y, PetscReal Lx, PetscReal Ly, PetscReal t);
PetscReal updateVy(PetscReal x, PetscReal y, PetscReal Lx, PetscReal Ly, PetscReal t);

/*
 * update the velocity field
 */
PetscErrorCode UpdateVelocityField(PetscReal t, MY_DATA *user)
{
  PetscInt       i,j;
  PetscInt       Nx_volume,Ny_volume,xs_volume,ys_volume,xm_volume,ym_volume;
  PetscErrorCode ierr;
  PetscReal      Lx, Ly;
  PetscReal      **array_vx_left, **array_vx_right, **array_vy_down, **array_vy_up;
  PetscReal      *array_x_vx, *array_y_vx, *array_x_vy, *array_y_vy;

  PetscFunctionBegin;

  Lx = X_LENGTH; Ly = Y_LENGTH;

  ierr = DMDAGetInfo(user->da_volume_dof_1, PETSC_IGNORE, &Nx_volume, &Ny_volume, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
                   PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE);CHKERRQ(ierr);

/*
 *      Get a pointer to vector data.
 *      - For default PETSc vectors, VecGetArray() returns a pointer to
 *        the data array.  Otherwise, the routine is implementation dependent.
 *      - You MUST call VecRestoreArray() when you no longer need access to
 *        the array.
 */
  ierr = DMDAVecGetArray(user->da_volume_dof_1, user->vx_left_global, &array_vx_left);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da_volume_dof_1, user->vx_right_global, &array_vx_right);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da_volume_dof_1, user->vy_down_global, &array_vy_down);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da_volume_dof_1, user->vy_up_global, &array_vy_up);CHKERRQ(ierr);

  ierr = VecGetArray(user->x_vx, &array_x_vx); CHKERRQ(ierr);
  ierr = VecGetArray(user->y_vx, &array_y_vx); CHKERRQ(ierr);
  ierr = VecGetArray(user->x_vy, &array_x_vy); CHKERRQ(ierr);
  ierr = VecGetArray(user->y_vy, &array_y_vy); CHKERRQ(ierr);

/*
 *      Get local grid boundaries (for 2-dimensional DMDA):
 *      xs, ys   - starting grid indices (no ghost points)
 *      xm, ym   - widths of local grid (no ghost points)
 *
 */
  ierr = DMDAGetCorners(user->da_volume_dof_1, &xs_volume, &ys_volume, PETSC_NULL, &xm_volume, &ym_volume, PETSC_NULL);CHKERRQ(ierr);

/*
 *      Update velocity field over the locally owned part of the grid
 */
  for (j=ys_volume; j<ys_volume+ym_volume; j++) {
    for (i=xs_volume; i<xs_volume+xm_volume; i++) {
      array_vx_left[j][i]  = updateVx(array_x_vx[i],array_y_vx[j],Lx,Ly,t);
      array_vx_right[j][i] = updateVx(array_x_vx[i+1],array_y_vx[j],Lx,Ly,t);
      array_vy_down[j][i]  = updateVy(array_x_vy[i],array_y_vy[j],Lx,Ly,t);
      array_vy_up[j][i]  = updateVy(array_x_vy[i],array_y_vy[j+1],Lx,Ly,t);
    }
  }

/*
 *      Restore vector
 */
  ierr = DMDAVecRestoreArray(user->da_volume_dof_1, user->vx_left_global, &array_vx_left);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da_volume_dof_1, user->vx_right_global, &array_vx_right);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da_volume_dof_1, user->vy_down_global, &array_vy_down);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da_volume_dof_1, user->vy_up_global, &array_vy_up);CHKERRQ(ierr);

  ierr = VecRestoreArray(user->x_vx, &array_x_vx); CHKERRQ(ierr);
  ierr = VecRestoreArray(user->y_vx, &array_y_vx); CHKERRQ(ierr);
  ierr = VecRestoreArray(user->x_vy, &array_x_vy); CHKERRQ(ierr);
  ierr = VecRestoreArray(user->y_vy, &array_y_vy); CHKERRQ(ierr);

  PetscFunctionReturn(0);
} 

/*
 * update vx based on (x,y) and total length (Lx,Ly) and time t
 */
PetscReal updateVx(PetscReal x, PetscReal y, PetscReal Lx, PetscReal Ly, PetscReal t)
{
    PetscReal Vx_max, Vx;

    PetscFunctionBegin;
    Vx_max = VX_MAX;
    Vx = 4.0*y*(Ly - y)*Vx_max/(Ly*Ly);
    PetscFunctionReturn(Vx);
}

/*
 * calculate vy based on (x,y) and total length (Lx,Ly) and time t
 */
PetscReal updateVy(PetscReal x, PetscReal y, PetscReal Lx, PetscReal Ly, PetscReal t)
{
    PetscFunctionBegin;
    PetscFunctionReturn(0.0);
}

#endif
