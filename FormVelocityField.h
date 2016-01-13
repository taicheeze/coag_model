#ifndef _FORM_VELOCITY_FIELD_
#define _FORM_VELOCITY_FIELD_

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
PetscErrorCode FormVelocityField(MY_DATA *user);
PetscReal findVx(PetscReal x, PetscReal y, PetscReal Lx, PetscReal Ly);
PetscReal findVy(PetscReal x, PetscReal y, PetscReal Lx, PetscReal Ly);

/*
 * set the velocity field
 */
PetscErrorCode FormVelocityField(MY_DATA *user)
{
  PetscInt       i,j;
  PetscInt       Nx_volume,Ny_volume,xs_volume,ys_volume,xm_volume,ym_volume;
  PetscReal      x,y,dx,dy;
  PetscErrorCode ierr;
  PetscReal      Lx, Ly;
  PetscReal      **array_vx_left, **array_vx_right, **array_vy_down, **array_vy_up;
  PetscReal      **array_x_volume, **array_y_volume, **array_dx_volume, **array_dy_volume;
  Vec            temp_x_vx, temp_y_vx, temp_x_vy, temp_y_vy;
  VecScatter     temp_ctx;

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

  ierr = DMDAVecGetArray(user->da_volume_dof_1, user->x_volume_global, &array_x_volume);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da_volume_dof_1, user->y_volume_global, &array_y_volume);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da_volume_dof_1, user->dx_volume_global, &array_dx_volume);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da_volume_dof_1, user->dy_volume_global, &array_dy_volume);CHKERRQ(ierr);

/*
 *      Get local grid boundaries (for 2-dimensional DMDA):
 *      xs, ys   - starting grid indices (no ghost points)
 *      xm, ym   - widths of local grid (no ghost points)
 *
 */
  ierr = DMDAGetCorners(user->da_volume_dof_1, &xs_volume, &ys_volume, PETSC_NULL, &xm_volume, &ym_volume, PETSC_NULL);CHKERRQ(ierr);

// create temp_x_vx, temp_y_vx, temp_x_vy, and temp_y_vy for creating x_vx, y_vx, x_vy, y_vy
  ierr = VecCreate(user->comm[0],&temp_x_vx); CHKERRQ(ierr);
  ierr = VecCreate(user->comm[0],&temp_y_vx); CHKERRQ(ierr);
  ierr = VecCreate(user->comm[0],&temp_x_vy); CHKERRQ(ierr);
  ierr = VecCreate(user->comm[0],&temp_y_vy); CHKERRQ(ierr);

// set the size of x_vx, y_vx, x_vy, and y_vy for interface with flow solver
  ierr = VecSetSizes(temp_x_vx,PETSC_DECIDE,Nx_volume+1); CHKERRQ(ierr); 
  ierr = VecSetSizes(temp_y_vx,PETSC_DECIDE,Ny_volume);   CHKERRQ(ierr); 
  ierr = VecSetSizes(temp_x_vy,PETSC_DECIDE,Nx_volume);   CHKERRQ(ierr); 
  ierr = VecSetSizes(temp_y_vy,PETSC_DECIDE,Ny_volume+1); CHKERRQ(ierr); 

  ierr = VecSetFromOptions(temp_x_vx); CHKERRQ(ierr);
  ierr = VecSetFromOptions(temp_y_vx); CHKERRQ(ierr);
  ierr = VecSetFromOptions(temp_x_vy); CHKERRQ(ierr);
  ierr = VecSetFromOptions(temp_y_vy); CHKERRQ(ierr);

/*
 * Initialize temp_x_vx, temp_y_vx, temp_x_vy, temp_y_vy to 0
 */
  ierr = VecSet(temp_x_vx, 0.0); CHKERRQ(ierr);
  ierr = VecSet(temp_y_vx, 0.0); CHKERRQ(ierr);
  ierr = VecSet(temp_x_vy, 0.0); CHKERRQ(ierr);
  ierr = VecSet(temp_y_vy, 0.0); CHKERRQ(ierr);

/*
 *      Form velocity field over the locally owned part of the grid
 */
  for (j=ys_volume; j<ys_volume+ym_volume; j++) {
    for (i=xs_volume; i<xs_volume+xm_volume; i++) {
      x = array_x_volume[j][i];
      y = array_y_volume[j][i];
      dx = array_dx_volume[j][i];
      dy = array_dy_volume[j][i];

      array_vx_left[j][i]  = findVx(x-dx*0.5,y,Lx,Ly);
      array_vx_right[j][i] = findVx(x+dx*0.5,y,Lx,Ly);
      array_vy_down[j][i]  = findVy(x,y-dy*0.5,Lx,Ly);
      array_vy_up[j][i]  = findVy(x,y+dy*0.5,Lx,Ly);

      // set x_vx, y_vx, x_vy, y_vy
      ierr = VecSetValue(temp_x_vx, i, x-dx*0.5, INSERT_VALUES); CHKERRQ(ierr);
      ierr = VecSetValue(temp_y_vx, j, y, INSERT_VALUES); CHKERRQ(ierr);

      ierr = VecSetValue(temp_x_vx, i+1, x+dx*0.5, INSERT_VALUES); CHKERRQ(ierr);
      ierr = VecSetValue(temp_y_vx, j, y, INSERT_VALUES); CHKERRQ(ierr);
      
      ierr = VecSetValue(temp_x_vy, i, x, INSERT_VALUES); CHKERRQ(ierr);
      ierr = VecSetValue(temp_y_vy, j, y-dy*0.5, INSERT_VALUES); CHKERRQ(ierr);

      ierr = VecSetValue(temp_x_vy, i, x, INSERT_VALUES); CHKERRQ(ierr);
      ierr = VecSetValue(temp_y_vy, j+1, y+dy*0.5, INSERT_VALUES); CHKERRQ(ierr);
    }
  }

/*
 *  Assembly temp_x_vx, temp_y_vx, temp_x_vy, temp_y_vy
 */
      ierr = VecAssemblyBegin(temp_x_vx); CHKERRQ(ierr); 
      ierr = VecAssemblyEnd(temp_x_vx); CHKERRQ(ierr);
      
      ierr = VecAssemblyBegin(temp_y_vx); CHKERRQ(ierr); 
      ierr = VecAssemblyEnd(temp_y_vx); CHKERRQ(ierr);

      ierr = VecAssemblyBegin(temp_x_vy); CHKERRQ(ierr); 
      ierr = VecAssemblyEnd(temp_x_vy); CHKERRQ(ierr);

      ierr = VecAssemblyBegin(temp_y_vy); CHKERRQ(ierr); 
      ierr = VecAssemblyEnd(temp_y_vy); CHKERRQ(ierr);

/*
 *  Create a copy of x_vx, y_vx, x_vy, y_vy on each processor
 *  for quicker access
 */
  ierr = VecScatterCreateToAll(temp_x_vx,&temp_ctx,&user->x_vx);CHKERRQ(ierr);
  ierr = VecScatterBegin(temp_ctx,temp_x_vx,user->x_vx,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
  ierr = VecScatterEnd(temp_ctx,temp_x_vx,user->x_vx,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
  ierr = VecScatterDestroy(&temp_ctx);CHKERRQ(ierr);
  ierr = VecDestroy(&temp_x_vx);CHKERRQ(ierr);

  ierr = VecScatterCreateToAll(temp_y_vx,&temp_ctx,&user->y_vx);CHKERRQ(ierr);
  ierr = VecScatterBegin(temp_ctx,temp_y_vx,user->y_vx,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
  ierr = VecScatterEnd(temp_ctx,temp_y_vx,user->y_vx,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
  ierr = VecScatterDestroy(&temp_ctx);CHKERRQ(ierr);
  ierr = VecDestroy(&temp_y_vx);CHKERRQ(ierr);

  ierr = VecScatterCreateToAll(temp_x_vy,&temp_ctx,&user->x_vy);CHKERRQ(ierr);
  ierr = VecScatterBegin(temp_ctx,temp_x_vy,user->x_vy,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
  ierr = VecScatterEnd(temp_ctx,temp_x_vy,user->x_vy,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
  ierr = VecScatterDestroy(&temp_ctx);CHKERRQ(ierr);
  ierr = VecDestroy(&temp_x_vy);CHKERRQ(ierr);

  ierr = VecScatterCreateToAll(temp_y_vy,&temp_ctx,&user->y_vy);CHKERRQ(ierr);
  ierr = VecScatterBegin(temp_ctx,temp_y_vy,user->y_vy,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
  ierr = VecScatterEnd(temp_ctx,temp_y_vy,user->y_vy,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
  ierr = VecScatterDestroy(&temp_ctx);CHKERRQ(ierr);
  ierr = VecDestroy(&temp_y_vy);CHKERRQ(ierr);


/*
 *      Restore vector
 */
  ierr = DMDAVecRestoreArray(user->da_volume_dof_1, user->vx_left_global, &array_vx_left);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da_volume_dof_1, user->vx_right_global, &array_vx_right);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da_volume_dof_1, user->vy_down_global, &array_vy_down);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da_volume_dof_1, user->vy_up_global, &array_vy_up);CHKERRQ(ierr);

  ierr = DMDAVecRestoreArray(user->da_volume_dof_1, user->x_volume_global, &array_x_volume);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da_volume_dof_1, user->y_volume_global, &array_y_volume);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da_volume_dof_1, user->dx_volume_global, &array_dx_volume);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da_volume_dof_1, user->dy_volume_global, &array_dy_volume);CHKERRQ(ierr);

  PetscFunctionReturn(0);
} 

/*
 * calculate vx based on (x,y) and total length (Lx,Ly)
 */
PetscReal findVx(PetscReal x, PetscReal y, PetscReal Lx, PetscReal Ly)
{
    PetscReal Vx_max, Vx;

    PetscFunctionBegin;
    Vx_max = VX_MAX;
    Vx = 4.0*y*(Ly - y)*Vx_max/(Ly*Ly);
    PetscFunctionReturn(Vx);
}

/*
 * calculate vy based on (x,y) and total length (Lx,Ly)
 */
PetscReal findVy(PetscReal x, PetscReal y, PetscReal Lx, PetscReal Ly)
{
    PetscFunctionBegin;
    PetscFunctionReturn(0.0);
}

#endif
