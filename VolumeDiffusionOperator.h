#ifndef _VOLUME_DIFFUSION_OPERATOR_H_
#define _VOLUME_DIFFUSION_OPERATOR_H_

#include "DataStructure.h"
#include "VolumeDiffusionSolverStruct.h"

//
//  Calculate diffusion operator
//
//  Case 1: Constant dx and geometrically enlarged dy:
//  dy_1/dy_0 = dy_0/dy_-1 = dy_-1/dy_-2 = r > 1
//
//            _____.____
//            |         |
//            *    x    *             dy_1
//   ____.____|____.____|____.____
//  |         |         |         |
//  *    x    *    x    *    x    *   dy_0
//  |____.____|____.____|____.____|
//            |         |
//            *    x    *             dy_-1
//            |____.____|
//            |         |
//            *    x    *             dy_-2
//            |____.____|
//      dx        dx         dx
//  
//
//     dx * dy_0 * dc/dt 
//  =  - 1/dx   * (c_0,0 - c_-1,0) * D * dy_0
//     + 1/dx   * (c_1,0 - c_0,0)  * D * dy_0
//     - 1/dy_0 * [(c_0,0 - c_0,-1) * p- + (c_0,-1 - c_0,-2) * q-] * D * dx
//     + 1/dy_0 * [(c_0,1 - c_0,0)  * p+ + (c_0,0  - c_0,-1) * q+] * D * dx
//
//  where D: diffusion coefficient
//  p-: 2*r*(3*r+1) / (r+1)^3
//  q-: 2*r^3*(r-1) / (r+1)^3
//  p+: 2*(3*r+1)   / (r+1)^3
//  q+: 2*r^2*(r-1) / (r+1)^3
//
//  Full equation:
//  p-: 2*dy0*(3*dy-1 + dy-2) / ( (dy-1 + dy0)  * (dy0 + 2*dy-1 + dy-2) )
//  q-: 2*dy0*(dy0 - dy-1)    / ( (dy-2 + dy-1) * (dy0 + 2*dy-1 + dy-2) )
//  p+: 2*dy0*(3*dy0  + dy-1) / ( (dy0  + dy1)  * (dy1 + 2*dy0  + dy-1) )
//  q+: 2*dy0*(dy1 - dy0)     / ( (dy-1 + dy0)  * (dy1 + 2*dy0  + dy-1) )
//
//
//  p-,q-,p+,q+ are derived using interpolation
//
//     dc/dt 
//  =  - 1/dx^2 * (c_0,0 - c_-1,0) * D       
//     + 1/dx^2 * (c_1,0 - c_0,0)  * D       
//     - 1/dy_0^2 * [(c_0,0 - c_0,-1) * p- + (c_0,-1 - c_0,-2) * q-] * D     
//     + 1/dy_0^2 * [(c_0,1 - c_0,0)  * p+ + (c_0,0  - c_0,-1) * q+] * D     
//
//  Case 2: Constant dx and constant dy
//
//     dc/dt 
//  =  - 1/dx^2 * (c_0,0 - c_-1,0) * D       
//     + 1/dx^2 * (c_1,0 - c_0,0)  * D       
//     - 1/dy_0^2 * (c_0,0 - c_0,-1)  * D     
//     + 1/dy_0^2 * (c_0,1 - c_0,0)   * D     
//
//  Case 3: boundary condition: no flux
//   (lower boundary:  flux from surface reactions, but here we treat as no flux,
//                    add flux from surface reactions each step)
//
//  this function calculate p-, q-, p+, q+ using r and store in VolumeDiffusionSolverStruct
//  meanwhile create DMDA Operator Vector (iwup1,iw0,iwdown1,iwdown2, iwleft1, iwright1) 
//  (implicit weight of neighboring concentrations in diffusion operator for each grid)
//  to save calculation time each time step
//

PetscErrorCode FormVolumeDiffusionOperator(VOLUME_DIFFUSION_SOLVER *solver_ptr, MY_DATA *user);
PetscErrorCode DestroyVolumeDiffusionOperator(VOLUME_DIFFUSION_SOLVER *solver_ptr);

/*
 * calculate the volume diffusion operator DMDA and vector
 */
PetscErrorCode FormVolumeDiffusionOperator(VOLUME_DIFFUSION_SOLVER *solver_ptr, MY_DATA *user)
{
  PetscInt          i,j;
  PetscInt          Nx_volume,Ny_volume,xs_volume,ys_volume,xm_volume,ym_volume, xe_volume, ye_volume, Ny_bottom;
  PetscErrorCode    ierr;
  PetscReal         dy_volume_ratio;
  PetscReal         p_plus, q_plus, p_minus, q_minus, temp_power, temp_constant;
  PetscReal         **array_dx_volume, **array_dy_volume;
  PetscReal         D_volume;
  IMPLICIT_WEIGHTS  **array_implicit_weight;


  PetscFunctionBegin;

  Ny_bottom          = NUMBER_OF_Y_GRIDS_BOTTOM;
  dy_volume_ratio    = user->dy_ratio;
  D_volume           = user->D_volume;

// calculate p+, q+, p-, q-
  temp_power = pow(dy_volume_ratio+1.0, 3);
  p_minus = 2.0 * dy_volume_ratio * (3.0 * dy_volume_ratio + 1.0) / temp_power;
  q_minus = 2.0 * pow(dy_volume_ratio, 3) * (dy_volume_ratio - 1.0) / temp_power;
  p_plus  = 2.0 * (3.0 * dy_volume_ratio + 1.0) / temp_power;
  q_plus  = 2.0 * pow(dy_volume_ratio, 2) * (dy_volume_ratio - 1.0) / temp_power;


// create da_diffusion_operator, DM for diffusion operator
  ierr = DMDACreate2d(user->comm[0], DMDA_BOUNDARY_NONE, DMDA_BOUNDARY_NONE,DMDA_STENCIL_STAR,-NUMBER_OF_X_GRIDS,-NUMBER_OF_Y_GRIDS,PETSC_DECIDE,1,VOLUME_DIFFUSION_STENCIL_SIZE,2,PETSC_NULL,PETSC_NULL,&solver_ptr->da_diffusion_operator);CHKERRQ(ierr);

  ierr = DMDAGetInfo(solver_ptr->da_diffusion_operator, PETSC_IGNORE, &Nx_volume, &Ny_volume, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
                   PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE);CHKERRQ(ierr);

// create diffusion_operator vector
  ierr = DMCreateGlobalVector(solver_ptr->da_diffusion_operator,&solver_ptr->diffusion_operator);CHKERRQ(ierr);

  // set diffusion operator to be a zero vector
  ierr = VecSet(solver_ptr->diffusion_operator, 0.0); CHKERRQ(ierr);

/*
 *      Get a pointer to vector data.
 *      - For default PETSc vectors, VecGetArray() returns a pointer to
 *        the data array.  Otherwise, the routine is implementation dependent.
 *      - You MUST call VecRestoreArray() when you no longer need access to
 *        the array.
 */
 // get array for volume coordinate vectors
  ierr = DMDAVecGetArray(user->da_volume_dof_1, user->dx_volume_global, &array_dx_volume);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da_volume_dof_1, user->dy_volume_global, &array_dy_volume);CHKERRQ(ierr);

 // get array for diffusion operator vectors
  ierr = DMDAVecGetArray(solver_ptr->da_diffusion_operator, solver_ptr->diffusion_operator, &array_implicit_weight);CHKERRQ(ierr);
  

/*
 *      Get local grid boundaries (for 2-dimensional DMDA):
 *      xs, ys   - starting grid indices (no ghost points)
 *      xm, ym   - widths of local grid (no ghost points)
 *      xe, ye   - ending grid indices (no ghost points)
 *
 */
  ierr = DMDAGetCorners(solver_ptr->da_diffusion_operator, &xs_volume, &ys_volume, PETSC_NULL, &xm_volume, &ym_volume, PETSC_NULL);CHKERRQ(ierr);
  xe_volume = xs_volume + xm_volume;
  ye_volume = ys_volume + ym_volume;

/*
 *      Form implicit weights over the locally owned part of the grid
 */
/*
 * Since we only parallel along x directions, ys_volume should be 0 and ym_volume should be Ny_volume
 * this will make upper/lower boundary conditions easy to handle
 */
  if(ys_volume != 0 || ye_volume != Ny_volume){
      SETERRQ(PETSC_COMM_SELF,1,"FormVolumeDiffusionOperator(): Parallel should be only along x directions, please make sure ys_volume==0 and ye_volume==Ny_volume.\n");
  }

  // remember diffusion operator is a zero vector now

  // first put left and right flux in
  for (j=ys_volume; j<ye_volume; j++) {
      for (i=xs_volume; i<xe_volume; i++) {
          temp_constant = D_volume / (array_dx_volume[j][i] * array_dx_volume[j][i]); // temp_constant = D/(dx^2)
          // left flux except for left boundary
          if(i!=0){
              array_implicit_weight[j][i].iw0      -= temp_constant;
              array_implicit_weight[j][i].iwleft1  += temp_constant;
          }
          // right flux except for right boundary
          if(i!=Nx_volume-1){
              array_implicit_weight[j][i].iw0      -= temp_constant;
              array_implicit_weight[j][i].iwright1 += temp_constant;
          }
      }
  }

  // then put upper flux in
  for (i=xs_volume; i<xe_volume; i++) {
      // 0 <= j < Ny_bottom-1:
      // absolute constant dy's
      for (j=0; j<Ny_bottom-1; j++) {
          temp_constant = D_volume / (array_dy_volume[j][i] * array_dy_volume[j][i]); // temp_constant = D/(dy^2)
          array_implicit_weight[j][i].iw0      -= temp_constant;
          array_implicit_weight[j][i].iwup1    += temp_constant;
      }

      // j = Ny_bottom-1:
      // variate dy from up but constant dy from down
      j = Ny_bottom-1;
          temp_constant = D_volume / (array_dy_volume[j][i] * array_dy_volume[j][i]); // temp_constant = D/(dy^2)
          array_implicit_weight[j][i].iwup1    += temp_constant * 8.0 / ((1.0 + dy_volume_ratio) * (3.0 + dy_volume_ratio));
          array_implicit_weight[j][i].iw0      -= temp_constant * (3.0 - dy_volume_ratio) / (1.0 + dy_volume_ratio);
          array_implicit_weight[j][i].iwdown1  -= temp_constant * (dy_volume_ratio - 1.0) / (dy_volume_ratio + 3.0);

      // Ny_bottom-1 < j < Ny_volume-1
      // variate dy
      for (j=Ny_bottom; j<Ny_volume-1; j++) {
          temp_constant = D_volume / (array_dy_volume[j][i] * array_dy_volume[j][i]); // temp_constant = D/(dy^2)
          array_implicit_weight[j][i].iwup1    += temp_constant * p_plus;
          array_implicit_weight[j][i].iw0      -= temp_constant * (p_plus - q_plus);
          array_implicit_weight[j][i].iwdown1  -= temp_constant * q_plus;
      }

      // j = Ny_volume-1:
      // no upper flux

  }

  // then put lower flux in
  for (i=xs_volume; i<xe_volume; i++) {
      // j=0:
      // no lower flux

      // 0 < j < Ny_bottom:
      // absolute constant dy's
      // j = Ny_bottom-1:
      // variate dy from up but constant dy from down
      for (j=1; j<Ny_bottom; j++) {
          temp_constant = D_volume / (array_dy_volume[j][i] * array_dy_volume[j][i]); // temp_constant = D/(dy^2)
          array_implicit_weight[j][i].iw0      -= temp_constant;
          array_implicit_weight[j][i].iwdown1  += temp_constant;
      }

      // j = Ny_bottom:
      // variate dy from up but the two lower cells have constant dy
      j = Ny_bottom;
          temp_constant = D_volume / (array_dy_volume[j][i] * array_dy_volume[j][i]); // temp_constant = D/(dy^2)
          array_implicit_weight[j][i].iw0      -= temp_constant * 8.0 * dy_volume_ratio / ((1.0 + dy_volume_ratio) * (3.0 + dy_volume_ratio));
          array_implicit_weight[j][i].iwdown1  += temp_constant * dy_volume_ratio * (3.0 - dy_volume_ratio) / (1.0 + dy_volume_ratio);
          array_implicit_weight[j][i].iwdown2  += temp_constant * dy_volume_ratio * (dy_volume_ratio - 1.0) / (dy_volume_ratio + 3.0);
      // Ny_bottom < j < Ny_volume
      // variate dy
      for (j=Ny_bottom+1; j<Ny_volume; j++) {
          temp_constant = D_volume / (array_dy_volume[j][i] * array_dy_volume[j][i]); // temp_constant = D/(dy^2)
          array_implicit_weight[j][i].iw0      -= temp_constant * p_minus;
          array_implicit_weight[j][i].iwdown1  += temp_constant * (p_minus - q_minus);
          array_implicit_weight[j][i].iwdown2  += temp_constant * q_minus;
      }
  }
  // j==0:
  // lower boundary, no flux from lower boundary
  // a stencil of size 4

  //  0<j<Ny_bottom-1:
  //  absolute constant dy's
  //  a stencil of size 5
  
  //  j=Ny_bottom-1:
  //  constant dy from lower boundary but variate dy from upper boundary
  //  a stencil of size 5

  //  j=Ny_bottom:
  //  the two lower cells have constant dy
  //  a stencil of size 6

  //  Ny_bottom < j < Ny_volume-1:
  //  variate dy
  //  a stencil of size 6
  
  //  j = Ny_volume-1:
  //  uppwer boundary, no flux from upper boundary
  //  a stencil of size 5

/*
 *      Restore vectors
 */
  ierr = DMDAVecRestoreArray(user->da_volume_dof_1, user->dx_volume_global, &array_dx_volume);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da_volume_dof_1, user->dy_volume_global, &array_dy_volume);CHKERRQ(ierr);

  ierr = DMDAVecRestoreArray(solver_ptr->da_diffusion_operator, solver_ptr->diffusion_operator, &array_implicit_weight);CHKERRQ(ierr);

  PetscFunctionReturn(0);
} 

PetscErrorCode DestroyVolumeDiffusionOperator(VOLUME_DIFFUSION_SOLVER *solver_ptr)
{
  PetscErrorCode  ierr;

  PetscFunctionBegin;

  ierr = VecDestroy(&solver_ptr->diffusion_operator);CHKERRQ(ierr);
  ierr = DMDestroy(&solver_ptr->da_diffusion_operator);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

#endif
