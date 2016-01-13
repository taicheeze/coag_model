#ifndef _VOLUME_CONVECTION_SOLVER_H_
#define _VOLUME_CONVECTION_SOLVER_H_

#include "DataStructure.h"

typedef struct {
   TS      ts;   /* time-stepping solver context */
} VOLUME_CONVECTION_SOLVER;

PetscErrorCode InitializeVolumeConvectionSolver(VOLUME_CONVECTION_SOLVER *solver_ptr, MY_DATA *user);
PetscErrorCode DestroyVolumeConvectionSolver(VOLUME_CONVECTION_SOLVER *solver_ptr);

PetscErrorCode VolumeConvectionSolverStep(VOLUME_CONVECTION_SOLVER *solver_ptr, PetscReal t, PetscReal dt, MY_DATA *user);

PetscErrorCode VolumeConvectionRHS(TS ts,PetscReal t,Vec x,Vec f,void *ptr);

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 *     Form and set volume convection TS solver
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
PetscErrorCode InitializeVolumeConvectionSolver(VOLUME_CONVECTION_SOLVER *solver_ptr, MY_DATA *user)
{
    PetscErrorCode ierr;
    PetscInt       number_of_volume_species;

    PetscFunctionBegin;

    number_of_volume_species = NUMBER_OF_SPECIES_IN_VOLUME;

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 *     Create TS to manage hierarchical solvers
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    ierr = TSCreate(user->comm[0],&solver_ptr->ts);CHKERRQ(ierr);

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 *     Set TS properties
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    ierr = TSSetProblemType(solver_ptr->ts,TS_LINEAR);CHKERRQ(ierr); // linear problem
    ierr = TSSetType(solver_ptr->ts,TSEULER);CHKERRQ(ierr); // forward euler

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 *     Set volume species concentration vector to be solution vector
 *     Also set the residual vector of volume species concentration
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    ierr = TSSetSolution(solver_ptr->ts, user->c_volume_global); CHKERRQ(ierr);

    ierr = TSSetRHSFunction(solver_ptr->ts,user->r_volume_global,VolumeConvectionRHS,user);CHKERRQ(ierr); // right-hand-side function

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 *    Set up TS
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    ierr = TSSetUp(solver_ptr->ts);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 *     Destroy volume reaction TS solver
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
PetscErrorCode DestroyVolumeConvectionSolver(VOLUME_CONVECTION_SOLVER *solver_ptr)
{
    PetscErrorCode ierr;

    PetscFunctionBegin;

    ierr = TSDestroy(&(solver_ptr->ts)); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 *     Execute one volume convection step
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
PetscErrorCode VolumeConvectionSolverStep(VOLUME_CONVECTION_SOLVER *solver_ptr, PetscReal t, PetscReal dt, MY_DATA *user)
{
    PetscErrorCode   ierr;
    PetscReal        vx_max, dx, sub_dt;
    PetscInt         n_sub, i; // number of sub steps

    PetscFunctionBegin;

    vx_max = VX_MAX;
    dx     = DX_LENGTH;

    if(dt*vx_max > dx){
       n_sub = dt*vx_max/dx + 1;
       sub_dt = dt/(PetscReal)n_sub;
       for(i=0;i<n_sub;++i){
          ierr = TSSetTime(solver_ptr->ts,t+sub_dt*i);CHKERRQ(ierr);
          ierr = TSSetTimeStep(solver_ptr->ts,sub_dt);CHKERRQ(ierr);
    
//          VolumeConvectionRHS(solver_ptr->ts, t+sub_dt*i, user->c_volume_global, user->r_volume_global, user);
//  PetscPrintf(PETSC_COMM_WORLD,"t=%G,sub_dt=%G,i=%d,current_t=%G",t,sub_dt,i,t+sub_dt*i);
//  VecView(user->c_volume_global,PETSC_VIEWER_STDOUT_WORLD);
//  VecView(user->r_volume_global,PETSC_VIEWER_STDOUT_WORLD);
//          VecAXPY(user->c_volume_global, sub_dt, user->r_volume_global);   
          ierr = TSStep(solver_ptr->ts);CHKERRQ(ierr);
       }
//       SETERRQ(PETSC_COMM_SELF,1,"VolumeConvectionSolverStep(): vx*dt>dx, cannot use Forward Euler any more.\n");
    } else {
          ierr = TSSetTime(solver_ptr->ts,t);CHKERRQ(ierr);
          ierr = TSSetTimeStep(solver_ptr->ts,dt);CHKERRQ(ierr);

          ierr = TSStep(solver_ptr->ts);CHKERRQ(ierr);
//          VolumeConvectionRHS(solver_ptr->ts, t, user->c_volume_global, user->r_volume_global, user);
//  PetscPrintf(PETSC_COMM_WORLD,"t=%G,dt=%G",t,dt);
//  VecView(user->c_volume_global,PETSC_VIEWER_STDOUT_WORLD);
//  VecView(user->r_volume_global,PETSC_VIEWER_STDOUT_WORLD);
//          VecAXPY(user->c_volume_global, dt, user->r_volume_global);   
    }

    PetscFunctionReturn(0);
}

PetscErrorCode VolumeConvectionRHS(TS ts,PetscReal t,Vec c_global,Vec f_global,void *ptr)
{
    MY_DATA          *user = (MY_DATA *)ptr;  /* user-defined application context */
    PetscErrorCode   ierr;
    PetscInt         ys,ym,ye,xs,xm,xe;
    PetscInt         xM, yM;
    PetscInt         i,j,k,number_of_volume_species;
//    PetscMPIInt      rank,size;
    Vec              c_local, f_local;
    Species_Volume   **array_c_local, **array_f_local;
    PetscReal        **array_dx, **array_dy, **array_vx_left, **array_vx_right, **array_vy_down, **array_vy_up;

    PetscFunctionBegin;

    number_of_volume_species = NUMBER_OF_SPECIES_IN_VOLUME;

/*
 *   Get pointers to global vector data for dx,dy,vx,vy
 */
    ierr = DMDAVecGetArray(user->da_volume_dof_1, user->dx_volume_global, &array_dx);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->da_volume_dof_1, user->dy_volume_global, &array_dy);CHKERRQ(ierr);

    ierr = DMDAVecGetArray(user->da_volume_dof_1, user->vx_left_global, &array_vx_left);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->da_volume_dof_1, user->vx_right_global, &array_vx_right);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->da_volume_dof_1, user->vy_down_global, &array_vy_down);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->da_volume_dof_1, user->vy_up_global, &array_vy_up);CHKERRQ(ierr);

/*
 * Get global grid boundaries to handle boundary conditions
 *  xM: right-most x grid
 *  yM: up-most y grid
 */
    ierr = DMDAGetInfo(user->da_c_volume, NULL, &xM, &yM, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);CHKERRQ(ierr);

//    ierr = PetscPrintf(PETSC_COMM_WORLD,"xM:%d, yM:%d\n",xM, yM);CHKERRQ(ierr);

/*
 *  Get local grid boundaries (for 2-dimensional DMDA):
 *   xs, ys   - starting grid indices (no ghost points)
 *   xm, ym   - widths of local grid (no ghost points)
 */
    ierr = DMDAGetCorners(user->da_c_volume, &xs, &ys, PETSC_NULL, &xm, &ym, PETSC_NULL);CHKERRQ(ierr);
    
    ye = ys + ym;
    xe = xs + xm;

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Get ready for local function computations
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/*  
 *  Set local vectors
 */
    c_local = user->c_volume_local;
    f_local = user->r_volume_local;

/*
    Scatter ghost points to local vector, using the 2-step process
        DMGlobalToLocalBegin(), DMGlobalToLocalEnd().
    By placing code between these two statements, computations can be
    done while messages are in transition.
*/
    ierr = DMGlobalToLocalBegin(user->da_c_volume,c_global,INSERT_VALUES,c_local);CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(user->da_c_volume,c_global,INSERT_VALUES,c_local);CHKERRQ(ierr);

/*
    Get array from local vectors so we can access ghost point
*/
    ierr = DMDAVecGetArray(user->da_c_volume, c_local, &array_c_local);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->da_c_volume, f_local, &array_f_local);CHKERRQ(ierr);

/* 
 * right-hand-side of each species
 * c_up: upstream relative to v
 *
 * dx*dy*dc/dt = vx_left*c_up_vx_left*dy + vx_right*c_up_vx_right*dy + vy_down*c_up_vy_down*dx + vy_up*c_up_vy_up*dx
 *
 * dc/dt = vx_left*c_up_vx_left/dx + vx_right*c_up_vx_right/dx + vy_down*c_up_vy_down/dy + vy_up*c_up_vy_up/dy
 *
 * 
 */

    for (i=xs; i<xe; i++) {
        for (j=ys; j<ye; j++) {
//            ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD,"CELL: i:%d, j:%d, k:0, c_center:%e, c_left:%e, c_right:%e, c_up:%e, c_down:%e\n",i,j,array_c_local[j][i].c[0],array_c_local[j][i-1].c[0],array_c_local[j][i+1].c[0],array_c_local[j+1][i].c[0],array_c_local[j-1][i].c[0]);CHKERRQ(ierr);
//            ierr = PetscSynchronizedFlush(PETSC_COMM_WORLD);CHKERRQ(ierr);
            for(k=0;k<number_of_volume_species;++k){
                array_f_local[j][i].c[k] = 0.0;

#ifdef LEFT_RIGHT_PERIODIC
		// flow from left cell
                if( array_vx_left[j][i] > 0.0 ){ // left cell is upstream
                    array_f_local[j][i].c[k] += array_vx_left[j][i] * array_c_local[j][i-1].c[k] / array_dx[j][i];
                } else if( array_vx_left[j][i] < 0.0){ // center cell is upstream
                    array_f_local[j][i].c[k] += array_vx_left[j][i] * array_c_local[j][i].c[k] / array_dx[j][i];
                }
#else
		// flow from left cell
                if( i == 0 ){  // Left boundary
                    if( array_vx_left[j][i] > 0.0){
                        if( array_c_local[j][i-1].c[k] == OUTLET_GHOST_CELL_VALUE ){
                            SETERRQ(PETSC_COMM_SELF,1,"VolumeConvectionRHS(): There is an outlet on the left boundary, but vx >= 0, boundary conditions need to be reviewed.\n");
                        } else if( array_c_local[j][i-1].c[k] != WALL_GHOST_CELL_VALUE ){
                            array_f_local[j][i].c[k] += array_vx_left[j][i] * array_c_local[j][i-1].c[k] / array_dx[j][i];
                        }
                    } else if( array_vx_left[j][i] < 0.0) {
                        if( array_c_local[j][i-1].c[k] != WALL_GHOST_CELL_VALUE ){
                            array_f_local[j][i].c[k] += array_vx_left[j][i] * array_c_local[j][i].c[k] / array_dx[j][i];
                        }
                    }
                } else { // others
                    if( array_vx_left[j][i] > 0.0 ){ // left cell is upstream
                        array_f_local[j][i].c[k] += array_vx_left[j][i] * array_c_local[j][i-1].c[k] / array_dx[j][i];
                    } else if(array_vx_left[j][i] < 0.0) { // center cell is upstream
                        array_f_local[j][i].c[k] += array_vx_left[j][i] * array_c_local[j][i].c[k] / array_dx[j][i];
                    }
		}
#endif
//            if(k==0){
//               ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD,"AFTER LEFT FLOW: dc:%e\n",array_f_local[j][i].c[0]);CHKERRQ(ierr);
//            }

#ifdef LEFT_RIGHT_PERIODIC
		// flow from right cell
                if( array_vx_right[j][i] > 0.0 ){ // center cell is upstream
                    array_f_local[j][i].c[k] -= array_vx_right[j][i] * array_c_local[j][i].c[k] / array_dx[j][i];
                } else if( array_vx_right[j][i] < 0.0 ){ // right cell is upstream
                    array_f_local[j][i].c[k] -= array_vx_right[j][i] * array_c_local[j][i+1].c[k] / array_dx[j][i];
                }
#else
		// flow from right cell
                if( i == (xM-1) ){  // Right boundary
                    if( array_vx_right[j][i] > 0.0){
                        if( array_c_local[j][i+1].c[k] != WALL_GHOST_CELL_VALUE ){
                            array_f_local[j][i].c[k] -= array_vx_right[j][i] * array_c_local[j][i].c[k] / array_dx[j][i];
                        }
                    } else if( array_vx_right[j][i] < 0.0 ){
                        if( array_c_local[j][i+1].c[k] == OUTLET_GHOST_CELL_VALUE ){
                            SETERRQ(PETSC_COMM_SELF,1,"VolumeConvectionRHS(): There is an outlet on the right boundary, but vx < 0, boundary conditions need to be reviewed.\n");
                        } else if( array_c_local[j][i+1].c[k] != WALL_GHOST_CELL_VALUE ){
                            array_f_local[j][i].c[k] -= array_vx_right[j][i] * array_c_local[j][i+1].c[k] / array_dx[j][i];
                        }
                    }
                } else { // others
                    if( array_vx_right[j][i] > 0.0 ){ // center cell is upstream
                        array_f_local[j][i].c[k] -= array_vx_right[j][i] * array_c_local[j][i].c[k] / array_dx[j][i];
                    } else if( array_vx_right[j][i] < 0.0 ){ // right cell is upstream
                        array_f_local[j][i].c[k] -= array_vx_right[j][i] * array_c_local[j][i+1].c[k] / array_dx[j][i];
                    }
		}
#endif
//            if(k==0){
//                ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD,"AFTER RIGHT FLOW: dc:%e\n",array_f_local[j][i].c[0]);CHKERRQ(ierr);
//            }
		// flow from down cell
                if( j == 0 ){  // Lower boundary
                    if( array_vy_down[j][i] > 0.0){
                        if( array_c_local[j-1][i].c[k] == OUTLET_GHOST_CELL_VALUE ){
                            SETERRQ(PETSC_COMM_SELF,1,"VolumeConvectionRHS(): There is an outlet on the lower boundary, but vy >= 0, boundary conditions need to be reviewed.\n");
                        } else if( array_c_local[j-1][i].c[k] != WALL_GHOST_CELL_VALUE ){
                            array_f_local[j][i].c[k] += array_vy_down[j][i] * array_c_local[j-1][i].c[k] / array_dy[j][i];
                        }
                    } else if( array_vy_down[j][i] < 0.0) {
                        if( array_c_local[j-1][i].c[k] != WALL_GHOST_CELL_VALUE ){
                            array_f_local[j][i].c[k] += array_vy_down[j][i] * array_c_local[j][i].c[k] / array_dy[j][i];
                        }
                    }
                } else { // others
                    if( array_vy_down[j][i] > 0.0 ){ // lower cell is upstream
                        array_f_local[j][i].c[k] += array_vy_down[j][i] * array_c_local[j-1][i].c[k] / array_dy[j][i];
                    } else if( array_vy_down[j][i] < 0.0){ // center cell is upstream
                        array_f_local[j][i].c[k] += array_vy_down[j][i] * array_c_local[j][i].c[k] / array_dy[j][i];
                    }
		}
//            if(k==0){
//                ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD,"AFTER DOWN FLOW: dc:%e\n",array_f_local[j][i].c[0]);CHKERRQ(ierr);
//            }

		// flow from up cell
                if( j == (yM-1) ){  // Upper boundary
                    if( array_vy_up[j][i] > 0.0){
                        if( array_c_local[j+1][i].c[k] != WALL_GHOST_CELL_VALUE ){
                            array_f_local[j][i].c[k] -= array_vy_up[j][i] * array_c_local[j][i].c[k] / array_dy[j][i];
                        }
                    } else if( array_vy_up[j][i] < 0.0){
                        if( array_c_local[j+1][i].c[k] == OUTLET_GHOST_CELL_VALUE ){
                            SETERRQ(PETSC_COMM_SELF,1,"VolumeConvectionRHS(): There is an outlet on the up boundary, but vy < 0, boundary conditions need to be reviewed.\n");
                        } else if( array_c_local[j+1][i].c[k] != WALL_GHOST_CELL_VALUE ){
                            array_f_local[j][i].c[k] -= array_vy_up[j][i] * array_c_local[j+1][i].c[k] / array_dy[j][i];
                        }
                    }
                } else { // others
                    if( array_vy_up[j][i] > 0.0 ){ // center cell is upstream
                        array_f_local[j][i].c[k] -= array_vy_up[j][i] * array_c_local[j][i].c[k] / array_dy[j][i];
                    } else if( array_vy_up[j][i] < 0.0){ // upper cell is upstream
                        array_f_local[j][i].c[k] -= array_vy_up[j][i] * array_c_local[j+1][i].c[k] / array_dy[j][i];
                    }
		}
//            if(k==0){
//                ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD,"AFTER UP FLOW: vy_up:%e, dc:%e\n",array_vy_up[j][i],array_f_local[j][i].c[0]);CHKERRQ(ierr);
//            }
            }
        }
    }

/*
 *   Restore array to vector
 */
    ierr = DMDAVecRestoreArray(user->da_volume_dof_1, user->dx_volume_global, &array_dx);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->da_volume_dof_1, user->dy_volume_global, &array_dy);CHKERRQ(ierr);

    ierr = DMDAVecRestoreArray(user->da_volume_dof_1, user->vx_left_global, &array_vx_left);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->da_volume_dof_1, user->vx_right_global, &array_vx_right);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->da_volume_dof_1, user->vy_down_global, &array_vy_down);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->da_volume_dof_1, user->vy_up_global, &array_vy_up);CHKERRQ(ierr);

    ierr = DMDAVecRestoreArray(user->da_c_volume, c_local, &array_c_local);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->da_c_volume, f_local, &array_f_local);CHKERRQ(ierr);

/*
     Insert values from the local OUTPUT vector into the global 
     output vector
*/
    ierr = DMLocalToGlobalBegin(user->da_c_volume,f_local,INSERT_VALUES,f_global);CHKERRQ(ierr);
    ierr = DMLocalToGlobalEnd(user->da_c_volume,f_local,INSERT_VALUES,f_global);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

#endif
