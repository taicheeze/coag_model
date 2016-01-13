#ifndef _DATA_STRUCTURE_H_
#define _DATA_STRUCTURE_H_

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
 *
 * Geometry
 *
 */

// x grids are equally spaced
// y grids are equally spaced near the bottom (with a small y), and are geometrically enlarged upwards out of the bottom region
#define NUMBER_OF_X_GRIDS         100 
#define NUMBER_OF_Y_GRIDS          25
#define NUMBER_OF_Y_GRIDS_BOTTOM   10 // this has to be at least 2 or it will cause error in diffusion solver

#define NUMBER_OF_GRIDS             NUMBER_OF_X_GRIDS*NUMBER_OF_Y_GRIDS

#define X_LENGTH          20000.0e-6
#define Y_LENGTH            250.0e-6
#define Y_BOTTOM_LENGTH      20.0e-6
#define BOTTOM_LAYER_DY      Y_BOTTOM_LENGTH/(PetscReal)NUMBER_OF_Y_GRIDS_BOTTOM

#define DX_LENGTH        X_LENGTH/(PetscReal)NUMBER_OF_X_GRIDS

#define NUMBER_OF_PATCHES 2

const PetscReal PATCH_LEFT_BOUNDARIES[NUMBER_OF_PATCHES]  = {3.0e-3, 5.0e-3};
const PetscReal PATCH_RIGHT_BOUNDARIES[NUMBER_OF_PATCHES] = {5.0e-3, 20.0e-3};

/*
 * boundary conditions
 */
// define LEFT_RIGHT_PERIODIC if it is periodic boundary condition on left/right wall
#define LEFT_RIGHT_PERIODIC //uncomment to add periodic boundary conditions

// if LEFT_RIGHT_PERIODIC is not defined there is an inlet and an outlet

// define OUTLET_GHOST_CELL_VALUE and WALL_GHOST_CELL_VALUE to be a negative number for control
#define OUTLET_GHOST_CELL_VALUE -2
#define WALL_GHOST_CELL_VALUE   -1

// not represented in absolute length but in the number of grids
const PetscInt INLET_X_BOUNDARIES[2] = {0, 0};
const PetscInt INLET_Y_BOUNDARIES[2] = {0, NUMBER_OF_Y_GRIDS-1};

const PetscInt OUTLET_X_BOUNDARIES[2] = {NUMBER_OF_X_GRIDS-1, NUMBER_OF_X_GRIDS-1};
const PetscInt OUTLET_Y_BOUNDARIES[2] = {0, NUMBER_OF_Y_GRIDS-1};

/*
 *
 *  Time stepping
 *
 */

#define T_START          0.0
#define T_END         2000.0
#define DT_MIN           0.5e-4
#define DT_MAX           0.150
//#define DT_STEP          0.01

#define OUTPUT_DT       10.0

// max change time step controller, define to use it, undefine to use others
#define USE_MAX_CHANGE_STEP_CONTROLLER

#define C_RELATIVE_TOL_MIN  0.5e-3
#define C_RELATIVE_TOL_MAX  1.0e-3

#define TIMESTEP_INCREASE_RATIO  1.01
#define TIMESTEP_DECREASE_RATIO  0.9

// double stepping time step controller, define to use it, undefine to use others
//#define USE_STEP_DOUBLING_STEP_CONTROLLER


/*
 *
 * Constants and tolerances
 *
 */

// absolute error tolerance 
#define ATOL 1.0e-9
// relative error tolerance
#define RTOL 0.5e-2

// safety fraction
#define SAFETY_FRACTION 0.9

//#define VX_MAX           1.25e-4
#define VX_MAX           1.25e-3

#define FLUX_CONVERSION_RATIO  1.0  // boundary flux constant due to unit conversion

//#define VOLUME_DIFFUSION_CONSTANT  1.0e-9
#define VOLUME_DIFFUSION_CONSTANT  5.0e-11

// error tolerance for any calculation
#define CALCULATION_ERROR_TOLERANCE   1.0e-14


#include "InitialConcentrations.h"

#include <stdbool.h>
#include <math.h>
#include <petscdmda.h>
#include <petscts.h>
#include <petscsys.h> 

typedef struct {
   PetscReal c[NUMBER_OF_SPECIES_IN_VOLUME];
} Species_Volume;

typedef struct {
   PetscReal c[NUMBER_OF_SPECIES_ON_SURFACE];
} Species_Surface;

typedef struct {
   MPI_Comm  comm[2];  // MPI communicator, comm[1] for volume, comm[2] for surface
                       // set while creating DM
                       // will use PETSC_COMM_WORLD as default for both now

   DM  da_c_volume;  // DM of species in volume (distributed data management)
   Vec c_volume_global; // species in volume
   Vec c_volume_local; // local vector with ghost cells 
   Vec r_volume_global; // residual vector for species in volume
   Vec r_volume_local; // local residual vector with ghost cells 

   DM  da_volume_dof_1; // DM of any vector in volume of dof(degree of freedom)=1
   Vec x_volume_global, y_volume_global, dx_volume_global, dy_volume_global; // coordinates in volume
   Vec vx_left_global, vx_right_global, vy_down_global, vy_up_global; // velocity in volume

   DM  da_c_surface; // DM of species on surface (distributed data management)
   Vec c_surface_global; // species on surface
   Vec c_surface_local; // local vector with ghost cells
   Vec r_surface_global; // residual vector for species on surface
   Vec r_surface_local; // local residual vector with ghost cells

   DM  da_surface_dof_1; // DM of any vector on surface of dof(degree of freedom)=1
   Vec x_surface_global, dx_surface_global; // coordinates on surface

   PetscReal  dy0;       // dy of the bottom layer grids
   PetscReal  dy_ratio;  // the ratio of  dy1/dy0 for geometrically enlarged grids, initiliazed in FormGrid() 
   PetscReal  D_volume;  // volume diffusion constant

   Vec c_volume_single; // global vector of a single species in volume, mainly for output, template, and error/change controls
   VecScatter volume_species_scatter[NUMBER_OF_SPECIES_IN_VOLUME]; // vecscatter for volume species

   Vec c_surface_single; // global vector of a single species on surface, mainly for output, template, and error/change controls
   VecScatter surface_species_scatter[NUMBER_OF_SPECIES_ON_SURFACE]; // vecscatter for surface species

   bool initialized; // flag to be true if initialized

#ifdef PETSC_USE_LOG
   PetscClassId SURFACE_REACTION_CLASSID, VOLUME_REACTION_CLASSID, VOLUME_CONVECTION_CLASSID, VOLUME_DIFFUSION_CLASSID, STEP_CONTROL_CLASSID;
#endif

} MY_DATA ;

#endif
