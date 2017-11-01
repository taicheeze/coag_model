#ifndef _SETTINGS_H_
#define _SETTINGS_H_
//#define DEBUG
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
#define NUMBER_OF_X_GRIDS         50 //reduced from 100 
#define NUMBER_OF_Y_GRIDS          25 //reduced from 25
#define NUMBER_OF_Y_GRIDS_BOTTOM   10 // this has to be at least 2 or it will cause error in diffusion solver ,reduced from 10

#define NUMBER_OF_GRIDS             NUMBER_OF_X_GRIDS*NUMBER_OF_Y_GRIDS

#define X_LENGTH          10000.0e-6 //reduced from 20000e-6
#define Y_LENGTH            10.0e-6
//reduced from 250e-6
#define Y_BOTTOM_LENGTH      1.0e-6
//reduced from 20e-6
#define BOTTOM_LAYER_DY      Y_BOTTOM_LENGTH/(PetscReal)NUMBER_OF_Y_GRIDS_BOTTOM

#define DX_LENGTH        X_LENGTH/(PetscReal)NUMBER_OF_X_GRIDS

#define NUMBER_OF_PATCHES 3

const PetscReal PATCH_LEFT_BOUNDARIES[NUMBER_OF_PATCHES]  = {0.0e-3, 3.0e-3, 5.0e-3};
const PetscReal PATCH_RIGHT_BOUNDARIES[NUMBER_OF_PATCHES] = {3.0e-3, 5.0e-3, 10.0e-3};

/*
 * boundary conditions
 */
// define LEFT_RIGHT_PERIODIC if it is periodic boundary condition on left/right wall
#define LEFT_RIGHT_PERIODIC

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
#define T_END         5000.0
#define DT_MIN           0.5e-4
#define DT_MAX           1.0
//#define DT_STEP          0.01

#define OUTPUT_DT       20.0

// max change time step controller, define to use it, undefine to use others
#define USE_MAX_CHANGE_STEP_CONTROLLER

#define C_RELATIVE_TOL_MIN  0.5e-3
#define C_RELATIVE_TOL_MAX  1.0e-3

#define TIMESTEP_INCREASE_RATIO  1.1
#define TIMESTEP_DECREASE_RATIO  0.9

// double stepping time step controller, define to use it, undefine to use others
//#define USE_STEP_DOUBLING_STEP_CONTROLLER


/*
 *
 * Constants and tolerances
 *
 */

// absolute error tolerance 
#define ATOL 1.0e-9 //relaxed from 1e-9
// relative error tolerance
#define RTOL 5.0e-3 //lowered from 5e-3

// safety fraction
#define SAFETY_FRACTION 0.9

//#define VX_MAX           5.0e-4 for 20 shear rate for 100um channel
#define VX_MAX           5e-3

//5.0e-4 is 20 shear for 100um and below

#define FLUX_CONVERSION_RATIO  1.0  // boundary flux constant due to unit conversion

//#define VOLUME_DIFFUSION_CONSTANT  1.0e-11
#define VOLUME_DIFFUSION_CONSTANT  5.0e-11

// error tolerance for any calculation
#define CALCULATION_ERROR_TOLERANCE   1.0e-14

#endif
