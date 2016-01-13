#ifndef _DATA_STRUCTURE_H_
#define _DATA_STRUCTURE_H_

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

#include "Settings.h"
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
   MPI_Comm  comm[2];  // MPI communicator, comm[0] for volume, comm[1] for surface
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

   Vec x_vx, y_vx; // vector to store the (x,y) coordinates of vx
   Vec x_vy, y_vy; // vector to store the (x,y) coordinates of vy
   
   bool initialized; // flag to be true if initialized

#ifdef PETSC_USE_LOG
   PetscClassId SURFACE_REACTION_CLASSID, VOLUME_REACTION_CLASSID, VOLUME_CONVECTION_CLASSID, VOLUME_DIFFUSION_CLASSID, STEP_CONTROL_CLASSID;
#endif

} MY_DATA ;

#endif
