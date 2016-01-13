
static char help[] = "Coagulation PDE in 2d.\n\n";

/*T
   Concepts: DMDA^using distributed arrays;
   Processors: n
T*/

/* 
   Include "petscdmda.h" so that we can use distributed arrays (DMDAs).
   Include "petscsnes.h" so that we can use SNES solvers. 
*/
#include <petscdmda.h>
#include <petscsnes.h>
#include "Headers.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
//  PetscInt                     its;                  /* iterations for convergence */
  PetscErrorCode               ierr;
  MY_DATA                      my_data;
  PetscReal                    t_start, t_end, t, dt, output_dt;
  PetscReal                    dt_min, dt_max;
  PetscReal                    min_volume, min_surface;
  PetscInt                     min_vol_index, min_sur_index;
  SURFACE_REACTION_SOLVER      surface_reaction_solver;
  VOLUME_REACTION_SOLVER       volume_reaction_solver;
  VOLUME_CONVECTION_SOLVER     volume_convection_solver;
  VOLUME_DIFFUSION_SOLVER      volume_diffusion_solver;
  TIME_STEP_CONTROLLER         time_step_controller;
  OUTPUT_HANDLER               output_handler;
//  PetscBool                    flg = PETSC_FALSE;
//  PetscBool                    matlab_function = PETSC_FALSE;
#ifdef PETSC_USE_LOG
  PetscLogEvent                surface_reaction, volume_reaction, volume_convection, volume_diffusion, step_control;
#endif

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Initialize program
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscInitialize(&argc,&argv,(char *)0,help);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Initialize my_data
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = InitializeMY_DATA(&my_data);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Set grids, velocity field, and initial solutions
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = FormGrid(&my_data);CHKERRQ(ierr);
  ierr = FormVelocityField(&my_data);CHKERRQ(ierr);
  ierr = FormInitialSolution(&my_data);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Set surface/volume reaction solver, volume convection solver,
         boundary flux factor solver
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = InitializeSurfaceReactionSolver(&surface_reaction_solver, &my_data);CHKERRQ(ierr);
  ierr = InitializeVolumeReactionSolver(&volume_reaction_solver, &my_data);CHKERRQ(ierr);
  ierr = InitializeVolumeConvectionSolver(&volume_convection_solver, &my_data);CHKERRQ(ierr);
  ierr = InitializeVolumeDiffusionSolver(&volume_diffusion_solver, &my_data);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Set time and time-step
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  t_start = T_START;
  t_end   = T_END;

  dt_min  = DT_MIN;
  dt_max  = DT_MAX;
  
  t  = t_start;
  dt = dt_min;

  output_dt = OUTPUT_DT;

  ierr = InitializeTimeStepController(t, dt, &time_step_controller, &my_data);CHKERRQ(ierr);
  ierr = InitializeOutputHandler(&output_handler, output_dt, "output_species_list.txt");CHKERRQ(ierr);

/*
  for(k=0;k<number_of_volume_species;++k){
      ierr = getSingleSpeciesVec(my_data.c_volume_global, my_data.c_volume_single, my_data.volume_species_scatter[k]); CHKERRQ(ierr);
      ierr = VecNorm(my_data.c_volume_global,NORM_MAX,&current_max_c_volume[k]); CHKERRQ(ierr);
      old_max_c_volume[k] = current_max_c_volume[k];      
  }
  for(k=0;k<number_of_surface_species;++k){
      ierr = getSingleSpeciesVec(my_data.c_surface_global, my_data.c_surface_single, my_data.surface_species_scatter[k]); CHKERRQ(ierr);
      ierr = VecNorm(my_data.c_surface_global,NORM_MAX,&current_max_c_surface[k]); CHKERRQ(ierr);
      old_max_c_surface[k] = current_max_c_surface[k];      
  }
*/
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Register event to track efficiency of solvers and step control
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
#ifdef PETSC_USE_LOG
  ierr = PetscLogEventRegister("Surface Reaction",my_data.SURFACE_REACTION_CLASSID,&surface_reaction); CHKERRQ(ierr);
  ierr = PetscLogEventRegister("Volume Reaction",my_data.VOLUME_REACTION_CLASSID,&volume_reaction); CHKERRQ(ierr);
  ierr = PetscLogEventRegister("Volume Convection",my_data.VOLUME_CONVECTION_CLASSID,&volume_convection); CHKERRQ(ierr);
  ierr = PetscLogEventRegister("Volume Diffusion",my_data.VOLUME_DIFFUSION_CLASSID,&volume_diffusion); CHKERRQ(ierr);
  ierr = PetscLogEventRegister("Step Control",my_data.STEP_CONTROL_CLASSID,&step_control); CHKERRQ(ierr);
#endif

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Solve system
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/*
  PetscViewer viewer;
  char     output_filename[50];
//  getSingleSpeciesVec(my_data.c_volume_global, my_data.c_volume_single, my_data.volume_species_scatter[0]);
//  PetscObjectSetName((PetscObject)my_data.c_volume_single,"volumeconcentration");
//  VecView(my_data.c_volume_single,viewer);
//  VecView(my_data.dy_volume_global,viewer);
  sprintf(output_filename, "coordinates.m");
  PetscViewerASCIIOpen(PETSC_COMM_WORLD, output_filename, &viewer);
  PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
  PetscObjectSetName((PetscObject)my_data.x_volume_global,"x_volume_global");
  VecView(my_data.x_volume_global,viewer);
  PetscObjectSetName((PetscObject)my_data.dx_volume_global,"dx_volume_global");
  VecView(my_data.dx_volume_global,viewer);
  PetscObjectSetName((PetscObject)my_data.y_volume_global,"y_volume_global");
  VecView(my_data.y_volume_global,viewer);
  PetscObjectSetName((PetscObject)my_data.dy_volume_global,"dy_volume_global");
  VecView(my_data.dy_volume_global,viewer);
  PetscViewerDestroy(&viewer);
*/
  ierr = OutputHandlerOtherOutput(&my_data); CHKERRQ(ierr);
  ierr = OutputHandlerSingleOutput(&output_handler, t, &my_data); CHKERRQ(ierr);

/*
  sprintf(output_filename, "concentration_volume_%d.m",(int)next_output);
  PetscViewerASCIIOpen(PETSC_COMM_WORLD, output_filename, &viewer);
  PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);

  getSingleSpeciesVec(my_data.c_volume_global, my_data.c_volume_single, my_data.volume_species_scatter[0]);
  PetscObjectSetName((PetscObject)my_data.c_volume_single,"VII_concentration");
  VecView(my_data.c_volume_single,viewer);
  getSingleSpeciesVec(my_data.c_volume_global, my_data.c_volume_single, my_data.volume_species_scatter[3]);
  PetscObjectSetName((PetscObject)my_data.c_volume_single,"IIa_concentration");
  VecView(my_data.c_volume_single,viewer);
  getSingleSpeciesVec(my_data.c_volume_global, my_data.c_volume_single, my_data.volume_species_scatter[7]);
  PetscObjectSetName((PetscObject)my_data.c_volume_single,"II_concentration");
  VecView(my_data.c_volume_single,viewer);
  PetscViewerDestroy(&viewer);
*/
/*
  PetscObjectSetName((PetscObject)my_data.c_volume_global,"volume_concentration");
  VecView(my_data.c_volume_global,viewer);
  PetscObjectSetName((PetscObject)my_data.c_surface_global,"surface_concentration");
  VecView(my_data.c_surface_global,viewer);

    PetscViewerASCIIOpen(PETSC_COMM_WORLD, output_filename, &viewer);
    PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
*/

/*
  sprintf(output_filename, "concentrations.m");
  PetscViewerASCIIOpen(PETSC_COMM_WORLD, output_filename, &viewer);
  PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
*/

#ifdef DEBUG
    ierr = PetscPrintf(PETSC_COMM_WORLD,"DEBUG mode\n");CHKERRQ(ierr);
    SNESConvergedReason reason;
    PetscInt            its;
    Vec                 tmp_x, tmp_dx;
    PetscReal           problem_t;

    problem_t = 10000.0;
#endif

  while(t<t_end){
    ierr = PetscPrintf(PETSC_COMM_WORLD,"t=%e\n",t);CHKERRQ(ierr);

#ifdef NON_STEADY_FLOW
    ierr = UpdateVelocityField(t, &my_data);CHKERRQ(ierr);
#endif

    while(time_step_controller.flg!=0){
#ifdef PETSC_USE_LOG
      ierr = PetscLogEventBegin(surface_reaction,0,0,0,0); CHKERRQ(ierr);
#endif
      ierr = SurfaceReactionSolverStep(&surface_reaction_solver, time_step_controller.tmp_t, time_step_controller.tmp_dt, &my_data); CHKERRQ(ierr);

      ierr = VecMin(my_data.c_volume_global,&min_vol_index,&min_volume); CHKERRQ(ierr);
      ierr = VecMin(my_data.c_surface_global,&min_sur_index,&min_surface); CHKERRQ(ierr);

      while(min_volume<0.0 && min_volume > -1.0e-300){
#ifdef DEBUG
          ierr = PetscPrintf(PETSC_COMM_WORLD,"DEBUG: corrected min_vol: %e, negative concentration at t=%e, after surface reaction\n", min_volume, t);CHKERRQ(ierr);
#endif
          ierr = VecSetValue(my_data.c_volume_global, min_vol_index, 0.0, INSERT_VALUES);CHKERRQ(ierr);
          ierr = VecAssemblyBegin(my_data.c_volume_global);CHKERRQ(ierr);
          ierr = VecAssemblyEnd(my_data.c_volume_global);CHKERRQ(ierr);
          ierr = VecMin(my_data.c_volume_global,&min_vol_index,&min_volume); CHKERRQ(ierr);
      }
      while(min_surface<0.0 && min_surface > -1.0e-300){
#ifdef DEBUG
          ierr = PetscPrintf(PETSC_COMM_WORLD,"DEBUG: corrected min_sur: %e, negative concentration at t=%e, after surface reaction\n", min_surface, t);CHKERRQ(ierr);
#endif
          ierr = VecSetValue(my_data.c_surface_global, min_sur_index, 0.0, INSERT_VALUES);CHKERRQ(ierr);
          ierr = VecAssemblyBegin(my_data.c_surface_global);CHKERRQ(ierr);
          ierr = VecAssemblyEnd(my_data.c_surface_global);CHKERRQ(ierr);
          ierr = VecMin(my_data.c_surface_global,&min_sur_index,&min_surface); CHKERRQ(ierr);
      }
      if(min_volume<0.0 || min_surface<0.0){
          ierr = PetscPrintf(PETSC_COMM_WORLD,"DEBUG: min_vol: %e, min_sur: %e, negative concentration at t=%e, after surface reaction\n", min_volume, min_surface, t);CHKERRQ(ierr);
/*
          ierr = SNESGetIterationNumber(surface_reaction_solver.snes,&its);CHKERRQ(ierr);
          ierr = SNESGetConvergedReason(surface_reaction_solver.snes, &reason);CHKERRQ(ierr);
          ierr = PetscPrintf(PETSC_COMM_WORLD,"Number of SNES iterations = %D, %s\n",its,SNESConvergedReasons[reason]);CHKERRQ(ierr);
          ierr = SNESGetSolution(surface_reaction_solver.snes,&tmp_x);CHKERRQ(ierr);
          ierr = VecView(tmp_x,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
          ierr = SNESGetSolutionUpdate(surface_reaction_solver.snes,&tmp_dx);CHKERRQ(ierr);
          ierr = VecView(tmp_dx,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
          ierr = OutputHandlerSingleOutput(&output_handler, t, &my_data); CHKERRQ(ierr);
*/
//          SETERRQ(PETSC_COMM_WORLD,-1,"Negative concentrations.\n");
      }
#ifdef DEBUG
      if(t>problem_t){
          ierr = OutputHandlerSingleOutput(&output_handler, t, &my_data); CHKERRQ(ierr);
      }
#endif

#ifdef PETSC_USE_LOG
      ierr = PetscLogEventEnd(surface_reaction,0,0,0,0); CHKERRQ(ierr);
#endif

#ifdef PETSC_USE_LOG
      ierr = PetscLogEventBegin(volume_reaction,0,0,0,0); CHKERRQ(ierr);
#endif
      ierr = VolumeReactionSolverStep(&volume_reaction_solver, time_step_controller.tmp_t, time_step_controller.tmp_dt, &my_data); CHKERRQ(ierr);

      ierr = VecMin(my_data.c_volume_global,&min_vol_index,&min_volume); CHKERRQ(ierr);
      ierr = VecMin(my_data.c_surface_global,&min_sur_index,&min_surface); CHKERRQ(ierr);

      while(min_volume<0.0 && min_volume > -1.0e-300){
#ifdef DEBUG
          ierr = PetscPrintf(PETSC_COMM_WORLD,"DEBUG: corrected min_vol: %e, negative concentration at t=%e, after volume reaction\n", min_volume, t);CHKERRQ(ierr);
#endif
          ierr = VecSetValue(my_data.c_volume_global, min_vol_index, 0.0, INSERT_VALUES);CHKERRQ(ierr);
          ierr = VecAssemblyBegin(my_data.c_volume_global);CHKERRQ(ierr);
          ierr = VecAssemblyEnd(my_data.c_volume_global);CHKERRQ(ierr);
          ierr = VecMin(my_data.c_volume_global,&min_vol_index,&min_volume); CHKERRQ(ierr);
      }
      while(min_surface<0.0 && min_surface > -1.0e-300){
#ifdef DEBUG
          ierr = PetscPrintf(PETSC_COMM_WORLD,"DEBUG: corrected min_sur: %e, negative concentration at t=%e, after volume reaction\n", min_surface, t);CHKERRQ(ierr);
#endif
          ierr = VecSetValue(my_data.c_surface_global, min_sur_index, 0.0, INSERT_VALUES);CHKERRQ(ierr);
          ierr = VecAssemblyBegin(my_data.c_surface_global);CHKERRQ(ierr);
          ierr = VecAssemblyEnd(my_data.c_surface_global);CHKERRQ(ierr);
          ierr = VecMin(my_data.c_surface_global,&min_sur_index,&min_surface); CHKERRQ(ierr);
      }
      if(min_volume<0.0 || min_surface<0.0){
          ierr = PetscPrintf(PETSC_COMM_WORLD,"DEBUG: min_vol: %e, min_sur: %e, negative concentration at t=%e, after volume reaction\n", min_volume, min_surface, t);CHKERRQ(ierr);
//          ierr = OutputHandlerSingleOutput(&output_handler, t, &my_data); CHKERRQ(ierr);
//          SETERRQ(PETSC_COMM_WORLD,-1,"Negative concentrations.\n");
      }
#ifdef DEBUG
      if(t>problem_t){
          ierr = OutputHandlerSingleOutput(&output_handler, t, &my_data); CHKERRQ(ierr);
      }
#endif

#ifdef PETSC_USE_LOG
      ierr = PetscLogEventEnd(volume_reaction,0,0,0,0); CHKERRQ(ierr);
#endif

#ifdef PETSC_USE_LOG
      ierr = PetscLogEventBegin(volume_convection,0,0,0,0); CHKERRQ(ierr);
#endif
      ierr = VolumeConvectionSolverStep(&volume_convection_solver, time_step_controller.tmp_t, time_step_controller.tmp_dt, &my_data); CHKERRQ(ierr);

      ierr = VecMin(my_data.c_volume_global,&min_vol_index,&min_volume); CHKERRQ(ierr);
      ierr = VecMin(my_data.c_surface_global,&min_sur_index,&min_surface); CHKERRQ(ierr);

      while(min_volume<0.0 && min_volume > -1.0e-300){
#ifdef DEBUG
          ierr = PetscPrintf(PETSC_COMM_WORLD,"DEBUG: corrected min_vol: %e, negative concentration at t=%e, after volume convection\n", min_volume, t);CHKERRQ(ierr);
#endif
          ierr = VecSetValue(my_data.c_volume_global, min_vol_index, 0.0, INSERT_VALUES);CHKERRQ(ierr);
          ierr = VecAssemblyBegin(my_data.c_volume_global);CHKERRQ(ierr);
          ierr = VecAssemblyEnd(my_data.c_volume_global);CHKERRQ(ierr);
          ierr = VecMin(my_data.c_volume_global,&min_vol_index,&min_volume); CHKERRQ(ierr);
      }
      while(min_surface<0.0 && min_surface > -1.0e-300){
#ifdef DEBUG
          ierr = PetscPrintf(PETSC_COMM_WORLD,"DEBUG: corrected min_sur: %e, negative concentration at t=%e, after volume convection\n", min_surface, t);CHKERRQ(ierr);
#endif
          ierr = VecSetValue(my_data.c_surface_global, min_sur_index, 0.0, INSERT_VALUES);CHKERRQ(ierr);
          ierr = VecAssemblyBegin(my_data.c_surface_global);CHKERRQ(ierr);
          ierr = VecAssemblyEnd(my_data.c_surface_global);CHKERRQ(ierr);
          ierr = VecMin(my_data.c_surface_global,&min_sur_index,&min_surface); CHKERRQ(ierr);
      }
      if(min_volume<0.0 || min_surface<0.0){
          ierr = PetscPrintf(PETSC_COMM_WORLD,"DEBUG: min_vol: %e, min_sur: %e, negative concentration at t=%e, after volume convection\n", min_volume, min_surface, t);CHKERRQ(ierr);
//          ierr = OutputHandlerSingleOutput(&output_handler, t, &my_data); CHKERRQ(ierr);
//          SETERRQ(PETSC_COMM_WORLD,-1,"Negative concentrations.\n");
      }
#ifdef DEBUG
      if(t>problem_t){
          ierr = OutputHandlerSingleOutput(&output_handler, t, &my_data); CHKERRQ(ierr);
      }
#endif

#ifdef PETSC_USE_LOG
      ierr = PetscLogEventEnd(volume_convection,0,0,0,0); CHKERRQ(ierr);
#endif

#ifdef PETSC_USE_LOG
      ierr = PetscLogEventBegin(volume_diffusion,0,0,0,0); CHKERRQ(ierr);
#endif
      ierr = VolumeDiffusionSolverStep(&volume_diffusion_solver, time_step_controller.tmp_t, time_step_controller.tmp_dt, &my_data); CHKERRQ(ierr);

      ierr = VecMin(my_data.c_volume_global,&min_vol_index,&min_volume); CHKERRQ(ierr);
      ierr = VecMin(my_data.c_surface_global,&min_sur_index,&min_surface); CHKERRQ(ierr);

      while(min_volume<0.0 && min_volume > -1.0e-300){
#ifdef DEBUG
          ierr = PetscPrintf(PETSC_COMM_WORLD,"DEBUG: corrected min_vol: %e, negative concentration at t=%e, after volume diffusion\n", min_volume, t);CHKERRQ(ierr);
#endif
          ierr = VecSetValue(my_data.c_volume_global, min_vol_index, 0.0, INSERT_VALUES);CHKERRQ(ierr);
          ierr = VecAssemblyBegin(my_data.c_volume_global);CHKERRQ(ierr);
          ierr = VecAssemblyEnd(my_data.c_volume_global);CHKERRQ(ierr);
          ierr = VecMin(my_data.c_volume_global,&min_vol_index,&min_volume); CHKERRQ(ierr);
      }
      while(min_surface<0.0 && min_surface > -1.0e-300){
#ifdef DEBUG
          ierr = PetscPrintf(PETSC_COMM_WORLD,"DEBUG: corrected min_sur: %e, negative concentration at t=%e, after volume diffusion\n", min_surface, t);CHKERRQ(ierr);
#endif
          ierr = VecSetValue(my_data.c_surface_global, min_sur_index, 0.0, INSERT_VALUES);CHKERRQ(ierr);
          ierr = VecAssemblyBegin(my_data.c_surface_global);CHKERRQ(ierr);
          ierr = VecAssemblyEnd(my_data.c_surface_global);CHKERRQ(ierr);
          ierr = VecMin(my_data.c_surface_global,&min_sur_index,&min_surface); CHKERRQ(ierr);
      }
      if(min_volume<0.0 || min_surface<0.0){
          ierr = PetscPrintf(PETSC_COMM_WORLD,"DEBUG: min_vol: %e, min_sur: %e, negative concentration at t=%e, after volume diffusion\n", min_volume, min_surface, t);CHKERRQ(ierr);
//          ierr = OutputHandlerSingleOutput(&output_handler, t, &my_data); CHKERRQ(ierr);
//          SETERRQ(PETSC_COMM_WORLD,-1,"Negative concentrations.\n");
      }
#ifdef DEBUG
      if(t>problem_t){
          ierr = OutputHandlerSingleOutput(&output_handler, t, &my_data); CHKERRQ(ierr);
      }
#endif

#ifdef PETSC_USE_LOG
      ierr = PetscLogEventEnd(volume_diffusion,0,0,0,0); CHKERRQ(ierr);
#endif

#ifdef PETSC_USE_LOG
      ierr = PetscLogEventBegin(step_control,0,0,0,0); CHKERRQ(ierr);
#endif

      ierr = SingleStepSizeControl(t, dt, &time_step_controller, &my_data); CHKERRQ(ierr);

#ifdef PETSC_USE_LOG
      ierr = PetscLogEventEnd(step_control,0,0,0,0); CHKERRQ(ierr);
#endif
    }

#ifdef DEBUG
//    ierr = OutputHandlerSingleOutput(&output_handler, t, &my_data); CHKERRQ(ierr);
#endif

    t    = time_step_controller.new_t;
    dt   = time_step_controller.new_dt;
    time_step_controller.flg = -1;

    if(t >= output_handler.next_output){
      ierr = OutputHandlerSingleOutput(&output_handler, t, &my_data); CHKERRQ(ierr);
/*
      sprintf(output_filename, "concentration_volume_%d.m",(int)t);
      PetscViewerASCIIOpen(PETSC_COMM_WORLD, output_filename, &viewer);
      PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
      getSingleSpeciesVec(my_data.c_volume_global, my_data.c_volume_single, my_data.volume_species_scatter[0]);
      PetscObjectSetName((PetscObject)my_data.c_volume_single,"VII_concentration");
      VecView(my_data.c_volume_single,viewer);
      getSingleSpeciesVec(my_data.c_volume_global, my_data.c_volume_single, my_data.volume_species_scatter[3]);
      PetscObjectSetName((PetscObject)my_data.c_volume_single,"IIa_concentration");
      VecView(my_data.c_volume_single,viewer);
      getSingleSpeciesVec(my_data.c_volume_global, my_data.c_volume_single, my_data.volume_species_scatter[7]);
      PetscObjectSetName((PetscObject)my_data.c_volume_single,"II_concentration");
      VecView(my_data.c_volume_single,viewer);

      PetscViewerDestroy(&viewer);
*/
    }
  }
#ifdef DEBUG
    ierr = OutputHandlerSingleOutput(&output_handler, t, &my_data); CHKERRQ(ierr);
#endif
//  PetscViewerDestroy(&viewer);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

/*
  PetscViewer viewer;
  PetscViewerASCIIOpen(PETSC_COMM_WORLD, "coordinate_volume.m", &viewer);
  PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
  PetscObjectSetName((PetscObject)my_data.x_volume_global,"volumecoordinate");
  VecView(my_data.x_volume_global,viewer);
  PetscObjectSetName((PetscObject)my_data.da_volume_dof_1,"da");
  DMView(my_data.da_volume_dof_1,viewer);
  PetscViewerDestroy(&viewer);
*/
/*
  PetscViewerBinaryMatlabOpen(PETSC_COMM_WORLD,"coordinate_surface",&viewer);
  DMDASetFieldName(da_coordinate_surface,0,"x");
  DMDASetFieldName(da_coordinate_surface,1,"dx");
  PetscViewerBinaryMatlabOutputVecDA(viewer,"coordinatesurface",vec_coordinate_surface,da_coordinate_surface);
  PetscViewerBinaryMatlabDestroy(&viewer);
*/


  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Free work space.  All PETSc objects should be destroyed when they
     are no longer needed.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Destroy step controller, solvers and my_data
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = DestroyOutputHandler(&output_handler);CHKERRQ(ierr);
  ierr = DestroyTimeStepController(&time_step_controller);CHKERRQ(ierr);

  ierr = DestroySurfaceReactionSolver(&surface_reaction_solver);CHKERRQ(ierr);
  ierr = DestroyVolumeReactionSolver(&volume_reaction_solver);CHKERRQ(ierr);
  ierr = DestroyVolumeConvectionSolver(&volume_convection_solver);CHKERRQ(ierr);
  ierr = DestroyVolumeDiffusionSolver(&volume_diffusion_solver);CHKERRQ(ierr);

  ierr = DestroyMY_DATA(&my_data);CHKERRQ(ierr);

  ierr = PetscFinalize();

  PetscFunctionReturn(0);
}

