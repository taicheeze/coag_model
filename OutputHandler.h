#ifndef _OUTPUT_HANDLER_H_
#define _OUTPUT_HANDLER_H_

#include "DataStructure.h"
#include "string.h"

typedef struct {
   PetscBool     volume_species_output_flag[NUMBER_OF_SPECIES_IN_VOLUME];   /* true: output, false: donot output */
   PetscBool     surface_species_output_flag[NUMBER_OF_SPECIES_ON_SURFACE];   /* true: output, false: donot output */
   PetscReal     output_dt;            /* time interval for consecutive outputs */
   PetscReal     next_output;          /* time for next output */
   PetscInt      output_counter;       /* how many output files have been generated */
} OUTPUT_HANDLER;


PetscErrorCode InitializeOutputHandler(OUTPUT_HANDLER *output_ptr, PetscReal output_dt, const char *filename);
PetscErrorCode DestroyOutputHandler(OUTPUT_HANDLER *output_ptr);

PetscErrorCode OutputHandlerOtherOutput(MY_DATA *user);
PetscErrorCode OutputHandlerSingleOutput(OUTPUT_HANDLER *output_ptr, PetscReal t, MY_DATA *user);

PetscErrorCode InitializeOutputHandler(OUTPUT_HANDLER *output_ptr, PetscReal output_dt, const char *filename)
{
    PetscErrorCode ierr;
    PetscInt       i, number_of_volume_species, number_of_surface_species;
    PetscInt       j, number_of_volume_output_species, number_of_surface_output_species;
    FILE           *fp;
    char           species_name[256];
//    PetscMPIInt    rank;
//    bool           notendoffile;

    PetscFunctionBegin;

    number_of_volume_species = NUMBER_OF_SPECIES_IN_VOLUME;
    number_of_surface_species = NUMBER_OF_SPECIES_ON_SURFACE;

    output_ptr->output_dt = output_dt;
    output_ptr->next_output = 0.0;

    for(i=0;i<number_of_volume_species;++i){
         output_ptr->volume_species_output_flag[i] = PETSC_FALSE;
    } 
    for(i=0;i<number_of_surface_species;++i){
         output_ptr->surface_species_output_flag[i] = PETSC_FALSE;
    } 

//    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);


    ierr = PetscFOpen(PETSC_COMM_WORLD,filename,"r",&fp); CHKERRQ(ierr);

//    notendoffile = true;
    ierr = PetscSynchronizedFGets(PETSC_COMM_WORLD,fp,256,species_name); CHKERRQ(ierr);
    number_of_volume_output_species = atoi(species_name);

    j=0;
//    while(notendoffile == true){
    while(j<number_of_volume_output_species){
//       strcpy(species_name,"NULL");
//       PetscPrintf(PETSC_COMM_WORLD,"%d,%s\n", strlen(species_name), species_name);
       ierr = PetscSynchronizedFGets(PETSC_COMM_WORLD,fp,256,species_name); CHKERRQ(ierr);
//       PetscPrintf(PETSC_COMM_WORLD,"%d,%s\n", strlen(species_name), species_name);
//       if(strcmp(species_name,"NULL") == 0){
//           PetscPrintf(PETSC_COMM_WORLD,"hello\n");
//           notendoffile = false;
//       } else {
           species_name[strlen(species_name)-1] = '\0'; // remove '\n' at the end
           for(i=0;i<number_of_volume_species;++i){
               if(strcmp(species_name,VOLUME_SPECIES_NAMES[i])==0){
                   output_ptr->volume_species_output_flag[i] = PETSC_TRUE;
               }
           }
//       }
       j++;
    }

//    notendoffile = true;
    ierr = PetscSynchronizedFGets(PETSC_COMM_WORLD,fp,256,species_name); CHKERRQ(ierr);
    number_of_surface_output_species = atoi(species_name);

    j=0;
//    while(notendoffile == true){
    while(j<number_of_surface_output_species){
//       strcpy(species_name,"NULL");
//       PetscPrintf(PETSC_COMM_WORLD,"%d,%s\n", strlen(species_name), species_name);
       ierr = PetscSynchronizedFGets(PETSC_COMM_WORLD,fp,256,species_name); CHKERRQ(ierr);
//       PetscPrintf(PETSC_COMM_WORLD,"%d,%s\n", strlen(species_name), species_name);
//       if(strcmp(species_name,"NULL") == 0){
//           PetscPrintf(PETSC_COMM_WORLD,"hello\n");
//           notendoffile = false;
//       } else {
           species_name[strlen(species_name)-1] = '\0'; // remove '\n' at the end
           for(i=0;i<number_of_surface_species;++i){
               if(strcmp(species_name,SURFACE_SPECIES_NAMES[i])==0){
                   output_ptr->surface_species_output_flag[i] = PETSC_TRUE;
               }
           }
//       }
       j++;
    }

    ierr = PetscFClose(PETSC_COMM_WORLD,fp); CHKERRQ(ierr);

    output_ptr->output_counter = 0;

    PetscFunctionReturn(0);

}

PetscErrorCode DestroyOutputHandler(OUTPUT_HANDLER *output_ptr)
{
    PetscFunctionBegin;

    PetscFunctionReturn(0);

}

PetscErrorCode OutputHandlerOtherOutput(MY_DATA *user)
{
    PetscErrorCode ierr;
    PetscViewer    viewer;
    char           output_filename[50];

    PetscFunctionBegin;

    sprintf(output_filename, "output/coordinates.m");

    ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, output_filename, &viewer); CHKERRQ(ierr);
    ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB); CHKERRQ(ierr);

    ierr = PetscViewerASCIIPrintf(viewer,"nx=[\n%d\n];\n",NUMBER_OF_X_GRIDS); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"ny=[\n%d\n];\n",NUMBER_OF_Y_GRIDS); CHKERRQ(ierr);

    ierr = PetscObjectSetName((PetscObject)user->x_volume_global,"x_volume_global"); CHKERRQ(ierr);
    ierr = VecView(user->x_volume_global,viewer); CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject)user->dx_volume_global,"dx_volume_global"); CHKERRQ(ierr);
    ierr = VecView(user->dx_volume_global,viewer); CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject)user->y_volume_global,"y_volume_global"); CHKERRQ(ierr);
    ierr = VecView(user->y_volume_global,viewer); CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject)user->dy_volume_global,"dy_volume_global"); CHKERRQ(ierr);
    ierr = VecView(user->dy_volume_global,viewer); CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);


    sprintf(output_filename, "output/velocity_coordinates.m");

    ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF, output_filename, &viewer); CHKERRQ(ierr);
    ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB); CHKERRQ(ierr);

    ierr = PetscViewerASCIIPrintf(viewer,"nx=[\n%d\n];\n",NUMBER_OF_X_GRIDS); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"ny=[\n%d\n];\n",NUMBER_OF_Y_GRIDS); CHKERRQ(ierr);

    ierr = PetscObjectSetName((PetscObject)user->x_vx,"x_vx"); CHKERRQ(ierr);
    ierr = VecView(user->x_vx,viewer); CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject)user->y_vx,"y_vx"); CHKERRQ(ierr);
    ierr = VecView(user->y_vx,viewer); CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject)user->x_vy,"x_vy"); CHKERRQ(ierr);
    ierr = VecView(user->x_vy,viewer); CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject)user->y_vy,"y_vy"); CHKERRQ(ierr);
    ierr = VecView(user->y_vy,viewer); CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

PetscErrorCode OutputHandlerSingleOutput(OUTPUT_HANDLER *output_ptr, PetscReal t, MY_DATA *user)
{
    PetscErrorCode ierr;
    PetscInt       i,number_of_volume_species,number_of_surface_species;
    PetscViewer    viewer;
    char           output_filename[50]="";

    PetscFunctionBegin;

    number_of_volume_species = NUMBER_OF_SPECIES_IN_VOLUME;
    number_of_surface_species = NUMBER_OF_SPECIES_ON_SURFACE;

    sprintf(output_filename, "output/concentrations_%d.m",output_ptr->output_counter);

    ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, output_filename, &viewer); CHKERRQ(ierr);
    ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB); CHKERRQ(ierr);

    ierr = PetscViewerASCIIPrintf(viewer,"t=[\n%e\n];\n",t); CHKERRQ(ierr);

    for(i=0;i<number_of_volume_species;++i){
#ifndef DEBUG
        if(output_ptr->volume_species_output_flag[i] == PETSC_TRUE){ //comment out to output everything
#endif
            ierr = getSingleSpeciesVec(user->c_volume_global, user->c_volume_single, user->volume_species_scatter[i]); CHKERRQ(ierr);
            ierr = PetscObjectSetName((PetscObject)user->c_volume_single,VOLUME_SPECIES_NAMES[i]); CHKERRQ(ierr);
            ierr = VecView(user->c_volume_single,viewer); CHKERRQ(ierr);
#ifndef DEBUG
        } //comment out to output everything
#endif
    }

    for(i=0;i<number_of_surface_species;++i){
#ifndef DEBUG
         if(output_ptr->surface_species_output_flag[i] == PETSC_TRUE){ //comment out to output everything
#endif
            ierr = getSingleSpeciesVec(user->c_surface_global, user->c_surface_single, user->surface_species_scatter[i]); CHKERRQ(ierr);
            ierr = PetscObjectSetName((PetscObject)user->c_surface_single,SURFACE_SPECIES_NAMES[i]); CHKERRQ(ierr);
            ierr = VecView(user->c_surface_single,viewer); CHKERRQ(ierr);
#ifndef DEBUG
         } //comment out to output everything
#endif
    }

    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

    output_ptr->output_counter++;
    output_ptr->next_output += output_ptr->output_dt;

    PetscFunctionReturn(0);
}

#endif
