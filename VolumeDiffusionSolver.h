#ifndef _VOLUME_DIFFUSION_SOLVER_H_
#define _VOLUME_DIFFUSION_SOLVER_H_

#include "DataStructure.h"
#include "SingleSpeciesScatter.h"
#include "VolumeDiffusionSolverStruct.h"
#include "VolumeDiffusionOperator.h"

PetscErrorCode InitializeVolumeDiffusionSolver(VOLUME_DIFFUSION_SOLVER *solver_ptr, MY_DATA *user);
PetscErrorCode DestroyVolumeDiffusionSolver(VOLUME_DIFFUSION_SOLVER *solver_ptr);

PetscErrorCode VolumeDiffusionSolverStep(VOLUME_DIFFUSION_SOLVER *solver_ptr, PetscReal t, PetscReal dt, MY_DATA *user);
PetscErrorCode VolumeDiffusionSingleSpeciesRHS(TS ts, PetscReal t, Vec x, Vec f, void *ptr);
PetscErrorCode VolumeDiffusionSingleSpeciesRHSMatrix(TS ts, PetscReal t, Vec x,Mat *J_ptr,Mat *Prec_ptr,MatStructure *flag,void *ptr);

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 *     Form and set volume diffusion TS solver
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
PetscErrorCode InitializeVolumeDiffusionSolver(VOLUME_DIFFUSION_SOLVER *solver_ptr, MY_DATA *user)
{
    PetscErrorCode ierr;
    PetscInt       number_of_volume_species, number_of_grids;

    PetscFunctionBegin;

    number_of_volume_species = NUMBER_OF_SPECIES_IN_VOLUME;
    number_of_grids          = NUMBER_OF_GRIDS;

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 *     Create TS to manage hierarchical solvers
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    ierr = TSCreate(user->comm[0],&solver_ptr->ts);CHKERRQ(ierr);

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 *     Set TS properties
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    ierr = TSSetProblemType(solver_ptr->ts,TS_NONLINEAR);CHKERRQ(ierr); // linear problem
    ierr = TSSetType(solver_ptr->ts,TSBEULER);CHKERRQ(ierr); // backward euler

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 *     Create x (solution vector) and r (residual vector) 
 *     by copying the structure of user->c_volume_single
 *     both should be of size number_of_grids 
 *     to hold only species concentrations in one grid on each cpu
 *     (species concentration in volume)
 *     Set solution vector to x 
 *     (remmeber x is just a "single grid window" with no internal data)
 *     store right hand side in r
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
//    ierr = DMCreateGlobalVector(user->da_volume_dof_1,&(solver_ptr->x_global));CHKERRQ(ierr);
    ierr = VecDuplicate(user->c_volume_single, &(solver_ptr->x_global)); CHKERRQ(ierr);
    ierr = VecDuplicate(user->c_volume_single, &(solver_ptr->x_global_backup)); CHKERRQ(ierr);
//    ierr = DMCreateLocalVector(user->da_volume_dof_1,&(solver_ptr->x_local));CHKERRQ(ierr);

    ierr = TSSetSolution(solver_ptr->ts, solver_ptr->x_global); CHKERRQ(ierr);

    ierr = VecDuplicate(solver_ptr->x_global, &(solver_ptr->r_global)); CHKERRQ(ierr);
//    ierr = VecDuplicate(solver_ptr->x_local, &(solver_ptr->r_local)); CHKERRQ(ierr);

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 *  Set right-hand-side function
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    ierr = TSSetRHSFunction(solver_ptr->ts,solver_ptr->r_global,VolumeDiffusionSingleSpeciesRHS,solver_ptr);CHKERRQ(ierr);

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 *   Create matrix data structure; set Jacobian evaluation routine
 *
 *    Set Jacobian matrix data structure and default Jacobian evaluation
 *    routine. User can override with:
 *     -snes_mf : matrix-free Newton-Krylov method with no preconditioning
 *               (unless user explicitly sets preconditioner)
 *     -snes_mf_operator : form preconditioning matrix as set by the user,
 *                         but use matrix-free approx for Jacobian-vector
 *                         products within Newton-Krylov method
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    ierr = TSSetDM(solver_ptr->ts,user->da_volume_dof_1); CHKERRQ(ierr);
    ierr = DMCreateMatrix(user->da_volume_dof_1,MATAIJ,&solver_ptr->J); CHKERRQ(ierr);
/*
    ierr = MatCreate(user->comm[0], &(solver_ptr->J));CHKERRQ(ierr);
    ierr = MatSetType(solver_ptr->J, MATAIJ);CHKERRQ(ierr);
    ierr = MatSetSizes(solver_ptr->J, PETSC_DECIDE, PETSC_DECIDE, number_of_grids, number_of_grids);CHKERRQ(ierr);
    ierr = MatSetUp(solver_ptr->J);
*/
    TSSetRHSJacobian(solver_ptr->ts, solver_ptr->J, solver_ptr->J, VolumeDiffusionSingleSpeciesRHSMatrix, solver_ptr); CHKERRQ(ierr);

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 *    Set up boundary flux factors solver
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    ierr = InitializeBoundaryFluxFactorStruct(&(solver_ptr->bffs), user);CHKERRQ(ierr);

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 *    Set up volume diffusion operator   
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    ierr = FormVolumeDiffusionOperator(solver_ptr, user); CHKERRQ(ierr);

//    MatStructure J_structure;
//    VolumeDiffusionSingleSpeciesRHSMatrix(solver_ptr->ts, 0.0, solver_ptr->x_global, &solver_ptr->J, &solver_ptr->J, &J_structure, solver_ptr); CHKERRQ(ierr);
//    TSSetRHSJacobian(solver_ptr->ts, solver_ptr->J, solver_ptr->J, TSComputeRHSJacobianConstant, solver_ptr); CHKERRQ(ierr);
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 *    Set up TS
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    ierr = TSSetUp(solver_ptr->ts);CHKERRQ(ierr);

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 *    Get the KSP in TS, PC in KSP
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    ierr = TSGetSNES(solver_ptr->ts,&(solver_ptr->snes));CHKERRQ(ierr);
    ierr = SNESGetKSP(solver_ptr->snes,&(solver_ptr->ksp));CHKERRQ(ierr);
//    ierr = TSGetKSP(solver_ptr->ts,&(solver_ptr->ksp));CHKERRQ(ierr);
    ierr = KSPGetPC(solver_ptr->ksp,&(solver_ptr->pc));CHKERRQ(ierr);


/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 * Set SNES Tolerances and Monitor
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
//    ierr = SNESSetTolerances(solver_ptr->snes,0.0,0.0,0.0,50,1000); CHKERRQ(ierr);
//    ierr = SNESSetConvergenceTest(solver_ptr->snes,SNESConvergenceTestFunction,NULL,NULL); CHKERRQ(ierr);
#ifdef DEBUG
//    ierr = SNESMonitorSet(solver_ptr->snes,CustomSNESMonitor,NULL,0);CHKERRQ(ierr);
#endif

/*
 * Set KSP method to bicg, PC method to ILU
 */
//    ierr = KSPSetType(solver_ptr->ksp, KSPBICG);CHKERRQ(ierr);
//    ierr = PCSetType(solver_ptr->pc, PCILU);CHKERRQ(ierr);
    ierr = KSPSetType(solver_ptr->ksp, KSPPREONLY);CHKERRQ(ierr);
    ierr = PCSetType(solver_ptr->pc, PCLU);CHKERRQ(ierr);
    ierr = PCFactorSetMatSolverPackage(solver_ptr->pc,MATSOLVERSUPERLU_DIST); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 *     Destroy volume diffusion TS solver
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
PetscErrorCode DestroyVolumeDiffusionSolver(VOLUME_DIFFUSION_SOLVER *solver_ptr)
{
    PetscErrorCode ierr;

    PetscFunctionBegin;

    ierr = VecDestroy(&(solver_ptr->x_global)); CHKERRQ(ierr);
    ierr = VecDestroy(&(solver_ptr->x_global_backup)); CHKERRQ(ierr);
    ierr = VecDestroy(&(solver_ptr->r_global)); CHKERRQ(ierr);
//    ierr = VecDestroy(&(solver_ptr->x_local)); CHKERRQ(ierr);
//    ierr = VecDestroy(&(solver_ptr->r_local)); CHKERRQ(ierr);
    ierr = MatDestroy(&(solver_ptr->J)); CHKERRQ(ierr);
    ierr = TSDestroy(&(solver_ptr->ts)); CHKERRQ(ierr);

    ierr = DestroyBoundaryFluxFactorStruct(&(solver_ptr->bffs));CHKERRQ(ierr);
    ierr = DestroyVolumeDiffusionOperator(solver_ptr);

    PetscFunctionReturn(0);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 *     Execute one volume diffusion step for each species
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
PetscErrorCode VolumeDiffusionSolverStep(VOLUME_DIFFUSION_SOLVER *solver_ptr, PetscReal t, PetscReal dt, MY_DATA *user)
{
    PetscErrorCode   ierr;
    PetscInt         k,number_of_volume_species;
    PetscReal        flux_conversion_ratio, ratiooverdy0;

    PetscFunctionBegin;

    number_of_volume_species = NUMBER_OF_SPECIES_IN_VOLUME;
    flux_conversion_ratio    = FLUX_CONVERSION_RATIO;

/*
 *  Calculate boundary flux factors
 */
    ierr = BoundaryFluxFactorsComputation(&(solver_ptr->bffs), t, user); CHKERRQ(ierr);

/*
 *  Divide boundary flux factos alpha and beta by dy0
 *  Because for the boundary layer,
 *  dx*dy* dc/dt = Diffusion(c) + flux_conversion_ratio * Boundary_flux(c) * dx
 *  =>
 *  dc/dt = Diffusion(c)/(dx*dy) + flux_conversion_ratio * Boundary_flux(c)/dy
 *
 */
    ratiooverdy0 = flux_conversion_ratio/user->dy0;
    ierr = VecScale(solver_ptr->bffs.alpha_global, ratiooverdy0); CHKERRQ(ierr);
    ierr = VecScale(solver_ptr->bffs.beta_global, ratiooverdy0); CHKERRQ(ierr);

/*
 *  Calculate diffusion equation for each volume species and restore to the whole vector
 */
    for(k=0;k<number_of_volume_species;++k){
    /*
     *  Get single species vector
     */
        ierr = getSingleSpeciesVec(user->c_volume_global, solver_ptr->x_global, user->volume_species_scatter[k]); CHKERRQ(ierr);
#ifdef DEBUG
        ierr = VecCopy(solver_ptr->x_global, solver_ptr->x_global_backup);
#endif

    /*
     *  Set current species in solver_ptr
     */
        solver_ptr->currentSpecies = k;

    /* 
     * Solve volume diffusion time-stepping over the locally owned part of the grid
     */
        ierr = TSSetTime(solver_ptr->ts,t);CHKERRQ(ierr);
        ierr = TSSetTimeStep(solver_ptr->ts,dt);CHKERRQ(ierr);

        ierr = TSStep(solver_ptr->ts);CHKERRQ(ierr);

#ifdef DEBUG_2
        PetscReal min;
        VecMin(solver_ptr->x_global, PETSC_NULL, &min);
/*
        if(min<0.0){
          ierr = PetscObjectSetName((PetscObject)solver_ptr->x_global, "x");CHKERRQ(ierr);
          ierr = VecView(solver_ptr->x_global,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
          ierr = PetscObjectSetName((PetscObject)solver_ptr->x_global_backup, "prev_x");CHKERRQ(ierr);
          ierr = VecView(solver_ptr->x_global_backup,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
          ierr = PetscObjectSetName((PetscObject)solver_ptr->J, "J");CHKERRQ(ierr);
          ierr = MatView(solver_ptr->J,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
          ierr = PetscObjectSetName((PetscObject)solver_ptr->bffs.alpha_global, "alpha");CHKERRQ(ierr);
          ierr = VecView(solver_ptr->bffs.alpha_global,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
          ierr = PetscObjectSetName((PetscObject)solver_ptr->bffs.beta_global, "beta");CHKERRQ(ierr);
          ierr = VecView(solver_ptr->bffs.beta_global,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
        }
*/
        PetscReal **array_c_volume, **array_prev_c_volume;
        Species_Volume *array_alpha_global, *array_beta_global;
        PetscInt xs, xe, xm, i;
        PetscReal    temp_c, prev_temp_c;

        ierr = DMDAVecGetArray(user->da_volume_dof_1, solver_ptr->x_global, &array_c_volume);CHKERRQ(ierr);
        ierr = DMDAVecGetArray(user->da_volume_dof_1, solver_ptr->x_global_backup, &array_prev_c_volume);CHKERRQ(ierr);
        ierr = DMDAVecGetArray(solver_ptr->bffs.da_alpha_beta, solver_ptr->bffs.alpha_global, &array_alpha_global);CHKERRQ(ierr);
        ierr = DMDAVecGetArray(solver_ptr->bffs.da_alpha_beta, solver_ptr->bffs.beta_global, &array_beta_global);CHKERRQ(ierr);

        ierr = DMDAGetCorners(user->da_c_surface, &xs, PETSC_NULL, PETSC_NULL, &xm, PETSC_NULL, PETSC_NULL);CHKERRQ(ierr);

        xe = xs + xm;

        for(i=xs;i<xe;++i){
          temp_c = solver_ptr->bffs.ratio_layer0 * array_c_volume[0][i] + solver_ptr->bffs.ratio_layer1 * array_c_volume[1][i];
          prev_temp_c = solver_ptr->bffs.ratio_layer0 * array_prev_c_volume[0][i] + solver_ptr->bffs.ratio_layer1 * array_prev_c_volume[1][i];
          if(temp_c < 0.0 && (array_alpha_global[i].c[k]!=0.0 || array_beta_global[i].c[k]!=0.0)){
//          if(temp_c < 0.0){

//          if(i==29 && k == 3){
              PetscPrintf(PETSC_COMM_SELF,"ah ha, t=%G, i=%d, k=%d, temp_c=%G, prev_temp_c=%G, ratio0=%G, ratio1=%G, alpha=%G, beta=%G\n",t,i,k,temp_c,prev_temp_c, solver_ptr->bffs.ratio_layer0, solver_ptr->bffs.ratio_layer1, array_alpha_global[i].c[k],array_beta_global[i].c[k]);
              PetscPrintf(PETSC_COMM_SELF,"prev:\n");
              PetscPrintf(PETSC_COMM_SELF,"%G, %G, %G\n", array_prev_c_volume[2][i-1], array_prev_c_volume[2][i],array_prev_c_volume[2][i+1]);
              PetscPrintf(PETSC_COMM_SELF,"%G, %G, %G\n", array_prev_c_volume[1][i-1], array_prev_c_volume[1][i],array_prev_c_volume[1][i+1]);
              PetscPrintf(PETSC_COMM_SELF,"%G, %G, %G\n", array_prev_c_volume[0][i-1], array_prev_c_volume[0][i],array_prev_c_volume[0][i+1]);
              PetscPrintf(PETSC_COMM_SELF,"new:\n");
              PetscPrintf(PETSC_COMM_SELF,"%G, %G, %G\n", array_c_volume[2][i-1], array_c_volume[2][i],array_c_volume[2][i+1]);
              PetscPrintf(PETSC_COMM_SELF,"%G, %G, %G\n", array_c_volume[1][i-1], array_c_volume[1][i],array_c_volume[1][i+1]);
              PetscPrintf(PETSC_COMM_SELF,"%G, %G, %G\n", array_c_volume[0][i-1], array_c_volume[0][i],array_c_volume[0][i+1]);
          }
        }

        ierr = DMDAVecRestoreArray(solver_ptr->bffs.da_alpha_beta, solver_ptr->bffs.alpha_global, &array_alpha_global);CHKERRQ(ierr);
        ierr = DMDAVecRestoreArray(solver_ptr->bffs.da_alpha_beta, solver_ptr->bffs.beta_global, &array_beta_global);CHKERRQ(ierr);
        ierr = DMDAVecRestoreArray(user->da_volume_dof_1, solver_ptr->x_global, &array_c_volume);CHKERRQ(ierr);
        ierr = DMDAVecRestoreArray(user->da_volume_dof_1, solver_ptr->x_global_backup, &array_prev_c_volume);CHKERRQ(ierr);
#endif

    /*
     *  Restore single species vector to the global vector
     */
        ierr = restoreSingleSpeciesVec(user->c_volume_global, solver_ptr->x_global, user->volume_species_scatter[k]); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}

PetscErrorCode VolumeDiffusionSingleSpeciesRHSMatrix(TS ts, PetscReal t,Vec x,Mat *J_ptr,Mat *Prec_ptr,MatStructure *flag,void *ptr)
{
   VOLUME_DIFFUSION_SOLVER  *solver_ptr = (VOLUME_DIFFUSION_SOLVER *)ptr;     /* user-defined application context */
   DM                       da;
   DMDALocalInfo            info;
   PetscInt                 i,j,k;
   PetscReal                val[6];
   PetscInt                 nc;
   MatStencil               row,col[6];
   PetscReal                ratio_layer0, ratio_layer1;
   PetscErrorCode           ierr;
   IMPLICIT_WEIGHTS         **array_implicit_weight;
   Species_Volume           *array_beta;

   PetscFunctionBegin;

   // get ratio_layer0 and ratio_layer1 for boundary flux for boundary layer
   ratio_layer0 = solver_ptr->bffs.ratio_layer0;
   ratio_layer1 = solver_ptr->bffs.ratio_layer1;

   // get array for diffusion operator vectors
   ierr = DMDAVecGetArray(solver_ptr->da_diffusion_operator, solver_ptr->diffusion_operator, &array_implicit_weight);CHKERRQ(ierr);
 
//   PetscPrintf(PETSC_COMM_WORLD,"diffusion operator:\n");
//   VecView(solver_ptr->diffusion_operator, PETSC_VIEWER_STDOUT_WORLD);

   // get array of beta from boundary flux
   ierr = DMDAVecGetArray(solver_ptr->bffs.da_alpha_beta, solver_ptr->bffs.beta_global, &array_beta);CHKERRQ(ierr);

//   PetscPrintf(PETSC_COMM_WORLD,"boundary flux beta:\n");
//   VecView(solver_ptr->bffs.beta_global, PETSC_VIEWER_STDOUT_WORLD);

   // get currentSpecies
   k = solver_ptr->currentSpecies;

   ierr = TSGetDM(ts,&da); CHKERRQ(ierr);
   DMDAGetLocalInfo(da,&info);
   for (j=info.ys; j<info.ys+info.ym; j++) {
     for (i=info.xs; i<info.xs+info.xm; i++) {
       nc = 0;
       row.i = i; row.j = j;
       
   /*
    * remember for boundary layer value we need to  include flux
    * remember c_boundary[i] = ratio_layer0 * c[0][i] + ratio_layer[1] * c[1][i]
    * therefore for grid (0,j), 
    * value of stencil (0,j) += ratio_layer0 * beta (beta is already divided by dy)
    * value of stencil (1,j) += ratio_layer1 * beta (beta is already divided by dy)
    */
       // center stencil 
       if(j != 0){
         col[nc].i = i; col[nc].j = j;   val[nc] = array_implicit_weight[j][i].iw0;
         nc++;
       } else {
         col[nc].i = i; col[nc].j = j;   val[nc] = array_implicit_weight[j][i].iw0 + array_beta[i].c[k] * ratio_layer0;
//         col[nc].i = i; col[nc].j = j;   val[nc] = array_implicit_weight[j][i].iw0;
         nc++;
       }

       // left stencil
       if(i != 0){
         col[nc].i = i-1; col[nc].j = j;   val[nc] = array_implicit_weight[j][i].iwleft1;
         nc++;
       }

       // right stencil
       if(i != info.mx-1){
         col[nc].i = i+1; col[nc].j = j;   val[nc] = array_implicit_weight[j][i].iwright1;
         nc++;
       }

       // down1 stencil
       if(j!=0){
         col[nc].i = i;   col[nc].j = j-1; val[nc] = array_implicit_weight[j][i].iwdown1;
         nc++;
       }

       // down2 stencil
       if(j>1){
         col[nc].i = i;   col[nc].j = j-2; val[nc] = array_implicit_weight[j][i].iwdown2;
         nc++;
       }

       // up stencil
       if (j == 0){
         col[nc].i = i;   col[nc].j = j+1; val[nc] = array_implicit_weight[j][i].iwup1 + array_beta[i].c[k] * ratio_layer1;
//         col[nc].i = i;   col[nc].j = j+1; val[nc] = array_implicit_weight[j][i].iwup1;
         nc++;
       } else if(j != info.my-1){
         col[nc].i = i;   col[nc].j = j+1; val[nc] = array_implicit_weight[j][i].iwup1;
         nc++;
       }

       MatSetValuesStencil(*Prec_ptr,1,&row,nc,col,val,INSERT_VALUES);
     }
   }


   /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Complete the matrix assembly process and set some options
      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
   /*
      Assemble matrix, using the 2-step process:
        MatAssemblyBegin(), MatAssemblyEnd()
      Computations can be done while messages are in transition
      by placing code between these two statements.
   */
   MatAssemblyBegin(*Prec_ptr,MAT_FINAL_ASSEMBLY);
   MatAssemblyEnd(*Prec_ptr,MAT_FINAL_ASSEMBLY);
   if (*J_ptr != *Prec_ptr) {
     MatAssemblyBegin(*J_ptr,MAT_FINAL_ASSEMBLY);
     MatAssemblyEnd(*J_ptr,MAT_FINAL_ASSEMBLY);
   }

   /*
      Set flag to indicate that the Jacobian matrix retains an identical
      nonzero structure throughout all timestepping iterations (although the
      values of the entries change). Thus, we can save some work in setting
      up the preconditioner (e.g., no need to redo symbolic factorization for
      ILU/ICC preconditioners).
       - If the nonzero structure of the matrix is different during
         successive linear solves, then the flag DIFFERENT_NONZERO_PATTERN
         must be used instead.  If you are unsure whether the matrix
         structure has changed or not, use the flag DIFFERENT_NONZERO_PATTERN.
       - Caution:  If you specify SAME_NONZERO_PATTERN, PETSc
         believes your assertion and does not check the structure
         of the matrix.  If you erroneously claim that the structure
         is the same when it actually is not, the new preconditioner
         will not function correctly.  Thus, use this optimization
         feature with caution!
   */
   *flag = SAME_NONZERO_PATTERN;

   /*
      Set and option to indicate that we will never add a new nonzero location 
      to the matrix. If we do, it will generate an error.
   */
   ierr = MatSetOption(*Prec_ptr,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE); CHKERRQ(ierr);

   // restore array for diffusion operator vectors
   ierr = DMDAVecRestoreArray(solver_ptr->da_diffusion_operator, solver_ptr->diffusion_operator, &array_implicit_weight);CHKERRQ(ierr);

   // get array of beta from boundary flux
   ierr = DMDAVecRestoreArray(solver_ptr->bffs.da_alpha_beta, solver_ptr->bffs.beta_global, &array_beta);CHKERRQ(ierr);

//   PetscPrintf(PETSC_COMM_WORLD,"Matrix:\n");
//   MatView(*Prec_ptr,PETSC_VIEWER_STDOUT_WORLD);

   PetscFunctionReturn(0);
}

PetscErrorCode VolumeDiffusionSingleSpeciesRHS(TS ts, PetscReal t, Vec x, Vec f, void *ptr)
{
   VOLUME_DIFFUSION_SOLVER  *solver_ptr = (VOLUME_DIFFUSION_SOLVER *)ptr;     /* user-defined application context */
   DM                       da;
   Mat                      A;
   MatStructure             A_structure;
   PetscErrorCode           ierr;
   PetscInt                 xs, xm, xe;
   PetscReal                **array_f;
   PetscInt                 i, k;
   Species_Volume           *array_alpha;

   PetscFunctionBegin;

   ierr = TSGetRHSJacobian(ts, &A, PETSC_NULL, PETSC_NULL, &ptr); CHKERRQ(ierr);
   ierr = VolumeDiffusionSingleSpeciesRHSMatrix(ts, t, x, &A, &A, &A_structure, solver_ptr); CHKERRQ(ierr);
#ifdef DEBUG
//   MatView(A, PETSC_VIEWER_STDOUT_WORLD);
#endif

   ierr = MatMult(A,x,f); CHKERRQ(ierr);

//   PetscPrintf(PETSC_COMM_WORLD,"boundary flux alpha:\n");
//   VecView(solver_ptr->bffs.alpha_global, PETSC_VIEWER_STDOUT_WORLD);

   // get array of f
   ierr = TSGetDM(ts,&da); CHKERRQ(ierr);
   ierr = DMDAVecGetArray(da, f, &array_f);CHKERRQ(ierr);
   // get array of alpha
   ierr = DMDAVecGetArray(solver_ptr->bffs.da_alpha_beta, solver_ptr->bffs.alpha_global, &array_alpha);CHKERRQ(ierr);

/*
 *  Get local grid boundaries (for 1-dimensional DMDA):
 *   xs   - starting grid indices (no ghost points)
 *   xm   - widths of local grid (no ghost points)
 */
   ierr = DMDAGetCorners(solver_ptr->bffs.da_alpha_beta, &xs, PETSC_NULL, PETSC_NULL, &xm, PETSC_NULL, PETSC_NULL);CHKERRQ(ierr);

   xe = xs + xm;

   k = solver_ptr->currentSpecies;

   for (i=xs; i<xe; i++) {
       array_f[0][i] += array_alpha[i].c[k];
   }

   // restore array of f and alpha
   ierr = DMDAVecRestoreArray(da, f, &array_f);CHKERRQ(ierr);
   ierr = DMDAVecRestoreArray(solver_ptr->bffs.da_alpha_beta, solver_ptr->bffs.alpha_global, &array_alpha);CHKERRQ(ierr);

//   PetscPrintf(PETSC_COMM_WORLD,"x:\n");
//   VecView(x, PETSC_VIEWER_STDOUT_WORLD);
//   PetscPrintf(PETSC_COMM_WORLD,"f:\n");
//   VecView(f, PETSC_VIEWER_STDOUT_WORLD);

   PetscFunctionReturn(0);
}

#endif
