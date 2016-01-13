#ifndef _CUSTOM_SNES_CONVERGENCE_TEST_H_
#define _CUSTOM_SNES_CONVERGENCE_TEST_H_

PetscErrorCode SNESConvergenceTestFunction(SNES snes, PetscInt its,PetscReal xnorm, PetscReal snorm, PetscReal fnorm, SNESConvergedReason *reason, void *ctx)
{
  PetscErrorCode ierr;
  Vec            x, dx, f; //, y, dy;
  PetscInt       pos;
  PetscReal      min, max1;
//  SNES           solver_snes = (SNES) ctx;
  Mat            J, J_pre;
  MatStructure   flg = DIFFERENT_NONZERO_PATTERN;

  PetscFunctionBegin;

  if(its==0){
     *reason = SNES_CONVERGED_ITERATING;
  } else if(its>10){
     *reason = SNES_CONVERGED_ITS;
#ifdef DEBUG
     ierr = SNESGetSolution(snes, &x); CHKERRQ(ierr);
     ierr = SNESGetSolutionUpdate(snes, &dx); CHKERRQ(ierr);
     ierr = VecMaxPointwiseDivide(dx,x,&max1); CHKERRQ(ierr);
     ierr = VecMin(x,&pos,&min); CHKERRQ(ierr);
#endif
 //    *reason = SNES_DIVERGED_MAX_IT;
  } else{
     ierr = SNESGetSolution(snes, &x); CHKERRQ(ierr);
     ierr = SNESGetSolutionUpdate(snes, &dx); CHKERRQ(ierr);
     ierr = SNESGetFunction(snes, &f, NULL, NULL); CHKERRQ(ierr);
/*
     ierr = VecDuplicate(x, &y); CHKERRQ(ierr);

     ierr = VecCopy(x, y); CHKERRQ(ierr);
     ierr = VecAbs(y); CHKERRQ(ierr);
     ierr = VecDuplicate(dx, &dy); CHKERRQ(ierr);
     ierr = VecCopy(dx, dy); CHKERRQ(ierr);
     ierr = VecAbs(dy); CHKERRQ(ierr);
*/
     ierr = VecMaxPointwiseDivide(dx,x,&max1); CHKERRQ(ierr);
/*
     ierr = VecDestroy(&y);CHKERRQ(ierr);
     ierr = VecDestroy(&dy);CHKERRQ(ierr);
*/
//     ierr = VecMaxPointwiseDivide(dx,x,&max2); CHKERRQ(ierr);
     ierr = VecMin(x,&pos,&min); CHKERRQ(ierr);
//     ierr = VecMaxPointwiseDivide(f,x,&max2); CHKERRQ(ierr);
//     PetscPrintf(PETSC_COMM_WORLD,"\nmax:%f\n",max);
//     if(max1 < 1.0e-5 && max2 < 1.0e-5){
//     if(max1 < 1.0e-5 && min >= 0.0){
//        ierr = PetscPrintf(PETSC_COMM_SELF,"Inside: Number of SNES iterations = %D, xnorm: %G, snorm: %G, fnorm: %G, max:%G, min:%G\n",its,xnorm,snorm,fnorm,max1, min);CHKERRQ(ierr);
     if(min<0){
        *reason =SNES_CONVERGED_ITERATING;
        while(min<0){
#ifdef DEBUG
          ierr = PetscPrintf(PETSC_COMM_SELF,"SNESConvergenceTest: perturb %G to 0.0\n",min);CHKERRQ(ierr);
#endif
          ierr = VecSetValue(x, pos, 0.0, INSERT_VALUES);CHKERRQ(ierr);
          ierr = VecAssemblyBegin(x);CHKERRQ(ierr);
          ierr = VecAssemblyEnd(x);CHKERRQ(ierr);
          ierr = VecMin(x,&pos,&min); CHKERRQ(ierr);
        }
        ierr = SNESComputeFunction(snes,x,f);CHKERRQ(ierr);
        ierr = SNESGetJacobian(snes,&J,&J_pre,NULL,NULL);CHKERRQ(ierr);
        ierr = SNESComputeJacobian(snes,x,&J,&J_pre,&flg);CHKERRQ(ierr);
     }else if(max1 < 1.0e-3 && snorm < 1.0e-10*xnorm){
        *reason = SNES_CONVERGED_TR_DELTA;
//     } else if(snorm < 1.0e-10*xnorm && min<0.0) {
//        *reason = SNES_CONVERGED_ITERATING;
     } else {
        *reason =SNES_CONVERGED_ITERATING;
//        while(min<0){
//          ierr = VecSetValue(x, pos, 1.0e-100, INSERT_VALUES);CHKERRQ(ierr);
//          ierr = VecAssemblyBegin(x);CHKERRQ(ierr);
//          ierr = VecAssemblyEnd(x);CHKERRQ(ierr);
//          ierr = VecMin(x,&pos,&min); CHKERRQ(ierr);
//        }
     }
  }
#ifdef DEBUG
  if(its>10){
     ierr = PetscPrintf(PETSC_COMM_SELF,"Number of SNES iterations = %D, xnorm: %G, snorm: %G, fnorm: %G, max:%G, min:%G\n",its,xnorm,snorm,fnorm,max1,min);CHKERRQ(ierr);
  }
#endif
  PetscFunctionReturn(0);
}

#endif
