#ifndef _CUSTOM_SNES_MONITOR_H_
#define _CUSTOM_SNES_MONITOR_H_

PetscErrorCode CustomSNESMonitor(SNES snes, PetscInt its, PetscReal fnorm, void *ctx)
{
  PetscErrorCode ierr;
  Vec            x, dx, f;
  SNESConvergedReason reason;

  PetscFunctionBegin;
if(its>30){
  ierr = SNESGetConvergedReason(snes, &reason);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_SELF,"Number of SNES iterations = %D, reason: %s, SNES Function norm %G\n",its,SNESConvergedReasons[reason],fnorm);CHKERRQ(ierr);
  ierr = SNESGetSolution(snes,&x);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)x, "x");CHKERRQ(ierr);
  ierr = VecView(x,PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr);
  ierr = SNESGetSolutionUpdate(snes,&dx);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)dx, "dx");CHKERRQ(ierr);
  ierr = VecView(dx,PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr);
  ierr = SNESGetFunction(snes,&f,NULL,NULL);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)f, "f");CHKERRQ(ierr);
  ierr = VecView(f,PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr);
}
  PetscFunctionReturn(0);
}

#endif
