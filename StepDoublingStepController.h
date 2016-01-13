#ifndef _STEP_DOUBLING_STEP_CONTROLLER_H_
#define _STEP_DOUBLING_STEP_CONTROLLER_H_

#include "DataStructure.h"

typedef struct {
    Vec         c_volume_single_step, c_surface_single_step; // result after single step, also to store errors
    Vec         c_volume_tmp, c_surface_tmp;  // result before the step in case the step fail and need reverting
    Vec         c_volume_err_tmp, c_surface_err_tmp;  // for comparison purpose
    Vec         c_volume_err_tol, c_surface_err_tol;  // error tolerance
    PetscReal   tmp_t, tmp_dt; // saved current_t and current_dt 
    PetscReal   new_t, new_dt; // intemediate t and dt, named new_t and new_dt fot consistency with other step controller
    PetscInt    flg; // flg to indicate which stage we are in
                     // flg =-1: reset by outer loop, before single step
                     // flg = 0: accept step after second of half step
                     // flg = 1: after single step
                     // flg = 2: after first of half step
                     // flg = 3: deny step after second of half step
} TIME_STEP_CONTROLLER; 

PetscErrorCode InitializeTimeStepController(PetscReal t, PetscReal dt, TIME_STEP_CONTROLLER *controller, MY_DATA *user);
PetscErrorCode SingleStepSizeControl(PetscReal t, PetscReal dt, TIME_STEP_CONTROLLER *controller, MY_DATA *user);
PetscErrorCode DestroyTimeStepController(TIME_STEP_CONTROLLER *controller);


PetscErrorCode InitializeTimeStepController(PetscReal t, PetscReal dt, TIME_STEP_CONTROLLER *controller, MY_DATA *user)
{
    PetscErrorCode ierr;
    PetscReal      atol, rtol;

    PetscFunctionBegin;

    atol = ATOL;
    rtol = RTOL;

    controller->tmp_t  = t;
    controller->tmp_dt = dt;
    controller->new_t  = t;
    controller->new_dt = dt;

    ierr = VecDuplicate(user->c_volume_global, &(controller->c_volume_single_step)); CHKERRQ(ierr);
    ierr = VecDuplicate(user->c_volume_global, &(controller->c_volume_tmp)); CHKERRQ(ierr);
    ierr = VecDuplicate(user->c_volume_global, &(controller->c_volume_err_tmp)); CHKERRQ(ierr);
    ierr = VecDuplicate(user->c_volume_global, &(controller->c_volume_err_tol)); CHKERRQ(ierr);

    ierr = VecDuplicate(user->c_surface_global, &(controller->c_surface_single_step)); CHKERRQ(ierr);
    ierr = VecDuplicate(user->c_surface_global, &(controller->c_surface_tmp)); CHKERRQ(ierr);
    ierr = VecDuplicate(user->c_surface_global, &(controller->c_surface_err_tmp)); CHKERRQ(ierr);
    ierr = VecDuplicate(user->c_surface_global, &(controller->c_surface_err_tol)); CHKERRQ(ierr);

    // copy into tmp in case step is rejected 
    ierr = VecCopy(user->c_volume_global, controller->c_volume_tmp); CHKERRQ(ierr);
    ierr = VecCopy(user->c_surface_global, controller->c_surface_tmp); CHKERRQ(ierr);

    // calculate error tolerance
    ierr = VecCopy(user->c_volume_global, controller->c_volume_err_tol); CHKERRQ(ierr);
    ierr = VecAbs(controller->c_volume_err_tol); CHKERRQ(ierr);
    ierr = VecScale(controller->c_volume_err_tol, rtol); CHKERRQ(ierr);
    ierr = VecShift(controller->c_volume_err_tol, atol); CHKERRQ(ierr);

    ierr = VecCopy(user->c_surface_global, controller->c_surface_err_tol); CHKERRQ(ierr);
    ierr = VecAbs(controller->c_surface_err_tol); CHKERRQ(ierr);
    ierr = VecScale(controller->c_surface_err_tol, rtol); CHKERRQ(ierr);
    ierr = VecShift(controller->c_surface_err_tol, atol); CHKERRQ(ierr);
   

    controller->flg = -1;

    PetscFunctionReturn(0);

}

/*
 *
 */
PetscErrorCode SingleStepSizeControl(PetscReal t, PetscReal dt, TIME_STEP_CONTROLLER *controller, MY_DATA *user)
{
    PetscErrorCode ierr;
    PetscReal      dt_increase_ratio, dt_decrease_ratio;
    PetscReal      c_relative_change_min, c_relative_change_max;
    PetscReal      dt_min, dt_max;
    PetscReal      atol, rtol;
    PetscReal      volume_err, surface_err;
    PetscReal      ratio, frac;
    PetscReal      volume_h, surface_h, new_h;

    PetscFunctionBegin;

    dt_min  = DT_MIN;
    dt_max  = DT_MAX;
    dt_increase_ratio = TIMESTEP_INCREASE_RATIO;
    dt_decrease_ratio = TIMESTEP_DECREASE_RATIO;
    c_relative_change_min = C_RELATIVE_TOL_MIN;
    c_relative_change_max = C_RELATIVE_TOL_MAX;
    atol = ATOL;
    rtol = RTOL;
    frac = SAFETY_FRACTION;

    /* flg to indicate which stage we are in, note it's different than 
     * in the main program as a operator splitting calculation has already been done at this point
     *   flg =-1: reset by outer loop, after single step
     *             tmp_t remains tmp_t (old_t), tmp_dt = tmp_dt/2 (h/2)
     *             copy c_global to c_single_step
     *             copy c_tmp to c_global
     *             set flg=1
     *   flg = 0: accept step after second of half step
     *            won't happen in this function
     *   flg = 1: after first of half step
     *             tmp_t = tmp_t + tmp_dt (old_t + h/2), tmp_dt remains tmp_dt (h/2)
     *             set flg=2
     *   flg = 2: after second of half step
     *             calculate error = |c_global - c_single_step|, store in c_single_step
     *             compare with error tolarance
     *             1: accept step: new_t  = tmp_t = tmp_t + tmp_dt (old_t + h)
     *                             calculate new_dt and set tmp_dt = new_dt
     *                             calculate new error tolerance
     *                             set flg=0
     *             2. reject step: tmp_t = tmp_t - tmp_dt (old_t)
     *                             calculate new tmp_dt
     *                             copy c_tmp to c_global
     *                             set flg = 3
     *   flg = 3: reset by inner loop, after single step
     *             the same as -1
     */
    if(controller->flg==-1 || controller->flg==3){
        // tmp_t remains tmp_t (old_t), tmp_dt = tmp_dt/2 (h/2)
        controller->tmp_dt /= 2.0;

        // copy c_global to c_single_step
        ierr = VecCopy(user->c_volume_global, controller->c_volume_single_step); CHKERRQ(ierr);
        ierr = VecCopy(user->c_surface_global, controller->c_surface_single_step); CHKERRQ(ierr);
        
        // copy c_tmp to c_global
        ierr = VecCopy(controller->c_volume_tmp, user->c_volume_global); CHKERRQ(ierr);
        ierr = VecCopy(controller->c_surface_tmp, user->c_surface_global); CHKERRQ(ierr);

        // set flg=1
        controller->flg = 1;

    } else if(controller->flg==1){
        // tmp_t = tmp_t + tmp_dt (old_t + h/2), tmp_dt remains tmp_dt (h/2)
        controller->tmp_t += controller->tmp_dt;

        // set flg=2
        controller->flg = 2;

    } else if(controller->flg==2){
        // calculate error 2|l_n| = |c_global - c_single_step|, store in c_single_step
        ierr = VecAXPY(controller->c_volume_single_step, -1.0, user->c_volume_global); CHKERRQ(ierr);
        ierr = VecAXPY(controller->c_surface_single_step, -1.0, user->c_surface_global); CHKERRQ(ierr);

        ierr = VecAbs(controller->c_volume_single_step); CHKERRQ(ierr);
        ierr = VecAbs(controller->c_surface_single_step); CHKERRQ(ierr);

        // compare |l_n| with error tolarance
        //  err_tmp = err_tol - err
        //  compare min(err_tmp) with 0
        ierr = VecWAXPY(controller->c_volume_err_tmp, -0.5, controller->c_volume_single_step, controller->c_volume_err_tol); CHKERRQ(ierr);
        ierr = VecWAXPY(controller->c_surface_err_tmp, -0.5, controller->c_surface_single_step, controller->c_surface_err_tol); CHKERRQ(ierr);

        ierr = VecMin(controller->c_volume_err_tmp, PETSC_NULL, &volume_err); CHKERRQ(ierr);
        ierr = VecMin(controller->c_surface_err_tmp, PETSC_NULL, &surface_err); CHKERRQ(ierr);

        if(volume_err >= 0.0 && surface_err >=0.0){
            ierr = PetscPrintf(PETSC_COMM_WORLD,"step accepted:%f\n",controller->tmp_dt);CHKERRQ(ierr);
            // 1: accept step: new_t  = tmp_t = tmp_t + tmp_dt (old_t + h)
            //     calculate new_dt and set tmp_dt = new_dt
            //     copy c_global to c_tmp
            //     calculate new error tolerance
            //     set flg=0
            controller->tmp_t += controller->tmp_dt;
            controller->new_t = controller->tmp_t;

            // (new_h/old_h)^2 = frac * err_tol / err
            // tmp_dt = h/2
            ratio = 8.0 * frac * controller->tmp_dt * controller->tmp_dt;

            ierr = VecPointwiseDivide(controller->c_volume_err_tmp, controller->c_volume_err_tol, controller->c_volume_single_step); CHKERRQ(ierr);
            ierr = VecPointwiseDivide(controller->c_surface_err_tmp, controller->c_surface_err_tol, controller->c_surface_single_step); CHKERRQ(ierr);

            ierr = VecScale(controller->c_volume_err_tmp, ratio);
            ierr = VecScale(controller->c_surface_err_tmp, ratio);

            // now it's actually new_h^2
            ierr = VecMin(controller->c_volume_err_tmp, PETSC_NULL, &volume_h); CHKERRQ(ierr);
            ierr = VecMin(controller->c_surface_err_tmp, PETSC_NULL, &surface_h); CHKERRQ(ierr);

            new_h = (volume_h<surface_h)?volume_h:surface_h;

            // now it's h, set new_dt and tmp_dt
            controller->tmp_dt = sqrt(new_h);
            controller->new_dt = controller->tmp_dt;
/*
            if(controller->tmp_dt < dt_max){
                 controller->tmp_dt *= (2.0*dt_increase_ratio);
            }
            controller->new_dt = controller->tmp_dt;
*/
            ierr = VecCopy(user->c_volume_global, controller->c_volume_err_tol); CHKERRQ(ierr);
            ierr = VecAbs(controller->c_volume_err_tol); CHKERRQ(ierr);
            ierr = VecScale(controller->c_volume_err_tol, rtol); CHKERRQ(ierr);
            ierr = VecShift(controller->c_volume_err_tol, atol); CHKERRQ(ierr);

            ierr = VecCopy(user->c_surface_global, controller->c_surface_err_tol); CHKERRQ(ierr);
            ierr = VecAbs(controller->c_surface_err_tol); CHKERRQ(ierr);
            ierr = VecScale(controller->c_surface_err_tol, rtol); CHKERRQ(ierr);
            ierr = VecShift(controller->c_surface_err_tol, atol); CHKERRQ(ierr);

            ierr = VecCopy(user->c_volume_global, controller->c_volume_tmp); CHKERRQ(ierr);
            ierr = VecCopy(user->c_surface_global, controller->c_surface_tmp); CHKERRQ(ierr);

            controller->flg = 0;
        } else {
            ierr = PetscPrintf(PETSC_COMM_WORLD,"step rejected:%f\n",controller->tmp_dt);CHKERRQ(ierr);
            // 2. reject step: tmp_t = tmp_t - tmp_dt (old_t)
            //     calculate new tmp_dt
            //     copy c_tmp to c_global
            //     set flg = 3
            controller->tmp_t -= controller->tmp_dt;

            // (new_h/old_h)^2 = frac * err_tol / err
            // tmp_dt = h/2
            ratio = 8.0 * frac * controller->tmp_dt * controller->tmp_dt;

            ierr = VecPointwiseDivide(controller->c_volume_err_tmp, controller->c_volume_err_tol, controller->c_volume_single_step); CHKERRQ(ierr);
            ierr = VecPointwiseDivide(controller->c_surface_err_tmp, controller->c_surface_err_tol, controller->c_surface_single_step); CHKERRQ(ierr);

            ierr = VecScale(controller->c_volume_err_tmp, ratio);
            ierr = VecScale(controller->c_surface_err_tmp, ratio);

            // now it's actually new_h^2
            ierr = VecMin(controller->c_volume_err_tmp, PETSC_NULL, &volume_h); CHKERRQ(ierr);
            ierr = VecMin(controller->c_surface_err_tmp, PETSC_NULL, &surface_h); CHKERRQ(ierr);
 
            new_h = (volume_h<surface_h)?volume_h:surface_h;

            // now it's h, set tmp_dt
            controller->tmp_dt = sqrt(new_h);
/*
            if(controller->tmp_dt > dt_min){
                 controller->tmp_dt *= (2.0*dt_decrease_ratio);
            } else{
                 controller->tmp_dt = dt_min;
            }
*/
            ierr = VecCopy(controller->c_volume_tmp, user->c_volume_global); CHKERRQ(ierr);
            ierr = VecCopy(controller->c_surface_tmp, user->c_surface_global); CHKERRQ(ierr);

            controller->flg = 3;
        }
    } else {
        SETERRQ1(PETSC_COMM_SELF,1,"SingleStepSizeControl() in StepDoublingController: controller flag %d shouldn't appear here.\n", controller->flg);
    }

    PetscFunctionReturn(0);
}

PetscErrorCode DestroyTimeStepController(TIME_STEP_CONTROLLER *controller)
{
    PetscErrorCode ierr;

    PetscFunctionBegin;

    ierr = VecDestroy(&controller->c_volume_single_step);CHKERRQ(ierr);
    ierr = VecDestroy(&controller->c_volume_tmp);CHKERRQ(ierr);
    ierr = VecDestroy(&controller->c_volume_err_tmp);CHKERRQ(ierr);
    ierr = VecDestroy(&controller->c_volume_err_tol);CHKERRQ(ierr);

    ierr = VecDestroy(&controller->c_surface_single_step);CHKERRQ(ierr);
    ierr = VecDestroy(&controller->c_surface_tmp);CHKERRQ(ierr);
    ierr = VecDestroy(&controller->c_surface_err_tmp);CHKERRQ(ierr);
    ierr = VecDestroy(&controller->c_surface_err_tol);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}
#endif
