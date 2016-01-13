#ifndef _HEADERS_H_
#define _HEADERS_H_

//#ifdef PETSC_USE_LOG
//#undef PETSC_USE_LOG
//#endif

#include <stdbool.h>
#include <math.h>
#include <petscdmda.h>
#include <petscts.h>
#include <petscsys.h>

#include "DataStructure.h"
#include "InitialConcentrations.h"
#include "SingleSpeciesScatter.h"
#include "InitializeMY_DATA.h"
#include "FormGrid.h"
#include "FormVelocityField.h"
#include "FormInitialSolution.h"
#include "VolumeReactionSolver.h"
#include "VolumeConvectionSolver.h"
#include "SurfaceReactionSolver.h"
#include "VolumeDiffusionSolver.h"
#include "OutputHandler.h"
//#include "BoundaryFluxFactorComputation.h"

#ifdef NON_STEADY_FLOW
#include "UpdateVelocityField.h"
#endif

#if defined USE_MAX_CHANGE_STEP_CONTROLLER
#include "MaxChangeTimeStepController.h"
#elif defined USE_STEP_DOUBLING_STEP_CONTROLLER
#include "StepDoublingStepController.h"
#endif

#include "DestroyMY_DATA.h"

#endif
