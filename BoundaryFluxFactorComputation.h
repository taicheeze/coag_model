#ifndef _BOUNDARY_FLUX_FACTORS_COMPUTATION_H_
#define _BOUNDARY_FLUX_FACTORS_COMPUTATION_H_

#include "DataStructure.h"
#include "SurfaceReactions.h"

typedef struct {
   Vec  c_s;             // surface species in one surface grid
   Vec  c_v;             // volume species in one surface grid

/*
 *     For each grid, for volume species k, the flux is
 *     F_k = stoich * surf_reactions
 *         = alpha_k + beta_k * c_v_k
 */
   Vec  alpha_single;    // boundary flux factor alpha for a single surface grid
   Vec  beta_single;    // boundary flux factor beta for a single surface grid

   Vec  alpha_global;    // boundary flux factor alpha for all surface grids
   Vec  beta_global;     // boundary flux factor beta for all surface grids

   PetscReal  ratio_layer0, ratio_layer1; // interpolation of boundary concentration
                                          // c_boundary[i] = ratio_layer0 * c[0][i] + ratio_layer1 * c[1][i]

   DM   da_alpha_beta;  // DMDA for alpha and beta
} BOUNDARY_FLUX_FACTOR_STRUCT;

PetscErrorCode InitializeBoundaryFluxFactorStruct(BOUNDARY_FLUX_FACTOR_STRUCT *bffs_ptr, MY_DATA *user);
PetscErrorCode DestroyBoundaryFluxFactorStruct(BOUNDARY_FLUX_FACTOR_STRUCT *bffs_ptr);

PetscErrorCode BoundaryFluxFactorsComputation(BOUNDARY_FLUX_FACTOR_STRUCT *bffs_ptr, PetscReal t, MY_DATA *user);
PetscErrorCode BoundaryFluxFactorsComputationSingleGrid(PetscReal t,Vec c_s,Vec c_v,Vec alpha_single, Vec beta_single);


/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 *     Form and set boundary flux factor struct
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
PetscErrorCode InitializeBoundaryFluxFactorStruct(BOUNDARY_FLUX_FACTOR_STRUCT *bffs_ptr, MY_DATA *user)
{
    PetscErrorCode ierr;
    PetscInt       number_of_volume_species, number_of_surface_species;

    PetscFunctionBegin;

    number_of_volume_species = NUMBER_OF_SPECIES_IN_VOLUME;
    number_of_surface_species = NUMBER_OF_SPECIES_ON_SURFACE;

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 *     Create c_s (vector for surface species on a single grid) without first assigning data array
 *     and create c_v (vector for volume species on a single surface grid) with assigning data array
 *     this is because c_v needs to be extrapolated from volume species concentration near the boundary
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    ierr = VecCreateSeqWithArray(MPI_COMM_SELF,1,number_of_surface_species,PETSC_NULL,&(bffs_ptr->c_s));CHKERRQ(ierr);

    ierr = VecCreate(MPI_COMM_SELF, &(bffs_ptr->c_v));CHKERRQ(ierr);
    ierr = VecSetType(bffs_ptr->c_v, VECSEQ);CHKERRQ(ierr);
    ierr = VecSetSizes(bffs_ptr->c_v, number_of_volume_species, number_of_volume_species);CHKERRQ(ierr);

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 *     Create alpha and beta for all the boundary grids
 *     For each grid, for volume species k, the flux is
 *     F_k = stoich * surf_reactions
 *         = alpha_k + beta_k * c_v_k
 *
 *     First create DMDA for it, then create global vector from DMDA
 *     create single grid alpha_single and beta_single without first
 *     assigning data array
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    ierr = DMDACreate1d(user->comm[1], DMDA_BOUNDARY_NONE, -NUMBER_OF_X_GRIDS, NUMBER_OF_SPECIES_IN_VOLUME,1,PETSC_NULL,&bffs_ptr->da_alpha_beta);CHKERRQ(ierr);

    ierr = DMCreateGlobalVector(bffs_ptr->da_alpha_beta,&(bffs_ptr->alpha_global));CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(bffs_ptr->da_alpha_beta,&(bffs_ptr->beta_global));CHKERRQ(ierr);

    ierr = VecCreateSeqWithArray(MPI_COMM_SELF,1,number_of_volume_species,PETSC_NULL,&(bffs_ptr->alpha_single));CHKERRQ(ierr);
    ierr = VecCreateSeqWithArray(MPI_COMM_SELF,1,number_of_volume_species,PETSC_NULL,&(bffs_ptr->beta_single));CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 *     Destroy boundary flux factor struct
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
PetscErrorCode DestroyBoundaryFluxFactorStruct(BOUNDARY_FLUX_FACTOR_STRUCT *bffs_ptr)
{
    PetscErrorCode ierr;

    PetscFunctionBegin;

    ierr = VecDestroy(&(bffs_ptr->c_s)); CHKERRQ(ierr);
    ierr = VecDestroy(&(bffs_ptr->c_v)); CHKERRQ(ierr);
    ierr = VecDestroy(&(bffs_ptr->alpha_single)); CHKERRQ(ierr);
    ierr = VecDestroy(&(bffs_ptr->beta_single)); CHKERRQ(ierr);
    ierr = VecDestroy(&(bffs_ptr->alpha_global)); CHKERRQ(ierr);
    ierr = VecDestroy(&(bffs_ptr->beta_global)); CHKERRQ(ierr);

    ierr = DMDestroy(&(bffs_ptr->da_alpha_beta)); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 *     Calculate boundary flux factos at all surface grids
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
PetscErrorCode BoundaryFluxFactorsComputation(BOUNDARY_FLUX_FACTOR_STRUCT *bffs_ptr, PetscReal t, MY_DATA *user)
{
    PetscErrorCode   ierr;
    PetscInt         xs,xm,xe;
    PetscInt         i,k,number_of_surface_species,number_of_volume_species;
    Species_Volume   **array_c_volume, *array_alpha_global, *array_beta_global;
    Species_Surface  *array_c_surface;
    PetscReal        **array_dy_volume, dy0, dy1, dy_ratio;
    PetscReal        *array_c_v_single;

    PetscFunctionBegin;

    number_of_volume_species  = NUMBER_OF_SPECIES_IN_VOLUME;
    number_of_surface_species = NUMBER_OF_SPECIES_ON_SURFACE;

/*
 *   Get a pointer to global vector data.
 */
    ierr = DMDAVecGetArray(user->da_c_volume, user->c_volume_global, &array_c_volume);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->da_c_surface, user->c_surface_global, &array_c_surface);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(bffs_ptr->da_alpha_beta, bffs_ptr->alpha_global, &array_alpha_global);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(bffs_ptr->da_alpha_beta, bffs_ptr->beta_global, &array_beta_global);CHKERRQ(ierr);

    ierr = DMDAVecGetArray(user->da_volume_dof_1, user->dy_volume_global, &array_dy_volume);CHKERRQ(ierr);

/*
 *  Get local grid boundaries (for 1-dimensional DMDA):
 *   xs   - starting grid indices (no ghost points)
 *   xm   - widths of local grid (no ghost points)
 */
    ierr = DMDAGetCorners(user->da_c_surface, &xs, PETSC_NULL, PETSC_NULL, &xm, PETSC_NULL, PETSC_NULL);CHKERRQ(ierr);

    xe = xs + xm;

    // put this out of the loop because it's the same for all the i's
    i = xs; 
        // calculate c_v with linear extrapolation
        dy0 = array_dy_volume[0][i]; // y=0, boundary layer
        dy1 = array_dy_volume[1][i]; // y=1, one layer up
        dy_ratio = dy0 / (dy0 + dy1);

        // set ratio_layer0 and ratio_layer1
        // c_boundary = ( 1 + dy0/(dy0+dy1) ) * c_0 - dy0/(dy0+dy1) * c_1
        bffs_ptr->ratio_layer0 = 1.0 + dy_ratio; 
        bffs_ptr->ratio_layer1 = 0.0 - dy_ratio;
       
/* 
 * Solve surface reaction time-stepping over the locally owned part of the grid
 * 
 */
    for (i=xs; i<xe; i++) {
        // place c_s to reflect surface species in single grid
        ierr = VecPlaceArray(bffs_ptr->c_s,array_c_surface[i].c); CHKERRQ(ierr);

        // calculate c_v in single surface grid
        ierr = VecGetArray(bffs_ptr->c_v,&array_c_v_single); CHKERRQ(ierr);
        // c_boundary = ( 1 + dy0/(dy0+dy1) ) * c_0 - dy0/(dy0+dy1) * c_1
        for (k=0; k<number_of_volume_species; ++k) {
            array_c_v_single[k] = ((PetscReal)1.0+dy_ratio) * array_c_volume[0][i].c[k] - dy_ratio * array_c_volume[1][i].c[k];
            if(array_c_v_single[k] < 0.0)
                array_c_v_single[k] = 0.0;
        }
        ierr = VecRestoreArray(bffs_ptr->c_v,&array_c_v_single); CHKERRQ(ierr);

        // place alpha and beta to reflect alpha and beta in single grid
        ierr = VecPlaceArray(bffs_ptr->alpha_single,array_alpha_global[i].c); CHKERRQ(ierr);
        ierr = VecPlaceArray(bffs_ptr->beta_single,array_beta_global[i].c); CHKERRQ(ierr);

        // calculate alpha and beta in the single grid
        ierr = BoundaryFluxFactorsComputationSingleGrid(t,bffs_ptr->c_s,bffs_ptr->c_v,bffs_ptr->alpha_single, bffs_ptr->beta_single); CHKERRQ(ierr);

        // reset c_s, alpha_single, and beta_single
        ierr = VecResetArray(bffs_ptr->c_s); CHKERRQ(ierr);
        ierr = VecResetArray(bffs_ptr->alpha_single); CHKERRQ(ierr);
        ierr = VecResetArray(bffs_ptr->beta_single); CHKERRQ(ierr);
    }

/*
 *  Restore array to vector
 */
    ierr = DMDAVecRestoreArray(user->da_c_volume, user->c_volume_global, &array_c_volume);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->da_c_surface, user->c_surface_global, &array_c_surface);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(bffs_ptr->da_alpha_beta, bffs_ptr->alpha_global, &array_alpha_global);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(bffs_ptr->da_alpha_beta, bffs_ptr->beta_global, &array_beta_global);CHKERRQ(ierr);

    ierr = DMDAVecRestoreArray(user->da_volume_dof_1, user->dy_volume_global, &array_dy_volume);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

/*
 *  Calculate boundary flux factors at a single surface grid
 */
PetscErrorCode BoundaryFluxFactorsComputationSingleGrid(PetscReal t,Vec c_s,Vec c_v,Vec alpha_single, Vec beta_single)
{
    PetscReal                reaction_rate[NUMBER_OF_REACTIONS_ON_SURFACE];
    PetscInt                 i, j, number_volume_species, number_surface_species, number_surface_reactions;
    PetscReal                *cs_ptr, *cv_ptr, *alpha_ptr, *beta_ptr, temp;
    PetscErrorCode           ierr;

    PetscFunctionBegin;

    number_volume_species   = NUMBER_OF_SPECIES_IN_VOLUME;
    number_surface_species   = NUMBER_OF_SPECIES_ON_SURFACE;
    number_surface_reactions = NUMBER_OF_REACTIONS_ON_SURFACE;

/*
 * Get the pointers pointing to data array of c_s, c_v, alpha_single, and beta_single respectively
 *  c_s: surface species concentrations
 *  c_v: volume species concentrations
 * For each grid, for volume species k, the flux is
 *   F_k = stoich * surf_reactions
 *       = alpha_k + beta_k * c_v_k
 */
    ierr = VecGetArray(c_s, &cs_ptr);CHKERRQ(ierr);
    ierr = VecGetArray(c_v, &cv_ptr);CHKERRQ(ierr);
    ierr = VecGetArray(alpha_single, &alpha_ptr);CHKERRQ(ierr);
    ierr = VecGetArray(beta_single, &beta_ptr);CHKERRQ(ierr);

    for(i=0;i<number_volume_species;++i){
        temp = cv_ptr[i];

        // set cv[i] to 0 to calculate alpha[i]
        cv_ptr[i] = 0.0;
        // Calling user-generated function from SurfaceReactions.h to calculate reation rate
        ierr = SurfaceReactionRateCalculation(t, cs_ptr, cv_ptr, reaction_rate); CHKERRQ(ierr);
        // calculate alpha[i]
        alpha_ptr[i] = 0.0;
        for(j=0;j<number_surface_reactions;++j){
            if(surface_stoichiometry_volume_species[i][j] != 0){
                 alpha_ptr[i] += (PetscReal)surface_stoichiometry_volume_species[i][j]*reaction_rate[j];
            }
        }

        // set cv[i] to 1 to calculate alpha[i]+beta[i]
        cv_ptr[i] = 1.0;
        // Calling user-generated function from SurfaceReactions.h to calculate reation rate
        ierr = SurfaceReactionRateCalculation(t, cs_ptr, cv_ptr, reaction_rate); CHKERRQ(ierr);
        // calculate alpha[i]+beta[i]
        beta_ptr[i] = 0.0;
        for(j=0;j<number_surface_reactions;++j){
            if(surface_stoichiometry_volume_species[i][j] != 0){
                 beta_ptr[i] += (PetscReal)surface_stoichiometry_volume_species[i][j]*reaction_rate[j];
            }
        }
        // calculate beta[i]
        beta_ptr[i] -= alpha_ptr[i];

        // restore cv[i]
        cv_ptr[i] = temp;
    }

/*
 * restore vectors c_s, c_v, alpha_single, and beta_single from data arrays
 */
    ierr = VecRestoreArray(c_s, &cs_ptr);CHKERRQ(ierr);
    ierr = VecRestoreArray(c_v, &cv_ptr);CHKERRQ(ierr);
    ierr = VecRestoreArray(alpha_single, &alpha_ptr);CHKERRQ(ierr);
    ierr = VecRestoreArray(beta_single, &beta_ptr);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

#endif
