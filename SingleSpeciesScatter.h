#ifndef _SINGLE_SPECIES_SCATTER_H_
#define _SINGLE_SPECIES_SCATTER_H_

PetscErrorCode setSingleSpeciesScatter(Vec vec_all_species, Vec vec_single_species, PetscInt number_of_species, PetscInt number_of_grids, PetscInt i, VecScatter *vscat_ptr);
PetscErrorCode destroySingleSpeciesScatter(VecScatter *vscat_ptr);
PetscErrorCode getSingleSpeciesVec(Vec vec_all_species, Vec vec_single_species, VecScatter vscat);
PetscErrorCode restoreSingleSpeciesVec(Vec vec_all_species, Vec vec_single_species, VecScatter vscat);

/* 
 * set scatter for vector of a single species for either volume or surface
 * vec_all_species: vector of all species
 * vec_single_species: vector of a single species
 * number_of_species: number of species in volume or on surface
 * number_of_grids: number of total grids in volume or on surface
 * i: index of the species in species vector c[number_of_species] in corresponding field
 *
 */
PetscErrorCode setSingleSpeciesScatter(Vec vec_all_species, Vec vec_single_species, PetscInt number_of_species, PetscInt number_of_grids, PetscInt i, VecScatter *vscat_ptr)
{
     PetscErrorCode ierr;
     IS is;
     PetscInt vec_all_size, vec_single_size, all_local_size, single_local_size,rstart;

     PetscFunctionBegin;

     // check if sizes match
     ierr = VecGetSize(vec_all_species,&vec_all_size); CHKERRQ(ierr);
     if(vec_all_size != number_of_species*number_of_grids){
        SETERRQ3(PETSC_COMM_SELF,1,"setSingleSpeciesScatter(): %d(number_of_species)*%d(number_of_grids) != %d(size(vec_all_species)).\n",number_of_species, number_of_grids, vec_all_size);
     }
     ierr = VecGetSize(vec_single_species,&vec_single_size); CHKERRQ(ierr);
     if(vec_single_size != number_of_grids){
        SETERRQ(PETSC_COMM_SELF,1,"setSingleSpeciesScatter(): number_of_grids != size(vec_single_species).\n");
     }

     // get local size of each vector
     ierr = VecGetLocalSize(vec_all_species,&all_local_size); CHKERRQ(ierr);
     ierr = VecGetLocalSize(vec_single_species,&single_local_size); CHKERRQ(ierr);
     ierr = VecGetOwnershipRange(vec_all_species,&rstart,PETSC_NULL); CHKERRQ(ierr);

     // create is
     ierr = ISCreateStride(PETSC_COMM_WORLD,single_local_size,rstart+i,number_of_species,&is); CHKERRQ(ierr);

     // creat scatter
     ierr = VecScatterCreate(vec_all_species,is,vec_single_species,PETSC_NULL, vscat_ptr); CHKERRQ(ierr);

     // free is
     ierr = ISDestroy(&is); CHKERRQ(ierr);

     PetscFunctionReturn(0);
}

/*
 * Destroy single species scatter
 */
PetscErrorCode destroySingleSpeciesScatter(VecScatter *vscat_ptr)
{
     PetscErrorCode ierr;

     PetscFunctionBegin;

     ierr = VecScatterDestroy(vscat_ptr); CHKERRQ(ierr);

     PetscFunctionReturn(0);
}

/*
 * Using a single species scatter to get a single species vector from the whole vector
 */
PetscErrorCode getSingleSpeciesVec(Vec vec_all_species, Vec vec_single_species, VecScatter vscat)
{
     PetscErrorCode  ierr;

     PetscFunctionBegin;

     // get single species's value
     ierr = VecScatterBegin(vscat,vec_all_species,vec_single_species,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
     ierr = VecScatterEnd(vscat,vec_all_species,vec_single_species,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

     PetscFunctionReturn(0);
}

/*
 * Using a single species scatter to restore a single species vector to the whole vector
 */
PetscErrorCode restoreSingleSpeciesVec(Vec vec_all_species, Vec vec_single_species, VecScatter vscat)
{
     PetscErrorCode  ierr;

     PetscFunctionBegin;

     // restore single species's value to the whole vector
     ierr = VecScatterBegin(vscat,vec_single_species,vec_all_species,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
     ierr = VecScatterEnd(vscat,vec_single_species,vec_all_species,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

     PetscFunctionReturn(0);
}

#endif
