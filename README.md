# Introduction
This is the repository for the 2D model for coagulation used by Petzold Group researchers at UCSB. The code was developed by Dr. Sheng Wu and Dr. Matthew Buoni during their time at UCSB.

## Download and install petsc
This model uses PETSc as its solver. This code was originally run with petsc version 3.4.0, and we cannot guarantee that it will work with the same efficiency with current or future versions of petsc. PETSc can be downloaded from the [official website](https://www.mcs.anl.gov/petsc/download/index.html)
Once downloaded and extracted, configure with the following commands:

1. $ ./configure PETSC_ARCH=arch-linux2-c-opt --download-superlu_dist --download-parmetis --download-metis --with-debugging=no

**DO NOT USE THIS ONE./configure PETSC_ARCH=arch-linux2-c-opt --with-cc=mpicc --with-fc=mpif90 --download-f-blas-lapack --download-mpich --download-superlu_dist --download-parmetis --download-metis --with-debugging=no**

2. $ make all test

3. Add $PETSC_DIR and $PETSC_ARCH into ~/.bash_profile or ~/.bashrc

## Running the Current Model

In order to run the model, run the following commands in the terminal:

1. $ make coagulationo3
2. $ make runcoagulation

## Analyzing data

After running the model, the data will be stored in the folder ./output unless specified elsewhere. The data is saved as a .m file to be parsed with MATLAB using the MATLAB function *custom_solver_data_analysis_general* The syntax for this function is as follows:
[t_list, conc_list, species_list, x, y, dx, dy] = custom_solver_data_analysis_general(number_of_intervals) where t_list is the time vector, conc_list is the species concentrations in a cell array, species list a cell array of output species names, x, y, dx and dy are coordinates and intervals of the mesh.
If you want average values over an area/line, rather than concentrations at every node, you can use the *coag_plot* MATLAB function to do this calculation. The syntax for the command is:
[t1, y] = coag_plot(C,t,species_list,dx,dy,cutoff,xrange,x,yrange,y)
where t1, and y are your output time and concentration arrays, C is the concentration cell array that is output from *custom_solver_data_analysis_general*, same with t, species_list, dx, x and y. xrange and yrange are two element vectors in the form [min_val max_val], cutoff is the index of the last volume species in the species list, since surface species have to be averaged over a line rather than an area.

# Creating your own model

We have a python script that helps us write up a template header for initial conditions, volume reactions and surface reactions. The python script is located in the folder ./orfeomodel_python/intrinsic_model. The command to run this script is: *python3 createPetscCode.py \'ModelFileName\'*

This script uses 3 files to write the template headers: model_file(any name), species.txt, initial_values.txt

To create a new model:
1. Edit species.txt with the list of all of your chemical species, list them in a numbered format and keep surface and bulk(fluid) species as two separate lists.

2. Edit initial_values.txt with the initial values for any species that starts off at a non-zero value. All other species will be assumed to start with initial concentration of 0.

3. Create or edit an existing model file, for an example, you can edit old files such as jan2017_model.txt. 
..* The model file first defines general constants, such as diffusivity, in this area, we also define the modifier which is the thickness of our boundary layer in order to calculate the boundary flux rates. The amplifier_on/off constants simply give the user an easy way to adjust all on and off binding reactions if they choose, these were left at 1 for our simulations. 
..* Next, the rate constants are defined, the default rates are defined in SI units, any other units must be defined so they can be adjusted to SI units for calculation. 
..* The next section defined the bulk reactions in the form: reaction #., reactants (separated by plus if more than 1), => (irreversible) or <=> (reversible), products(separated by plus if more than 1), space, kf# = *constant for forward reaction* space kr = *constant for reverse reaction*.
..* The last section surface reactions follows the same format except the key string to identify forward and reverse reactions are sf and sr for surface forward and reverse constants.
..* Once these files are set, run $*python3 createPetscCode.py \'ModelFileName\'* in the terminal to generate three header files InitialConcentrations.h, SurfaceReactions.h, VolumeReactions.h

These files are ready to be used if:

1. All stoichiometry coefficients are 1.
2. You want the entire bottom surface to be the injury site. (You may need to adjust the settings file to account for this)
3. All reactions are mass action.
..* Once they are ready to be used, copy these files and overwrite the files already in the main directory. You may adjust the settings or you can just run the model at this point.

### Adjusting the Stoichiometry Matrix
If you have a reaction in which there is a non-1 coefficient in the stoichiometry matrix, you must open the appropriate header file (if it's a volume reaction, open VolumeReactions.h, if it's a surface reaction, open SurfaceReactions.h), find where the stoichiometry is defined, it should be after the line "PetscInt  volume_stoichiometry[NUMBER_OF_SPECIES_IN_VOLUME][NUMBER_OF_REACTIONS_IN_VOLUME]"

In the stoichiometry matrix, the columns refer to the reaction number and the rows refer to the species number, adjust the constant to match your reaction.

### Creating a Injury Site "Patch"

If you wish for your injury site to have differing concentrations along the bottom boundary, you must define the concentrations of these sections separately. For example, the default situation has 3 patches, a non-reactive patch, the injury site followed by another non-reactive patch. Your InitialConcentrations.h file must be adjusted to reflect this. To do this, two arrays need to be added to the array after "const PetscReal SURFACE_INITIAL_CONCENTRATIONS_PATCHES[NUMBER_OF_PATCHES][NUMBER_OF_SPECIES_ON_SURFACE]"

Currently there is only 1 array after this line. If we call this array P2, to be used in the middle patch, then this portion needs to be modified as an array of 3 arrays in the format
={
  {P1
  },
  {P2
  },
  {P3
  }
 };
 
### Using a Non-Mass Action Reaction

When using the python script to add reaction that isn't in a standard mass action form, define the reaction considering only the reactants and products of the reaction so that the generated stoichiometry matrix will be correct (if you need to adjust this matrix, see the above section on modifying the stoichiometry matrix).

1. Run the python script to generate the header files
2. Open the header file with the non-standard reaction (if it's a volume reaction, open VolumeReactions.h, if it's a surface reaction, open SurfaceReactions.h)
3. Scroll down to where the reaction rates are defined, there will be a list of constants defined followed by the reaction rates. If you need help finding a specific reaction, scroll to the top where each reaction is listed, note that reversible reactions have been turned into two separate reactions.
4. Re-define the reaction rate using r_constants for the different constants, x_ptr for the different species (these are volume species in the volume header and surface species in the surface header), for the correct index number of species can be figured out with the list of reactions at the top as well and cv_ptr to refer to volume species in SurfaceReactions.h. Additional constants may be defined by adding more values to the r_constant array. 
..* For example, we can have a Michaelis-Menten reaction in which A is converted to B with a catalyst C. The rate at which A is converted to B is given by r1[C][A]/(r2+[A]), we can write that as: reaction_rate[48] = r_constant[48] * x_ptr[50] * x_ptr[39] / (r_constant[49] + x_ptr[39]); if 48 was the reaction number, 39 is the index of the catalyst, and 50 is the index for the substrate.
5. Re-define the jacobian matrix. The jacobian matrix in our solver has to be explicitly defined. Copy any additional constants that were used in the constants array above the jacobian matrix, these two arrays should be identical to avoid confusion. The jacobian is the first partial derivative of each species with respect to all species. For example, the entry J[0][0] refers to the partial derivative of the first species with respect to the first species. So for any concentrations that are affected by the reaction (involved in the stoichiometry matrix), we need to consider its partial derivative with respect to all other species involved in the reaction. For mass-action reactions, this part is done automatically, for non-mass action reactions, this part is also done automatically, but incorrectly, you will need to not only add the correct value, but you must also remove the incorrect value.
..*Using the same example: r1[C][A]/(r2+[A]), if the index values for A, B and C are 39, 54 and 50 respectively, then J[39][39], J[39][50], J[39][54], J[54][39], J[54][50], J[54][54] must be edited. Since the converion of A to B is 1 to 1, then J[39][39] = -J[54][39]. The partial derivative obtained from this reaction is ADDED to the end of the entry if an entry already exists, and a new entry is created if one does not exist. The automatically created, incorrect term may already be there, this must be deleted for correct calculations.

