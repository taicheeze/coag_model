# Introduction
This is the repository for the 2D model for coagulation used by Petzold Group researchers at UCSB. The code was developed by Dr. Sheng Wu and Dr. Matthew Buoni during their time at UCSB. The current model for coagulation was created by Tie Bo Wu. This README was written by Tie Bo Wu and can be contacted via e-mail: tiebo(at)outlook(dot)com. 

## Download and install petsc
This model uses PETSc as its solver. This code was originally run with petsc version 3.4.0, and we cannot guarantee that it will work with the same efficiency with current or future versions of petsc. PETSc can be downloaded from the [official website](https://www.mcs.anl.gov/petsc/download/index.html)
Once downloaded and extracted, configure with the following commands:

1. $ ./configure PETSC_ARCH=arch-linux2-c-opt --download-superlu_dist --download-parmetis --download-metis --with-debugging=no

**DO NOT USE THIS ONE./configure PETSC_ARCH=arch-linux2-c-opt --with-cc=mpicc --with-fc=mpif90 --download-f-blas-lapack --download-mpich --download-superlu_dist --download-parmetis --download-metis --with-debugging=no**

2. $ make all test

3. Add $PETSC_DIR and $PETSC_ARCH into ~/.bash_profile or ~/.bashrc

## Running the Current Model

In order to run the model, run the following commands in the terminal, this should be the last steps **after** all modifications have been made:

1. $ make coagulationo3
2. $ make runcoagulation

## Analyzing data

After running the model, the data will be stored in the folder ./output unless specified elsewhere. The data is saved as a .m file to be parsed with MATLAB using the MATLAB function *custom_solver_data_analysis_general*. 

The syntax for this function is as follows:
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

Once they are ready to be used, copy these files and overwrite the files already in the main directory. You may adjust the settings or you can just run the model at this point.

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

For example, we can have a Michaelis-Menten reaction in which A is converted to B with a catalyst C. The rate at which A is converted to B is given by r1[C][A]/(r2+[A]), we can write that as: reaction_rate[48] = r_constant[48] * x_ptr[50] * x_ptr[39] / (r_constant[49] + x_ptr[39]); if 48 was the reaction number, 39 is the index of the catalyst, and 50 is the index for the substrate.

5. Re-define the jacobian matrix. The jacobian matrix in our solver has to be explicitly defined. Copy any additional constants that were used in the constants array above the jacobian matrix, these two arrays should be identical to avoid confusion. The jacobian is the first partial derivative of each species with respect to all species. For example, the entry J[0][0] refers to the partial derivative of the first species with respect to the first species. So for any concentrations that are affected by the reaction (involved in the stoichiometry matrix), we need to consider its partial derivative with respect to all other species involved in the reaction. For mass-action reactions, this part is done automatically, for non-mass action reactions, this part is also done automatically, but incorrectly, you will need to not only add the correct value, but you must also remove the incorrect value.

Using the same example: r1[C][A]/(r2+[A]), if the index values for A, B and C are 39, 54 and 50 respectively, then J[39][39], J[39][50], J[39][54], J[54][39], J[54][50], J[54][54] must be edited. Since the converion of A to B is 1 to 1, then J[39][39] = -J[54][39]. The partial derivative obtained from this reaction is ADDED to the end of the entry if an entry already exists, and a new entry is created if one does not exist. The automatically created, incorrect term may already be there, this must be deleted for correct calculations.

# Modifying the output

There are two main components to modify about the output, the outputed species and the output interval. Edit these files before compiling the solver with the $make coagulationo3 command.

## Modifying the output species
Since this is a PDE solver that outputs every concentration at every node, if we decide to output every species, a computer can struggle to parse such a large file. Instead we must specify which species we want to output.

1. In the main model directory, locate the file 'output_species_list.txt' and open it
2. The format of the file goes: number of volume species and below it the name of the volume species to track, followed by number of surface species, and the number of surface species it should track. Adjust this to reflect the desired output.

## Modifying the Output Interval
The max timestep in the solver is set to 1 second, however, outputting the concentrations at every second is often unnecessary and would produce large output files that may be difficult to parse, so a larger output interval is often used.

1. In the main model directory, locate the file 'Settings.h' and open it.
2. Locate the line that begins: "#define OUTPUT_DT" and edit the value on the same line to reflect the frequency at which to store output data.

# Modifying Initial Conditions

If you simply want to run the same model again, but with different initial conditions, you do not need to generate a new model. Simply modify the InitialConcentrations.h file and change the values in the concentration arrays. There will be a string array above the number array that specifies which element corresponds to which species. There is a value array for the initial concentrations of fluid species, as well as arrays for each patch of surface species. All values are in SI units, so mol/m^3.

# Modifying Other Model Parameters

Most of the other modifications to the model are done from the Settings.h file. We will only cover the lines we recommend editing in this readme file. For the modifications in this section, open the 'Settings.h' file.

## Modifying Grid Refinement

#define NUMBER_OF_X_GRIDS         50
#define NUMBER_OF_Y_GRIDS          25
#define NUMBER_OF_Y_GRIDS_BOTTOM   10 

The first two lines define the number of nodes along the x and y direction, the third line specifies how concentrated the nodes are to the bottom surface, the higher this value, the finer the mesh near the bottom surface.

## Modifying the Geometry

#define X_LENGTH          10000.0e-6
#define Y_LENGTH            10.0e-6
#define Y_BOTTOM_LENGTH      1.0e-6

The first two lines defines the size of the channel in meters, the third line defines the size of the bottom layer which will have disproportionately higher density of nodes.

## Modifying the Patches

#define NUMBER_OF_PATCHES 3

const PetscReal PATCH_LEFT_BOUNDARIES[NUMBER_OF_PATCHES]  = {0.0e-3, 3.0e-3, 5.0e-3};
const PetscReal PATCH_RIGHT_BOUNDARIES[NUMBER_OF_PATCHES] = {3.0e-3, 5.0e-3, 10.0e-3};

The first line defines the number of patches, in most of our simulations we use 3. The next two lines define the boundaries of the patches. The first array defines the left boundaries, the second array defines the right boundaries. In this example, we have 3 patches, the first patch from 0 to 3mm, then from 3mm to 5mm then 5mm to 10mm. The initial concentrations of these patches are defined in the InitialConcentrations.h header, with the first array in the initial surface concentrations representing the initial concentrations of the first patch etc... For more information see the **Creating a Injury Site "Patch"** portion of this Readme.

## Modifying the Boundary Conditions

#define LEFT_RIGHT_PERIODIC

If the above line is uncommented, the simulation will run with periodic boundary conditions. However, if the above line is deleted or commented, it will be run with an inlet and outlet and the concentrations at the inlet will be set to the initial concentrations of the fluid.

## Modifying the Time Stepping 

#define T_START          0.0
#define T_END         5000.0
#define DT_MIN           0.5e-4
#define DT_MAX           1.0
//#define DT_STEP          0.01 (leave commenetd)

#define OUTPUT_DT       20.0

This portion defines the start, end and max/min time steps for the solver as well as the output frequency. It not recommended that you do not change the dt_min or dt_max unless you are familiar with numerical solvers.

## Modifying the Tolerances

// absolute error tolerance 
#define ATOL 1.0e-9
// relative error tolerance
#define RTOL 5.0e-3

The absolute and relative tolerances can be adjusted here, we do not recommend adjusting these values unless the user is familiar with numerical solvers.

## Modifying Flow Velocity

#define VX_MAX           0

This creates a parabolic flow profile with max velocity defined above. The flow will go from the left boundary to the right boundary. The values are defined in meters per second.

## Modifying Diffusion Constant

#define VOLUME_DIFFUSION_CONSTANT  5.0e-11

This simulation uses a single value for diffusion constant defined here. The units are in meters^2/sec

# Instantaneous Concentration Changes (treatments)

If you wish to instantaneously add/remove a species to the system, you must adjust Coagulation.c. We use this to introduce drug interventions. In the variable definitions, notice the line: 

"PetscInt		       TXA=1,AP=1;" 

To turn on TXA or AP treatment, change their respective value to 0. 

To adjust modify the treatment modify this line: 

"PetscScalar		       scalefac=0.1,shiftfac=1.0e-3;"

scalefac represents the percent of plasminogen leftover after TXA is applied, shiftfac refers to the amount of antiplasmin added to the system in mols/m^3. 

To adjust the timing at which these treatments are applied, scroll toward the bottom and look for the lines that read:

if(t >= 4000 && TXA==0)

or

if(t >= 4000 && AP==0)

To adjust when the intervention is applied, change the value after "t>=", in this example, the intervention is applied after 4000 seconds.

# Additional Information
The file Coagulation.tmp can be opened to see what which timesteps have been solved, scroll to the bottom to see the last timestep that was solved.
