
#ifndef _VOLUME_REACTIONS_H_
#define _VOLUME_REACTIONS_H_

#include "DataStructure.h"

/* VOLUME REACTIONS:
 *
 *  0.   Xa + VII -> Xa + VIIa                                                           ([2] + [0] -> [2] + [1])
 *  1.   IIa + VII -> IIa + VIIa                                                         ([3] + [0] -> [3] + [1])
 *  2.   Xa + II -> Xa + IIa                                                             ([2] + [7] -> [2] + [3])
 *  3.   IIa + VIII -> IIa + VIIIa                                                       ([3] + [8] -> [3] + [9])
 *  4.   VIIIa -> VIIIa1L + VIIIa2                                                       ([9] -> [10] + [11])
 *  5.   VIIIa1L + VIIIa2 -> VIIIa                                                       ([10] + [11] -> [9])
 *  6.   IIa + V -> IIa + Va                                                             ([3] + [12] -> [3] + [13])
 *  7.   Xa + TFPI -> XaTFPI                                                             ([2] + [15] -> [16])
 *  8.   XaTFPI -> Xa + TFPI                                                             ([16] -> [2] + [15])
 *  9.   Xa + ATIII -> XaATIII                                                           ([2] + [17] -> [18])
 *  10.  mIIa + ATIII -> mIIaATIII                                                       ([14] + [17] -> [19])
 *  11.  IXa + ATIII -> IXaATIII                                                         ([6] + [17] -> [20])
 *  12.  IIa + ATIII -> IIaATIII                                                         ([3] + [17] -> [21])
 *  13.  IIa + substrate -> IIa + fluorescence                                           ([3] + [22] -> [3] + [23])
 *  14.  IXa + X -> IXa + Xa                                                             ([6] + [4] -> [6] + [2])
 *  15.  mIIa + V -> mIIa + Va                                                           ([14] + [12] -> [14] + [13])
 *  16.  APC + Va -> APCVa                                                               ([24] + [13] -> [25])
 *  17.  APCVa -> APC + Va                                                               ([25] -> [24] + [13])
 *  18.  APCVa -> APC + Va5                                                              ([25] -> [24] + [26])
 *  19.  APCVa -> APC + Va3                                                              ([25] -> [24] + [27])
 *  20.  APC + Va5 -> APCVa5                                                             ([24] + [26] -> [28])
 *  21.  APCVa5 -> APC + Va5                                                             ([28] -> [24] + [26])
 *  22.  APC + Va3 -> APCVa3                                                             ([24] + [27] -> [29])
 *  23.  APCVa3 -> APC + Va3                                                             ([29] -> [24] + [27])
 *  24.  APCVa5 -> APC + Va53                                                            ([28] -> [24] + [30])
 *  25.  APCVa3 -> APC + Va53                                                            ([29] -> [24] + [30])
 *  26.  Va3 -> HCF + LCA1                                                               ([27] -> [31] + [32])
 *  27.  Va53 -> HCF + LCA1                                                              ([30] -> [31] + [32])
 *  28.  APC + LCA1 -> APCLCA1                                                           ([24] + [32] -> [33])
 *  29.  APCLCA1 -> APC + LCA1                                                           ([33] -> [24] + [32])
 */


/*
 * volume_stoichiometry: stoichiometry matrix of volume species for volume reactions
 */
/*
 * Begin of user-generated volume reaction stoichiometries
 */
PetscInt  volume_stoichiometry[NUMBER_OF_SPECIES_IN_VOLUME][NUMBER_OF_REACTIONS_IN_VOLUME]
          ={
            {-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0,-1, 1,-1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 1,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 1, 1,-1, 1,-1, 1, 1, 1, 0, 0,-1, 1},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,-1, 1, 0, 0,-1, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0,-1, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0,-1, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0,-1, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,-1, 1},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
           };


/*
 * End of user-generated volume reactions stoichiometries
 */

/*
 * User-generated routine to calculate reaction rates of surface reactions in a single grid
 * x_ptr: pointer to volume species array 
 * reaction_rate: pointer to reaction rates array
 */
PetscErrorCode VolumeReactionRateCalculation(PetscReal t,PetscReal x_ptr[NUMBER_OF_SPECIES_IN_VOLUME],PetscReal reaction_rate[NUMBER_OF_REACTIONS_IN_VOLUME])
{
    PetscReal r_constant[NUMBER_OF_REACTIONS_IN_VOLUME];
    PetscFunctionBegin;

/*
 * Begin of user-generated volume reaction rates
 */
    r_constant[0] = 1.3000e+04;
    r_constant[1] = 2.3000e+01;
    r_constant[2] = 7.5000e+00;
    r_constant[3] = 2.0000e+04;
    r_constant[4] = 6.0000e-03;
    r_constant[5] = 2.2000e+01;
    r_constant[6] = 2.0000e+04;
    r_constant[7] = 9.0000e+02;
    r_constant[8] = 3.6000e-04;
    r_constant[9] = 1.5000e+00;
    r_constant[10] = 7.1000e+00;
    r_constant[11] = 4.9000e-01;
    r_constant[12] = 7.1000e+00;
    r_constant[13] = 1.0000e+00;
    r_constant[14] = 5.7000e+00;
    r_constant[15] = 3.0000e+03;
    r_constant[16] = 1.0000e+05;
    r_constant[17] = 7.0000e-01;
    r_constant[18] = 1.0000e+00;
    r_constant[19] = 1.9200e-01;
    r_constant[20] = 1.0000e+05;
    r_constant[21] = 7.0000e-01;
    r_constant[22] = 1.0000e+05;
    r_constant[23] = 7.0000e-01;
    r_constant[24] = 1.0000e+00;
    r_constant[25] = 1.0000e+00;
    r_constant[26] = 2.8000e-02;
    r_constant[27] = 2.8000e-02;
    r_constant[28] = 1.0000e+05;
    r_constant[29] = 7.0000e-01;

    reaction_rate[0] = r_constant[0] * x_ptr[0] * x_ptr[2];
    reaction_rate[1] = r_constant[1] * x_ptr[0] * x_ptr[3];
    reaction_rate[2] = r_constant[2] * x_ptr[2] * x_ptr[7];
    reaction_rate[3] = r_constant[3] * x_ptr[3] * x_ptr[8];
    reaction_rate[4] = r_constant[4] * x_ptr[9];
    reaction_rate[5] = r_constant[5] * x_ptr[10] * x_ptr[11];
    reaction_rate[6] = r_constant[6] * x_ptr[3] * x_ptr[12];
    reaction_rate[7] = r_constant[7] * x_ptr[2] * x_ptr[15];
    reaction_rate[8] = r_constant[8] * x_ptr[16];
    reaction_rate[9] = r_constant[9] * x_ptr[2] * x_ptr[17];
    reaction_rate[10] = r_constant[10] * x_ptr[14] * x_ptr[17];
    reaction_rate[11] = r_constant[11] * x_ptr[6] * x_ptr[17];
    reaction_rate[12] = r_constant[12] * x_ptr[3] * x_ptr[17];
    reaction_rate[13] = r_constant[13] * x_ptr[3] * x_ptr[22];
    reaction_rate[14] = r_constant[14] * x_ptr[4] * x_ptr[6];
    reaction_rate[15] = r_constant[15] * x_ptr[12] * x_ptr[14];
    reaction_rate[16] = r_constant[16] * x_ptr[13] * x_ptr[24];
    reaction_rate[17] = r_constant[17] * x_ptr[25];
    reaction_rate[18] = r_constant[18] * x_ptr[25];
    reaction_rate[19] = r_constant[19] * x_ptr[25];
    reaction_rate[20] = r_constant[20] * x_ptr[24] * x_ptr[26];
    reaction_rate[21] = r_constant[21] * x_ptr[28];
    reaction_rate[22] = r_constant[22] * x_ptr[24] * x_ptr[27];
    reaction_rate[23] = r_constant[23] * x_ptr[29];
    reaction_rate[24] = r_constant[24] * x_ptr[28];
    reaction_rate[25] = r_constant[25] * x_ptr[29];
    reaction_rate[26] = r_constant[26] * x_ptr[27];
    reaction_rate[27] = r_constant[27] * x_ptr[30];
    reaction_rate[28] = r_constant[28] * x_ptr[24] * x_ptr[32];
    reaction_rate[29] = r_constant[29] * x_ptr[33];

/*
 * End of user-generated volume reaction rates
 */

    PetscFunctionReturn(0);
}

/*
 * User-generated routine to calculate jacobian matrix of volume reactions for volume species in a single grid
 * x_ptr: pointer to volume species array 
 * J: pointer to jacobian matrix for volume species
 */
PetscErrorCode VolumeReactionJacobianCalculation(PetscReal t,PetscReal x_ptr[NUMBER_OF_SPECIES_IN_VOLUME], PetscReal J[NUMBER_OF_SPECIES_IN_VOLUME][NUMBER_OF_SPECIES_IN_VOLUME])
{
    PetscReal r_constant[NUMBER_OF_REACTIONS_IN_VOLUME];

    PetscFunctionBegin;

/*
 * Begin of user-generated jacobian function
 */
    r_constant[0] = 1.3000e+04;
    r_constant[1] = 2.3000e+01;
    r_constant[2] = 7.5000e+00;
    r_constant[3] = 2.0000e+04;
    r_constant[4] = 6.0000e-03;
    r_constant[5] = 2.2000e+01;
    r_constant[6] = 2.0000e+04;
    r_constant[7] = 9.0000e+02;
    r_constant[8] = 3.6000e-04;
    r_constant[9] = 1.5000e+00;
    r_constant[10] = 7.1000e+00;
    r_constant[11] = 4.9000e-01;
    r_constant[12] = 7.1000e+00;
    r_constant[13] = 1.0000e+00;
    r_constant[14] = 5.7000e+00;
    r_constant[15] = 3.0000e+03;
    r_constant[16] = 1.0000e+05;
    r_constant[17] = 7.0000e-01;
    r_constant[18] = 1.0000e+00;
    r_constant[19] = 1.9200e-01;
    r_constant[20] = 1.0000e+05;
    r_constant[21] = 7.0000e-01;
    r_constant[22] = 1.0000e+05;
    r_constant[23] = 7.0000e-01;
    r_constant[24] = 1.0000e+00;
    r_constant[25] = 1.0000e+00;
    r_constant[26] = 2.8000e-02;
    r_constant[27] = 2.8000e-02;
    r_constant[28] = 1.0000e+05;
    r_constant[29] = 7.0000e-01;

    J[0][0] =  -1*r_constant[0]*x_ptr[2] -1*r_constant[1]*x_ptr[3];
    J[0][2] =  -1*r_constant[0]*x_ptr[0];
    J[0][3] =  -1*r_constant[1]*x_ptr[0];
    J[1][0] =  +1*r_constant[0]*x_ptr[2] +1*r_constant[1]*x_ptr[3];
    J[1][2] =  +1*r_constant[0]*x_ptr[0];
    J[1][3] =  +1*r_constant[1]*x_ptr[0];
    J[2][2] =  -1*r_constant[7]*x_ptr[15] -1*r_constant[9]*x_ptr[17];
    J[2][4] =  +1*r_constant[14]*x_ptr[6];
    J[2][6] =  +1*r_constant[14]*x_ptr[4];
    J[2][15] =  -1*r_constant[7]*x_ptr[2];
    J[2][16] =  +1*r_constant[8];
    J[2][17] =  -1*r_constant[9]*x_ptr[2];
    J[3][2] =  +1*r_constant[2]*x_ptr[7];
    J[3][3] =  -1*r_constant[12]*x_ptr[17];
    J[3][7] =  +1*r_constant[2]*x_ptr[2];
    J[3][17] =  -1*r_constant[12]*x_ptr[3];
    J[4][4] =  -1*r_constant[14]*x_ptr[6];
    J[4][6] =  -1*r_constant[14]*x_ptr[4];
    J[6][6] =  -1*r_constant[11]*x_ptr[17];
    J[6][17] =  -1*r_constant[11]*x_ptr[6];
    J[7][2] =  -1*r_constant[2]*x_ptr[7];
    J[7][7] =  -1*r_constant[2]*x_ptr[2];
    J[8][3] =  -1*r_constant[3]*x_ptr[8];
    J[8][8] =  -1*r_constant[3]*x_ptr[3];
    J[9][3] =  +1*r_constant[3]*x_ptr[8];
    J[9][8] =  +1*r_constant[3]*x_ptr[3];
    J[9][9] =  -1*r_constant[4];
    J[9][10] =  +1*r_constant[5]*x_ptr[11];
    J[9][11] =  +1*r_constant[5]*x_ptr[10];
    J[10][9] =  +1*r_constant[4];
    J[10][10] =  -1*r_constant[5]*x_ptr[11];
    J[10][11] =  -1*r_constant[5]*x_ptr[10];
    J[11][9] =  +1*r_constant[4];
    J[11][10] =  -1*r_constant[5]*x_ptr[11];
    J[11][11] =  -1*r_constant[5]*x_ptr[10];
    J[12][3] =  -1*r_constant[6]*x_ptr[12];
    J[12][12] =  -1*r_constant[6]*x_ptr[3] -1*r_constant[15]*x_ptr[14];
    J[12][14] =  -1*r_constant[15]*x_ptr[12];
    J[13][3] =  +1*r_constant[6]*x_ptr[12];
    J[13][12] =  +1*r_constant[6]*x_ptr[3] +1*r_constant[15]*x_ptr[14];
    J[13][13] =  -1*r_constant[16]*x_ptr[24];
    J[13][14] =  +1*r_constant[15]*x_ptr[12];
    J[13][24] =  -1*r_constant[16]*x_ptr[13];
    J[13][25] =  +1*r_constant[17];
    J[14][14] =  -1*r_constant[10]*x_ptr[17];
    J[14][17] =  -1*r_constant[10]*x_ptr[14];
    J[15][2] =  -1*r_constant[7]*x_ptr[15];
    J[15][15] =  -1*r_constant[7]*x_ptr[2];
    J[15][16] =  +1*r_constant[8];
    J[16][2] =  +1*r_constant[7]*x_ptr[15];
    J[16][15] =  +1*r_constant[7]*x_ptr[2];
    J[16][16] =  -1*r_constant[8];
    J[17][2] =  -1*r_constant[9]*x_ptr[17];
    J[17][3] =  -1*r_constant[12]*x_ptr[17];
    J[17][6] =  -1*r_constant[11]*x_ptr[17];
    J[17][14] =  -1*r_constant[10]*x_ptr[17];
    J[17][17] =  -1*r_constant[9]*x_ptr[2] -1*r_constant[10]*x_ptr[14] -1*r_constant[11]*x_ptr[6] -1*r_constant[12]*x_ptr[3];
    J[18][2] =  +1*r_constant[9]*x_ptr[17];
    J[18][17] =  +1*r_constant[9]*x_ptr[2];
    J[19][14] =  +1*r_constant[10]*x_ptr[17];
    J[19][17] =  +1*r_constant[10]*x_ptr[14];
    J[20][6] =  +1*r_constant[11]*x_ptr[17];
    J[20][17] =  +1*r_constant[11]*x_ptr[6];
    J[21][3] =  +1*r_constant[12]*x_ptr[17];
    J[21][17] =  +1*r_constant[12]*x_ptr[3];
    J[22][3] =  -1*r_constant[13]*x_ptr[22];
    J[22][22] =  -1*r_constant[13]*x_ptr[3];
    J[23][3] =  +1*r_constant[13]*x_ptr[22];
    J[23][22] =  +1*r_constant[13]*x_ptr[3];
    J[24][13] =  -1*r_constant[16]*x_ptr[24];
    J[24][24] =  -1*r_constant[16]*x_ptr[13] -1*r_constant[20]*x_ptr[26] -1*r_constant[22]*x_ptr[27] -1*r_constant[28]*x_ptr[32];
    J[24][25] =  +1*r_constant[17] +1*r_constant[18] +1*r_constant[19];
    J[24][26] =  -1*r_constant[20]*x_ptr[24];
    J[24][27] =  -1*r_constant[22]*x_ptr[24];
    J[24][28] =  +1*r_constant[21] +1*r_constant[24];
    J[24][29] =  +1*r_constant[23] +1*r_constant[25];
    J[24][32] =  -1*r_constant[28]*x_ptr[24];
    J[24][33] =  +1*r_constant[29];
    J[25][13] =  +1*r_constant[16]*x_ptr[24];
    J[25][24] =  +1*r_constant[16]*x_ptr[13];
    J[25][25] =  -1*r_constant[17] -1*r_constant[18] -1*r_constant[19];
    J[26][24] =  -1*r_constant[20]*x_ptr[26];
    J[26][25] =  +1*r_constant[18];
    J[26][26] =  -1*r_constant[20]*x_ptr[24];
    J[26][28] =  +1*r_constant[21];
    J[27][24] =  -1*r_constant[22]*x_ptr[27];
    J[27][25] =  +1*r_constant[19];
    J[27][27] =  -1*r_constant[22]*x_ptr[24] -1*r_constant[26];
    J[27][29] =  +1*r_constant[23];
    J[28][24] =  +1*r_constant[20]*x_ptr[26];
    J[28][26] =  +1*r_constant[20]*x_ptr[24];
    J[28][28] =  -1*r_constant[21] -1*r_constant[24];
    J[29][24] =  +1*r_constant[22]*x_ptr[27];
    J[29][27] =  +1*r_constant[22]*x_ptr[24];
    J[29][29] =  -1*r_constant[23] -1*r_constant[25];
    J[30][28] =  +1*r_constant[24];
    J[30][29] =  +1*r_constant[25];
    J[30][30] =  -1*r_constant[27];
    J[31][27] =  +1*r_constant[26];
    J[31][30] =  +1*r_constant[27];
    J[32][24] =  -1*r_constant[28]*x_ptr[32];
    J[32][27] =  +1*r_constant[26];
    J[32][30] =  +1*r_constant[27];
    J[32][32] =  -1*r_constant[28]*x_ptr[24];
    J[32][33] =  +1*r_constant[29];
    J[33][24] =  +1*r_constant[28]*x_ptr[32];
    J[33][32] =  +1*r_constant[28]*x_ptr[24];
    J[33][33] =  -1*r_constant[29];

/*
 * End of user-generated jacobian function
 */

    PetscFunctionReturn(0);
}

#endif
