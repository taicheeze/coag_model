
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
 *  30.  APC -> NULL                                                                     ([24] -> )
 *  31.  Fg + IIa -> Fg_IIa                                                              ([35] + [3] -> [37])
 *  32.  Fg_IIa -> Fg + IIa                                                              ([37] -> [35] + [3])
 *  33.  Fg_IIa -> FnI + IIa + FPA                                                       ([37] -> [36] + [3] + [40])
 *  34.  FnI + IIa -> FnI_IIa                                                            ([36] + [3] -> [38])
 *  35.  FnI_IIa -> FnI + IIa                                                            ([38] -> [36] + [3])
 *  36.  FnI_IIa -> FnII + IIa + FPB                                                     ([38] -> [39] + [3] + [41])
 *  37.  FnI -> FnI2                                                                     ([36] -> [42])
 *  38.  FnI2 -> FnI                                                                     ([42] -> [36])
 *  39.  FnI2 + IIa -> FnI2_IIa                                                          ([42] + [3] -> [43])
 *  40.  FnI2_IIa -> FnI2 + IIa                                                          ([43] -> [42] + [3])
 *  41.  FnI2_IIa -> FnII2 + IIa + FPB                                                   ([43] -> [44] + [3] + [41])
 *  42.  FnII + IIa -> FnII_IIa                                                          ([39] + [3] -> [45])
 *  43.  FnII_IIa -> FnII + IIa                                                          ([45] -> [39] + [3])
 *  44.  FnI2_IIa + ATIII -> FnI2_IIa_ATIII                                              ([43] + [17] -> [46])
 *  45.  FnI_IIa + ATIII -> FnI_IIa_ATIII                                                ([38] + [17] -> [47])
 *  46.  FnII_IIa + ATIII -> FnII_IIa_ATIII                                              ([45] + [17] -> [48])
 *  47.  Pn + AP -> Pn_AP                                                                ([50] + [51] -> [52])
 *  48.  FnII -> FDP                                                                     ([39] -> [54])
 *  49.  FnII2 -> FDP                                                                    ([44] -> [54])
 *  50.  stPA + PAI -> stPA_PAI                                                          ([55] + [53] -> [56])
 *  51.  stPA_PAI -> stPA + PAI                                                          ([56] -> [55] + [53])
 *  52.  Pg -> Pn                                                                        ([49] -> [50])
 *  53.  stPA -> NULL                                                                    ([55] -> )
 *  54.  FnI -> FDP                                                                      ([36] -> [54])
 *  55.  FnI2 -> FDP                                                                     ([42] -> [54])
 *  56.  PAI + APC -> NULL                                                               ([53] + [24] -> )
 */


/*
 * volume_stoichiometry: stoichiometry matrix of volume species for volume reactions
 */
/*
 * Begin of user-generated volume reaction stoichiometries
 */
PetscInt  volume_stoichiometry[NUMBER_OF_SPECIES_IN_VOLUME][NUMBER_OF_REACTIONS_IN_VOLUME]
          ={
            {-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0,-1, 1,-1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 1,-1, 1, 1, 0, 0,-1, 1, 1,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 1,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 1, 1,-1, 1,-1, 1, 1, 1, 0, 0,-1, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,-1, 1, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 1, 0,-2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1,-1, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 1, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0,-1},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 1, 2, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 0,-1, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0}
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
    r_constant[30] = 1.1111e-03;
    r_constant[31] = 1.0000e+05;
    r_constant[32] = 7.2000e+02;
    r_constant[33] = 8.4000e+01;
    r_constant[34] = 1.0000e+05;
    r_constant[35] = 7.5000e+02;
    r_constant[36] = 7.4000e+00;
    r_constant[37] = 1.0000e+03;
    r_constant[38] = 6.4000e-02;
    r_constant[39] = 1.0000e+05;
    r_constant[40] = 7.5000e+02;
    r_constant[41] = 4.9000e+01;
    r_constant[42] = 1.0000e+05;
    r_constant[43] = 1.0000e+03;
    r_constant[44] = 1.6000e+01;
    r_constant[45] = 1.6000e+01;
    r_constant[46] = 1.0000e+01;
    r_constant[47] = 4.0000e+02; //AP to Pn binding changed from 3e3
    r_constant[48] = 4.7000e-01; //fibrin degradation to FDP formation, increased from 4.7e-1
    r_constant[49] = 2.1000e-03; //fibrin degradation to FDP formation
    r_constant[50] = 4.0000e+04; //pai stpa binding
    r_constant[51] = 1.6000e-02; //increased from 1.6e-2, stpa_pai unbinding
    r_constant[52] = 1.0000e+00; //pai stpa binding
    r_constant[53] = 1.0000e+00; //stpa_pai unbinding
    r_constant[54] = 7.7000e-05; //Pg to pn conversion
    r_constant[55] = 4.1000e-04; //Pg to pn conversion
    r_constant[56] = 3.0000e-04; //Pg to pn conversion
    r_constant[57] = 9.0000e-02; //should be 9e-2
    r_constant[58] = 1.0000e-03; //stpa_pai unbinding
    r_constant[59] = 2.8900e-03; //stpa hepatic clearance
    r_constant[60] = 3.1300e-01; //protofibril degradation to FDP formation
    r_constant[61] = 1.8000e+02; //APC to PAI binding default 1.8e2

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
    reaction_rate[30] = r_constant[30] * x_ptr[24];
    reaction_rate[31] = r_constant[31] * x_ptr[3] * x_ptr[35];
    reaction_rate[32] = r_constant[32] * x_ptr[37];
    reaction_rate[33] = r_constant[33] * x_ptr[37];
    reaction_rate[34] = r_constant[34] * x_ptr[3] * x_ptr[36];
    reaction_rate[35] = r_constant[35] * x_ptr[38];
    reaction_rate[36] = r_constant[36] * x_ptr[38];
    reaction_rate[37] = r_constant[37] * pow(x_ptr[36],2.0);
    reaction_rate[38] = r_constant[38] * x_ptr[42];
    reaction_rate[39] = r_constant[39] * x_ptr[3] * x_ptr[42];
    reaction_rate[40] = r_constant[40] * x_ptr[43];
    reaction_rate[41] = r_constant[41] * x_ptr[43];
    reaction_rate[42] = r_constant[42] * x_ptr[3] * x_ptr[39];
    reaction_rate[43] = r_constant[43] * x_ptr[45];
    reaction_rate[44] = r_constant[44] * x_ptr[17] * x_ptr[43];
    reaction_rate[45] = r_constant[45] * x_ptr[17] * x_ptr[38];
    reaction_rate[46] = r_constant[46] * x_ptr[17] * x_ptr[45];
    reaction_rate[47] = r_constant[47] * x_ptr[50] * x_ptr[51];
    reaction_rate[48] = r_constant[48] * x_ptr[50] * x_ptr[39] / (r_constant[49] + x_ptr[39]);
    reaction_rate[49] = r_constant[48] * x_ptr[50] * x_ptr[44] / (r_constant[49] + x_ptr[44]);
    reaction_rate[50] = r_constant[52] * r_constant[50] * x_ptr[53] * x_ptr[55];
    reaction_rate[51] = r_constant[53] * r_constant[51] * x_ptr[56] * (x_ptr[39]+2*x_ptr[44])/(r_constant[58] + (x_ptr[39]+2*x_ptr[44]));
    reaction_rate[52] = r_constant[57] * x_ptr[55] * x_ptr[49]*(x_ptr[39]+2*x_ptr[44])/(r_constant[55]*((x_ptr[39]+2*x_ptr[44])+r_constant[56])+((x_ptr[39]+2*x_ptr[44])*x_ptr[50]+r_constant[54])*x_ptr[49]);
    reaction_rate[53] = r_constant[59] * x_ptr[55];
    reaction_rate[54] = r_constant[60] * x_ptr[50] * x_ptr[36] / (r_constant[49] + x_ptr[36]);
    reaction_rate[55] = r_constant[60] * x_ptr[50] * x_ptr[42] / (r_constant[49] + x_ptr[42]);
    reaction_rate[56] = r_constant[61] * x_ptr[24] * x_ptr[53];

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
    r_constant[30] = 1.1111e-03;
    r_constant[31] = 1.0000e+05;
    r_constant[32] = 7.2000e+02;
    r_constant[33] = 8.4000e+01;
    r_constant[34] = 1.0000e+05;
    r_constant[35] = 7.5000e+02;
    r_constant[36] = 7.4000e+00;
    r_constant[37] = 1.0000e+03;
    r_constant[38] = 6.4000e-02;
    r_constant[39] = 1.0000e+05;
    r_constant[40] = 7.5000e+02;
    r_constant[41] = 4.9000e+01;
    r_constant[42] = 1.0000e+05;
    r_constant[43] = 1.0000e+03;
    r_constant[44] = 1.6000e+01;
    r_constant[45] = 1.6000e+01;
    r_constant[46] = 1.0000e+01;
    r_constant[47] = 4.0000e+02; //AP to Pn binding changed from 3e3
    r_constant[48] = 4.7000e-01; //fibrin degradation to FDP formation, increased from 4.7e-1
    r_constant[49] = 2.1000e-03; //fibrin degradation to FDP formation
    r_constant[50] = 4.0000e+04; //pai stpa binding
    r_constant[51] = 1.6000e-02; //increased from 1.6e-2, stpa_pai unbinding
    r_constant[52] = 1.0000e+00; //pai stpa binding
    r_constant[53] = 1.0000e+00; //stpa_pai unbinding
    r_constant[54] = 7.7000e-05; //Pg to pn conversion
    r_constant[55] = 4.1000e-04; //Pg to pn conversion
    r_constant[56] = 3.0000e-04; //Pg to pn conversion
    r_constant[57] = 9.0000e-02; //should be 9e-2
    r_constant[58] = 1.0000e-03; //stpa_pai unbinding
    r_constant[59] = 2.8900e-03; //stpa hepatic clearance
    r_constant[60] = 3.1300e-01; //protofibril degradation to FDP formation
    r_constant[61] = 1.8000e+02; //APC to PAI binding default 1.8e2

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
    J[3][3] =  -1*r_constant[12]*x_ptr[17] -1*r_constant[31]*x_ptr[35] -1*r_constant[34]*x_ptr[36] -1*r_constant[39]*x_ptr[42] -1*r_constant[42]*x_ptr[39];
    J[3][7] =  +1*r_constant[2]*x_ptr[2];
    J[3][17] =  -1*r_constant[12]*x_ptr[3];
    J[3][35] =  -1*r_constant[31]*x_ptr[3];
    J[3][36] =  -1*r_constant[34]*x_ptr[3];
    J[3][37] =  +1*r_constant[32] +1*r_constant[33];
    J[3][38] =  +1*r_constant[35] +1*r_constant[36];
    J[3][39] =  -1*r_constant[42]*x_ptr[3];
    J[3][42] =  -1*r_constant[39]*x_ptr[3];
    J[3][43] =  +1*r_constant[40] +1*r_constant[41];
    J[3][45] =  +1*r_constant[43];
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
    J[17][17] =  -1*r_constant[9]*x_ptr[2] -1*r_constant[10]*x_ptr[14] -1*r_constant[11]*x_ptr[6] -1*r_constant[12]*x_ptr[3] -1*r_constant[44]*x_ptr[43] -1*r_constant[45]*x_ptr[38] -1*r_constant[46]*x_ptr[45];
    J[17][38] =  -1*r_constant[45]*x_ptr[17];
    J[17][43] =  -1*r_constant[44]*x_ptr[17];
    J[17][45] =  -1*r_constant[46]*x_ptr[17];
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
    J[24][24] =  -1*r_constant[16]*x_ptr[13] -1*r_constant[20]*x_ptr[26] -1*r_constant[22]*x_ptr[27] -1*r_constant[28]*x_ptr[32] -1*r_constant[30] -1*r_constant[61]*x_ptr[53];
    J[24][25] =  +1*r_constant[17] +1*r_constant[18] +1*r_constant[19];
    J[24][26] =  -1*r_constant[20]*x_ptr[24];
    J[24][27] =  -1*r_constant[22]*x_ptr[24];
    J[24][28] =  +1*r_constant[21] +1*r_constant[24];
    J[24][29] =  +1*r_constant[23] +1*r_constant[25];
    J[24][32] =  -1*r_constant[28]*x_ptr[24];
    J[24][33] =  +1*r_constant[29];
    J[24][53] =  -1*r_constant[61]*x_ptr[24];
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
    J[35][3] =  -1*r_constant[31]*x_ptr[35];
    J[35][35] =  -1*r_constant[31]*x_ptr[3];
    J[35][37] =  +1*r_constant[32];
    J[36][3] =  -1*r_constant[34]*x_ptr[36];
    J[36][36] =  -1*r_constant[34]*x_ptr[3] -2*r_constant[37]*x_ptr[36] -1*r_constant[60] * r_constant[49] * x_ptr[50]/ pow((r_constant[49] + x_ptr[36]),2.0);
    J[36][37] =  +1*r_constant[33];
    J[36][38] =  +1*r_constant[35];
    J[36][42] =  +1*r_constant[38];
    J[36][50] =  -1*r_constant[60] * x_ptr[36] / (r_constant[49] + x_ptr[36]);
    J[37][3] =  +1*r_constant[31]*x_ptr[35];
    J[37][35] =  +1*r_constant[31]*x_ptr[3];
    J[37][37] =  -1*r_constant[32] -1*r_constant[33];
    J[38][3] =  +1*r_constant[34]*x_ptr[36];
    J[38][17] =  -1*r_constant[45]*x_ptr[38];
    J[38][36] =  +1*r_constant[34]*x_ptr[3];
    J[38][38] =  -1*r_constant[35] -1*r_constant[36] -1*r_constant[45]*x_ptr[17];
    J[39][3] =  -1*r_constant[42]*x_ptr[39];
    J[39][38] =  +1*r_constant[36];
    J[39][39] =  -1*r_constant[42]*x_ptr[3] -1*r_constant[48]  * r_constant[49] * x_ptr[50]/ pow((r_constant[49] + x_ptr[39]),2.0);
    J[39][45] =  +1*r_constant[43];
    J[39][50] =  -1*r_constant[48] * x_ptr[39] / (r_constant[49] + x_ptr[39]);
    J[40][37] =  +1*r_constant[33];
    J[41][38] =  +1*r_constant[36];
    J[41][43] =  +1*r_constant[41];
    J[42][3] =  -1*r_constant[39]*x_ptr[42];
    J[42][36] =  +2*r_constant[37]*x_ptr[36];
    J[42][42] =  -1*r_constant[38] -1*r_constant[39]*x_ptr[3] -1*r_constant[55]-1*r_constant[60] * r_constant[49] * x_ptr[50]/ pow((r_constant[49] + x_ptr[42]),2.0);
    J[42][43] =  +1*r_constant[40];
    J[42][50] =  -1*r_constant[60] * x_ptr[42] / (r_constant[49] + x_ptr[42]);
    J[43][3] =  +1*r_constant[39]*x_ptr[42];
    J[43][17] =  -1*r_constant[44]*x_ptr[43];
    J[43][42] =  +1*r_constant[39]*x_ptr[3];
    J[43][43] =  -1*r_constant[40] -1*r_constant[41] -1*r_constant[44]*x_ptr[17];
    J[44][43] =  +1*r_constant[41];
    J[44][44] =  -1*r_constant[49]-1*r_constant[48]  * r_constant[49] * x_ptr[50]/ pow((r_constant[49] + x_ptr[44]),2.0);
    J[44][50] =  -1*r_constant[48] * x_ptr[44] / (r_constant[49] + x_ptr[39]);
    J[45][3] =  +1*r_constant[42]*x_ptr[39];
    J[45][17] =  -1*r_constant[46]*x_ptr[45];
    J[45][39] =  +1*r_constant[42]*x_ptr[3];
    J[45][45] =  -1*r_constant[43] -1*r_constant[46]*x_ptr[17];
    J[46][17] =  +1*r_constant[44]*x_ptr[43];
    J[46][43] =  +1*r_constant[44]*x_ptr[17];
    J[47][17] =  +1*r_constant[45]*x_ptr[38];
    J[47][38] =  +1*r_constant[45]*x_ptr[17];
    J[48][17] =  +1*r_constant[46]*x_ptr[45];
    J[48][45] =  +1*r_constant[46]*x_ptr[17];
    J[49][39] =  -1*r_constant[57]*x_ptr[55]*x_ptr[49]*(r_constant[54]*x_ptr[49]+r_constant[55]*r_constant[56])/pow(((x_ptr[39]+2*x_ptr[44])*r_constant[55]+r_constant[55]*r_constant[56]+(x_ptr[39]+2*x_ptr[44])*x_ptr[49]+r_constant[54]*x_ptr[49]),2.0);
    J[49][44] =  -2*r_constant[57]*x_ptr[55]*x_ptr[49]*(r_constant[54]*x_ptr[49]+r_constant[55]*r_constant[56])/pow(((x_ptr[39]+2*x_ptr[44])*r_constant[55]+r_constant[55]*r_constant[56]+(x_ptr[39]+2*x_ptr[44])*x_ptr[49]+r_constant[54]*x_ptr[49]),2.0);
    J[49][49] =  -1*r_constant[57]*x_ptr[55]*(x_ptr[39]+2*x_ptr[44])*r_constant[55]*(r_constant[56]+(x_ptr[39]+2*x_ptr[44]))/pow(((x_ptr[39]+2*x_ptr[44])*r_constant[55]+r_constant[55]*r_constant[56]+(x_ptr[39]+2*x_ptr[44])*x_ptr[49]+r_constant[54]*x_ptr[49]),2.0);
    J[49][55] =  -1*r_constant[57] * x_ptr[49]*(x_ptr[39]+2*x_ptr[44])/(r_constant[55]*((x_ptr[39]+2*x_ptr[44])+r_constant[56])+((x_ptr[39]+2*x_ptr[44])*x_ptr[50]+r_constant[54])*x_ptr[49]);
    J[50][39] =  +1*r_constant[57]*x_ptr[55]*x_ptr[49]*(r_constant[54]*x_ptr[49]+r_constant[55]*r_constant[56])/pow(((x_ptr[39]+2*x_ptr[44])*r_constant[55]+r_constant[55]*r_constant[56]+(x_ptr[39]+2*x_ptr[44])*x_ptr[49]+r_constant[54]*x_ptr[49]),2.0);
    J[50][44] =  +2*r_constant[57]*x_ptr[55]*x_ptr[49]*(r_constant[54]*x_ptr[49]+r_constant[55]*r_constant[56])/pow(((x_ptr[39]+2*x_ptr[44])*r_constant[55]+r_constant[55]*r_constant[56]+(x_ptr[39]+2*x_ptr[44])*x_ptr[49]+r_constant[54]*x_ptr[49]),2.0);
    J[50][49] =  -1*r_constant[57]*x_ptr[55]*(x_ptr[39]+2*x_ptr[44])*r_constant[55]*(r_constant[56]+(x_ptr[39]+2*x_ptr[44]))/pow(((x_ptr[39]+2*x_ptr[44])*r_constant[55]+r_constant[55]*r_constant[56]+(x_ptr[39]+2*x_ptr[44])*x_ptr[49]+r_constant[54]*x_ptr[49]),2.0);
    J[50][50] =  -1*r_constant[47]*x_ptr[51];
    J[50][51] =  -1*r_constant[47]*x_ptr[50];
    J[50][55] =  +1*r_constant[57] * x_ptr[49]*(x_ptr[39]+2*x_ptr[44])/(r_constant[55]*((x_ptr[39]+2*x_ptr[44])+r_constant[56])+((x_ptr[39]+2*x_ptr[44])*x_ptr[50]+r_constant[54])*x_ptr[49]);
    J[51][50] =  -1*r_constant[47]*x_ptr[51];
    J[51][51] =  -1*r_constant[47]*x_ptr[50];
    J[52][50] =  +1*r_constant[47]*x_ptr[51];
    J[52][51] =  +1*r_constant[47]*x_ptr[50];
    J[53][24] =  -1*r_constant[56]*x_ptr[53];
    J[53][39] =  +1*r_constant[53]*r_constant[51]*x_ptr[56]*r_constant[58]/pow(((x_ptr[39]+2*x_ptr[44])+r_constant[58]),2.0);
    J[53][44] =  +2*r_constant[53]*r_constant[51]*x_ptr[56]*r_constant[58]/pow(((x_ptr[39]+2*x_ptr[44])+r_constant[58]),2.0);
    J[53][53] =  -1*r_constant[52] * r_constant[50] * x_ptr[55] -1*r_constant[61] * x_ptr[24];
    J[53][55] =  -1*r_constant[52] * r_constant[50] * x_ptr[53];
    J[53][56] =  +1*r_constant[53] * r_constant[51] * (x_ptr[39]+2*x_ptr[44])/(r_constant[58] + (x_ptr[39]+2*x_ptr[44]));
    J[54][36] =  +1*r_constant[60] * r_constant[49] * x_ptr[50]/ pow((r_constant[49] + x_ptr[36]),2.0);
    J[54][39] =  +1*r_constant[48]  * r_constant[49] * x_ptr[50]/ pow((r_constant[49] + x_ptr[39]),2.0);
    J[54][42] =  +2*r_constant[60] * r_constant[49] * x_ptr[50]/ pow((r_constant[49] + x_ptr[42]),2.0);
    J[54][44] =  +2*r_constant[48]  * r_constant[49] * x_ptr[50]/ pow((r_constant[49] + x_ptr[44]),2.0);
    J[54][50] =  +1*r_constant[60] * x_ptr[36] / (r_constant[49] + x_ptr[36]) + 1*r_constant[48] * x_ptr[39] / (r_constant[49] + x_ptr[39]) + 1*r_constant[60] * x_ptr[42] / (r_constant[49] + x_ptr[42]) + 1*r_constant[48] * x_ptr[44] / (r_constant[49] + x_ptr[39]);
    J[55][39] =  +1*r_constant[53]*r_constant[51]*x_ptr[56]*r_constant[58]/pow(((x_ptr[39]+2*x_ptr[44])+r_constant[58]),2.0);
    J[55][44] =  +2*r_constant[53]*r_constant[51]*x_ptr[56]*r_constant[58]/pow(((x_ptr[39]+2*x_ptr[44])+r_constant[58]),2.0);
    J[55][53] =  -1*r_constant[52] * r_constant[50] * x_ptr[55];
    J[55][55] =  -1*r_constant[52] * r_constant[50] * x_ptr[53] -1*r_constant[59];
    J[55][56] =  +1*r_constant[53] * r_constant[51] * (x_ptr[39]+2*x_ptr[44])/(r_constant[58] + (x_ptr[39]+2*x_ptr[44]));
    J[56][39] =  -1*r_constant[53]*r_constant[51]*x_ptr[56]*r_constant[58]/pow(((x_ptr[39]+2*x_ptr[44])+r_constant[58]),2.0);
    J[56][44] =  -2*r_constant[53]*r_constant[51]*x_ptr[56]*r_constant[58]/pow(((x_ptr[39]+2*x_ptr[44])+r_constant[58]),2.0);
    J[56][53] =  +1*r_constant[52] * r_constant[50] * x_ptr[55];
    J[56][55] =  +1*r_constant[52] * r_constant[50] * x_ptr[53];
    J[56][56] =  -1*r_constant[53] * r_constant[51] * (x_ptr[39]+2*x_ptr[44])/(r_constant[58] + (x_ptr[39]+2*x_ptr[44]));

/*
 * End of user-generated jacobian function
 */

    PetscFunctionReturn(0);
}

#endif
