
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
 *  30.  XI + IIa -> IIaXI                                                               ([35] + [3] -> [37])
 *  31.  IIaXI -> XI + IIa                                                               ([37] -> [35] + [3])
 *  32.  IIaXI -> IIa + XIa                                                              ([37] -> [3] + [36])
 *  33.  APC -> NULL                                                                     ([24] -> )
 *  34.  AI + LPS -> AI                                                                  ([43] + [38] -> [43])
 *  35.  IL1b + AI -> AI                                                                 ([40] + [43] -> [43])
 *  36.  TNFa + AI -> AI                                                                 ([41] + [43] -> [43])
 *  37.  IL6 + AI -> AI                                                                  ([42] + [43] -> [43])
 *  38.  IL1a + AI -> AI                                                                 ([39] + [43] -> [43])
 *  39.  AI -> NULL                                                                      ([43] -> )
 *  40.  Fg + IIa -> Fg_IIa                                                              ([45] + [3] -> [47])
 *  41.  Fg_IIa -> Fg + IIa                                                              ([47] -> [45] + [3])
 *  42.  Fg_IIa -> FnI + IIa + FPA                                                       ([47] -> [46] + [3] + [50])
 *  43.  FnI + IIa -> FnI_IIa                                                            ([46] + [3] -> [48])
 *  44.  FnI_IIa -> FnI + IIa                                                            ([48] -> [46] + [3])
 *  45.  FnI_IIa -> FnII + IIa + FPB                                                     ([48] -> [49] + [3] + [51])
 *  46.  FnI -> FnI2                                                                     ([46] -> [52])
 *  47.  FnI2 -> FnI                                                                     ([52] -> [46])
 *  48.  FnI2 + IIa -> FnI2_IIa                                                          ([52] + [3] -> [53])
 *  49.  FnI2_IIa -> FnI2 + IIa                                                          ([53] -> [52] + [3])
 *  50.  FnI2_IIa -> FnII2 + IIa + FPB                                                   ([53] -> [54] + [3] + [51])
 *  51.  FnII + IIa -> FnII_IIa                                                          ([49] + [3] -> [55])
 *  52.  FnII_IIa -> FnII + IIa                                                          ([55] -> [49] + [3])
 *  53.  FnI2_IIa + ATIII -> FnI2_IIa_ATIII                                              ([53] + [17] -> [56])
 *  54.  FnI_IIa + ATIII -> FnI_IIa_ATIII                                                ([48] + [17] -> [57])
 *  55.  FnII_IIa + ATIII -> FnII_IIa_ATIII                                              ([55] + [17] -> [58])
 *  56.  Pn + AP -> Pn_AP                                                                ([60] + [61] -> [62])
 *  57.  FnII -> FDP                                                                     ([49] -> [64])
 *  58.  FnII2 -> FDP                                                                    ([54] -> [64])
 *  59.  stPA + PAI -> stPA_PAI                                                          ([65] + [63] -> [66])
 *  60.  stPA_PAI -> stPA + PAI                                                          ([66] -> [65] + [63])
 *  61.  Pg -> Pn                                                                        ([59] -> [60])
 *  62.  stPA -> NULL                                                                    ([65] -> )
 */


/*
 * volume_stoichiometry: stoichiometry matrix of volume species for volume reactions
 */
/*
 * Begin of user-generated volume reaction stoichiometries
 */
PetscInt  volume_stoichiometry[NUMBER_OF_SPECIES_IN_VOLUME][NUMBER_OF_REACTIONS_IN_VOLUME]
          ={
            {-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0,-1, 1,-1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 1, 0, 0, 0, 0, 0, 0, 0,-1, 1, 1,-1, 1, 1, 0, 0,-1, 1, 1,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 1,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 1, 1,-1, 1,-1, 1, 1, 1, 0, 0,-1, 1, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,-1, 1, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 1, 0,-2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1,-1, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 1, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 0,-1},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0}
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
    r_constant[30] = 1.0000e+05;
    r_constant[31] = 1.0000e+01;
    r_constant[32] = 1.4300e+00;
    r_constant[33] = 1.1111e-03;
    r_constant[34] = 2.2219e+01;
    r_constant[35] = 4.1526e+01;
    r_constant[36] = 1.6975e+02;
    r_constant[37] = 8.2750e+02;
    r_constant[38] = 7.8502e+02;
    r_constant[39] = 1.0286e-03;
    r_constant[40] = 1.0000e+05;
    r_constant[41] = 7.2000e+02;
    r_constant[42] = 8.4000e+01;
    r_constant[43] = 1.0000e+05;
    r_constant[44] = 7.5000e+02;
    r_constant[45] = 7.4000e+00;
    r_constant[46] = 1.0000e+03;
    r_constant[47] = 6.4000e-02;
    r_constant[48] = 1.0000e+05;
    r_constant[49] = 7.5000e+02;
    r_constant[50] = 4.9000e+01;
    r_constant[51] = 1.0000e+05;
    r_constant[52] = 1.0000e+03;
    r_constant[53] = 1.6000e+01;
    r_constant[54] = 1.6000e+01;
    r_constant[55] = 1.0000e+01; 
    r_constant[56] = 3.0000e+03; //AP to Pn binding
    r_constant[57] = 4.7000e-01; //fibrin degradation to FDP formation, increased from 4.7e-1
    r_constant[58] = 2.1000e-03; //fibrin degradation to FDP formation
    r_constant[59] = 4.0000e+04; //pai stpa binding
    r_constant[60] = 1.6000e-02; //increased from 1.6e-2, stpa_pai unbinding
    r_constant[61] = 1.0000e+00; //pai stpa binding
    r_constant[62] = 1.0000e+00; //stpa_pai unbinding
    r_constant[63] = 7.7000e-05; //Pg to pn conversion
    r_constant[64] = 4.1000e-04; //Pg to pn conversion
    r_constant[65] = 3.0000e-04; //Pg to pn conversion
    r_constant[66] = 9.0000e-02; //should be 9e-2
    r_constant[67] = 1.0000e-03; //stpa_pai unbinding
    r_constant[68] = 2.8900e-03; //stpa hepatic clearance

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
    reaction_rate[30] = r_constant[30] * x_ptr[3] * x_ptr[35];
    reaction_rate[31] = r_constant[31] * x_ptr[37];
    reaction_rate[32] = r_constant[32] * x_ptr[37];
    reaction_rate[33] = r_constant[33] * x_ptr[24];
    reaction_rate[34] = r_constant[34] * x_ptr[38] * x_ptr[43];
    reaction_rate[35] = r_constant[35] * x_ptr[40] * x_ptr[43];
    reaction_rate[36] = r_constant[36] * x_ptr[41] * x_ptr[43];
    reaction_rate[37] = r_constant[37] * x_ptr[42] * x_ptr[43];
    reaction_rate[38] = r_constant[38] * x_ptr[39] * x_ptr[43];
    reaction_rate[39] = r_constant[39] * x_ptr[43];
    reaction_rate[40] = r_constant[40] * x_ptr[3] * x_ptr[45];
    reaction_rate[41] = r_constant[41] * x_ptr[47];
    reaction_rate[42] = r_constant[42] * x_ptr[47];
    reaction_rate[43] = r_constant[43] * x_ptr[3] * x_ptr[46];
    reaction_rate[44] = r_constant[44] * x_ptr[48];
    reaction_rate[45] = r_constant[45] * x_ptr[48];
    reaction_rate[46] = r_constant[46] * x_ptr[46] * x_ptr[46];
    reaction_rate[47] = r_constant[47] * x_ptr[52];
    reaction_rate[48] = r_constant[48] * x_ptr[3] * x_ptr[52];
    reaction_rate[49] = r_constant[49] * x_ptr[53];
    reaction_rate[50] = r_constant[50] * x_ptr[53];
    reaction_rate[51] = r_constant[51] * x_ptr[3] * x_ptr[49];
    reaction_rate[52] = r_constant[52] * x_ptr[55];
    reaction_rate[53] = r_constant[53] * x_ptr[17] * x_ptr[53];
    reaction_rate[54] = r_constant[54] * x_ptr[17] * x_ptr[48];
    reaction_rate[55] = r_constant[55] * x_ptr[17] * x_ptr[55];
    reaction_rate[56] = r_constant[56] * x_ptr[60] * x_ptr[61];
    reaction_rate[57] = r_constant[57] * x_ptr[60] * x_ptr[49] / (r_constant[58] + x_ptr[49]);
    reaction_rate[58] = r_constant[57] * x_ptr[60] * x_ptr[54] / (r_constant[58] + x_ptr[54]);
    reaction_rate[59] = r_constant[61]*r_constant[59] * x_ptr[63] * x_ptr[65];
    reaction_rate[60] = r_constant[62]*r_constant[60] * x_ptr[66]* (x_ptr[49]+2*x_ptr[54])/(r_constant[67] + (x_ptr[49]+2*x_ptr[54]));
    reaction_rate[61] = r_constant[66] *x_ptr[65]* x_ptr[59]*(x_ptr[49]+2*x_ptr[54])/ (r_constant[64]*(x_ptr[49]+2*x_ptr[54])+r_constant[64]*r_constant[65]+(x_ptr[49]+2*x_ptr[54])*x_ptr[59]+r_constant[63]*x_ptr[59]);
    reaction_rate[62] = r_constant[68] * x_ptr[65];

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
    r_constant[30] = 1.0000e+05;
    r_constant[31] = 1.0000e+01;
    r_constant[32] = 1.4300e+00;
    r_constant[33] = 1.1111e-03;
    r_constant[34] = 2.2219e+01;
    r_constant[35] = 4.1526e+01;
    r_constant[36] = 1.6975e+02;
    r_constant[37] = 8.2750e+02;
    r_constant[38] = 7.8502e+02;
    r_constant[39] = 1.0286e-03;
    r_constant[40] = 1.0000e+05;
    r_constant[41] = 7.2000e+02;
    r_constant[42] = 8.4000e+01;
    r_constant[43] = 1.0000e+05;
    r_constant[44] = 7.5000e+02;
    r_constant[45] = 7.4000e+00;
    r_constant[46] = 1.0000e+03;
    r_constant[47] = 6.4000e-02;
    r_constant[48] = 1.0000e+05;
    r_constant[49] = 7.5000e+02;
    r_constant[50] = 4.9000e+01;
    r_constant[51] = 1.0000e+05;
    r_constant[52] = 1.0000e+03;
    r_constant[53] = 1.6000e+01;
    r_constant[54] = 1.6000e+01;
    r_constant[55] = 1.0000e+01; 
    r_constant[56] = 3.0000e+03; //AP to Pn binding
    r_constant[57] = 4.7000e-01; //fibrin degradation to FDP formation, increased from 4.7e-1
    r_constant[58] = 2.1000e-03; //fibrin degradation to FDP formation
    r_constant[59] = 4.0000e+04; //pai stpa binding
    r_constant[60] = 1.6000e-02; //increased from 1.6e-2, stpa_pai unbinding
    r_constant[61] = 1.0000e+00; //pai stpa binding
    r_constant[62] = 1.0000e+00; //stpa_pai unbinding
    r_constant[63] = 7.7000e-05; //Pg to pn conversion
    r_constant[64] = 4.1000e-04; //Pg to pn conversion
    r_constant[65] = 3.0000e-04; //Pg to pn conversion
    r_constant[66] = 9.0000e-02; //should be 9e-2
    r_constant[67] = 1.0000e-03; //stpa_pai unbinding
    r_constant[68] = 2.8900e-03; //stpa hepatic clearance

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
    J[3][3] =  -1*r_constant[12]*x_ptr[17] -1*r_constant[30]*x_ptr[35] -1*r_constant[40]*x_ptr[45] -1*r_constant[43]*x_ptr[46] -1*r_constant[48]*x_ptr[52] -1*r_constant[51]*x_ptr[49];
    J[3][7] =  +1*r_constant[2]*x_ptr[2];
    J[3][17] =  -1*r_constant[12]*x_ptr[3];
    J[3][35] =  -1*r_constant[30]*x_ptr[3];
    J[3][37] =  +1*r_constant[31] +1*r_constant[32];
    J[3][45] =  -1*r_constant[40]*x_ptr[3];
    J[3][46] =  -1*r_constant[43]*x_ptr[3];
    J[3][47] =  +1*r_constant[41] +1*r_constant[42];
    J[3][48] =  +1*r_constant[44] +1*r_constant[45];
    J[3][49] =  -1*r_constant[51]*x_ptr[3];
    J[3][52] =  -1*r_constant[48]*x_ptr[3];
    J[3][53] =  +1*r_constant[49] +1*r_constant[50];
    J[3][55] =  +1*r_constant[52];
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
    J[17][17] =  -1*r_constant[9]*x_ptr[2] -1*r_constant[10]*x_ptr[14] -1*r_constant[11]*x_ptr[6] -1*r_constant[12]*x_ptr[3] -1*r_constant[53]*x_ptr[53] -1*r_constant[54]*x_ptr[48] -1*r_constant[55]*x_ptr[55];
    J[17][48] =  -1*r_constant[54]*x_ptr[17];
    J[17][53] =  -1*r_constant[53]*x_ptr[17];
    J[17][55] =  -1*r_constant[55]*x_ptr[17];
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
    J[24][24] =  -1*r_constant[16]*x_ptr[13] -1*r_constant[20]*x_ptr[26] -1*r_constant[22]*x_ptr[27] -1*r_constant[28]*x_ptr[32] -1*r_constant[33];
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
    J[35][3] =  -1*r_constant[30]*x_ptr[35];
    J[35][35] =  -1*r_constant[30]*x_ptr[3];
    J[35][37] =  +1*r_constant[31];
    J[36][37] =  +1*r_constant[32];
    J[37][3] =  +1*r_constant[30]*x_ptr[35];
    J[37][35] =  +1*r_constant[30]*x_ptr[3];
    J[37][37] =  -1*r_constant[31] -1*r_constant[32];
    J[38][38] =  -1*r_constant[34]*x_ptr[43];
    J[38][43] =  -1*r_constant[34]*x_ptr[38];
    J[39][39] =  -1*r_constant[38]*x_ptr[43];
    J[39][43] =  -1*r_constant[38]*x_ptr[39];
    J[40][40] =  -1*r_constant[35]*x_ptr[43];
    J[40][43] =  -1*r_constant[35]*x_ptr[40];
    J[41][41] =  -1*r_constant[36]*x_ptr[43];
    J[41][43] =  -1*r_constant[36]*x_ptr[41];
    J[42][42] =  -1*r_constant[37]*x_ptr[43];
    J[42][43] =  -1*r_constant[37]*x_ptr[42];
    J[43][43] =  -1*r_constant[39];
    J[45][3] =  -1*r_constant[40]*x_ptr[45];
    J[45][45] =  -1*r_constant[40]*x_ptr[3];
    J[45][47] =  +1*r_constant[41];
    J[46][3] =  -1*r_constant[43]*x_ptr[46];
    J[46][46] =  -1*r_constant[43]*x_ptr[3] -2*r_constant[46] * x_ptr[46];
    J[46][47] =  +1*r_constant[42];
    J[46][48] =  +1*r_constant[44];
    J[46][52] =  +1*r_constant[47];
    J[47][3] =  +1*r_constant[40]*x_ptr[45];
    J[47][45] =  +1*r_constant[40]*x_ptr[3];
    J[47][47] =  -1*r_constant[41] -1*r_constant[42];
    J[48][3] =  +1*r_constant[43]*x_ptr[46];
    J[48][17] =  -1*r_constant[54]*x_ptr[48];
    J[48][46] =  +1*r_constant[43]*x_ptr[3];
    J[48][48] =  -1*r_constant[44] -1*r_constant[45] -1*r_constant[54]*x_ptr[17];
    J[49][3] =  -1*r_constant[51]*x_ptr[49];
    J[49][48] =  +1*r_constant[45];
    J[49][49] =  -1*r_constant[57]*r_constant[58]*x_ptr[60] / ((r_constant[58] + x_ptr[49]) * (r_constant[58] + x_ptr[49]))-1*r_constant[51] * x_ptr[3];
    J[49][55] =  +1*r_constant[52];
    J[49][60] =  -1*r_constant[57]*x_ptr[49]/ (r_constant[58] + x_ptr[49]);
    J[50][47] =  +1*r_constant[42];
    J[51][48] =  +1*r_constant[45];
    J[51][53] =  +1*r_constant[50];
    J[52][3] =  -1*r_constant[48]*x_ptr[52];
    J[52][46] =  +2*r_constant[46] * x_ptr[46];
    J[52][52] =  -1*r_constant[47] -1*r_constant[48]*x_ptr[3];
    J[52][53] =  +1*r_constant[49];
    J[53][3] =  +1*r_constant[48]*x_ptr[52];
    J[53][17] =  -1*r_constant[53]*x_ptr[53];
    J[53][52] =  +1*r_constant[48]*x_ptr[3];
    J[53][53] =  -1*r_constant[49] -1*r_constant[50] -1*r_constant[53]*x_ptr[17];
    J[54][53] =  +1*r_constant[50];
    J[54][54] =  -1*r_constant[57]*r_constant[58]*x_ptr[60] / ((r_constant[58] + x_ptr[54]) * (r_constant[58] + x_ptr[54]));
    J[54][60] =  -1*r_constant[57]*x_ptr[54]/ (r_constant[58] + x_ptr[54]);
    J[55][3] =  +1*r_constant[51]*x_ptr[49];
    J[55][17] =  -1*r_constant[55]*x_ptr[55];
    J[55][49] =  +1*r_constant[51]*x_ptr[3];
    J[55][55] =  -1*r_constant[52] -1*r_constant[55]*x_ptr[17];
    J[56][17] =  +1*r_constant[53]*x_ptr[53];
    J[56][53] =  +1*r_constant[53]*x_ptr[17];
    J[57][17] =  +1*r_constant[54]*x_ptr[48];
    J[57][48] =  +1*r_constant[54]*x_ptr[17];
    J[58][17] =  +1*r_constant[55]*x_ptr[55];
    J[58][55] =  +1*r_constant[55]*x_ptr[17];
    J[59][49] =  -1*r_constant[66]*x_ptr[65]*x_ptr[59]*(r_constant[63]*x_ptr[59]+r_constant[64]*r_constant[65])/(((x_ptr[49]+2*x_ptr[54])*r_constant[64]+r_constant[64]*r_constant[65]+(x_ptr[49]+2*x_ptr[54])*x_ptr[59]+r_constant[63]*x_ptr[59])*((x_ptr[49]+2*x_ptr[54])*r_constant[64]+r_constant[64]*r_constant[65]+(x_ptr[49]+2*x_ptr[54])*x_ptr[59]+r_constant[63]*x_ptr[59]));
    J[59][54] =  -2*r_constant[66]*x_ptr[65]*x_ptr[59]*(r_constant[63]*x_ptr[59]+r_constant[64]*r_constant[65])/(((x_ptr[49]+2*x_ptr[54])*r_constant[64]+r_constant[64]*r_constant[65]+(x_ptr[49]+2*x_ptr[54])*x_ptr[59]+r_constant[63]*x_ptr[59])*((x_ptr[49]+2*x_ptr[54])*r_constant[64]+r_constant[64]*r_constant[65]+(x_ptr[49]+2*x_ptr[54])*x_ptr[59]+r_constant[63]*x_ptr[59]));
    J[59][59] =  -1*r_constant[66]*x_ptr[65]*(x_ptr[49]+2*x_ptr[54])*r_constant[64]*(r_constant[65]+(x_ptr[49]+2*x_ptr[54]))/(((x_ptr[49]+2*x_ptr[54])*r_constant[64]+r_constant[64]*r_constant[65]+(x_ptr[49]+2*x_ptr[54])*x_ptr[59]+r_constant[63]*x_ptr[59])*((x_ptr[49]+2*x_ptr[54])*r_constant[64]+r_constant[64]*r_constant[65]+(x_ptr[49]+2*x_ptr[54])*x_ptr[59]+r_constant[63]*x_ptr[59]));
    J[59][65] =  -1*r_constant[66]*x_ptr[59]*(x_ptr[49]+2*x_ptr[54])/(((x_ptr[49]+2*x_ptr[54])+r_constant[63])*(x_ptr[59]+r_constant[64]*((x_ptr[49]+2*x_ptr[54])+r_constant[65]))/((x_ptr[49]+2*x_ptr[54])+r_constant[63]));
    J[60][49] =  +1*r_constant[66]*x_ptr[65]*x_ptr[59]*(r_constant[63]*x_ptr[59]+r_constant[64]*r_constant[65])/(((x_ptr[49]+2*x_ptr[54])*r_constant[64]+r_constant[64]*r_constant[65]+(x_ptr[49]+2*x_ptr[54])*x_ptr[59]+r_constant[63]*x_ptr[59])*((x_ptr[49]+2*x_ptr[54])*r_constant[64]+r_constant[64]*r_constant[65]+(x_ptr[49]+2*x_ptr[54])*x_ptr[59]+r_constant[63]*x_ptr[59]));
    J[60][54] =  +2*r_constant[66]*x_ptr[65]*x_ptr[59]*(r_constant[63]*x_ptr[59]+r_constant[64]*r_constant[65])/(((x_ptr[49]+2*x_ptr[54])*r_constant[64]+r_constant[64]*r_constant[65]+(x_ptr[49]+2*x_ptr[54])*x_ptr[59]+r_constant[63]*x_ptr[59])*((x_ptr[49]+2*x_ptr[54])*r_constant[64]+r_constant[64]*r_constant[65]+(x_ptr[49]+2*x_ptr[54])*x_ptr[59]+r_constant[63]*x_ptr[59]));
    J[60][59] =  +1*r_constant[66]*x_ptr[65]*(x_ptr[49]+2*x_ptr[54])*r_constant[64]*(r_constant[65]+(x_ptr[49]+2*x_ptr[54]))/(((x_ptr[49]+2*x_ptr[54])*r_constant[64]+r_constant[64]*r_constant[65]+(x_ptr[49]+2*x_ptr[54])*x_ptr[59]+r_constant[63]*x_ptr[59])*((x_ptr[49]+2*x_ptr[54])*r_constant[64]+r_constant[64]*r_constant[65]+(x_ptr[49]+2*x_ptr[54])*x_ptr[59]+r_constant[63]*x_ptr[59]));
    J[60][60] =  -1*r_constant[56]*x_ptr[61];
    J[60][61] =  -1*r_constant[56]*x_ptr[60];
    J[60][65] =  +1*r_constant[66]*x_ptr[59]*(x_ptr[49]+2*x_ptr[54])/(((x_ptr[49]+2*x_ptr[54])+r_constant[63])*(x_ptr[59]+r_constant[64]*((x_ptr[49]+2*x_ptr[54])+r_constant[65]))/((x_ptr[49]+2*x_ptr[54])+r_constant[63]));
    J[61][60] =  -1*r_constant[56]*x_ptr[61];
    J[61][61] =  -1*r_constant[56]*x_ptr[60];
    J[62][60] =  +1*r_constant[56]*x_ptr[61];
    J[62][61] =  +1*r_constant[56]*x_ptr[60];
    J[63][49] =  +1*r_constant[62]*r_constant[60]*x_ptr[66]*r_constant[67]/(((x_ptr[49]+2*x_ptr[54])+r_constant[67])*((x_ptr[49]+2*x_ptr[54])+r_constant[67]));
    J[63][54] =  +2*r_constant[62]*r_constant[60]*x_ptr[66]*r_constant[67]/(((x_ptr[49]+2*x_ptr[54])+r_constant[67])*((x_ptr[49]+2*x_ptr[54])+r_constant[67]));
    J[63][63] =  -1*r_constant[61]*r_constant[59]*x_ptr[65];
    J[63][65] =  -1*r_constant[61]*r_constant[59]*x_ptr[63]-1*r_constant[68];
    J[63][66] =  +1*r_constant[62]*r_constant[60]*(x_ptr[49]+2*x_ptr[54])/(r_constant[67]+(x_ptr[49]+2*x_ptr[54]));
    J[64][49] =  +1*r_constant[57]*r_constant[58]*x_ptr[60] / ((r_constant[58] + x_ptr[49]) * (r_constant[58] + x_ptr[49]));
    J[64][54] =  +2*r_constant[57]*r_constant[58]*x_ptr[60] / ((r_constant[58] + x_ptr[54]) * (r_constant[58] + x_ptr[54]));
    J[64][60] =  +1*r_constant[57]*x_ptr[49]/ (r_constant[58] + x_ptr[49]) +2*r_constant[57]*x_ptr[54]/ (r_constant[58] + x_ptr[54]);
    J[65][49] =  +1*r_constant[62]*r_constant[60]*x_ptr[66]*r_constant[67]/(((x_ptr[49]+2*x_ptr[54])+r_constant[67])*((x_ptr[49]+2*x_ptr[54])+r_constant[67]));
    J[65][54] =  +2*r_constant[62]*r_constant[60]*x_ptr[66]*r_constant[67]/(((x_ptr[49]+2*x_ptr[54])+r_constant[67])*((x_ptr[49]+2*x_ptr[54])+r_constant[67]));
    J[65][63] =  -1*r_constant[61]*r_constant[59]*x_ptr[65];
    J[65][65] =  -1*r_constant[61]*r_constant[59]*x_ptr[63];
    J[65][66] =  +1*r_constant[62]*r_constant[60]*(x_ptr[49]+2*x_ptr[54])/(r_constant[67]+(x_ptr[49]+2*x_ptr[54]));
    J[66][49] =  -1*r_constant[62]*r_constant[60]*x_ptr[66]*r_constant[67]/(((x_ptr[49]+2*x_ptr[54])+r_constant[67])*((x_ptr[49]+2*x_ptr[54])+r_constant[67]));
    J[66][54] =  -2*r_constant[62]*r_constant[60]*x_ptr[66]*r_constant[67]/(((x_ptr[49]+2*x_ptr[54])+r_constant[67])*((x_ptr[49]+2*x_ptr[54])+r_constant[67]));
    J[66][63] =  +1*r_constant[61]*r_constant[59]*x_ptr[65];
    J[66][65] =  +1*r_constant[61]*r_constant[59]*x_ptr[63];
    J[66][66] =  -1*r_constant[62]*r_constant[60]*(x_ptr[49]+2*x_ptr[54])/(r_constant[67]+(x_ptr[49]+2*x_ptr[54]));

/*
 * End of user-generated jacobian function
 */

    PetscFunctionReturn(0);
}

#endif
