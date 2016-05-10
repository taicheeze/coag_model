#ifndef _INITIAL_CONCENTRATIONS_H_
#define _INITIAL_CONCENTRATIONS_H_

#include <petscsys.h>

#define NUMBER_OF_SPECIES_IN_VOLUME     67
#define NUMBER_OF_SPECIES_ON_SURFACE    69

const char VOLUME_SPECIES_NAMES[NUMBER_OF_SPECIES_IN_VOLUME][256]
        = {"VII"                   ,"VIIa"                  ,"Xa"                    ,"IIa"                   ,"X"                     ,
           "IX"                    ,"IXa"                   ,"II"                    ,"VIII"                  ,"VIIIa"                 ,
           "VIIIa1L"               ,"VIIIa2"                ,"V"                     ,"Va"                    ,"mIIa"                  ,
           "TFPI"                  ,"XaTFPI"                ,"ATIII"                 ,"XaATIII"               ,"mIIaATIII"             ,
           "IXaATIII"              ,"IIaATIII"              ,"substrate"             ,"fluorescence"          ,"APC"                   ,
           "APCVa"                 ,"Va5"                   ,"Va3"                   ,"APCVa5"                ,"APCVa3"                ,
           "Va53"                  ,"HCF"                   ,"LCA1"                  ,"APCLCA1"               ,"PC"                    ,
           "XI"                    ,"XIa"                   ,"IIaXI"                 ,"LPS"                   ,"IL1a"                  ,
           "IL1b"                  ,"TNFa"                  ,"IL6"                   ,"AI"                    ,"ATF"                   ,
           "Fg"                    ,"FnI"                   ,"Fg_IIa"                ,"FnI_IIa"               ,"FnII"                  ,
           "FPA"                   ,"FPB"                   ,"FnI2"                  ,"FnI2_IIa"              ,"FnII2"                 ,
           "FnII_IIa"              ,"FnI2_IIa_ATIII"        ,"FnI_IIa_ATIII"         ,"FnII_IIa_ATIII"        ,"Pg"                    ,
           "Pn"                    ,"AP"                    ,"Pn_AP"                 ,"PAI"                   ,"FDP"                   ,
           "stPA"                  ,"stPA_PAI"              
          };


const PetscReal VOLUME_INITIAL_CONCENTRATIONS[NUMBER_OF_SPECIES_IN_VOLUME]
        = {  1.0000e-05,  1.0000e-07,  0.0000e+00,  0.0000e+00,  1.6000e-04,
             9.0000e-05,  0.0000e+00,  1.4000e-03,  7.0000e-07,  0.0000e+00,
             0.0000e+00,  0.0000e+00,  2.0000e-05,  0.0000e+00,  0.0000e+00,
             2.5000e-06,  0.0000e+00,  3.4000e-03,  0.0000e+00,  0.0000e+00,
             0.0000e+00,  0.0000e+00,  1.0000e-03,  0.0000e+00,  0.0000e+00,
             0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,
             0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  6.5000e-05,
             0.0000e+00,  0.0000e+00,  0.0000e+00,  6.6000e-04,  0.0000e+00,
             0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,
             9.0000e-03,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,
             0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,
             0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  2.0000e-03,
             0.0000e+00,  1.0000e-03,  0.0000e+00,  4.0000e-07,  0.0000e+00,
             0.0000e+00,  7.0000e-08
           };


const char SURFACE_SPECIES_NAMES[NUMBER_OF_SPECIES_ON_SURFACE][256]
        = {"TF"                    ,"TFVII"                 ,"TFVIIa"                ,"TFVIIaX"               ,"TFVIIaXa"              ,
           "TFVIIaIX"              ,"TFVIIaXaTFPI"          ,"TFVIIaATIII"           ,"P9"                    ,"IXP9"                  ,
           "IXaP9"                 ,"P9s"                   ,"IXaP9s"                ,"P10"                   ,"XP10"                  ,
           "XaP10"                 ,"P5"                    ,"VP5"                   ,"VaP5"                  ,"P8"                    ,
           "VIIIP8"                ,"VIIIaP8"               ,"P2"                    ,"IIP2"                  ,"IIaP2"                 ,
           "VP5XaP10"              ,"VIIIP8XaP10"           ,"IXaP9VIIIaP8"          ,"IXaP9VIIIaP8XP10"      ,"IXaP9sVIIIaP8"         ,
           "IXaP9sVIIIaP8XP10"     ,"XaP10VaP5"             ,"XaP10VaP5IIP2"         ,"mIIaP2"                ,"VIIIa1LP8"             ,
           "IXaP9VIIIaP8X"         ,"XaP10VaP5II"           ,"IXaP9sVIIIaP8X"        ,"Va3P5"                 ,"Va5P5"                 ,
           "Va53P5"                ,"XaP10Va3P5"            ,"XaP10Va5P5"            ,"XaP10Va5P5IIP2"        ,"XaP10Va3P5IIP2"        ,
           "XaP10Va53P5"           ,"XaP10Va53P5IIP2"       ,"IIP2VaP5"              ,"TM"                    ,"TMIIa"                 ,
           "TMmIIa"                ,"TMIIaPC"               ,"TMmIIaPC"              ,"TMIIaAPC"              ,"IXP9XIa"               ,
           "MRL"                   ,"MR1a"                  ,"MR1b"                  ,"MRT"                   ,"MR6"                   ,
           "MRL_LPS"               ,"MR1a_IL1a"             ,"MR1b_IL1b"             ,"MRT_TNFa"              ,"MR6_IL6"               ,
           "tPA"                   ,"tPA_PAI"               ,"ECPR"                  ,"ECPR_Pg"               
          };


const PetscReal SURFACE_INITIAL_CONCENTRATIONS_PATCHES[NUMBER_OF_PATCHES][NUMBER_OF_SPECIES_ON_SURFACE]
        = {
           {  6.2500e-13,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,
              0.0000e+00,  0.0000e+00,  0.0000e+00,  1.0379e-09,  0.0000e+00,
              0.0000e+00,  1.0379e-09,  0.0000e+00,  1.1209e-08,  0.0000e+00,
              0.0000e+00,  1.2454e-08,  0.0000e+00,  0.0000e+00,  1.8682e-09,
              0.0000e+00,  0.0000e+00,  8.3029e-09,  0.0000e+00,  0.0000e+00,
              0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,
              0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,
              0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,
              0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,
              0.0000e+00,  0.0000e+00,  0.0000e+00,  1.0000e-12,  0.0000e+00,
              0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,
              2.0000e-06,  2.0000e-06,  2.0000e-06,  2.0000e-06,  2.0000e-06,
              0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,
              7.0000e-11,  0.0000e+00,  2.7000e-09,  0.0000e+00
           }
          };


#define NUMBER_OF_REACTIONS_IN_VOLUME    66
#define NUMBER_OF_REACTIONS_ON_SURFACE   172

#endif
