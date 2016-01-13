#ifndef _INITIAL_CONCENTRATIONS_H_
#define _INITIAL_CONCENTRATIONS_H_

#include <petscsys.h>

#define NUMBER_OF_SPECIES_IN_VOLUME     35
#define NUMBER_OF_SPECIES_ON_SURFACE    54

const char VOLUME_SPECIES_NAMES[NUMBER_OF_SPECIES_IN_VOLUME][256]
        = {"VII"                   ,"VIIa"                  ,"Xa"                    ,"IIa"                   ,"X"                     ,
           "IX"                    ,"IXa"                   ,"II"                    ,"VIII"                  ,"VIIIa"                 ,
           "VIIIa1L"               ,"VIIIa2"                ,"V"                     ,"Va"                    ,"mIIa"                  ,
           "TFPI"                  ,"XaTFPI"                ,"ATIII"                 ,"XaATIII"               ,"mIIaATIII"             ,
           "IXaATIII"              ,"IIaATIII"              ,"substrate"             ,"fluorescence"          ,"APC"                   ,
           "APCVa"                 ,"Va5"                   ,"Va3"                   ,"APCVa5"                ,"APCVa3"                ,
           "Va53"                  ,"HCF"                   ,"LCA1"                  ,"APCLCA1"               ,"PC"                    
          };


const PetscReal VOLUME_INITIAL_CONCENTRATIONS[NUMBER_OF_SPECIES_IN_VOLUME]
        = {  1.0000e-05,  1.0000e-07,  0.0000e+00,  0.0000e+00,  1.6000e-04,
             9.0000e-05,  0.0000e+00,  1.4000e-03,  7.0000e-07,  0.0000e+00,
             0.0000e+00,  0.0000e+00,  2.0000e-05,  0.0000e+00,  0.0000e+00,
             2.5000e-06,  0.0000e+00,  3.4000e-03,  0.0000e+00,  0.0000e+00,
             0.0000e+00,  0.0000e+00,  1.0000e-03,  0.0000e+00,  0.0000e+00,
             0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,
             0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  6.5000e-05
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
           "TMmIIa"                ,"TMIIaPC"               ,"TMmIIaPC"              ,"TMIIaAPC"              
          };


const PetscReal SURFACE_INITIAL_CONCENTRATIONS_PATCHES[NUMBER_OF_PATCHES][NUMBER_OF_SPECIES_ON_SURFACE]
        = {  5.0000e-08,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,
             0.0000e+00,  0.0000e+00,  0.0000e+00,  1.0379e-11,  0.0000e+00,
             0.0000e+00,  1.0379e-11,  0.0000e+00,  1.1209e-10,  0.0000e+00,
             0.0000e+00,  1.2454e-10,  0.0000e+00,  0.0000e+00,  1.8682e-11,
             0.0000e+00,  0.0000e+00,  8.3029e-11,  0.0000e+00,  0.0000e+00,
             0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,
             0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,
             0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,
             0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,
             0.0000e+00,  0.0000e+00,  0.0000e+00,  1.0000e-09,  0.0000e+00,
             0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00
           };


#define NUMBER_OF_REACTIONS_IN_VOLUME    30
#define NUMBER_OF_REACTIONS_ON_SURFACE   137

#endif
