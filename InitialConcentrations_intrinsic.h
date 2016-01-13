#ifndef _INITIAL_CONCENTRATIONS_H_
#define _INITIAL_CONCENTRATIONS_H_

#include <petscsys.h>

#define NUMBER_OF_SPECIES_IN_VOLUME     38
#define NUMBER_OF_SPECIES_ON_SURFACE    55

const char VOLUME_SPECIES_NAMES[NUMBER_OF_SPECIES_IN_VOLUME][256]
        = {"VII"                   ,"VIIa"                  ,"Xa"                    ,"IIa"                   ,"X"                     ,
           "IX"                    ,"IXa"                   ,"II"                    ,"VIII"                  ,"VIIIa"                 ,
           "VIIIa1L"               ,"VIIIa2"                ,"V"                     ,"Va"                    ,"mIIa"                  ,
           "TFPI"                  ,"XaTFPI"                ,"ATIII"                 ,"XaATIII"               ,"mIIaATIII"             ,
           "IXaATIII"              ,"IIaATIII"              ,"substrate"             ,"fluorescence"          ,"APC"                   ,
           "APCVa"                 ,"Va5"                   ,"Va3"                   ,"APCVa5"                ,"APCVa3"                ,
           "Va53"                  ,"HCF"                   ,"LCA1"                  ,"APCLCA1"               ,"PC"                    ,
           "XI"                    ,"XIa"                   ,"IIaXI"                 
          };


const PetscReal VOLUME_INITIAL_CONCENTRATIONS[NUMBER_OF_SPECIES_IN_VOLUME]
        = {  0.0000e-05,  1.0000e-05,  3.4700e-06,  4.2900e-06,  1.6200e-06,
             3.2600e-05,  1.6700e-05,  1.2400e-05,  0.0000e-07,  1.2000e-09,
             6.9900e-07,  6.9900e-07,  0.0000e-05,  1.7200e-06,  4.4200e-09,
             2.6300e-07,  1.4000e-06,  1.6600e-03,  1.4000e-04,  4.7400e-06,
             3.7500e-05,  1.5600e-03,  9.0100e-04,  9.9100e-05,  1.4000e-08,
             1.2600e-09,  3.5700e-06,  8.4500e-09,  2.9200e-09,  6.9300e-12,
             1.0400e-07,  1.2900e-05,  1.2900e-05,  2.5500e-08,  6.4900e-05,
             0.0000e-05,  3.0000e-05,  0.0000e+00
           };
;


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
           "TMmIIa"                ,"TMIIaPC"               ,"TMmIIaPC"              ,"TMIIaAPC"              ,"IXP9XIa"               
          };


const PetscReal SURFACE_INITIAL_CONCENTRATIONS_PATCHES[NUMBER_OF_PATCHES][NUMBER_OF_SPECIES_ON_SURFACE]
        = {
           {  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,
              0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,
              0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,
              0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,
              0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,
              0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,
              0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,
              0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,
              0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,
              0.0000e+00,  0.0000e+00,  0.0000e+00,  2.5000e-13,  0.0000e+00,
              0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00
           },
           {  1.2500e-12,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,
              0.0000e+00,  0.0000e+00,  0.0000e+00,  2.0757e-09,  0.0000e+00,
              0.0000e+00,  2.0757e-09,  0.0000e+00,  2.2418e-08,  0.0000e+00,
              0.0000e+00,  2.4909e-08,  0.0000e+00,  0.0000e+00,  3.7363e-09,
              0.0000e+00,  0.0000e+00,  1.6606e-08,  0.0000e+00,  0.0000e+00,
              0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,
              0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,
              0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,
              0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,
              0.0000e+00,  0.0000e+00,  0.0000e+00,  2.5000e-13,  0.0000e+00,
              0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00
           },
           {  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,
              0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,
              0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,
              0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,
              0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,
              0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,
              0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,
              0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,
              0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,
              0.0000e+00,  0.0000e+00,  0.0000e+00,  2.5000e-13,  0.0000e+00,
              0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00
           }
          };


#define NUMBER_OF_REACTIONS_IN_VOLUME    34
#define NUMBER_OF_REACTIONS_ON_SURFACE   140

#endif