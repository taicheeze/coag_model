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
        = {  2.2200e-15,  1.0052e-05,  8.9182e-07,  7.3817e-07,  1.2600e-04,
             8.9200e-05,  8.1900e-10,  1.5700e-04,  3.2800e-46,  1.3900e-09,
             7.0600e-07,  7.0600e-07,  5.8100e-46,  1.9600e-06,  1.3800e-06,
             6.9600e-07,  1.8000e-06,  2.1100e-03,  2.2300e-05,  1.3200e-04,
             3.0900e-09,  1.1400e-03,  9.3700e-04,  6.2500e-05,  1.3100e-08,
             1.3500e-09,  4.1200e-06,  9.1500e-09,  3.1600e-09,  7.0100e-12,
             1.1300e-07,  1.3900e-05,  1.3900e-05,  2.5800e-08,  6.4900e-05,
             0.0000e+00,  0.0000e+00,  0.0000e+00
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
              0.0000e+00,  0.0000e+00,  0.0000e+00,  1.0379e-09,  0.0000e+00,
              0.0000e+00,  1.0379e-09,  0.0000e+00,  1.1209e-08,  0.0000e+00,
              0.0000e+00,  1.2454e-08,  0.0000e+00,  0.0000e+00,  1.8682e-09,
              0.0000e+00,  0.0000e+00,  8.3029e-09,  0.0000e+00,  0.0000e+00,
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
