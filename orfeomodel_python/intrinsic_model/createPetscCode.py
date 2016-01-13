# returns list of surf and bulk species
def createListOfSpecies(filename):
    file = open(filename, 'r')
    s = file.read()
    file.close()

    LLines = s.split('\n')
    
    for i in range(10):
        for line in LLines:
            if "#" not in line and "." not in line:
                LLines.remove(line)

    for line in LLines:
        if "#surface species" in line:
            surface_species_begin = LLines.index(line)
        if "#bulk species" in line:
            bulk_species_begin = LLines.index(line)

    end_of_species = len(LLines)

    if surface_species_begin<bulk_species_begin:
        LSurf = LLines[(surface_species_begin+1):bulk_species_begin]
        LBulk = LLines[(bulk_species_begin+1):end_of_species]
    else:
        LSurf = LLines[(surface_species_begin+1):end_of_species]
        LBulk = LLines[(bulk_species_begin+1):surface_species_begin]


    LSurf2 = [line.split()[1] for line in LSurf]
    LBulk2 = [line.split()[1] for line in LBulk]    

    return (LSurf2, LBulk2)

# returns list of initial surf and bulk concentrations
def createListOfInitConc(filename, LSurfSpecies, LBulkSpecies):

    file = open(filename, 'r')
    s = file.read()
    file.close()

    LLines = s.split('\n')

    for i in range(10):
        for line in LLines:
            if "#" not in line and "=" not in line:
                LLines.remove(line)

    for line in LLines:
        if "#initial values" in line:
            initial_values_begin = LLines.index(line)

    initial_values_end = len(LLines)

    for line in LLines[:initial_values_begin]:
        if "mol/l" in line:
            line = line + "*1.0e3"
        if '[' in line:
            line = line[:line.index('[')] + line[line.index(']')+1:]
        exec(line)

    LSurfInit = [0.0]*len(LSurfSpecies)
    LBulkInit = [0.0]*len(LBulkSpecies)
    for line in LLines[(initial_values_begin+1):initial_values_end]:
        if "mol/l" in line:
            line = line + "*1.0e3"
        if '[' in line:
            line = line[:line.index('[')] + line[line.index(']')+1:]
        K = line.split(' = ')
        spName = K[0][:-5]  # slice off '_init'
        spConc = eval(K[1])

        if spName in LSurfSpecies:
            idx = LSurfSpecies.index(spName)
            LSurfInit[idx] = spConc
        else:
            idx = LBulkSpecies.index(spName)
            LBulkInit[idx] = spConc

    return (LSurfInit, LBulkInit)


# reads entire model file and returns
# lists of reactions, locations, and rates
def readModelFile(filename):
    file = open(filename, 'r')
    s = file.read()
    file.close()

    LLines = s.split('\n')

    for i in range(10):
        for line in LLines:
            if "#" not in line and "=" not in line:
                LLines.remove(line)
 
    initial_values_begin = -1

    for line in LLines:
        if "#Constants" in line:
            constants_begin = LLines.index(line)
        if "#initial values" in line:
            initial_values_begin = LLines.index(line)
        if "#Rate Constants" in line:
            rate_constants_begin = LLines.index(line)
        if "#Bulk Reactions" in line:
            bulk_reactions_begin = LLines.index(line)
        if "#Surface Reactions" in line:
            surface_reactions_begin = LLines.index(line)

    end_of_reactions = len(LLines)

    if initial_values_begin == -1:
        constants_end = rate_constants_begin
    else:
        constants_end = initial_values_begin

    DDefs = {}
    for line in (LLines[(constants_begin+1):constants_end] + LLines[(rate_constants_begin+1):bulk_reactions_begin]):
        if "l/(mol*s)" in line:
            line = line + "*1.0e-3"
        if '[' in line:
            line = line[:line.index('[')] + line[line.index(']')+1:]

        #print(line)
        exec(line)  # define variables
        K = line.split(' = ')
        DDefs[K[0]] = eval(K[0])
        
    # DDefs is a dictionary of names and float values      

    # bulk reactions
    LAllRxns = []
    LRxnLoc = []
    KList = []
    i = 0
    for line in LLines[(bulk_reactions_begin+1):surface_reactions_begin]:
        K = line.split()
        rxn = K[1]
        i += 1
        #print(i)
        if '<=>' in rxn:  # two reactions
            LHS = rxn.split('<=>')[0]
            RHS = rxn.split('<=>')[1]

            kf = eval(K[4])
            kb = eval(K[7])

        else:  # one reaction
            LHS = rxn.split('=>')[0]
            RHS = rxn.split('=>')[1]

            kf = eval(K[4])

            
        LHS_sp_ord = LHS.split('+')  # includes some '2*sp' type strings
        RHS_sp_ord = RHS.split('+')

        # split up reacting species from their order (number)
        LHS_sp = [];  LHS_ord = []
        for sp_ord in LHS_sp_ord:
            if '*' in sp_ord:
                LHS_ord.append(sp_ord.split('*')[0])
                LHS_sp.append(sp_ord.split('*')[1])
            else:
                LHS_ord.append('1')
                LHS_sp.append(sp_ord)

        RHS_sp = [];  RHS_ord = []
        for sp_ord in RHS_sp_ord:
            if '*' in sp_ord:
                RHS_ord.append(sp_ord.split('*')[0])
                RHS_sp.append(sp_ord.split('*')[1])
            else:
                RHS_ord.append('1')
                RHS_sp.append(sp_ord)
            

        if '<=>' in rxn:  # two reactions
            LAllRxns.append( (LHS_sp_ord, LHS_ord, LHS_sp, RHS_sp_ord, RHS_ord, RHS_sp) )
            LAllRxns.append( (RHS_sp_ord, RHS_ord, RHS_sp, LHS_sp_ord, LHS_ord, LHS_sp) )
            KList.append(kf)
            KList.append(kb)
            LRxnLoc.append('b')
            LRxnLoc.append('b')
        else:
            LAllRxns.append( (LHS_sp_ord, LHS_ord, LHS_sp, RHS_sp_ord, RHS_ord, RHS_sp) )
            KList.append(kf)
            LRxnLoc.append('b')
            
        

    # surface reactions
    i = 0
    for line in LLines[(surface_reactions_begin+1):end_of_reactions]:
        K = line.split()
        rxn = K[1]
        i += 1
        #print(i)
        if '<=>' in rxn:  # two reactions
            LHS = rxn.split('<=>')[0]
            RHS = rxn.split('<=>')[1]

            kf = eval(K[4])
            kb = eval(K[7])

        else:  # one reaction
            LHS = rxn.split('=>')[0]
            RHS = rxn.split('=>')[1]

            kf = eval(K[4])

            
        LHS_sp_ord = LHS.split('+')  # includes some '2*sp' type strings
        RHS_sp_ord = RHS.split('+')

        # split up reacting species from their order (number)
        LHS_sp = [];  LHS_ord = []
        for sp_ord in LHS_sp_ord:
            if '*' in sp_ord:
                LHS_ord.append(sp_ord.split('*')[0])
                LHS_sp.append(sp_ord.split('*')[1])
            else:
                LHS_ord.append('1')
                LHS_sp.append(sp_ord)

        RHS_sp = [];  RHS_ord = []
        for sp_ord in RHS_sp_ord:
            if '*' in sp_ord:
                RHS_ord.append(sp_ord.split('*')[0])
                RHS_sp.append(sp_ord.split('*')[1])
            else:
                RHS_ord.append('1')
                RHS_sp.append(sp_ord)
            

        if '<=>' in rxn:  # two reactions
            LAllRxns.append( (LHS_sp_ord, LHS_ord, LHS_sp, RHS_sp_ord, RHS_ord, RHS_sp) )
            LAllRxns.append( (RHS_sp_ord, RHS_ord, RHS_sp, LHS_sp_ord, LHS_ord, LHS_sp) )
            KList.append(kf)
            KList.append(kb)
            LRxnLoc.append('s')
            LRxnLoc.append('s')
        else:
            LAllRxns.append( (LHS_sp_ord, LHS_ord, LHS_sp, RHS_sp_ord, RHS_ord, RHS_sp) )
            KList.append(kf)
            LRxnLoc.append('s')
            
        


    return (LAllRxns, LRxnLoc, KList)


def createBulkChemicalReactionModelTextFile(LAllRxns, LRxnLoc, KList, LBulkSpecies, LSurfSpecies):

    NBulkSpecies = len(LBulkSpecies)
    NBulkRxns    = LRxnLoc.count('b')

    # first, create the documentation of the reactions at the top of the file
    s = \
'''
#ifndef _VOLUME_REACTIONS_H_
#define _VOLUME_REACTIONS_H_

#include "DataStructure.h"

/* VOLUME REACTIONS:
 *
'''

    # create stoichiometry matrix
    stoich = []
    for i in range(NBulkSpecies):
        stoich_singleline = [0]*NBulkRxns
        stoich.append(stoich_singleline)

    reactants_coef = []
    for i in range(NBulkRxns):
        reactants_coef_singleline = [0]*NBulkSpecies
        reactants_coef.append(reactants_coef_singleline)

    reactants = []

    count = 0
    for i in range(len(LAllRxns)):
        LHS_sp_ord, LHS_ord, LHS_sp, RHS_sp_ord, RHS_ord, RHS_sp = LAllRxns[i]
        if LRxnLoc[i] == 'b':
            s += '{}{:5}'.format(' *  ', str(count) + '.')

            # form the reactions description
            t = ' + '.join(LHS_sp_ord)
            t += ' -> '
            t += ' + '.join(RHS_sp_ord)
            s += '{:80}'.format(t)
            
            if LHS_sp[0] != 'NULL':
                LHS_idx = [LBulkSpecies.index(sp) for sp in LHS_sp]
                reactants.append(LHS_idx)
            else:
                LHS_idx = []
                    
            if RHS_sp[0] != 'NULL':
                RHS_idx = [LBulkSpecies.index(sp) for sp in RHS_sp]
            else:
                RHS_idx = []


            t = '('
            for j in range(len(LHS_idx)):
                stoich[LHS_idx[j]][count] -= int(LHS_ord[j])
                reactants_coef[count][LHS_idx[j]] += int(LHS_ord[j])
                if LHS_ord[j] != '1':
                    t += '{}*[{}]'.format(LHS_ord[j], LHS_idx[j])
                else:
                    t += '[{}]'.format(LHS_idx[j])

                if j != len(LHS_ord)-1:
                    t += ' + '
                    
            t += ' -> '

            for j in range(len(RHS_idx)):
                stoich[RHS_idx[j]][count] += int(RHS_ord[j])
                if RHS_ord[j] != '1':
                    t += '{}*[{}]'.format(RHS_ord[j], RHS_idx[j])
                else:
                    t += '[{}]'.format(RHS_idx[j])

                if j != len(RHS_ord)-1:
                    t += ' + '
            
            t += ')'
            s += t + '\n'
            count += 1

    s += ' */\n\n'

    s += \
'''
/*
 * volume_stoichiometry: stoichiometry matrix of volume species for volume reactions
 */
/*
 * Begin of user-generated volume reaction stoichiometries
 */
PetscInt  volume_stoichiometry[NUMBER_OF_SPECIES_IN_VOLUME][NUMBER_OF_REACTIONS_IN_VOLUME]
          ={
'''

    # now create the stoichiometry matrix
    for i in range(NBulkSpecies):
        s += ' '*12+'{';
        for j in range(NBulkRxns):
            t = '{:2}'.format(stoich[i][j])
            if j != NBulkRxns-1:
               t += ','
            s += t
        if i != NBulkSpecies-1:
            s += '},\n'
        else:
            s += '}\n'

    s += '           };\n\n'

    s += \
'''
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
'''
    
    # now the reaction rate constants
    reaction_rate_constants_string = ''
    count = 0
    for i in range(len(LAllRxns)):
        LHS_sp_ord, LHS_ord, LHS_sp, RHS_sp_ord, RHS_ord, RHS_sp = LAllRxns[i]
        if LRxnLoc[i] == 'b':
            reaction_rate_constants_string += '{}{}{}{}'.format(' '*4, 'r_constant[', count,'] = ')
            reaction_rate_constants_string += '{:10.4e};\n'.format(KList[i])
            count += 1

    s += reaction_rate_constants_string + '\n'

    # reaction rate equations
    for i in range(NBulkRxns):
        s += '    reaction_rate[' + str(i) + '] = r_constant[' + str(i) + ']'
        for j in range(len(reactants_coef[i])):
            if(reactants_coef[i][j]!=0):
                 if(reactants_coef[i][j]==1):
                      s += ' * x_ptr[' + str(j) + ']'
                 else:
                      s += ' * pow(x_ptr[' + str(j) + '],' + reactants_coef[i][j] + ')'
        s += ';\n'

    s += \
'''
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
'''

    # reaction rate constants again
    s += reaction_rate_constants_string + '\n'

    # jac_flag to indicate if J[i][j]==0
    jac_flag = []
    for i in range(NBulkSpecies):
        jac_flag_singleline = [0]*NBulkSpecies
        jac_flag.append(jac_flag_singleline)

    for i in range(NBulkSpecies):
        for k in range(NBulkRxns):
            if stoich[i][k]!=0:
                for j in range(len(reactants_coef[k])):
                    if reactants_coef[k][j]!=0 :
                         jac_flag[i][j] = 1

    # reaction rate jacobians
    for i in range(NBulkSpecies):
        for j in range(NBulkSpecies):
            if jac_flag[i][j] != 0:
                s += '    J[' + str(i) + '][' + str(j) + '] = '
                output_jacobian = ''
                for k in range(NBulkRxns):
                    if stoich[i][k]!=0 and reactants_coef[k][j]!=0:
                        if len(reactants[k])>2:
                            print("Volume reaction more than 2 reactants, need to rewrite python code\n")
                        elif len(reactants[k])==2:
                            if reactants_coef[k][reactants[k][0]] + reactants_coef[k][reactants[k][1]] > 2:
                                print("Volume reaction more than second order, need to rewrite python code\n")
                        if len(reactants[k])==1:
                            if reactants_coef[k][j]==1:
                               if(stoich[i][k]>0):
                                    output_jacobian += ' +' + str(stoich[i][k]) + '*r_constant['  + str(k) + ']'
                               else:
                                    output_jacobian += ' ' + str(stoich[i][k]) + '*r_constant['  + str(k) + ']'
                            else:
                               if(stoich[i][k]>0):
                                   output_jacobian += ' +' + str(stoich[i][k]) + '*r_constant['  + str(k) + ']*' + str(reactants_coef[k][j]) + '*x_ptr[' + str(j) + ']'
                               else:
                                   output_jacobian += ' ' + str(stoich[i][k]) + '*r_constant['  + str(k) + ']*' + str(reactants_coef[k][j]) + '*x_ptr[' + str(j) + ']'
                        elif reactants[k][0]==j:
                            if(stoich[i][k]>0):
                                output_jacobian += ' +' + str(stoich[i][k]) + '*r_constant['  + str(k) + ']*x_ptr[' + str(reactants[k][1]) + ']'
                            else:
                                output_jacobian += ' ' + str(stoich[i][k]) + '*r_constant['  + str(k) + ']*x_ptr[' + str(reactants[k][1]) + ']'
                        elif reactants[k][1]==j:
                            if(stoich[i][k]>0):
                                output_jacobian += ' +' + str(stoich[i][k]) + '*r_constant['  + str(k) + ']*x_ptr[' + str(reactants[k][0]) + ']'
                            else:
                                output_jacobian += ' ' + str(stoich[i][k]) + '*r_constant['  + str(k) + ']*x_ptr[' + str(reactants[k][0]) + ']'
                output_jacobian += ';\n'
                s += output_jacobian

    s += \
'''
/*
 * End of user-generated jacobian function
 */

    PetscFunctionReturn(0);
}

#endif
'''

    file = open('VolumeReactions.h', 'w')
    file.write(s)
    file.close()

def createSurfChemicalReactionModelTextFile(LAllRxns, LRxnLoc, KList, LBulkSpecies, LSurfSpecies):

    NSurfSpecies = len(LSurfSpecies)
    NBulkSpecies = len(LBulkSpecies)
    NSurfRxns    = LRxnLoc.count('s')
    
    # first, create the documentation of the reactions at the top of the file
    s = \
'''
#ifndef _SURFACE_REACTIONS_H_
#define _SURFACE_REACTIONS_H_

#include "DataStructure.h"

/* SURFACE REACTIONS:
 *
'''

    # create surface stoichiometry matrix for surface species and bulk species
    stoich_surf = []
    for i in range(NSurfSpecies):
        stoich_singleline = [0]*NSurfRxns
        stoich_surf.append(stoich_singleline)
    
    stoich_bulk = []
    for i in range(NBulkSpecies):
        stoich_singleline = [0]*NSurfRxns
        stoich_bulk.append(stoich_singleline)

    reactants_surf_coef = []
    for i in range(NSurfRxns):
        reactants_coef_singleline = [0]*NSurfSpecies
        reactants_surf_coef.append(reactants_coef_singleline)

    reactants_bulk_coef = []
    for i in range(NSurfRxns):
        reactants_coef_singleline = [0]*NBulkSpecies
        reactants_bulk_coef.append(reactants_coef_singleline)

    reactants_surf = []
    reactants_bulk = []

    count = 0
    for i in range(len(LAllRxns)):
        LHS_sp_ord, LHS_ord, LHS_sp, RHS_sp_ord, RHS_ord, RHS_sp = LAllRxns[i]
        if LRxnLoc[i] == 's':            
            s += '{}{:5}'.format(' *  ', str(count) + '.')

            # form the reactions description
            t = ' + '.join(LHS_sp_ord)
            t += ' -> '
            t += ' + '.join(RHS_sp_ord)
            s += '{:80}'.format(t)

            LHS_idx_b = [];  LHS_idx_s = []
            for sp in LHS_sp:
                if sp != 'NULL':
                    if sp in LBulkSpecies:
                        LHS_idx_b.append(LBulkSpecies.index(sp))
                        LHS_idx_s.append('')
                    else:
                        LHS_idx_b.append('')
                        LHS_idx_s.append(LSurfSpecies.index(sp))
                    
            RHS_idx_b = [];  RHS_idx_s = []
            for sp in RHS_sp:
                if sp != 'NULL':
                    if sp in LBulkSpecies:
                        RHS_idx_b.append(LBulkSpecies.index(sp))
                        RHS_idx_s.append('')
                    else:
                        RHS_idx_b.append('')
                        RHS_idx_s.append(LSurfSpecies.index(sp))

##            print('LHS_idx_b = ', LHS_idx_b)
##            print('LHS_idx_s = ', LHS_idx_s)
##            print('RHS_idx_b = ', RHS_idx_b)
##            print('RHS_idx_s = ', RHS_idx_s)

            LHS_idx_b_stripped = [s for s in LHS_idx_b if s!='']
            LHS_idx_s_stripped = [s for s in LHS_idx_s if s!='']
            reactants_bulk.append(LHS_idx_b_stripped)
            reactants_surf.append(LHS_idx_s_stripped)

            t = '('
            for j in range(len(LHS_idx_b)):
                if LHS_idx_b[j] != '':
                    stoich_bulk[LHS_idx_b[j]][count] -= int(LHS_ord[j])
                    reactants_bulk_coef[count][LHS_idx_b[j]] += int(LHS_ord[j])
                else:
                    stoich_surf[LHS_idx_s[j]][count] -= int(LHS_ord[j])
                    reactants_surf_coef[count][LHS_idx_s[j]] += int(LHS_ord[j])
                    
                if LHS_ord[j] != '1':
                    if LHS_idx_b[j] != '':
                        t += '{}*[{}]'.format(LHS_ord[j], LHS_idx_b[j])
                    else:
                        t += '{}*<{}>'.format(LHS_ord[j], LHS_idx_s[j])
                else:
                    if LHS_idx_b[j] != '':
                        t += '[{}]'.format(LHS_idx_b[j])
                    else:
                        t += '<{}>'.format(LHS_idx_s[j])

                if j != len(LHS_ord)-1:
                    t += ' + '
                    
            t += ' -> '

            for j in range(len(RHS_idx_b)):
                if RHS_idx_b[j] != '':
                    stoich_bulk[RHS_idx_b[j]][count] += int(RHS_ord[j])
                else:
                    stoich_surf[RHS_idx_s[j]][count] += int(RHS_ord[j])
                    
                if RHS_ord[j] != '1':
                    if RHS_idx_b[j] != '':
                        t += '{}*[{}]'.format(RHS_ord[j], RHS_idx_b[j])
                    else:
                        t += '{}*<{}>'.format(RHS_ord[j], RHS_idx_s[j])
                else:
                    if RHS_idx_b[j] != '':
                        t += '[{}]'.format(RHS_idx_b[j])
                    else:
                        t += '<{}>'.format(RHS_idx_s[j])

                if j != len(RHS_ord)-1:
                    t += ' + '
            
            t += ')'
            s += t + '\n'
            count += 1

    s += ' */\n\n'
    
    s += \
'''
/*
 * surface_stoichiometry_volume_species: stoichiometry matrix of volume species in surface reactions
 * surface_stoichiometry_surface_species: stoichiometry matrix of surface species in surface reactions
 */
/*
 * Begin of user-generated surface reaction stoichiometries
 */
PetscInt  surface_stoichiometry_volume_species[NUMBER_OF_SPECIES_IN_VOLUME][NUMBER_OF_REACTIONS_ON_SURFACE]
          ={
'''

    # now create the surface stoichiometry matrix for bulk species
    for i in range(NBulkSpecies):
        s += ' '*12+'{';
        for j in range(NSurfRxns):
            t = '{:2}'.format(stoich_bulk[i][j])
            if j != NSurfRxns-1:
               t += ','
            s += t
        if i != NBulkSpecies-1:
            s += '},\n'
        else:
            s += '}\n'

    s += '           };\n\n'

    s += \
'''
PetscInt  surface_stoichiometry_surface_species[NUMBER_OF_SPECIES_ON_SURFACE][NUMBER_OF_REACTIONS_ON_SURFACE]
          ={
'''

    # now create the surface stoichiometry matrix for surface species
    for i in range(NSurfSpecies):
        s += ' '*12+'{';
        for j in range(NSurfRxns):
            t = '{:2}'.format(stoich_surf[i][j])
            if j != NSurfRxns-1:
               t += ','
            s += t
        if i != NSurfSpecies-1:
            s += '},\n'
        else:
            s += '}\n'

    s += '           };\n\n'

    s += \
'''
/*
 * End of user-generated surface reactions stoichiometries
 */


/*
 * User-generated routine to calculate reaction rates of surface reactions in a single grid
 * x_ptr: pointer to surface species array 
 * cv_ptr: pointer to volume species array
 * reaction_rate: pointer to reaction rates array
 */
PetscErrorCode SurfaceReactionRateCalculation(PetscReal t,PetscReal x_ptr[NUMBER_OF_SPECIES_ON_SURFACE],PetscReal cv_ptr[NUMBER_OF_SPECIES_IN_VOLUME],PetscReal reaction_rate[NUMBER_OF_REACTIONS_ON_SURFACE])
{
    PetscReal  r_constant[NUMBER_OF_REACTIONS_ON_SURFACE];

    PetscFunctionBegin;

/*
 * Begin of auto-generated surface reaction rates
 */
'''
    
    # now the reaction rate constants
    reaction_rate_constants_string = ''
    count = 0
    for i in range(len(LAllRxns)):
        LHS_sp_ord, LHS_ord, LHS_sp, RHS_sp_ord, RHS_ord, RHS_sp = LAllRxns[i]
        if LRxnLoc[i] == 's':
            reaction_rate_constants_string += '{}{}{}{}'.format(' '*4, 'r_constant[', count,'] = ')
            reaction_rate_constants_string += '{:10.4e};\n'.format(KList[i])
            count += 1

    s += reaction_rate_constants_string + '\n'

    # reaction rate equations
    for i in range(NSurfRxns):
        s += '    reaction_rate[' + str(i) + '] = r_constant[' + str(i) + ']'
        for j in range(len(reactants_surf_coef[i])):
            if(reactants_surf_coef[i][j]!=0):
                 if(reactants_surf_coef[i][j]==1):
                      s += ' * x_ptr[' + str(j) + ']'
                 else:
                      s += ' * pow(x_ptr[' + str(j) + '],' + reactants_surf_coef[i][j] + ')'
        for j in range(len(reactants_bulk_coef[i])):
            if(reactants_bulk_coef[i][j]!=0):
                 if(reactants_bulk_coef[i][j]==1):
                      s += ' * cv_ptr[' + str(j) + ']'
                 else:
                      s += ' * pow(cv_ptr[' + str(j) + '],' + reactants_bulk_coef[i][j] + ')'
        s += ';\n'

    s += \
'''
/*
 * End of auto-generated surface reaction rates
 */

    PetscFunctionReturn(0);
}

/*
 * User-generated routine to calculate jacobian matrix of surface reactions for surface species in a single grid
 * x_ptr: pointer to surface species array 
 * cv_ptr: pointer to volume species array
 * J: pointer to jacobian matrix for surface species
 */
PetscErrorCode SurfaceReactionJacobianCalculation(PetscReal t,PetscReal x_ptr[NUMBER_OF_SPECIES_ON_SURFACE],PetscReal cv_ptr[NUMBER_OF_SPECIES_IN_VOLUME], PetscReal J[NUMBER_OF_SPECIES_ON_SURFACE][NUMBER_OF_SPECIES_ON_SURFACE])
{
    PetscReal  r_constant[NUMBER_OF_REACTIONS_ON_SURFACE];

    PetscFunctionBegin;

/*
 * Begin of user-generated jacobian function
 */
'''

    # reaction rate constants again
    s += reaction_rate_constants_string + '\n'

    # jac_flag to indicate if J[i][j]==0
    jac_flag = []
    for i in range(NSurfSpecies):
        jac_flag_singleline = [0]*NSurfSpecies
        jac_flag.append(jac_flag_singleline)

    for i in range(NSurfSpecies):
        for k in range(NSurfRxns):
            if stoich_surf[i][k]!=0:
                for j in range(len(reactants_surf_coef[k])):
                    if reactants_surf_coef[k][j]!=0 :
                         jac_flag[i][j] = 1

    # reaction rate jacobians
    for i in range(NSurfSpecies):
        for j in range(NSurfSpecies):
            if jac_flag[i][j] != 0:
                s += '    J[' + str(i) + '][' + str(j) + '] = '
                output_jacobian = ''
                for k in range(NSurfRxns):
                    if stoich_surf[i][k]!=0 and reactants_surf_coef[k][j]!=0:
                        if len(reactants_surf[k])+len(reactants_bulk[k])==0:
                            print("Surface reaction zero-th order, need to rewrite python code\n")
                        elif len(reactants_surf[k])+len(reactants_bulk[k])>2:
                            print("Surface reaction more than 2 reactants, need to rewrite python code\n")
                        elif len(reactants_surf[k])==2:
                            if reactants_surf_coef[k][reactants_surf[k][0]] + reactants_surf_coef[k][reactants_surf[k][1]] > 2:
                                print("Surface reaction more than second order, need to rewrite python code\n")
                                
                        if len(reactants_surf[k])==1:
                            if reactants_surf_coef[k][j]==1:
                               if(stoich_surf[i][k]>0):
                                    output_jacobian += ' +' + str(stoich_surf[i][k]) + '*r_constant['  + str(k) + ']'
                               else:
                                    output_jacobian += ' ' + str(stoich_surf[i][k]) + '*r_constant['  + str(k) + ']'
                            else:
                               if(stoich_surf[i][k]>0):
                                   output_jacobian += ' +' + str(stoich_surf[i][k]) + '*r_constant['  + str(k) + ']*' + str(reactants_surf_coef[k][j]) + '*x_ptr[' + str(j) + ']'
                               else:
                                   output_jacobian += ' ' + str(stoich_surf[i][k]) + '*r_constant['  + str(k) + ']*' + str(reactants_surf_coef[k][j]) + '*x_ptr[' + str(j) + ']'
                        elif reactants_surf[k][0]==j:
                            if(stoich_surf[i][k]>0):
                                output_jacobian += ' +' + str(stoich_surf[i][k]) + '*r_constant['  + str(k) + ']*x_ptr[' + str(reactants_surf[k][1]) + ']'
                            else:
                                output_jacobian += ' ' + str(stoich_surf[i][k]) + '*r_constant['  + str(k) + ']*x_ptr[' + str(reactants_surf[k][1]) + ']'
                        elif reactants_surf[k][1]==j:
                            if(stoich_surf[i][k]>0):
                                output_jacobian += ' +' + str(stoich_surf[i][k]) + '*r_constant['  + str(k) + ']*x_ptr[' + str(reactants_surf[k][0]) + ']'
                            else:
                                output_jacobian += ' ' + str(stoich_surf[i][k]) + '*r_constant['  + str(k) + ']*x_ptr[' + str(reactants_surf[k][0]) + ']'
                                                            
                        for l in range(len(reactants_bulk_coef[k])):
                            if(reactants_bulk_coef[k][l]!=0):
                                if(reactants_bulk_coef[k][l]==1):
                                    output_jacobian += '*cv_ptr[' + str(l) + ']'
                                else:
                                    output_jacobian += '*pow(cv_ptr[' + str(l) + '],' + reactants_bulk_coef[k][l] + ')'

                output_jacobian += ';\n'
                s += output_jacobian

    s += \
'''
/*
 * End of user-generated jacobian function
 */

    PetscFunctionReturn(0);
}

#endif
'''

    file = open('SurfaceReactions.h', 'w')
    file.write(s)
    file.close()

        
def createInitConcFile(LSurfSpecies, LBulkSpecies, LSurfInit, LBulkInit, NBulkSpecies, NBulkRxns, NSurfSpecies, NSurfRxns):
    s = \
'''#ifndef _INITIAL_CONCENTRATIONS_H_
#define _INITIAL_CONCENTRATIONS_H_

#include <petscsys.h>

'''

    t = '#define NUMBER_OF_SPECIES_IN_VOLUME     {}\n'.format(NBulkSpecies)
    s += t

    t = '#define NUMBER_OF_SPECIES_ON_SURFACE    {}\n'.format(NSurfSpecies)
    s += t

    s += '\n'

    s += 'const char VOLUME_SPECIES_NAMES[NUMBER_OF_SPECIES_IN_VOLUME][256]\n'
    t = '        = {'
    s += t

    for i in range(0, len(LBulkInit), 5):
        iMax = min(i+5, len(LBulkInit))
        if i!=0:
            s += ' '*len(t)
        for j in range(i, iMax):
            if j < len(LBulkInit)-1:
               s += '{:24}'.format('\"'+LBulkSpecies[j]+'\"') + ','
            else:
               s += '{:24}'.format('\"'+LBulkSpecies[j]+'\"')
        s += '\n'
    s += ' '*(len(t)-1) + '};\n'

    s += '\n\n'

    s += 'const PetscReal VOLUME_INITIAL_CONCENTRATIONS[NUMBER_OF_SPECIES_IN_VOLUME]\n'

    t = '        = {'
    s += t
    for i in range(0, len(LBulkInit), 5):
        iMax = min(i+5, len(LBulkInit))
        if i != 0:
            s += ' '*len(t)
        for j in range(i, iMax):
            if j < len(LBulkInit)-1:
                s += '{:12.4e},'.format(LBulkInit[j])
            else:
                s += '{:12.4e}'.format(LBulkInit[j])
        s += '\n'
    s += ' '*len(t) + '};\n'

    s += '\n'*2

    s += 'const char SURFACE_SPECIES_NAMES[NUMBER_OF_SPECIES_ON_SURFACE][256]\n'
    t = '        = {'
    s += t

    for i in range(0, len(LSurfInit), 5):
        iMax = min(i+5, len(LSurfInit))
        if i!=0:
            s += ' '*len(t)
        for j in range(i, iMax):
            if j < len(LSurfInit)-1:
                s += '{:24}'.format('\"'+LSurfSpecies[j]+'\"') + ','
            else:
                s += '{:24}'.format('\"'+LSurfSpecies[j]+'\"')
        s += '\n'
    s += ' '*(len(t)-1) + '};\n'

    s += '\n\n'
    
    s += 'const PetscReal SURFACE_INITIAL_CONCENTRATIONS_PATCHES[NUMBER_OF_PATCHES][NUMBER_OF_SPECIES_ON_SURFACE]\n'

    t = '        = {\n'
    s += t
    t = '           {'
    s += t
    for i in range(0, len(LSurfInit), 5):
        iMax = min(i+5, len(LSurfInit))
        if i != 0:
            s += ' '*len(t)
        for j in range(i, iMax):
            if j < len(LSurfInit)-1:
                s += '{:12.4e},'.format(LSurfInit[j])
            else:
                s += '{:12.4e}'.format(LSurfInit[j])
        s += '\n'
    s += ' '*(len(t)-1) + '}\n'
    s += ' '*(len(t)-2) + '};\n'

    s += '\n'*2    

    t = '#define NUMBER_OF_REACTIONS_IN_VOLUME    {}\n'.format(NBulkRxns)
    s += t

    t = '#define NUMBER_OF_REACTIONS_ON_SURFACE   {}\n'.format(NSurfRxns)
    s += t

    s += '\n'

    s += '#endif\n'
    
    file = open('InitialConcentrations.h', 'w')
    file.write(s)
    file.close()
   
    
def main():
    import sys

    if len(sys.argv) == 2:
        LAllRxns, LRxnLoc, KList = readModelFile(sys.argv[1])
    else:
        print("Please use format: python3 createPetscCode.py \'ModelFileName\'")
        return (-1)

    LSurfSpecies, LBulkSpecies = createListOfSpecies('species.txt')
    
    LSurfInit, LBulkInit = createListOfInitConc('initial_values.txt', LSurfSpecies, LBulkSpecies)

    createBulkChemicalReactionModelTextFile(LAllRxns, LRxnLoc, KList, LBulkSpecies, LSurfSpecies)
    createSurfChemicalReactionModelTextFile(LAllRxns, LRxnLoc, KList, LBulkSpecies, LSurfSpecies)
    createInitConcFile(LSurfSpecies, LBulkSpecies, LSurfInit, LBulkInit, len(LBulkSpecies), LRxnLoc.count('b'), len(LSurfSpecies), LRxnLoc.count('s'))


main()

