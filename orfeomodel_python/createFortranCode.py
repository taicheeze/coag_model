# returns list of surf and bulk species
def createListOfSpecies(filename):
    file = open(filename, 'r')
    s = file.read()
    file.close()

    LLines = s.split('\n')

    LSurf = LLines[2:56]
    LBulk = LLines[59:94]

    LSurf2 = [line.split()[1] for line in LSurf]
    LBulk2 = [line.split()[1] for line in LBulk]    
    

    return (LSurf2, LBulk2)

# returns list of initial surf and bulk concentrations
def createListOfInitConc(filename, LSurfSpecies, LBulkSpecies):

    file = open(filename, 'r')
    s = file.read()
    file.close()

    LLines = s.split('\n')

    for line in LLines[:2]:
        if '[' in line:
            line = line[:line.index('[')]
        exec(line)

    LSurfInit = [0.0]*len(LSurfSpecies)
    LBulkInit = [0.0]*len(LBulkSpecies)
    for line in LLines[4:24]:
        if '[' in line:
            line = line[:line.index('[')]
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
def readModelFile():
    file = open('model_dec4_v2.txt', 'r')
    s = file.read()
    file.close()

    LLines = s.split('\n')
    DDefs = {}
    for line in (LLines[1:8] + LLines[10:132]):
        if '[' in line:
            line = line[:line.index('[')]

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
    for line in LLines[135:159]:
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
    for line in LLines[161:253]:
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


# returns a list of items the same length as the no. of
# reactions, with 's', 'b' indicating surface and bulk
# location of reactions
def createListOfReactionLocations():
    file = open('BulkRxnRates.txt', 'r');
    s = file.read()
    file.close()

    LLines = s.split('\n')
    while '' in LLines:  LLines.remove('')

    D = {}
    for line in LLines:
        i = line.split()[0]
        k = int(i[1:])
        v = 'b'  # bulk
        D[k] = v

    file = open('SurfRxnRates.txt', 'r');
    s = file.read()
    file.close()

    LLines = s.split('\n')
    while '' in LLines:  LLines.remove('')

    for line in LLines:
        i = line.split()[0]
        k = int(i[2:])
        v = 's'  # surface
        D[k] = v

    L = list(range(len(D)))
    for k, v in D.items():
        L[k-1] = v

    return L

# return a list of all the rate constants
def createListOfRateConstants():
    file = open('KList.txt', 'r')
    s = file.read()
    file.close()

    LLines = s.split('\n')
    while '' in LLines:  LLines.remove('')

    L = []
    for line in LLines:
        rate = float(line.split()[1])
        L.append(rate)
        
    return L    
    

# output is a list of tuples with the following structure:
# (LHS_species, RHS_species), where LHS_species and
# RHS_species are lists of species
def createListOfAllReactions():
    file = open('networkOutput.txt', 'r')
    s = file.read()
    file.close()

    LLines = s.split('\n')
    while '' in LLines:  LLines.remove('')

    
    LAllRxns = []
    for line in LLines:
        line = line.split(',')[1]
        LHS=line.split(' = ')[0]
        RHS=line.split(' = ')[1]
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
            

        LAllRxns.append( (LHS_sp_ord, LHS_ord, LHS_sp, RHS_sp_ord, RHS_ord, RHS_sp) )
        
        
    return LAllRxns


def createChemicalReactionModelTextFile(LAllRxns, LRxnLoc, KList, LBulkSpecies, LSurfSpecies):

    # first, create the documentation of the reactions at the top of the file
    s = ''
    s += ' '*4 + '!'*6 + ' BULK REACTIONS ' + '!'*55 + '\n'
    s += ' '*4 + '!\n'

    count = 0
    for i in range(len(LAllRxns)):
        LHS_sp_ord, LHS_ord, LHS_sp, RHS_sp_ord, RHS_ord, RHS_sp = LAllRxns[i]
        if LRxnLoc[i] == 'b':
            count += 1
            s += '{}{:5}'.format(' '*4 + '! ', str(count) + '.')

            # form the reactions description
            t = ' + '.join(LHS_sp_ord)
            t += ' -> '
            t += ' + '.join(RHS_sp_ord)
            s += '{:80}'.format(t)

            LHS_idx = [str(LBulkSpecies.index(sp) + 1) for sp in LHS_sp]
            RHS_idx = [str(LBulkSpecies.index(sp) + 1) for sp in RHS_sp]


            t = '('
            for j in range(len(LHS_ord)):
                if LHS_ord[j] != '1':
                    t += '{}*[{}]'.format(LHS_ord[j], LHS_idx[j])
                else:
                    t += '[{}]'.format(LHS_idx[j])

                if j != len(LHS_ord)-1:
                    t += ' + '
                    
            t += ' -> '

            for j in range(len(RHS_ord)):
                if RHS_ord[j] != '1':
                    t += '{}*[{}]'.format(RHS_ord[j], RHS_idx[j])
                else:
                    t += '[{}]'.format(RHS_idx[j])

                if j != len(RHS_ord)-1:
                    t += ' + '
            
            t += ')'
            s += t + '\n'

    s += '\n'*2
    
    # now the Fortran Reaction code
    # bulk reactions
    count = 0
    for i in range(len(LAllRxns)):
        LHS_sp_ord, LHS_ord, LHS_sp, RHS_sp_ord, RHS_ord, RHS_sp = LAllRxns[i]
        if LRxnLoc[i] == 'b':
            count += 1
            s += '{}{}{}  '.format(' '*4 + '! ', 'reaction ', count)

            LHS_idx = [str(LBulkSpecies.index(sp) + 1) for sp in LHS_sp]
            RHS_idx = [str(LBulkSpecies.index(sp) + 1) for sp in RHS_sp]
            
            t = '('
            for j in range(len(LHS_ord)):
                if LHS_ord[j] != '1':
                    t += '{}*[{}]'.format(LHS_ord[j], LHS_idx[j])
                else:
                    t += '[{}]'.format(LHS_idx[j])

                if j != len(LHS_ord)-1:
                    t += ' + '
                    
            t += ' -> '

            for j in range(len(RHS_ord)):
                if RHS_ord[j] != '1':
                    t += '{}*[{}]'.format(RHS_ord[j], RHS_idx[j])
                else:
                    t += '[{}]'.format(RHS_idx[j])

                if j != len(RHS_ord)-1:
                    t += ' + '
            
            t += ')'
            s += t + '\n'
            for j in range(len(LHS_ord)):
                s += '    {}({}, {}) = {}.0_DP\n'.format('rxn_coeffs_LHS', count, LHS_idx[j], LHS_ord[j])
            s += '\n'
            for j in range(len(RHS_ord)):
                s += '    {}({}, {}) = {}.0_DP\n'.format('rxn_coeffs_RHS', count, RHS_idx[j], RHS_ord[j])
            s += '\n'
            
            s += '    kF({}) = {}_DP\n'.format(count, KList[i])
            s += '    kB({}) = {}_DP\n'.format(count, 0.0)

            s += '\n'*2




    s += '\n'*2
    
    # surface reactions
    s += ' '*4 + '!'*6 + ' SURFACE REACTIONS ' + '!'*55 + '\n'
    s += ' '*4 + '!\n'

    count = 0
    for i in range(len(LAllRxns)):
        LHS_sp_ord, LHS_ord, LHS_sp, RHS_sp_ord, RHS_ord, RHS_sp = LAllRxns[i]
        if LRxnLoc[i] == 's':
            count += 1
            s += '{}{:5}'.format(' '*4 + '! ', str(count) + '.')

            # form the reactions description
            t = ' + '.join(LHS_sp_ord)
            t += ' -> '
            t += ' + '.join(RHS_sp_ord)
            s += '{:80}'.format(t)

            LHS_idx_b = [];  LHS_idx_s = []
            for sp in LHS_sp:
                if sp in LBulkSpecies:
                    LHS_idx_b.append(str(LBulkSpecies.index(sp) + 1))
                    LHS_idx_s.append('')
                else:
                    LHS_idx_b.append('')
                    LHS_idx_s.append(str(LSurfSpecies.index(sp) + 1))
                    
            RHS_idx_b = [];  RHS_idx_s = []
            for sp in RHS_sp:
                if sp in LBulkSpecies:
                    RHS_idx_b.append(str(LBulkSpecies.index(sp) + 1))
                    RHS_idx_s.append('')
                else:
                    RHS_idx_b.append('')
                    RHS_idx_s.append(str(LSurfSpecies.index(sp) + 1))

##            print('LHS_idx_b = ', LHS_idx_b)
##            print('LHS_idx_s = ', LHS_idx_s)
##            print('RHS_idx_b = ', RHS_idx_b)
##            print('RHS_idx_s = ', RHS_idx_s)
            
            


            t = '('
            for j in range(len(LHS_ord)):
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

            for j in range(len(RHS_ord)):
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

    s += '\n'*2

    s += \
'''
    ! initialize all coefficients to 0.0 so that only 
    ! those that are non-zero need specified below
    surface_rxn_bulk_coeffs_LHS(:, :) = 0.0d0
    surface_rxn_bulk_coeffs_RHS(:, :) = 0.0d0 
    surface_rxn_surface_coeffs_LHS(:, :) = 0.0d0 
    surface_rxn_surface_coeffs_RHS(:, :) = 0.0d0

    
'''
    
    # now the Fortran Reaction code
    # surface reactions
    count = 0
    for i in range(len(LAllRxns)):
        LHS_sp_ord, LHS_ord, LHS_sp, RHS_sp_ord, RHS_ord, RHS_sp = LAllRxns[i]
        if LRxnLoc[i] == 's':
            count += 1
            s += '{}{}{}  '.format(' '*4 + '! ', 'reaction ', count)


            LHS_idx_b = [];  LHS_idx_s = []
            for sp in LHS_sp:
                if sp in LBulkSpecies:
                    LHS_idx_b.append(str(LBulkSpecies.index(sp) + 1))
                    LHS_idx_s.append('')
                else:
                    LHS_idx_b.append('')
                    LHS_idx_s.append(str(LSurfSpecies.index(sp) + 1))
                    
            RHS_idx_b = [];  RHS_idx_s = []
            for sp in RHS_sp:
                if sp in LBulkSpecies:
                    RHS_idx_b.append(str(LBulkSpecies.index(sp) + 1))
                    RHS_idx_s.append('')
                else:
                    RHS_idx_b.append('')
                    RHS_idx_s.append(str(LSurfSpecies.index(sp) + 1))

##            print('LHS_idx_b = ', LHS_idx_b)
##            print('LHS_idx_s = ', LHS_idx_s)
##            print('RHS_idx_b = ', RHS_idx_b)
##            print('RHS_idx_s = ', RHS_idx_s)
            

                    
            
            t = '('
            for j in range(len(LHS_ord)):
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

            for j in range(len(RHS_ord)):
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
            for j in range(len(LHS_ord)):
                if LHS_idx_b[j] != '':
                    s += '    {}({}, {}) = {}.0_DP\n'.format('surface_rxn_bulk_coeffs_LHS', count, LHS_idx_b[j], LHS_ord[j])
                else:
                    s += '    {}({}, {}) = {}.0_DP\n'.format('surface_rxn_surface_coeffs_LHS', count, LHS_idx_s[j], LHS_ord[j])
            s += '\n'
            
            for j in range(len(RHS_ord)):
                if RHS_idx_b[j] != '':
                    s += '    {}({}, {}) = {}.0_DP\n'.format('surface_rxn_bulk_coeffs_RHS', count, RHS_idx_b[j], RHS_ord[j])
                else:
                    s += '    {}({}, {}) = {}.0_DP\n'.format('surface_rxn_surface_coeffs_RHS', count, RHS_idx_s[j], RHS_ord[j])
            s += '\n'
            
            s += '    k_surf({}) = {}_DP\n'.format(count, KList[i])

            s += '\n'*2
            


    file = open('setup_chemical_reaction_model_orfeomodel.txt', 'w')
    file.write(s)
    file.close()
        
def createGlobalDimensions3DModule(NBulkSpecies, NBulkRxns, NSurfSpecies, NSurfRxns):
    s = \
'''module global_dimensions_3D

    
    ! Spatial domain definitions
    integer, parameter :: NX = 100  ! number of spatial cells in x-direction
    integer, parameter :: NZ = 1  ! number of spatial cells in z-direction (must be >=2 - set NZ=2 for 2D problem)
    integer, parameter :: NY = 50  ! number of spatial cells in y-direction
    integer, parameter :: J1 = 30  ! number of cells in y-direction over uniform region
    integer, parameter :: J2 = 49  ! non-uniform region (J1+1 - J2)

    integer, parameter :: MX = 10  ! number of subcells in x-direction for seeding interface
    integer, parameter :: MZ = 10  ! number of subcells in z-direction for seeding interface
    integer, parameter :: MY = 10  ! number of subcells in y-direction for seeding interface

'''
    t = '    integer, parameter :: N_species = {}'.format(NBulkSpecies)
    s += '{:52}! number of bulk chemical species\n'.format(t)

    t = '    integer, parameter :: N_rxns = {}'.format(NBulkRxns)
    s += '{:52}! no. of chemical reaction pathways (in bulk solution)\n'.format(t)

    t = '    integer, parameter :: N_species_surface = {}'.format(NSurfSpecies)
    s += '{:52}! number of surface chemical species\n'.format(t)

    t = '    integer, parameter :: N_surface_rxns = {}'.format(NSurfRxns)
    s += '{:52}! no. of chemical reaction pathways (on active surface involving both adsorbed and bulk species)\n'.format(t)

    s += \
'''    
    

end module global_dimensions_3D
'''
    file = open('global_dimensions_3D_module_orfeomodel.txt', 'w')
    file.write(s)
    file.close()

def createInitConcFile(LSurfSpecies, LBulkSpecies, LSurfInit, LBulkInit):
    s = \
'''    ! define the number of chemical species (in solution) and their physical properties
                                                                   
    ! BULK SPECIES 
    !
'''
    for i in range(0, len(LBulkInit), 5):
        iMax = min(i+5, len(LBulkInit))
        t = '    !  ' + '{}-{}.'.format(i+1, iMax)
        s += '{:15}'.format(t)
        for j in range(i, iMax):
            s += '{:20}'.format(LBulkSpecies[j] + ',')
        s += '\n'

    s += '\n'*2
    s += '''    ! SURFACE SPECIES 
    !
'''
    for i in range(0, len(LSurfInit), 5):
        iMax = min(i+5, len(LSurfInit))
        t = '    !  ' + '{}-{}.'.format(i+1, iMax)
        s += '{:15}'.format(t)
        for j in range(i, iMax):
            s += '{:20}'.format(LSurfSpecies[j] + ',')
        s += '\n'

    
        
    s += '\n'*3
    t = '    double precision, dimension(N_species), parameter :: D_species = (/ '
    s += t
    for i in range(0, len(LBulkInit), 5):
        iMax = min(i+5, len(LBulkInit))
        if i != 0:
            s += ' '*len(t)
        for j in range(i, iMax):
            if j < len(LBulkInit)-1:
                s += '{:8}'.format('1.0_DP,')
            else:
                s += '{:8}'.format('1.0_DP')
        s += ' &\n'
    s += ' '*len(t) + '/)*1.0e-9_DP/(L_scale**2) ! Diffusivity'


    s += '\n'*3
    s += '    ! not relevant for blood model\n'
    t = '    double precision, dimension(N_species), parameter :: z_species = (/ '
    s += t
    for i in range(0, len(LBulkInit), 5):
        iMax = min(i+5, len(LBulkInit))
        if i != 0:
            s += ' '*len(t)
        for j in range(i, iMax):
            if j < len(LBulkInit)-1:
                s += '{:8}'.format('1.0_DP,')
            else:
                s += '{:8}'.format('1.0_DP')
        s += ' &\n'
    s += ' '*len(t) + '/)  ! Charge'
    
    
    s += '\n'*3
    s += '    ! not relevant for blood model\n'
    t = '    double precision, dimension(N_species), parameter :: u_species = (/ '
    s += t
    for i in range(0, len(LBulkInit), 5):
        iMax = min(i+5, len(LBulkInit))
        if i != 0:
            s += ' '*len(t)
        for j in range(i, iMax):
            if j < len(LBulkInit)-1:
                s += '{:8}'.format('1.0_DP,')
            else:
                s += '{:8}'.format('1.0_DP')
        s += ' &\n'
    s += ' '*len(t) + '/)*(1.0e-9_DP/(8.31_DP*298.0_DP))/(L_scale**2)  ! Mobility'


    s += '\n'*3
    s += '    ! chemical species concentrations at Dirichlet (left) boundary\n'
    t = '    double precision, dimension(N_species), parameter :: c_top = (/ '
    s += t
    for i in range(0, len(LBulkInit), 5):
        iMax = min(i+5, len(LBulkInit))
        if i != 0:
            s += ' '*len(t)
        for j in range(i, iMax):
            if j < len(LBulkInit)-1:
                s += '{:10.2e}_DP,'.format(LBulkInit[j])
            else:
                s += '{:10.2e}_DP'.format(LBulkInit[j])
        s += ' &\n'
    s += ' '*len(t) + '/)'


    s += '\n'*3
    s += '    double precision, parameter :: potential_top = 0.0_DP ! Electric potential at top boundary\n'
    t = '    double precision, dimension(N_species_surface), parameter :: c0_surface = (/'
    s += t
    for i in range(0, len(LSurfInit), 5):
        iMax = min(i+5, len(LSurfInit))
        if i != 0:
            s += ' '*len(t)
        for j in range(i, iMax):
            if j < len(LSurfInit)-1:
                s += '{:10.2e}_DP,'.format(LSurfInit[j])
            else:
                s += '{:10.2e}_DP'.format(LSurfInit[j])
        s += ' &\n'
    s += ' '*len(t) + '/)'
    
    

    file = open('global_data_3D_module_orfeomodel.txt', 'w')
    file.write(s)
    file.close()


def createSurfaceDensityFile(N_species_surface):
    s = '    ! define the surface species densities (mol/L_scale^2)\n'
    t = '    double precision, dimension(N_species_surface), parameter :: surface_density = (/ '
    s += t
    for i in range(0, N_species_surface, 5):
        iMax = min(i+5, N_species_surface)
        if i != 0:
            s += ' '*len(t)
        for j in range(i, iMax):
            if j < N_species_surface-1:
                s += '{:8}'.format('1.0_DP,')
            else:
                s += '{:8}'.format('1.0_DP')
        s += ' &\n'
    s += ' '*len(t) + '/)'

    file = open('surface_rxns_module_orfeomodel.txt', 'w')
    file.write(s)
    file.close()
    
    
def main():
    LAllRxns, LRxnLoc, KList = readModelFile()

    LSurfSpecies, LBulkSpecies = createListOfSpecies('species.txt')
    LSurfInit, LBulkInit = createListOfInitConc('initial_values.txt', LSurfSpecies, LBulkSpecies)


    createChemicalReactionModelTextFile(LAllRxns, LRxnLoc, KList, LBulkSpecies, LSurfSpecies)
    createGlobalDimensions3DModule(len(LBulkSpecies), LRxnLoc.count('b'), len(LSurfSpecies), LRxnLoc.count('s'))
    createInitConcFile(LSurfSpecies, LBulkSpecies, LSurfInit, LBulkInit)
    createSurfaceDensityFile(len(LSurfSpecies))

main()

