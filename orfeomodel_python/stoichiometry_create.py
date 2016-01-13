stoich = []
for i in range(len(LBulkSpecies)):
    stoich_singleline = [0]*LRxnLoc.count('b')
    stoich.append(stoich_singleline)

count = 0
for i in range(len(LAllRxns)):
     LHS_sp_ord, LHS_ord, LHS_sp, RHS_sp_ord, RHS_ord, RHS_sp = LAllRxns[i]
     if LRxnLoc[i] == 'b':
        LHS_idx = [LBulkSpecies.index(sp) for sp in LHS_sp]
        RHS_idx = [LBulkSpecies.index(sp) for sp in RHS_sp]
        for j in range(len(LHS_ord)):
           stoich[LHS_idx[j]][count] -= int(LHS_ord[j])
        for j in range(len(RHS_ord)):
           stoich[RHS_idx[j]][count] += int(RHS_ord[j])
        count += 1
 

