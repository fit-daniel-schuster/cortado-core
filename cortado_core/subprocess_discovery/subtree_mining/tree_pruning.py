from collections import defaultdict
from cortado_core.subprocess_discovery.concurrency_trees.cTrees import cTreeOperator
from cortado_core.subprocess_discovery.subtree_mining.freq_counting import PruningSets


def compute_f2_pruned_set(F):
    
    fA = set()
    fDF = set()
    fCC = set()  
    
    for f in F: 
        child = f.tree.children[0]
        
        if f.tree.op == cTreeOperator.Sequential:
            if child.label:
                fDF.add(child.label)
    
        elif f.tree.op == cTreeOperator.Concurrent:
            if child.label:
                fCC.add(child.label)
                
        elif f.tree.op == cTreeOperator.Fallthrough:
            if child.label:
                fA.add(child.label)
                
    return fA, fDF, fCC

def compute_f3_pruned_set(pSets : PruningSets, F):
    
    # TODO Delete the invalid f3 patterns 
    
    fA = defaultdict(set)
    fDF = defaultdict(set)
    fCC = defaultdict(set)
    fOP = defaultdict(set)
    fOPOP = defaultdict(set)
    fOPACT = defaultdict(set)
    
    valid_F = set()
    
    
    for f in F:
        if len(f.tree.children) > 1: 
            
            lChild = f.tree.children[0]
            rChild = f.tree.children[1]
            
            if f.tree.op == cTreeOperator.Sequential:
                if lChild.label:
                    
                    if rChild.label:
                        fDF[lChild.label].add(rChild.label)
                    else: 
                        fOP[lChild.label].add(rChild.op)
                    
                    valid_F.add(f)
                    
            elif f.tree.op == cTreeOperator.Concurrent:
                if lChild.label:
                    
                    if rChild.label:
                        fCC[lChild.label].add(rChild.label)
                    else: 
                        fOP[lChild.label].add(rChild.op)
                        
                    valid_F.add(f)
                
            elif f.tree.op == cTreeOperator.Fallthrough:
                if lChild.label:
                    
                    if rChild.label:
                        fA[lChild.label].add(rChild.label)
                    else: 
                        fOP[lChild.label].add(rChild.op)
                        
                    valid_F.add(f)
                        
            if lChild.op and rChild.op: 
                fOPOP[lChild.op].add(rChild.op)
                
            if lChild.op and rChild.label: 
                fOPACT[lChild.op].add(rChild.label)
                
        else: 
            valid_F.add(f)
                
    pSets.ftLabelPrune = fA
    pSets.dfLabelPrune = fDF
    pSets.ccLabelPrune = fCC
    pSets.operatorPrune = fOP
    pSets.operatorActivityPrune = fOPACT
    pSets.operatorOperatorPrune = fOPOP

    return pSets, valid_F


