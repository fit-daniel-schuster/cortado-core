from cortado_core.subprocess_discovery.concurrency_trees.cTrees import BARROW, cTreeOperator
from cortado_core.subprocess_discovery.subtree_mining.tree_pattern import TreePattern

def __compute_size_of_tree(tree) -> int: 
    
    size = 1
    
    if tree.op: 
        
        for child in tree.children: 
            size += __compute_size_of_tree(child)
    
    return size


def compute_valid_leaf_eliminated_children(tree):
    
    valid_substrings = []
    full_rep = [(repr(child) + BARROW) for child in tree.children]
    
    if len(tree.children) > 2:
        
        if tree.op != cTreeOperator.Sequential: 
              
            for index, child in enumerate(tree.children): 
                
                if child.label:
            
                    cString = "".join(full_rep[:index] + full_rep[index + 1:])
                    valid_substrings.append(tree.op.value + cString)
        
    for index, child in enumerate(tree.children): 
        
        if child.op:
            
            for sString in compute_valid_leaf_eliminated_children(child): 

                cString = "".join(full_rep[:index] + [sString + BARROW] + full_rep[index + 1:])
                valid_substrings.append(tree.op.value + cString)
            
        
    return valid_substrings

    
    
def _compute_subtree_eliminated_children(tree): 
    
    valid_substrings = []
    tSize = __compute_size_of_tree(tree)
    full_rep = [(repr(child) + BARROW) for child in tree.children]
    
    if len(tree.children) > 2:
        
        if tree.op != cTreeOperator.Sequential: 
            
            for index, child in enumerate(tree.children): 
                
                # All op nodes in the left tree have at least 2 children
                if child.op:
            
                    # Add the string of the child left-out
                    cString = "".join(full_rep[:index] + full_rep[index + 1:])
                    cSize  =  tSize - __compute_size_of_tree(child)
                    
                    valid_substrings.append((cSize, tree.op.value + cString))
          
    for index, child in enumerate(tree.children): 

        cSize  = __compute_size_of_tree(child)

        for (iSize, sString) in _compute_subtree_eliminated_children(child):  
            sSize = tSize - (cSize - iSize)
            cString = "".join(full_rep[:index] + [sString + BARROW] + full_rep[index + 1:])
            valid_substrings.append((sSize, tree.op.value + cString))
            
    return valid_substrings     
    
# Computes valid subtree nested under the Root
def _get_root_enclosed_subtrees(root, full): 
    
    valid_substrings = []
    
    if full or len(root.children) == 2: 
    
        for child in root.children: 
            
            if child.op: 
                
                cSize = __compute_size_of_tree(child)
                valid_substrings.append((cSize, repr(child)))
            
    return valid_substrings



def _compute_left_out_subtree_strings(tree): 
    
    valid_substrings = []
    tSize = __compute_size_of_tree(tree)
    tmpSize = 0
    
    full_rep = [(repr(child) + BARROW) for child in tree.children]
    
    for cIdx, child in enumerate(tree.children): 
        cSize = __compute_size_of_tree(child)
        
        if child.op: 
            
            sStrings = _compute_left_out_subtree_strings(child)
            
            if 0 < cIdx < len(tree.children) - 1:
                cString = ''.join(full_rep[cIdx + 1:])
                cString = str(tree.op.value) + cString
                valid_substrings.append((tSize - tmpSize - cSize, cString))
                
                for sSize, sString in sStrings: 
                    cString = ''.join(full_rep[cIdx + 1:])
                    cString = str(tree.op.value) + sString + BARROW + cString
                    valid_substrings.append((tSize - tmpSize - (cSize - sSize), cString)) 

            for sSize, sString in sStrings: 
                cString = ''.join(full_rep[:cIdx] + [sString + BARROW] + full_rep[cIdx + 1:])
                cString = str(tree.op.value) + cString
                valid_substrings.append((tSize - (cSize - sSize), cString))
            
        tmpSize += cSize
                    
    return valid_substrings 



# Move along the left-most path and eleiminate left-most activity leafs in sequential groups
def _compute_left_most_path_eliminated_leafs(tree): 
    
    valid_substrings = []

    
    if len(tree.children) > 0: 
        
        full_rep = [(repr(child) + BARROW) for child in tree.children]
        lmc = tree.children[0]

        if len(tree.children) > 2:
            
            if tree.op == cTreeOperator.Sequential: 
                
                if lmc.label:   
                    
                    cString = "".join(full_rep[1:])
                    valid_substrings.append(tree.op.value  + cString)
        
        if tree.op == cTreeOperator.Sequential: 
        
            if lmc.op: 
            
                for sString in _compute_left_most_path_eliminated_leafs(lmc): 

                        cString = "".join([sString + BARROW] + full_rep[1:])
                        valid_substrings.append(tree.op.value + cString)
        else: 
            
            for index, child in enumerate(tree.children): 
                    
                if child.op:
                    
                    for sString in _compute_left_most_path_eliminated_leafs(child): 

                        cString = "".join(full_rep[:index] + [sString + BARROW] + full_rep[index + 1:])
                        valid_substrings.append(tree.op.value + cString)
            
    return valid_substrings

# Move along the right-most path and eleminate right-most activity leafs in sequential groups
def _compute_right_most_path_eliminated_leafs(tree): 
    
    valid_substrings = []
    
    if len(tree.children) > 0: 
        
        full_rep = [(repr(child) + BARROW) for child in tree.children]
        rmc = tree.children[-1]
        
        if len(tree.children) > 2:
            
            if tree.op == cTreeOperator.Sequential: 
                
                if rmc.label: 
                    
                    cString = "".join(full_rep[:-1])
                    valid_substrings.append(tree.op.value  + cString)
                
        
        if tree.op == cTreeOperator.Sequential: 
                        
            if rmc.op: 
            
                for sString in _compute_right_most_path_eliminated_leafs(rmc): 

                        cString = "".join(full_rep[:-1] + [sString + BARROW])
                        valid_substrings.append(tree.op.value + cString)
        else: 
            
            for index, child in enumerate(tree.children): 
                    
                if child.op:
                    
                    for sString in _compute_right_most_path_eliminated_leafs(child): 

                        cString = "".join(full_rep[:index] + [sString + BARROW] + full_rep[index + 1:])
                        valid_substrings.append(tree.op.value + cString)
            
    return valid_substrings

# Move along the right-most path and eleminate right-most subtrees
def _compute_right_most_path_eliminated_subtree(tree): 
    
    valid_substrings = []
    
    if len(tree.children) > 0: 
        tSize = __compute_size_of_tree(tree)
        full_rep = [(repr(child) + BARROW) for child in tree.children]
        
        rmc = tree.children[-1]
        rSize = __compute_size_of_tree(rmc)
        if len(tree.children) > 2:
            if tree.op == cTreeOperator.Sequential: 
                if rmc.op: 
                    cString = "".join(full_rep[:-1])
                    cSize  =  tSize - rSize
                    valid_substrings.append((cSize, tree.op.value  + cString))
            
        if tree.op == cTreeOperator.Sequential: 
            if rmc.op: 
            
                for (iSize, sString) in _compute_right_most_path_eliminated_subtree(rmc): 
                    sSize = tSize - (rSize - iSize)
                    cString = "".join(full_rep[:-1] + [sString + BARROW])
                    valid_substrings.append((sSize, tree.op.value + cString))
                    
        else: 
            
            for index, child in enumerate(tree.children): 
                cSize  = __compute_size_of_tree(child)
                for (iSize, sString) in _compute_right_most_path_eliminated_subtree(child): 
                    sSize = tSize - (cSize - iSize)
                    cString = "".join(full_rep[:index] + [sString + BARROW] + full_rep[index + 1:])
                    valid_substrings.append((sSize, tree.op.value + cString))
   
            
    return valid_substrings

# Move along the left-most path and eliminate left-most subtrees
def _compute_left_most_path_eliminated_subtree(tree): 
    
    valid_substrings = []
    
    if len(tree.children) > 0: 
        tSize = __compute_size_of_tree(tree)
        full_rep = [(repr(child) + BARROW) for child in tree.children]
        
        lmc = tree.children[0]
        lSize = __compute_size_of_tree(lmc)
        
        if len(tree.children) > 2:
            if tree.op == cTreeOperator.Sequential: 
                if lmc.op:   
                        
                    cString = "".join(full_rep[1:])
                    cSize  =  tSize - lSize
                    valid_substrings.append((cSize, tree.op.value  + cString))
                

        if tree.op == cTreeOperator.Sequential: 
            if lmc .op: 
                for (iSize, sString) in _compute_left_most_path_eliminated_subtree(lmc): 
                    sSize = tSize - (lSize - iSize)
                    cString = "".join([sString + BARROW] + full_rep[1:])
                    valid_substrings.append((sSize, tree.op.value + cString))
                    
        else: 
            for index, child in enumerate(tree.children): 
                cSize  = __compute_size_of_tree(child)
                for (iSize, sString) in _compute_left_most_path_eliminated_subtree(child): 
                    sSize = tSize - (cSize - iSize)
                    cString = "".join(full_rep[:index] + [sString + BARROW] + full_rep[index + 1:])
                    valid_substrings.append((sSize, tree.op.value + cString))
   
            
    return valid_substrings


if __name__ == '__main__':
    
    import timeit
    from cortado_core.subprocess_discovery.concurrency_trees.cTrees import ConcurrencyTree, cTreeOperator
    from cortado_core.subprocess_discovery.subtree_mining.treebank import create_treebank_from_cv_variants
    from cortado_core.tests.pattern_mining.example_log import create_example_log_1, create_example_log_2
    from cortado_core.utils.cvariants import get_concurrency_variants
    from cortado_core.subprocess_discovery.subtree_mining.right_most_path_extension.min_sub_mining import min_sub_mining
    from pm4py.objects.log.importer.xes.importer import apply as xes_import
    from cortado_core.subprocess_discovery.subtree_mining.treebank import (
        create_treebank_from_cv_variants,
    )
    from cortado_core.subprocess_discovery.subtree_mining.maximal_connected_components.maximal_connected_check import (
        set_maximaly_closed_patterns,
    )
    from cortado_core.subprocess_discovery.concurrency_trees.cTrees import cTreeOperator
    from cortado_core.subprocess_discovery.subtree_mining.freq_counting import (
        FrequencyCountingStrategy,
        _get_prune_sets,
        compute_frequent_activity_sets,
    )
    from cortado_core.subprocess_discovery.subtree_mining.right_most_path_extension.initial_candidate_generation import (
        generate_initial_candidates,
    )
    
    from cortado_core.subprocess_discovery.subtree_mining.tree_pruning import (
        compute_f3_pruned_set,
    )
    
    from cortado_core.utils.timestamp_utils import TimeUnit

    freq_strat = FrequencyCountingStrategy.TraceTransaction

    l = xes_import("C:\\Users\\Michael\\Desktop\\VQL_Evaluation\\Datasets\\BPI_Challenge_2012.xes")
    art_start = False
    variants = get_concurrency_variants(l, False, TimeUnit.MS)
    treebank = create_treebank_from_cv_variants(variants, artifical_start=art_start)

    min_sup = 500
    k = 100 

    start_time = timeit.default_timer()
    cm_k_patterns = min_sub_mining(
        treebank,
        variants,
        frequency_counting_strat=freq_strat,
        k_it=100,
        min_sup=min_sup,
        artifical_start=art_start,
    )
    
    print('CM Time', timeit.default_timer() - start_time)
    print()
    print("Setting CM Closed")
    set_maximaly_closed_patterns(cm_k_patterns)
    
    fSets = compute_frequent_activity_sets(
        variants, freq_strat, min_sup
    )
    F = generate_initial_candidates(
        treebank, min_sup, freq_strat, fSets, False
    )
    
    # Store the results
    k_pattern = {2: F}
    pSets = _get_prune_sets(fSets, F)
    skipPrune = True
    has_fallthroughs = False

    # For every k > 2 create the k pattern from the frequent k-1 pattern
    for k in range(6):

        newF = []

        for tp in F:

            # Compute the right most path extension of all k-1 pattern
            tps = tp.right_most_path_extension(pSets, skipPrune, has_fallthroughs)

            if(str(tp) == '∧(→(A_FINALIZED, O_CREATED)) with support: 726'): 
                print()
                print('Generated Extensions:')
                for p in tps:
                    print(p)
        
            sup_to_gain = tp.support

            for c in tps:
                if f := c.update_rmo_list(
                    treebank, min_sup, freq_strat, sup_to_gain
                ):
                    newF.append(f)

        # For each candidate update the rmo and through this compute the support
        F = newF

        # Break early, if there is no frequent pattern lefts
        if len(F) > 0:
            k_pattern[3 + k] = F
        else:
            break

        if k == 0:
            pSets, F = compute_f3_pruned_set(pSets, F)
            skipPrune = False
    
   
    
    for k in range(2, 6):
        print()
        print('K: ', k)
        for p in cm_k_patterns[k]:
            
            print(p)
                        