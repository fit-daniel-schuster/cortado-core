from typing import Mapping
from cortado_core.subprocess_discovery.concurrency_trees.cTrees import cTreeOperator
from cortado_core.subprocess_discovery.subtree_mining.ct_frequency_counting import ct_compute_frequent_activity_sets
from cortado_core.subprocess_discovery.subtree_mining.freq_counting import (
    FrequencyCountingStrategy,
    _get_prune_sets,
    _update_sets_artifical_start,
    compute_frequent_activity_sets,
)
from cortado_core.subprocess_discovery.subtree_mining.right_most_path_extension.initial_candidate_generation import (
    generate_initial_candidates,
)
from cortado_core.subprocess_discovery.subtree_mining.tree_pruning import (
    compute_f3_pruned_set,
)
from cortado_core.subprocess_discovery.subtree_mining.treebank import TreeBankEntry
from cortado_core.subprocess_discovery.subtree_mining.utilities import (
    _contains_fallthrough,
)
from cortado_core.subprocess_discovery.subtree_mining.maximal_connected_components.maximal_connected_check import (
    check_if_valid_tree,
)

from cortado_core.subprocess_discovery.subtree_mining.folding_label import fold_loops

def min_sub_mining(
    treebank: Mapping[int, TreeBankEntry],
    frequency_counting_strat: FrequencyCountingStrategy,
    k_it,
    min_sup,
    loop=False,
):

    """ """

    if loop:
        fold_loops(treebank, loop, fSets)

    fSets = ct_compute_frequent_activity_sets(
        treebank, frequency_counting_strat, min_sup
    )

    F = generate_initial_candidates(
        treebank, min_sup, frequency_counting_strat, fSets
    )

    # No Pattern with fallthroughs of size 2 => No patterns with Fallthrougs
    has_fallthroughs = any([_contains_fallthrough(f.tree) for f in F])

    # Store the results
    k_pattern = {2: F}

    pSets = _get_prune_sets(fSets, F)

    skipPrune = True

    # For every k > 2 create the k pattern from the frequent k-1 pattern
    for k in range(k_it):

        newF = []

        for tp in F:

            # Compute the right most path extension of all k-1 pattern
            tps = tp.right_most_path_extension(pSets, skipPrune, has_fallthroughs)
            
            sup_to_gain = tp.support

            for c in tps:
                if f := c.update_rmo_list(
                    treebank, min_sup, frequency_counting_strat, sup_to_gain
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
            
    return k_pattern


if __name__ == "__main__":
    from pm4py.objects.log.importer.xes.importer import apply as xes_import
    from cortado_core.utils.cvariants import get_concurrency_variants
    from cortado_core.subprocess_discovery.subtree_mining.treebank import (
        create_treebank_from_cv_variants,
    )
    
    from cortado_core.subprocess_discovery.subtree_mining.maximal_connected_components.maximal_connected_check import set_maximaly_closed_patterns
    from cortado_core.experiments.subpattern_eval.Algos.asai_performance import min_sub_mining_asai
    import timeit
    from cortado_core.utils.timestamp_utils import TimeUnit

    freq_strat = FrequencyCountingStrategy.TraceTransaction
    #l = xes_import("C:\\Users\\Michael\\Desktop\\VQL_Evaluation\\Datasets\\BPI_2017_Only_Fallthrough.xes")
    l = xes_import("C:\\Users\\Michael\\Desktop\\reducedBPI2012.xes")
    art_start = False
    variants = get_concurrency_variants(l, False, TimeUnit.MS)
    
    treebank = create_treebank_from_cv_variants(variants, artifical_start=art_start)
    min_sup = 200
    k = 100
    
    start = timeit.default_timer() 
    
    rmo_k_patterns = min_sub_mining(
        treebank,
        frequency_counting_strat=freq_strat,
        k_it=k,
        min_sup=min_sup,
    )
    
    print('Runtime Valid', timeit.default_timer() - start)
    
    print()
    print("Setting RMO Closed")
    
    set_maximaly_closed_patterns(rmo_k_patterns)
    
    start = timeit.default_timer() 
    
    cm_k_patterns, _ = min_sub_mining_asai(
        treebank,
        frequency_counting_strat=freq_strat,
        k_it=k,
        min_sup=min_sup,
    )
    
    print('Runtime Asai', timeit.default_timer() - start)
    
    print()
    print("Setting CM Closed")
    set_maximaly_closed_patterns(cm_k_patterns) 
    
    print()
    print("Setting RMO Closed")
    set_maximaly_closed_patterns(rmo_k_patterns)
        
    print()
    print(
        "Closed RMO",
        sum(
            [
                len(
                    [
                        pattern
                        for pattern in patterns
                        if pattern.closed and check_if_valid_tree(pattern.tree)
                    ]
                )
                for patterns in rmo_k_patterns.values()
            ]
        ),
    )
    print(
        "Maximal RMO",
        sum(
            [
                len(
                    [
                        pattern
                        for pattern in patterns
                        if pattern.maximal and check_if_valid_tree(pattern.tree)
                    ]
                )
                for patterns in rmo_k_patterns.values()
            ]
        ),
    )
    print()
    print(
        "Closed CM",
        sum(
            [
                len(
                    [
                        pattern
                        for pattern in patterns
                        if pattern.closed and check_if_valid_tree(pattern.tree)
                    ]
                )
                for patterns in cm_k_patterns.values()
            ]
        ),
    )
    print(
        "Maximal CM",
        sum(
            [
                len(
                    [
                        pattern
                        for pattern in patterns
                        if pattern.maximal and check_if_valid_tree(pattern.tree)
                    ]
                )
                for patterns in cm_k_patterns.values()
            ]
        ),
    )

    print()

    rm_k_patterns_nested = {
        k: {str(pattern): pattern for pattern in patterns}
        for k, patterns in rmo_k_patterns.items()
    }
    cm_k_patterns_nested = {
        k: {str(pattern): pattern for pattern in patterns}
        for k, patterns in cm_k_patterns.items()
    }

    for k in cm_k_patterns_nested.keys():

        print()
        print("K:", k)
        print("Total:", "RMO:", len(rmo_k_patterns[k]), "CM:", len(cm_k_patterns[k]))
        print(
            "Valid:",
            "RMO:",
            len(
                [
                    pattern
                    for pattern in rmo_k_patterns[k]
                    if check_if_valid_tree(pattern.tree)
                ]
            ),
            "CM:",
            len(
                [
                    pattern
                    for pattern in cm_k_patterns[k]
                    if check_if_valid_tree(pattern.tree)
                ]
            ),
        )
        print(
            "Closed:",
            "RMO:",
            len(
                [
                    pattern
                    for pattern in rmo_k_patterns[k]
                    if pattern.closed and check_if_valid_tree(pattern.tree)
                ]
            ),
            "CM:",
            len(
                [
                    pattern
                    for pattern in cm_k_patterns[k]
                    if pattern.closed and check_if_valid_tree(pattern.tree)
                ]
            ),
        )
        print(
            "Maxmial:",
            "RMO:",
            len(
                [
                    pattern
                    for pattern in rmo_k_patterns[k]
                    if pattern.maximal and check_if_valid_tree(pattern.tree)
                ]
            ),
            "CM:",
            len(
                [
                    pattern
                    for pattern in cm_k_patterns[k]
                    if pattern.maximal and check_if_valid_tree(pattern.tree)
                ]
            ),
        )

        patterns_rmo = set(rm_k_patterns_nested[k].keys())
        patterns_cm = set(cm_k_patterns_nested[k].keys())

        print()
        print("Intersection")
        for pattern in patterns_cm.intersection(patterns_rmo):

            if check_if_valid_tree(rm_k_patterns_nested[k][pattern].tree):
                if (
                    rm_k_patterns_nested[k][pattern].closed
                    != cm_k_patterns_nested[k][pattern].closed
                ):
                    print(
                        pattern,
                        "Closed RMO",
                        rm_k_patterns_nested[k][pattern].closed,
                        "Closed CM",
                        cm_k_patterns_nested[k][pattern].closed,
                    )

                    print(pattern)
                    print(repr(rm_k_patterns_nested[k][pattern].tree))

                if (
                    rm_k_patterns_nested[k][pattern].maximal
                    != cm_k_patterns_nested[k][pattern].maximal
                ):
                    print(
                        pattern,
                        "Maximal RMO",
                        rm_k_patterns_nested[k][pattern].maximal,
                        "Maximal CM",
                        cm_k_patterns_nested[k][pattern].maximal,
                    )

        print()
        print("Only in RMO")
        for pattern in patterns_rmo.difference(patterns_cm):
            if check_if_valid_tree(rm_k_patterns_nested[k][pattern].tree):
                if (
                    rm_k_patterns_nested[k][pattern].closed
                    or rm_k_patterns_nested[k][pattern].maximal
                ):
                    print(
                        pattern,
                        "Closed RMO",
                        rm_k_patterns_nested[k][pattern].closed,
                        "Maximal RMO",
                        rm_k_patterns_nested[k][pattern].maximal,
                    )

    if patterns_cm.difference(patterns_rmo):
        print("CM has pattern not in RMO")