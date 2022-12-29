from typing import List, Mapping
from cortado_core.subprocess_discovery.concurrency_trees.cTrees import cTreeOperator
from cortado_core.subprocess_discovery.subtree_mining.blanket_mining.compute_root_occurence_blanket import (
    check_root_occurence_blanket,
)
from cortado_core.subprocess_discovery.subtree_mining.ct_frequency_counting import ct_compute_frequent_activity_sets
from cortado_core.subprocess_discovery.subtree_mining.obj import PruningSets
from cortado_core.subprocess_discovery.subtree_mining.tree_pruning import (
    compute_f3_pruned_set,
)
from cortado_core.subprocess_discovery.subtree_mining.treebank import TreeBankEntry
from cortado_core.subprocess_discovery.subtree_mining.utilities import (
    _contains_fallthrough,
)   

from cortado_core.subprocess_discovery.subtree_mining.maximal_connected_components.maximal_connected_check import (
    check_if_valid_tree
)

from cortado_core.subprocess_discovery.subtree_mining.blanket_mining.cm_tree_pattern import (
    CMTreePattern,
)
from cortado_core.subprocess_discovery.subtree_mining.blanket_mining.compute_frequency_blanket import (
    check_frequency_blanket,
)
from cortado_core.subprocess_discovery.subtree_mining.blanket_mining.compute_occurence_blanket import (
    check_occ_blanket,
)
from cortado_core.subprocess_discovery.subtree_mining.blanket_mining.compute_transaction_blanket import (
    check_transaction_blanket,
)
from cortado_core.subprocess_discovery.subtree_mining.blanket_mining.create_initial_candidates import (
    generate_initial_candidates,
)
from cortado_core.subprocess_discovery.subtree_mining.freq_counting import (
    FrequencyCountingStrategy,
    _get_prune_sets,
)
from cortado_core.utils.split_graph import Group
from cortado_core.subprocess_discovery.subtree_mining.folding_label import fold_loops


def cm_min_sub_mining(
    treebank : Mapping[int, TreeBankEntry],
    frequency_counting_strat: FrequencyCountingStrategy,
    k_it : int,
    min_sup : int,
    loop : int = 0,
):

    if loop:
        fold_loops(treebank, loop, fSets)
    
    fSets = ct_compute_frequent_activity_sets(treebank, frequency_counting_strat, min_sup)
    
    C = generate_initial_candidates(
        treebank, min_sup, frequency_counting_strat, fSets
    )
    
    k_pattern = {2: C}

    has_fallthroughs = any([_contains_fallthrough(f.tree) for f in C])

    # Define the initial Prune Sets
    pSets = _get_prune_sets(fSets, C)

    # Skip the initial pruning step to compute the full set of F3 Patterns too properly update the Pruning Sets
    skipPrune = True

    for k in range(k_it):
        E = []

        for c in C:

            if res := cm_grow(
                c,
                treebank,
                min_sup,
                frequency_counting_strat,
                skipPrune,
                pSets,
                has_fallthroughs,
            ):
                E.extend(res)

            c.tree

        if len(E) == 0:
            break

        C = E

        if len(E) > 0:
            k_pattern[k + 3] = E
        else:
            break

        if k == 0: 
            pSets, C = compute_f3_pruned_set(pSets, C)
            skipPrune = False

    return k_pattern


def cm_grow(
    tp: CMTreePattern,
    treebank,
    min_sup: int,
    frequency_counting_strat: FrequencyCountingStrategy,
    skipPrune: bool,
    pSets: PruningSets,
    has_fallthroughs: bool,
):
    
    occurenceBased = (
        frequency_counting_strat == FrequencyCountingStrategy.TraceOccurence
        or frequency_counting_strat == FrequencyCountingStrategy.VariantOccurence
    )
    
    E = []

    B_left_occ_not_empty, B_occ_not_empty = check_occ_blanket(tp)

    if (not skipPrune) and B_left_occ_not_empty:
        return None

    else:

        patterns = tp.right_most_path_extension(pSets, skipPrune, has_fallthroughs)
        sup_to_gain = tp.support

        for e in patterns:
            if p := e.update_rmo_list(
                treebank, min_sup, frequency_counting_strat, sup_to_gain
            ):
                E.append(p)
        

    if not B_occ_not_empty and check_if_valid_tree(tp.tree):

        if (occurenceBased and not check_root_occurence_blanket(tp)) or (
            not occurenceBased and not check_transaction_blanket(tp)
        ):

            tp.closed = True

            # NO VALID EXTENSION EXISTS, thus we can check if it is maximal
            # if len(E) == 0:  Check if any Extension exists, that isn't based on an operator Node => Easier Exclusion
            if not any([p.rml.label for p in E]):

                B_freq_not_empty = check_frequency_blanket(
                    tp=tp,
                    min_sup=min_sup,
                    treeBank=treebank,
                    strategy=frequency_counting_strat,
                )

                if not B_freq_not_empty:
                   tp.maximal = True               

    return E

if __name__ == "__main__":

    from pm4py.objects.log.importer.xes.importer import apply as xes_import
    from cortado_core.utils.cvariants import get_concurrency_variants
    from cortado_core.subprocess_discovery.subtree_mining.treebank import (
        create_treebank_from_cv_variants,
    )
    from cortado_core.subprocess_discovery.subtree_mining.blanket_mining.cm_grow import (
        cm_min_sub_mining,
    )
    from cortado_core.subprocess_discovery.subtree_mining.maximal_connected_components.maximal_connected_check import (
        check_if_valid_tree,
        set_maximaly_closed_patterns,
    )
    from cortado_core.experiments.subpattern_eval.Algos.asai_performance import (
        min_sub_mining_asai,
    )
    
    from cortado_core.utils.timestamp_utils import TimeUnit

    freq_strat = FrequencyCountingStrategy.TraceTransaction

    l = xes_import("C:\\Users\\Michael\\Desktop\\Sepsis Cases - Event Log.xes")
    art_start = False
    variants = get_concurrency_variants(l, False, TimeUnit.MS)
    
    treebank = create_treebank_from_cv_variants(variants, artifical_start=art_start)
    min_sup = 20
    k = 100 
    
    print(min_sup)
    
    cm_k_patterns = cm_min_sub_mining(
        treebank,
        frequency_counting_strat=freq_strat,
        k_it=k,
        min_sup=min_sup,
    )

    print()
    print("Setting CM Closed")
    set_maximaly_closed_patterns(cm_k_patterns)

    rmo_k_patterns, _ = min_sub_mining_asai(
        treebank,
        frequency_counting_strat=freq_strat,
        k_it=k,
        min_sup=min_sup,
    )

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