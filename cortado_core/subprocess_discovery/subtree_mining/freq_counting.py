from collections import Counter, defaultdict

from typing import List, Mapping
from pm4py.objects.log.obj import Trace
from cortado_core.subprocess_discovery.concurrency_trees.cTrees import cTreeOperator
from cortado_core.subprocess_discovery.subtree_mining.obj import (
    FrequencyCountingStrategy,
    FrequentActivitySets,
    PruningSets,
)
from cortado_core.subprocess_discovery.subtree_mining.tree_pruning import (
    compute_f2_pruned_set,
)

from cortado_core.utils.constants import ARTIFICAL_END_NAME, ARTIFICAL_START_NAME
from cortado_core.utils.split_graph import Group


def compute_frequent_activity_sets(
    variants: Mapping[Group, List[Trace]],
    freq_strat: FrequencyCountingStrategy,
    min_sup: int,
):

    """
    Computes the sets of frequent activites in a single pass over the variants, takes into account the different frequency counting strategies

    Args:
        variants (Mapping[Group, List[Trace]]): The variants object as created in cvariants
        freq_strat FrequencyCountingStrategy: Frequency Counting Strategy
        min_sup int: The minimal support

    Returns:
        _type_: _description_
    """

    directly_follows_counter = Counter()
    eventually_follows_counter = Counter()
    concurrent_counter = Counter()
    start_activities = Counter()
    end_activities = Counter()
    activities_counter = Counter()

    for variant, ts in variants.items():
        
        df_C = defaultdict(list)
        c_C = defaultdict(list)
        ef_C = defaultdict(list)
        sa_C = defaultdict(list)
        ea_C = defaultdict(list)
        act_C = defaultdict(list)
       
            
        for graph, _ in variant.graphs.items(): 

            # Count the number of Edges, Concurrency Pairs are already sorted
            if (
                freq_strat == FrequencyCountingStrategy.VariantTransaction or 
                freq_strat == FrequencyCountingStrategy.TraceTransaction
            ):
                
                for k in graph.directly_follows.keys(): 
                    df_C[k].append(1)
                for k in graph.concurrency_pairs.keys(): 
                    c_C[k].append(1)
                for k in graph.follows.keys(): 
                    ef_C[k].append(1)
                for k in graph.start_activities.keys(): 
                    sa_C[k].append(1)
                for k in graph.end_activities.keys(): 
                    ea_C[k].append(1)
                for k in graph.events.keys(): 
                    act_C[k].append(1)
                
            # Count the number of Edges, Concurrency Pairs are already sorted
            if (
                freq_strat == FrequencyCountingStrategy.TraceOccurence or 
                freq_strat == FrequencyCountingStrategy.VariantOccurence 
            ):
                
                for k, v in graph.directly_follows.items(): 
                    df_C[k].append(len(set([e[0] for e in v])))
                for k, v in graph.concurrency_pairs.items():
                    c_C[k].append(len(set([e[0] for e in v])))
                for k, v in graph.follows.items(): 
                    ef_C[k].append(len(set([e[0] for e in v])))
                for k, v in graph.start_activities.items(): 
                    sa_C[k].append(len(v))
                for k, v in graph.end_activities.items(): 
                    ea_C[k].append(len(v))
                for k, v in graph.events.items(): 
                    act_C[k].append(len(v))
                
        if (
            freq_strat == FrequencyCountingStrategy.VariantTransaction or 
            freq_strat == FrequencyCountingStrategy.TraceTransaction
        ): 
            df_C = { k : 1 for k, v in df_C.items() if len(v) == len(variant.graphs)}
            c_C = { k : 1 for k, v in c_C.items() if len(v) == len(variant.graphs)}
            ef_C = { k : 1 for k, v in ef_C.items() if len(v) == len(variant.graphs)}
            
            sa_C = { k : 1 for k, v in sa_C.items() if len(v) == len(variant.graphs)}
            ea_C = { k : 1 for k, v in ea_C.items() if len(v) == len(variant.graphs)}
            act_C = { k : 1 for k, v in act_C.items() if len(v) == len(variant.graphs)}
        
        if (
            freq_strat == FrequencyCountingStrategy.VariantOccurence or 
            freq_strat == FrequencyCountingStrategy.TraceOccurence
        ): 
            
            df_C = { k : min(v) for k, v in df_C.items() if len(v) == len(variant.graphs)}
            c_C = { k : min(v) for k, v in c_C.items() if len(v) == len(variant.graphs)}
            ef_C = { k : min(v) for k, v in ef_C.items() if len(v) == len(variant.graphs)}
            
            sa_C = { k : min(v) for k, v in sa_C.items() if len(v) == len(variant.graphs)}
            ea_C = { k : min(v) for k, v in ea_C.items() if len(v) == len(variant.graphs)}
            act_C = { k : min(v) for k, v in act_C.items() if len(v) == len(variant.graphs)}
        
        if (
            freq_strat == FrequencyCountingStrategy.TraceTransaction
            or freq_strat == FrequencyCountingStrategy.TraceOccurence
            ):
            nT = len(ts)

        else:
            nT = 1
                
        directly_follows_counter.update(
            {key:  count * nT for key, count in df_C.items()}
        )
        
        eventually_follows_counter.update(
            {key: count * nT for key, count in ef_C.items()}
        )
        
        concurrent_counter.update({key: count * nT for key, count in c_C.items()})

        activities_counter.update({key: count * nT for key, count in act_C.items()})
        start_activities.update({key: count *  nT for key, count in sa_C.items()})
        end_activities.update({key: count * nT for key, count in ea_C.items()})

    # Check if the Activites are above a certain support
    frequent_df_pairs = set(
        [
            pair
            for pair in directly_follows_counter
            if directly_follows_counter[pair] > min_sup
        ]
    )
    frequent_ef_pairs = set(
        [
            pair
            for pair in eventually_follows_counter
            if eventually_follows_counter[pair] > min_sup
        ]
    )
    frequent_cc_pairs = set(
        [pair for pair in concurrent_counter if concurrent_counter[pair] > min_sup]
    )

    frequent_activities = set(
        [act for act in activities_counter if activities_counter[act] > min_sup]
    )
    frequent_end_activity = set(
        [act for act in end_activities if end_activities[act] > min_sup]
    )
    frequent_start_activity = set(
        [act for act in start_activities if start_activities[act] > min_sup]
    )

    def flatten_pairs(pairs, both_sides=False):

        """
        Flatten a pair into a {l : rs} dict
        """

        freq_dict = defaultdict(set)

        for l, r in pairs:

            freq_dict[l].add(r)

            if both_sides:
                freq_dict[r].add(l)

        return freq_dict

    frequent_df_relations = flatten_pairs(frequent_df_pairs)
    frequent_cc_relations = flatten_pairs(frequent_cc_pairs, both_sides=True)
    frequent_ef_relations = flatten_pairs(frequent_ef_pairs)

    return FrequentActivitySets(
        fA=frequent_activities,
        fStart=frequent_start_activity,
        fEnd=frequent_end_activity,
        dfR=frequent_df_relations,
        efR=frequent_ef_relations,
        ccR=frequent_cc_relations,
    )


def _update_sets_artifical_start(fSets: FrequentActivitySets):
    fSets.dfR[ARTIFICAL_START_NAME] = fSets.fStart
    fSets.efR[ARTIFICAL_START_NAME] = fSets.fA
    
    fSets.fA.add(ARTIFICAL_START_NAME)
    fSets.fA.add(ARTIFICAL_END_NAME)

    for activity in fSets.fEnd:
        fSets.dfR[activity].add(ARTIFICAL_END_NAME)

    for activity in fSets.fA:
        fSets.efR[activity].add(ARTIFICAL_END_NAME)

    return fSets


def _get_prune_sets(fSets: FrequentActivitySets, F) -> PruningSets:

    # Define the initial Prune Sets
    ftNestPrune, dfNestPrune, ccNestPrune = compute_f2_pruned_set(F)

    ccLabelPrune = fSets.ccR
    ccFollowsPrune = fSets.ccR

    dfLabelPrune = fSets.dfR
    dfFollowsPrune = fSets.dfR

    efFollowsPrune = fSets.efR

    ftLabelPrune = {a: ftNestPrune for a in ftNestPrune}
    operatorPrune = {
        a: set(
            [
                cTreeOperator.Sequential,
                cTreeOperator.Concurrent,
                cTreeOperator.Fallthrough,
            ]
        )
        for a in fSets.fA 
    }
    
    operatorActivityPrune = {}
    operatorActivityPrune[cTreeOperator.Sequential] = fSets.fA # Use fA as dfNest does not take into account leafs only following
    operatorActivityPrune[cTreeOperator.Fallthrough] = fSets.fA
    operatorActivityPrune[cTreeOperator.Concurrent] = fSets.fA

    operatorOperatorPrune = {}
    operatorOperatorPrune[cTreeOperator.Sequential] = set([cTreeOperator.Fallthrough])
    operatorOperatorPrune[cTreeOperator.Fallthrough] = set([cTreeOperator.Sequential, cTreeOperator.Concurrent, cTreeOperator.Fallthrough])
    operatorOperatorPrune[cTreeOperator.Concurrent] = set([cTreeOperator.Concurrent, cTreeOperator.Fallthrough])
    
    
    return PruningSets(
        ftNestPrune=ftNestPrune,
        dfNestPrune=dfNestPrune,
        ccNestPrune=ccNestPrune,
        ccLabelPrune=ccLabelPrune,
        ccFollowsPrune=ccFollowsPrune,
        dfLabelPrune=dfLabelPrune,
        dfFollowsPrune=dfFollowsPrune,
        efFollowsPrune=efFollowsPrune,
        ftLabelPrune=ftLabelPrune,
        operatorPrune=operatorPrune,
        operatorOperatorPrune=operatorOperatorPrune,
        operatorActivityPrune = operatorActivityPrune,
    )
    
    

