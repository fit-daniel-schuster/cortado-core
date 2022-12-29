from collections import Counter

from typing import Mapping
import unittest

from pm4py.objects.log.obj import Trace
from pm4py.objects.log.util.interval_lifecycle import to_interval
from cortado_core.subprocess_discovery.subtree_mining.freq_counting import FrequencyCountingStrategy, compute_frequent_activity_sets
from cortado_core.subprocess_discovery.subtree_mining.obj import FrequentActivitySets
from cortado_core.subprocess_discovery.subtree_mining.right_most_path_extension.initial_candidate_generation import generate_initial_candidates
from cortado_core.subprocess_discovery.subtree_mining.treebank import create_treebank_from_cv_variants
from cortado_core.tests.pattern_mining.example_log import create_example_log_1

from cortado_core.utils.cgroups_graph import ConcurrencyGroup, cgroups_graph
from cortado_core.utils.cvariants import get_concurrency_variants, unique_activities
from cortado_core.utils.timestamp_utils import TimeUnit

l = create_example_log_1()

class FrequencyCounting(unittest.TestCase):

    def test_generate_graphs(self): 
        
        interval_log = to_interval(l)
        log_renamed, names = unique_activities(interval_log)

        graphs : Mapping[ConcurrencyGroup, Trace] = {}        

        for trace, original_trace in zip(log_renamed, interval_log):
            variant = cgroups_graph(trace, TimeUnit.SEC)
            graphs[variant] = graphs.get(variant, []) + [original_trace]
            
        g1 = list(graphs.keys())[0]
        
        self.assertTrue('A0' in g1.events)
        self.assertTrue('B0' in g1.events)

        g2 = list(graphs.keys())[-1]
        
        self.assertTrue(('H0', 'C0') in g2.follows)
        self.assertTrue(('B0', 'C0') in g2.directly_follows)
        self.assertFalse(('H0', 'C0') in g2.concurrency_pairs) # Not in Log
        self.assertTrue(('B0', 'G0') in g2.concurrency_pairs)

    def test_count_frequencies_trace_occurence(self):
        
        ### Counting Strat
        
        strat = FrequencyCountingStrategy.TraceOccurence    
        
        variants = get_concurrency_variants(l)
        
        frequentActivitySet : FrequentActivitySets = compute_frequent_activity_sets(variants, strat, min_sup = 0)
        
        self.assertTrue(frequentActivitySet.fStart == set(['A', 'I', 'H']))
        self.assertTrue(set(frequentActivitySet.dfR['A']) == set(['B', 'D', 'A', 'C']))
        self.assertTrue(set(frequentActivitySet.ccR.keys()) == set(['A', 'D', 'I', 'B', 'G', 'C']))
        self.assertTrue('I' in frequentActivitySet.ccR['B'])
        self.assertTrue(set(frequentActivitySet.dfR.keys()) == set(['A', 'B', 'C', 'G', 'H']))
        self.assertTrue(set(frequentActivitySet.ccR['C']) == set(['I', 'G', 'D', 'B']))
        
        self.assertFalse('B' in frequentActivitySet.ccR['A']) # Not in Log 
        self.assertFalse('A' in frequentActivitySet.fEnd) # Not in Log

        frequentActivitySet : FrequentActivitySets  = compute_frequent_activity_sets(variants, strat, min_sup = 5)
        
        self.assertTrue(frequentActivitySet.fStart == set(['A']))
        self.assertTrue(set(frequentActivitySet.dfR['A']) == set(['B']) and set(frequentActivitySet.dfR['B']) == set(['C']))
        self.assertTrue(set(frequentActivitySet.ccR.keys()) == set())
        self.assertTrue(set(frequentActivitySet.dfR.keys())  == set(['A', 'B']))
        
        self.assertFalse('D' in frequentActivitySet.dfR['A']) # sup too low
        self.assertFalse('A' in frequentActivitySet.fEnd) # Not in Log
        
        
    def test_count_frequencies_trace_transaction(self):
        
        ### Counting Strat
        
        strat = FrequencyCountingStrategy.TraceTransaction    
        
        variants = get_concurrency_variants(l)
        
        frequentActivitySet : FrequentActivitySets = compute_frequent_activity_sets(variants, strat, min_sup = 3)
        
        self.assertTrue(frequentActivitySet.fStart == set(['A', 'H']))
        self.assertTrue(set(frequentActivitySet.dfR['A']) == set(['B']))
        self.assertTrue(set(frequentActivitySet.dfR['H']) == set(['B']))
        self.assertFalse('A' in frequentActivitySet.fEnd) 
        self.assertFalse('B' in set(frequentActivitySet.ccR.keys())) # Sup 3 
        self.assertFalse('G' in set(frequentActivitySet.ccR.keys())) # Sup 3 
        self.assertFalse('G' in frequentActivitySet.dfR['H']) 
        self.assertFalse('C' in frequentActivitySet.dfR['A']) 
        
        
    def test_count_frequencies_variant_transaction(self):
            
        ### Counting Strat
        
        strat = FrequencyCountingStrategy.VariantTransaction    
        
        variants = get_concurrency_variants(l)
        
        frequentActivitySet : FrequentActivitySets =  compute_frequent_activity_sets(variants, strat, min_sup = 2)
        
        self.assertTrue(set(frequentActivitySet.dfR['A']) == set(['B', 'C']))
        self.assertTrue(set(frequentActivitySet.ccR.keys()) == set(['B', 'C']))
        self.assertTrue(set(frequentActivitySet.dfR['B']) == set('C'))
        self.assertTrue('H' in frequentActivitySet.dfR) # Appears in 3 variants
        
        frequentActivitySet : FrequentActivitySets = compute_frequent_activity_sets(variants, strat, min_sup = 7)

        self.assertFalse('A' in frequentActivitySet.dfR) # Sup A -> B 7 
        
        
        
    def test_count_frequencies_variant_occurence(self):
            
        ### Counting Strat
        
        strat = FrequencyCountingStrategy.VariantOccurence    
        variants = get_concurrency_variants(l)
        
        frequentActivitySet : FrequentActivitySets = compute_frequent_activity_sets(variants, strat, min_sup = 2)
        
        self.assertTrue(set(frequentActivitySet.dfR['A']) == set(['B', 'C']))
        self.assertTrue(set(frequentActivitySet.ccR.keys()) == set(['B', 'C']))
        self.assertTrue(set(frequentActivitySet.dfR['B']) == set('C'))
        self.assertTrue('H' in  frequentActivitySet.dfR) # Appears in 3 Variants
        
        frequentActivitySet : FrequentActivitySets = compute_frequent_activity_sets(variants, strat, min_sup = 7)

        self.assertTrue('A' in frequentActivitySet.dfR) # Sup A -> B 8
        
class InitialCandidateGeneration(unittest.TestCase):
    
    def test_create_candidates_rmo(self):
        
        strat = FrequencyCountingStrategy.TraceTransaction    
        variants = get_concurrency_variants(l)
        
        frequentActivitySet : FrequentActivitySets = compute_frequent_activity_sets(variants, strat, min_sup = 2)
        
        treebank = create_treebank_from_cv_variants(variants, artifical_start = False)
    
        initial_candidates = generate_initial_candidates(treebank,
                2,
                strat, 
                frequentActivitySet,
                )
        
        candidates_dict = {}
        
        for candidate in initial_candidates: 
            candidates_dict[str(candidate.tree)] = candidate.support
            
        self.assertTrue(candidates_dict['→(∧( ))'] == 6)  
        self.assertTrue(candidates_dict['∧(→( ))'] == 4)
        self.assertTrue(candidates_dict['→(A)'] == 8)
        self.assertTrue(candidates_dict['→(H)'] == 5)
        self.assertTrue(candidates_dict['→(B)'] == 10)
        self.assertTrue(candidates_dict['∧(G)'] == 3)
        self.assertTrue(candidates_dict['∧(B)'] == 3)    
