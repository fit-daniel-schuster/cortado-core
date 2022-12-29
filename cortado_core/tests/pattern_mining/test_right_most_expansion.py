from collections import Counter
import timeit
import unittest
from cortado_core.subprocess_discovery.concurrency_trees.cTrees import ConcurrencyTree, cTreeOperator
from cortado_core.subprocess_discovery.subtree_mining.freq_counting import FrequencyCountingStrategy, compute_frequent_activity_sets
from cortado_core.subprocess_discovery.subtree_mining.treebank import create_treebank_from_cv_variants
from cortado_core.tests.pattern_mining.example_log import create_example_log_1, create_example_log_2
from cortado_core.utils.cvariants import get_concurrency_variants
from cortado_core.subprocess_discovery.subtree_mining.right_most_path_extension.min_sub_mining import min_sub_mining


l = create_example_log_1()
variants = get_concurrency_variants(l)
treebank = create_treebank_from_cv_variants(variants, artifical_start = False)
   
class RightMostPathMining(unittest.TestCase):
    
    def test_TraceTransaction_Mining(self):
        strat = FrequencyCountingStrategy.TraceTransaction
        
        k_patterns = min_sub_mining(treebank, strat, 20, 2)
    
        self.assertEqual(len(k_patterns.keys()), 6, 'Max Pattern Size')
        
        self.assertEqual(len(k_patterns[2]), 8, 'Patterns of Size 2')
        self.assertEqual(len(k_patterns[3]), 11, 'Patterns of Size 3')
        self.assertEqual(len(k_patterns[4]), 9, 'Patterns of Size 4')
        self.assertEqual(len(k_patterns[5]), 6, 'Patterns of Size 5')
        self.assertEqual(len(k_patterns[6]), 3, 'Patterns of Size 6')
        self.assertEqual(len(k_patterns[7]), 1, 'Patterns of Size 7')
        
        patterns_3 = [str(x) for x in k_patterns[3]]
        
        self.assertIn('→(B, C) with support: 6', patterns_3, '→(B, C) with support: 6')
        self.assertIn('∧(B, C) with support: 3', patterns_3, '∧(B, C) with support: 3')
        self.assertIn('→(A, B) with support: 5', patterns_3, '→(A, B) with support: 5')
        self.assertIn('∧(G, →( )) with support: 3', patterns_3, '∧(G, →( )) with support: 3')
        
        patterns_4 = [str(x) for x in k_patterns[4]]
        
        self.assertIn('→(∧(B, C)) with support: 3', patterns_4)
        self.assertIn('→(A, ∧(B)) with support: 3', patterns_4)
        self.assertIn('→(H, ∧(G)) with support: 3', patterns_4)
        self.assertIn('∧(G, →(B)) with support: 3', patterns_4)
        
        patterns_5 =  [str(x) for x in k_patterns[5]]
        
        self.assertIn('→(A, ∧(B, C)) with support: 3', patterns_5)
        self.assertIn('→(∧(G, →(B))) with support: 3', patterns_5)
        
        patterns_6 =  [str(x) for x in k_patterns[6]]
        
        self.assertIn('→(H, ∧(G, →(B))) with support: 3', patterns_6)
        self.assertIn('→(∧(G, →(B, C))) with support: 3', patterns_6)
          
    
        self.assertEqual(str(k_patterns[7][0].tree), '→(H, ∧(G, →(B, C)))')
                
                
    def test_TraceOccurences_Mining(self):
        strat = FrequencyCountingStrategy.TraceOccurence
        
        k_patterns = min_sub_mining(treebank, strat, 20, 2)
         
        patterns_3 = [str(x) for x in k_patterns[3]]
        
        self.assertIn('→(A, B) with support: 5', patterns_3, '→(B, C) with support: 6')
        self.assertIn('∧(B, C) with support: 3', patterns_3, '∧(B, C) with support: 3')
        self.assertIn('∧(G, →( )) with support: 3', patterns_3, '∧(G, →( )) with support: 3')
        
        
    def test_VariantTransaction_Mining(self):
        strat = FrequencyCountingStrategy.VariantTransaction
        
        k_patterns = min_sub_mining(treebank, strat, 20, 2)
         
        self.assertEqual(len(k_patterns.keys()), 4, 'Max Pattern Size')
        patterns_3 = [str(x) for x in k_patterns[3]]
        
        self.assertIn('→(A, B) with support: 4', patterns_3, '→(B, C) with support: 6')
        self.assertIn('→(B, C) with support: 4', patterns_3, '∧(B, C) with support: 3')
        
        patterns_4 = [str(x) for x in k_patterns[4]]
        
        self.assertIn('→(∧(B, C)) with support: 3', patterns_4, '→(B, C) with support: 6')

        patterns_5 = [str(x) for x in k_patterns[5]]
        
        self.assertIn('→(A, ∧(B, C)) with support: 3' , patterns_5, '→(A, ∧(B, C)) with support: 3')

    
    def test_VariantOccurence_Mining(self):
        
        strat = FrequencyCountingStrategy.VariantOccurence
        k_patterns = min_sub_mining(treebank, strat, 20, 2)
         
        self.assertEqual(len(k_patterns.keys()), 4, 'Max Pattern Size')
        patterns_3 = [str(x) for x in k_patterns[3]]
        
        self.assertIn('→(A, B) with support: 4', patterns_3, '→(B, C) with support: 6')
        self.assertIn('→(B, C) with support: 4', patterns_3, '∧(B, C) with support: 3')
        
        patterns_4 = [str(x) for x in k_patterns[4]]
        
        self.assertIn('→(∧(B, C)) with support: 3', patterns_4, '→(B, C) with support: 3')

        patterns_5 = [str(x) for x in k_patterns[5]]
        
        self.assertIn('→(A, ∧(B, C)) with support: 3' , patterns_5, '→(A, ∧(B, C)) with support: 3')

if __name__ == '__main__':
    from pm4py.objects.log.importer.xes.importer import apply as xes_import
    min_sup = 0
    
    freq_strat = FrequencyCountingStrategy.TraceTransaction
    
    l = create_example_log_2()
    variants = get_concurrency_variants(l)
    treebank = create_treebank_from_cv_variants(variants, artifical_start = False)
    
    strat = FrequencyCountingStrategy.TraceTransaction

    k_patterns = min_sub_mining(treebank, strat, 20, 1)
    
    def _check_seq_predecessor(t): 
        
        node = t.rml
        while node.parent != None: 
            if node.parent.op == cTreeOperator.Sequential: 
                if node.parent.children[0] != node:
                        return False

            node = node.parent 
        
        return True

    for k in k_patterns: 
        print()
        print('K: ', k)
        
        for p in k_patterns[k]:
            print(p, _check_seq_predecessor(p))
                        
                        
                        
                        
    l = create_example_log_1()
   
 
    unittest.main()
    
