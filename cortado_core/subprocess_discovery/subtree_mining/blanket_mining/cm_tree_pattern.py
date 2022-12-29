from typing import Mapping, Tuple
from cortado_core.subprocess_discovery.concurrency_trees.cTrees import (
    ConcurrencyTree,
    cTreeOperator,
)
from cortado_core.subprocess_discovery.subtree_mining.blanket_mining.cm_tree import (
    CMConcurrencyTree,
)
from cortado_core.subprocess_discovery.subtree_mining.freq_counting import (
    FrequencyCountingStrategy,
)
from cortado_core.subprocess_discovery.subtree_mining.obj import PruningSets
from cortado_core.subprocess_discovery.subtree_mining.tree_pattern import (
    TreePattern,
    extend_node,
)
from cortado_core.subprocess_discovery.subtree_mining.treebank import TreeBankEntry
from cortado_core.subprocess_discovery.subtree_mining.utilities import (
    _compute_unique_roots,
    _update_current_sup,
    _update_subToGain,
    flatten,
)
from itertools import repeat


class CMTreePattern(TreePattern):
    def __init__(self, T: CMConcurrencyTree, rml: CMConcurrencyTree, heightDiff: int):
        """ """

        self.tree: CMConcurrencyTree = T
        self.parent: CMTreePattern = None
        self.support: int = None
        self.maximal: bool = False
        self.closed: bool = False
        self.nInner: int = None
        self.nLeaf: int = None
        self.rml: CMConcurrencyTree = rml
        self.heightDiff: int = heightDiff
        self.rightBlanket: CMConcurrencyTree = None
        self.rmo: Mapping[int, Tuple[int, ConcurrencyTree]] = {}

    def __str__(self):
        """ """
        return str(self.tree) + " with support: " + str(self.support)

    def right_most_path_extension(
        self,
        pSets : PruningSets,
        skipPrune : bool,
        has_fallthroughs : bool,
    ):
        """ """

        # print('Right Blanket', self.rightBlanket)

        extended_motifes = []

        # The extension_position as an offset in height from the rml/rmo Node in the pattern / tree
        currentNode = self.rml
        extensionOffset = 0

        # If the node is not closed, we can terminate early as we will only extend at the node until it has at least 2 children
        if currentNode.op and (not skipPrune) and len(currentNode.children) < 2:

            extended_motifes = extend_node(
                self,
                currentNode,
                extensionOffset,
                pSets,
                skipPrune,
                has_fallthroughs,
            )

            # Backtrack from the rmo to the root along the right-most path, but only if the bottom node is closed, i.e it has at least 2 operators
        else:

            while currentNode is not None:

                # Only append to operator Nodes
                if currentNode.op:
                    extended_motifes.extend(
                        extend_node(
                            self,
                            currentNode,
                            extensionOffset,
                            pSets,
                            skipPrune,
                            has_fallthroughs,
                        )
                    )

                    # If the current Operator Node isn't closed, break early
                    if (not skipPrune) and len(currentNode.children) < 2:
                        break

                if currentNode == self.rightBlanket:
                    break

                currentNode = currentNode.parent
                extensionOffset += 1

        return extended_motifes

    def update_rmo_list(
        self,
        treeBank: Mapping[int, TreeBankEntry],
        min_sup: int,
        freq_strat: FrequencyCountingStrategy,
        supToGain: int,
    ):
        """
        Updates the List of Right Most Occurences in the Treebank
        """

        rmo = {}
        cur_sup = 0
        mask_map = {}

        # Check if there is no sequential predecessor, allow then matching all children of a nested Sequential Operator Node
        no_seq_predecessor = _check_no_seq_predecessor(self)

        for oc in self.rmo:

            occurences = []
            mask = []

            for entry in self.rmo[oc]:

                offsetHeight = self.heightDiff
                rid, eNode = entry
                nExpansions = 0
                skipped_right = False 
                
                # Go back to find the parent node
                while offsetHeight != 0:

                    eChild = eNode
                    eNode = eNode.parent
                    offsetHeight -= 1

                    # Remaining offset too walk is larger 1, and we have a sequential operator that does not cover all siblings
                    if offsetHeight != 0 and eNode.op == cTreeOperator.Sequential and eChild != eNode.children[-1]:
                        skipped_right = True
                        break
                    
                if skipped_right:
                    mask.append(nExpansions)
                    continue 

                # We extended somewhere along the right-most path
                if self.heightDiff > 0:

                    # There exists a right sibling
                    if eChild.rSib:

                        if eNode.op == cTreeOperator.Sequential:

                            right_sibling = eChild.rSib

                            # Check if the siblings label is the same as the label of the rml
                            if (
                                self.rml.label is not None
                                and right_sibling.label == self.rml.label
                            ) or (
                                self.rml.op is not None
                                and right_sibling.op == self.rml.op
                            ):
                                occurences.append((rid, right_sibling))
                                nExpansions += 1

                        else:

                            right_sibling = eChild.rSib

                            while right_sibling is not None:

                                if (
                                    self.rml.op is not None
                                    and self.rml.op == right_sibling.op
                                ):
                                    occurences.append((rid, right_sibling))
                                    nExpansions += 1

                                elif (
                                    self.rml.label is not None
                                    and self.rml.label == right_sibling.label
                                ):
                                    occurences.append((rid, right_sibling))
                                    nExpansions += 1
                                    break  # Stop early on the Left most child, TODO Consider the Case of an S Label

                                right_sibling = right_sibling.rSib

                # We extened the "leaf node", thus we need to check all the children of the rmo in the tree, if they match the
                # newly added node under the "leaf node"
                else:
                     
                    if no_seq_predecessor or eNode.op != cTreeOperator.Sequential: 
                        for child in eNode.children:
                            if (
                                self.rml.label is not None and child.label == self.rml.label
                            ) or (self.rml.op is not None and child.op == self.rml.op):
                                occurences.append((rid, child))
                                nExpansions += 1
                            else:
                                continue
                    else: 
                        lmc = eNode.children[0] 
                        if (
                                self.rml.label is not None and lmc.label == self.rml.label
                            ) or (self.rml.op is not None and lmc.op == self.rml.op):
                                occurences.append((rid, lmc))
                                nExpansions += 1

                mask.append(nExpansions)

            assert sum(mask) == len(occurences)
            assert len(mask) == len(self.rmo[oc])

            # Add to the rmo dict the occurences:
            if len(occurences) > 0:
                mask_map[oc] = mask

                rmo[oc] = occurences
                cur_sup = _update_current_sup(
                    treeBank, freq_strat, cur_sup, oc, _compute_unique_roots(occurences)
                )  # TODO, check such that the latter computation is only done for Occurence Settings

            # If the support to gain is smaller than the gap between current support and min_sup we can stop early
            supToGain = _update_subToGain(
                treeBank, freq_strat, supToGain, oc, self.rmo[oc]
            )

            if min_sup > (supToGain + cur_sup):
                del self
                return False

        # Keep the pattern and update the parents closed/maximal state if needed
        if cur_sup > min_sup:
            self.rmo = rmo
            update_occurences_list(self.tree, mask_map)
            self.rml.occList = {k: [v[1] for v in vs] for k, vs in rmo.items()}
            self.support = cur_sup

            return self

        else:
            del self
            return False



def _check_no_seq_predecessor(t : TreePattern): 
        
    node = t.rml
    
    while node.parent != None: 
        
        if node.parent.op == cTreeOperator.Sequential: 
            if node.parent.children[0] != node:
                     return False

        node = node.parent 
        
    return True


def update_occurences_list(tree : CMConcurrencyTree, maskMap):

    occList = {}

    if tree.occList:
    
        for occ in tree.occList:

            if occ in maskMap:

                occList[occ] = flatten(
                    [list(repeat(x, y)) for x, y in zip(tree.occList[occ], maskMap[occ])]
                )
                
        tree.occList = occList

        for child in tree.children:
            update_occurences_list(child, maskMap)
