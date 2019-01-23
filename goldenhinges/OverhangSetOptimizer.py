from .biotools import memo_reverse_complement, sequences_differences
from numpy.random import choice
import numpy as np

class OverhangSetOptimizer:

    def __init__(self, set_size, possible_overhangs, external_overhangs,
                 initial_set=None, mutations=1):
        self.set_size = set_size
        self.external_overhangs = external_overhangs
        self.possible_overhangs = possible_overhangs
        if initial_set is not None:
            self.selected_overhangs = initial_set
        else:
            self.selected_overhangs = np.random.choice(
                possible_overhangs, set_size, replace=False)
        self.selected_overhangs = set(self.selected_overhangs)
        self.leftover_overhangs = (set(possible_overhangs)
                                   .difference(self.selected_overhangs))
        self.current_score = self.score(self.selected_overhangs,
                                        self.external_overhangs)
        self.mutations = mutations

    def score(self, overhangs, external_overhangs):
        if not hasattr(self, '_overhangs_scores'):
            all_overhangs = self.possible_overhangs + self.external_overhangs
            self._overhangs_scores = {}
            for i, o1 in enumerate(all_overhangs[:-1]):
                ro1 = memo_reverse_complement(o1)
                for o2 in all_overhangs[i + 1:]:
                    diff1 = sequences_differences(o1, o2)
                    diff2 = sequences_differences(ro1, o2)
                    score = min(diff1, diff2)
                    self._overhangs_scores[(o1, o2)] = score
                    self._overhangs_scores[(o2, o1)] = score
        overhangs = list(overhangs) + list(external_overhangs)
        return 0 if (len(overhangs) == 0) else np.mean([
        self._overhangs_scores[(o1, o2)]
        for i, o1 in enumerate(overhangs)
        for o2 in overhangs[i+1:]
    ])

    def optimize(self, iterations=100):
        for i in range(iterations):
            removed = np.random.choice(list(self.selected_overhangs),
                                       self.mutations, replace=False)
            added = np.random.choice(list(self.leftover_overhangs),
                                     self.mutations, replace=False)
            new_overhangs = (self.selected_overhangs
                                 .difference(removed).union(added))
            new_score = self.score(new_overhangs, self.external_overhangs)
            if new_score > self.current_score:
                self.selected_overhangs.update(added)
                self.selected_overhangs.difference_update(removed)
                self.leftover_overhangs.update(removed)
                self.leftover_overhangs.difference_update(added)
                self.current_score = new_score
