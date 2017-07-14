import itertools as itt
import Numberjack as nj
import numpy as np
from .clique_methods import find_compatible_overhangs
from .biotools import (reverse_complement, sequences_differences, gc_content,
                       list_overhangs)


class OverhangsSelector:
    """A selector of comppatible overhangs for Golden-Gate assembly and others.

    Parameters
    ----------

    gc_min
      Minimal amount of GC allowed in the overhangs, e.g. 0.25 for 25%
      of GC in the overhang

    gc_max
      Maximal amount of GC allowed in the overhangs, e.g. 0.75 for 75%
      of GC or less in the overhang.

    differences
      Number of nucleotides by which all the selected overhangs and their
      reverse complement should differ.

    overhang_size
      Number of nucleotides for the overhangs, e.g. 4 for golden gate assembly.

    forbidden_overhang
      List of all forbidden overhangs.

    time_limit
      Time in seconds after which the solvers should stop if no solution
      was found yet.
    """

    def __init__(self, gc_min=0, gc_max=1, differences=1, overhangs_size=4,
                 forbidden_overhangs=(), time_limit=None):
        """Initialize the object (see class description)."""
        self.overhangs_size = overhangs_size
        self.gc_min = gc_min
        self.gc_max = gc_max
        self.forbidden_overhangs = set(forbidden_overhangs)
        self.differences = differences
        self._compatible_overhang_pairs_memo = None
        self._precompute_standard_overhangs()
        self.time_limit = time_limit

    def _precompute_standard_overhangs(self):
        """Precompute the standard form of the allowed overhangs for the selector.

        The allowed overhangs are the ones with the right GC content.
        The standard form of an overhang is the overhang or its complement,
        whichever is lexically smaller.
        """
        self.all_overhangs = list_overhangs(self.overhangs_size)
        self.standard_overhangs = set()
        self.standard_overhangs_list = []
        self.overhang_to_number = {
            overhang: index
            for index, overhang in enumerate(self.all_overhangs)
        }
        for overhang in self.all_overhangs:
            reverse = reverse_complement(overhang)
            if (reverse not in self.standard_overhangs) and \
               (self.gc_min <= gc_content(overhang) <= self.gc_max) and \
               (overhang != reverse) and \
               (overhang not in self.forbidden_overhangs):
                self.standard_overhangs.add(overhang)
                self.standard_overhangs_list.append(overhang)

    def _standardize_overhang(self, overhang):
        """Return the standard form of the overhang.

        The result is either: the overhang itself if it is already in
        self.standard_overhangs. The reverse complement of the overhang if the
        reverse is in self.standard_overhangs. If neither the overhang nor the
        reverse are standard (i.e. if the overhang is forbidden by its GC
        content) these None is returned.
        """
        if overhang in self.standard_overhangs:
            return overhang
        else:
            reverse = reverse_complement(overhang)
            if reverse in self.standard_overhangs:
                return reverse
            else:
                return None

    def _compatible_overhangs_pairs(self, two_sided=False):
        """Return the list of all compatible standard overhangs pairs, i.e.
        overhangs whith at least self.differences differences between them.

        This least is only computed once then kept in memor for subsequent
        calls.
        If two_sided is True, the least contains each pair twice, in both
        directions, i.e. (o1, o2) and (o2, o1). This is used when the order
        of the overhangs in the results matters.
        """
        if self._compatible_overhang_pairs_memo is None:
            self._compatible_overhang_pairs_memo = [
                (self.overhang_to_number[o1], self.overhang_to_number[o2])
                for (o1, o2) in itt.combinations(self.standard_overhangs, 2)
                if (sequences_differences(o1, o2) >= self.differences) and
                (sequences_differences(o1, reverse_complement(o2)) >= self.differences)
            ]
        result = self._compatible_overhang_pairs_memo
        if two_sided:
            result = result + [(o2, o1) for o1, o2 in result]
        return result

    def _list_overhangs_in_sequence(self, sequence, annotations=()):
        """Return the list all subsequences of size ``self.overhang_size``."""
        return [
            sequence[i:i + self.overhangs_size]
            for i in range(len(sequence) - self.overhangs_size)
        ]

    def _find_overhang_in_interval(self, sequence, interval, overhang):
        """Find the position in the subsequence (sequence interval) that
        contain the overhang.

        If several positions are found, the position nearest to the interval's
        center is selected.

        Parameters
        ----------

        sequence
          An ATGC string.

        interval
          A tuple (start, stop) representing the segment of the sequence where
          to look for the overhang

        overhang
          An ATGC sequence representing an overhang
        """
        start, end = interval
        region = sequence[start:end]
        overhangs_locations = {}
        center = 0.5 * (start + end)
        for i, ovh in enumerate(self._list_overhangs_in_sequence(region)):
            std_ovh = self._standardize_overhang(ovh)
            if (std_ovh not in overhangs_locations) or \
                (abs(overhangs_locations[std_ovh] - center) > abs(i - center)):
                overhangs_locations[std_ovh] = i
        return overhangs_locations[overhang] + start

    def select_from_sets(self, sets_list, solutions=1):
        """Find compatible overhangs, picking one from each provided set.

        Parameters
        ----------

        sets_list
          A list of either sets or lists of overhangs.

        solutions
          Either 1 for a unique solution, a number k for a list of solutions,
          or "iter" which returns an iterator over all solutions.
        """
        variables = [nj.Variable([self.overhang_to_number[o] for o in _set])
                     for _set in sets_list]
        if self.differences == 1:
            constraints = [nj.AllDiff(variables)]
        else:
            # if all sets are equal then the variables are interchangeable,
            # otherwise the order matters. This is reflected in the symmetry
            # of the allowed pairs of overhangs.
            two_sided = any(
                set1 != set2
                for set1, set2 in zip(sets_list, sets_list[1:])
            )
            possible_pairs = self._compatible_overhangs_pairs(two_sided)
            constraints = [
                nj.Table((v1, v2), possible_pairs)
                for v1, v2 in itt.combinations(variables, 2)
            ]
        model = nj.Model(*constraints)
        solver = model.load("Mistral", variables)
        if self.time_limit is not None:
            solver.setTimeLimit(self.time_limit)
        solver.solve()

        def get_solution():
            if any([v.get_value() is None for v in variables]):
                return None
            else:
                return [self.all_overhangs[v.get_value()] for v in variables]

        if solutions == 1:
            returned = get_solution()
            solver.delete()
            return returned
        else:
            def solutions_iterator():
                while True:
                    solution = get_solution()
                    yield solution
                    if solution is None:
                        break
                    solver.getNextSolution()
            iterator = solutions_iterator()
            if isinstance(solutions, int):
                returned = [next(iterator) for i in range(solutions)]
                solver.delete()
                return returned
            else:
                return iterator

    def cut_sequence_at_intervals(self, sequence, intervals, solutions=1,
                                  allow_sequence_optimization=False):
        """Select compatible-overhangs cut locations, one in each interval.

        Parameters
        ----------

        sequence
          An ATGC string or a Biopython record

        intervals
          solutions=1

        solutions
          If equal to 1, one solution is returned (i.e. a list of cuts).
          If larger than 1, a list of solutions is returned.
          If equal to "iter", an iterator over solutions is returned
        """
        sets_list = [
            set([
                self._standardize_overhang(o)
                for o in self._list_overhangs_in_sequence(sequence[start:end])
                if self._standardize_overhang(o) is not None
            ])
            for (start, end) in intervals
        ]
        if any([len(s) == 0 for s in sets_list]):
            return None

        choices = self.select_from_sets(sets_list, solutions=solutions)

        def get_solution(choices):
            if choices is None:
                return None
            return [
                self._find_overhang_in_interval(sequence, interval, overhang)
                for interval, overhang in zip(intervals, choices)
            ]

        if solutions == 1:
            return get_solution(choices)
        elif isinstance(solutions, int):
            return [get_solution(c) for c in choices]
        else:
            return (get_solution(c) for c in choices)

    def cut_sequence_nearest_from_indices(self, sequence, indices,
                                          max_radius=10,
                                          allow_sequence_optimization=False):
        """Cut a sequence at locations forming intercompatible overhangs,
        each location as close as possible as an index in the provided list.
        """
        if isinstance(max_radius, int):
            max_radius = [max_radius for i in indices]
        largest_max_radius = max(max_radius)
        radius = 0
        while radius < largest_max_radius:
            radius += 1
            intervals = [
                (max(0, i - min(local_max_radius, radius)),
                 i + max(1, min(local_max_radius, radius)))
                for i, local_max_radius in zip(indices, max_radius)
            ]
            solution = self.cut_sequence_at_intervals(sequence, intervals)
            if solution is not None:
                return solution
        return None

    def cut_sequence_into_similar_lengths(self, sequence, nsegments,
                                          max_radius=10,
                                          include_extremities=False):
        """Return a list of compatible cut locations forming segments of
        similar length.

        Parameters
        ----------

        sequence
          An ATGC string of the sequence to assemble.


        nsegments
          Number of segments desired for the assembly

        max_radius
          Maximal radius around each location for the search.

        include_extremities
          If True, the extremities of the sequence are considered cut, i.e.
          the first and last basepairs are considered overhangs.

        """
        cuts_indices = np.linspace(0, len(sequence), nsegments + 1).astype(int)
        if include_extremities:
            max_radius = [0] + [max_radius for i in range(nsegments - 1)] + [0]
        else:
            cuts_indices = cuts_indices[1:-1]
        return self.cut_sequence_nearest_from_indices(sequence, cuts_indices,
                                                      max_radius=max_radius)

    def generate_overhangs_set(self, n_overhangs=None, mandatory_overhangs=(),
                               start_at=2, step=2, message_callback=None,
                               n_cliques=None):
        """Generate a set of compatible overhangs, eg ``{"ATTC", "ATCG", ...}``

        Parameters
        ----------

        n_overhangs
          Size of the desired overhang set. If left to None, the algorithm will
          return the largest set it can find.

        mandatory_overhangs
          Overhangs which must be in the final set.

        step
          Increment to use for the set size when looking for the larget
          possible set (case ``n_overhangs=None``). Note that this should not
          change the final result, but a well-chosen step can improve the
          computations speed several fold

        """
        if message_callback is None:
            message_callback = lambda *a: None

        if n_cliques is not None:
            return find_compatible_overhangs(
                overhangs_size=self.overhangs_size,
                mandatory_overhangs=mandatory_overhangs,
                forbidden_overhangs=self.forbidden_overhangs,
                min_gc_content=self.gc_min, max_gc_content=self.gc_max,
                min_overhangs_differences=self.differences,
                min_reverse_overhangs_differences=self.differences,
                n_solutions_considered=500, score='subset_size',
                progress_bar=False
            )

        if n_overhangs is None:
            n_overhangs = max([start_at, 2, len(mandatory_overhangs) + 1])
            result = None
            while True:
                message_callback("Attempting to find %d overhangs" %
                                 n_overhangs)
                solution = self.generate_overhangs_set(n_overhangs,
                                                       mandatory_overhangs)
                if solution is None:
                    # the step increment went too far, there was no solution,
                    # conduct a finer search from the last increment.
                    for n in range(n_overhangs - step + 1, n_overhangs):
                        message_callback("Attempting refining to %d overhangs"
                                         % n)
                        solution = self.generate_overhangs_set(
                            n, mandatory_overhangs)
                        if solution is None:
                            break
                        result = solution
                    break
                result = solution
                n_overhangs += step

            return result

        L = len(mandatory_overhangs)
        if L > 0:
            standard_mandatory_overhangs = [
                self._standardize_overhang(o)
                if (self._standardize_overhang(o) is not None)
                else o
                for o in mandatory_overhangs
            ]
            new_standard_overhangs = \
                self.standard_overhangs.\
                difference(set(standard_mandatory_overhangs))
            sets_list = ([[o] for o in standard_mandatory_overhangs] +
                         [new_standard_overhangs
                          for i in range(n_overhangs - L)])
            solution = self.select_from_sets(sets_list)
            if solution is not None:
                solution = list(mandatory_overhangs) + solution[L:]
            return solution
        else:
            sets_list = [self.standard_overhangs for i in range(n_overhangs)]
            return self.select_from_sets(sets_list)
