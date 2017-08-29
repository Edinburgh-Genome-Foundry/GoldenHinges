import itertools as itt
import Numberjack as nj
import numpy as np
from .clique_methods import find_compatible_overhangs
from .biotools import (reverse_complement, sequences_differences, gc_content,
                       list_overhangs)
from tqdm import tqdm
try:
    import dnachisel as dc
    DNACHISEL_AVAILABLE = True
except:
    DNACHISEL_AVAILABLE = False

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

    overhangs_size
      Number of nucleotides for the overhangs, e.g. 4 for golden gate assembly.

    forbidden_overhang
      List of all forbidden overhangs.

    time_limit
      Time in seconds after which the solvers should stop if no solution
      was found yet.
    """

    def __init__(self, gc_min=0, gc_max=1, differences=1, overhangs_size=4,
                 forbidden_overhangs=(), time_limit=None,
                 external_overhangs=(), progress_logger=None):
        """Initialize the object (see class description)."""
        self.overhangs_size = overhangs_size
        self.gc_min = gc_min
        self.gc_max = gc_max
        self.differences = differences
        self._compatible_overhang_pairs_memo = None
        self.time_limit = time_limit
        if progress_logger is None:
            def progress_logger(**k):
                pass
        self.progress_logger = progress_logger
        self.external_overhangs = external_overhangs
        self.all_overhangs = list_overhangs(self.overhangs_size)
        for o1 in self.all_overhangs:
            if any(((sequences_differences(o1, o2) < self.differences)
                    for o2 in external_overhangs)):
                forbidden_overhangs.append(o1)
        self.forbidden_overhangs = set(forbidden_overhangs)
        self._precompute_standard_overhangs()


    def _precompute_standard_overhangs(self):
        """Precompute the standard form of the allowed overhangs for the selector.

        The allowed overhangs are the ones with the right GC content.
        The standard form of an overhang is the overhang or its complement,
        whichever is lexically smaller.
        """
        self.standard_overhangs = set()
        self.standard_overhangs_list = []
        self.overhang_to_number = {
            overhang: index
            for index, overhang in enumerate(self.all_overhangs)
        }
        for o in self.all_overhangs:
            reverse = reverse_complement(o)
            if (reverse not in self.standard_overhangs) and \
               (self.gc_min <= gc_content(o) <= self.gc_max) and \
               sequences_differences(o, reverse) >= self.differences and \
               (o not in self.forbidden_overhangs):
                self.standard_overhangs.add(o)
                self.standard_overhangs_list.append(o)

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
                (sequences_differences(o1, reverse_complement(o2)) >=
                 self.differences)
            ]
        result = self._compatible_overhang_pairs_memo
        if two_sided:
            result = result + [(o2, o1) for o1, o2 in result]
        return result

    def _overhangs_are_compatible(self, o1, o2):
        return tuple(sorted([o1, o2])) in self._compatible_overhangs_pairs()


    def select_from_sets(self, sets_list, solutions=1, optimize_score=True):
        """Find compatible overhangs, picking one from each provided set.

        This is the central solver for methods cut_sequence,
        cut_sequence_into_similar_lengths,

        Parameters
        ----------

        sets_list
          A list of either sets or lists of overhangs.

        solutions
          Either 1 for a unique solution, a number k for a list of solutions,
          or "iter" which returns an iterator over all solutions.
        """

        # Todo:
        # Inspect each set, for sets with only one solution, remove the set,
        # eliminate all incompatible overhangs in other sets



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

        if optimize_score:
            scores_variables = []
            for _set, variable in zip(sets_list, variables):
                max_score = max([o['score'] for o in _set.values()])
                score_variable = nj.Variable(0, max_score + 1)
                scores_variables.append(score_variable)

                model.add(nj.Table((variable, score_variable), [
                    (self.overhang_to_number[o], int(_set[o]['score']))
                    for o in _set
                ]))
            total_score = sum(scores_variables)
            model.add(nj.Minimize(total_score))

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

    def _list_overhangs_in_sequence(self, zone, sequence=None,
                                    constraints_problem=None):
        """Return the list all subsequences of size ``self.overhangs_size``."""
        # local_region_size = min(local_region_size, zone[1] - zone[0])
        start, end = zone
        region_size = self.overhangs_size
        if constraints_problem is not None:
            sequence = constraints_problem.sequence
            for i in range(start, max(start + 1, end - region_size)):
                mutation_space = constraints_problem.mutation_space.localized(
                    (i, i + region_size)
                )
                if len(mutation_space.choices) == 0:
                    yield dict(
                        sequence=sequence[i:i + self.overhangs_size],
                        location=i
                    )
                else:
                    location = dc.Location(*mutation_space.choices_span)
                    subsequence = location.extract_sequence(sequence)
                    localized_constraints = [
                        _constraint.localized(location)
                        for _constraint in constraints_problem.constraints
                        if not _constraint.enforced_by_nucleotide_restrictions
                    ]
                    local_problem = dc.DnaOptimizationProblem(
                        sequence=sequence,
                        constraints=localized_constraints,
                        mutation_space=mutation_space
                    )
                    for variant in mutation_space.all_variants(sequence):
                        local_problem.sequence = variant
                        if local_problem.all_constraints_pass():
                            mutated_region = variant[location.start: location.end]
                            j_end = len(location) - self.overhangs_size
                            for j in range(0, j_end + 1):
                                seq = mutated_region[j: j + region_size]
                                yield dict(
                                    sequence=seq,
                                    location=i + j,
                                    n_mutations=sequences_differences(
                                        subsequence, mutated_region),
                                    mutated_region=(i, mutated_region)
                                )
        else:
            for i in range(zone[0], zone[1] - self.overhangs_size):
                yield dict(
                    sequence=sequence[i:i + self.overhangs_size],
                    location=i
                )

    def cut_sequence(self, sequence, intervals=None, solutions=1,
                     allow_edits=False, include_extremities=True,
                     optimize_score=True, edit_penalty=10,
                     equal_segments=None, max_radius=10,
                     target_indices=None):
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

        # FIRST SUBSCENARIO: CUT THE SEQUENCE INTO EQUAL LENGTHS

        if equal_segments is not None:
            target_indices = np.linspace(0, len(sequence), equal_segments + 1)
            target_indices = target_indices.astype(int)[1:-1]

        # SECOND SUBSCENARIO: CUT THE SEQUENCE NEAREST FROM TARGET INDICES

        if target_indices is not None:
            if isinstance(max_radius, int):
                max_radius = [max_radius for i in target_indices]
            largest_max_radius = max(max_radius)
            radius = 0
            self.progress_logger(max_radius=max_radius)
            while radius < largest_max_radius:
                radius += 1
                self.progress_logger(radius=radius)
                intervals = [
                    (max(0, i - min(local_max_radius, radius)),
                     i + max(1, min(local_max_radius, radius)))
                    for i, local_max_radius in zip(target_indices, max_radius)
                ]
                solution = self.cut_sequence(
                    sequence, intervals, optimize_score=optimize_score,
                    include_extremities=include_extremities,
                    allow_edits=allow_edits, edit_penalty=edit_penalty)
                if solution is not None:
                    return solution
            return None

        # MAIN SUBSCENARIO: CUT THE SEQUENCE IN PREDEFINED INTERVALS

        # either these intervals are provided or they are extracted from the
        # sequence-record

        if intervals is None:
            if not hasattr(sequence, 'features'):
                raise ValueError("Provide either intervals or a record with "
                                 "intervals marked as !cut")
            intervals = [
                (int(f.location.start), int(f.location.end))
                for f in sequence.features
                if ''.join(f.qualifiers.get('label', '')) == "!cut"
            ]

        if include_extremities:
            intervals = [(0, self.overhangs_size)] + intervals + [
                (len(sequence) - self.overhangs_size, len(sequence))
            ]

        cst_problem = None
        if hasattr(sequence, 'features'):
            if allow_edits:
                if not DNACHISEL_AVAILABLE:
                    raise ImportError(
                        "Looks like you are trying to use a GoldenHinges"
                        " feature which requires DnaChisel installed")
                cst_problem = dc.DnaOptimizationProblem.from_record(sequence)
                if len(cst_problem.constraints) == 0:
                    cst_problem = None
            sequence = str(sequence.seq)
        sets_list = []
        self.progress_logger(interval_ind=0, n_intervals=len(intervals),
                             message="Now exploring possibilities...")
        for i, (start, end) in enumerate(intervals):

            overhangs_dict = {}
            middle_location = int(0.5*(end + start))
            all_possible_overhangs = self._list_overhangs_in_sequence(
                sequence=sequence, zone=(start, end),
                constraints_problem=cst_problem)
            for o in all_possible_overhangs:
                std_o = self._standardize_overhang(o['sequence'])
                if std_o is None:
                    continue
                o['score'] = abs(o['location'] - middle_location)
                if 'n_mutations' in o:
                    o['score'] += edit_penalty * o['n_mutations']
                is_new = (std_o not in overhangs_dict)
                if is_new or (o['score'] < overhangs_dict[std_o]['score']):
                    overhangs_dict[std_o] = o
            self.progress_logger(interval_ind=i + 1)
            sets_list.append(overhangs_dict)

        if any([len(s) == 0 for s in sets_list]):
            return None

        self.progress_logger(message="Now thinking...")
        choices = self.select_from_sets(sets_list, solutions=solutions,
                                        optimize_score=optimize_score)

        def get_solution(choices):
            if choices is None:
                return None
            return [
                _set[overhang]
                for _set, overhang in zip(sets_list, choices)
            ]

        if solutions == 1:
            return get_solution(choices)
        elif isinstance(solutions, int):
            return [get_solution(c) for c in choices]
        else:
            return (get_solution(c) for c in choices)

    def generate_overhangs_set(self, n_overhangs=None, mandatory_overhangs=(),
                               start_at=2, step=2, n_cliques=None):
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
                self.progress_logger(n_overhangs=n_overhangs)
                solution = self.generate_overhangs_set(n_overhangs,
                                                       mandatory_overhangs)
                if solution is None:
                    # the step increment went too far, there was no solution,
                    # conduct a finer search from the last increment.
                    for n in range(n_overhangs - step + 1, n_overhangs):
                        self.progress_logger(n_overhangs=n)
                        solution = self.generate_overhangs_set(
                            n, mandatory_overhangs)
                        if solution is None:
                            break
                        result = solution
                    break
                result = solution
                n_overhangs += step

            return result

        # Case where a number of overhangs was specified.
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
            solution = self.select_from_sets(sets_list, optimize_score=False)
            if solution is not None:
                solution = list(mandatory_overhangs) + solution[L:]
            return solution
        else:
            sets_list = [self.standard_overhangs for i in range(n_overhangs)]
            return self.select_from_sets(sets_list, optimize_score=False)
