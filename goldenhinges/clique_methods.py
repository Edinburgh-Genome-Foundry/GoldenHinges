import networkx
import numpy as np
from tqdm import tqdm
import itertools
from .biotools import (sequences_differences, gc_content, list_overhangs,
                        memo_reverse_complement)

def find_large_compatible_subset(all_elements, mandatory_elements=(),
                                 compatibility_conditions=(),
                                 elements_filters=(),
                                 solution_validity_conditions=(),
                                 score="subset_size",
                                 n_solutions_considered=5000,
                                 progress_bar=False,
                                 randomize=False):
    """Return a maximal subset of `all_elements` where all elements are valid
    and inter-compatibles.

    The algorithm takes in a set of elements `all_elements` and filters out
    some elements using the filters in `elements_filters`. Then it creates
    a graph whose nodes are the remaining elements. Edges are added between
    all pairs of elements which are "compatible" as defined by the
    `compatibility_conditions`. Finally we look for the "cliques" i.e. subsets
    of the graph made of elements that are all inter-compatibles. We consider
    a number `n_solutions_considered` of these, and return the one which scored
    highest as defined by the `score` function.

    Parameters
    ----------

    all_elements
      A tuple or list of all possible elements.

    mandatory_elements
      A tuple or list of elements contained in `all_elements` that must be
      included in the final solution.

    elements_filters
      Functions used to pre-filter the list of `all_elements`. Must be a
      list or tuple of functions `fun(element)->True/False`. Only elements
      of `all_elements` such that `fun(element)` is True for all filters
      are considered.

    compatibility_conditions
      Functions used to determine whether two elements are compatible, i.e.
      whether a solution can feature these two elements at the same time.
      Must be a list or tuple of functions `fun(e1, e2)->True/False`.
      Two elements of `all_elements` are compatible when `fun(e1, e2)`
      is True for all functions in `compatibility_conditions`.

    solution_validity_conditions
      Additional validity conditions used to filter out some solutions.
      Must be a list or tuple of functions `fun(solution)->True/False`, where
      the `solution` is a list of elements of `all_elements`.
      A solution is considered valid when `fun(solution)` is True for all
      functions in `solution_validity_conditions`.

    score
      Of all the solutions explored, the which scores highest is returned.
      Can be either a function `fun(solution)->float` where `solution` is
      a list of elements, or it can be the default "subset_size" which means
      the score will be the length of the solution found (the solution returned
      will have as many elements as could be found)

    n_solutions_considered
      Number of cliques of the graph that are itered through (set to None if
      you want to consider all cliques in the graph, which may take a very
      long time)

    progress_bar
      If True, progress bars are displayed as the edges are computed and the
      graph cliques are explored.

    """

    def progress(iterator, title=None, total_loops=None):
        """Wrap the iterator into a tqdm progress bar"""
        if progress_bar:
            return tqdm(iterator, desc=title, total=total_loops)
        else:
            return iterator

    if score == "subset_size":
        score = lambda subset: len(subset)

    def elements_are_compatible(e1, e2):
        """Gather all compatibility conditions into a single function"""
        return all(
            test(e1, e2)
            for test in compatibility_conditions
        )

    def compatible_with_all_mandatory(element):
        """Test whether an element is compatible with all (other) mandatory
        elements"""
        return all(
            elements_are_compatible(element, m)
            for m in mandatory_elements
            if (m != element)
        )

    elements_filters += [compatible_with_all_mandatory]

    def element_is_valid(element):
        """Gather all elements pre-filters into a single function."""
        return all(fl(element) for fl in elements_filters)

    def solution_is_valid(elements_set):
        """Gather all validity conditions into a single function."""
        return all(test(elements_set) for test in solution_validity_conditions)

    filtered_elements = list(filter(element_is_valid, all_elements))

    wrong_mandatories = [
        m for m in mandatory_elements
        if m not in filtered_elements
    ]
    if wrong_mandatories != []:
        raise ValueError("These mandatory elements don't pass your tests: " +
                         str(wrong_mandatories))

    graph = networkx.Graph()
    L = len(filtered_elements)
    for (e1, e2) in progress(itertools.combinations(filtered_elements, 2),
                             title="Computing graph edges",
                             total_loops=L * (L - 1) / 2):
        if elements_are_compatible(e1, e2):
            graph.add_edge(e1, e2)

    best_solution = []
    best_score = None
    if randomize:

        ind_to_oh = {i: oh for i, oh in enumerate(graph.nodes())}
        oh_to_ind = {oh: i for i, oh in ind_to_oh.items()}
        graph = networkx.Graph([
            (oh_to_ind[n1], oh_to_ind[n2])
            for (n1, n2) in graph.edges()
        ])
        def cliques_generator():
            """Return random cliques from the connections graph.

            The randomness is introduced by permuting the nodes names,
            running `networkx.circular_paths` once, permuting the nodes
            names again, etc.
            """
            # print (original_adj)
            while True:
                permutation = np.arange(len(graph.nodes()))
                np.random.shuffle(permutation)
                antipermutation = np.argsort(permutation)
                new_graph = networkx.Graph([
                    (permutation[n1], permutation[n2])
                    for (n1, n2) in graph.edges()
                ])
                for clique in networkx.find_cliques(new_graph):
                    clique = [antipermutation[i] for i in clique]
                    clique = [ind_to_oh[i] for i in clique]
                    yield clique
                    break
        cliques_iterator = cliques_generator()
    else:
        cliques_iterator = networkx.find_cliques(graph)
    for i, clique in progress(enumerate(cliques_iterator),
                              title="Exploring graph cliques",
                              total_loops=n_solutions_considered):
        if solution_is_valid(clique):
            clique_score = score(clique)
            if (best_score is None) or (clique_score > best_score):
                best_solution = clique
                best_score = clique_score
        if ((n_solutions_considered is not None) and
                (i > n_solutions_considered)):
            break

    return best_solution


def find_compatible_overhangs(overhangs_size=4,
                              mandatory_overhangs=(),
                              forbidden_overhangs=(),
                              min_gc_content=0, max_gc_content=1,
                              min_overhangs_differences=2,
                              min_reverse_overhangs_differences=2,
                              elements_filters=(),
                              compatibility_conditions=(),
                              solution_validity_conditions=(),
                              n_solutions_considered=5000,
                              score="subset_size",
                              progress_bar=False,
                              randomize=False):
    """Return a list of compatible overhangs for (Golden Gate) assembly,
    satisfying all the specified conditions.

    Parameters
    ----------

    overhangs_size
      The size of the overhangs. Four is the default and the most common case.

    mandatory_overhangs
      A list `["ATGC", "TTGC"...]` of the overhangs that must be part of the
      final solution. If these do not respect the other conditions or if they
      are not compatible between themselves, an error will be raised.

    forbidden_overhangs
      A list `["ATGC", "TTGC"...]` of overhangs that should NOT be part of the
      final solution.

    min_gc_content
      Float between 0.0 and 1.0 indicating the minimum proportion of G and C
      that valid overhangs should contain.

    max_gc_content
      Float between 0.0 and 1.0 indicating the maximum proportion of G and C
      that valid overhangs should contain.

    min_overhangs_differences
      Minimal number of different basepairs between two overhangs for them
      to be compatible (1 is an acceptable value but 2 is advised to really
      ensure the specificity of the assembly).

    min_reverse_overhangs_differences=2
      Minimal number of different basepairs between an overhang and the
      reverse-complement of a second overhang for these two overhangs to be
      to be compatible (1 is an acceptable value but 2 is advised to really
      ensure the specificity of the assembly).

    elements_filters
      Additional filters to narrow down the possible overhangs. Must be a
      list or tuple of functions `fun(element)->True/False`. Only overhangs
      such that `fun(element)` is True for all filters are kept.

    compatibility_conditions
      Additional conditions to determine whether two overhands are compatible,
      i.e. whether a solution can feature these two elements at the same time.
      Must be a list or tuple of functions `fun(e1, e2)->True/False`.
      Two overhangs are compatible when `fun(e1, e2)` is True for all functions
      in `compatibility_conditions`.

    solution_validity_conditions
      Additional validity conditions used to filter out some solutions.
      Must be a list or tuple of functions `fun(solution)->True/False`, where
      the `solution` is a list of elements of `all_elements`.
      A solution is considered valid when `fun(solution)` is True for all
      functions in `solution_validity_conditions`.

    score
      Of all the solutions explored, the which scores highest is returned.
      Can be either a function `fun(solution)->float` where `solution` is
      a list of overhangs, or it can be the default "subset_size" which means
      the score will be the length of the solution found (the solution returned
      will have as many different overhangs as could be found)

    n_solutions_considered
      Number of solutions considered (set to None if you want to consider all
      possible solutions which may take a very long time)

    progress_bar
      If True, progress bars are displayed as the edges are computed and the
      graph cliques are explored. (see `find_best_compatible_subset`)
    """
    all_elements = list_overhangs(overhangs_size)
    elements_filters = list(elements_filters) + [
        lambda e: e not in forbidden_overhangs,
        lambda e: min_gc_content <= gc_content(e) <= max_gc_content,
        lambda e: (sequences_differences(e, memo_reverse_complement(e)) >=
                   min_reverse_overhangs_differences)
    ]

    compatibility_conditions = list(compatibility_conditions) + [
        lambda e1, e2: (sequences_differences(e1, e2) >=
                        min_overhangs_differences),
        lambda e1, e2: (
            sequences_differences(e1, memo_reverse_complement(e2)) >=
            min_reverse_overhangs_differences)
    ]

    return find_large_compatible_subset(
        all_elements, mandatory_elements=mandatory_overhangs,
        elements_filters=elements_filters,
        compatibility_conditions=compatibility_conditions,
        solution_validity_conditions=solution_validity_conditions,
        n_solutions_considered=n_solutions_considered,
        score=score,
        progress_bar=progress_bar,
        randomize=randomize)
