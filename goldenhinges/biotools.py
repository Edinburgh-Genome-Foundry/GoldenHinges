import itertools as itt

complements = {"A": "T", "T": "A", "C": "G", "G": "C"}


def reverse_complement(sequence):
    """Return the reverse-complement of the DNA sequence.
    For instance ``complement("ATGC")`` returns ``"GCAT"``.

    The sequence must be an ATGC string.
    """
    return "".join([complements[c] for c in sequence[::-1]])


def gc_content(sequence):
    """Return the proportion of G and C in the sequence (between 0 and 1).

    The sequence must be an ATGC string.
    """
    return 1.0 * len([c for c in sequence if c in "GC"]) / len(sequence)


def sequences_differences(seq1, seq2):
    """Return the number of different basepairs between sequences ``seq1``
    and ``seq2`` (which must be ATGC strings)
    """
    return len([c1 for c1, c2 in zip(seq1, seq2) if c1 != c2])


def list_overhangs(overhang_size=4, filters=()):
    """Return the list of all possible ATGC overhangs of the given size, such
    that ``fl(overhang)`` is true for every function ``fl`` in ``filters``.
    """
    return [
        "".join(overhang)
        for overhang in itt.product(* overhang_size * ("ATGC",))
        if all((fl(overhang) for fl in filters))
    ]
