import itertools as itt

complements = {"A": "T", "T": "A", "C": "G", "G": "T"}


def rev_complement(sequence):
    """Return the reverse-complement of the DNA sequence.
    For instance ``complement("ATGC")`` returns ``"GCAT"``.
    """
    return "".join([complements[c] for c in sequence[::-1]])


def gc_content(sequence):
    return 1.0 * len([c for c in sequence if c in "GC"]) / len(sequence)


def seq_differences(seq1, seq2):
    return len([c1 for c1, c2 in zip(seq1, seq2) if c1 != c2])


def list_overhangs(overhang_size=4, filters=()):
    return [
        "".join(overhang)
        for overhang in itt.product(* overhang_size * ("ATGC",))
        if all((fl(overhang) for fl in filters))
    ]
