import itertools as itt
from copy import deepcopy
from functools import lru_cache
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import DNAAlphabet
import numpy as np


complements = {"A": "T", "T": "A", "C": "G", "G": "C"}
def reverse_complement(sequence):
    """Return the reverse-complement of the DNA sequence.
    For instance ``complement("ATGC")`` returns ``"GCAT"``.

    The sequence must be an ATGC string.
    """
    return "".join([complements[c] for c in sequence[::-1]])

@lru_cache(maxsize=4096)
def memo_reverse_complement(sequence):
    return reverse_complement(sequence)

def gc_content(sequence):
    """Return the proportion of G and C in the sequence (between 0 and 1).

    The sequence must be an ATGC string.
    """
    return 1.0 * len([c for c in sequence if c in "GC"]) / len(sequence)


# def sequences_differences(seq1, seq2):
#     """Return the number of different basepairs between sequences ``seq1``
#     and ``seq2`` (which must be ATGC strings)
#     """
#     return len([c1 for c1, c2 in zip(seq1, seq2) if c1 != c2])

def sequences_differences_array(seq1, seq2):
    """Return an array [0, 0, 1, 0, ...] with 1s for sequence differences.

    seq1, seq2 should both be ATGC strings.
    """
    if len(seq1) != len(seq2):
        raise ValueError("Only use on same-size sequences (%d, %d)" %
                         (len(seq1), len(seq2)))
    arr1 = np.fromstring(seq1, dtype="uint8")
    arr2 = np.fromstring(seq2, dtype="uint8")
    return arr1 != arr2


def sequences_differences(seq1, seq2):
    """Return the number of nucleotides that differ in the two sequences.

    seq1, seq2 should be strings of DNA sequences e.g. "ATGCTGTGC"
    """
    return int(sequences_differences_array(seq1, seq2).sum())


def list_overhangs(overhang_size=4, filters=()):
    """Return the list of all possible ATGC overhangs of the given size, such
    that ``fl(overhang)`` is true for every function ``fl`` in ``filters``.
    """
    return [
        "".join(overhang)
        for overhang in itt.product(* overhang_size * ("ATGC",))
        if all((fl(overhang) for fl in filters))
    ]


def crop_record(record, crop_start, crop_end, features_suffix=" (part)"):
    """Return the cropped record with possibly cropped features.

    Note that this differs from ``record[start:end]`` in that in the latter
    expression, cropped features are discarded.

    Parameters
    ----------

    record
      A Biopython record

    crop_start, crop_end
      Start and end of the segment to be cropped.

    features_suffix
      All cropped features will have their label appended with this suffix.
    """
    features = []
    for feature in record.features:
        start, end = sorted([feature.location.start, feature.location.end])
        new_start, new_end = max(start, crop_start), min(end, crop_end)
        if new_end <= new_start:
            continue
        new_start, new_end = new_start - crop_start, new_end - crop_start

        feature = deepcopy(feature)
        feature.location = FeatureLocation(new_start, new_end,
                                           feature.location.strand)
        label = "".join(feature.qualifiers.get("label", ""))
        feature.qualifiers["label"] = label + features_suffix
        features.append(feature)

    new_record = record[crop_start: crop_end]
    new_record.features = features
    return new_record


def sequences_differences_segments(seq1, seq2):
    """Return the list of segments on which sequence seq1 differs from seq2.

    The list is of the form [(start1, end1), (start2, end2), etc.]

    Parameters
    ----------

    seq1, seq2
      ATGC sequences to be compared
    """
    arr1 = np.fromstring(seq1, dtype="uint8")
    arr2 = np.fromstring(seq2, dtype="uint8")
    arr = 1 * (arr1 != arr2)
    diffs = np.diff([0] + list(arr) + [0]).nonzero()[0]
    half = int(len(diffs)/2)
    return [(diffs[2*i], diffs[2*i+1]) for i in range(half)]


def annotate_record(seqrecord, location="full", feature_type="misc_feature",
                    margin=0, **qualifiers):
    """Add a feature to a Biopython SeqRecord.

    Parameters
    ----------

    seqrecord
      The biopython seqrecord to be annotated.

    location
      Either (start, end) or (start, end, strand). (strand defaults to +1)

    feature_type
      The type associated with the feature

    margin
      Number of extra bases added on each side of the given location.

    qualifiers
      Dictionnary that will be the Biopython feature's `qualifiers` attribute.
    """
    if location == "full":
        location = (margin, len(seqrecord)-margin)

    strand = location[2] if len(location) == 3 else 1
    seqrecord.features.append(
        SeqFeature(
            FeatureLocation(location[0], location[1], strand),
            qualifiers=qualifiers,
            type=feature_type
        )
    )

def sequence_to_biopython_record(sequence, id='<unknown id>',
                                 name='<unknown name>', features=()):
    return SeqRecord(Seq(sequence, alphabet=DNAAlphabet()),
                     id=id, name=name, features=list(features))
