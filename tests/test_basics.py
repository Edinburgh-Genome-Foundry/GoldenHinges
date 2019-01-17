"""
Basic tests to check that the core functionalities are at least running.
"""
import os
import itertools
import numpy as np
from goldenhinges import (OverhangsSelector, list_overhangs, gc_content,
                          sequences_differences, reverse_complement)
from dnachisel import (random_dna_sequence, sequence_to_biopython_record,
                       annotate_record)
import pytest


@pytest.fixture
def data():
    data_path = os.path.join('tests', 'test_data')
    with open(os.path.join(data_path, "phage_sequence.txt"), "r") as f:
        phage_sequence = f.read()
    return {
        "phage_sequence": phage_sequence
    }


def test_generate_overhangs_collection():
    selector = OverhangsSelector(gc_min=0.5, gc_max=0.5,
                                 differences=2, time_limit=2)
    collection = selector.generate_overhangs_set(n_overhangs=18, n_cliques=100)
    collection = selector.generate_overhangs_set(start_at=len(collection))
    assert len(collection) == 24
    for o1, o2 in itertools.combinations(collection, 2):
        assert sequences_differences(o1, o2) >= 2
        assert sequences_differences(o1, reverse_complement(o2)) >= 2


def test_generate_overhangs_collection2():
    selector = OverhangsSelector(gc_min=0.25, gc_max=0.75,
                                 differences=2, time_limit=2)
    collection = selector.generate_overhangs_set()
    assert len(collection) >= 24
    for o1, o2 in itertools.combinations(collection, 2):
        assert sequences_differences(o1, o2) >= 2
        assert sequences_differences(o1, reverse_complement(o2)) >= 2

def test_generate_overhangs_collection_with_possible():
    selector = OverhangsSelector(gc_min=0.25, gc_max=0.75,
                                 differences=1,
                                 possible_overhangs=['ATTC', 'AAAA', 'GAAT',
                                                     'CTCA'],
                                 time_limit=2)
    collection = selector.generate_overhangs_set()
    assert len(collection) == 2

def test_cut_sequence_into_similar_lengths(data):
    def invalid_overhang(overhang):
        gc = gc_content(overhang)
        three_gc = max([gc_content(overhang[:-1]),
                        gc_content(overhang[1:])]) == 1
        return (gc != 0.5) and (three_gc or (gc != 0.75))

    forbidden_overhangs = list_overhangs(filters=[invalid_overhang])

    selector = OverhangsSelector(forbidden_overhangs=forbidden_overhangs,
                                 differences=1, time_limit=2)
    sequence = data["phage_sequence"]
    solution = selector.cut_sequence(
        sequence, equal_segments=50, max_radius=20, include_extremities=False)
    indices = [o['location'] for o in solution]
    diffs = np.diff([0] + indices + [len(sequence)])
    assert len(diffs) == 50
    assert int(diffs.mean()) == 970

def test_from_record():
    seq = random_dna_sequence(7202, seed=123)
    record = sequence_to_biopython_record(seq)
    zone = (1900, len(seq) - 1900)
    annotate_record(record, location=zone,
                    label="Gene: acs", color='#8edfff')
    annotate_record(record, location=zone,
                    label="@EnforceTranslation")
    annotate_record(record, location=(zone[0]-1800, zone[0], 0),
                    label="@AvoidChanges")
    annotate_record(record, location=(zone[1], 1800 + zone[1], 0),
                    label="@AvoidChanges")

    # ADD SEMI-RANDOM CUTTING ZONES
    cut_region_size = 70
    zones = [
        (x + int(200*np.sin(x)),
         x + cut_region_size + int(200*np.sin(x) - 50*np.cos(x)),
         0)
        for x in range(50, len(seq), 1030)[1:]
    ]
    for zone in zones:
        annotate_record(record, location=zone, label="!cut")

    # SOLVE PROBLEM

    selector = OverhangsSelector(gc_min=0.25, gc_max=0.75, differences=2)
    solution = selector.cut_sequence(record, allow_edits=True,
                                     include_extremities=True)
    assert (solution is not None)
