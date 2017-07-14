"""
Basic tests to check that the core functionalities are at least running.
"""

import os
import numpy as np
from goldenhinges import OverhangsSelector, list_overhangs, gc_content
import pytest


@pytest.fixture
def data():
    data_path = os.path.join('tests', 'test_data')
    with open(os.path.join(data_path, "phage_sequence.txt"), "r") as f:
        phage_sequence = f.read()
    return {
        "phage_sequence": phage_sequence
    }


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
    indices = selector.cut_sequence_into_similar_lengths(
        sequence, nsegments=50, max_radius=20, include_extremities=False)
    diffs = np.diff([0] + indices + [len(sequence)])
    assert len(diffs) == 50
    assert int(diffs.mean()) == 970


def test_generate_overhangs_collection():
    selector = OverhangsSelector(gc_min=0.25, gc_max=0.75,
                                 differences=2, time_limit=2)
    collection = selector.generate_overhangs_set(n_overhangs=18)

    assert len(collection) == 18
