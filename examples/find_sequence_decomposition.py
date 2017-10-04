"""Decompose a 50kb viral sequence into 50 assemblable 1kb fragments.

The decomposition produces 50 compatible 4bp protusions.

Notes that the fragments contain various type-2S sites (BsaI, BsmBI...)
so the whole assembly cannot be done using standard Golden-Gate-like methods,
see this paper for the associated protocol:

https://www.nature.com/articles/srep10655

"""
from goldenhinges import OverhangsSelector
import numpy
import os

with open(os.path.join("data", "phage_sequence.txt"), "r") as f:
    sequence = f.read()

selector = OverhangsSelector(gc_min=0.25, gc_max=0.75, differences=1,
                             time_limit=2)
solution = selector.cut_sequence(sequence, equal_segments=50,
                                 max_radius=20,
                                 include_extremities=False)
indices = [o['location'] for o in solution]
print ("solution:", indices)

diffs = numpy.diff([0] + indices + [len(sequence)])
print len(diffs), diffs.min(), diffs.max(), numpy.mean(diffs), numpy.std(diffs)
