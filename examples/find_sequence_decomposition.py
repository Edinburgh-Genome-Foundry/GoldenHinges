from goldenhinges import OverhangsSelector, list_overhangs

with open("data/phage_sequence.txt", "r") as f:
    sequence = f.read()

forbidden_overhangs = 

selector = OverhangsSelector(gc_min=0.25, gc_max=0.75, differences=1,
                             time_limit=2)
indices = selector.cut_sequence_into_similar_lengths(sequence, nsegments=50,
                                                     max_radius=20,
                                                     include_extremities=False)
print "solution", indices
import numpy
diffs = numpy.diff([0] + indices + [len(sequence)])
print len(diffs), diffs.min(), diffs.max(), numpy.mean(diffs), numpy.std(diffs)