from goldenhinges import OverhangsSelector, list_overhangs
import numpy

with open("data/phage_sequence.txt", "r") as f:
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
