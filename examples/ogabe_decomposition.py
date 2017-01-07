from goldenhinges import OverhangsSelector, list_overhangs, gc_content


with open("data/phage_sequence.txt", "r") as f:
    sequence = f.read()

def invalid_overhang(overhang):
	gc = gc_content(overhang)
	three_gc = max([gc_content(overhang[:-1]), gc_content(overhang[1:])]) == 1
	return (gc != 0.5) and (three_gc or (gc != 0.75))
forbidden_overhangs = list_overhangs(filters=[invalid_overhang])

selector = OverhangsSelector(forbidden_overhangs=forbidden_overhangs,
	                         differences=1, time_limit=2)
indices = selector.cut_sequence_into_similar_lengths(sequence, nsegments=50,
                                                     max_radius=20,
                                                     include_extremities=False)
print "solution", indices
import numpy
diffs = numpy.diff([0] + indices + [len(sequence)])
print len(diffs), diffs.min(), diffs.max(), numpy.mean(diffs), numpy.std(diffs)