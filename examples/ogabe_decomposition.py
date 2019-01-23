import os
import numpy as np
from goldenhinges import (OverhangsSelector, list_overhangs, gc_content,
                          write_report_for_cutting_solution)


with open(os.path.join("data", "phage_sequence.txt"), "r") as f:
    sequence = f.read()

def invalid_overhang(overhang):
    gc = gc_content(overhang)
    three_gc = max([gc_content(overhang[:-1]), gc_content(overhang[1:])]) == 1
    return (gc != 0.5) and (three_gc or (gc != 0.75))


forbidden_overhangs = list_overhangs(filters=[invalid_overhang])

selector = OverhangsSelector(forbidden_overhangs=forbidden_overhangs,
                             differences=1, time_limit=2)
solution = selector.cut_sequence(
    sequence, equal_segments=50, max_radius=20, include_extremities=False)

indices = [o['location'] for o in solution]
diffs = np.diff([0] + indices + [len(sequence)])
print ("""
Cut indices: %s

Number of segments: %d
Segments length: %d - %d
Segments mean and std: %.02f +/- %.02f
""" % (indices, len(diffs), diffs.min(), diffs.max(), diffs.mean(),
       np.std(diffs)))

print ("Writing a report...")
write_report_for_cutting_solution(
    solution=solution, sequence=sequence, left_flank='', right_flank='',
    target=os.path.join('results', 'ogabe_decomposition'))
print ("Done! See report in results/ogabe_decomposition/")