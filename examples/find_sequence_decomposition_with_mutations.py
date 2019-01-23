"""Example of cutting problem defined by a record.

This example has a long
"""

import os
from Bio import SeqIO
from goldenhinges import OverhangsSelector
from goldenhinges.reports import write_report_for_cutting_solution


genbank_file = os.path.join('data', 'sequence_with_constraints.gb')
record = SeqIO.read(genbank_file, 'genbank')
selector = OverhangsSelector(gc_min=0.25, gc_max=0.75, differences=2)
solution = selector.cut_sequence(record, allow_edits=True,
                                 equal_segments=20,
                                 include_extremities=True)

print ("Writing the report...")
write_report_for_cutting_solution(
    solution=solution, target=os.path.join('results', 'with_mutations'),
    sequence=record, left_flank='CGTCTCA', right_flank='TGAGACG',
)
print ("Done! See the report in results/with_mutations/")
