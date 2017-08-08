"""Example of cutting problem defined by a record."""

import numpy as np
from Bio import SeqIOfrom goldenhinges import OverhangsSelector
from dnachisel.biotools import (annotate_record, random_dna_sequence,
                                sequence_to_biopython_record)
from dnachisel.reports import SpecAnnotationsTranslator


# DESIGN A SEQUENCE

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

cut_region_size = 600
zones = [
    (x+ int(200*np.sin(x)),
     x + cut_region_size + int(200*np.sin(x) - 50*np.cos(x)),
     0)
    for x in range(50, len(seq), 1030)[1:]
]
for zone in zones:
    annotate_record(record, location=zone, label="!cut")

# MAKE A PLOT AND A GENBANK FROM THE PROBLEM
translator = SpecAnnotationsTranslator()
grecord = translator.translate_record(record)
ax, _ = grecord.plot(figure_width=7.2, annotate_inline=True)
ax.figure.savefig(os.path.join('results', 'problem.svg'), bbox_inches="tight")
SeqIO.write(record, os.path.join('results', 'problem.gb'), "genbank")

# SOLVE PROBLEM

selector = OverhangsSelector(gc_min=0.25, gc_max=0.75, differences=2)
solution = selector.cut_sequence(record, allow_edits=True,
                                 include_extremities=True)
print ("solution", solution)
if solution is not None:
    print ([o['score'] for o in solution])
