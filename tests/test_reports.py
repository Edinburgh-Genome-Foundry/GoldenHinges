import matplotlib
matplotlib.use("Agg")

from goldenhinges.reports import write_report_for_cutting_solution
from dnachisel import (random_dna_sequence, sequence_to_biopython_record)
from goldenhinges import OverhangsSelector

def test_sequence_cut_with_report():
    seq = random_dna_sequence(8000, seed=123)
    record = sequence_to_biopython_record(seq)
    selector = OverhangsSelector(differences=1, time_limit=2)
    solution = selector.cut_sequence(record, equal_segments=12,
                                     include_extremities=True)
    zip_data = write_report_for_cutting_solution(
        solution=solution, target="@memory", sequence=record,
        left_flank="CGTCTCA",
        right_flank="TGAGACG"
    )
