from dna_features_viewer import BiopythonTranslator
import numpy as np
from copy import deepcopy
from Bio import SeqIO
import flametree
import matplotlib.pyplot as plt
from geneblocks import DiffBlocks

from .biotools import (annotate_record, sequence_to_biopython_record,
                       sequences_differences_segments, crop_record)

def new_sequence_from_cutting_solution(solution, sequence):
    """Return a new sequence will all mutations in the cutting solution.
    
    Input sequence and returned sequences are ATGC strings.
    """
    new_sequence = np.array(list(sequence))
    for o in solution:
        if o.get('n_mutations', 0) == 0:
            continue

        start, mutated_seq = o['mutated_region']
        end = start + len(mutated_seq)
        new_sequence[start:end] = np.array(list(mutated_seq))
    return "".join(new_sequence)

def write_report_for_cutting_solution(solution, target, sequence,
                                      left_flank='', right_flank='',
                                      display_positions=False):
    """Write a complete report for Type IIS arbitrary sequence assembly.

    Parameters
    -----------

    solution
      The solution returned by an OverhangsSelector's ``cut_sequence`` method.

    target
      Either a path to a folder, a zip, or "@memory" to return raw ZIP file
      data instead of writing files. If ``target`` poitns to an exsisting
      folder/zip, it will be completely overwritten.

    sequence
      Sequence to be cut (can be a record)

    left_flank
      Left flank to be added to every fragment

    right_flank
      Right flank to be added to every fragment

    display_positions
      If True, the exact coordinate of each cut will be reported in the plot.
    """

    root = flametree.file_tree(target, replace=True)
    if isinstance(left_flank, str):
        left_flank = sequence_to_biopython_record(left_flank)
        annotate_record(left_flank, label='left_flank')

    if isinstance(right_flank, str):
        right_flank = sequence_to_biopython_record(right_flank)
        annotate_record(right_flank, label='right_flank')

    if hasattr(sequence, "seq"):
        record = sequence
        sequence = str(record.seq)
    else:
        record = sequence_to_biopython_record(sequence)
    # COMPUTE THE EDITED SEQUENCE (MAY BE EQUAL TO ORIGINAL IF NO EDITS)
    new_sequence = new_sequence_from_cutting_solution(solution, sequence)
    edited_segments = sequences_differences_segments(sequence, new_sequence)
    blocks = DiffBlocks.from_sequences(sequence, new_sequence)
    if hasattr(sequence, 'features'):
        ax, _ = blocks.plot(separate_axes=True)
    else:
        ax = blocks.plot(separate_axes=False)
    ax.set_title("Edits in new sequence vs. original")
    ax.figure.savefig(root._file('edits.pdf').open('wb'), format='pdf',
                      bbox_inches='tight')
    plt.close(ax.figure)
    # PLOT SUMMARY FIGURE

    plot_record = sequence_to_biopython_record(sequence)
    display_positions = False
    for o in solution:
        start, end = o['location'], o['location'] + len(o['sequence'])
        label = ("%s\n(%d)" % (o['sequence'], o['location'])
                 if display_positions else o['sequence'])
        annotate_record(plot_record, (start, end, 0), label=label)

    translator = BiopythonTranslator()
    gr = translator.translate_record(plot_record)
    ax, _ = gr.plot(with_ruler=False, figure_width=max(8, len(solution) / 2))
    ax.set_title("Selected overhangs", loc="left",
                 fontdict=dict(weight='bold', fontsize=13))
    # ax.figure.set_size_inches((max(8, 0.7*len(o)), 2))
    ax.set_ylim(top=ax.get_ylim()[1] + 2)

    xx = [x for (a, b) in edited_segments for x in range(a, b)]
    ax.plot(xx, [0 for x in xx], marker='o', c='r', lw=0,
            label='sequence edits')
    L = len(sequence)
    ax.set_xlim(-.1 * L, 1.1 * L)
    ax.legend(loc=2, fontsize=12)
    locs = sorted([o['location'] for o in solution])
    diffs = np.diff(locs)
    text = "Segment size: %d +/- %d bp. (mean +/- 1std)" % (diffs.mean(),
                                                            diffs.std())
    ax.text(L / 2, -1, text, horizontalalignment="center",
            verticalalignment="top", fontsize=14)
    ax.figure.savefig(root._file('summary_plot.pdf').open('wb'), format='pdf',
                      bbox_inches='tight')
    plt.close(ax.figure)

    #  WRITE GENBANK RECORD OF FINAL SEQUENCE

    report_record = deepcopy(record)
    report_record.seq = sequence_to_biopython_record(new_sequence).seq
    for (start, end) in edited_segments:
        annotate_record(report_record, (int(start), int(end), 0),
                        label="!edited")
    for o in solution:
        start = int(o['location'])
        end = int(o['location'] + len(o['sequence']))
        annotate_record(report_record, (start, end, 0), label='overhang')
        annotate_record(report_record, (start, end, 0), label='@DoNotModify')
    SeqIO.write(report_record, root._file('final_sequence.gb'), 'genbank')

    #  WRITE GENBANK RECORDS OF ALL FRAGMENTS
    sequences = []
    fragments_records_dir = root._dir("fragments_records")
    for i, (o1, o2) in enumerate(zip(solution, solution[1:])):
        seqname = "fragment_%02d" % (i + 1)
        start, end = o1['location'], o2['location'] + len(o2['sequence'])
        fragment = crop_record(report_record, start, end)
        seqrecord = left_flank + fragment + right_flank
        SeqIO.write(seqrecord, fragments_records_dir._file(seqname + ".gb"),
                    'genbank')
        sequences.append(";".join([seqname, str(seqrecord.seq)]))
    root._file("fragments_sequences.csv").write("\n".join(sequences))

    root._file('overhangs_list.csv').write(", ".join([
        o['sequence'] for o in solution
    ]))

    return root._close()
