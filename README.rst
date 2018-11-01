.. raw:: html

    <p align="center">
    <img alt="Golden Hinges Logo" title="Golden Hinges Logo" src="https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/GoldenHinges/master/docs/_static/images/title.png"
   :alt: [logo]" width="600">
    <br /><br />
    </p>

.. image:: https://travis-ci.org/Edinburgh-Genome-Foundry/GoldenHinges.svg?branch=master
   :target: https://travis-ci.org/Edinburgh-Genome-Foundry/GoldenHinges
   :alt: Travis CI build status

Golden Hinges (full documentation `here <https://edinburgh-genome-foundry.github.io/GoldenHinges/>`_) is a Python library to find sets
of overhangs (also called junctions, or protusions) for multipart DNA assembly
such as Golden Gate assembly.

Given a set of constraints (GC content bounds, differences between overhangs,
mandatory and forbidden overhangs) Golden Hinges enables to find:

- Maximal sets of valid and inter-compatible overhangs.
- Sequence decompositions (i.e. position of cuts) which produce valid and
  inter-compatible overhangs, for type-2S DNA assembly.
- Sequence mutations (subject to constraints) which enable the sequence
  decomposition, in exterme cases where the original sequence does not allow
  for such decomposition.

You can see Golden Hinges in action in this
`web demo <http://cuba.genomefoundry.org/design-overhangs>`_:

Examples of use
----------------

Finding maximal overhangs sets
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Let us compute a collection of overhangs, as large as possible, where

- All overhangs have 25-75 GC%
- There is a 2-basepair difference between any two overhangs (and their reverse-complement)
- The overhangs ``ATGC`` and ``CCGA`` are forbidden

Here is the code

.. code:: python

  from goldenhinges import OverhangsSelector
  selector = OverhangsSelector(
      gc_min=0.25,
      gc_max=0.5,
      differences=2,
      forbidden_overhangs=['ATGC', 'CCGA']
  )
  overhangs = selector.generate_overhangs_set()
  print (overhangs)

Result:

.. code:: python

    >>> ['AACG', 'CAAG', 'ACAC', 'TGAC', 'ACGA', 'AGGT',
         'TGTG', 'ATCC', 'AAGC', 'AGTC', 'TCTC', 'TAGG',
         'AGCA', 'GTAG', 'TGGA', 'ACTG', 'GAAC', 'TCAG',
         'ATGG', 'TTGC', 'TTCG', 'GATG', 'AGAG', 'TACC']

In some cases this may take some time to complete, as the algorithm slowly builds
collections of increasing sizes. An alternative algorithm consisting in finding
random maximal sets of compatible overhangs is much faster, but gives suboptimal
solutions:

.. code:: python

    overhangs = selector.generate_overhangs_set(n_cliques=5000)

Result:

.. code:: python

    >>> ['CAAA', 'GTAA', 'ATTC', 'AATG', 'ACAT', 'ATCA',
         'AGAG', 'GCTT', 'AGTT', 'TCGT', 'CTGA', 'TGGA',
         'TAGG', 'GGTA', 'GACA']

The two approaches can be combined to first find an approximate solution, then
attempt to find larger sets:

.. code:: python

    test_overhangs = selector.generate_overhangs_set(n_cliques=5000)
    overhangs = selector.generate_overhangs_set(start_at=len(test_overhangs))


Using experimental annealing data from Potapov 2018
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`This study by Potapov *et al.* <https://www.biorxiv.org/content/early/2018/05/15/322297>`_
provides insightful data on overhangs annealing, in particular which overhangs
have weak general annealing power, and which pairs of overhangs have significant
"cross-talk". You can use the data in this paper via the Python
`tatapov <https://github.com/Edinburgh-Genome-Foundry/tatapov>`_ library
to identify which overhangs or overhang pairs you want the GoldenHinges
``OverhangSelector`` to exclude:


.. code:: python

    import tatapov
    from goldenhinges import OverhangsSelector

    annealing_data = tatapov.annealing_data['37C']['01h']

    self_annealings = tatapov.relative_self_annealings(annealing_data)
    weak_self_annealing_overhangs = [
        overhang
        for overhang, self_annealing in self_annealings.items()
        if self_annealing < 0.05
    ]

    cross_annealings = tatapov.cross_annealings(annealing_data)
    high_cross_annealing_pairs = [
        overhang_pair
        for overhang_pair, cross_annealing in cross_annealings.items()
        if cross_annealing > 0.005
    ]

    selector = OverhangsSelector(
        forbidden_overhangs=weak_self_annealing_overhangs,
        forbidden_pairs=high_cross_annealing_pairs
    )


Finding a sequence decomposition
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


In this example, we find where to cut a 50-kilobasepair sequence to create
assemblable fragments with 4-basepair overhangs. We indicate that:

- There should be 50 fragments, with a minimum of variance in their sizes.
- The fragments overhangs should have 25-75 GC% with a 1-basepair difference
  between any two overhangs (and their reverse-complement). They should also be
  compatible with the 4-basepair extremities of the sequence.

.. code:: python

    from Bio import SeqIO
    from goldenhinges import OverhangsSelector

    sequence = SeqIO.read
    selector = OverhangsSelector(gc_min=0.25, gc_max=0.75, differences=1)
    solution = selector.cut_sequence(
        sequence, equal_segments=50, max_radius=20,
        include_extremities=True
    )

This returns a list of dictionnaries, each representing one overhang with
properties ``o['location']`` (coordinate of the overhang in the sequence)
and ``o['sequence']`` (sequence of the overhang).

This solution can be turned into a full report featuring all sequences to order
(with restriction sites added on the left and right flanks), and a graphic of
the overhang's positions, using the following function:


.. code:: python

    from goldenhinges.reports import write_report_for_cutting_solution

    write_report_for_cutting_solution(
        solution, 'full_report.zip', sequence,
        left_flank='CGTCTCA', right_flank='TGAGACG',
        display_positions=False
    )

Sequence mutation and decomposition from a Genbank file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If the input sequence is a Genbank record (or a Biopython record) has locations
annotated vy features feature labeled ``!cut``, GoldenHinges will attempt to
find a decomposition with exactly one cut in each of these locations (favoring
cuts located near the middle of each region).

GoldenHinges also allows to modify the sequence to enable some decomposition.
Note that solutions involving base changes are penalized and solutions involving
the original solution will always be prefered, so no base change will be
suggested unless strictly necessary.

If the input record has `DnaChisel <https://github.com/Edinburgh-Genome-Foundry/DnaChisel>`_
annotations such as ``@AvoidChanges`` or ``@EnforceTranslation``, these will be
enforced to forbid some mutations.

Here is an example of such a record:

.. image:: https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/GoldenHinges/master/examples/data/sequence_with_constraints.png
   :alt: [sequence with constraints]
   :align: center
   :width: 672px

And here is the code to optimize and decompose it:

.. code:: python

    record = SeqIO.read(genbank_file, 'genbank')
    selector = OverhangsSelector(gc_min=0.25, gc_max=0.75,
                                 differences=2)
    solution = selector.cut_sequence(record, allow_edits=True,
                                     include_extremities=True)

Installation
--------------

Install Numberjack's dependencies first:

.. code:: python

    sudo apt install python-dev swig libxml2-dev zlib1g-dev libgmp-dev

If you have PIP installed, just type in a terminal:

.. code:: python

    (sudo) pip install goldenhinges

Golden Hinges can be installed by unzipping the source code in one directory and
using this command:

.. code:: python

    sudo python setup.py install



Contribute !
--------------

Golden Hinges is an open-source software originally written at the
`Edinburgh Genome Foundry <http://edinburgh-genome-foundry.github.io/home.html>`_
by `Zulko <https://github.com/Zulko>`_ and
`released on Github <https://github.com/Edinburgh-Genome-Foundry/GoldenHinges>`_
under the MIT licence. Everyone is welcome to contribute !
