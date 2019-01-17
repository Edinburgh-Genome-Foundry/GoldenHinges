Golden Hinges
==============
Golden Hinges (full documentation `here <https://edinburgh-genome-foundry.github.io/GoldenHinges/>`_) is a Python library to find sets of overhangs (also called junctions, or protusions) for multipart DNA assembly such as Golden Gate assembly.

Given a set of constraints (GC content bounds, differences between overhangs,
mandatory and forbidden overhangs) Golden Hinges enables to find:

- Maximal sets of valid and inter-compatible overhangs.
- Sequence decompositions (i.e. position of cuts) which produce valid and
  inter-compatible overhangs, for type-2S DNA assembly.
- Sequence mutations (subject to constraints) which enable the sequence
  decomposition, in extreme cases where the original sequence does not allow
  for such decomposition.

Infos
-----

**PIP installation:**

.. code:: bash

  pip install goldenhinges

**Web documentation:**

`<https://edinburgh-genome-foundry.github.io/GoldenHinges/>`_

**Github Page (with examples)**

`<https://github.com/Edinburgh-Genome-Foundry/GoldenHinges>`_

**Live demo**

`<http://cuba.genomefoundry.org/design_overhangs>`_

**License:** MIT, Copyright Edinburgh Genome Foundry

More biology software
-----------------------

.. image:: https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/Edinburgh-Genome-Foundry.github.io/master/static/imgs/logos/egf-codon-horizontal.png
  :target: https://edinburgh-genome-foundry.github.io/

Golden Hinges is part of the `EGF Codons <https://edinburgh-genome-foundry.github.io/>`_ synthetic biology software suite for DNA design, manufacturing and validation.
