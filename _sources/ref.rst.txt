Golden Hinges Reference Manual
==============================

Golden Hinges has one main class called ``OverhangsSelector``. One instance of
this class determines which Overhangs are accepted and what "intercompatibility"
of overhangs means. The different methods the class then enable to solve
different scenarii.

OverhangsSelector
------------------

.. autoclass:: goldenhinges.OverhangsSelector
   :members:

Clique methods
--------------

Methods based on graph cliques enable to quickly find large albeit non necessary optimal sets of compatible overhangs

.. automodule:: goldenhinges.clique_methods
   :members:

Biotools
---------

.. automodule:: goldenhinges.biotools
   :members:
