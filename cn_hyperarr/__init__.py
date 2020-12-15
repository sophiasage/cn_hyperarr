"""
Congruence Normality Package
===============================================================

This is a package to compute congruence normality of simplicial
hyperplane arrangements of rank three and make the database of 
known arrangements of rank three usable/available.

With this package, you can compute the following:

    - vector configurations
    - covectors and shard covectors
    - forcing oriented graph on shards
    - congruence normality of the poset of regions
    - hyperplane arrangements from vector configurations
    - the vector configurations of the three infinite families of simplicial 
      rank-three arrangments

With this package, you can load:

    - normals of all known simplicial arrangements of rank-three
    - invariants of all known simplicial arrangements of rank-three
    - wiring diagrams of all known simplicial arrangements of rank-three
"""
from __future__ import absolute_import
from .main import *
del absolute_import
