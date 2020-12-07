=======================================================================================
Welcome to cn_hyperarr's documentation!
=======================================================================================

cn_hyperarr is a package for the computation of congruence normality of
three-dimensional hyperplane arrangements.

A dictionary containing the normal vectors to all currently known simplicial 
arrangements of rank 3 is available.

There is a module for vector configurations which could be useful independently.
A vector configuration can be seen as the set of normals to a hyperplane
arrangement. A simplicial hyperplane arrangement has a lattice of regions 
associated to each chamber. To check whether this lattice is obtainable through
a sequence of doublings of convex sets, i.e. is congruence normal, we use the 
theory developed in [1]. For each choice of chamber, there is an associated 
acyclic vector configuration. We compute the shard covectors
and the forcing oriented graph on the shard covectors. If this oriented 
graph is acyclic, then the arrangement is congruence normal with 
respect to the chosen chamber; this property is tested in the method 
:func:`~vector_classes.VectorConfiguration.is_congruence_normal`.

In the computations module, we test congruence normality for each region 
of a hyperplane arrangement.

Finally there is a module for constructing the three infinite families of 
rank-3 simplicial arrangements.

REFERENCES:

    - [1] Michael Cuntz, Sophia Elia, and Jean-Philippe Labbé. Congruence normality of simplicial hyperplane arrangements via oriented matroids, 2020. arXiv:2009.14152.

AUTHORS:

- Jean-Philippe Labbé (2020): Initial version
- Sophia Elia (2020): Initial version

Congruence Normality for Hyperplane Arrangements
=======================================================================================

.. toctree::
   :maxdepth: 2

   installation
   computations
   vector_conf
   inf_fam

.. automodule:: cn_hyperarr
   :members:
   :undoc-members:
   :show-inheritance:

Indices and Tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
