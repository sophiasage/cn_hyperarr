=======================================================================================
Welcome to cn_hyperarr's documentation!
=======================================================================================

cn_hyperarr is a package for the computation of congruence normality of
three-dimensional hyperplane arrangements.

A :ref:`database<Database>` containing a dictionary with the normal vectors to 
all known simplicial arrangements of rank 3 is available. 
There is a module for :ref:`vector configurations<vector configurations>`,
which could be used independently.
A vector configuration can be seen as the set of normals to a hyperplane
arrangement. A simplicial hyperplane arrangement has a lattice of regions
associated to each chamber. To check whether this lattice is obtainable through
a sequence of doublings of convex sets, i.e. is congruence normal, we use the
theory developed in [CEL]_. For each choice of chamber, there is an associated
acyclic vector configuration. 
Using this ayclic configuration, we compute the shard covectors
and the forcing oriented graph on the shard covectors. 
The forcing oriented graph is acyclic if and only if the arrangement is 
congruence normal with respect to the chosen chamber.
This property is tested in the method
:func:`~vector_classes.VectorConfiguration.is_congruence_normal`.
In the :ref:`computations module<computations module>`, we test congruence 
normality for each region of a hyperplane arrangement.
Finally there is a module for constructing the 
:ref:`three infinite families<inf_fam>` of rank-3 simplicial arrangements.

To summarize, with this package, you can compute the following:

    - vector configurations
    - covectors and shard covectors
    - forcing oriented graph on shards
    - congruence normality of the poset of regions
    - hyperplane arrangements from vector configurations
    - the vector configurations of the three infinite families of simplicial
      rank-three arrangements [Gru]_

With this package, you can load:

    - normals of all known simplicial arrangements of rank-three
    - invariants of all known simplicial arrangements of rank-three
    - wiring diagrams of all known simplicial arrangements of rank-three

EXAMPLES:

Here is an example of the arrangement `A(10,60)_3` that is sometimes
congruence normal.
First we load the normals of the arrangement from the database.
The entries of the database are labeled in the same way as in [CEL]_::

    sage: from cn_hyperarr import *
    sage: s_normals = db_normals_CEL[(10, 60, 3)]

Now we make the normals into a vector configuration::

    sage: s_vc = VectorConfiguration([vector(x) for x in s_normals])
     
To test congruence normality, use the :func:`~main.RegionsCongruenceNormal` method::

    sage: somet_check = RegionsCongruenceNormality(s_vc)
    sage: somet_vals_list = list(somet_check.values())
    sage: [somet_vals_list.count(True), somet_vals_list.count(False)]
    [40,20]

REFERENCES:

    .. [CEL] Michael Cuntz, Sophia Elia, and Jean-Philippe Labbé. Congruence normality of simplicial hyperplane arrangements via oriented matroids, 2020. arXiv:2009.14152.

    .. [Gru] Branko Grunbaum. A catalogue of simplicial arrangements in the real projective plane, 2009. Ars Math. Contemp. 2, no. 1, 1-25.

AUTHORS:

    - Michael Cuntz (2020): Initial version
    - Sophia Elia (2020): Initial version
    - Jean-Philippe Labbé (2020): Initial version


Congruence Normality for Hyperplane Arrangements
=======================================================================================

.. toctree::
   :maxdepth: 2

   installation
   computations
   inf_fam
   database
   vector_conf

.. automodule:: cn_hyperarr
   :members:
   :undoc-members:
   :show-inheritance:

Indices and Tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
