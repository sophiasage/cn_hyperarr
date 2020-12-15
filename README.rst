===================================================
Congruence Normality for Hyperplane Arrangements
===================================================
.. image:: https://travis-ci.org/sophiasage/cn_hyperarr.svg?branch=master
    :target: https://travis-ci.org/sophiasage/cn_hyperarr

This package is a `SageMath <http://www.sagemath.org>`_ package for computing congruence normality of rank-three, simplicial hyperplane arrangements.

This package includes a database of known rank-three simplicial hyperplane
arrangements. It also includes modules for creating vector configurations and 
the three infinite families of simplicial rank-three arrangements.
A vector configuration can be seen as the set of normals to a hyperplane
arrangement. A simplicial hyperplane arrangement has a lattice of regions 
associated to each chamber. This lattice is congruence normal if it is 
obtainable through a sequence of doublings of convex sets. 
A hyperplane arrangement can be always or sometimes or never congruence normal,
depending on whether its lattices of regions are congruence normal.

Here are examples of arrangements that are always, sometimes, and never
congruence normal. 
First we load the normals of the three arrangements from the database. 
The entries of the database are labeled in the same way as in [CEL]_ ::

    sage: from cn_hyperarr import *
    sage: always_normals = db_normals_CEL[(6,24,1)] 
    sage: somet_normals = db_normals_CEL[(10,60,3)]
    sage: never_normals = db_normals_CEL[(22,288,1)]

Now we make them into vector configurations::

    sage: always_vc = VectorConfiguration((vector(x) for x in always_normals)) 
    sage: somet_vc = VectorConfiguration((vector(x) for x in somet_normals)) 
    sage: never_vc = VectorConfiguration((vector(x) for x in never_normals))
     
To test congruence normality, use the :func:`~main.RegionsCongruenceNormal`
method::

    sage: always_check = RegionsCongruenceNormality(always_vc)
    sage: always_vals_list = list(always_check.values())
    sage: [always_vals_list.count(True), always_vals_list.count(False)]
    [24,0]
    sage: somet_check = RegionsCongruenceNormality(somet_vc)
    sage: somet_vals_list = list(somet_check.values())
    sage: [somet_vals_list.count(True), somet_vals_list.count(False)]
    [40,20]
    sage: never_check = RegionsCongruenceNormality(never_vc)
    sage: never_vals_list = list(never_check.values())
    sage: [never_vals_list.count(True), never_vals_list.count(False)]
    [0,288]

The full documentation for the package can be found at https://sophiasage.github.io/cn_hyperarr/doc/html/


Installation
------------

Local install from source
^^^^^^^^^^^^^^^^^^^^^^^^^

Download the source from the git repository::

    $ git clone https://github.com/sophiasage/cn_hyperarr.git

Change to the root directory and run::

    $ sage -pip install --upgrade --no-index -v .

For convenience this package contains a [makefile](makefile) with this
and other often used commands. Should you wish too, you can use the
shorthand::

    $ make install


Usage
-----


Once the package is installed, you can use it in Sage. To do so you have to import it with::

    sage: from cn_hyperarr import *
    
You can test whether a hyperplane arrangement is congruence normal as follows::

    sage: vc = VectorConfiguration([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[0,1,1],[1,1,1]]);
    sage: check_result = RegionsCongruenceNormality(vc);
    sage: check_result.values()
    dict_values([True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True])

To load the list of normals in sage::

    sage: from cn_hyperarr.database import *
    sage: db_normals_CEL[(7,32,1)]
    ((0, 0, 1), (0, 1, 0), (0, 1, 1), (1, 0, 0), (1, 0, 1), (1, 1, 0), (1, 1, 1))
