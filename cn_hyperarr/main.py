# -*- coding: utf-8 -*-
r"""
Module for the computations related to congruence normality/uniformity of
posets of regions of hyperplane arrangements.

This module has methods to test simplicial hyperplane arrangements of rank 3
for congruence normality.  A hyperplane arrangement is congruence normal if
each of its posets of regions is obtainable from the one element lattice
through a finite sequence of doublings of convex sets.  Equivalently, one can
create a directed graph on shards of the arrangement for each base region and
check if it is acyclic, see the module 
:class:`~vector_classes.VectorConfiguration` and its method
:func:`~vector_classes.VectorConfiguration.is_congruence_normal`.  The method
:func:`~main.vectorconf_to_hyperplane_arrangement` creates a hyperplane
arrangement from a three-dimensional vector configuration.  Every region of
each arrangement is then tested to determine if the corresponding poset is
congruence normal in the :func:`~main.RegionsCongruenceNormal` method.

EXAMPLES:

This vector configuration is congruence normal for any choice of base region in
the corresponding hyperplane arrangement. We load it from the database::

    sage: from cn_hyperarr.main import *
    sage: always_normals = db_normals_CEL[(6,24,1)]
    sage: always_normals
    ((0, 0, 1), (0, 1, 0), (0, 1, 1), (1, 0, 0), (1, 0, 1), (1, 1, 1))
    sage: vc = VectorConfiguration(list(vector(x) for x in always_normals));
    sage: h = vectorconf_to_hyperplane_arrangement(vc);
    sage: h.n_regions()
    24
    sage: check_result = RegionsCongruenceNormality(vc);
    sage: check_result.values()
    dict_values([True, True, True, True, True, True, True, True, True, True, True, 
    True, True, True, True, True, True, True, True, True, True, True, True, True])
    sage: len(check_result.values())
    24

This vector configuration of ten vectors is sometimes congruence normal.
It is the smallest simplicial arrangement of rank three that is not always
congruence normal and is referred to as A(10,60)_3 in [CEL]_. We also
load it from the database::

    sage: from cn_hyperarr.main import *
    sage: sometimes_normals = db_normals_CEL[(10,60,3)]
    sage: sometimes_normals
    ((2*tau + 1, 2*tau, tau),
     (2*tau + 2, 2*tau + 1, tau + 1),
     (1, 1, 1),
     (tau + 1, tau + 1, tau),
     (2*tau, 2*tau, tau),
     (tau + 1, tau + 1, 1),
     (1, 1, 0),
     (0, 1, 0),
     (1, 0, 0),
     (tau + 1, tau, tau))
    sage: ncn_conf = VectorConfiguration(list(vector(x) for x in sometimes_normals));
    sage: check_result = RegionsCongruenceNormality(ncn_conf)
    sage: check_result.values()
    dict_values([False, False, True, True, True, True, False, False, True, True, False, True, True, False, True, False, True, True, True, True, True, True, False, False, False, True, True, True, True, True, False, False, True, True, True, True, False, False, True, True, False, True, True, False, True, False, True, True, True, True, True, True, False, False, False, True, True, True, True, True])

AUTHORS:

- Jean-Philippe Labbé (2020): Initial version
- Sophia Elia (2020): Initial version
"""

##############################################################################
#     Copyright (C) 2020 Jean-Philippe Labbe <labbe at math.fu-berlin.de>
#                   2020 Sophia Elia         <sophiae56 at math.fu-berlin.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
##############################################################################

import time
import multiprocessing
from .database import *
from .vector_classes import *
from .infinite_families import *
from sage.rings.integer_ring import ZZ
from sage.geometry.hyperplane_arrangement.arrangement import HyperplaneArrangements

##############################################################################
# Helper functions
##############################################################################


def vectorconf_to_hyperplane_arrangement(vector_conf, backend=None):
    r"""
    Return the hyperplane arrangement associated to the vector configuration
    ``vector_conf``.

    INPUT:

    - ``vector_conf`` -- a vector configuration

    - ``backend`` -- string (default = ``None``). The backend to use.

    OUTPUT:

    A hyperplane arrangement

    EXAMPLES:

    The arrangement A(10,60)`_3` with 10 hyperplanes is the smallest rank-three 
    simplicial arrangement that is not congruence normal::

        sage: from cn_hyperarr.main import *
        sage: tau = AA((1+sqrt(5))/2)
        sage: ncn = [[2*tau + 1, 2*tau, tau], [2*tau + 2, 2*tau + 1, tau + 1], 
        ....:        [1, 1, 1], [tau + 1, tau + 1, tau], [2*tau, 2*tau, tau], 
        ....:        [tau + 1, tau + 1, 1], [1, 1, 0], [0, 1, 0], [1, 0, 0], 
        ....:        [tau + 1, tau, tau]]
        sage: ncn_conf = VectorConfiguration(ncn);
        sage: ncn_arr = vectorconf_to_hyperplane_arrangement(ncn_conf); ncn_arr
        Arrangement of 10 hyperplanes of dimension 3 and rank 3

    Using normaliz as backend::

        sage: ncn_conf_norm = VectorConfiguration(ncn, 'normaliz')   # optional - pynormaliz
        sage: ncn_conf_norm.backend()                                # optional - pynormaliz
        'normaliz'
    """
    if not isinstance(vector_conf, VectorConfiguration):
        vector_conf = VectorConfiguration(vector_conf, backend=backend)
    if vector_conf.base_ring() == ZZ:
        H = HyperplaneArrangements(QQ, names='xyz')
    else:
        H = HyperplaneArrangements(vector_conf.base_ring(), names='xyz')
    x, y, z = H.gens()
    A = H(backend=vector_conf.backend())
    for v in vector_conf:
        a, b, c = v
        A = A.add_hyperplane(a*x + b*y + c*z)
    return A


def wrapper_forcing_acyclic(vectorconf):
    r"""
    Return whether the vector configuration coming from a specific choice
    of base region gives rise to a congruence normal poset of regions.

    INPUT:

    - ``vectorconf`` -- a vector configuration

    OUTPUT:

    A tuple. The first entry says if the forcing_orient_graph is acyclic. The
    second says the number of vertices.

    EXAMPLES:

    The arrangement A(10,60)_3 with 10 hyperplanes is not congruence 
    normal for the implicit choice of base region::

        sage: from cn_hyperarr.main import *
        sage: tau = AA((1+sqrt(5))/2)
        sage: ncn = [[2*tau+1,2*tau,tau],[2*tau+2,2*tau+1,tau+1],[1,1,1],[tau+1,tau+1,tau],[2*tau,2*tau,tau],[tau+1,tau+1,1],[1,1,0],[0,1,0],[1,0,0],[tau+1,tau,tau]]
        sage: ncn_conf = VectorConfiguration(ncn);
        sage: wrapper_forcing_acyclic(ncn_conf);
        (False, 29)

    A congruence normal arrangement, with an acyclic forcing_oriented_graph::

        sage: vc = VectorConfiguration([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[0,1,1],[1,1,1]])
        sage: wrapper_forcing_acyclic(vc)
        (True, 11)
    """
    forcing = vectorconf.forcing_oriented_graph()
    return forcing.is_directed_acyclic(), forcing.num_verts()

##############################################################################
# Verification functions
##############################################################################


def RegionsCongruenceNormality(vector_conf, backend=None, verbose=False, nb_proc=4):
    r"""
    Return whether regions in the hyperplane arrangement generated by the
    vector configuration ``vector_conf`` lead to congruence normal lattice of
    regions.

    INPUT:

    - ``vectorconf`` -- a vector configuration

    - ``backend`` -- string (default = ``None``).

    - ``verbose`` -- (default = False). If ``True`` return number of
      hyperplanes, normal vectors, number of regions.

    - ``nb_proc`` -- (default = ``4``). The number of processors to use to
      verify. 

    OUTPUT:

    A dictionary whose keys are tuple (vector,vc) where vector is a
    representative vector of a region of the hyperplane arrangement and vc is
    the corresponding acyclic vector configuration which is a reorientation of
    ``vector_conf``.

    EXAMPLES:

    An example of an arrangement that is always congruence normal::

        sage: from cn_hyperarr.main import *
        sage: vc = VectorConfiguration([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[0,1,1],[1,1,1]])
        sage: h = vectorconf_to_hyperplane_arrangement(vc);
        sage: h.n_regions()
        24
        sage: check_result = RegionsCongruenceNormality(vc);
        sage: check_result.values()
        dict_values([True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True])
        sage: len(check_result.values())
        24

    This vector configuration of ten vectors is sometimes congruence normal.
    It is the smallest simplicial arrangement of rank three that is not always
    congruence normal and is referred to as A(10,60)_3 in [CEL]_ ::

        sage: tau = AA((1+sqrt(5))/2)
        sage: ncn = [[2*tau+1,2*tau,tau],[2*tau+2,2*tau+1,tau+1],[1,1,1],[tau+1,tau+1,tau],[2*tau,2*tau,tau],[tau+1,tau+1,1],[1,1,0],[0,1,0],[1,0,0],[tau+1,tau,tau]]
        sage: ncn_conf = VectorConfiguration(ncn);
        sage: check = RegionsCongruenceNormality(ncn_conf)
        sage: vals_list = list(check.values())
        sage: [vals_list.count(True), vals_list.count(False)]
        [40, 20]
    """
    A = vectorconf_to_hyperplane_arrangement(vector_conf, backend)

    if verbose:
        print("The arrangement has {} hyperplanes.".format(len(A)))
        print("The normal vectors are:\n {}".format(vector_conf))
        print("Getting regions", flush=True)
    R = A.regions()
    if verbose:
        print("Found {} regions".format(len(R)), flush=True)
    repr_vector = [c.representative_point() for c in R]

    the_vector_configs = []
    for rc in repr_vector:
        list_of_normals = []
        for h in A.hyperplanes():
            normal = h.A()
            if normal.dot_product(rc) > 0:
                list_of_normals += [tuple(list(-normal))]
            else:
                list_of_normals += [tuple(list(normal))]
        rc.set_immutable()
        # The following line is A REAL PROBLEM
        # the_vector_configs += [VectorConfiguration(list_of_normals, backend='a_backend')]
        # Do we want normaliz here? It seems not necessary, but might give some
        # speed within the class VectorConf ... let's leave it as it for now.
        the_vector_configs += [VectorConfiguration(list_of_normals)]
    if verbose:
        print("Done creating the vector configs. Now checking all regions", flush=True)

    acyclic_dict = {}
    with multiprocessing.Pool(processes=nb_proc) as pool:
        forcings = pool.map(wrapper_forcing_acyclic, the_vector_configs)
    if verbose:
        print("Done checking all regions. Cleaning up", flush=True)
    for i in range(len(R)):
        rc = repr_vector[i]
        vc = the_vector_configs[i]
        acyclic_dict[tuple([rc, vc])] = forcings[i][0]
    return acyclic_dict
