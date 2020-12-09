# -*- coding: utf-8 -*-
r"""
Library containing the three infinite families of simplicial arrangements of rank three.

EXAMPLES:

The first infinite family is the family of near pencils. All lines except for
one share a common intersection point, as may be seen in the hypergraph::

    sage: from cn_hyperarr import *
    sage: np4 = near_pencil_family(4,'normaliz'); np4             # optional - pynormaliz
    Vector configuration of 4 vectors in dimension 3
    sage: np4_hypergraph = near_pencil_hypergraph(4, 'normaliz'); # optional - pynormaliz
    sage: np4_hypergraph.block_sizes()                              # optional - pynormaliz
    [3, 2, 2, 2]

The second infinite family comes from regular polygons and their lines of
mirror symmetry. They are defined by the total number of hyperplanes, equal to
twice the number of edges of the polygon::

    sage: ftwo_6 = family_two(6,'normaliz'); ftwo_6               # optional - pynormaliz
    Vector configuration of 6 vectors in dimension 3
    sage: ftwo_10 = family_two(10,'normaliz'); ftwo_10            # optional - pynormaliz
    Vector configuration of 10 vectors in dimension 3

The third infinite family comes from regular polygons with an even number of
edges and their lines of mirror symmetry along with the line at infinity::

    sage: fthree_17 = family_three(17, 'normaliz'); fthree_17     # optional - pynormaliz
    Vector configuration of 17 vectors in dimension 3

REFERENCES:

    - [1] Michael Cuntz, Sophia Elia, and Jean-Philippe Labbé. Congruence
      normality of simplicial hyperplane arrangements via oriented matroids,
      2020. arXiv:2009.14152.
    - [2] Branko Grunbaum. A catalogue of simplicial arrangements in the real
      projective plane, 2009. Ars Math. Contemp. 2, no. 1, 1-25.

AUTHORS:

- Jean-Philippe Labbé (2020): Initial version
- Sophia Elia (2020): Initial version
"""

##############################################################################
#     Copyright (C) 2020 Jean-Philippe Labbe <labbe at math.fu-berlin.de>
#                   2020 Sophia Elia         < sophiae56 at math.fu-berlin.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
##############################################################################

from sage.geometry.polyhedron.library import polytopes
from sage.geometry.polyhedron.constructor import Polyhedron
from sage.rings.qqbar import QQbar
from .vector_classes import *
from sage.misc.persist import load, save
from sage.rings.rational_field import QQ
from sage.combinat.designs.incidence_structures import IncidenceStructure


def near_pencil_family(n, backend=None):
    r"""
    Return the vector configuration of the near pencil with ``n`` lines.

    INPUT:

    - ``n`` -- integer. ``n`` `\geq 3`. The number of lines in the near pencil.
    
    - ``backend`` -- string (default = ``None``). The backend to use.

    OUTPUT:

    A vector configuration.

    EXAMPLES:

    The near pencil with three hyperplanes::

        sage: from cn_hyperarr import *
        sage: np = near_pencil_family(3,'normaliz'); np         # optional - pynormaliz
        Vector configuration of 3 vectors in dimension 3

    The near pencil arrangements are always congruence normal. Note, this test
    just shows there exists a region such that they are CN::

        sage: near_pencils = [near_pencil_family(n,'normaliz') for n in range(3,6)] # optional - pynormaliz
        sage: [ np.is_congruence_normal() for np in near_pencils] # optional - pynormaliz
        [True, True, True]

    TESTS:

    Test that the backend is normaliz::

        sage: np = near_pencil_family(3,'normaliz');            # optional - pynormaliz
        sage: np.backend()                                        # optional - pynormaliz
        'normaliz'
    """
    z = QQbar.zeta(2*n)
    vecs = [[(z**k).real(), (z**k).imag(), 0] for k in range(1, n)]
    vecs += [vector([0, -1, 1])]
    return VectorConfiguration(vecs, backend=backend)

def near_pencil_matroid(n, backend=None):
    r"""
    Return the matroid of the near pencil with ``n`` lines.

    INPUT:

    - ``n`` -- integer.  ``n`` `\geq 3`. The number of lines in the near pencil.

    - ``backend`` -- string (default = ``None``). The backend to use.

    OUTPUT:

    a matroid.


    EXAMPLES:

    The matroid of the near pencil with three hyperplanes::

        sage: from cn_hyperarr import *
        sage: np_mat = near_pencil_matroid(3, 'normaliz'); np_mat   # optional - pynormaliz
        Matroid of rank 3 on 3 elements with 1 bases

    TESTS::

        sage: list(np_mat.cocircuits())                             # optional - pynormaliz
        [frozenset({0}), frozenset({1}), frozenset({2})]
    """
    try:
        # print("Try loading near pencil matroid {} from library".format(n))
        mat = load('pencil_matroid_{}'.format(n))
        # print("Found it.")
    except FileNotFoundError:
        # print("Failed. Compute {}.".format(n))
        mat = near_pencil_family(n, backend).underlying_matroid()
        # print("Finished computing np {}: {}. Saving...".format(n, mat))
        save(mat, "pencil_matroid_{}.sobj".format(n))
        # print("Finished saving {}.".format(n))
    return mat


def near_pencil_hypergraph(n,backend=None):
    r"""
    Return the hypergraph of the near pencil with ``n`` lines.

    This hypergraph has ``n`` vertices corresponding to the ``n`` lines of the
    arrangement. The blocks correspond to the intersections.

    INPUT:

    - ``n`` -- integer. ` ``n`` \geq 3`. The number of lines in the near pencil.

    - ``backend`` -- string (default = ``None``). The backend to use.

    OUTPUT:

    a hypergraph.

    EXAMPLES::

        sage: from cn_hyperarr import *
        sage: np_hyp = near_pencil_hypergraph(3,'normaliz'); np_hyp   # optional - pynormaliz
        Incidence structure with 3 points and 3 blocks

    TESTS::

        sage: near_pencil_hypergraph(6, None).blocks()                  # optional - pynormaliz
        [[0, 1, 2, 3, 4], [0, 5], [1, 5], [2, 5], [3, 5], [4, 5]]
    """
    try:
        # print("Try loading near pencil hg {} from library".format(n))
        hg = load('pencil_hg_{}'.format(n))
        # print("Found it.")
    except FileNotFoundError:
        # print("Failed. Compute {}.".format(n))
        hg = IncidenceStructure([(n-1, i) for i in range(n-1)] + [tuple(range(n-1))])
        # hg = near_pencil_family(n, backend).underlying_hypergraph()
        # print("Finished computing np {}: {}. Saving...".format(n, hg))
        save(hg, "pencil_hg_{}.sobj".format(n))
        # print("Finished saving {}.".format(n))
    return hg


def family_two(n, backend=None):
    r"""
    Return the vector configuration of the simplicial arrangement
    `\A(n,1)` from the family `\mathcal R(1)` in Grunbaum's list.

    The arrangement will have an even number of hyperplanes consisting
    of the edges of the regular ` ``n/2`-gon and the ` ``n``/2` lines of 
    mirror symmetry.

    INPUT:

    - ``n`` -- integer. ` ``n`` \geq 6`. The number of lines in the arrangement.
    
    - ``backend`` -- string (default = ``None``). The backend to use.

    OUTPUT:

    A vector configuration.

    EXAMPLES::

        sage: from cn_hyperarr import *
        sage: pf = family_two(8,'normaliz'); pf   # optional - pynormaliz
        Vector configuration of 8 vectors in dimension 3

    The number of lines must be even::

        sage: pf3 = family_two(3,'normaliz');     # optional - pynormaliz
        Traceback (most recent call last):
        ...
        AssertionError: n must be even

    The number of lines must be at least 6::

        sage: pf4 = family_two(4,'normaliz')      # optional - pynormaliz
        Traceback (most recent call last):
        ...
        ValueError: n (=2) must be an integer greater than 2
    """
    assert n % 2 == 0, "n must be even"
    reg_poly = polytopes.regular_polygon(n/QQ(2), backend='normaliz')
    reg_cone = Polyhedron(rays=[list(v.vector()) + [1] for v in reg_poly.vertices()], backend=backend)
    vecs = [h.A() for h in reg_cone.Hrepresentation()]

    z = QQbar.zeta(n)
    vecs += [[(z**k).real(), (z**k).imag(), 0] for k in range(n/QQ(2))]
    return VectorConfiguration(vecs, backend=backend)


def family_two_matroid(n,backend=None):
    r"""
    Return the matroid of the simplicial arrangement
    `\A(n,1)` from the family `\mathcal R(1)` in Grunbaum's list.

    The arrangement will have an even number of hyperplanes consisting
    of the edges of the regular ``n``/2-gon and the ``n``/2 lines of mirror symmetry.

    INPUT:

    - ``n`` -- integer. ``n`` `\geq 6`. The number of lines in the arrangement.
 
    - ``backend`` -- string (default = ``None``). The backend to use.

    OUTPUT:

    A matroid.

    EXAMPLES::

        sage: from cn_hyperarr import *
        sage: pm = family_two_matroid(6,'normaliz'); pm        # optional - pynormaliz
        Matroid of rank 3 on 6 elements with 16 bases

    TESTS::

        sage: pm10 = family_two_matroid(10, 'normaliz'); pm10  # optional - pynormaliz
        Matroid of rank 3 on 10 elements with 100 bases
    """
    try:
        # print("Try loading polygon matroid {} from library".format(n))
        mat = load('polygon_matroid_{}'.format(n))
        # print("Found it.")
    except FileNotFoundError:
        # print("Failed. Compute {}.".format(n))
        mat = family_two(n,backend).underlying_matroid()
        # print("Finished computing poly {}: {}. Saving...".format(n, mat))
        save(mat, "polygon_matroid_{}.sobj".format(n))
        # print("Finished saving {}.".format(n))
    return mat


def family_two_hypergraph(n,backend=None):
    r"""
    Return the hypergraph of the simplicial arrangement
    `\A(n,1)` from the family `\mathcal R(1)` in Grunbaum's list.

    The arrangement will have an even number of hyperplanes consisting
    of the edges of the regular ``n``/2-gon and the ``n``/2 lines of mirror symmetry.

    The blocks of the hypergraph correspond to intersections of the lines.

    INPUT:

    - ``n`` -- integer. ``n`` `\geq 6`. The number of lines in the arrangement.

    - ``backend`` -- string (default = ``None``). The backend to use.

    OUTPUT:

    A hypergraph.

    EXAMPLES::

        sage: from cn_hyperarr import *
        sage: ph6 = family_two_hypergraph(6, 'normaliz'); ph6   # optional - pynormaliz
        Incidence structure with 6 points and 7 blocks
        sage: ph6.blocks()                                        # optional - pynormaliz
        [[0, 1, 4], [0, 2, 3], [0, 5], [1, 2, 5], [1, 3], [2, 4], [3, 4, 5]]

    TESTS::

        sage: ph6.ground_set()                                    # optional - pynormaliz
        [0, 1, 2, 3, 4, 5]
    """
    try:
        # print("Try loading polygon hg {} from library".format(n))
        hg = load('polygon_hg_{}'.format(n))
        # print("Found it.")
    except FileNotFoundError:
        # print("Failed. Compute {}.".format(n))
        hg = family_two(n,backend).underlying_hypergraph()
        # print("Finished computing poly {}: {}. Saving...".format(n, hg))
        save(hg, "polygon_hg_{}.sobj".format(n))
        # print("Finished saving {}.".format(n))

    return hg


def family_three(n, backend=None):
    r"""
    Return the vector configuration of the arrangement `A`(``n``,2) from the
    family `\mathcal R(2)` in Grunbaum's list.

    The arrangement will have ``n`` hyperplanes consisting of the edges of the
    regular even (``n``-1)/2-gon and the (``n``-1)/2 lines of mirror symmetry
    together with the line at infinity.
    Grunbaum's notation is `A`(``n``,2) for the resulting arrangement.

    INPUT:

    - ``n`` -- integer. (``n``-1)/2 must be even. and ``n`` `\geq` 9.
    
    - ``backend`` -- string (default = ``None``). The backend to use.

    OUTPUT:

    A vector configuration.

    EXAMPLES::

        sage: from cn_hyperarr import *
        sage: f3_9 = family_three(9,'normaliz'); f3_9    # optional - pynormaliz
        Vector configuration of 9 vectors in dimension 3

    The number of lines minus one should be zero mod 4::

        sage: f3_8 = family_three(8, 'normaliz')         # optional - pynormaliz
        Traceback (most recent call last):
        ...
        AssertionError: (n-1) must be 0 mod 4
    """
    assert (n-1) % 4 == 0, "(n-1) must be 0 mod 4"
    li = family_two(n-1, backend).vectors()
    li += tuple([vector([0, 0, 1])])
    return VectorConfiguration(li, backend=backend)


def family_three_matroid(n, backend=None):
    r"""
    Return the matroid of the arrangement `A`(``n``,2) from the
    family `\mathcal R`(2) in Grunbaum's list.

    The arrangement will have ``n`` hyperplanes consisting of the edges of the
    regular even (``n``-1)/2-gon and the (``n``-1)/2 lines of mirror symmetry
    together with the line at infinity.
    Grunbaum's notation is `A`(``n``,2) for the resulting arrangement.

    INPUT:

    - ``n`` -- integer. (``n``-1) must be 0 mod 4, and ``n``>= 9.

    - ``backend`` -- string (default = ``None``). The backend to use.

    OUTPUT:

    A matroid.

    EXAMPLES::

        sage: from cn_hyperarr import *
        sage: f3_9_mat = family_three_matroid(9, 'normaliz'); f3_9_mat  # optional - pynormaliz
        Matroid of rank 3 on 9 elements with 68 bases

        sage: f13 = family_three_matroid(13, None); f13
        Matroid of rank 3 on 13 elements with 242 bases
    """
    try:
        # print("Try loading polygon inf matroid {} from library".format(n))
        mat = load('polygon_inf_matroid_{}'.format(n))
        # print("Found it.")
    except FileNotFoundError:
        # print("Failed. Compute {}.".format(n))
        mat = family_three(n,backend).underlying_matroid()
        # print("Finished computing poly inf {}: {}. Saving...".format(n, mat))
        save(mat, "polygon_inf_matroid_{}.sobj".format(n))
        # print("Finished saving {}.".format(n))
    return mat


def family_three_hypergraph(n, backend=None):
    r"""
    Return the hypergraph of the arrangement `A`(``n``,2) from the
    family `\mathcal R(2)` in Grunbaum's list.

    The arrangement will have ``n`` hyperplanes consisting of the edges of the
    regular even (``n``-1)/2-gon and the (``n``-1)/2 lines of mirror symmetry
    together with the line at infinity.
    Grunbaum's notation is `A`(``n``,2) for the resulting arrangement.

    The hypergraph has blocks corresponding to the intersections of the lines.

    INPUT:

    - ``n`` -- integer. (``n``-1) must be 0 mod 4 and ``n`` `\geq` 9.
      
    - ``backend`` -- string (default = ``None``). The backend to use.

    OUTPUT:

    a hypergraph.

    EXAMPLES::

        sage: from cn_hyperarr import *
        sage: f3_9 = family_three_hypergraph(9, 'normaliz'); f3_9   # optional - pynormaliz
        Incidence structure with 9 points and 13 blocks

    The square has four vertices where three lines intersect. The center of the
    square and two points at infinity also have 4 lines intersecting. All other
    instersections involve two lines::

        sage: f3_9.block_sizes()                                      # optional - pynormaliz
        [3, 3, 4, 2, 4, 3, 2, 3, 2, 2, 4, 2, 2]
    """
    try:
        # print("Try loading polygon inf hg {} from library".format(n))
        hg = load('polygon_inf_hg_{}'.format(n))
        # print("Found it.")
    except FileNotFoundError:
        # print("Failed. Compute {}.".format(n))
        hg = family_three(n,backend).underlying_hypergraph()
        # print("Finished computing poly inf {}: {}. Saving...".format(n, hg))
        save(hg, "polygon_inf_hg_{}.sobj".format(n))
        # print("Finished saving {}.".format(n))

    return hg
