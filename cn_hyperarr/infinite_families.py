# -*- coding: utf-8 -*-
r"""
Library containing the three infinite families of simplicial arrangements of rank three.

EXAMPLES:

The first infinite family is the family of near pencils. All lines except for
one share a common intersection point::

    sage: from cn_hyperarr.infinite_families import *
    sage: np4 = near_pencil_family(4,'normaliz'); np4             # optional - pynormaliz
    Vector configuration of 4 vectors in dimension 3

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
from sage.rings.rational_field import QQ


def near_pencil_family(n, backend=None):
    r"""
    Return the vector configuration of the near pencil with ``n`` lines.

    All lines except for one share a common intersection point.

    INPUT:

    - ``n`` -- integer. ``n`` `\geq 3`. The number of lines in the near pencil.
    
    - ``backend`` -- string (default = ``None``). The backend to use.

    OUTPUT:

    A vector configuration.

    EXAMPLES:

    The near pencil with three hyperplanes::

        sage: from cn_hyperarr.infinite_families import *
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


def family_two(n, backend=None):
    r"""
    Return the vector configuration of the simplicial arrangement
    `A(n,1)` from the family `\mathcal R(1)` in Grunbaum's list.

    The arrangement will have an ``n`` hyperplanes, with ``n`` even, consisting
    of the edges of the regular `n/2`-gon and the `n/2` lines of 
    mirror symmetry.

    INPUT:

    - ``n`` -- integer. ``n`` `\geq 6`. The number of lines in the arrangement.
    
    - ``backend`` -- string (default = ``None``). The backend to use.

    OUTPUT:

    A vector configuration.

    EXAMPLES::

        sage: from cn_hyperarr.infinite_families import *
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


def family_three(n, backend=None):
    r"""
    Return the vector configuration of the arrangement `A(n,2)` from the
    family `\mathcal R(2)` in Grunbaum's list.

    The arrangement will have ``n`` hyperplanes consisting of the edges of the
    regular, even `(n-1)/2`-gon and the `(n-1)/2` lines of mirror symmetry
    together with the line at infinity.
    Grunbaum's notation is `A(n,2)` for the resulting arrangement.

    INPUT:

    - ``n`` -- integer. (``n``-1)/2 must be even, and ``n`` `\geq` 9.
    
    - ``backend`` -- string (default = ``None``). The backend to use.

    OUTPUT:

    A vector configuration.

    EXAMPLES::

        sage: from cn_hyperarr.infinite_families import *
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
