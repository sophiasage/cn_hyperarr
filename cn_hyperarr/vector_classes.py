# -*- coding: utf-8 -*-
r"""
Module for vector configurations and covectors. 
A vector configuration is a labeled set of vectors :math:`\{v_i~:~i\in I\}` in a
common vector space `V`. In the :class:`~vector_classes.VectorConfiguration` 
class, there are methods to test if a vector configuration is acyclic or totally 
acyclic, and also to return the underlying matroid and hypergraph. 
When the vector configuration is three dimensional, 
one can compute its cocircuits, i.e. maximal covectors.

A covector is a vector of signs -1,0,+1 (3=*, for restricted covectors)
obtained by an affine linear map on a vector configuration.

.. math::

    C_{c,a}:= (sign(c\cdot p_i+a))_{i\in[m]}

where `c` is a vector, `a` is a scalar, and `\{p_i\}_{i\in[m]}` is a vector
configuration. Members of the :class:`~vector_classes.Covector` class can be 
constructed either by declaring their entries or by passing a vector, scalar, 
and vector configuration. This class also has the 
:func:`~vector_classes.Covector.intersection` method as defined in [1].

There are also methods for examining congruence uniformity/normality of 
posets of regions of hyperplane arrangements.
A vector configuration can be seen as the set of normals to a hyperplane
arrangement. A simplicial hyperplane arrangement has a lattice of regions 
associated to each chamber. To check whether this lattice is obtainable through
a sequence of doublings of convex sets, i.e. is congruence normal, we use the 
theory developed in [1]. For each choice of chamber, there is an associated 
acyclic vector configuration. In this module, we compute the shard covectors
and the forcing oriented graph on the shard covectors. If this oriented 
graph is acyclic, then the arrangement is congruence normal with 
respect to the chosen chamber; this property is tested in the method 
:func:`~vector_classes.VectorConfiguration.is_congruence_normal`.

EXAMPLES::

    sage: from cn_hyperarr import *
    sage: tau = AA((1+sqrt(5))/2)
    sage: ncn = [[2*tau+1,2*tau,tau],[2*tau+2,2*tau+1,tau+1],[1,1,1],[tau+1,tau+1,tau],[2*tau,2*tau,tau],[tau+1,tau+1,1],[1,1,0],[0,1,0],[1,0,0],[tau+1,tau,tau]]
    sage: ncn_conf = VectorConfiguration(ncn);
    sage: len(ncn_conf.shard_covectors())
    29
    sage: ncn_conf.forcing_oriented_graph()
    Digraph on 29 vertices
    sage: ncn_conf.is_congruence_normal()
    False

REFERENCES:

    - [1] Michael Cuntz, Sophia Elia, and Jean-Philippe Labbé. Congruence normality of simplicial hyperplane arrangements via oriented matroids, 2020. arXiv:2009.14152.

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

from sage.misc.cachefunc import cached_method, cached_function
from itertools import combinations
from sage.modules.free_module_element import vector
from sage.misc.flatten import flatten
from sage.structure.sequence import Sequence
from sage.geometry.polyhedron.constructor import Polyhedron
from sage.functions.generalized import sign
from sage.sets.set import Set
from sage.graphs.digraph import DiGraph
from sage.matroids.constructor import Matroid
from sage.combinat.designs.incidence_structures import IncidenceStructure
from sage.matrix.constructor import matrix
from sage.matrix.special import identity_matrix
from sage.functions.other import sqrt
from sage.misc.functional import n
from sage.categories.cartesian_product import cartesian_product


class VectorConfiguration():
    r"""
    A vector configuration.

    A vector configuration is a labeled set of vectors `\{v_i~:~i\in I\}` in a
    common vector space `V`.

    EXAMPLES::

        sage: from cn_hyperarr import *
        sage: v = [[1, -1, -1], [1, 1, -1], [1, 1, 1], [1, -1, 1], [-1, -1, 1],
        ....:  [-1, -1, -1], [-1, 1, -1], [-1, 1, 1]]
        sage: vc = VectorConfiguration(v); vc
        Vector configuration of 8 vectors in dimension 3

    You can ask different things::

        sage: vc.vectors()
        ((1, -1, -1),
         (1, 1, -1),
         (1, 1, 1),
         (1, -1, 1),
         (-1, -1, 1),
         (-1, -1, -1),
         (-1, 1, -1),
         (-1, 1, 1))
        sage: vc.ambient_dimension()
        3
        sage: vc.three_dim_cocircuits()
        (((0, 1), (2, 3, 5, 6), (4, 7)),
         ((0, 3), (1, 2, 4, 5), (6, 7)),
         ((0, 5), (1, 3, 4, 6), (2, 7)),
         ((1, 2), (0, 3, 6, 7), (4, 5)),
         ((1, 6), (0, 2, 5, 7), (3, 4)),
         ((2, 3), (0, 1, 4, 7), (5, 6)),
         ((2, 7), (1, 3, 4, 6), (0, 5)),
         ((3, 4), (0, 2, 5, 7), (1, 6)),
         ((4, 5), (0, 3, 6, 7), (1, 2)),
         ((4, 7), (2, 3, 5, 6), (0, 1)),
         ((5, 6), (0, 1, 4, 7), (2, 3)),
         ((6, 7), (1, 2, 4, 5), (0, 3)))

    The vector configuration of an arrangement of hyperplanes that is not
    congruence normal. It is the arrangement A(10,60)_3 in [1]::

        sage: tau = AA((1+sqrt(5))/2)
        sage: ncn = [[2*tau+1,2*tau,tau],[2*tau+2,2*tau+1,tau+1],[1,1,1],[tau+1,tau+1,tau],[2*tau,2*tau,tau],[tau+1,tau+1,1],[1,1,0],[0,1,0],[1,0,0],[tau+1,tau,tau]]
        sage: ncn_conf = VectorConfiguration(ncn); ncn_conf
        Vector configuration of 10 vectors in dimension 3
        sage: ncn_conf.base_ring()
        Algebraic Real Field

    TESTS:

    The vectors should all have the same dimension::

        sage: VectorConfiguration([[0,1],[0,0,1]])
        Traceback (most recent call last):
        ...
        AssertionError: The vectors are not all of the same dimension

    The list of vectors can have repeats::

        sage: vc = VectorConfiguration([[0,0],[0,0]])
        sage: vc
        Vector configuration of 2 vectors in dimension 2
    """
    def __init__(self, vector_list, backend=None):
        r"""
        Construct a vector configuration.

        INPUT:

        - ``vector_list`` -- a list of vectors

        - ``backend`` -- a string, a polyhedral backend or``None`` (default)

        OUTPUT:

        A vector configuration.


        EXAMPLES::

            sage: from cn_hyperarr import *
            sage: vc = VectorConfiguration([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[0,1,1],[1,1,1]])
            sage: vc
            Vector configuration of 6 vectors in dimension 3
            sage: vc.ambient_dimension()
            3
            sage: vc.vectors()
            ((1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (0, 1, 1), (1, 1, 1))

        The same vector configuration with backend ``'normaliz'``::

            sage: vc = VectorConfiguration([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[0,1,1],[1,1,1]], backend = 'normaliz')  # optional - pynormaliz
            sage: vc                                    # optional - pynormaliz
            Vector configuration of 6 vectors in dimension 3
            sage: vc.backend()                          # optional - pynormaliz
            'normaliz'

        The vectors should all have the same dimension::

            sage: vc = VectorConfiguration([[1,0,0],[0,1,0,0]])
            Traceback (most recent call last):
            ...
            AssertionError: The vectors are not all of the same dimension

        The list of vectors can have repeats::

            sage: vc = VectorConfiguration([[0,0],[0,0]])
            sage: vc
            Vector configuration of 2 vectors in dimension 2

        TESTS:

        An empty vector configuration::

            sage: vc = VectorConfiguration([])
            sage: vc.vectors()
            ()
            sage: vc.ambient_dimension()
            -1
        """
        values = flatten([list(v) for v in vector_list])
        br = Sequence(values).universe()
        self._base_ring = br
        vector_list = tuple(vector(br, v, immutable=True) for v in vector_list)
        self._nb_vectors = len(vector_list)
        self._vectors = vector_list
        if self._nb_vectors == 0:
            self._dimension = -1
        else:
            self._dimension = len(self._vectors[0])
        # add that empty configuration has dim -1
        self._backend = backend
        if self._nb_vectors != 0:
            assert set([len(v) for v in self._vectors]) == set([self._dimension]), "The vectors are not all of the same dimension"

    def __getitem__(self, key):
        r"""
        Return the key-th element in the vector configuration.

        INPUT:

        - key -- an index

        OUTPUT:

        A vector.

        EXAMPLES::

            sage: from cn_hyperarr import *
            sage: vc = VectorConfiguration([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[0,1,1],[1,1,1]])
            sage: vc[3]
            (1, 1, 0)
        """
        return self._vectors.__getitem__(key)

    def __repr__(self):
        r"""
        String representation of the vector configuration

        OUTPUT:

        a string that describes the vector configuration.

        EXAMPLES::

            sage: from cn_hyperarr import *
            sage: vc = VectorConfiguration([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[0,1,1],[1,1,1]])
            sage: vc
            Vector configuration of 6 vectors in dimension 3
        """
        return "Vector configuration of {} vectors in dimension {}".format(self._nb_vectors, self._dimension)

    def ambient_dimension(self):
        r"""
        Return the ambient dimension of the vector configuration

        OUTPUT:

        An integer.

        EXAMPLES::

            sage: from cn_hyperarr import *
            sage: vc = VectorConfiguration([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[0,1,1],[1,1,1]])
            sage: vc.ambient_dimension()
            3
            sage: VectorConfiguration([]).ambient_dimension()
            -1
        """
        return self._dimension

    def base_ring(self):
        r"""
        Return the base ring of the vector configuration

        OUTPUT:

        A ring.

        EXAMPLES::

            sage: from cn_hyperarr import *
            sage: vc1 = VectorConfiguration([[1,0,0],[0,1,0],[3,0,0]])
            sage: vc1.base_ring()
            Integer Ring

            sage: vc2 = VectorConfiguration([[1/2,0,0],[0,1/3,0],[3,0,0]])
            sage: vc2.base_ring()
            Rational Field

            sage: vc3 = VectorConfiguration([[sqrt(2),1,0],[0,1,0],[0,0,1]])
            sage: vc3.base_ring()
            Symbolic Ring

            sage: K.<a> = QuadraticField(-5)
            sage: vc4 = VectorConfiguration([[a,1/2],[2,0],[4,5/6]])
            sage: vc4.base_ring()
            Number Field in a with defining polynomial x^2 + 5 with a = 2.236067977499790?*I

            sage: vc5 = VectorConfiguration([[1,0,0],[0,1,0],[3,0,1e-10]])
            sage: vc5.base_ring()
            Real Field with 53 bits of precision

            sage: vc6 = VectorConfiguration([])
            sage: vc6.base_ring()
            Category of objects
        """
        return self._base_ring

    def backend(self):
        r"""
        Return the backend used for the computations

        OUTPUT:

        A string stating the backend.

        EXAMPLES:

        The default backend is ``None``::

            sage: from cn_hyperarr import *
            sage: vc = VectorConfiguration([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[0,1,1],[1,1,1]])
            sage: vc.backend()

        You can specify any backend that Polyhedron objects accept::

            sage: vc = VectorConfiguration([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[0,1,1],[1,1,1]], backend='normaliz')    # optional - pynormaliz
            sage: vc.backend()   # optional - pynormaliz
            'normaliz'
        """
        return self._backend

    def vectors(self):
        r"""
        Return the list of vectors

        OUTPUT:

        A tuple.

        EXAMPLES::

            sage: from cn_hyperarr import *
            sage: vc = VectorConfiguration([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[0,1,1],[1,1,1]])
            sage: vc.vectors()
            ((1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (0, 1, 1), (1, 1, 1))
        """
        return self._vectors

    def n_vectors(self):
        r"""
        Return the number of vectors in the configuration

        OUTPUT:

        An integer.

        EXAMPLES::

            sage: from cn_hyperarr import *
            sage: vc = VectorConfiguration([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[0,1,1],[1,1,1]])
            sage: vc.n_vectors()
            6
        """
        return self._nb_vectors

    @cached_method
    def is_acyclic(self):
        r"""
        Return whether ``self`` is acyclic.

        A vector configuration is acyclic if the origin is not in the relative
        interior of the cone spanned by the vectors.

        OUTPUT:

        Boolean. Whether ``self`` is acyclic.

        .. SEEALSO::

            :meth:`is_totally_cyclic`.

        EXAMPLES::

            sage: from cn_hyperarr import *
            sage: vc = VectorConfiguration([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[0,1,1],[1,1,1]])
            sage: vc.is_acyclic()
            True
            sage: v = [[1, -1, -1],
            ....:      [1, 1, -1],
            ....:      [1, 1, 1],
            ....:      [1, -1, 1],
            ....:      [-1, -1, 1],
            ....:      [-1, -1, -1],
            ....:      [-1, 1, -1],
            ....:      [-1, 1, 1]]
            sage: vc = VectorConfiguration(v)
            sage: vc.is_acyclic()
            False

        TESTS::

            sage: vc = VectorConfiguration([])
            sage: vc.is_acyclic()
            True

            sage: vc2 = VectorConfiguration([[1,0,0],[-1,0,0]])
            sage: vc2.is_acyclic()
            False

            sage: vc3 = VectorConfiguration([[0,0,0]])
            sage: vc3.is_acyclic()
            Traceback (most recent call last):
            ...
            ValueError: PPL::ray(e):
            e == 0, but the origin cannot be a ray.

            sage: vc = VectorConfiguration([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[0,1,1],[1,1,1],[-1,0,0]])
            sage: vc.is_acyclic()
            False
        """
        acyc_vc = Polyhedron(rays=self, backend=self.backend())
        if acyc_vc.relative_interior_contains(vector([0] * acyc_vc.ambient_dim())):
            return False
        elif not acyc_vc.lines() == ():
            return False
        else:
            return True

    def is_totally_cyclic(self):
        r"""
        Return whether ``self`` is totally cyclic.

        A vector configuration is totally cyclic if the origin is in the relative
        interior of the cone spanned by the vectors.

        OUTPUT:

        Boolian. Whether ``self`` is totally cyclic.

        .. SEEALSO::

            :meth:`is_acyclic`.

        EXAMPLES::

            sage: from cn_hyperarr import *
            sage: vc = VectorConfiguration([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[0,1,1],[1,1,1]])
            sage: vc.is_totally_cyclic()
            False
            sage: v = [[1, -1, -1],
            ....:      [1, 1, -1],
            ....:      [1, 1, 1],
            ....:      [1, -1, 1],
            ....:      [-1, -1, 1],
            ....:      [-1, -1, -1],
            ....:      [-1, 1, -1],
            ....:      [-1, 1, 1]]
            sage: vc = VectorConfiguration(v)
            sage: vc.is_totally_cyclic()
            True

        TESTS::

            sage: vc = VectorConfiguration([])
            sage: vc.is_totally_cyclic()
            False

            sage: vc2 = VectorConfiguration([[1,0,0],[-1,0,0]])
            sage: vc2.is_totally_cyclic()
            True

            sage: vc3 = VectorConfiguration([[0,0,0]])
            sage: vc3.is_totally_cyclic()
            Traceback (most recent call last):
            ...
            ValueError: PPL::ray(e):
            e == 0, but the origin cannot be a ray.
        """
        return not self.is_acyclic()

    @cached_method
    def three_dim_cocircuits(self):
        r"""
        Return the cocircuits of a 3-dimensional vector configuration.

        OUTPUT:

        A tuple of cocircuits in the form (-,0,+).

        EXAMPLES::

            sage: from cn_hyperarr import *
            sage: vc = VectorConfiguration([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[0,1,1],[1,1,1]])
            sage: vc.three_dim_cocircuits()
            (((), (0, 1, 3), (2, 4, 5)),
             ((), (0, 2), (1, 3, 4, 5)),
             ((), (1, 2, 4), (0, 3, 5)),
             ((0,), (2, 3, 5), (1, 4)),
             ((0, 2, 5), (3, 4), (1,)),
             ((0, 3), (1, 5), (2, 4)),
             ((0, 3, 5), (1, 2, 4), ()),
             ((1,), (3, 4), (0, 2, 5)),
             ((1, 3), (0, 4, 5), (2,)),
             ((1, 3, 4, 5), (0, 2), ()),
             ((1, 4), (2, 3, 5), (0,)),
             ((2,), (0, 4, 5), (1, 3)),
             ((2, 4), (1, 5), (0, 3)),
             ((2, 4, 5), (0, 1, 3), ()))

        The vectors should be 3-dimensional::

            sage: vc = VectorConfiguration([[1,0,0,0],[0,1,0,0]])
            sage: vc.three_dim_cocircuits()
            Traceback (most recent call last):
            ...
            AssertionError: The ambient dimension is not 3

        Three linearly dependent vectors in dimension 3 form a cocircuit::

            sage: vc = VectorConfiguration([[1,0,0],[0,1,0],[2,-1,0]])
            sage: vc.three_dim_cocircuits()
            (((), (0, 1, 2), ()),)
        """
        assert self.ambient_dimension() == 3, "The ambient dimension is not 3"
        cocircuits = set()
        nb = self.n_vectors()
        for pair in combinations(self, 2):
            pt1, pt2 = pair
            Eqns = Polyhedron(rays=[pt1, pt2], backend=self._backend).equations()
            # if nb >= 28:
            #     Eqns = Polyhedron(rays=[pt1, pt2], backend=self._backend).equations()
            # else:
            #     Eqns = Polyhedron(rays=[pt1, pt2]).equations()
            if len(Eqns) == 1:  # Cheap way to check if they span a 2d-space
                Eq = Eqns[0]  # Get the normal vector of that 2d-space
                affine_eval = [sign(Eq.eval(i)) for i in self]
                minus = tuple([i for i in range(nb) if affine_eval[i] < 0])
                zero = tuple([i for i in range(nb) if affine_eval[i] == 0])
                plus = tuple([i for i in range(nb) if affine_eval[i] > 0])
                coc1 = tuple([minus, zero, plus])
                coc2 = tuple([plus, zero, minus])
                if coc1 not in cocircuits:
                    cocircuits.add(coc1)
                    cocircuits.add(coc2)
        return tuple(sorted(list(cocircuits)))

    @cached_method
    def underlying_matroid(self):
        r"""
        Return the underlying matroid of the vector configuration.

        OUTPUT:

        A matroid

        NOTE:

        Doesn't work if the matrix is not full dimensional
        why not use Matroid(matrix = mat) to define?

        EXAMPLES::

            sage: from cn_hyperarr import *
            sage: vc = VectorConfiguration([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[0,1,1],[1,1,1]])
            sage: vc.underlying_matroid()
            Matroid of rank 3 on 6 elements with 16 bases

        For a set of three linearly independent vectors in dimension three,
        there is just one basis::

            sage: vc = VectorConfiguration([[1,0,0],[0,1,0],[0,0,1]])
            sage: vc.underlying_matroid()
            Matroid of rank 3 on 3 elements with 1 bases
        """
        mat = matrix(self.vectors()).transpose()
        dim = mat.nrows()
        bases = [cols for cols in combinations(range(self.n_vectors()), dim) if mat.matrix_from_columns(cols).det() != 0]
        return Matroid(bases=bases)

    @cached_method
    def underlying_hypergraph(self):
        r"""
        Return the underlying hypergraph of cocircuits of the vector configuration.

        Edges of the hypergraph correspond to zeros in a cocircuit.

        OUTPUT:

        A hypergraph

        EXAMPLES::

            sage: from cn_hyperarr import *
            sage: vc = VectorConfiguration([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[0,1,1],[1,1,1]])
            sage: vc.underlying_hypergraph()
            Incidence structure with 6 points and 7 blocks

        For a set of three linearly dependent vectors in dimension three there
        is just one cocicuit and thus one block::

            sage: vc = VectorConfiguration([[1,0,0],[0,1,0],[2,-1,0]])
            sage: vc.underlying_hypergraph()
            Incidence structure with 3 points and 1 blocks
        """
        coc = set([v[1] for v in self.three_dim_cocircuits()])
        return IncidenceStructure(coc)

    @cached_method
    def dominating_pairs(self):
        r"""
        Return the dominating pairs for each point in a 3-d acyclic vector_configuration.

        A point p has a dominating pair (a,b) when p lies in the positive
        linear span of a and b.

        OUTPUT:

        A dictionary where keys are the indices of the vectors and the values
        are the associated dominating pairs.

        EXAMPLES::

            sage: from cn_hyperarr import *
            sage: vc = VectorConfiguration([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[0,1,1],[1,1,1]])
            sage: vc.dominating_pairs()
            {0: set(), 1: set(), 2: set(), 3: {(0, 1)}, 4: {(1, 2)}, 5: {(0, 4), (2, 3)}}
            sage: tau = AA((1+sqrt(5))/2)
            sage: ncn = [[2*tau+1,2*tau,tau],[2*tau+2,2*tau+1,tau+1],[1,1,1],[tau+1,tau+1,tau],[2*tau,2*tau,tau],[tau+1,tau+1,1],[1,1,0],[0,1,0],[1,0,0],[tau+1,tau,tau]]
            sage: ncn_conf = VectorConfiguration(ncn);
            sage: ncn_conf.dominating_pairs()
            {0: {(4, 8), (6, 9)},
             1: {(0, 2), (3, 8), (5, 9)},
             2: set(),
             3: {(2, 6), (7, 9)},
             4: {(1, 7), (2, 6)},
             5: {(0, 7), (2, 6)},
             6: {(7, 8)},
             7: set(),
             8: set(),
             9: {(2, 8)}}

        The vector [0,0,1] is in the positive linear span of [-1,0,1] and
        [1,0,1]::

            sage: vc = VectorConfiguration([[0,0,1],[-1,0,1],[1,0,1]])
            sage: vc.dominating_pairs()
            {0: {(1, 2)},
             1: set(),
             2: set()}
        """
        assert self.is_acyclic(), "The vector configuration should be acyclic"
        nb_pts = self.n_vectors()
        cocircuits = self.three_dim_cocircuits()
        dom_dict = {}

        for index in range(nb_pts):
            the_dominating_pairs = set()
            my_v = self[index]
            # We only consider the cocircuits that are zero at given index and contains at least 3 zero elements.
            relevant_cocircuits = (coc for coc in cocircuits if index in coc[1] and len(coc[1]) >= 3)
            for r_coc in relevant_cocircuits:
                the_cone = Polyhedron(rays=[self[j] for j in r_coc[1]], backend=self._backend)
                if the_cone.relative_interior_contains(my_v):
                    # This means that my_v is dominated.
                    # Now we catch the dominating rays:
                    # We are interested in the index of
                    # the rays and not just the ray (which
                    # may not be equal to the vector we started with
                    dom_rays = [j for j in r_coc[1] if not the_cone.relative_interior_contains(self[j])]
                    assert len(dom_rays) == 2, "problem with the cocircuits"
                    the_dominating_pairs.add(tuple(dom_rays))

            dom_dict[index] = the_dominating_pairs

        return dom_dict

    @cached_method
    def line_vertices(self):
        r"""
        Return the vertices for each line in a 3-d acyclic vector_configuration.

        NOTE:

        An acyclic vector configuration corresponds to the normal vectors of
        a hyperplane arrangement with a selected base region.

        OUTPUT:

        A list of pairs.

        EXAMPLES::

            sage: from cn_hyperarr import *
            sage: vc = VectorConfiguration([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[0,1,1],[1,1,1]])
            sage: vc.line_vertices()
            [(0, 1), (1, 2), (2, 3), (0, 4)]

        The vectors [-1,0,1] and [1,0,1] are the vertices of the line containing
        [0,0,1]::

            sage: vc = VectorConfiguration([[0,0,1],[-1,0,1],[1,0,1]])
            sage: vc.line_vertices()
            [(1, 2)]
        """
        assert self.is_acyclic(), "The vector configuration should be acyclic"
        nb_pts = self.n_vectors()
        cocircuits = self.three_dim_cocircuits()
        list_verts = []

        for coc in cocircuits:
            zero_indices = coc[1]
            if len(zero_indices) > 2:
                the_cone = Polyhedron(rays=[self[j] for j in zero_indices], backend=self._backend)
                the_verts = tuple([_ for _ in zero_indices if not the_cone.relative_interior_contains(self[_])])
                assert len(the_verts) == 2, "problem with the cocircuits"
                if the_verts not in list_verts:
                    list_verts += [the_verts]

        return list_verts

    def shard_covectors(self):
        r"""
        Return the shard covectors for the given planar vector configuraton.

        NOTE:

        An acyclic vector configuration corresponds to a hyperplane arrangement
        with selected base region.

        OUTPUT:

        A tuple of covectors as sequences of -1,0,1,'3'.

        EXAMPLES::

            sage: from cn_hyperarr import *
            sage: vc = VectorConfiguration([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[0,1,1],[1,1,1]])
            sage: len(vc.shard_covectors())
            11

        The following vector configuration is not congruence normal and has 29
        shards::

            sage: tau = AA((1+sqrt(5))/2)
            sage: ncn = [[2*tau+1,2*tau,tau],[2*tau+2,2*tau+1,tau+1],[1,1,1],[tau+1,tau+1,tau],[2*tau,2*tau,tau],[tau+1,tau+1,1],[1,1,0],[0,1,0],[1,0,0],[tau+1,tau,tau]]
            sage: ncn_conf = VectorConfiguration(ncn);
            sage: len(ncn_conf.shard_covectors())
            29

        The vector configuration [[-1,0,1],[1,0,1],[0,0,1],[0,1,0]] has five shards::

            sage: vc = VectorConfiguration([[-1,0,1],[1,0,1],[0,0,1],[0,1,0]])
            sage: len(vc.shard_covectors())
            5
        """
        nb_pts = self.n_vectors()
        dom_dict = self.dominating_pairs()
        cocircuits = self.three_dim_cocircuits()
        rms_covectors = []
        for index in range(nb_pts):
            the_pairs = dom_dict[index]
            if len(the_pairs) == 0:
                # The vector is not dominated
                # We put a joker everywhere
                # and set to 0 the entry at the index
                cov = [3] * nb_pts
                cov[index] = 0
                rms_covectors += [Covector(cov)]
            else:
                # The vector is dominated at least once
                the_restricted_set = flatten([list(p) for p in the_pairs])
                for pair in the_pairs:
                    p1, p2 = pair
                    the_indices = list(pair) + [index]
                    # Get the two cocircuits that contain the triple
                    the_coc_l = [coc for coc in cocircuits if all(i in coc[1] for i in the_indices)]
                    assert len(the_coc_l) == 2, "There is a problem with the cocircuits"
                    coc1, coc2 = the_coc_l
                    # Setup the two restricted covectors
                    restricted_coc1 = [[i for i in s if i in the_restricted_set] for s in coc1]
                    restricted_coc2 = [[i for i in s if i in the_restricted_set] for s in coc2]

                    # Creating the four possible fillings
                    # By placing p1 and p2 in the four possible ways
                    l_minres_coc = []
                    l_minres_coc += [[restricted_coc1[0] + [p1], [index], restricted_coc1[2] + [p2]]]
                    l_minres_coc += [[restricted_coc1[0] + [p2], [index], restricted_coc1[2] + [p1]]]
                    l_minres_coc += [[restricted_coc2[0] + [p1], [index], restricted_coc2[2] + [p2]]]
                    l_minres_coc += [[restricted_coc2[0] + [p2], [index], restricted_coc2[2] + [p1]]]

                    # From the indices, get the covectors
                    for r_coc in l_minres_coc:
                        cov = Covector(r_coc, nb_pts)
                        rms_covectors += [cov]
            rms_covectors = list(set(rms_covectors))
        return tuple(rms_covectors)

    def line_covectors(self):
        r"""
        Return a dictionary with pairs of normal vectors as keys and a list
        containing the line covectors as value (up to multiplication by -1).

        The line covectors encode the sign pattern of all the vectors in
        ``self`` on the intersection of two hyperplanes for vector
        configurations in dimension 3.

        OUTPUT

        a dictionary.

        EXAMPLES::

            sage: from cn_hyperarr import *
            sage: vc = VectorConfiguration([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[0,1,1],[1,1,1]])
            sage: vc.line_covectors()
            {{4, 5}: ((0, 1, -1, 1, 0, 0), (0, -1, 1, -1, 0, 0)),
             {3, 5}: ((1, -1, 0, 0, -1, 0), (-1, 1, 0, 0, 1, 0)),
             {3, 4}: ((1, -1, 1, 0, 0, 1), (-1, 1, -1, 0, 0, -1)),
             {2, 5}: ((-1, 1, 0, 0, 1, 0), (1, -1, 0, 0, -1, 0)),
             {2, 4}: ((-1, 0, 0, -1, 0, -1), (1, 0, 0, 1, 0, 1)),
             {2, 3}: ((-1, 1, 0, 0, 1, 0), (1, -1, 0, 0, -1, 0)),
             {1, 5}: ((1, 0, -1, 1, -1, 0), (-1, 0, 1, -1, 1, 0)),
             {1, 4}: ((1, 0, 0, 1, 0, 1), (-1, 0, 0, -1, 0, -1)),
             {1, 3}: ((0, 0, -1, 0, -1, -1), (0, 0, 1, 0, 1, 1)),
             {1, 2}: ((1, 0, 0, 1, 0, 1), (-1, 0, 0, -1, 0, -1)),
             {0, 5}: ((0, -1, 1, -1, 0, 0), (0, 1, -1, 1, 0, 0)),
             {0, 4}: ((0, -1, 1, -1, 0, 0), (0, 1, -1, 1, 0, 0)),
             {0, 3}: ((0, 0, 1, 0, 1, 1), (0, 0, -1, 0, -1, -1)),
             {0, 2}: ((0, -1, 0, -1, -1, -1), (0, 1, 0, 1, 1, 1)),
             {0, 1}: ((0, 0, 1, 0, 1, 1), (0, 0, -1, 0, -1, -1))}
            sage: tau = AA((1+sqrt(5))/2)
            sage: ncn = [[2*tau+1,2*tau,tau],[2*tau+2,2*tau+1,tau+1],[1,1,1],[tau+1,tau+1,tau],[2*tau,2*tau,tau],[tau+1,tau+1,1],[1,1,0],[0,1,0],[1,0,0],[tau+1,tau,tau]]

        TESTS:

        The ambient dimension should be 3::

            sage: vc = VectorConfiguration([[1], [2], [3]])
            sage: vc.line_covectors()
            Traceback (most recent call last):
            ...
            AssertionError: The ambient dimension should be 3

        TODO:

            Make it work in higher dimension.
        """
        assert self.ambient_dimension() == 3, "The ambient dimension should be 3"
        nb_hyperpl = self.n_vectors()
        hic_dict = {}
        for index_1, index_2 in combinations(range(nb_hyperpl), 2):
            new_key = Set([index_1, index_2])
            vector_1 = self[index_1]
            vector_2 = self[index_2]
            intersection = vector_1.cross_product(vector_2)
            new_value = list()
            for index in range(nb_hyperpl):
                dot_prod = intersection.dot_product(self[index])
                if dot_prod == 0:
                    new_value.append(0)
                elif dot_prod < 0:
                    new_value.append(-1)
                else:
                    new_value.append(1)
            hic1 = Covector(new_value)
            hic2 = Covector([-v for v in new_value])
            hic_dict[new_key] = tuple([hic1, hic2])
        return hic_dict

    @cached_method
    def forcing_oriented_graph(self):
        r"""
        Return the forcing oriented graph for a 3-dimensional vector
        configuration.

        This is the directed graph on shards. The hyperplane arrangement with
        base region corresponding to the acyclic vector configuration is
        congruence normal if and only if the directed graph on shards
        is acyclic.

        OUTPUT:

        An oriented graph

        EXAMPLES::

            sage: from cn_hyperarr import *
            sage: vc = VectorConfiguration([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[0,1,1],[1,1,1]])
            sage: fog = vc.forcing_oriented_graph(); fog
            Digraph on 11 vertices
            sage: fog.num_edges()
            24
            sage: tau = AA((1+sqrt(5))/2)
            sage: ncn = [[2*tau+1,2*tau,tau],[2*tau+2,2*tau+1,tau+1],[1,1,1],[tau+1,tau+1,tau],[2*tau,2*tau,tau],[tau+1,tau+1,1],[1,1,0],[0,1,0],[1,0,0],[tau+1,tau,tau]]
            sage: ncn_conf = VectorConfiguration(ncn);
            sage: ncn_fog = ncn_conf.forcing_oriented_graph(); ncn_fog
            Digraph on 29 vertices
            sage: ncn_fog.num_edges()
            96
            sage: ncn_fog.is_directed_acyclic()
            False

        TESTS:

        The ambient dimension should be 3::

            sage: vc = VectorConfiguration([[1],[2],[3]])
            sage: vc.forcing_oriented_graph()
            Traceback (most recent call last):
            ...
            AssertionError: The ambient dimension is not 3
        """
        nb_pts = self.n_vectors()
        shard_covectors = self.shard_covectors()
        hic = self.line_covectors()

        G = DiGraph()
        for shard_cov in shard_covectors:
            G.add_vertex(shard_cov)
            shard_index = shard_cov.index(0)
            potential_shards = (v for v in G.vertices() if v.index(0) != shard_index)
            for other_shard_cov in potential_shards:
                other_shard_index = other_shard_cov.index(0)
                # Get the potential ordering of forcing:
                forcing = True
                if other_shard_cov[shard_index] in [1, -1]:
                    # shard_cov -> other_shard_cov
                    first_cov, second_cov = shard_cov, other_shard_cov
                elif shard_cov[other_shard_index] in [1, -1]:
                    # other_shard_cov -> shard_cov
                    first_cov, second_cov = other_shard_cov, shard_cov
                else:
                    forcing = False
                if forcing:
                    # A forcing takes place
                    hic1, hic2 = hic[Set([shard_index, other_shard_index])]
                    sic1 = hic1 & shard_cov & other_shard_cov
                    sic2 = hic2 & shard_cov & other_shard_cov
                    if sic1 == hic1 or sic2 == hic2:
                        G.add_edge(first_cov, second_cov)
        return G.copy(immutable=True)

    def is_congruence_normal(self):
        r"""
        Return whether the poset of regions of the corresponding hyperplane
        arrangement is congruence normal.

        An acyclic vector configuration corresponds to a hyperplane arrangement
        with selected base region.

        OUTPUT:

        Boolean. Whether the poset of regions is congruence normal with respect
        to implicitly chosen base region.

        EXAMPLES::

            sage: from cn_hyperarr import *
            sage: vc = VectorConfiguration([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[0,1,1],[1,1,1]])
            sage: vc.is_congruence_normal()
            True

        The first arrangement that is not congruence normal::

            sage: tau = AA((1+sqrt(5))/2)
            sage: ncn = [[2*tau+1,2*tau,tau],[2*tau+2,2*tau+1,tau+1],[1,1,1],[tau+1,tau+1,tau],[2*tau,2*tau,tau],[tau+1,tau+1,1],[1,1,0],[0,1,0],[1,0,0],[tau+1,tau,tau]]
            sage: ncn_conf = VectorConfiguration(ncn);
            sage: ncn_conf.is_congruence_normal()
            False

        TESTS:

        The arrangement should be acyclic::

            sage: vc = VectorConfiguration([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[0,1,1],[1,1,1],[-1,0,0]])
            sage: vc.is_congruence_normal()
            Traceback (most recent call last):
            ...
            AssertionError: The vector configuration should be acyclic
        """
        assert self.is_acyclic(), "The vector configuration should be acyclic"
        forcing_graph = self.forcing_oriented_graph()
        return forcing_graph.is_directed_acyclic()

    def affine_basis(self):
        r"""
        Return the basis of an acylic vector configuration.

        These are the rays of the cone spanned by the vectors.
        This only works for simplicial vector configurations.

        OUTPUT:

        A set. The indices of the vectors that form the rays of the cone spanned
        by the vectors..

        EXAMPLES:

        The standard basis vectors are the rays of the following acyclic
        configuration::

            sage: from cn_hyperarr import *
            sage: vc = VectorConfiguration([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[0,1,1],[1,1,1]])
            sage: vc.affine_basis()
            {0, 1, 2}

        The cone spanned by the vectors should be simplicial::

            sage: vc = VectorConfiguration([[0,0,1],[1,0,1],[1,1,1],[0,1,1]])
            sage: vc.affine_basis()
            Traceback (most recent call last):
            ...
            AssertionError: The cone is not simplicial

        TESTS:

        The arrangement should be acyclic::

            sage: vc = VectorConfiguration([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[0,1,1],[1,1,1],[-1,0,0]])
            sage: vc.affine_basis()
            Traceback (most recent call last):
            ...
            AssertionError: The vector configuration should be acyclic

        The ambient dimension should be three::

            sage: vc = VectorConfiguration([[2,0,0,0],[0,2,0,0],[0,0,1,0],[0,0,0,1]])
            sage: vc.affine_basis()
            Traceback (most recent call last):
            ...
            AssertionError: The ambient dimension is not 3
        """
        # the arrangement should be acyclic and three dimensional
        assert self.is_acyclic(), "The vector configuration should be acyclic"
        assert self.ambient_dimension() == 3, "The ambient dimension is not 3"

        tdc = self.three_dim_cocircuits()
        border_lines = [set(c[1]) for c in tdc if len(c[0]) == 0]
        # check that the cone is simplicial
        assert len(border_lines) == 3, "The cone is not simplicial"
        return (border_lines[0] & border_lines[1]).union(border_lines[0] & border_lines[2]).union(border_lines[1] & border_lines[2])


class Covector(tuple):
    r"""
    A covector

    A covector is a vector of signs -1,0,+1 (3=*, for restricted covectors)
    obtained by an affine linear map on a vector configuration.

    .. MATH::

        C_{c,a}:= (sign(c\cdot p_i+a))_{i\in[m]}

    where `c` is a vector, `a` is a scalar, and `\{p_i\}_{i\in[m]}` is a vector
    configuration.

    EXAMPLES:

    It is possible to create a covector from various input. The first
    possibility is a list of signs and stars::

        sage: from cn_hyperarr import *
        sage: c = Covector([1,3,1,0,-1,3,0,-1,1]); c
        (1, 3, 1, 0, -1, 3, 0, -1, 1)
        sage: c.stars()
        (1, 5)

    The second possibility is giving the signs of the indices [-,0,+] and the
    length. The unspecified indices will have the star symbol::

        sage: d = Covector([[1,3,5],[2,8],[4,7]],10); d
        (3, -1, 0, -1, 1, -1, 3, 1, 0, 3)

    Finally, one may specify a vector configuration, a vector, a constant, and
    the starred indices::

        sage: vc = VectorConfiguration([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[0,1,1],[1,1,1],[0,0,0]])
        sage: c = vector([0,1,2])
        sage: a = -2
        sage: stars = [4]
        sage: Covector(vc,c,a,stars)
        (-1, -1, 0, -1, 3, 1, -1)
    """
    def __new__(cls, *data):
        r"""
        Create a covector.

        INPUT:

        - ``data`` -- either a tuple (vector_conf,c,a,starred), a vector of
          signs, or a tuple ((-,0,+),size) describing the sign vector.

        OUTPUT:

        a covector.

        EXAMPLES::

            sage: from cn_hyperarr import *
            sage: Covector([1,3,1,0,-1])
            (1, 3, 1, 0, -1)

            sage: Covector([1, 3, 1, 0, -1]) == Covector([[4], [3], [0, 2]],5)
            True
        """
        if len(data) == 1:  # From a vector
            covector = vector(*data)
            instance = super(Covector, cls).__new__(cls, covector)
            instance._vectorconf = None
            instance._normal = None
            instance._constant = None
            instance._stars = (i for i in range(len(covector)) if covector[i] == 3)
            instance._covector = vector(covector)
            return instance
        elif len(data) == 4:  # From vector configuration
            return cls._from_vector_conf(data)
        elif len(data) == 2:  # From a tuple ((-,0,+),size)
            return cls._from_signs(data)
        else:
            raise ValueError("the input data is not valid")

    @classmethod
    def _from_vector_conf(cls, data):
        r"""
        Create from a vector configuration

        INPUT:

        - ``data`` -- tuple of length 4 containing a vector configuration, a
          vector, a scalar constant and a set of starred values

        OUTPUT:

        A covector.

        EXAMPLES::

            sage: from cn_hyperarr import *
            sage: vc = VectorConfiguration([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[0,1,1],[1,1,1]])
            sage: c = vector([1,1,1])
            sage: a = -3/2
            sage: stars = [4]
            sage: Covector(vc,c,a,stars)
            (-1, -1, -1, 1, 3, 1)

        TESTS:

        If the vector is 0, the inner product will be 0::

            sage: vc = VectorConfiguration([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[0,1,1],[1,1,1]])
            sage: c = vector([0,0,0])
            sage: a = 0
            sage: stars = [4]
            sage: Covector(vc,c,a,stars) # indirect doctest
            (0, 0, 0, 0, 3, 0)
        """
        vectorconf = VectorConfiguration(data[0])
        normal = vector(data[1])
        constant = data[2]
        stars = tuple(data[3])
        covector = [sign(normal * v + constant) for v in vectorconf]
        for index in stars:
            covector[index] = 3
        instance = super(Covector, cls).__new__(cls, covector)
        instance._vectorconf = vectorconf
        instance._normal = normal
        instance._constant = constant
        instance._stars = stars
        instance._covector = vector(covector)
        return instance

    @classmethod
    def _from_signs(cls, data):
        r"""
        Create from signed indices

        INPUT:

        - ``data`` -- tuple of length 2 containing a triple (-,0,+) and the
          size of the vector

        OUTPUT:

        A covector.

        EXAMPLES::

            sage: from cn_hyperarr import *
            sage: Covector([[0, 2], [5], [4, 3]], 6)
            (-1, 3, -1, 1, 1, 0)

            sage: Covector([1,-1,0,3]) == Covector([[1],[2],[0]],4)  # indirect doctest
            True
        """
        signs, size = data
        covector = []
        for j in range(size):
            if j in signs[0]:
                covector += [-1]
            elif j in signs[1]:
                covector += [0]
            elif j in signs[2]:
                covector += [1]
            else:
                covector += [3]
        instance = super(Covector, cls).__new__(cls, covector)
        instance._vectorconf = None
        instance._normal = None
        instance._constant = None
        instance._stars = (i for i in range(len(covector)) if covector[i] == 3)
        instance._covector = vector(covector)
        return instance

    def stars(self):
        r"""
        Return the indices with stars

        OUTPUT:

        A tuple

        EXAMPLES::

            sage: from cn_hyperarr import *
            sage: Covector([0,3,3,0,-1,1]).stars()
            (1, 2)

        A covector with no stars::

            sage: Covector([1,-1,0,0,1]).stars()
            ()
        """
        return tuple(self._stars)

    def as_vector(self):
        r"""
        Return the covector as a vector object

        OUTPUT:

        A vector.

        EXAMPLES::

            sage: from cn_hyperarr import *
            sage: v = Covector([0,3,3,0,-1,1]).as_vector();v
            (0, 3, 3, 0, -1, 1)
            sage: 3*v
            (0, 9, 9, 0, -3, 3)
        """
        return self._covector

    @cached_method
    def intersection(self, other):
        r"""
        Return the intersection of ``self`` with ``other``.

        INPUT:

        - ``other`` -- a :class:`Covector`

        OUTPUT:

        A covector. The intersection of ``self`` with ``other``.

        EXAMPLES::

            sage: from cn_hyperarr import *
            sage: a = Covector([-1]*5)
            sage: b = Covector([1]*5)
            sage: a.intersection(b)
            (0, 0, 0, 0, 0)
            sage: a & b
            (0, 0, 0, 0, 0)
            sage: c = Covector([-1,0,1,0,-1])
            sage: a & c
            (-1, 0, 0, 0, -1)
            sage: b & c
            (0, 0, 1, 0, 0)
            sage: c & b
            (0, 0, 1, 0, 0)
            sage: c & a
            (-1, 0, 0, 0, -1)
            sage: a & b & c
            (0, 0, 0, 0, 0)
            sage: b & c & a
            (0, 0, 0, 0, 0)
            sage: d = Covector([3,3,-1,1,0])
            sage: a & d
            (-1, -1, -1, 0, 0)
            sage: b & d
            (1, 1, 0, 1, 0)
            sage: c & d
            (-1, 0, 0, 0, 0)

        TESTS::

            sage: a = Covector([1]*5)
            sage: b = Covector([-1]*4)
            sage: a & b
            Traceback (most recent call last):
            ...
            ValueError: The covectors don't have the same length
        """
        if len(self) != len(other):
            raise ValueError("The covectors don't have the same length")
        return Covector([inter_binary(self[i], other[i]) for i in range(len(self))])

    __and__ = intersection

##############################################################################
# Helper functions


@cached_function
def inter_binary(left, right):
    r"""
    Intersection operation on -,0,+,3, the 3 plays the role of a star.

    The commutative intersection operation is defined as follows:
    0 & {0, +, -, 3} = 0,
    * & {+, -, 3} = {+, -, 3}
    + & {+, -} = {+, 0}

    INPUT:

    - left, right -- elements of `\{-,0,+,3\}`

    OUTPUT:

    - element of `\{-,0,+,3\}`

    EXAMPLES::

        sage: from cn_hyperarr import *
        sage: from itertools import product
        sage: for i,j in product([-1,0,1,3],repeat=2):
        ....:     print(i,j,inter_binary(i,j))
        ....:
        -1 -1 -1
        -1 0 0
        -1 1 0
        -1 3 -1
        0 -1 0
        0 0 0
        0 1 0
        0 3 0
        1 -1 0
        1 0 0
        1 1 1
        1 3 1
        3 -1 -1
        3 0 0
        3 1 1
        3 3 3

    TESTS::

        sage: inter_binary(3,0)
        0
        sage: inter_binary(1,-1)
        0
        sage: inter_binary(3,1)
        1
    """
    if left == 0 or right == 0:  # Absorbant
        return 0
    elif left == 3:  # Star is the identity
        return right
    elif right == 3:  # Star is the identity
        return left
    elif left == right:  # Idempotent
        return left
    else:  # The only cases left are +,- and -,+
        return 0
