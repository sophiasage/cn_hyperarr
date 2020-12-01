# -*- coding: utf-8 -*-
r"""
Module for the computations related to congruence normality/uniformity of posets of
regions of hyperplane arrangements.

This module has methods to test simplicial hyperplane arrangements of rank 3 
for congruence normality. Every region of each arrangement is tested to see
if the corresponding poset is congruence normal. 
There are also methods to double check the output and sort the arrangements into
categories of always, sometimes, and never congruence normal.

EXAMPLES:

This vector configuration is congruence normal for any choice of base region in
the corresponding hyperplane arrangement::

    sage: from cn_hyperarr import *
    sage: vc = VectorConfiguration([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[0,1,1],[1,1,1]]);
    sage: h = vectorconf_to_hyperplane_arrangement(vc);
    sage: h.n_regions()
    24
    sage: check_result = check_all_regions(vc);
    sage: check_result[0].values()
    dict_values([True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True])
    sage: len(check_result[0].values())
    24

REFERENCES:

    - [1] Michael Cuntz, Sophia Elia, and Jean-Philippe Labbé. Congruence normality of simplicial hyperplane arrangements via oriented matroids, 2020. arXiv:2009.14152.

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

import time, multiprocessing
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
    
    - ``backend`` -- string (default = None). The backend to use.

    OUTPUT:

    A hyperplane arrangement

    EXAMPLES:

    This arrangement with 10 hyperplanes is the smallest rank-three simplicial 
    arrangement that is not congruence normal::
    
        sage: from cn_hyperarr import *
        sage: tau = AA((1+sqrt(5))/2)
        sage: ncn = [[2*tau+1,2*tau,tau],[2*tau+2,2*tau+1,tau+1],[1,1,1],[tau+1,tau+1,tau],[2*tau,2*tau,tau],[tau+1,tau+1,1],[1,1,0],[0,1,0],[1,0,0],[tau+1,tau,tau]]
        sage: ncn_conf = VectorConfiguration(ncn);
        sage: ncn_arr = vectorconf_to_hyperplane_arrangement(ncn_conf); ncn_arr
        Arrangement of 10 hyperplanes of dimension 3 and rank 3

    Using normaliz as backend::

        sage: ncn_conf_norm = VectorConfiguration(ncn, 'normaliz')   # optional - pynormaliz
        sage: ncn_conf_norm.backend()                                # optional - pynormaliz
        'normaliz'
    """
    if not isinstance(vector_conf,VectorConfiguration):
        vector_conf = VectorConfiguration(vector_conf, backend=backend)
    if vector_conf.base_ring() == ZZ:
        H = HyperplaneArrangements(QQ,names='xyz')
    else:
        H = HyperplaneArrangements(vector_conf.base_ring(),names='xyz')
    x,y,z = H.gens()
    A = H(backend=vector_conf.backend())
    for v in vector_conf:
        a,b,c = v
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

    This arrangement with 10 hyperplanes is not congruence normal for the
    implicit choice of base region::

        sage: from cn_hyperarr import *
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
    return forcing.is_directed_acyclic(),forcing.num_verts()

##############################################################################
# Verification functions
##############################################################################

def RegionsCongruenceNormality(vector_conf, backend=None, verbose=False):
    r"""
    Return whether regions in the hyperplane arrangement generated by the
    vector configuration lead to congruence normal lattice of regions.
    
    INPUT:
    
    - ``vectorconf`` -- a vector configuration

    - ``backend`` -- string (default = ``None``).

    - ``verbose`` -- (default = False). If ``True`` return number of
      hyperplanes, normal vectors, number of regions.

    OUTPUT:

    A dictionary whose keys are tuple (vector,vc) where vector is a
    representative vector of a region of the hyperplane arrangement and vc is
    the corresponding acyclic vector configuration which is a reorientation of
    ``vector_conf``.

    EXAMPLES:

    An example of an arrangement that is always congruence normal::

        sage: from cn_hyperarr import *
        sage: vc = VectorConfiguration([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[0,1,1],[1,1,1]])
        sage: h = vectorconf_to_hyperplane_arrangement(vc);
        sage: h.n_regions()
        24
        sage: check_result = check_all_regions(vc);
        sage: check_result[0].values()
        dict_values([True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True])
        sage: len(check_result[0].values())
        24
    
    An arrangement that is not always congruence normal::

        sage: tau = AA((1+sqrt(5))/2)
        sage: ncn = [[2*tau+1,2*tau,tau],[2*tau+2,2*tau+1,tau+1],[1,1,1],[tau+1,tau+1,tau],[2*tau,2*tau,tau],[tau+1,tau+1,1],[1,1,0],[0,1,0],[1,0,0],[tau+1,tau,tau]]
        sage: ncn_conf = VectorConfiguration(ncn);
        sage: check = check_all_regions(ncn_conf)
        sage: vals_list = list(check[0].values())
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
    with multiprocessing.Pool(processes=8) as pool:
        forcings = pool.map(wrapper_forcing_acyclic, the_vector_configs)
    if verbose:
        print("Done checking all regions. Cleaning up", flush=True)
    for i in range(len(R)):
        rc = repr_vector[i]
        vc = the_vector_configs[i]
        acyclic_dict[tuple([rc,vc])] = forcings[i][0]
    return acyclic_dict

def report_on_congnorm(vector_conf, cuntz_index, folder):
    r"""
    This function computes the congruence normality for every choice of base
    region of the vector configuration ``vector_conf`` which has index
    ``cuntz_index`` in Cuntz's list of hyperplane arrangements.

    The standard output of the computations are written in a file named
    ``filename`` with Cuntz's index, and the results are saved in a Sage object

    INPUT:

    - ``vector_conf`` -- a vector configuration

    - ``cuntz_index`` -- an integer, the index of the vector configuration in
      Cuntz's list

    - ``folder`` -- a string, the name of the folder to write the output

    OUTPUT:

    A tuple (res_dict,nb_true,nb_false)

    TODO:

    Change cuntz_index and references to Cuntz' list

    EXAMPLES:

    An example of an arrangement that is always congruence normal::

        sage: from cn_hyperarr import *
        sage: vc = VectorConfiguration([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[0,1,1],[1,1,1]])
    """
    sys.stdout = open('{}/arr_{}.txt'.format(folder,cuntz_index), 'w')
    verbose = True
    backend = 'normaliz'
    nb_hyp = len(vector_conf)
    the_res_dict,nb_shards = check_all_regions(vector_conf, backend=backend, verbose=verbose)
    save(the_res_dict,'{}/dict_result_{}.sobj'.format(folder,cuntz_index))
    save(nb_shards,'{}/nb_reg_result_{}.sobj'.format(folder,cuntz_index))
    print("Saved the results.")

    nb_reg = len(the_res_dict.values())
    nb_true = sum(the_res_dict.values())
    all_true = all(the_res_dict.values())
    if verbose:
        print("The number of shards are: {}".format(nb_shards), flush=True)        
        if all_true:
            print("The lattice is always congruence normal: {}".format(all_true), flush=True)
        else:
            print("The lattice is congruence normal in {} cases.\n".format(nb_true) +
                  "The lattice is NOT congruence normal in {} cases.".format(nb_reg-nb_true), flush=True)
    sys.stdout.close()
    return True

def launch_verification(indices=None, folder=None, main_out=None):
    r"""
    This check congruence normality for all of the arrangements in
    the normal vectors list.

    The function computes the given indices (by default all indices) and 
    writes the output in the given ``folder``. The main standard output 
    is written in ``main_out``

    TODO:

        Rename files like recent_cuntz_list

    INPUT:

    - ``indices`` - a list (default=None). the indices in the file to check.

    - ``folder`` - a string (default=None). The folder to save the output. By 
      default the folder is named ``'output'``.

    - ``main_out`` - a string (default=None). The name of the file where the 
      main standard output is written. By defualt the file is 
      called ``'main_out'``.
    """
    if not folder:
        folder = 'output'
    if not main_out:
        main_out = '{}/main_out.txt'.format(folder)
    else:
        main_out = '{}/{}.txt'.format(folder,main_out)
    sys.stdout = open(main_out, 'w')
    print("Started computations on {}.".format(time.ctime()))
    print("Loading Cuntz's list...")
    load("recent_cuntz_list.sage")
    print("Done!")

    if not indices:
        nb = len(arrangements_list)
        indices = range(nb)
    else:
        nb = len(indices)
    print("There are {} arrangements to compute".format(nb))
    for new_index in range(nb):
        cuntz_index = indices[new_index]
        vector_conf = arrangements_list[cuntz_index]
        print("Shipping index {} (Cuntz's index: {}) on {}".format(new_index,cuntz_index,time.ctime()))
        sys.stdout.close()
        report_on_congnorm(vector_conf, cuntz_index, folder)
        sys.stdout = open(main_out, 'a')
        print("Finished computations on {}.".format(time.ctime()), flush=True)

    sys.stdout = sys.__stdout__
    print("Finished the FULL computations on {}.".format(time.ctime()), flush=True)
    return True

def check_output_dict(the_dict):
    r"""
    Checks whether the output fits Sage's results.

    Another way to check for congruence normality is to see if the poset is 
    obtainable through doublings of intervals.

    INPUT:

    - ``the_dict`` -- a dictionary. keys are tuples consisting of a vector
      and a vector configuration. 
    """
    print("Start a check")
    out_list = []
    prob_list = []
    for the_v,vc in the_dict.keys():
        the_ha = vectorconf_to_hyperplane_arrangement(vc,backend='normaliz')
        the_base_region = the_ha.region_containing_point(the_v)
        try:
            the_lattice = LatticePoset(the_ha.poset_of_regions(the_base_region))
            if the_lattice.is_semidistributive():
                our_res = the_dict[the_v,vc]
                print("Our computation gives: {}".format(our_res))
                sage_res = the_lattice.is_constructible_by_doublings('interval')
                print("Sage gives: {}".format(sage_res))
                the_prob = the_v,vc
                if our_res != sage_res:
                    print("There is a problem!")
                    prob_list += [the_prob]
                if not our_res:
                    print("We got our prize with {} hyperplanes".format(vc.n_vectors()))
                    out_list += [the_prob]
            else:
                print("For {}, the lattice is not SD".format(tuple([the_v,vc])))
        except:
            print("For {}, the poset is not a lattice.".format(tuple([the_v,vc])))

    if len(prob_list) == 0:
        print("Finished: all good.")
        if len(out_list) == 0:
            return True,[]
        else:
            print("Finished with not CU and tight")
            return True,out_list
    else:
        print("Finished with problems.")
        return False,prob_list

def filter_arrangements_into_families(indices=None,excluded_nb=[],parallel_bound=28):
    r"""
    Filters the arrangements into the three infinite families or "others".
    """
    print("Loading Cuntz's list...")
    load("recent_cuntz_list.sage")
    print("Done!")

    if not indices:
        nb = len(arrangements_list)
        indices = range(nb)
    else:
        nb = len(indices)
    print("There are {} arrangements to compute".format(nb))

    list_np = []
    list_poly = []
    list_poly_inf = []
    other = []

    # Prepare the computation
    vc_list = []
    needed_ns_para = []
    needed_ns_series = []
    for new_index in range(nb):
        cuntz_index = indices[new_index]
        vl = arrangements_list[cuntz_index]
        nb_v = len(vl) 
        if nb_v not in excluded_nb:
            if nb_v <= parallel_bound:
                if tuple([nb_v,None]) not in needed_ns_para:
                    vector_conf = VectorConfiguration(vl)
                    needed_ns_para += [tuple([nb_v,None])]
            else:
                if tuple([nb_v,"normaliz"]) not in needed_ns_series:
                    vector_conf = VectorConfiguration(vl,backend='normaliz')
                    needed_ns_series += [tuple([nb_v,"normaliz"])]
        vc_list += [vector_conf]

    print("The needed numbers are {}".format(needed_ns_para + needed_ns_series))

    # To be used without normaliz!
    with multiprocessing.Pool(processes=8) as pool:
        print("Computing pencils in parallel")
        np_mat = pool.map(near_pencil_hypergraph,needed_ns_para)
        near_pencils = {needed_ns_para[index][0]:np_mat[index] for index in range(len(needed_ns_para))}
        odd_para = [_ for _ in needed_ns_para if _[0] % 2 == 1]
        even_para = [_ for _ in needed_ns_para if _[0] % 2 == 0]
        print("Computing polygons in parallel")
        p_mat = pool.map(polygons_hypergraph,even_para)
        polygons = {even_para[index][0]:p_mat[index] for index in range(len(even_para))}
        print("Computing polygons inf in parallel")        
        p_inf_mat = pool.map(polygons_infinity_hypergraph,odd_para)
        polygons_inf = {odd_para[index][0]:p_inf_mat[index] for index in range(len(odd_para))}
    
    print("Computing pencils in series")
    np_mat = list(map(near_pencil_hypergraph,needed_ns_series))
    for index in range(len(needed_ns_series)):
        near_pencils[needed_ns_series[index][0]] = np_mat[index]
    odd_series = [_ for _ in needed_ns_series if _[0] % 2 == 1]
    even_series = [_ for _ in needed_ns_series if _[0] % 2 == 0]
    print("Computing polygons in series")
    p_mat = list(map(polygons_hypergraph,even_series))
    for index in range(len(even_series)):
        polygons[even_series[index][0]] = p_mat[index]
    print("Computing polygons inf in series")
    p_inf_mat = list(map(polygons_infinity_hypergraph,odd_series))
    for index in range(len(odd_series)):
        polygons_inf[odd_series[index][0]] = p_inf_mat[index]
    print("Done precomputing matroids")

    print("Creating inputs")
    the_inputs = []
    for new_index in range(nb):
        vector_conf = vc_list[new_index]
        nb_v = vector_conf.n_vectors()
        if nb_v not in excluded_nb:
            if nb_v % 2 == 0:
                input_data = tuple([vector_conf,near_pencils[nb_v],polygons[nb_v],0])
            else:
                input_data = tuple([vector_conf,near_pencils[nb_v],0,polygons_inf[nb_v]])
            the_inputs += [input_data]
        else:
            input_data = tuple([0,0,0,0])
            the_inputs += [input_data]

    with multiprocessing.Pool(processes=8) as pool:
        print("Launching the types fetching")
        types = pool.map(get_hg_type,the_inputs)

    for new_index in range(nb):
        cuntz_index = indices[new_index]
        the_type = types[new_index] 
        if the_type == 0:
            # Found a near pencil
            print("It is a near pencil!")
            list_np += [cuntz_index]
        elif the_type == 1:
            # Found a polygon
            print("It is a polygon!")
            list_poly += [cuntz_index]
        elif the_type == 2:
            # Found a polygon inf
            print("It is a polygon infinity!")
            list_poly_inf += [cuntz_index]
        else:
            # It does not fit 
            print("It is something else!")
            other += [cuntz_index]

    return list_np,list_poly,list_poly_inf,other

def get_hg_type(input_data):
    r"""
    Return the type of the vector configuration.
    """
    if input_data == tuple([0,0,0,0]):
        return 3
    else:
        vector_conf,np,p,p_inf = input_data
        nb_v = vector_conf.n_vectors()
        print("Computing the matroid with {} vectors".format(nb_v))
        new_hypergraph = vector_conf.underlying_hypergraph()
        print("Finished computing the matroid")
    
        if new_hypergraph.is_isomorphic(np):
            print("Found a near pencil with {} vectors".format(nb_v))
            return 0
        elif nb_v % 2 == 0:
            if new_hypergraph.is_isomorphic(p):
                print("Found a polygon with {} vectors".format(nb_v))
                return 1
            else:
                return 3
        else:
            if new_hypergraph.is_isomorphic(p_inf):
                print("Found a polygon inf with {} vectors".format(nb_v))
                return 2
            else:
                return 3

def get_mat_type(input_data):
    r"""
    Return the type of the vector configuration.
    """
    if input_data == tuple([0,0,0]):
        return 2
    else:
        vector_conf,p,p_inf = input_data
        nb_v = vector_conf.n_vectors()
        print("Computing the matroid with {} vectors".format(nb_v))
        new_matroid = vector_conf.underlying_matroid()
        print("Finished computing the matroid: {}".format(new_matroid))
    
        if nb_v % 2 == 0:
            if new_matroid.is_isomorphic(p):
                print("Found a polygon with {} vectors".format(nb_v))
                return 0
            else:
                return 2
        elif nb_v % 4 == 1:
            if new_matroid.is_isomorphic(p_inf):
                print("Found a polygon inf with {} vectors".format(nb_v))
                return 1
            else:
                return 2
        else:
            return 2
