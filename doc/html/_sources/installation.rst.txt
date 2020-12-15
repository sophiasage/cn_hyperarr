.. nodoctest

Installation & Usage
====================

**Local install from source**


Download the source from the git repository::

    $ git clone https://github.com/sophiasage/cn_hyperarr.git

Change to the root directory and run::

    $ sage -pip install --upgrade --no-index -v .

For convenience this package contains a [makefile](makefile) with this
and other often used commands. Should you wish too, you can use the
shorthand::

    $ make install
    
**Usage**


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
