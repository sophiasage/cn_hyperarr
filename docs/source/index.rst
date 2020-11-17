=======================================================================================
Welcome to the cn_hyperarr's documentation!
=======================================================================================

cn_hyperarr is a package that Computation of congruence normality for hyperplane arrangements.

Installation
============

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
    sage: answer_to_ultimate_question()
    42


Congruence Normality for Hyperplane Arrangements
=======================================================================================

.. toctree::
   :maxdepth: 2

   main
   infinite_families
   vector_classes


Indices and Tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
