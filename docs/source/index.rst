.. ecape documentation master file, created by
   sphinx-quickstart on Sat Jun  3 23:46:52 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.
   role :: raw-tex(raw)
   :format: latex html

.. toctree::
   :hidden:

   home <self>
   ecape_py <ecape_py>

ecape
====================

ecape is a simple module that contains an entraining CAPE, or ECAPE, calculation described by :cite:t:`peters2023analytic`.
Additionally, :cite:t:`Peters2023` -provided MatLab scripts serve as a reference and test verification data.
The module leans heavily on Metpy [:cite:t:`metpy`] for meteorological calculations.

Installation
--------------

To use ecape, install it with pip:

.. code-block:: console

    pip install ecape


Contact
---------
 - Robert Capella
 - bob.capella@gmail.com
 - https://twitter.com/minusthebob

Future Work
-------------
 - if useful, incorporate into MetPy
 - provide .nc, .csv, & aws support
 - provide just-in-time compiling support for faster 2D work

References
------------
.. bibliography::
