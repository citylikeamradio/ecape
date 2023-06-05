.. ecape documentation master file, created by
   sphinx-quickstart on Sat Jun  3 23:46:52 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.
   role :: raw-tex(raw)
   :format: latex html

.. toctree::
   :hidden:

   home <self>
   ecape <ecape>
   example

ecape
====================

ecape is a simple module that contains an entraining CAPE, or ECAPE, calculation described by :cite:t:`peters2023analytic`.
Additionally, :cite:t:`Peters2023` -provided MatLab scripts serve as a reference and test verification data.
The module leans heavily on MetPy [:cite:t:`metpy`] for meteorological calculations.

Package published via `Hatch <https://hatch.pypa.io/latest/>`_.

Installation & Use
-------------------

To use ecape, install it with pip:

.. code-block:: console

    pip install ecape

See the :ref:`example` page.

.. code-block:: python

   from ecape.calc import calc_ecape
   ...
   ecape = calc_ecape(height, pressure, temperature, specific_humidity, u_wind, v_wind, cape_type)

Source
-------------
https://github.com/citylikeamradio/ecape

Contact
---------
 - Robert Capella |:cloud_tornado:|
 - bob.capella@gmail.com |:cloud_lightning:|
 - https://twitter.com/minusthebob |:cloud_snow:|

Questions, comments, and feedback are certainly welcome. This project is a personal exercise
in learning how to publish packages to Github & PyPI, so excuse the excessive documentation for
one function |:smile:|.

Future Work
-------------
 - if useful, incorporate into MetPy
 - provide cli .nc, .csv, & aws support

Disclaimer
-------------
There is a ~10% difference in ECAPE between calc_ecape and Peters' published matlab scripts.
This is primarily due to a difference in calculated CAPE. The tests describe other sources of error.

Since:
 - the methods here are within ~1% of Peters' calculations when CAPE is equivalent in the sample data
 - Peters et. al. specifically mention MetPy for determining CAPE
 - MetPy is a reliable, open-source, and frequently used meteorological calculation package

MetPy's CAPE calculations were chosen for ease of readability and implementation.

References
------------
.. bibliography::
