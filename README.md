# ecape-py

ecape-py is a simple module that contains an entraining CAPE, or ECAPE, calculation described by Peters et. al. 2023.
Additionally, Peters et. al. 2023 -provided MatLab scripts serve as a reference and test verification data.
The module leans heavily on Metpy for meteorological calculations.

[![PyPI - Version](https://img.shields.io/pypi/v/ecape-py.svg)](https://pypi.org/project/ecape-py)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/ecape-py.svg)](https://pypi.org/project/ecape-py)

Installation
------------

```console
pip install ecape-py
```

Documentation, Source
-------------
https://github.com/citylikeamradio/ecape-py
https://citylikeamradio.github.io/ecape-py/html/index.html

Contact
-------------
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
Ryan M. May, Sean C. Arms, Patrick Marsh, Eric Bruning, John R. Leeman, Kevin Goebbert, Jonathan E. Thielen, Zachary S Bruick, and M. Drew. Camron. Metpy: a Python package for meteorological data. 2023. URL: Unidata/MetPy, doi:10.5065/D6WW7G29.

John Peters. ECAPE scripts. 2 2023. URL: https://figshare.com/articles/software/ECAPE_scripts/21859818, doi:10.6084/m9.figshare.21859818.v4.

John M. Peters, Daniel R. Chavas, Hugh Morrison, Chun-Yian Su, and Brice E. Coffer. An analytic formula for entraining cape in mid-latitude storm environments. 2023. arXiv:2301.04712.

Licence
------------

`ecape-py` is distributed under the terms of the [MIT](https://spdx.org/licenses/MIT.html) license.
