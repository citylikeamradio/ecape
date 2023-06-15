# ecape

ecape is a simple module that contains an entraining CAPE, or ECAPE, calculation described by Peters et. al. 2023.
Peters-provided MatLab scripts serve as a reference and test verification data.
The module leans heavily on MetPy for meteorological calculations.

[![PyPI - Version](https://img.shields.io/pypi/v/ecape.svg)](https://pypi.org/project/ecape)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/ecape.svg)](https://pypi.org/project/ecape)

## Installation & Use

In console:

```console
pip install ecape
```

See the example in linked documentation.

```python
   from ecape.calc import calc_ecape
   ...
   ecape = calc_ecape(height, pressure, temperature, specific_humidity, u_wind, v_wind)
```

### Documentation & Source

https://github.com/citylikeamradio/ecape

https://citylikeamradio.github.io/ecape

### Contact

 - Robert Capella
 - bob.capella@gmail.com
 - https://twitter.com/minusthebob

Questions, comments, and feedback are certainly welcome. This project is a personal exercise
in learning how to publish packages to Github & PyPI, so excuse the excessive documentation for
one function.

### Future Work
 - add support for other water content variables
 - if useful, incorporate into MetPy
 - provide cli .nc, .csv, & aws support

### A note on undiluted CAPE & calculation accuracy
**If users prefer their own CAPE calculations, use the `undiluted_cape` parameter:**

When comparing calc_ecape.py & COMPUTE_ECAPE.m run on Peters 2023 sample data,
there is a ~10% difference in the resultant ECAPE. This is almost entirely due to a difference in calculated MUCAPE.
The tests describe other sources of variation (~1%).

Given:
 - the methods here are within ~1% of Peters' calculations when undiluted CAPE is equivalent
 - Peters et. al. specifically mention MetPy for determining undiluted CAPE
 - MetPy is a reliable, open-source, and frequently used meteorological calculation package

MetPy's undiluted CAPE calculations were chosen for ease of readability and implementation.


### References
Ryan M. May, Sean C. Arms, Patrick Marsh, Eric Bruning, John R. Leeman, Kevin Goebbert, Jonathan E. Thielen, Zachary S Bruick, and M. Drew. Camron. Metpy: a Python package for meteorological data. 2023. URL: Unidata/MetPy, doi:10.5065/D6WW7G29.

John Peters. ECAPE scripts. 2 2023. URL: https://figshare.com/articles/software/ECAPE_scripts/21859818, doi:10.6084/m9.figshare.21859818.v4.

Peters, J. M., D. R. Chavas, C. Su, H. Morrison, and B. E. Coffer, 2023: An analytic formula for entraining CAPE in mid-latitude storm environments. J. Atmos. Sci., https://doi.org/10.1175/JAS-D-23-0003.1, in press.

### Licence

`ecape` is distributed under the terms of the [MIT](https://spdx.org/licenses/MIT.html) license.
