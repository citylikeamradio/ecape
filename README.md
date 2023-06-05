# ecape

ecape is a simple module that contains an entraining CAPE, or ECAPE, calculation described by Peters et. al. 2023.
Additionally, Peters et. al. 2023 -provided MatLab scripts serve as a reference and test verification data.
The module leans heavily on Metpy for meteorological calculations.

[![PyPI - Version](https://img.shields.io/pypi/v/ecape.svg)](https://pypi.org/project/ecape)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/ecape.svg)](https://pypi.org/project/ecape)

## Installation & Use


```console
pip install ecape
```

See the example in linked documentation:

```python
   from ecape.calc import calc_ecape
   ...
   ecape = calc_ecape(height, pressure, temperature, specific_humidity, u_wind, v_wind, cape_type)
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
 - if useful, incorporate into MetPy
 - provide cli .nc, .csv, & aws support

### Disclaimer
There is a ~10% difference in ECAPE between calc_ecape and Peters' published matlab scripts. 
This is primarily due to a difference in calculated CAPE. The tests describe other sources of error.

Since:
 - the methods here are within ~1% of Peters' calculations when CAPE is equivalent in the sample data
 - Peters et. al. specifically mention MetPy for determining CAPE
 - MetPy is a reliable, open-source, and frequently used meteorological calculation package

MetPy's CAPE calculations were chosen for ease of readability and implementation.

### References
Ryan M. May, Sean C. Arms, Patrick Marsh, Eric Bruning, John R. Leeman, Kevin Goebbert, Jonathan E. Thielen, Zachary S Bruick, and M. Drew. Camron. Metpy: a Python package for meteorological data. 2023. URL: Unidata/MetPy, doi:10.5065/D6WW7G29.

John Peters. ECAPE scripts. 2 2023. URL: https://figshare.com/articles/software/ECAPE_scripts/21859818, doi:10.6084/m9.figshare.21859818.v4.

John M. Peters, Daniel R. Chavas, Hugh Morrison, Chun-Yian Su, and Brice E. Coffer. An analytic formula for entraining cape in mid-latitude storm environments. 2023. arXiv:2301.04712.

### Licence

`ecape` is distributed under the terms of the [MIT](https://spdx.org/licenses/MIT.html) license.
