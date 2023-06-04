from pathlib import Path

import numpy as np
from metpy.calc import dewpoint_from_specific_humidity, most_unstable_parcel
from metpy.units import units
from pytest import approx

from src.ecape_py.ecape import (
    calc_ecape,
    calc_ecape_a,
    calc_el_height,
    calc_integral_arg,
    calc_lfc_height,
    calc_mse,
    calc_ncape,
    calc_psi,
    calc_sr_wind,
)

"""w/in 1% of authors calculations, given errors described below, when mucape is equivalent"""

def sample_data():
    """
    values via author's published matlab scripts (sigma = 1.6)
    https://figshare.com/articles/software/ECAPE_scripts/21859818
    """

    sounding_loc = Path("./sounding.txt")
    data = np.genfromtxt(sounding_loc, delimiter=",")

    height = data[:, 0] * units("m")  # height, m
    pressure = data[:, 1] * units("Pa")  # pressure, Pa
    temperature = data[:, 2] * units("K")  # temp, K
    specific_humidity = data[:, 3] * units("kg/kg")  # water vapor mass fraction (specific humidity), kg / kg
    u_wind = data[:, 4] * units("m/s")  # u wind, m / s
    v_wind = data[:, 5] * units("m/s")  # v wind, m / s

    dew_point_temperature = dewpoint_from_specific_humidity(pressure, temperature, specific_humidity)

    return height, pressure, temperature, specific_humidity, u_wind, v_wind, dew_point_temperature


def test_end_to_end_ecape():
    """
    values via author's published matlab scripts (sigma = 1.6)
    https://figshare.com/articles/software/ECAPE_scripts/21859818
    """

    height, pressure, temperature, specific_humidity, u_wind, v_wind, dew_point_temperature = sample_data()

    for cape_type in ['most_unstable', 'surface_based', 'mixed_layer']:
        ecape = calc_ecape(height, pressure, temperature, specific_humidity, u_wind, v_wind, cape_type)
        print(f"{cape_type}: {ecape}")

    assert ecape


def test_calc_psi():
    """
    values via author's published matlab scripts (sigma = 1.6)
    https://figshare.com/articles/software/ECAPE_scripts/21859818
    """

    el_z = 11750.0 * units("m")

    psi = calc_psi(el_z)
    assert psi.magnitude == approx(0.0034, rel=0.001)


def test_calc_ecape_a():
    """
    values via author's published matlab scripts (sigma = 1.6)
    https://figshare.com/articles/software/ECAPE_scripts/21859818
    """

    sr_wind = 16.662798431352986 * units('m/s')
    psi = 0.003401863644631 * units('dimensionless')
    ncape = 7.604878130037112e02 * units('m**2/s**2')
    cape = 3.530029673046427e03 * units('J/kg')

    ecape_a = calc_ecape_a(sr_wind, psi, ncape, cape)
    assert ecape_a.magnitude == approx(3.343908138651551e03 , rel=0.0001)


def test_calc_integral_arg():
    """
    values via author's published matlab scripts (sigma = 1.6)
    https://figshare.com/articles/software/ECAPE_scripts/21859818
    """

    intarg_loc = Path("./intarg.txt")
    data = np.genfromtxt(intarg_loc, delimiter=",")
    mseo_bar = data[:, 0] * units("J/kg")
    mseo_star = data[:, 1] * units("J/kg")
    t0 = data[:, 2] * units('K')
    int_arg = data[:, 3]

    integral_arg = calc_integral_arg(mseo_bar, mseo_star, t0)

    for test, verify in zip(integral_arg, int_arg):
        assert test.magnitude == approx(verify, rel=0.0001)


def test_calc_ncape():
    """
    values via author's published matlab scripts (sigma = 1.6)
    https://figshare.com/articles/software/ECAPE_scripts/21859818
    """

    ncape_loc = Path("./ncape.txt")
    data = np.genfromtxt(ncape_loc, delimiter=",")
    integral_arg = data[:, 0] * units('m/s**2')
    height = data[:, 1] * units('m')
    lfc_idx = 16
    el_idx = 117

    ncape = calc_ncape(integral_arg, height, lfc_idx, el_idx)
    assert ncape.magnitude == 7.604878130037112e02


def test_calc_mse():
    """
    values via author's published matlab scripts (sigma = 1.6)
    https://figshare.com/articles/software/ECAPE_scripts/21859818
    """

    height, pressure, temperature, specific_humidity, u_wind, v_wind, dew_point_temperature = sample_data()

    mseo_loc = Path("./mseo.txt")
    data = np.genfromtxt(mseo_loc, delimiter=",")
    mse0_bar = data[:, 1] * units('J/kg')
    mse0_star = data[:, 2] * units('J/kg')

    mseo_bar_sounding, mseo_star_sounding = calc_mse(pressure, height, temperature, specific_humidity)

    for test, verify in zip(mseo_bar_sounding, mse0_bar):
        assert test.magnitude == approx(verify, rel=.005)
    for test, verify in zip(mseo_star_sounding, mse0_star):
        assert test.magnitude == approx(verify, rel=.005)


def test_calc_sr_wind():
    """
    values via author's published matlab scripts (sigma = 1.6)
    https://figshare.com/articles/software/ECAPE_scripts/21859818

    Note: there is a difference in Bunkers right mover compenents
    between the author's compute_VSR.m and MetPy's metpy.calc.bunkers_storm_motion

    For the example sounding:
    compute_VSR.m: C_x: 15.634223483535706, C_y: 4.74162399790177, V_sr: 16.662798431352986
    metpy.calc.bunkers_storm_motion: C_x: 14.635206387655197, C_y: 4.668955528029753, V_sr: 15.929554223814558
    ~4.6 % lower V_sr -> ~1.1 % lower ECAPE result in this case
    Given modest effect, no further interrogation conducted, although I expect use of WCA vs. mean to be the cause

    Given equivalent bunkers-right components, the remainder of the function produces equivalent results
    MetPy's bunkers_storm_motion was used for code readability and consistency.
    """

    height, pressure, temperature, specific_humidity, u_wind, v_wind, dew_point_temperature = sample_data()

    sr_wind = calc_sr_wind(pressure, u_wind, v_wind, height)

    assert sr_wind.magnitude == 15.929554223814558


def test_calc_el_height():
    """
    Note: there is a difference in el & lfc idx/heights between the authors published scripts and
    the following calculations. This may be due to the selection of indexes or the el/lfc calculations.

    COMPUTE_ECAPE.m: MU_EL_idx: 117, MU_EL: 11750 -> ECAPE 3087
    calc_el_height.py: el_idx: 115, el_z: 11500 -> ECAPE: 3089

    Given a modest (~.1 %) impact on the resultant ECAPE (via the ncape & psi calculations)
    a process leveraging metpy.calc.el was used for code readability and consistency.
    """

    height, pressure, temperature, specific_humidity, u_wind, v_wind, dew_point_temperature = sample_data()

    el_idx, el_z = calc_el_height(pressure, height, temperature, dew_point_temperature, most_unstable_parcel)

    assert el_idx == 115
    assert el_z.magnitude == 11500


def test_calc_lfc_height():
    """
    Note: there is a difference in el & lfc idx/heights between the authors published scripts and
    the following calculations. This may be due to the selection of indexes or the el/lfc calculations.

    COMPUTE_ECAPE.m: MU_LFC_idx: 16, MU_EL: 1650 -> ECAPE 3086
    calc_el_height.py: lfc_idx: 18, lfc_z: 1800 -> ECAPE: 3089

    Given a modest (~.1 %) impact on the resultant ECAPE (via the ncape calculation)
    a process leveraging metpy.calc.lfc was used for code readability and consistency.
    """

    height, pressure, temperature, specific_humidity, u_wind, v_wind, dew_point_temperature = sample_data()

    lfc_idx, lfc_z = calc_lfc_height(pressure, height, temperature, dew_point_temperature, most_unstable_parcel)

    assert lfc_idx == 18
    assert lfc_z.magnitude == 1800
