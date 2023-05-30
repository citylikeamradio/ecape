from metpy.units import units
import numpy as np
from pathlib import Path
from pytest import approx
from src.ecape_py.ecape import calc_ecape, calc_pitchfork, calc_ecape_a, calc_integral_arg, calc_ncape


def sample_data():
    sounding_loc = Path('./sounding.txt')
    data = np.genfromtxt(sounding_loc, delimiter=',')
    height = data[:, 0] * units('m')  # height, m
    pressure = data[:, 1] * units('Pa')  # pressure, Pa
    temperature = data[:, 2] * units('K')  # temp, K
    specific_humidity = data[:, 3] * units('kg/kg')  # water vapor mass fraction (specific humidity), kg / kg
    u_wind = data[:, 4] * units('m/s')  # u wind, m / s
    v_wind = data[:, 5] * units('m/s')  # v wind, m / s
    return height, pressure, temperature, specific_humidity, u_wind, v_wind


def test_end_to_end_ecape():
    height, pressure, temperature, specific_humidity, u_wind, v_wind = sample_data()
    ecape = calc_ecape(height, pressure, temperature, specific_humidity, u_wind, v_wind)
    assert ecape


def test_calc_pitchfork():
    # values via author's published matlab scripts (sigma = 1.6)
    # https://figshare.com/articles/software/ECAPE_scripts/21859818

    el_z = 11750.0 * units('m')
    pitchfork = calc_pitchfork(el_z)
    assert round(pitchfork, 4) == 0.0034 * units('dimensionless')


def test_calc_ecape_a():
    # values via author's published matlab scripts (sigma = 1.6)
    # https://figshare.com/articles/software/ECAPE_scripts/21859818

    sr_wind = 16.662798431352986
    pitchfork = 0.003401863644631
    ncape = 7.604878130037112e02
    cape = 3.530029673046427e03
    ecape_a = calc_ecape_a(sr_wind, pitchfork, ncape, cape)
    assert ecape_a == approx(3.343908138651551e03)


def test_calc_integral_arg():
    sounding_loc = Path('./intarg.txt')
    data = np.genfromtxt(sounding_loc, delimiter=',')
    mseo_bar = data[:, 0]
    mseo_star = data[:, 1]
    t0 = data[:, 2]
    int_arg = data[:, 3]

    integral_arg = calc_integral_arg(mseo_bar, mseo_star, t0)
    for i, x in enumerate(integral_arg.magnitude):
        assert x == approx(int_arg[i], abs=0.0001)


def test_calc_ncape():
    # values via author's published matlab scripts (sigma = 1.6)
    # https://figshare.com/articles/software/ECAPE_scripts/21859818

    sounding_loc = Path('./ncape.txt')
    data = np.genfromtxt(sounding_loc, delimiter=',')
    integral_arg = data[:, 0]
    height = data[:, 1]
    lfc_idx = 16
    el_idx = 117

    ncape = calc_ncape(integral_arg, height, lfc_idx, el_idx)
    assert ncape == 7.604878130037112e02
