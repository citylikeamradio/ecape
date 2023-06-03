# SPDX-FileCopyrightText: 2023-present Robert Capella <bob.capella@gmail.com>
# SPDX-License-Identifier: MIT

"""Calculate the entraining CAPE (ECAPE) of a parcel"""
from typing import List

import metpy.calc as mpcalc
import numpy as np
import pint
from metpy.constants import dry_air_spec_heat_press, earth_gravity
from metpy.units import check_units, units


@check_units("[pressure]", "[temperature]", "[temperature]")
def _get_parcel_profile(pressure, temperature, dew_point_temperature, parcel_func=None):
    """

    :param pressure:
    :param temperature:
    :param dew_point_temperature:
    :param parcel_func:
    :return:
    """

    if parcel_func:
        parcel_p, parcel_t, parcel_td, *parcel_i = parcel_func(pressure, temperature, dew_point_temperature)
        parcel_profile = mpcalc.parcel_profile(pressure, parcel_t, parcel_td)
    else:
        parcel_profile = None

    return parcel_profile


@check_units("[pressure]", "[length]", "[temperature]", "[temperature]")
def calc_lfc_height(pressure, height, temperature, dew_point_temperature, parcel_func):
    """

    :param pressure:
    :param height:
    :param temperature:
    :param dew_point_temperature:
    :param parcel_func:
    :return:
    """

    parcel_profile = _get_parcel_profile(pressure, temperature, dew_point_temperature, parcel_func)

    lfc_p, lfc_t = mpcalc.lfc(pressure, temperature, dew_point_temperature, parcel_temperature_profile=parcel_profile)
    lfc_idx = (pressure - lfc_p > 0).nonzero()[0][-1]
    lfc_z = height[lfc_idx]
    return lfc_idx, lfc_z


@check_units("[pressure]", "[length]", "[temperature]", "[temperature]")
def calc_el_height(pressure, height, temperature, dew_point_temperature, parcel_func):
    """

    :param pressure:
    :param height:
    :param temperature:
    :param dew_point_temperature:
    :param parcel_func:
    :return:
    """
    parcel_profile = _get_parcel_profile(pressure, temperature, dew_point_temperature, parcel_func)

    el_p, el_t = mpcalc.el(pressure, temperature, dew_point_temperature, parcel_temperature_profile=parcel_profile)
    el_idx = (pressure - el_p > 0).nonzero()[0][-1]
    el_z = height[el_idx]
    return el_idx, el_z


@check_units("[pressure]", "[speed]", "[speed]", "[length]")
def calc_sr_wind(pressure, u_wind, v_wind, height):
    """

    :param pressure:
    :param u_wind:
    :param v_wind:
    :param height:
    :return:
    """

    bunkers_right = mpcalc.bunkers_storm_motion(pressure, u_wind, v_wind, height)[0]

    u_sr = u_wind - bunkers_right[0]
    v_sr = v_wind - bunkers_right[1]

    u_sr_1km = u_sr[np.nonzero(height <= 1000 * units("m"))]
    v_sr_1km = v_sr[np.nonzero(height <= 1000 * units("m"))]

    sr_wind = np.mean(mpcalc.wind_speed(u_sr_1km, v_sr_1km))
    return sr_wind


@check_units("[pressure]", "[length]", "[temperature]", "[mass]/[mass]")
def calc_mse(pressure, height, temperature, specific_humidity):
    """

    :param pressure:
    :param height:
    :param temperature:
    :param specific_humidity:
    :return:
    """

    moist_static_energy = mpcalc.moist_static_energy(height, temperature, specific_humidity)
    saturation_mixing_ratio = mpcalc.saturation_mixing_ratio(pressure, temperature)
    moist_static_energy_star = mpcalc.moist_static_energy(height, temperature, saturation_mixing_ratio)

    moist_static_energy_bar = np.cumsum(moist_static_energy) / np.arange(1, len(moist_static_energy) + 1)
    moist_static_energy_bar = moist_static_energy_bar.to("J/kg")
    moist_static_energy_star = moist_static_energy_star.to("J/kg")
    return moist_static_energy_bar, moist_static_energy_star


@check_units("[energy]/[mass]", "[energy]/[mass]", "[temperature]")
def calc_integral_arg(moist_static_energy_bar, moist_static_energy_star, temperature):
    """

    :param moist_static_energy_bar:
    :param moist_static_energy_star:
    :param temperature:
    :return:
    """

    integral_arg = -(earth_gravity / (dry_air_spec_heat_press * temperature)) * (
        moist_static_energy_bar - moist_static_energy_star
    )
    return integral_arg


@check_units("[length]/[time]**2", "[length]","[dimensionless]", "[dimensionless]")
def calc_ncape(integral_arg, height, lfc_idx, el_idx):
    """

    :param integral_arg:
    :param height:
    :param lfc_idx:
    :param el_idx:
    :return:
    """

    ncape = np.sum(
        (0.5 * integral_arg[lfc_idx:el_idx] + 0.5 * integral_arg[lfc_idx + 1 : el_idx + 1])
        * (height[lfc_idx + 1 : el_idx + 1] - height[lfc_idx:el_idx])
    )

    return ncape


@check_units("[speed]", "[dimensionless]", "[length]**2/[time]**2", "[energy]/[mass]")
def calc_ecape_a(sr_wind, pitchfork, ncape, cape):
    """

    :param sr_wind:
    :param pitchfork:
    :param ncape:
    :param cape:
    :return:
    """

    term_a = sr_wind**2 / 2.0
    term_b = (-1 - pitchfork - (2 * pitchfork / sr_wind**2) * ncape) / (4 * pitchfork / sr_wind**2)
    term_c = (
        np.sqrt(
            (1 + pitchfork + (2 * pitchfork / sr_wind**2) * ncape) ** 2
            + 8 * (pitchfork / sr_wind**2) * (cape - (pitchfork * ncape))
        )
    ) / (4 * pitchfork / sr_wind**2)
    ecape_a = term_a + term_b + term_c
    return ecape_a.to("J/kg") if ecape_a >= 0 else 0


@check_units("[length]")
def calc_pitchfork(el_z):
    """

    :param el_z:
    :return:
    """

    # additional constants
    sigma = 1.6 * units("dimensionless")
    alpha = 0.8 * units("dimensionless")
    l_mix = 120.0 * units("m")
    pr = (1.0 / 3.0) * units("dimensionless")  # prandtl number
    ksq = 0.18 * units("dimensionless")  # von karman constant

    pitchfork = (ksq * alpha**2 * np.pi**2 * l_mix) / (4.0 * pr * sigma**2 * el_z)
    return pitchfork


@check_units("[length]", "[pressure]", "[temperature]", "[mass]/[mass]", "[speed]", "[speed]")
def calc_ecape(
    height: List[pint.Quantity],
    pressure: List[pint.Quantity],
    temperature: List[pint.Quantity],
    specific_humidity: List[pint.Quantity],
    u_wind: List[pint.Quantity],
    v_wind: List[pint.Quantity],
    cape_type: str = "most_unstable",
) -> pint.Quantity:
    """
    Calculate the entraining CAPE (ECAPE) of a parcel

    Parameters:
    ----------
        height : list of 'pint.Quantity', length units
            Atmospheric heights at the levels given by 'pressure'.
        pressure : list of 'pint.Quantity, pressure units
            Total atmospheric pressure
        temperature : list of 'pint.Quantity', temperature units
            Air temperature
        specific humidity : list of 'pint.Quantity', mass/mass (or 'dimensionless') units
            Specific humidity
        u_wind : list of 'pint.Quantity', speed units
            X component of the wind
        v_wind : list of 'pint.Quantity', speed units
            Y component of the wind
        cape_type: str
            Variation of CAPE desired. 'most_unstable' (default), 'surface_based', or 'mixed_layer'

    Returns:
    ----------
        ecape : 'pint.Quantity'
            Entraining CAPE
    """

    cape_func = {
        "most_unstable": mpcalc.most_unstable_cape_cin,
        "surface_based": mpcalc.surface_based_cape_cin,
        "mixed_layer": mpcalc.mixed_layer_cape_cin,
    }

    parcel_func = {
        "most_unstable": mpcalc.most_unstable_parcel,
        "surface_based": None,
        "mixed_layer": mpcalc.mixed_parcel,
    }

    # calculate cape
    dew_point_temperature = mpcalc.dewpoint_from_specific_humidity(pressure, temperature, specific_humidity)
    cape, _ = cape_func[cape_type](pressure, temperature, dew_point_temperature)

    # calculate the level of free convection (lfc) and equilibrium level (el) indexes
    lfc_idx, _ = calc_lfc_height(pressure, height, temperature, dew_point_temperature, parcel_func[cape_type])
    el_idx, el_z = calc_el_height(pressure, height, temperature, dew_point_temperature, parcel_func[cape_type])

    # calculate the buoyancy dilution potential (ncape)
    moist_static_energy_bar, moist_static_energy_star = calc_mse(pressure, height, temperature, specific_humidity)
    integral_arg = calc_integral_arg(moist_static_energy_bar, moist_static_energy_star, temperature)
    ncape = calc_ncape(integral_arg, height, lfc_idx, el_idx)

    # calculate the storm relative (sr) wind
    sr_wind = calc_sr_wind(pressure, u_wind, v_wind, height)

    # calculate the entraining cape (ecape)
    pitchfork = calc_pitchfork(el_z)
    ecape_a = calc_ecape_a(sr_wind, pitchfork, ncape, cape)

    return ecape_a
