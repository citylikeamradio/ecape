# SPDX-FileCopyrightText: 2023-present Robert Capella <bob.capella@gmail.com>
# SPDX-License-Identifier: MIT

"""Calculate the entraining CAPE (ECAPE) of a parcel"""
import metpy.calc as mpcalc
import numpy as np
import pint
from metpy.constants import dry_air_spec_heat_press, earth_gravity
from metpy.units import units


def lfc_height(pressure, height, temperature, dew_point_temperature):
    lfc_p, lfc_t = mpcalc.lfc(pressure, temperature, dew_point_temperature)
    lfc_idx = np.abs(pressure - lfc_p).argmin()
    lfc_z = height[lfc_idx]
    return lfc_idx, lfc_z


def el_height(pressure, height, temperature, dew_point_temperature):
    el_p, el_t = mpcalc.el(pressure, temperature, dew_point_temperature)
    el_idx = np.abs(pressure - el_p).argmin()
    el_z = height[el_idx]
    return el_idx, el_z


def sr_wind_1km_bunkers_right(pressure, u_wind, v_wind, height):
    bunkers_right = mpcalc.bunkers_storm_motion(pressure, u_wind, v_wind, height)[0]
    u_sr = u_wind - bunkers_right[0]
    v_sr = v_wind - bunkers_right[1]

    avg_u_sr_1km = np.mean(u_sr[np.nonzero(height <= 1000 * units("m"))])
    avg_v_sr_1km = np.mean(v_sr[np.nonzero(height <= 1000 * units("m"))])

    sr_wind = mpcalc.wind_speed(avg_u_sr_1km, avg_v_sr_1km)

    return sr_wind


def calc_mse(pressure, height, temperature, specific_humidity):
    moist_static_energy = mpcalc.moist_static_energy(height, temperature, specific_humidity)
    saturation_mixing_ratio = mpcalc.saturation_mixing_ratio(pressure, temperature)
    moist_static_energy_star = mpcalc.moist_static_energy(height, temperature, saturation_mixing_ratio)

    moist_static_energy_bar = np.cumsum(moist_static_energy) / np.arange(1, len(moist_static_energy) + 1)
    moist_static_energy_bar = moist_static_energy_bar.to("J/kg")
    moist_static_energy_star = moist_static_energy_star.to("J/kg")
    return moist_static_energy_bar, moist_static_energy_star


def calc_integral_arg(moist_static_energy_bar, moist_static_energy_star, temperature):
    integral_arg = -(earth_gravity / (dry_air_spec_heat_press * temperature)) * (
        moist_static_energy_bar - moist_static_energy_star
    )
    return integral_arg


def calc_ncape(integral_arg, height, lfc_idx, el_idx):

    ncape = np.sum(
        (0.5 * integral_arg[lfc_idx:el_idx] + 0.5 * integral_arg[lfc_idx + 1 : el_idx + 1])
        * (height[lfc_idx + 1 : el_idx + 1] - height[lfc_idx:el_idx])
    )

    return ncape


def calc_ecape_a(sr_wind, pitchfork, ncape, cape):
    term_a = sr_wind**2 / 2.0
    term_b = (-1 - pitchfork - (2 * pitchfork / sr_wind**2) * ncape) / (4 * pitchfork / sr_wind**2)
    term_c = (
        np.sqrt(
            (1 + pitchfork + (2 * pitchfork / sr_wind**2) * ncape) ** 2
            + 8 * (pitchfork / sr_wind**2) * (cape - (pitchfork * ncape))
        )
    ) / (4 * pitchfork / sr_wind**2)
    ecape_a = term_a + term_b + term_c
    return ecape_a


def calc_pitchfork(el_z):
    # additional constants
    sigma = 1.6 * units("dimensionless")
    alpha = 0.8 * units("dimensionless")
    l_mix = 120.0 * units("m")
    pr = (1.0 / 3.0) * units("dimensionless")  # prandtl number
    ksq = 0.18 * units("dimensionless")  # von karman constant
    pitchfork = (ksq * alpha**2 * np.pi**2 * l_mix) / (4.0 * pr * sigma**2 * el_z)
    return pitchfork


def calc_ecape(
    height: pint.Quantity,
    pressure: pint.Quantity,
    temperature: pint.Quantity,
    specific_humidity: pint.Quantity,
    u_wind: pint.Quantity,
    v_wind: pint.Quantity,
    cape_type: str = "most unstable",
) -> pint.Quantity:
    """
    Calculate the entraining CAPE (ECAPE) of a parcel

    Parameters:
    ----------
        height : 'pint.Quantity'
            Atmospheric heights at the levels given by 'pressure'.
        pressure : 'pint.Quantity'
            Total atmospheric pressure
        temperature : 'pint.Quantity'
            Air temperature
        specific humidity : 'pint.Quantity'
            Specific humidity
        u_wind : 'pint.Quantity'
            X component of the wind
        v_wind : 'pint.Quantity'
            Y component of the wind
        cape_type: str
            Variation of CAPE desired. 'most unstable' (default), 'surface based', or 'mixed layer'

    Returns:
    ----------
        ecape : 'pint.Quantity'
            Entraining CAPE
    """

    cape_func = {
        "most unstable": mpcalc.most_unstable_cape_cin,
        "surface based": mpcalc.surface_based_cape_cin,
        "mixed layer": mpcalc.mixed_layer_cape_cin,
    }

    # calculate cape
    dew_point_temperature = mpcalc.dewpoint_from_specific_humidity(pressure, temperature, specific_humidity)
    cape, _ = cape_func[cape_type](pressure, temperature, dew_point_temperature)

    # calculate the level of free convection (lfc) and equilibrium level (el) indexes
    lfc_idx, _ = lfc_height(pressure, height, temperature, dew_point_temperature)
    el_idx, el_z = el_height(pressure, height, temperature, dew_point_temperature)

    # calculate the buoyancy dilution potential (ncape)
    moist_static_energy_bar, moist_static_energy_star = calc_mse(pressure, height, temperature, specific_humidity)
    integral_arg = calc_integral_arg(moist_static_energy_bar, moist_static_energy_star, temperature)
    ncape = calc_ncape(integral_arg, height, lfc_idx, el_idx)

    # calculate the storm relative (sr) wind
    sr_wind = sr_wind_1km_bunkers_right(pressure, u_wind, v_wind, height)

    # calculate the entraining cape (ecape)
    pitchfork = calc_pitchfork(el_z)
    ecape_a = calc_ecape_a(sr_wind, pitchfork, ncape, cape)

    return ecape_a
