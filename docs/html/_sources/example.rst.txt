.. _example:

example
--------------

For use in a script:

.. code-block:: python

    from ecape.calc import calc_ecape

    from pathlib import Path
    import numpy as np
    from metpy.units import units

    sounding_loc = Path("./sounding.txt")
    data = np.genfromtxt(sounding_loc, delimiter=",")

    height = data[:, 0] * units("m")
    pressure = data[:, 1] * units("Pa")
    temperature = data[:, 2] * units("K")
    specific_humidity = data[:, 3] * units("kg/kg")
    u_wind = data[:, 4] * units("m/s")
    v_wind = data[:, 5] * units("m/s")

    cape_type = 'most_unstable'

    ecape = calc_ecape(height, pressure, temperature, specific_humidity, u_wind, v_wind, cape_type)
    print(f"{cape_type} ECAPE: {ecape}")
