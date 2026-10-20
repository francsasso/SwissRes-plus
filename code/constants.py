"""Surface properties and reference coefficients for the SwissRes+ thermal model.

Element arrays follow this order: north, south, east and west opaque walls;
ground floor; roof; north, south, east and west windows. OP denotes opaque
construction, GR ground contact and W glazing. NV/SV/EV/WV/HOR identify the
surface-irradiance columns used throughout the model.

Heat-transfer coefficients are in W/(m2 K), areal heat capacities in J/(m2 K),
and solar absorptance and sky-view factors are dimensionless. Glazing angular
reduction and Perez coefficients support surface-irradiance calculations.
"""

import numpy as np

_ISO52016_AOI_DEG = np.array([0, 30, 45, 60, 70, 75, 80, 85, 90], dtype=float)


_ISO52016_F = np.array([1.0, 0.99, 0.97, 0.86, 0.65, 0.48, 0.29, 0.12, 0.0], dtype=float)


typology_elements = np.array(["OP", "OP", "OP", "OP", "GR", "OP", "W", "W", "W", "W"], dtype=object)


solar_abs_elements = np.array([0.6, 0.6, 0.6, 0.6, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0], dtype=np.float64)


heat_convective_elements_internal = np.array(
    [2.5, 2.5, 2.5, 2.5, 0.7, 5.0, 2.5, 2.5, 2.5, 2.5], dtype=np.float64
)


heat_radiative_elements_internal = np.array(
    [5.13, 5.13, 5.13, 5.13, 5.13, 5.13, 5.13, 5.13, 5.13, 5.13], dtype=np.float64
)


orientation_elements = np.array(
    ["NV", "SV", "EV", "WV", "HOR", "HOR", "NV", "SV", "EV", "WV"], dtype=object
)


heat_convective_elements_external = np.array(
    [20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0], dtype=np.float64
)


heat_radiative_elements_external = np.array(
    [4.68, 4.68, 4.68, 4.68, 0, 4.89, 4.63, 4.63, 4.63, 4.63], dtype=np.float64
)


sky_factor_elements = np.array([0.5, 0.5, 0.5, 0.5, 0.0, 1.0, 0.5, 0.5, 0.5, 0.5], dtype=np.float64)


thermal_capacity_elements_old = np.array(
    [425000, 425000, 425000, 425000, 400000, 400000, 0, 0, 0, 0], dtype=np.float64
)


thermal_capacity_elements_new = np.array(
    [230000, 230000, 230000, 230000, 400000, 400000, 0, 0, 0, 0], dtype=np.float64
)


PEREZ_COEFFICIENTS = np.array(
    [
        [1.065, -0.008, 0.588, -0.062, -0.06, 0.072, -0.022],
        [1.23, 0.13, 0.683, -0.151, -0.019, 0.066, -0.029],
        [1.5, 0.33, 0.487, -0.221, 0.055, -0.064, -0.026],
        [1.95, 0.568, 0.187, -0.295, 0.109, -0.152, -0.014],
        [2.8, 0.873, -0.392, -0.362, 0.226, -0.462, 0.001],
        [4.5, 1.132, -1.237, -0.412, 0.288, -0.823, 0.056],
        [6.2, 1.06, -1.6, -0.359, 0.264, -1.127, 0.131],
        [99999, 0.678, -0.327, -0.25, 0.156, -1.377, 0.251],
    ],
    dtype=float,
)
