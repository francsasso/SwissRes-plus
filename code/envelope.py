"""Thermal-node representation of building envelopes and ground contact.

Assign solar absorption, areal heat capacities and inter-node conductances to
walls, roofs, floors and windows. Ground-contact calculations provide monthly
boundary temperatures and an equivalent thermal resistance.
"""

import numpy as np
import pandas as pd
from constants import typology_elements


def solar_absorption_of_element(building_object):
    """Return node solar-absorptance factors with shape (5, 10).

    Copy building_object.solar_abs_elements into the first/external node row;
    other rows are zero. The element count follows the module typology array."""
    el_list = len(typology_elements)
    solar_abs_elements = building_object.solar_abs_elements
    node_solar_absorptance = np.zeros((5, el_list))
    node_solar_absorptance[0, :] = solar_abs_elements
    return node_solar_absorptance


def temp_calculation_of_ground(
    archetype_object, weather, lambda_gr=2.0, R_si=0.17, R_se=0.04, psi_k=0.05
):
    """Calculate the equivalent ground resistance and monthly temperatures.

    Parameters
    ----------
    archetype_object : archetype_object
        Envelope areas, setpoints, wall thickness, floor resistance and bridges.
    weather : pandas.DataFrame
        DatetimeIndex and temperature [degC], including the warmup period.
    lambda_gr : float, optional
        Ground conductivity [W/(m K)].
    R_si, R_se : float, optional
        Internal and external surface resistances [m2 K/W].
    psi_k : float, optional
        Ground perimeter thermal-bridge coefficient [W/(m K)].

    Returns
    -------
    r_gr_ve : float
        Equivalent ground-layer resistance [m2 K/W].
    theta_gr_ve : numpy.ndarray
        Monthly ground boundary temperatures [degC], January through December.
    thermal_bridge_heat : float
        Whole-building thermal-bridge coefficient [W/K].

    Notes
    -----
    Floor geometry is derived from south/east facade and ground areas. Zero wall
    thickness and floor resistance are assigned 0.35 m and 5.3 m2 K/W in place.
    The floor resistance is an independent ground-model parameter. A perimeter
    bridge is assigned when the area-based value is zero and no explicit override
    exists. Climate statistics include every supplied row, including warmup.
    """
    R_gr = 0.5 / lambda_gr
    try:
        external_temperature_monthly_averages = weather["temperature"].resample("ME").mean()
        external_temperature_monthly_minima = weather["temperature"].resample("ME").min()
        external_temperature_monthly_maxima = weather["temperature"].resample("ME").max()
    except:
        external_temperature_monthly_averages = weather["temperature"].resample("M").mean()
        external_temperature_monthly_minima = weather["temperature"].resample("M").min()
        external_temperature_monthly_maxima = weather["temperature"].resample("M").max()
    amplitude_of_external_temperature_variations = (
        external_temperature_monthly_maxima - external_temperature_monthly_minima
    ).mean() / 2
    annual_mean_external_temperature = external_temperature_monthly_averages.mean()
    heating_setpoint = archetype_object.heating_setpoint
    cooling_setpoint = archetype_object.cooling_setpoint
    annual_mean_internal_temperature = (heating_setpoint + cooling_setpoint) / 2
    amplitude_of_internal_temperature_variations = (cooling_setpoint - heating_setpoint) / 2
    coldest_month = 1
    months = np.arange(1, 13)
    internal_temperature_by_month = (
        annual_mean_internal_temperature
        - amplitude_of_internal_temperature_variations
        * np.cos(2 * np.pi * (months - coldest_month) / 12)
    )
    area = archetype_object.area_elements
    sog_area = area[4]
    south_facade_area = area[1] + area[7]
    east_facade_area = area[2] + area[8]
    facade_height = np.sqrt(east_facade_area * south_facade_area / sog_area)
    sog_width = south_facade_area / facade_height
    sog_length = sog_area / sog_width
    exposed_perimeter = 2 * (sog_length + sog_width)
    characteristic_floor_dimension = sog_area / (0.5 * exposed_perimeter)
    if archetype_object.wall_thickness == 0:
        archetype_object.__setattr__("wall_thickness", 0.35)
    if archetype_object.thermal_resistance_floor == 0:
        archetype_object.__setattr__("thermal_resistance_floor", 5.3)
    equivalent_ground_thickness = archetype_object.wall_thickness + lambda_gr * (
        archetype_object.thermal_resistance_floor + R_se
    )
    if equivalent_ground_thickness < characteristic_floor_dimension:
        U_sog = (
            2
            * lambda_gr
            / (np.pi * characteristic_floor_dimension + equivalent_ground_thickness)
            * np.log(np.pi * characteristic_floor_dimension / equivalent_ground_thickness + 1)
        )
    else:
        U_sog = lambda_gr / (0.457 * characteristic_floor_dimension + equivalent_ground_thickness)
    r_gr_ve = 1 / U_sog - R_si - archetype_object.thermal_resistance_floor - R_gr
    if archetype_object._thermal_bridge_heat is None and archetype_object.thermal_bridge_heat == 0:
        archetype_object.thermal_bridge_heat = exposed_perimeter * psi_k
    steady_state_heat_transfer_coefficient = sog_area * U_sog + exposed_perimeter * psi_k
    periodic_penetration_depth = 3.2
    h_pi = (
        sog_area
        * lambda_gr
        / equivalent_ground_thickness
        * np.sqrt(
            2
            / (np.float_power(1 + periodic_penetration_depth / equivalent_ground_thickness, 2) + 1)
        )
    )
    h_pe = (
        0.37
        * exposed_perimeter
        * lambda_gr
        * np.log(periodic_penetration_depth / equivalent_ground_thickness + 1)
    )
    annual_average_heat_flow_rate = steady_state_heat_transfer_coefficient * (
        annual_mean_internal_temperature - annual_mean_external_temperature
    )
    months = np.arange(1, 13)
    periodic_heat_flow_due_to_internal_temperature_variation = (
        -h_pi
        * amplitude_of_internal_temperature_variations
        * np.cos(2 * np.pi * (months - coldest_month) / 12)
    )
    periodic_heat_flow_due_to_external_temperature_variation = (
        h_pe
        * amplitude_of_external_temperature_variations
        * np.cos(2 * np.pi * (months - coldest_month - 1) / 12)
    )
    average_heat_flow_rate = (
        annual_average_heat_flow_rate
        + periodic_heat_flow_due_to_internal_temperature_variation
        + periodic_heat_flow_due_to_external_temperature_variation
    )
    theta_gr_ve = internal_temperature_by_month - (
        average_heat_flow_rate
        - exposed_perimeter
        * psi_k
        * (annual_mean_internal_temperature - annual_mean_external_temperature)
    ) / (sog_area * U_sog)
    thermal_bridge_heat = archetype_object.thermal_bridge_heat
    return (r_gr_ve, theta_gr_ve, thermal_bridge_heat)


def number_of_nodes_element(archetype_object):
    """Return total node count, per-element counts and cumulative offsets.

    Opaque and ground elements have five nodes; windows have two. Zone air is
    global node zero. The remaining nodes follow the element order in constants.py.
    Return (n_total_nodes, n_nodes_per_element, element_node_offsets).
    """
    el_list = len(archetype_object.typology_elements)
    n_nodes_per_element = np.full(el_list, 5)
    n_nodes_per_element[archetype_object.typology_elements == "W"] = 2
    element_node_offsets = np.array([0] * el_list)
    for eli in range(1, el_list):
        element_node_offsets[eli] = element_node_offsets[eli - 1] + n_nodes_per_element[eli - 1]
    n_total_nodes = element_node_offsets[-1] + n_nodes_per_element[-1] + 1
    return (n_total_nodes, n_nodes_per_element, element_node_offsets)


def conductance_node_of_element(archetype_object, lambda_gr=2.0):
    """Calculate inter-node conductances per unit area [W/(m2 K)].

    Return a (4, 10) array of links between consecutive element nodes. Subtract
    internal and external film resistances from 1/U and distribute the resulting
    construction resistance. lambda_gr is ground conductivity [W/(m K)].

    Write construction resistances back to the archetype. Nonfinite or nonpositive
    net resistance uses max(0.1 / U, 1e-6) [m2 K/W].
    """
    ground_layer_resistance = 0.5 / lambda_gr
    element_typology = archetype_object.typology_elements
    inter_node_conductance = np.zeros((4, len(element_typology)))
    u_eli = archetype_object.transmittance_u_elements
    # Update the archetype resistance array in place.
    element_construction_resistance = archetype_object.thermal_resistance_r_elements
    internal_convective_coefficient = archetype_object.heat_convective_elements_internal
    internal_radiative_coefficient = archetype_object.heat_radiative_elements_internal
    external_convective_coefficient = archetype_object.heat_convective_elements_external
    external_radiative_coefficient = archetype_object.heat_radiative_elements_external
    for i in range(0, len(element_typology)):
        r_raw = 1.0 / u_eli[i]
        r_int = 1.0 / (internal_convective_coefficient[i] + internal_radiative_coefficient[i])
        r_ext = 1.0 / (external_convective_coefficient[i] + external_radiative_coefficient[i])
        r_net = r_raw - r_int - r_ext
        if not np.isfinite(r_net) or r_net <= 0.0:
            r_net = max(r_raw * 0.1, 1e-06)
        element_construction_resistance[i] = r_net
    archetype_object.thermal_resistance_r_elements = element_construction_resistance
    layer_no = 0
    for i in range(len(element_typology)):
        if archetype_object.thermal_resistance_r_elements[i] != 0:
            if element_typology[i] == "OP":
                inter_node_conductance[0, i] = 6 / element_construction_resistance[i]
            elif element_typology[i] == "W":
                inter_node_conductance[0, i] = 1 / element_construction_resistance[i]
            elif element_typology[i] == "GR":
                inter_node_conductance[0, i] = 2 / ground_layer_resistance
    layer_no = 1
    for i in range(len(element_typology)):
        if archetype_object.thermal_resistance_r_elements[i] != 0:
            if element_typology[i] == "OP":
                inter_node_conductance[layer_no, i] = 3 / element_construction_resistance[i]
            elif element_typology[i] == "GR":
                inter_node_conductance[layer_no, i] = 1 / (
                    element_construction_resistance[i] / 4 + ground_layer_resistance / 2
                )
    layer_no = 2
    for i in range(len(element_typology)):
        if archetype_object.thermal_resistance_r_elements[i] != 0:
            if element_typology[i] == "OP":
                inter_node_conductance[layer_no, i] = 3 / element_construction_resistance[i]
            elif element_typology[i] == "GR":
                inter_node_conductance[layer_no, i] = 2 / element_construction_resistance[i]
    layer_no = 3
    for i in range(len(element_typology)):
        if archetype_object.thermal_resistance_r_elements[i] != 0:
            if element_typology[i] == "OP":
                inter_node_conductance[layer_no, i] = 6 / element_construction_resistance[i]
            elif element_typology[i] == "GR":
                inter_node_conductance[layer_no, i] = 4 / element_construction_resistance[i]
    return inter_node_conductance


def areal_heat_capacity_of_element(archetype_object):
    """Distribute element areal heat capacity among thermal nodes [J/(m2 K)].

    Return shape (5, 10). construction_class selects class_i, class_e, class_ie,
    class_d or class_m mass placement. Window capacities are zero. Ground nodes
    combine floor capacity with the prescribed soil contribution; other element
    capacities follow the building-age category.
    """
    element_typology = archetype_object.typology_elements
    element_areal_heat_capacity = archetype_object.thermal_capacity_elements
    node_areal_heat_capacity = np.zeros((5, len(element_typology)))
    if archetype_object.construction_class == "class_i":
        node = 1
        for i in range(len(element_typology)):
            if archetype_object.typology_elements[i] == "GR":
                node_areal_heat_capacity[node, i] = 1000000.0
        node = 4
        for i in range(len(element_typology)):
            if archetype_object.typology_elements[i] == "OP":
                node_areal_heat_capacity[node, i] = element_areal_heat_capacity[i]
            elif archetype_object.typology_elements[i] == "GR":
                node_areal_heat_capacity[node, i] = element_areal_heat_capacity[i]
    elif archetype_object.construction_class == "class_e":
        node = 0
        for i in range(len(element_typology)):
            if archetype_object.typology_elements[i] == "OP":
                node_areal_heat_capacity[node, i] = element_areal_heat_capacity[i]
            elif archetype_object.typology_elements[i] == "GR":
                node = 2
                node_areal_heat_capacity[node, i] = element_areal_heat_capacity[i]
    elif archetype_object.construction_class == "class_ie":
        node = 0
        for i in range(len(element_typology)):
            if archetype_object.typology_elements[i] == "OP":
                node_areal_heat_capacity[node, i] = element_areal_heat_capacity[i] / 2
            elif archetype_object.typology_elements[i] == "GR":
                node_areal_heat_capacity[node, i] = element_areal_heat_capacity[i] / 2
        node = 4
        for i in range(len(element_typology)):
            if archetype_object.typology_elements[i] == "OP":
                node_areal_heat_capacity[node, i] = element_areal_heat_capacity[i] / 2
            elif archetype_object.typology_elements[i] == "GR":
                node_areal_heat_capacity[node, i] = element_areal_heat_capacity[i] / 2
    elif archetype_object.construction_class == "class_d":
        node_list_1 = [1, 2, 3]
        for node in node_list_1:
            for i in range(len(element_typology)):
                if archetype_object.typology_elements[i] == "OP":
                    node_areal_heat_capacity[node, i] = element_areal_heat_capacity[i] / 4
                if archetype_object.typology_elements[i] == "GR":
                    if node == 2:
                        node_areal_heat_capacity[node, i] = element_areal_heat_capacity[i] / 4
                    if node == 3:
                        node_areal_heat_capacity[node, i] = element_areal_heat_capacity[i] / 2
        node_list_2 = [0, 4]
        for node in node_list_2:
            for i in range(len(element_typology)):
                if archetype_object.typology_elements[i] == "OP":
                    node_areal_heat_capacity[node, i] = element_areal_heat_capacity[i] / 8
                if archetype_object.typology_elements[i] == "GR":
                    if node == 4:
                        node_areal_heat_capacity[node, i] = element_areal_heat_capacity[i] / 4
    elif archetype_object.construction_class == "class_m":
        node = 2
        for i in range(len(element_typology)):
            if archetype_object.typology_elements[i] == "OP":
                node_areal_heat_capacity[node, i] = element_areal_heat_capacity[i]
            if archetype_object.typology_elements[i] == "GR":
                node = 3
                node_areal_heat_capacity[node, i] = element_areal_heat_capacity[i]
    return node_areal_heat_capacity
