"""Hourly thermal balances, operative-temperature control and energy demand.

The Numba solver uses nominal glazing g-factors. Separate postprocessing
functions calculate heat-pump electricity and heating-system final energy.
"""

import numpy as np
import pandas as pd
from numba import njit
from hourly_inputs import (
    irradiation_direction_ISO_fast,
    vent_heat_transf_coef_and_Int_gains,
)
from envelope import (
    solar_absorption_of_element,
    temp_calculation_of_ground,
    number_of_nodes_element,
    conductance_node_of_element,
    areal_heat_capacity_of_element,
)


@njit
def build_matrices_numba(
    n_total_nodes,
    n_nodes_per_element,
    element_node_offsets,
    area_elements,
    area_elements_tot,
    irradiance_values,
    Phi_sol_zi,
    node_solar_absorptance,
    node_areal_heat_capacity,
    inter_node_conductance,
    element_type,
    theta_ext,
    theta_sup,
    phi_int,
    r_gr_ve,
    theta_gr_ve,
    Ah_ci,
    heat_conv_int,
    heat_rad_int,
    heat_conv_ext,
    heat_rad_ext,
    sky_factor_elements,
    delta_Theta_er,
    time_step_seconds,
    thermal_bridge_heat,
    internal_heat_capacity,
    n_load_cases,
    Phi_HC_nd_calc,
    internal_conv_frac,
    solar_conv_frac,
    H_conv_frac,
    C_conv_frac,
    Theta_old,
    h_ve,
):
    """Assemble one time step of the zone-air and envelope thermal balances.

    Node layout follows n_total_nodes, n_nodes_per_element and element_node_offsets.
    element_type codes outdoor/ground exposure as 0/1; Theta_old holds previous
    temperatures [degC]. time_step_seconds is in s, areas in m2, node capacities
    in J/(m2 K), surface conductances in W/(m2 K), internal_heat_capacity in J/K.
    Air/ground temperatures are degC, delta_Theta_er K and irradiance W/m2.

    Whole-zone gains and trial loads Phi_HC_nd_calc are W (heating positive,
    cooling negative); h_ve and thermal_bridge_heat are W/K. Convective fractions
    and sky factors are dimensionless. Bridges use theta_ext; ventilation uses
    theta_sup. Backward Euler discretizes storage. Air balance is in W;
    envelope balances are in W/m2.

    Returns
    -------
    MatA : numpy.ndarray
        Coefficient matrix (n_total_nodes, n_total_nodes).
    VecB : numpy.ndarray
        Right-hand sides (n_total_nodes, 3); first n_load_cases columns active.

    References
    ----------
    ISO 52016-1:2017, 6.5.6.1, Equation (37): matrix system;
    6.5.6.2, Equation (38): zone-air balance; 6.5.6.3, Equations (39)--(41):
    surface and internal node balances. Node coefficients are prepared using 6.5.7.
    The array layout and linearized surface-coupling implementation are described
    in MODEL_FORMULATION.md.
    """
    n_elements = len(n_nodes_per_element)
    MatA = np.zeros((n_total_nodes, n_total_nodes), dtype=np.float64)
    VecB = np.zeros((n_total_nodes, 3), dtype=np.float64)
    heat_radiative_elements_internal_mn = 0.0
    for i in range(n_elements):
        heat_radiative_elements_internal_mn += area_elements[i] * heat_rad_int[i]
    heat_radiative_elements_internal_mn /= area_elements_tot
    # Zone-air balance [W]: storage, surface convection, bridges and ventilation.
    row_idx = 0
    # Backward-Euler history and prescribed gains form the right-hand side.
    rhs_contribution = (
        thermal_bridge_heat * theta_ext
        + theta_sup * h_ve
        + internal_conv_frac * phi_int
        + solar_conv_frac * Phi_sol_zi
        + internal_heat_capacity / time_step_seconds * Theta_old[row_idx]
    )
    # Add the convective share of each heating/cooling trial load to the air node.
    for load_case_idx in range(n_load_cases):
        HC_conv_frac = H_conv_frac if Phi_HC_nd_calc[load_case_idx] > 0 else C_conv_frac
        x_temp_val = rhs_contribution.item()
        phi_hc = float(Phi_HC_nd_calc[load_case_idx])
        VecB[row_idx, load_case_idx] += x_temp_val + HC_conv_frac * phi_hc
    MatA[row_idx, row_idx] += (
        internal_heat_capacity / time_step_seconds + Ah_ci + thermal_bridge_heat + h_ve
    )
    # Couple zone air to each element's innermost surface through A*h_c.
    for element_idx in range(n_elements):
        element_node_count = n_nodes_per_element[element_idx]
        col_idx = element_node_offsets[element_idx] + element_node_count
        MatA[row_idx, col_idx] -= area_elements[element_idx] * heat_conv_int[element_idx]
    row_idx = 0
    for element_idx in range(n_elements):
        # Envelope-node balances [W/m2]: areal storage and absorbed short-wave radiation.
        for node_idx in range(n_nodes_per_element[element_idx]):
            row_idx += 1
            col_idx = element_node_offsets[element_idx] + node_idx + 1
            rhs_contribution = (
                node_solar_absorptance[node_idx, element_idx] * irradiance_values[element_idx]
                + node_areal_heat_capacity[node_idx, element_idx]
                / time_step_seconds
                * Theta_old[row_idx]
            )
            for load_case_idx in range(n_load_cases):
                VecB[row_idx, load_case_idx] += rhs_contribution
            if node_idx == n_nodes_per_element[element_idx] - 1:
                # Distribute nonconvective internal, solar and HVAC gains uniformly per internal surface area.
                rhs_contribution = (1 - internal_conv_frac) * phi_int + (
                    1 - solar_conv_frac
                ) * Phi_sol_zi
                for load_case_idx in range(n_load_cases):
                    HC_conv_frac = H_conv_frac if Phi_HC_nd_calc[load_case_idx] > 0 else C_conv_frac
                    VecB[row_idx, load_case_idx] += (
                        rhs_contribution + (1 - HC_conv_frac) * Phi_HC_nd_calc[load_case_idx]
                    ) / area_elements_tot
            MatA[row_idx, col_idx] += (
                node_areal_heat_capacity[node_idx, element_idx] / time_step_seconds
            )
            # Adjacent-node conductance contributes a diagonal loss and an off-diagonal coupling.
            if node_idx > 0:
                MatA[row_idx, col_idx] += inter_node_conductance[node_idx - 1, element_idx]
                MatA[row_idx, col_idx - 1] -= inter_node_conductance[node_idx - 1, element_idx]
            if node_idx < n_nodes_per_element[element_idx] - 1:
                MatA[row_idx, col_idx] += inter_node_conductance[node_idx, element_idx]
                MatA[row_idx, col_idx + 1] -= inter_node_conductance[node_idx, element_idx]
            if node_idx == 0:
                # Outdoor surface: convective/radiative exchange with air plus the sky-temperature correction.
                if element_type[element_idx] == 0:
                    for load_case_idx in range(n_load_cases):
                        VecB[row_idx, load_case_idx] += (
                            heat_conv_ext[element_idx] + heat_rad_ext[element_idx]
                        ) * theta_ext - sky_factor_elements[element_idx] * heat_rad_ext[
                            element_idx
                        ] * delta_Theta_er
                    MatA[row_idx, col_idx] += heat_conv_ext[element_idx] + heat_rad_ext[element_idx]
                # Ground-contact surface: exchange with the monthly virtual-ground temperature.
                elif element_type[element_idx] == 1:
                    for load_case_idx in range(n_load_cases):
                        VecB[row_idx, load_case_idx] += 1 / r_gr_ve * theta_gr_ve
                    MatA[row_idx, col_idx] += 1 / r_gr_ve
            if node_idx == n_nodes_per_element[element_idx] - 1:
                MatA[row_idx, col_idx] += (
                    heat_conv_int[element_idx] + heat_radiative_elements_internal_mn
                )
                # Internal surface convection couples back to the zone-air temperature.
                MatA[row_idx, 0] -= heat_conv_int[element_idx]
                # Long-wave surface coupling uses area-weighted temperatures and radiative coefficients.
                for other_element_idx in range(n_elements):
                    other_surface_node_idx = n_nodes_per_element[other_element_idx] - 1
                    other_col_idx = (
                        1 + element_node_offsets[other_element_idx] + other_surface_node_idx
                    )
                    MatA[row_idx, other_col_idx] -= (
                        area_elements[other_element_idx]
                        / area_elements_tot
                        * heat_rad_int[other_element_idx]
                    )
    return (MatA, VecB)


def temperature_and_energy_needs_calculation_numba(
    position,
    weather_file,
    archetype_object,
    presence_month,
    n_load_cases=2,
    internal_areal_heat_capacity=10000,
    internal_conv_frac=0.4,
    solar_conv_frac=0.1,
    H_conv_frac=1,
    C_conv_frac=1,
    delta_Theta_er=11,
    *,
    schedule=None,
    c_air=1006,
    rho_air=1.204,
    **kwargs,
):
    """Simulate hourly operative temperature and heating/cooling demand.

    Parameters
    ----------
    position : pandas.DataFrame
        canton, latitude, longitude [degrees].
    weather_file : pandas.DataFrame
        Complete hourly calendar year: date, ghi/dhi/dni [W/m2],
        temperature [degC], wspeed_avg for the building's canton.
    archetype_object : archetype_object
        Populated building, modified by ground/node preparation.
        Use separate instances for independent runs.
    presence_month : pandas.DataFrame or mapping
        Occupancy factors keyed by buildingtype/month.
    n_load_cases : int, optional
        Simultaneous load cases; updated by control logic each time step.
    internal_areal_heat_capacity : float, optional
        Internal capacity per reference area [J/(m2 K)], default 10000.
    internal_conv_frac, solar_conv_frac, H_conv_frac, C_conv_frac : float, optional
        Convective fractions of internal, solar, heating and cooling gains.
    delta_Theta_er : float, optional
        Outdoor-air to sky temperature difference [K], default 11.
    schedule : pandas.DataFrame
        Weekday/weekend occupancy, appliance, lighting and control profiles.
    c_air, rho_air : float, optional
        Air heat capacity [J/(kg K)] and density [kg/m3]; defaults 1006, 1.204.
    **kwargs
        Accepted without affecting the calculation.

    Returns
    -------
    hourly_results : pandas.DataFrame
        Date-indexed, excluding 744 warmup hours. Q_HC [W]: heating positive,
        cooling negative. T_op/T_ext [degC]; Q_H/Q_C: nonnegative W,
        truncated to int64.
    annual_results_df : pandas.DataFrame
        Q_H_annual/Q_C_annual [MWh], weighted by number_buildings;
        Q_H_annual_per_sqm/Q_C_annual_per_sqm [kWh/m2], representative building.

    Notes
    -----
    Validate shared inputs once with hourly_inputs.validate_hourly_inputs.
    Control uses operative temperature: mean of air and mean radiant temperature.
    Initial node temperature is 20 degC; time step is 3600 s. Glazing gains use
    nominal g_factor and glass_share. Window irradiance is reduced to 10 percent
    above 300 W/m2 outside steps 3960..6912, above 500 W/m2 inside that interval.
    Step indices and ground climate statistics include warmup; annual aggregation
    uses integer hourly loads.

    References
    ----------
    CEN ISO/TR 52016-2, 6.5.5: load determination by interpolation.
    ISO 52016-1:2017, 6.5.5.3, Equations (35)--(36): operative and mean radiant temperatures;
    6.5.6: thermal balances; 6.5.11: zone internal capacity. The use profiles,
    seasonal controls, shading activation, initialization choice and integer-load
    reporting are SwissRes+ model conventions.
    """
    # Hourly boundary conditions: surface irradiance, warmup, use profiles and ventilation.
    weather = irradiation_direction_ISO_fast(
        archetype_object, weather_file, position, glass_share=archetype_object.glass_share
    )
    weather.index = weather.index.tz_localize(None)
    h_ve, phi_int, sim_df = vent_heat_transf_coef_and_Int_gains(
        archetype_object,
        weather,
        presence_month,
        c_air=c_air,
        rho_air=rho_air,
        schedule=schedule,
    )
    h_ve = h_ve.to_numpy()
    n_steps = len(sim_df)
    temperature_array = sim_df["temperature"].to_numpy()
    t_supply_array = sim_df["t_supply"].to_numpy()
    heating_t_set_array = sim_df["heating_t_set"].to_numpy()
    cooling_t_set_array = sim_df["cooling_t_set"].to_numpy()
    phi_int_array = phi_int.to_numpy()
    month_array = sim_df.index.month.to_numpy()
    # Wind-dependent exterior convection [W/(m2 K)] for the hourly outdoor balance.
    hc_ext = 4.0 + 4.0 * sim_df["wspeed_avg"].to_numpy()
    hc_ext_array = np.empty(10, dtype=np.float64)
    # Initialize trial loads, node temperatures and the hourly state arrays.
    Phi_HC_nd_calc = np.zeros(3, dtype=np.float64)
    Phi_HC_nd_act = np.zeros(n_steps, dtype=np.float64)
    Theta_int_air = np.zeros((n_steps, 3), dtype=np.float64)
    Theta_int_r_mn = np.zeros((n_steps, 3), dtype=np.float64)
    Theta_int_op = np.zeros((n_steps, 3), dtype=np.float64)
    Theta_op_act = np.zeros(n_steps, dtype=np.float64)
    time_step_seconds = 3600.0
    active_col_idx = 0
    n_elements = len(archetype_object.typology_elements)
    typology_elements = archetype_object.typology_elements
    element_type = n_elements * ["EXT"]
    element_type[np.where(typology_elements == "GR")[0][0]] = "GR"
    area_elements = np.array(archetype_object.area_elements)
    area_elements_tot = np.sum(area_elements)
    orientation_elements = np.array(archetype_object.orientation_elements)
    element_g_factor = np.asarray(archetype_object.g_factor, dtype=np.float64)
    # Combine nominal solar gain coefficients with glazed window area.
    tau_area = element_g_factor * area_elements * archetype_object.glass_share
    irradiance_matrix = sim_df[orientation_elements].to_numpy()
    irradiance_matrix_shaded = irradiance_matrix.copy()
    internal_heat_capacity = internal_areal_heat_capacity * archetype_object.era_specific
    Ah_ci = np.dot(area_elements, archetype_object.heat_convective_elements_internal)
    heat_radiative_elements_internal_mn = (
        np.dot(area_elements, archetype_object.heat_radiative_elements_internal) / area_elements_tot
    )
    n_total_nodes, n_nodes_per_element, element_node_offsets = number_of_nodes_element(
        archetype_object
    )
    final_node_indices = element_node_offsets + n_nodes_per_element
    Theta_old = 20 * np.ones(n_total_nodes)
    VecB = 20 * np.ones((n_total_nodes, 3))
    type_map = {"EXT": 0, "GR": 1}
    element_type_codes = np.array([type_map[t] for t in element_type])
    # Prepare ground boundary conditions and the envelope thermal network.
    r_gr_ve, theta_gr_ve, thermal_bridge_heat = temp_calculation_of_ground(
        archetype_object, weather
    )
    sim_df.index = pd.to_datetime(sim_df.index)
    inter_node_conductance = conductance_node_of_element(archetype_object)
    node_areal_heat_capacity = areal_heat_capacity_of_element(archetype_object)
    node_solar_absorptance = solar_absorption_of_element(archetype_object)
    # Advance one hour using the previously selected node-temperature solution.
    for step_idx in range(n_steps):
        Theta_H_set = heating_t_set_array[step_idx]
        Theta_C_set = cooling_t_set_array[step_idx]
        Theta_old = np.asarray(VecB[:, active_col_idx])
        hc_ext_array[:] = hc_ext[step_idx]
        irradiance_values = irradiance_matrix[step_idx]
        # Apply shading thresholds using step indices including warmup.
        if step_idx < 3960 or step_idx > 6912:
            mask = irradiance_matrix[step_idx] >= 300
        else:
            mask = irradiance_matrix[step_idx] >= 500
        irradiance_matrix_shaded[step_idx][mask] *= 0.1
        irradiance_shaded = irradiance_matrix_shaded[step_idx]
        # Transmitted solar gain [W] = shaded irradiance * nominal g * glazed window area.
        Phi_sol_components = irradiance_shaded * tau_area
        Phi_sol_zi = Phi_sol_components.sum()
        theta_ext = temperature_array[step_idx]
        theta_sup = t_supply_array[step_idx]
        phi_int_value = phi_int_array[step_idx]
        h_ve_value = h_ve[step_idx]
        theta_gr_ve_month = theta_gr_ve[month_array[step_idx] - 1]
        if Theta_H_set < -99:
            power_heating_max_act = 0
        else:
            power_heating_max_act = archetype_object.power_heating_max
        if Theta_C_set > 99:
            power_cooling_max_act = 0
        else:
            power_cooling_max_act = archetype_object.power_cooling_max
        # Start with zero load and configure the heating/cooling trials enabled by the setpoints.
        Phi_HC_nd_calc[0] = 0
        if power_heating_max_act == 0 and power_cooling_max_act == 0:
            n_load_cases = 1
        elif power_cooling_max_act == 0:
            H_col_idx = 1
            n_load_cases = 2
            Phi_HC_nd_calc[H_col_idx] = power_heating_max_act
        elif power_heating_max_act == 0:
            C_col_idx = 1
            n_load_cases = 2
            Phi_HC_nd_calc[C_col_idx] = power_cooling_max_act
        else:
            n_load_cases = 3
            H_col_idx = 1
            C_col_idx = 2
            Phi_HC_nd_calc[H_col_idx] = power_heating_max_act
            Phi_HC_nd_calc[C_col_idx] = power_cooling_max_act
        # Solve the trial cases, then repeat with an accepted interpolated load.
        iterate = True
        while iterate:
            iterate = False
            MatA, VecB = build_matrices_numba(
                n_total_nodes,
                n_nodes_per_element,
                element_node_offsets,
                area_elements,
                area_elements_tot,
                irradiance_values,
                Phi_sol_zi,
                node_solar_absorptance,
                node_areal_heat_capacity,
                inter_node_conductance,
                element_type_codes,
                theta_ext,
                theta_sup,
                phi_int_value,
                r_gr_ve,
                theta_gr_ve_month,
                Ah_ci,
                archetype_object.heat_convective_elements_internal,
                archetype_object.heat_radiative_elements_internal,
                hc_ext_array,
                archetype_object.heat_radiative_elements_external,
                archetype_object.sky_factor_elements,
                delta_Theta_er,
                time_step_seconds,
                archetype_object.thermal_bridge_heat,
                internal_heat_capacity,
                n_load_cases,
                Phi_HC_nd_calc,
                internal_conv_frac,
                solar_conv_frac,
                H_conv_frac,
                C_conv_frac,
                Theta_old,
                h_ve_value,
            )
            # Solve all active right-hand-side columns against the same coefficient matrix.
            theta = np.linalg.solve(MatA, VecB)
            VecB = theta
            Theta_int_air[step_idx, :] = VecB[0, :]
            Theta_int_r_mn[step_idx, :] = 0
            # Mean radiant temperature is represented by the area-weighted inner-surface mean.
            Theta_int_r_mn[step_idx, :] = (
                np.dot(area_elements, VecB[final_node_indices, :]) / area_elements_tot
            )
            Theta_int_op[step_idx, :] = (
                Theta_int_air[step_idx, :] + Theta_int_r_mn[step_idx, :]
            ) / 2
            # Interpolate load from the trial operative temperatures; enforce prescribed power limits.
            if n_load_cases > 1:
                if Theta_int_op[step_idx, 0] < Theta_H_set:
                    Theta_op_set = Theta_H_set
                    Phi_HC_nd_act[step_idx] = (
                        archetype_object.power_heating_max
                        * (Theta_op_set - Theta_int_op[step_idx, 0])
                        / (Theta_int_op[step_idx, H_col_idx] - Theta_int_op[step_idx, 0])
                    )
                    if Phi_HC_nd_act[step_idx] > archetype_object.power_heating_max:
                        Phi_HC_nd_act[step_idx] = archetype_object.power_heating_max
                        Theta_op_act[step_idx] = Theta_int_op[step_idx, H_col_idx]
                        active_col_idx = H_col_idx
                    else:
                        Phi_HC_nd_calc[0] = Phi_HC_nd_act[step_idx]
                        Theta_op_act[step_idx] = Theta_op_set
                        active_col_idx = 0
                        n_load_cases = 1
                        # Re-solve the node temperatures using the interpolated heating load.
                        iterate = True
                elif Theta_int_op[step_idx, 0] > Theta_C_set:
                    Theta_op_set = Theta_C_set
                    Phi_HC_nd_act[step_idx] = (
                        archetype_object.power_cooling_max
                        * (Theta_op_set - Theta_int_op[step_idx, 0])
                        / (Theta_int_op[step_idx, C_col_idx] - Theta_int_op[step_idx, 0])
                    )
                    if Phi_HC_nd_act[step_idx] < archetype_object.power_cooling_max:
                        Phi_HC_nd_act[step_idx] = archetype_object.power_cooling_max
                        Theta_op_act[step_idx] = Theta_int_op[step_idx, C_col_idx]
                        active_col_idx = C_col_idx
                    else:
                        Phi_HC_nd_calc[0] = Phi_HC_nd_act[step_idx]
                        Theta_op_act[step_idx] = Theta_op_set
                        active_col_idx = 0
                        n_load_cases = 1
                        # Re-solve the node temperatures using the interpolated cooling load.
                        iterate = True
                else:
                    Phi_HC_nd_act[step_idx] = 0
                    Theta_op_act[step_idx] = Theta_int_op[step_idx, 0]
                    active_col_idx = 0
            else:
                Phi_HC_nd_act[step_idx] = Phi_HC_nd_calc[0]
                Theta_op_act[step_idx] = Theta_int_op[step_idx, 0]
                active_col_idx = 0
    # Remove warmup before aggregating integer hourly loads.
    Tstep_first_act = 744
    hourly_results = pd.DataFrame(
        data=np.vstack(
            (
                Phi_HC_nd_act[Tstep_first_act:],
                Theta_op_act[Tstep_first_act:],
                sim_df["temperature"][Tstep_first_act:],
            )
        ).T,
        index=sim_df[Tstep_first_act:].index,
        columns=["Q_HC", "T_op", "T_ext"],
    )
    hourly_results["Q_H"] = 0
    mask = hourly_results["Q_HC"] > 0
    hourly_results.loc[mask, "Q_H"] = hourly_results.loc[mask, "Q_HC"].astype("int64")
    hourly_results["Q_C"] = 0
    mask = hourly_results["Q_HC"] < 0
    hourly_results.loc[mask, "Q_C"] = -hourly_results.loc[mask, "Q_HC"].astype("int64")
    # One-hour load sums are numerically Wh; stock totals use buildings as the aggregation weight.
    Q_H_annual = hourly_results["Q_H"].sum() * archetype_object.number_buildings / 1000000
    Q_C_annual = hourly_results["Q_C"].sum() * archetype_object.number_buildings / 1000000
    Q_H_annual_per_sqm = hourly_results["Q_H"].sum() / archetype_object.era_specific / 1000
    Q_C_annual_per_sqm = hourly_results["Q_C"].sum() / archetype_object.era_specific / 1000
    annual_results_dic = {
        "Q_H_annual": Q_H_annual,
        "Q_C_annual": Q_C_annual,
        "Q_H_annual_per_sqm": Q_H_annual_per_sqm,
        "Q_C_annual_per_sqm": Q_C_annual_per_sqm,
    }
    annual_results_df = pd.DataFrame([annual_results_dic])
    return (hourly_results, annual_results_df)


# Optional electricity and final-energy postprocessing
def calculate_cop(
    df: pd.DataFrame,
    archetype,
    efficiency_heating: float = 0.45,
    efficiency_cooling: float = 0.35,
    cop_max_heat: float = 10.0,
    cop_max_cool: float = 8.0,
    dT_cond_heat: float = 0.0,
    dT_evap_heat: float = 0.0,
    dT_cond_cool: float = 7.0,
    dT_evap_cool: float = 5.0,
    min_dt_K: float = 1.0,
):
    """Estimate hourly heat-pump COP from cycle temperatures and efficiencies.

    Parameters
    ----------
    df : pandas.DataFrame
        Hourly DatetimeIndex and T_ext [degC].
    archetype : archetype_object
        Heating supply: 35 degC for retrofit or age > 2000, otherwise 60 degC.
        Cooling supply: 10 degC.
    efficiency_heating, efficiency_cooling : float, optional
        Ideal-cycle efficiency fractions, defaults 0.45 and 0.35.
    cop_max_heat, cop_max_cool : float, optional
        COP upper bounds, defaults 10 and 8.
    dT_cond_heat, dT_evap_heat, dT_cond_cool, dT_evap_cool : float, optional
        Fluid-to-cycle approach temperature differences [K].
    min_dt_K : float, optional
        Lower bound on cycle temperature difference [K].

    Returns
    -------
    pandas.DataFrame
        Copy of df with dimensionless COP. Heating: September 15--May 15
        inclusive; cooling otherwise.

    References
    ----------
    This is an energy-system postprocessing model. ISO 52016-1:2017 provides
    the useful thermal loads; this function's cycle efficiencies, supply temperatures
    and COP bounds are separate heat-pump modelling assumptions.
    """
    if "T_ext" not in df.columns:
        raise ValueError("DataFrame must contain 'T_ext' with outdoor air temperature in °C.")
    if not hasattr(archetype, "retrofit") or not hasattr(archetype, "age"):
        raise AttributeError("archetype must expose 'retrofit' (bool) and 'age' (int).")
    out = df.copy()
    # Set heating and cooling supply temperatures [degC].
    heating_temp_C = 35.0 if bool(archetype.retrofit) or archetype.age > 2000 else 60.0
    cooling_temp_C = 10.0
    month = out.index.month
    day = out.index.day
    heating_mask = (
        (month < 5) | (month == 5) & (day <= 15) | (month == 9) & (day >= 15) | (month > 9)
    )
    T_air_K = out["T_ext"].astype(float) + 273.15
    T_heat_supply_K = heating_temp_C + 273.15
    T_cool_supply_K = cooling_temp_C + 273.15
    T_cond_heat_K = T_heat_supply_K + dT_cond_heat
    T_evap_heat_K = T_air_K - dT_evap_heat
    T_evap_cool_K = T_cool_supply_K - dT_evap_cool
    T_cond_cool_K = T_air_K + dT_cond_cool
    denom_heat = (T_cond_heat_K - T_evap_heat_K).clip(lower=min_dt_K)
    cop_heat_carnot = T_cond_heat_K / denom_heat
    denom_cool = (T_cond_cool_K - T_evap_cool_K).clip(lower=min_dt_K)
    cop_cool_carnot = T_evap_cool_K / denom_cool
    cop_heat_real = (efficiency_heating * cop_heat_carnot).clip(lower=0, upper=cop_max_heat)
    cop_cool_real = (efficiency_cooling * cop_cool_carnot).clip(lower=0, upper=cop_max_cool)
    out["COP"] = np.where(heating_mask, cop_heat_real, cop_cool_real)
    return out


def add_electric_heating_cooling_demand(hourly_results, annual_results_df, archetype):
    """Convert Q_H/Q_C [W] to heat-pump electricity with hourly COP.

    hourly_results requires a DatetimeIndex and T_ext [degC]; a copy receives
    COP and Q_H_el/Q_C_el [W]. Update annual_results_df with Q_H_el_annual/
    Q_C_el_annual [MWh], weighted by number_buildings.
    Return (annual_results_df, hourly_results).
    """
    hourly_results = calculate_cop(hourly_results, archetype)
    hourly_results["Q_H_el"] = np.where(
        hourly_results["COP"] > 0, hourly_results["Q_H"] / hourly_results["COP"], 0
    )
    hourly_results["Q_C_el"] = np.where(
        hourly_results["COP"] > 0, hourly_results["Q_C"] / hourly_results["COP"], 0
    )
    Q_H_el_annual = hourly_results["Q_H_el"].sum()
    Q_C_el_annual = hourly_results["Q_C_el"].sum()
    annual_results_df["Q_H_el_annual"] = Q_H_el_annual * archetype.number_buildings / 1000000
    annual_results_df["Q_C_el_annual"] = Q_C_el_annual * archetype.number_buildings / 1000000
    return (annual_results_df, hourly_results)


def add_final_energy_to_annual_results(
    archetype, annual_results_df, distribution_df, efficiency_df
):
    """Aggregate heating final energy by system share and efficiency.

    annual_results_df supplies Q_H_annual/Q_H_el_annual [MWh]. distribution_df
    is keyed by buildingtype, with fractional share_gas/share_oil/share_wood/
    share_dh/share_hp. efficiency_df supplies nec_gas/nec_oil/nec_wood/nec_dh
    from its first row. Divide useful demand by efficiency and weight by share;
    add the heat-pump share of heating electricity.

    Update and return annual_results_df with Q_H_annual_final [MWh] and,
    when a heat-pump share is defined, Q_H_el_real [MWh].
    """
    buildingtype = archetype.building_type
    energy_need = annual_results_df["Q_H_annual"].values[0]
    hp_energy_final = annual_results_df["Q_H_el_annual"].values[0]
    distribution_mapping = {
        "share_gas": "gas",
        "share_oil": "oil",
        "share_wood": "wood",
        "share_dh": "dh",
    }
    efficiency_mapping = {"gas": "nec_gas", "oil": "nec_oil", "wood": "nec_wood", "dh": "nec_dh"}
    distribution_df = distribution_df.set_index("buildingtype").apply(
        pd.to_numeric, errors="coerce"
    )
    efficiency_df = efficiency_df.apply(pd.to_numeric, errors="coerce")
    efficiencies = {
        system: efficiency_df[eff_col].values[0] for system, eff_col in efficiency_mapping.items()
    }
    final_energy_total = 0
    if buildingtype in distribution_df.index:
        for dist_col, system in distribution_mapping.items():
            if dist_col in distribution_df.columns and system in efficiencies:
                distribution = distribution_df.at[buildingtype, dist_col]
                efficiency = efficiencies[system]
                if pd.notna(distribution) and pd.notna(efficiency) and (efficiency != 0):
                    final_energy_total += energy_need * distribution / efficiency
        share_hp = distribution_df.at[buildingtype, "share_hp"]
        if pd.notna(share_hp):
            final_energy_total += hp_energy_final * share_hp
            annual_results_df["Q_H_el_real"] = hp_energy_final * share_hp
    annual_results_df["Q_H_annual_final"] = final_energy_total
    return annual_results_df
