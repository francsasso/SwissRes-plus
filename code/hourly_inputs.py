"""Hourly boundary conditions and building-use profiles for SwissRes+.

Prepare weather, surface irradiance, occupancy, appliances, lighting, thermal
controls, internal gains and ventilation. Schedules are supplied explicitly.
December warmup provides the initial thermal history for annual simulations.
"""

import numpy as np
import pandas as pd
from numba import njit
from constants import _ISO52016_AOI_DEG, _ISO52016_F, PEREZ_COEFFICIENTS
from archetype_inputs import validate_table


# Weather input
def weather_file(weather_path, archetype):
    """Read hourly weather for the archetype's canton from a CSV file.

    Parameters
    ----------
    weather_path : path-like
        CSV containing date, canton, ghi, dhi, dni, temperature and wspeed_avg.
    archetype : archetype_object
        Building whose canton selects the weather rows.

    Returns
    -------
    pandas.DataFrame
        Weather with date as a column, irradiance in W/m2 and temperature in
        degC. Numeric parsing converts invalid weather values to NaN, which
        validate_hourly_inputs rejects.
    """
    df = pd.read_csv(weather_path, parse_dates=["date"])
    df.columns = [c.strip() for c in df.columns]
    cols_keep = ["date", "canton", "ghi", "dhi", "dni", "temperature", "wspeed_avg"]
    weather = df.loc[df["canton"] == archetype.canton, cols_keep].drop(columns=["canton"])
    numeric_cols = ["ghi", "dhi", "dni", "temperature", "wspeed_avg"]
    weather[numeric_cols] = weather[numeric_cols].apply(pd.to_numeric, errors="coerce")
    return weather


# Occupancy, appliance and lighting profiles


def merge_hourly_profiles(df, archetype, *, schedule=None):
    """Merge weekday/weekend use profiles into an hourly weather frame.

    df has a DatetimeIndex. schedule supplies buildingtype, day, hour, occupancy,
    appliances, lighting and heating/cooling control labels. The merge selects
    archetype.building_type and matches day/hour. The input index and day/hour
    columns are updated in place. Return the merged DataFrame.
    """
    df.index = pd.to_datetime(df.index)
    df["day"] = np.where(df.index.dayofweek < 5, "weekday", "weekend")
    df["hour"] = df.index.hour
    schedule_filtered = schedule[schedule["buildingtype"] == archetype.building_type]
    occupancy = pd.merge(
        df, schedule_filtered, on=["day", "hour"], how="left", validate="many_to_one"
    )
    occupancy.index = df.index
    return occupancy


# Setpoints, setbacks and seasonal controls
def apply_temperature_setpoints(occupancy, archetype):
    """Convert control labels to temperatures and apply seasonal operation.

    Map setpoint/setback/off using the archetype's settings. Heating is off
    May 15--September 15 inclusive; cooling is off outside that interval and
    for School during July. Off is represented by -100 degC for heating and
    100 degC for cooling. Return (occupancy, summer_mask).
    """
    heating_map = {
        "setback": archetype.heating_setback,
        "setpoint": archetype.heating_setpoint,
        "off": -100,
    }
    cooling_map = {
        "setback": archetype.cooling_setback,
        "setpoint": archetype.cooling_setpoint,
        "off": 100,
    }
    occupancy["heating_t_set"] = occupancy["heating_t_set"].map(heating_map)
    occupancy["cooling_t_set"] = occupancy["cooling_t_set"].map(cooling_map)
    month = occupancy.index.month
    day = occupancy.index.day
    summer_mask = ((month > 5) | (month == 5) & (day >= 15)) & (
        (month < 9) | (month == 9) & (day <= 15)
    )
    occupancy.loc[summer_mask, "heating_t_set"] = -100
    occupancy.loc[~summer_mask, "cooling_t_set"] = 100
    july_mask = month == 7
    if archetype.building_type == "School":
        occupancy.loc[july_mask, "cooling_t_set"] = 100
    return (occupancy, summer_mask)


# Non-residential optimal start and residential ramp
def apply_non_residential_optimal_start(occupancy, archetype):
    """Advance non-residential weekday heating according to outdoor temperature.

    A three-hour rolling temperature sets the ramp duration before 08:00.
    Duration is rounded and bounded to one through four hours; Monday can add
    one hour. Update heating_t_set in place and return the profile DataFrame.
    """
    occupancy_start_hour = 8
    base_warmup = 1
    k_per_K = 0.25
    max_warmup = 4
    monday_extra = 1
    cool_filter_h = 3
    is_non_res = archetype.building_type not in ["SFH", "MFH"]
    if is_non_res:
        idx = occupancy.index
        Text = occupancy["temperature"].rolling(f"{cool_filter_h}h").mean().fillna(method="bfill")
        Tset = archetype.heating_setpoint
        dT = (Tset - Text).clip(lower=0)
        warmup_hours = (
            (base_warmup + k_per_K * dT)
            .round()
            .clip(lower=base_warmup, upper=max_warmup)
            .astype(int)
        )
        is_weekday = occupancy["day"].eq("weekday")
        is_start_hour = occupancy["hour"].eq(occupancy_start_hour)
        is_monday = idx.dayofweek == 0
        monday_boost_vec = np.where(is_monday & is_start_hour & is_weekday, monday_extra, 0)
        starts = occupancy.index[is_weekday & is_start_hour]
        for t0 in starts:
            w0 = int(warmup_hours.loc[t0]) + int(
                monday_boost_vec[is_weekday & is_start_hour & (idx == t0)].sum()
            )
            if w0 <= 0:
                continue
            w0 = min(w0, max_warmup + monday_extra)
            hours_back = pd.date_range(end=t0 - pd.Timedelta(hours=1), periods=w0, freq="H")
            mask_hours = occupancy.index.isin(hours_back)
            mask_apply = mask_hours & (occupancy["heating_t_set"] != -100)
            if mask_apply.any():
                steps = np.arange(1, w0 + 1, dtype=float)
                weights = steps / w0
                weights_series = pd.Series(weights, index=hours_back)
                weights_series = weights_series.reindex(occupancy.index).fillna(0.0)
                delta_h = archetype.heating_setpoint - archetype.heating_setback
                occupancy.loc[mask_apply, "heating_t_set"] = (
                    archetype.heating_setback + weights_series[mask_apply] * delta_h
                )
    return occupancy


def apply_physical_residential_heating_ramp(occupancy, archetype, summer_mask):
    """Apply morning and evening heating ramps to SFH/MFH profiles in place.

    occupancy contains hour and heating_t_set; summer_mask marks heating-off
    rows. At 05:00 use setback, 06:00 halfway to setpoint and 07:00 setpoint.
    At 22:00 use the halfway value and at 23:00 setback. Summer rows keep the
    heating-off marker. Return occupancy.
    """
    if archetype.building_type not in ("SFH", "MFH"):
        return occupancy
    heating_setpoint = float(archetype.heating_setpoint)
    heating_setback = float(archetype.heating_setback)
    delta_t = heating_setpoint - heating_setback
    hour = occupancy["hour"].to_numpy(dtype=np.int8, copy=False)
    summer_mask_array = np.asarray(summer_mask, dtype=bool)
    heating_active = ~summer_mask_array
    heating_t_set = occupancy["heating_t_set"].to_numpy(dtype=np.float64, copy=True)
    mask_05 = heating_active & (hour == 5)
    mask_06 = heating_active & (hour == 6)
    mask_07 = heating_active & (hour == 7)
    heating_t_set[mask_05] = heating_setback
    heating_t_set[mask_06] = heating_setback + delta_t / 2.0
    heating_t_set[mask_07] = heating_setpoint
    mask_22 = heating_active & (hour == 22)
    mask_23 = heating_active & (hour == 23)
    heating_t_set[mask_22] = heating_setback + 0.5 * delta_t
    heating_t_set[mask_23] = heating_setback
    heating_t_set[summer_mask_array] = -100.0
    occupancy["heating_t_set"] = heating_t_set
    return occupancy


# Ventilation supply temperature and profile assembly
def apply_supply_air_temperature(occupancy, archetype, summer_mask):
    """Calculate ventilation supply temperature with heat recovery and bypass.

    Mechanical ventilation uses recovery effectiveness 0.7 and a 0.3 K bypass
    deadband, with heating/cooling setpoints as exhaust-temperature proxies.
    Natural ventilation uses outdoor temperature. Add t_supply [degC] to
    occupancy and return the updated frame.
    """
    if archetype.ventilation:
        eta_rec = 0.7
        deadband = 0.3
        T_ext = occupancy["temperature"].to_numpy()
        T_hset = occupancy["heating_t_set"].to_numpy()
        T_cset = occupancy["cooling_t_set"].to_numpy()
        T_target = np.where(summer_mask, T_cset, T_hset)
        t_supply_hr = T_ext + eta_rec * (T_target - T_ext)
        bypass_winter = ~summer_mask & (T_ext > T_hset + deadband)
        bypass_summer = summer_mask & (T_ext < T_cset - deadband)
        bypass = bypass_winter | bypass_summer
        t_supply = np.where(bypass, T_ext, t_supply_hr)
        occupancy["t_supply"] = t_supply
    else:
        occupancy["t_supply"] = occupancy["temperature"]
    return occupancy


def occupancy_hourly(df, archetype, *, schedule=None):
    """Assemble hourly use profiles, setpoints and ventilation supply temperature.

    Parameters
    ----------
    df : pandas.DataFrame
        Hourly weather with a DatetimeIndex and temperature [degC]. Its index
        and day/hour columns are updated during profile preparation.
    archetype : archetype_object
        Building type, setpoint/setback temperatures and ventilation status.
    schedule : pandas.DataFrame
        Explicit profile table with buildingtype, day, hour, occupancy,
        appliances, lighting, heating_t_set and cooling_t_set. Supply all
        48 weekday/weekend-hour combinations. Controls use setpoint/setback/off.

    Returns
    -------
    pandas.DataFrame
        Weather and use profiles with numeric heating/cooling setpoints and
        t_supply [degC], including seasonal controls and heating ramps.
    """
    occupancy = merge_hourly_profiles(df, archetype, schedule=schedule)
    occupancy, summer_mask = apply_temperature_setpoints(occupancy, archetype)
    occupancy = apply_non_residential_optimal_start(occupancy, archetype)
    occupancy = apply_physical_residential_heating_ramp(
        occupancy=occupancy, archetype=archetype, summer_mask=summer_mask
    )
    return apply_supply_air_temperature(occupancy, archetype, summer_mask)


# Monthly presence factors
def apply_monthly_occupancy_factor(df, archetype, presence_month, *, schedule=None):
    """Apply monthly presence factors to the hourly occupant profile.

    presence_month is a buildingtype/month/factor DataFrame or a mapping keyed
    by (buildingtype, month). Missing combinations use 1.0. Return the prepared
    hourly frame with occupancy scaled; appliance, lighting and control profiles
    are independent of this scaling. Weather fields are updated by occupancy_hourly.
    """
    occupancy_final = occupancy_hourly(df, archetype, schedule=schedule)
    occupancy_final.index = pd.to_datetime(occupancy_final.index)
    occupancy_final["month"] = occupancy_final.index.month
    if isinstance(presence_month, pd.DataFrame):
        presence_dict = presence_month.set_index(["buildingtype", "month"])["factor"].to_dict()
    else:
        presence_dict = presence_month
    occupancy_final["key"] = list(zip(occupancy_final["buildingtype"], occupancy_final["month"]))
    # Presence scales occupants only; electrical profiles and controls stay as prepared.
    occupancy_final["occupancy"] *= occupancy_final["key"].map(presence_dict).fillna(1.0)
    occupancy_final.drop(columns=["month", "key"], inplace=True)
    return occupancy_final


# Internal gains and ventilation/infiltration
@njit
def compute_h_ve_and_phi_int(
    buildingtype,
    occupancy,
    appliances,
    lighting,
    hour,
    era_specific,
    air_flow_day,
    air_flow_night,
    air_change_rate_extra,
    gain_appliance,
    gain_light,
    gain_occupants,
    c_air=1006.0,
    rho_air=1.204,
):
    """Compute total airflow, ventilation conductance and internal gains.

    All hourly arrays share one length. occupancy/appliances/lighting are
    schedule multipliers; hour is 0..23; era_specific is representative ERA
    [m2]. air_flow_day/night and air_change_rate_extra are specific flows
    [m3/(m2 h)]. gain_appliance/light/occupants are specific gains [W/m2].
    c_air is J/(kg K), rho_air is kg/m3; divide hourly flow by 3600.

    Returns
    -------
    h_ve : numpy.ndarray
        Ventilation plus infiltration heat transfer coefficient [W/K].
    phi_int : numpy.ndarray
        Total occupant, appliance and lighting gains [W].
    air_flow_rate : numpy.ndarray
        Total ventilation plus infiltration flow [m3/h].

    Notes
    -----
    Day flow is used at hours 07..20 inclusive. Administration/School have
    zero intentional ventilation when occupancy equals zero; infiltration
    remains active. Positive fractional occupancy does not scale the airflow."""
    n_steps = len(occupancy)
    h_ve = np.zeros(n_steps)
    phi_int = np.zeros(n_steps)
    air_flow_rate = np.zeros(n_steps)
    for i in range(n_steps):
        if buildingtype in ["Administration", "School"] and occupancy[i] == 0:
            air_flow = 0
        elif 6 < hour[i] < 21:
            air_flow = air_flow_day[i] * era_specific
        else:
            air_flow = air_flow_night[i] * era_specific
        # Convert both intentional ventilation and infiltration to total m3/h.
        air_flow += air_change_rate_extra * era_specific
        air_flow_rate[i] = air_flow
        h_ve[i] = c_air * rho_air / 3600.0 * air_flow
        phi_int[i] = (
            gain_appliance * appliances[i]
            + gain_light * lighting[i]
            + gain_occupants * occupancy[i]
        ) * era_specific
    return (h_ve, phi_int, air_flow_rate)


def vent_heat_transf_coef_and_Int_gains(
    archetype_object, weather_file, presence_month, c_air=1006, rho_air=1.204, *, schedule=None
):
    """Calculate hourly ventilation conductance and internal heat gains.

    Parameters
    ----------
    archetype_object : archetype_object
        Building use, reference area, ventilation rates and specific heat gains.
    weather_file : pandas.DataFrame
        Hourly weather with DatetimeIndex and temperature [degC].
    presence_month : pandas.DataFrame or mapping
        Monthly occupancy factors keyed by buildingtype and month.
    c_air : float, optional
        Air specific heat capacity [J/(kg K)], default 1006.
    rho_air : float, optional
        Air density [kg/m3], default 1.204.
    schedule : pandas.DataFrame
        Explicit weekday/weekend use profiles and control labels.

    Returns
    -------
    h_ve : pandas.Series
        Ventilation and infiltration heat-transfer coefficient [W/K].
    phi_int : pandas.Series
        Occupant, appliance and lighting heat gains [W].
    sim_df : pandas.DataFrame
        Weather, use profiles, controls, air flow rate [m3/h], h_ve and phi_int.
    """
    sim_df = apply_monthly_occupancy_factor(
        weather_file, archetype_object, presence_month, schedule=schedule
    )
    sim_df.index = pd.to_datetime(sim_df.index)
    occupancy = sim_df["occupancy"].to_numpy()
    appliances = sim_df["appliances"].to_numpy()
    lighting = sim_df["lighting"].to_numpy()
    hour = sim_df["hour"].to_numpy()
    era_specific = archetype_object.era_specific
    air_flow_day = np.full(len(occupancy), archetype_object.air_flow_rate_day)
    air_flow_night = np.full(len(occupancy), archetype_object.air_flow_rate_night)
    air_change_rate_extra = archetype_object.air_change_rate_extra
    gain_appliance = archetype_object.internal_gains_appliance
    gain_light = archetype_object.internal_gains_light
    gain_occupants = archetype_object.internal_gains_occupants
    buildingtype = archetype_object.building_type
    h_ve, phi_int, air_flow_rate = compute_h_ve_and_phi_int(
        buildingtype,
        occupancy,
        appliances,
        lighting,
        hour,
        era_specific,
        air_flow_day,
        air_flow_night,
        air_change_rate_extra,
        gain_appliance,
        gain_light,
        gain_occupants,
        c_air=c_air,
        rho_air=rho_air,
    )
    sim_df["air flow rate"] = air_flow_rate
    sim_df["h_ve"] = h_ve
    sim_df["phi_int"] = phi_int
    return (sim_df["h_ve"], sim_df["phi_int"], sim_df)


# Solar geometry and glazing angular reduction
def Equation_of_time(day_of_year: pd.Series) -> pd.Series:
    """Calculate the equation-of-time correction [minutes] from calendar day.

    day_of_year is a numeric pandas Series. The correction uses piecewise
    seasonal expressions and returns a Series with the same index.
    """
    t_eq = 0.45 * (day_of_year - 359)
    mask = day_of_year < 21
    t_eq[mask] = 2.6 + 0.44 * day_of_year[mask]
    mask = (day_of_year >= 21) & (day_of_year < 136)
    t_eq[mask] = 5.2 + 9 * np.cos((day_of_year[mask] - 43) * 0.0357)
    mask = (day_of_year >= 136) & (day_of_year < 241)
    t_eq[mask] = 1.4 - 5 * np.cos((day_of_year[mask] - 135) * 0.0449)
    mask = (day_of_year >= 241) & (day_of_year < 336)
    t_eq[mask] = -6.3 - 10 * np.cos((day_of_year[mask] - 306) * 0.036)
    return t_eq


def Air_mass_calc(solar_altitude_angle: pd.Series) -> pd.Series:
    """Calculate relative air mass from solar altitude [radians].

    Return a dimensionless pandas Series. Altitudes below 10 degrees use the
    low-altitude expression; higher altitudes use the reciprocal sine.
    """
    solar_altitude_angle_deg = np.degrees(solar_altitude_angle)
    air_mass = 1 / (
        np.sin(solar_altitude_angle) + 0.15 * np.power(solar_altitude_angle_deg + 3.885, -1.253)
    )
    mask = solar_altitude_angle_deg >= 10
    air_mass[mask] = 1 / np.sin(solar_altitude_angle[mask])
    return air_mass


def iso52016_reduction(theta_deg_array: np.ndarray) -> np.ndarray:
    """Calculate the incidence-angle reduction factor for glazing solar gains.

    theta_deg_array contains angles in degrees. Interpolate the angular curve
    after clipping to 0--90 degrees. Return dimensionless factors of the same shape.
    """
    theta = np.clip(np.asarray(theta_deg_array, dtype=float), 0.0, 90.0)
    return np.interp(theta, _ISO52016_AOI_DEG, _ISO52016_F)


def Hour_angle_calc(solar_time: pd.Series) -> pd.Series:
    """Calculate solar hour angle [degrees] from solar time [hours].

    Use a 12.5-hour origin and wrap angles to the -180..180 degree interval.
    """
    hour_angle_deg = 180 / 12 * (12.5 - solar_time)
    hour_angle_deg[hour_angle_deg > 180] -= 360
    hour_angle_deg[hour_angle_deg < -180] += 360
    return hour_angle_deg


# December warmup and surface irradiance
def prepare_december_warmup(sim_df):
    """Prepend December weather and irradiance as the annual warmup period.

    Copy December rows to the preceding year, keeping their solar geometry and
    all columns. Profile preparation uses the shifted dates. Return the extended
    frame; the thermal solver removes its first 744 hours from reported results.
    """
    dec_mask = sim_df.index.month == 12
    december = sim_df.loc[dec_mask].copy()
    if not december.empty:
        december.index = december.index.map(lambda dt: dt.replace(year=dt.year - 1))
        sim_df = pd.concat([december, sim_df])
    return sim_df


def irradiation_direction_ISO_fast(
    archetype_object,
    sim_df,
    position,
    include_circumsolar=True,
    add_dynamic_g=True,
    add_transmitted_columns=True,
    glass_share=0.75,
):
    """Calculate oriented surface irradiance and prepare the December warmup.

    Parameters
    ----------
    archetype_object : archetype_object
        Canton, orientation shading factors and nominal glazing g-factors.
    sim_df : pandas.DataFrame
        Weather with date, ghi/dhi/dni [W/m2] and temperature [degC].
    position : pandas.DataFrame
        Canton coordinates in latitude/longitude [degrees].
    include_circumsolar : bool, optional
        Include circumsolar irradiance in the shaded direct component.
    add_dynamic_g : bool, optional
        Add angle-dependent g_NV_eff/g_SV_eff/g_EV_eff/g_WV_eff columns.
    add_transmitted_columns : bool, optional
        Add I_win_NV/I_win_SV/I_win_EV/I_win_WV [W/m2]. Requires add_dynamic_g=True.
    glass_share : float, optional
        Glazed fraction applied to I_win_* columns, default 0.75.

    Returns
    -------
    pandas.DataFrame
        Weather indexed by timezone-naive date with NV/SV/EV/WV/HOR irradiance
        [W/m2], optional glazing columns and a prepended December warmup.

    Notes
    -----
    Solar geometry uses fixed UTC+1, hour+1 and a 365-day orbital approximation.
    Timezone labels are removed without converting the clock time. Copied warmup
    rows keep their calculated irradiance. The thermal solver uses nominal g;
    the angle-dependent glazing columns are additional diagnostics.
    """
    sim_df = sim_df.copy()
    sim_df["date"] = pd.to_datetime(sim_df["date"])
    sim_df["date"] = pd.to_datetime(sim_df["date"]).dt.tz_localize(None)
    sim_df = sim_df.sort_values("date").set_index("date")
    dni = sim_df["dni"].to_numpy(dtype=float)
    dhi = sim_df["dhi"].to_numpy(dtype=float)
    ghi = sim_df["ghi"].to_numpy(dtype=float)
    T = len(sim_df)
    coords = position.loc[position["canton"] == archetype_object.canton]
    if coords.empty:
        raise ValueError(f"Canton '{archetype_object.canton}' not found in position data.")
    lat_deg = float(coords.iloc[0]["latitude"])
    lon_deg = float(coords.iloc[0]["longitude"])
    tz = 1
    doy = sim_df.index.dayofyear.to_numpy(int)
    hod = (sim_df.index.hour + 1).to_numpy(int)
    n_days = 365
    earth_orbit = 2 * np.pi * (doy / n_days)
    decl_deg = (
        0.33281
        - 22.984 * np.cos(earth_orbit)
        - 0.3499 * np.cos(2 * earth_orbit)
        - 0.1398 * np.cos(3 * earth_orbit)
        + 3.7872 * np.sin(earth_orbit)
        + 0.03205 * np.sin(2 * earth_orbit)
        + 0.07187 * np.sin(3 * earth_orbit)
    )
    decl = np.radians(decl_deg)
    t_eq = Equation_of_time(pd.Series(doy)).to_numpy(float)
    time_shift = tz - lon_deg / 15.0
    solar_time = hod - t_eq / 60.0 - time_shift
    ha_deg = Hour_angle_calc(pd.Series(solar_time)).to_numpy(float)
    ha = np.radians(ha_deg)
    lat = np.radians(lat_deg)
    sin_alt = np.sin(decl) * np.sin(lat) + np.cos(decl) * np.cos(lat) * np.cos(ha)
    sin_alt = np.clip(sin_alt, 0.0, None)
    alt = np.arcsin(sin_alt)
    zen = np.pi / 2 - alt
    cos_zen = np.cos(zen)
    I_ext = 1370.0 * (1.0 + 0.033 * np.cos(earth_orbit))
    b_perez = np.maximum(np.cos(np.radians(85.0)), cos_zen)
    K_eps = 1.104
    mask_dhi = dhi > 0
    clearness = np.full(T, 999.0)
    clearness[mask_dhi] = (
        (dhi[mask_dhi] + dni[mask_dhi]) / dhi[mask_dhi] + K_eps * alt[mask_dhi] ** 3
    ) / (1.0 + K_eps * alt[mask_dhi] ** 3)
    air_mass = Air_mass_calc(pd.Series(alt)).to_numpy(float)
    sky_brightness = np.zeros(T)
    sky_brightness[mask_dhi] = air_mass[mask_dhi] * dhi[mask_dhi] / I_ext[mask_dhi]
    Perez = PEREZ_COEFFICIENTS
    breaks = Perez[:, 0]
    idx = np.digitize(clearness, breaks, right=False)
    idx[idx >= len(breaks)] = len(breaks) - 1
    f11 = Perez[idx, 1]
    f12 = Perez[idx, 2]
    f13 = Perez[idx, 3]
    f21 = Perez[idx, 4]
    f22 = Perez[idx, 5]
    f23 = Perez[idx, 6]
    F1 = np.maximum(0.0, f11 + f12 * sky_brightness + f13 * zen)
    F2 = f21 + f22 * sky_brightness + f23 * zen
    beta = np.radians(np.array([90.0, 90.0, 90.0, 90.0]))[None, :]
    gamma = np.radians(np.array([180.0, 0.0, 90.0, -90.0]))[None, :]
    decl_ = decl[:, None]
    ha_ = ha[:, None]
    cos_beta = np.cos(beta)
    sin_beta = np.sin(beta)
    cos_gamma = np.cos(gamma)
    sin_gamma = np.sin(gamma)
    term = (
        np.sin(decl_) * np.sin(lat) * cos_beta
        - np.sin(decl_) * np.cos(lat) * sin_beta * cos_gamma
        + np.cos(decl_) * np.cos(lat) * cos_beta * np.cos(ha_)
        + np.cos(decl_) * np.sin(lat) * sin_beta * cos_gamma * np.cos(ha_)
        + np.cos(decl_) * sin_beta * sin_gamma * np.sin(ha_)
    )
    cos_inc = np.clip(term, 0.0, None)
    I_dir = dni[:, None] * cos_inc
    I_circ = dhi[:, None] * F1[:, None] * cos_inc / (b_perez[:, None] + 1e-12)
    I_dif = dhi[:, None] * (
        (1.0 - F1)[:, None] * (1.0 + cos_beta) / 2.0
        + F1[:, None] * cos_inc / (b_perez[:, None] + 1e-12)
        + F2[:, None] * sin_beta
    )
    I_dif_ground = (dhi + dni * np.sin(alt))[:, None] * 0.2 * (1.0 - cos_beta) / 2.0
    I_dif_tot = np.maximum(I_dif - I_circ + I_dif_ground, 0.0)
    if include_circumsolar:
        I_dir_tot = np.maximum(I_dir + I_circ, 0.0)
    else:
        I_dir_tot = np.maximum(I_dir, 0.0)
    shade = np.array(
        [
            archetype_object.shading[0],
            archetype_object.shading[1],
            archetype_object.shading[2],
            archetype_object.shading[3],
        ],
        dtype=float,
    )[None, :]
    POA = shade * I_dir_tot + I_dif_tot
    sim_df["NV"] = POA[:, 0]
    sim_df["SV"] = POA[:, 1]
    sim_df["EV"] = POA[:, 2]
    sim_df["WV"] = POA[:, 3]
    sim_df["HOR"] = ghi
    if add_dynamic_g:
        aoi_rad = np.arccos(np.clip(cos_inc, 0.0, 1.0))
        aoi_deg = np.degrees(aoi_rad)
        F = iso52016_reduction(aoi_deg)
        g_nom = np.asarray(archetype_object.g_factor, dtype=float)
        g_nom_by_ori = np.array([g_nom[6], g_nom[7], g_nom[8], g_nom[9]], dtype=float)
        g_eff = F * g_nom_by_ori[None, :]
        sim_df["g_NV_eff"] = g_eff[:, 0]
        sim_df["g_SV_eff"] = g_eff[:, 1]
        sim_df["g_EV_eff"] = g_eff[:, 2]
        sim_df["g_WV_eff"] = g_eff[:, 3]
    if add_transmitted_columns:
        I_win = POA * g_eff * glass_share
        sim_df["I_win_NV"] = I_win[:, 0]
        sim_df["I_win_SV"] = I_win[:, 1]
        sim_df["I_win_EV"] = I_win[:, 2]
        sim_df["I_win_WV"] = I_win[:, 3]
    sim_df = prepare_december_warmup(sim_df)
    fill_cols = [
        "NV",
        "SV",
        "EV",
        "WV",
        "HOR",
        "g_NV_eff",
        "g_SV_eff",
        "g_EV_eff",
        "g_WV_eff",
        "I_win_NV",
        "I_win_SV",
        "I_win_EV",
        "I_win_WV",
    ]
    fill_cols = [c for c in fill_cols if c in sim_df.columns]
    sim_df[fill_cols] = sim_df[fill_cols].fillna(0.0)
    return sim_df


# Technical input checks
def validate_hourly_inputs(weather, archetype, position, schedule, presence_month):
    """Check hourly inputs once before simulations sharing the same setup.

    Call after loading/parsing, once per weather year, canton and building type.
    Require unique consecutive hours from January 1 00:00 to December 31 23:00,
    finite weather values, unique position keys and all 48 weekday/weekend-hour
    schedule rows. Accept leap years. Clock times are checked without converting
    timezone labels. Raise ValueError for invalid inputs.
    """
    for name, table in [("schedule", schedule), ("position", position)]:
        validate_table(name, table)
    if isinstance(presence_month, pd.DataFrame):
        validate_table("presence_month", presence_month)
    needed = ["date", "ghi", "dhi", "dni", "temperature", "wspeed_avg"]
    dates = pd.DatetimeIndex(pd.to_datetime(weather["date"], errors="raise")).tz_localize(None)
    year = int(dates[0].year)
    expected = pd.date_range(f"{year}-01-01", f"{year + 1}-01-01", freq="h", inclusive="left")
    if not dates.sort_values().equals(expected):
        raise ValueError("Provide one complete calendar year with consecutive hourly timestamps")
    for col in needed[1:]:
        if not np.isfinite(weather[col].to_numpy()).all():
            raise ValueError(f"Weather {col} must be finite")
    if "canton" in weather and (not weather["canton"].eq(archetype.canton).all()):
        raise ValueError("Weather must already be filtered to the selected canton")
    chosen = schedule.loc[schedule["buildingtype"].eq(archetype.building_type)]
    required = {(day, hour) for day in ["weekday", "weekend"] for hour in range(24)}
    if set(zip(chosen["day"], chosen["hour"])) != required:
        raise ValueError("Supply all 48 weekday/weekend schedule rows for the building type")
