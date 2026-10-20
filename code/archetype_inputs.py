"""Building-archetype inputs and assembly for SwissRes+.

Load preassembled CSV tables, validate their keys and columns, and construct
representative buildings with geometry, thermal properties and use profiles.
TABLE_COLUMNS and TABLE_KEYS define the supported table schemas. All element
arrays use the orientation order defined in constants.py.
"""

import numpy as np
import pandas as pd
from constants import (
    typology_elements,
    solar_abs_elements,
    heat_convective_elements_internal,
    heat_radiative_elements_internal,
    orientation_elements,
    heat_convective_elements_external,
    heat_radiative_elements_external,
    sky_factor_elements,
    thermal_capacity_elements_old,
    thermal_capacity_elements_new,
)

TABLE_COLUMNS = {
    "envelope_df": [
        "archetype_id",
        "buildingtype",
        "age",
        "typology",
        "era_specific",
        "g_value",
        "u_wall",
        "u_ground",
        "u_roof",
        "u_wind",
        "a_wall_nord",
        "a_wall_south",
        "a_wall_est",
        "a_wall_west",
        "a_ground",
        "a_roof",
        "a_window_nord",
        "a_window_south",
        "a_window_est",
        "a_window_west",
    ],
    "archetype_info_df": [
        "archetype_id",
        "buildingtype",
        "canton",
        "buildings",
        "era",
        "ventilation_presence",
        "retrofit",
    ],
    "archetype_prop_df": ["buildingtype", "Qs_Wp", "Occ_m2p", "Ea_Wm2", "El_Wm2"],
    "shading_df": ["buildingtype", "typology", "f_nord", "f_south", "f_est", "f_west"],
    "thermal_bridge_df": ["buildingtype", "renovation_period", "utb"],
    "temperature_correction": ["buildingtype", "age", "temperature_correction"],
    "schedule": [
        "buildingtype",
        "day",
        "hour",
        "occupancy",
        "appliances",
        "lighting",
        "heating_t_set",
        "cooling_t_set",
    ],
    "presence_month": ["buildingtype", "month", "factor"],
    "position": ["canton", "latitude", "longitude"],
    "weather": ["date", "canton", "ghi", "dhi", "dni", "temperature", "wspeed_avg"],
    "distribution_df": [
        "buildingtype",
        "share_gas",
        "share_oil",
        "share_wood",
        "share_dh",
        "share_hp",
    ],
    "efficiency_df": ["nec_gas", "nec_oil", "nec_wood", "nec_dh"],
}

TABLE_KEYS = {
    "envelope_df": ["buildingtype", "archetype_id"],
    "archetype_info_df": ["buildingtype", "archetype_id"],
    "archetype_prop_df": ["buildingtype"],
    "shading_df": ["buildingtype", "typology"],
    "thermal_bridge_df": ["buildingtype", "renovation_period"],
    "temperature_correction": ["buildingtype", "age"],
    "schedule": ["buildingtype", "day", "hour"],
    "presence_month": ["buildingtype", "month"],
    "position": ["canton"],
    "weather": ["canton", "date"],
    "distribution_df": ["buildingtype"],
    "efficiency_df": [],
}

TEXT_COLUMNS = {
    "archetype_id",
    "buildingtype",
    "typology",
    "canton",
    "day",
    "heating_t_set",
    "cooling_t_set",
}


def parse_csv_bool(value):
    """Parse a boolean from True/False or 1/0, including CSV strings.

    Raise ValueError for missing or ambiguous values.
    """
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, str):
        value = value.strip().lower()
        if value in ("true", "1"):
            return True
        if value in ("false", "0"):
            return False
    elif isinstance(value, (int, float, np.integer, np.floating)):
        if value in (0, 1):
            return bool(value)
    raise ValueError(f"A boolean True/False or 1/0 is required, got {value!r}")


def validate_table(name, df):
    """Check required columns, unique keys and valid values for a table.

    Call after numeric parsing, once per table. Return df without modifying it;
    invalid values and duplicate keys raise ValueError.
    """
    columns = TABLE_COLUMNS[name]
    if df[columns].isna().any().any():
        raise ValueError(f"{name}: required values cannot be missing")
    if TABLE_KEYS[name] and df.duplicated(TABLE_KEYS[name]).any():
        raise ValueError(f"{name}: duplicate keys {TABLE_KEYS[name]}")
    for col in TABLE_COLUMNS[name]:
        if col not in TEXT_COLUMNS | {"date", "ventilation_presence", "retrofit"}:
            if not np.isfinite(df[col].to_numpy()).all():
                raise ValueError(f"{name}.{col}: non-finite values")
    if name == "schedule":
        if not df["day"].isin(["weekday", "weekend"]).all():
            raise ValueError("schedule.day must be weekday/weekend")
        if not df["hour"].isin(range(24)).all():
            raise ValueError("schedule.hour must be an integer in 0..23")
        for col in ["heating_t_set", "cooling_t_set"]:
            if not df[col].isin(["setpoint", "setback", "off"]).all():
                raise ValueError(f"schedule.{col}: unknown control label")
    if name == "presence_month" and (not df["month"].isin(range(1, 13)).all()):
        raise ValueError("presence_month.month must be an integer in 1..12")
    return df


def load_csv_tables(paths, **read_csv_kwargs):
    """Load and validate the CSV tables specified by name.

    Parameters
    ----------
    paths : mapping
        Table names from TABLE_COLUMNS mapped to file paths or readable buffers.
    read_csv_kwargs
        Options passed to pandas.read_csv, such as delimiter, encoding or dtype.

    Returns
    -------
    dict
        Requested DataFrames keyed by table name. Boolean fields, numeric fields
        and weather dates are parsed according to their table schemas.
    """
    tables = {}
    for name, path in paths.items():
        columns = TABLE_COLUMNS[name]
        df = pd.read_csv(path, **read_csv_kwargs)
        if name == "archetype_info_df":
            for col in ["ventilation_presence", "retrofit"]:
                df[col] = df[col].map(parse_csv_bool)
        if name == "weather":
            df["date"] = pd.to_datetime(df["date"], errors="raise")
        for col in columns:
            if col not in TEXT_COLUMNS | {"date", "ventilation_presence", "retrofit"}:
                df[col] = pd.to_numeric(df[col], errors="raise")
        tables[name] = validate_table(name, df)
    return tables


class archetype_object(object):
    """Representative building geometry, thermal properties and control settings.

    Create a populated instance with archetype_table_from_df. era_specific is the
    representative energy reference area [m2]; number_buildings weights annual
    stock demand. Element arrays use the order defined in constants.py.

    thermal_bridge_heat is derived from utb and wall area unless assigned an
    explicit coefficient [W/K]. Assign None to restore the area-based value.
    Ground and node preparation update this object's thermal properties; use
    separate instances for independent simulations. Surface coefficient arrays
    are shared with constants.py.
    """

    def __init__(self, building_type: str, age: int):
        """Initialize metadata, shared coefficient arrays and unpopulated geometry.

        building_type selects control/ventilation rules; age selects capacity
        and infiltration categories. Populate remaining inputs before simulation."""
        self.archetype_id = 0
        self.typology = ""
        self.building_type = building_type
        self.age = age
        self.era = 0.0
        self.era_specific = 0.0
        self.canton = 0.0
        self.number_buildings = 0.0
        self.volume = 0.0
        self.ventilation = False
        self.heating_setpoint = 21.0
        self.cooling_setpoint = 26.0
        self.occupancy = 0.0
        self.internal_gains_occupants = 0.0
        self.internal_gains_appliance = 0.0
        self.internal_gains_light = 0.0
        self.wall_thickness = 0.0
        self.thermal_resistance_floor = 0.0
        self.area_elements = []
        self.transmittance_u_elements = []
        self.thermal_resistance_r_elements = []
        self.heat_convective_elements_internal = heat_convective_elements_internal
        self.heat_radiative_elements_internal = heat_radiative_elements_internal
        self.heat_convective_elements_external = heat_convective_elements_external
        self.heat_radiative_elements_external = heat_radiative_elements_external
        self.sky_factor_elements = sky_factor_elements
        self.solar_abs_elements = solar_abs_elements
        self.g_factor = []
        self.glass_share = 0.75
        self._thermal_bridge_heat = None
        self.retrofit = False
        self.construction_class = "class_i"
        self.typology_elements = typology_elements
        self.orientation_elements = orientation_elements
        self.shading = []
        self.cooling_setback = 26.0
        self.utb = 0.0

    @property
    def heating_setback(self):
        """Return heating setback [degC]: residential setpoint-2, otherwise 16."""
        if self.building_type in ["SFH", "MFH"]:
            return self.heating_setpoint - 2
        else:
            return 16.0

    @property
    def air_flow_rate_day(self):
        """Return specific day ventilation 36/occupancy [m3/(m2 h)], hours 07..20."""
        if self.building_type in ["SFH", "MFH"]:
            return 36 / self.occupancy
        elif self.building_type == "School":
            return 36 / self.occupancy
        else:
            return 36 / self.occupancy

    @property
    def air_flow_rate_night(self):
        """Return specific night ventilation 15/occupancy for SFH/MFH, else zero."""
        if self.building_type in ["SFH", "MFH"]:
            return 15 / self.occupancy
        else:
            return 0.0

    @property
    def thermal_capacity_elements(self):
        """Return shared element capacities [J/(m2 K)], split at construction year 1980."""
        if self.age < 1980:
            return thermal_capacity_elements_old
        else:
            return thermal_capacity_elements_new

    @property
    def air_change_rate_extra(self):
        """Return specific infiltration flow [m3/(m2 h)] by age and retrofit status.

        Retrofit buildings use 0.125. Other buildings use 0.4 before 1951,
        0.3 for 1951--1990, 0.25 for 1991--2000 and 0.125 from 2001 onward.
        """
        if self.retrofit:
            return 0.125
        if self.age < 1951:
            return 0.4
        elif self.age < 1991:
            return 0.3
        elif self.age < 2001:
            return 0.25
        else:
            return 0.125

    @property
    def power_heating_max(self):
        """Return maximum heating power [W]: 100 W/m2 times representative ERA."""
        return 100 * self.era_specific

    @property
    def power_cooling_max(self):
        """Return maximum cooling power [W]: -100 W/m2 times representative ERA."""
        return -100 * self.era_specific

    @property
    def thermal_bridge_heat(self):
        """Return an explicit bridge coefficient or utb times opaque-wall area [W/K].

        None means no override. An explicit zero is distinct from no override."""
        if self._thermal_bridge_heat is not None:
            return self._thermal_bridge_heat
        return self.utb * sum(self.area_elements[:4])

    @thermal_bridge_heat.setter
    def thermal_bridge_heat(self, value):
        """Set a finite bridge coefficient [W/K], or None to restore area derivation."""
        if value is None:
            self._thermal_bridge_heat = None
            return
        value = float(value)
        if not np.isfinite(value):
            raise ValueError("thermal_bridge_heat must be finite.")
        self._thermal_bridge_heat = value


def archetype_table_from_df(
    b_type,
    n_row,
    envelope_df,
    archetype_info_df,
    archetype_prop_df,
    shading_df,
    thermal_bridge_df,
    temperature_correction,
    *,
    glass_share=0.75,
):
    """Construct one representative building from six archetype input tables.

    Parameters
    ----------
    b_type : str
        Building-type key shared by the input tables.
    n_row : int
        Zero-based row in archetype_info_df after filtering by buildingtype.
    envelope_df : pandas.DataFrame
        Geometry and fabric properties keyed by buildingtype/archetype_id:
        age, typology, era_specific [m2], g_value [-], u_wall/u_ground/u_roof/
        u_wind [W/(m2 K)], and a_wall_*/a_window_*/a_ground/a_roof [m2].
        Orientation suffixes are nord, south, est and west.
    archetype_info_df : pandas.DataFrame
        buildingtype/archetype_id, canton, buildings, era, ventilation_presence
        and retrofit. buildings weights annual stock totals; era is metadata.
    archetype_prop_df : pandas.DataFrame
        buildingtype, Qs_Wp [W/person], Occ_m2p [m2/person], Ea_Wm2 and El_Wm2
        [W/m2] for occupants, appliances and lighting.
    shading_df : pandas.DataFrame
        buildingtype/typology and dimensionless f_nord/f_south/f_est/f_west.
    thermal_bridge_df : pandas.DataFrame
        buildingtype/renovation_period and utb [W/(m2 K)] per opaque wall area.
        renovation_period matches the envelope age.
    temperature_correction : pandas.DataFrame
        buildingtype/age and temperature_correction [K].
    glass_share : float, optional
        Glazed fraction of window area, default 0.75.

    Returns
    -------
    archetype_object
        Populated building with nominal g_value assigned to all four windows.

    Notes
    -----
    Existing-building wall U includes utb; the construction U is u_wall - utb.
    Retrofit wall U describes the construction alone and is used directly.
    Heating setpoint is 23 degC for retrofit, 19 degC for SFH age 1919, and
    21 degC plus the tabulated correction otherwise.
    """
    envelope = envelope_df[envelope_df["buildingtype"] == b_type].reset_index(drop=True)
    archetype_info = archetype_info_df[archetype_info_df["buildingtype"] == b_type].reset_index(
        drop=True
    )
    archetype_info_row = archetype_info.iloc[n_row].squeeze()
    archetype_id = archetype_info_row["archetype_id"]
    envelope_archetype = envelope[envelope["archetype_id"] == archetype_id]
    if len(envelope_archetype) != 1:
        raise ValueError("Expected one envelope row for the selected buildingtype/archetype_id.")
    envelope_archetype = envelope_archetype.iloc[0]
    archetype = archetype_object(building_type=b_type, age=envelope_archetype["age"])
    archetype.archetype_id = archetype_id
    archetype.typology = envelope_archetype["typology"]
    archetype.era_specific = envelope_archetype["era_specific"]
    g_value = float(envelope_archetype["g_value"])
    if not np.isfinite(g_value) or not 0 <= g_value <= 1:
        raise ValueError("g_value must be finite and between 0 and 1.")
    archetype.g_factor = np.array([0.0] * 6 + [g_value] * 4, dtype=np.float64)
    glass_share = float(glass_share)
    if not np.isfinite(glass_share) or not 0 < glass_share <= 1:
        raise ValueError("glass_share must be greater than 0 and at most 1.")
    archetype.glass_share = glass_share
    archetype.area_elements = np.array(
        [
            envelope_archetype["a_wall_nord"],
            envelope_archetype["a_wall_south"],
            envelope_archetype["a_wall_est"],
            envelope_archetype["a_wall_west"],
            envelope_archetype["a_ground"],
            envelope_archetype["a_roof"],
            envelope_archetype["a_window_nord"],
            envelope_archetype["a_window_south"],
            envelope_archetype["a_window_est"],
            envelope_archetype["a_window_west"],
        ],
        dtype=np.float64,
    )
    archetype.canton = archetype_info_row["canton"]
    archetype.number_buildings = archetype_info_row["buildings"]
    archetype.era = archetype_info_row["era"]
    archetype.ventilation = parse_csv_bool(archetype_info_row["ventilation_presence"])
    archetype.retrofit = parse_csv_bool(archetype_info_row["retrofit"])
    prop = archetype_prop_df[archetype_prop_df["buildingtype"] == b_type]
    if len(prop) != 1:
        raise ValueError("Expected one property row for the selected buildingtype.")
    prop = prop.iloc[0]
    archetype.internal_gains_occupants = prop["Qs_Wp"] / prop["Occ_m2p"]
    archetype.internal_gains_appliance = prop["Ea_Wm2"]
    archetype.internal_gains_light = prop["El_Wm2"]
    archetype.occupancy = prop["Occ_m2p"]
    shading = shading_df[
        (shading_df["buildingtype"] == b_type)
        & (shading_df["typology"].values == archetype.typology)
    ]
    if len(shading) != 1:
        raise ValueError("Expected one shading row for the selected buildingtype/typology.")
    shading = shading.iloc[0]
    archetype.shading = np.array(
        [
            shading["f_nord"],
            shading["f_south"],
            shading["f_est"],
            shading["f_west"],
            0,
            1,
            shading["f_nord"],
            shading["f_south"],
            shading["f_est"],
            shading["f_west"],
        ],
        dtype=np.float64,
    )
    tb = thermal_bridge_df[
        (thermal_bridge_df["buildingtype"] == b_type)
        & (thermal_bridge_df["renovation_period"] == archetype.age)
    ]
    if len(tb) != 1:
        raise ValueError("Expected one thermal-bridge row for the selected buildingtype/age.")
    archetype.utb = tb["utb"].iloc[0]
    # Separate the bridge addition only for existing, non-retrofitted walls.
    u_wall = float(envelope_archetype["u_wall"])
    if not archetype.retrofit:
        u_wall -= archetype.utb
    if not np.isfinite(u_wall) or u_wall <= 0:
        raise ValueError("The separated wall U-value must be finite and positive.")
    archetype.transmittance_u_elements = np.array(
        [u_wall] * 4
        + [envelope_archetype["u_ground"], envelope_archetype["u_roof"]]
        + [envelope_archetype["u_wind"]] * 4,
        dtype=np.float64,
    )
    archetype.thermal_resistance_r_elements = 1 / archetype.transmittance_u_elements
    temperature_correction = temperature_correction.loc[
        (temperature_correction["buildingtype"] == b_type)
        & (temperature_correction["age"] == archetype.age)
    ]
    if archetype.retrofit:
        archetype.heating_setpoint = 23.0
    elif b_type == "SFH" and archetype.age == 1919:
        archetype.heating_setpoint = 19.0
    else:
        if len(temperature_correction) != 1:
            raise ValueError(
                "Expected one temperature-correction row for the selected buildingtype/age."
            )
        archetype.heating_setpoint = 21.0 + temperature_correction["temperature_correction"].iloc[0]
    return archetype
