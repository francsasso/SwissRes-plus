# Input data

This document describes the CSV inputs supplied with the [Geneva case study](../case_study/geneva/). The files contain preassembled building archetypes and supporting parameters for the SwissRes+ hourly calculation method.

The case study includes 109 archetypes: 27 offices (`Administration`), 30 multi-family houses (`MFH`), 30 single-family houses (`SFH`) and 22 schools (`School`). All represent the canton of Geneva (`GE`).

The underlying RegBL and CECB records and the procedure used to derive the archetypes are not included.

For the calculation sequence and output conventions, see the [model guide](../code/MODEL_GUIDE.md).

## Sources

| Input group | Source |
|---|---|
| Hourly weather | MeteoSwiss, Genève-Cointrin station, 2020 |
| Building-use parameters, including internal loads and use schedules | SIA 2024:2021 [2] |
| Archetype characterisation, shading factors, thermal bridges and temperature corrections | Sasso and Patel (2026) [1] and the previous studies cited therein |


## Files and table relationships

The CSV files are located directly in `case_study/geneva/`. Column names retain those used by the research model.

| File | Matching keys | Role |
|---|---|---|
| `archetype_info_df.csv` | `buildingtype`, `archetype_id` | Archetype metadata and stock weights |
| `envelope_df.csv` | `buildingtype`, `archetype_id` | Geometry and envelope properties |
| `archetype_prop_df.csv` | `buildingtype` | Building-use parameters |
| `shading_df.csv` | `buildingtype`, `typology` | Orientation-specific shading |
| `thermal_bridge_df.csv` | `buildingtype`, `renovation_period` | Thermal-bridge coefficients; period matches archetype `age` |
| `temperature_correction.csv` | `buildingtype`, `age` | Heating-setpoint corrections |
| `schedule.csv` | `buildingtype`, `day`, `hour` | Hourly weekday and weekend profiles |
| `presence_month.csv` | `buildingtype`, `month` | Monthly occupancy factors |
| `position.csv` | `canton` | Coordinates for solar calculations |
| `weather_GE_2020.csv` | `canton`, `date` | Hourly outdoor conditions |

The archetype builder selects a row of `archetype_info_df` within a building category and matches its envelope by identifier. Supporting properties are matched by category, construction period or urban context.

Hourly loads describe one representative building. The number of buildings represented by each archetype is used separately to aggregate annual stock demand.

## 1. Archetype information

### `archetype_info_df.csv`

| Column | Meaning | Unit or values |
|---|---|---|
| `archetype_id` | Archetype identifier | Identifier |
| `buildingtype` | Building category | `Administration`, `MFH`, `SFH`, `School` |
| `age` | Construction-period label used by the model | Year label |
| `typology` | Urban-context category | `Rural`, `Suburban`, `Urban` |
| `canton` | Canton identifier | `GE` |
| `ventilation_presence` | Mechanical ventilation flag | Boolean |
| `buildings` | Number of buildings represented; annual stock aggregation weight | Buildings |
| `era_specific` | Energy reference area of the representative building | m² |
| `era` | Total energy reference area represented by the archetype | m² |
| `retrofit` | Retrofit flag used to select thermal and control assumptions | Boolean |
| `n_dh`, `n_electro`, `n_oil`, `n_hp`, `n_gas`, `n_wood` | Building counts by heating system: district heating, electric resistance, oil, heat pump, gas and wood | Buildings |

In the supplied dataset:

`era = era_specific × buildings`

Construction-period labels are category keys, not necessarily the exact construction year of an individual building.

Heating-system counts are retained as stock metadata. The selected thermal solver does not use them to calculate useful heating and cooling loads.

The builder takes the representative geometry and its construction-period label from the envelope table.

## 2. Envelope properties

### `envelope_df.csv`

| Column or group | Meaning | Unit |
|---|---|---|
| `archetype_id`, `buildingtype` | Keys linking the envelope to the archetype information | Identifier |
| `age`, `typology` | Construction-period and urban-context labels | Category |
| `era_specific` | Representative building energy reference area | m² |
| `a_wall_nord`, `a_wall_south`, `a_wall_est`, `a_wall_west` | Opaque wall areas by orientation | m² |
| `a_window_nord`, `a_window_south`, `a_window_est`, `a_window_west` | Window areas by orientation | m² |
| `a_roof`, `a_ground` | Roof and ground-contact floor areas | m² |
| `u_wall`, `u_wind`, `u_roof`, `u_ground` | Thermal transmittances of walls, windows, roof and ground floor | W/(m² K) |
| `g_value` | Nominal total solar energy transmittance of glazing | Dimensionless |

Orientation suffixes preserve the original table naming: `nord` = north, `south` = south, `est` = east and `west` = west.

The file also retains `a_wall_ext`, `a_wall_int`, `a_window`, `abs_wall`, `emi_wall`, `abs_roof`, `emi_roof` and `emi_wind`. These columns are not read by the selected DataFrame archetype builder. It uses the orientation-specific areas, while surface absorption and heat-transfer parameters are assigned through the model constants.

### Element order and thermal conventions

The model arrays follow this order:

1. North, south, east and west opaque walls.
2. Ground-contact floor and roof.
3. North, south, east and west windows.

For non-retrofitted walls, the input `u_wall` includes the thermal-bridge addition `utb`. The builder subtracts `utb` to obtain the wall construction U-value. For retrofitted walls, `u_wall` is used directly. Thermal bridges enter the zone balance separately.

The builder assigns `g_value` to the four window orientations. `glass_share`, defaulting to 0.75 in the code, is the glazed fraction of window area and is distinct from `g_value`.

The selected thermal solver uses nominal glazing g-factors. Angle-dependent glazing columns generated during solar preparation are diagnostic for this calculation path.

## 3. Building-use parameters

### `archetype_prop_df.csv`

| Column | Meaning | Unit |
|---|---|---|
| `buildingtype` | Building category | Category |
| `Occ_m2p` | Reference area per occupant | m²/person |
| `Qs_Wp` | Sensible heat gain per occupant | W/person |
| `Ea_Wm2` | Specific appliance heat gains | W/m² |
| `El_Wm2` | Specific lighting heat gains | W/m² |

The specific occupant gain is calculated as `Qs_Wp / Occ_m2p`. Hourly occupancy, appliance and lighting multipliers are then applied before multiplication by the representative reference area.

`X_ghp`, `ventilation_rate_m3_night` and `ventilation_rate_m3_day` are retained from the input tables but are not read by the selected builder. In this implementation, ventilation rates are defined by the archetype class rather than those CSV columns.

Building-use parameters are based on SIA 2024:2021 [2].

### `schedule.csv`

Each building category has 48 rows: 24 hours for `weekday` and 24 for `weekend`.

| Column | Meaning | Unit or values |
|---|---|---|
| `buildingtype` | Building category | Category |
| `day` | Day type | `weekday`, `weekend` |
| `hour` | Hour used to match the weather timestamp | 0–23 |
| `occupancy` | Occupant-presence multiplier | Dimensionless |
| `appliances` | Appliance-use multiplier | Dimensionless |
| `lighting` | Lighting-use multiplier | Dimensionless |
| `heating_t_set` | Heating control label | `setpoint`, `setback`, `off` |
| `cooling_t_set` | Cooling control label | `setpoint`, `setback`, `off` |
| `water` | Retained water-use profile; not used by the selected thermal workflow | Not used in this calculation |

Control labels are converted to temperatures using the archetype properties. Seasonal masks and heating-start rules subsequently modify these profiles; the CSV labels alone do not describe the final hourly controls.

Numerical off markers in the code are control conventions rather than physical indoor temperatures.

### `presence_month.csv`

| Column | Meaning | Unit or values |
|---|---|---|
| `buildingtype` | Building category | Category |
| `month` | Calendar month | 1–12 |
| `factor` | Monthly factor applied to the occupant-presence profile | Dimensionless |

The supplied table covers all 12 months for each building category.

Monthly factors directly scale occupancy only, not the appliance or lighting profiles. Occupancy can also affect ventilation through the non-residential zero-occupancy rule in the code.

## 4. Shading, thermal bridges and temperature corrections


### `shading_df.csv`

| Column | Meaning | Unit |
|---|---|---|
| `buildingtype`, `typology` | Building-category and urban-context keys | Category |
| `f_nord`, `f_south`, `f_est`, `f_west` | Orientation-specific shading factors | Dimensionless |

In the selected solar routine, these factors act on the direct/circumsolar component. Diffuse and ground-reflected contributions are added separately.

### `thermal_bridge_df.csv`

| Column | Meaning | Unit |
|---|---|---|
| `buildingtype` | Building category | Category |
| `renovation_period` | Period key matched to the envelope `age` | Year label |
| `utb` | Thermal-bridge addition per unit opaque wall area | W/(m² K) |

The default whole-building thermal-bridge coefficient is `utb` multiplied by the sum of the four opaque wall areas, giving W/K.

Ground-boundary preparation and explicit overrides are described in the function docstrings.

### `temperature_correction.csv`

| Column | Meaning | Unit |
|---|---|---|
| `buildingtype`, `age` | Building-category and construction-period keys | Category |
| `temperature_correction` | Correction added to the base heating setpoint | K |

The current builder uses 21 °C plus this correction for ordinary non-retrofitted archetypes. Its explicit exceptions are 23 °C for retrofitted buildings and 19 °C for SFH with period label 1919.

These rules are model assumptions, not a general prescription of SIA 2024.

## 5. Coordinates and weather

### `position.csv`

| Column | Meaning | Unit or values |
|---|---|---|
| `canton` | Key matching the archetype canton | `GE` |
| `latitude` | Latitude used for solar geometry | Degrees north |
| `longitude` | Longitude used for solar geometry | Degrees east |

### `weather_GE_2020.csv`


| Column | Meaning | Unit or format |
|---|---|---|
| `date` | Hourly timestamp | `YYYY-MM-DD HH:MM:SS` |
| `canton` | Canton identifier | `GE` |
| `ghi` | Global horizontal irradiance | W/m² |
| `dhi` | Diffuse horizontal irradiance | W/m² |
| `dni` | Direct normal irradiance | W/m² |
| `temperature` | Outdoor air temperature | °C |
| `wspeed_avg` | Wind speed | m/s |

The file contains 8,784 consecutive hourly records from 1 January to 31 December 2020, including 29 February.

The solar calculation uses fixed UTC+1. Timezone labels, if present, are removed without converting clock times. Timestamps must therefore be interpreted consistently with that convention.

December is copied internally to provide a 744-hour warmup preceding the reported year. No additional warmup rows are required in the input CSV. Warmup hours are excluded from reported hourly and annual results.

## Scope of the supplied inputs

The Geneva tables illustrate the structure and content of the model
inputs through a cantonal subset. They do not constitute a standalone
simulation example and do not include the full national input dataset,
the original energy-certificate records or all future climate and
building-stock scenarios used in the research.

Optional final-energy postprocessing requires additional heating-system share and efficiency tables (`distribution_df` and `efficiency_df`), which are not part of these ten case-study input files.

## References

[1] Sasso, F., & Patel, M. K. (2026). *The effect of climate change, building retrofit and cooling system adoption on future thermal energy demand in the Swiss building stock.* Energy & Buildings, 363, 117544. https://doi.org/10.1016/j.enbuild.2026.117544

[2] Swiss Society of Engineers and Architects (SIA). (2021). *Données d’utilisation des locaux pour l’énergie et les installations du bâtiment.* Norme SIA 2024.

[3] MeteoSwiss. Hourly weather data for 2020, Genève-Cointrin station. https://www.meteoswiss.admin.ch/
