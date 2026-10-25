# SwissRes+ thermal simulation modules

SwissRes+ calculates hourly operative temperatures and heating/cooling demand
for representative buildings. The envelope is represented by thermal nodes,
with separate zone-air, surface, ventilation and ground-contact balances.

## Module overview

| Module | Purpose |
|---|---|
| `archetype_inputs.py` | Load CSV tables, check their schemas and assemble representative buildings |
| `constants.py` | Define element order, surface properties and solar reference coefficients |
| `hourly_inputs.py` | Prepare weather, irradiance, use profiles, controls, internal gains and ventilation |
| `envelope.py` | Build node layouts, conductances, capacities and ground boundary conditions |
| `thermal_simulation.py` | Solve hourly heat balances and calculate thermal demand and downstream energy use |

The modules use NumPy, pandas and Numba. Non-residential optimal-start processing
uses pandas 2.x APIs. Keep the Python files in the same directory and import them
from the research script or notebook.

## Assemble and simulate a building

```python
from archetype_inputs import archetype_table_from_df
from hourly_inputs import validate_hourly_inputs
from thermal_simulation.py import temperature_and_energy_needs_calculation_numba

building = archetype_table_from_df(
    b_type, n_row,
    envelope_df, archetype_info_df, archetype_prop_df,
    shading_df, thermal_bridge_df, temperature_correction,
)

validate_hourly_inputs(weather, building, position, schedule, presence_month)

hourly, annual = temperature_and_energy_needs_calculation_numba(
    position, weather, building, presence_month, schedule=schedule,
)
```

The six archetype tables describe geometry, fabric, use, shading, bridges and
temperature corrections. `n_row` is the zero-based position within the information
table filtered by building type. Other rows are matched by their documented keys.
`TABLE_COLUMNS` and `TABLE_KEYS` in `archetype_inputs.py` specify the CSV schemas.

CSV tables are parsed and checked once during loading. The builder selects each
matching row once and checks its uniqueness. For DataFrames already loaded in
memory, `validate_table(name, frame)` can check each parsed table before assembly.

Call `validate_hourly_inputs` once for each weather year, canton and building
type, before a loop of simulations sharing those inputs. The solver does not
repeat the full-year checks. Profile merges enforce unique day/hour matches.

`load_csv_tables` loads only the requested tables:

```python
from archetype_inputs import load_csv_tables

tables = load_csv_tables({
    "envelope_df": "envelope_df.csv",
    "archetype_info_df": "archetype_info_df.csv",
    "archetype_prop_df": "archetype_prop_df.csv",
    "shading_df": "shading_df.csv",
    "thermal_bridge_df": "thermal_bridge_df.csv",
    "temperature_correction": "temperature_correction.csv",
})
building = archetype_table_from_df(b_type, n_row, **tables)
```

Use a separate building instance for each independent run. Ground and node
preparation update its thermal properties. Surface coefficient arrays are shared
with `constants.py`.

## Hourly inputs

- `weather`: one complete hourly calendar year with a `date` column, `ghi`, `dhi`,
  `dni`, `temperature` and `wspeed_avg`, filtered to the building's canton.
  Irradiance is in W/m2 and temperature in degC. A combined CSV also contains
  `canton`; `weather_file(path, building)` reads and filters such a file.
- `position`: `canton`, `latitude` and `longitude` in degrees.
- `schedule`: `buildingtype`, `day`, `hour`, `occupancy`, `appliances`, `lighting`,
  `heating_t_set` and `cooling_t_set`. Provide all 48 weekday/weekend-hour pairs
  for the building type. Controls use the labels `setpoint`, `setback` and `off`.
- `presence_month`: `buildingtype`, `month` and `factor`, or a mapping keyed by
  `(buildingtype, month)`. Factors scale occupant gains; missing pairs use 1.0.

The calculation uses one-hour steps and prepends December as a 744-hour warmup.
Input dates retain their clock time when timezone labels are removed. Solar
geometry uses fixed UTC+1. Supply a consistent hourly time basis, including for
leap years.

## Physical quantities and conventions

Element arrays contain opaque walls N/S/E/W, ground floor, roof and windows
N/S/E/W. Areas refer to the representative building, in m2. U-values are in
W/(m2 K), node capacities in J/(m2 K) and inter-node conductances in W/(m2 K).

Existing wall input U includes `utb`; construction U is obtained by subtraction.
Retrofit wall U describes the construction alone. Glazing gains use nominal
`g_value` and the glazed fraction `glass_share`. Additional angle-dependent
columns from irradiance preparation are diagnostics.

Ventilation and infiltration are supplied as specific flows in m3/(m2 h) and
scaled by the representative reference area. Air specific heat capacity and
density default to 1006 J/(kg K) and 1.204 kg/m3. Set `c_air` and `rho_air` in
the simulation call or gains wrapper to use other values. Heating and cooling control is based
on operative temperature, the mean of air and mean radiant temperatures.

## Results

`hourly` has a timezone-naive date index after warmup removal:

| Column | Quantity | Unit |
|---|---|---|
| `Q_HC` | Signed thermal load: positive heating, negative cooling | W |
| `T_op` | Operative temperature | degC |
| `T_ext` | Outdoor temperature | degC |
| `Q_H`, `Q_C` | Heating and cooling magnitudes, truncated to integers | W |

`annual` contains weighted stock heating/cooling demand in MWh and representative
specific demand in kWh/m2. `number_buildings` weights stock totals; `era_specific`
is the representative reference area. Annual sums use the integer hourly loads.

Electricity and final-energy calculations are explicit follow-up calls:

```python
from simulation_numba import (
    add_electric_heating_cooling_demand,
    add_final_energy_to_annual_results,
)

annual, hourly = add_electric_heating_cooling_demand(hourly, annual, building)
annual = add_final_energy_to_annual_results(
    building, annual, distribution_df, efficiency_df,
)
```

The electricity helper returns annual results first. Final heating energy uses
system shares from `distribution_df` and efficiencies from `efficiency_df`.
