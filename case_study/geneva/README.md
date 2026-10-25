# Geneva case study

This case study provides a cantonal subset of SwissRes+ input data
for Geneva, Switzerland. It illustrates the structure, format and
parameters required by the model, covering residential and
non-residential building archetypes and hourly weather data.

The files complement the [calculation modules](../../code/) and
[model guide](../../code/MODEL_GUIDE.md). No case-specific execution
script is included.

## Building archetypes

The dataset contains 109 archetypes covering four building categories:

| Building category | Input identifier | Number of archetypes |
|---|---|---:|
| Offices | Administration | 27 |
| Multi-family houses | MFH | 30 |
| Single-family houses | SFH | 30 |
| Schools | School | 22 |

Each archetype describes a representative building and the number of
buildings it represents. Geometry and thermal properties are supplied
as preassembled inputs.

The underlying RegBL and CECB records and the procedure used to derive
these archetypes are not included.

## Input files

The CSV files are located directly in this directory.

| File | Contents |
|---|---|
| `archetype_info_df.csv` | Archetype identifiers, building categories, canton, stock weights and retrofit status |
| `envelope_df.csv` | Representative building geometry, surface areas, thermal transmittances and glazing properties |
| `archetype_prop_df.csv` | Building-use parameters, including area per occupant and specific internal gains |
| `shading_df.csv` | Orientation-specific shading factors by building category and urban context |
| `thermal_bridge_df.csv` | Thermal-bridge coefficients by building category and construction period |
| `temperature_correction.csv` | Heating-setpoint corrections by building category and construction period |
| `schedule.csv` | Weekday and weekend profiles for occupancy, appliances, lighting and temperature controls |
| `presence_month.csv` | Monthly occupancy adjustment factors |
| `position.csv` | Latitude and longitude used for the Geneva solar calculations |
| `weather_GE_2020.csv` | Hourly outdoor temperature, solar irradiance and wind speed |

Some CSV columns are retained from the research input tables but are
not used by the selected calculation workflow. The model guide and
function docstrings describe how the inputs enter the calculation.

## Weather and calculation period

The weather file contains 8,784 consecutive hourly records for the
calendar year 2020, including 29 February. The `canton` column is `GE`.

The calculation prepends December as a warmup period and excludes
these additional 744 hours from the reported results. Input timestamps
must follow the time convention described in the model guide; the
solar calculation uses fixed UTC+1.

## Scope

The Geneva dataset documents the types of input data required by
SwissRes+. It is not intended as a standalone simulation example
or as a reproduction of the national analysis and scenarios
presented in the associated publication.

Complete research results will be deposited on Yareta, the University
of Geneva research data repository. A link will be added once available.

## Associated publication

Sasso, F., & Patel, M. K. (2026).
*The effect of climate change, building retrofit and cooling system
adoption on future thermal energy demand in the Swiss building stock.*
Energy & Buildings, 363, 117544.

https://doi.org/10.1016/j.enbuild.2026.117544
