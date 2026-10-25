# Methodology

SwissRes+ combines an archetype-based representation of the Swiss building stock with an hourly dynamic thermal model. It calculates indoor operative temperatures and useful energy demand for space heating and cooling.

The associated [publication](https://doi.org/10.1016/j.enbuild.2026.117544) describes the model, input characterisation, research scenarios and validation. The governing equations and calculation procedures are documented in EN ISO 52016-1 and the accompanying CEN ISO/TR 52016-2.

## Building-stock representation

The building stock is represented by archetypes describing residential buildings, offices and schools. Each archetype combines representative geometry, envelope properties, building-use parameters and operating conditions.

The repository uses preassembled archetype tables. The derivation of archetypes from the underlying RegBL and CECB records is outside the scope of the public workflow.

Each representative building is simulated individually. Its annual demand is weighted by the number of buildings it represents to obtain stock-level totals.

## Calculation workflow

The calculation follows five stages:

1. **Archetype assembly:** combine geometry, thermal properties, shading factors and building-use parameters from the input tables.
2. **Hourly boundary conditions:** prepare weather, surface irradiance, occupancy, internal gains, ventilation and temperature-control profiles.
3. **Thermal-network preparation:** assign envelope nodes, conductances, heat capacities and ground-contact boundary conditions.
4. **Hourly simulation:** solve the thermal balances and determine heating or cooling loads according to operative-temperature controls and available system capacity.
5. **Aggregation:** remove the initialisation period and calculate annual representative-building and stock-weighted demand.

The thermal calculation uses the resistance-capacitance representation described in EN ISO 52016-1. Readers should refer to the standard and its accompanying Technical Report for the detailed equations, node definitions and solution procedures. Function docstrings identify the relevant references.

## Model assumptions and conventions

SwissRes+ combines the standard-based thermal formulation with model-specific parameterisations. These include archetype properties, occupancy and equipment profiles, ventilation assumptions, temperature settings, seasonal operation and heating-start rules.

The calculation uses one-hour time steps. A copy of December is prepended as a 744-hour initialisation period and excluded from reported results.

Heating and cooling control is based on indoor operative temperature. Hourly thermal loads describe the representative building; annual stock totals apply the corresponding building-count weight. Specific annual demand is expressed per unit of representative energy reference area.

Optional postprocessing converts useful thermal demand into heat-pump electricity use or final energy. These calculations use separate technology assumptions and are distinct from the building thermal balances.

## Geneva case study

The [Geneva case study](../case_study/geneva/) provides input tables
for 109 residential and non-residential archetypes and hourly weather
data for 2020.

Its purpose is to illustrate the structure, format and parameters
required by SwissRes+. No case-specific execution script is included,
and the dataset does not cover the full national analysis or all
research scenarios.

For input descriptions and sources, see [Input data](input_data.md).
For module responsibilities, function calls, units and output
conventions, see the [model guide](../code/MODEL_GUIDE.md).

## References

- Sasso, F., & Patel, M. K. (2026). *The effect of climate change, building retrofit and cooling system adoption on future thermal energy demand in the Swiss building stock.* Energy & Buildings, 363, 117544. https://doi.org/10.1016/j.enbuild.2026.117544

- EN ISO 52016-1:2017. *Energy performance of buildings — Energy needs for heating and cooling, internal temperatures and sensible and latent heat loads — Part 1: Calculation procedures.*

- CEN ISO/TR 52016-2. *Energy performance of buildings — Energy needs for heating and cooling, internal temperatures and sensible and latent heat loads — Part 2: Explanation and justification of ISO 52016-1 and ISO 52017-1.*
