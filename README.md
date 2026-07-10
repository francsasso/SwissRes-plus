# SwissRes+

**SwissRes+** is a dynamic bottom-up building stock model developed to analyse the evolution of space heating and cooling demand in the Swiss building sector under different climate, retrofit, and technology adoption scenarios. The model combines a detailed archetype-based representation of the Swiss building stock with an hourly dynamic thermal simulation framework based on **EN ISO 52016**. 

SwissRes+ represents residential and non-residential buildings through archetypes derived from the Swiss Building and Dwelling Register (RegBL), energy performance certificates (CECB), and national statistics. The model simulates thermal energy demand at hourly resolution, accounting for building characteristics, climate conditions, retrofit measures, heating technologies, and cooling system adoption. 

## Main Features

- Dynamic hourly simulation of space heating and cooling demand.
- Building physics based on the EN ISO 52016 resistance-capacitance (R-C) methodology.
- Detailed representation of the Swiss building stock through archetypes.
- Residential and non-residential sectors (single-family houses, multi-family houses, schools, and offices).
- Assessment of future climate scenarios (CH2018 RCP 2.6 and RCP 8.5).
- Simulation of building stock evolution, demographic growth, retrofit policies, and heating system transitions.
- Estimation of both annual energy demand and hourly load profiles. 

## Research Applications

SwissRes+ can be used to:

- Assess long-term building decarbonization pathways.
- Quantify the impact of retrofit policies.
- Evaluate heat pump deployment scenarios.
- Estimate future cooling demand under climate change.
- Analyse impacts on electricity demand and peak loads.
- Support national and cantonal energy planning. 

## Model Workflow

```text
Building stock and energy certificate data
        │
        ▼
 Archetype generation
        │
        ▼
 Climate and scenario data
        │
        ▼
 Dynamic simulation (ISO 52016)
        │
        ▼
 Heating & cooling demand
        │
        ▼
 National aggregation
        │
        ▼
 Scenario analysis and results
```

## Repository Contents

This repository currently provides:

- Model outputs used in published research.
- Selected input datasets.

The simulation code is not currently included in the repository. Future releases may provide additional model components depending on data licensing and project developments.

## Associated Publication

Sasso, F., & Patel, M. K. (2026).

> **The effect of climate change, building retrofit and cooling system adoption on future thermal energy demand in the Swiss building stock**
>
> *Energy & Buildings*, 363, 117544.
>
> https://doi.org/10.1016/j.enbuild.2026.117544

## Citation

If you use SwissRes+ results, input datasets, or methodology in your work, please cite the associated publication:

Sasso, F., & Patel, M. K. (2026). *The effect of climate change, building retrofit and cooling system adoption on future thermal energy demand in the Swiss building stock*. Energy & Buildings, 363, 117544. https://doi.org/10.1016/j.enbuild.2026.117544


## Contact

**Francesco Sasso**  
Energy Efficiency Group  
Institute for Environmental Sciences  
University of Geneva

📧 francesco.sasso@unige.ch

## License


The publication associated with SwissRes+ is distributed under the
Creative Commons Attribution 4.0 International (CC BY 4.0) License.

Unless otherwise stated, the documentation and datasets contained in this
repository are made available under the same license.

