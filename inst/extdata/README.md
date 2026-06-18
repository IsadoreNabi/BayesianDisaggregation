
# Example Data Files for BayesianDisaggregation

This directory contains example data files for testing and demonstration:

## CPI.xlsx
- Annual Consumer Price Index data (2019-2023)
- Contains aggregate index (Total) and component indices
- Categories: Food, Housing, Transport, Healthcare, Education, Recreation, Other
- Base year: 2019 (index = 100)

## WEIGHTS.xlsx
- Industry weights matrix for CPI components
- Rows: Industries/Categories
- Columns: Years (2019-2023)
- Each year's weights sum to 1.0
- Format: Industry | 2019 | 2020 | 2021 | 2022 | 2023

These are minimal example files for package testing and documentation.
For real analysis, users should provide their own data files with appropriate structure.

