# Nutrient Composition & Cognitive Function Analysis (NHANES)

## Overview
This project investigates how nutrients co-occur in broader dietary patterns and whether those dietary patterns are associated with cognitive performance in adults.

Using NHANES 2013–2014 dietary and cognitive function data, factor analysis was applied to reduce high-dimensional nutrient intake variables into interpretable latent dietary factors. Regression modeling was then used to evaluate associations between extracted dietary patterns and composite cognitive performance scores.

## Data
- NHANES 2013–2014 dietary intake data
- Cognitive function assessment data
- 26 selected nutrient variables including macronutrients, vitamins, minerals, cholesterol, and caffeine

## Methods
- Principal Factor Analysis (PFA)
- Varimax rotation
- Scree plot and eigenvalue-based factor selection
- Composite cognitive score standardization
- Multivariate regression modeling

The final model identified 6 interpretable dietary patterns representing broader nutritional behaviors.

## Key Findings
Several extracted dietary patterns showed modest but statistically significant associations with cognitive performance outcomes. Patterns associated with fiber-rich foods, fruits & vegetables, and protein-dense dietary structures demonstrated positive associations with cognitive function scores.

## Tools Used
- SAS
- SAS Studio
- PROC FACTOR
- PROC REG

## Current Enhancements
Current improvements in progress include:
- Incorporating 2-day averaged dietary recall data
- Adding total caloric intake covariates
- Adding demographic control variables
- Expanding regression diagnostics and model robustness checks

## Repository Contents
- SAS analysis code
- Factor analysis output
- Regression diagnostics
- Project write-up
