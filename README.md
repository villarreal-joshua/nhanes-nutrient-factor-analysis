# Nutrient Composition Covariance & Cognitive Function (NHANES)

## Overview
This project investigates how nutrients co-occur in dietary patterns and whether those latent dietary patterns are associated with cognitive performance.

Using NHANES 2013–2014 data, we applied factor analysis to nutrient intake variables to extract interpretable dietary factors and then evaluated associations with a composite cognitive score.

## Data
- NHANES 2013–2014 dietary intake data
- Cognitive function measures (standardized tests combined into a composite score)
- Selected 26 nutrient variables (filtered from a larger set for interpretability and to reduce redundancy)

## Methodology
- Principal factor extraction
- Varimax rotation for interpretability
- Factor selection using scree plot, eigenvalues (> 1), and residual diagnostics (RMS off-diagonal residuals)
- Regression modeling of cognitive composite score on extracted dietary factors

## Key Findings
Several extracted dietary patterns showed statistically significant (but modest) associations with cognitive performance. Patterns related to processed/animal protein, fiber-rich foods, and fruits/vegetables were positively associated, while sugars and caffeine patterns were not significant.

## Tools Used
- SAS (PROC FACTOR, regression procedures)
- Factor analysis and diagnostics
- Regression modeling
