# Contouring Groundwater Information: The Search For Truth
This repository contains code and example data accompanying the paper "Contouring Groundwater Information: The Search For Truth."

# Abstract
Contour maps of groundwater information, such as water levels, contamination, and physical properties, are essential tools for decision-making in a wide range of hydrological applications. Methods for generating these maps vary from expert-drawn contours combined with simple geometric calculations to advanced spatio-temporal statistical methods and complex multi-dimensional numerical groundwater models. In this study, we evaluate several methods from a decision-making perspective, focusing on reproducibility, quality assurance, and the implications of decisions informed by these methods. Specifically, we compare seven different approaches for mapping groundwater contamination: (a) the three-point method (TPM), (b) thin-plate spline smoothing (TPS), (c) thin-plate splines with a time component (TPS-T), (d) regularized canonical correlation regression with spatial smoothing (rCCR-S), (e) a Bayesian method, (f) space-time kriging using integrated nested Laplace approximation (INLA), and (g) random forest spatial interpolation (RFSI). While this study is not exhaustive, it aims to illustrate the strengths and limitations of these methods through a simulation study and an application to chloroform concentrations in groundwater in Henderson, Nevada, over more than 40 years. Data for the simulation are collected from multiple well locations over time, providing insights into the performance of these methods under realistic conditions. The evaluation considers varying conditions, such as sample size (number of wells) and time step intervals. Results from the simulation study, which is based on a two-dimensional simplification of a groundwater plume, are presented and compared using root mean square prediction error.

## Data
Three datasets are available:
- Five time steps of a groundwater plume simulated in ModFlow for example implementation with the chosen methods (complex_ts_(01-40).rds)
- Data frame with concentration predictions for each method for plotting contours (predictions_complex_plume.RDS)
- MSPE results for one simulated well network across all time points for 15, 30, 50, 75, and 100 wells (sim_complex_df_25_itergroup_4.RDS)

## Code
The script simulation_code.R runs the methods defined in model_functions.R. 

The script plot_sim_results.R shows how results are plotted as seen in the manuscript using the example data provided. 
