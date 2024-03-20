# Data-driven simulations to assess the impact of study imperfections in time-to-event analyses

R scripts and data to reproduce the two examples of data-driven simulations presented in the manuscript:

Abrahamowicz M, Beauchamp ME, Boulesteix AL, Morris TP, Sauerbrei W, Kaufman JS, on behalf of the STRATOS Simulation Panel (SP). Data-Driven Simulations to Assess the Impact of Study Imperfections in Time-to-Event Analyses. *Am J Epidemiol* 2024.

## Content

Code written by Marie-Eve Beauchamp. For questions or comments about the code contact marie-eve.beauchamp (at) rimuhc.ca.

The code was written in R with the following software versions:
- R version 4.2.0 (2022-04-22 ucrt)
- Platform: x86_64-w64-mingw32/x64 (64-bit)
- Running under: Windows 10 x64 (build 17134)
- Using R packages:
  - survival version 3.3-1 
  - PermAlgo version 1.2  

Relevant files for Example 1 and Example 2 of the manuscript are included in the folder with the matching name:  

## /Example1

Example 1 of data-driven simulations is run by the script `Code_ex1_simulations.R`, which:
- Calls the functions from the script `Code_ex1_functions.R`
- Saves detailed simulation results in a separate file for each scenario (`SimResults_sc1.RData`, ..., `SimResults_sc8.RData`)
- Reproduces Table 1 of the manuscript (saved in `Table1.RData`).

## /Example2

Example 2 of data-driven simulations is run by the script `Code_ex2_simulations.R`, which:
- Calls the functions from the scripts `Code_ex2_functions.R` and `Code_permalgorithm.realDat.R`
- Loads the original synthetic data from `dataOriginal_ex2.RData` 
- Saves detailed simulation results in a separate file for each scenario (`SimResults_sc1.RData`, ..., `SimResults_sc4.RData` and `SimResults_QBA.RData`) 
- Reproduces Table 2 of the manuscript and Web Table 3 (saved in `Table2.RData` and `WebTable3.RData`).

The script `Code_permalgorithm.realdat.R` contains the documentation and code for the function 'permalgorithm.realdat', which implements the adaptation of the permutational algorithm to input data including time-dependent covariates with values known up to different follow-up times across subjects. This version of the function 'permalgorithm.realdat' was used for the simulations for Example 2. For the latest version of this function, see the function 'permalgorithm.realdat' that will be implemented in the next version of the PermAlgo package on CRAN (after version 1.2). The function 'permalgorithm.realdat' is complementary to the function 'permalgorithm' already included in the PermAlgo R package. See the documentation for 'permalgorithm.realdat' for further details.
