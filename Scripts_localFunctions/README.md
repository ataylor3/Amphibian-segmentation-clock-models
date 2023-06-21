# Scripts and local functions

All script and local function files include comments with step-by-step explanations of the code. 

## Scripts 

__DDE_solns.m__ 
Gives the code and general flow used to solve the DDE systems that different parameter combinations yield; we include the local functions 
called to solved DDE systems and assess oscillatory behavior (period and amplitude of gene expression).
We also use this general method of solving to show that for the same DDE system, the period of mRNA and protein expression are very nearly
equivalent (Supplement 1); we show this for a DDE system corresponding to X. laevis and A. mexicanum Brownian Motion models. 

__amplitude.m__
Script used to generate amplitude plots. Resulting amplitude data matrices are given for every model shown in Figures S6 and S7 in Amplitude of expression.xlsx 

__analytical_results.m__
Script that uses steady state solution to Lewis' DDE system to generate model specific plots of critical delay vale Tcrit across increasing 
protein half-life. Can be used to generate plots in Figures S1 and S2 (when model-specific parameters are plugged in).

__inc_genestability.m__
Script used to test if additional increases in mRNA and protein stability can act to recapitulate the 155 minute period of gene expression/
rate of somite segmentation seen in A. mexicanum. Can be used to generate plots seen in Figure S6. Corresponding period data matrices are given in Amex_BM_hmTexp_incGeneProductStability.xlsx

__nuclear_export_simulations.m__
Script used for nuclear export simulations. Used to generate species- and diffusion-specific estimates for export delay (Texp), and used to 
generate plots in Figures 3, S3 and S4. Data matrices for each model type (Brownian Motion with initial position at the origin; fractional Brownian 
Motion with initial position at the origin; Brownian Motion with initial position drawn from a uniform distribution; and fractional Brownian 
Motion with initial position drawn from a uniform distribution) are given in nuclear export simulations.xlsx. Matrices are 10,000 by 26 where rows correspondto the 10,000 iterations and columns correspond to increasing nuclear radii (0.5 to 13 µm with increments of 0.5 µm); this data is used to generate export time disributions, and the mean export time for each nuclear radius. Statistics (mean, standard deviation, and coefficient of variation for each radius length for each model type) are also given in nuclear export simulations.xlsx

__period.m__
Script used to generate period plots in Figures 1 and 2. Resulting period data matrices for each model are given in Period of expression.xlsx

__sensitivity analysis.m__
Script used to carry out sensitivty analyses shown in Tables 3 and S2. 

## Local functions

__Amp.m__ 
Gives average amplitude of a DDE solution

__P.m__ 
Gives average periof of a DDE solution 

__ddefun_nested.m__ 
Uses ddesd to solve DDE systems corresponding to the different models and parameter combinations

__history.m__ 
Sets intial mRNA and protein molecule values to 100 

__osc_behavior.m__ 
Stores local extrema and time stamps in vectors, cuts first 5 cycles of oscillation, and uses local functions Amp.m and P.m
to retrieve period and amplitude of gene expression for each solution
