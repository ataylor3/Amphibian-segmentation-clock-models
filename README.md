# Amphibian-segmentation-clock-models

## Overview

This repository includes simluation code and raw data used to generate figures and results in the manuscript *Amphibian segmentation clock models suggest mechanisms of slowed development across increasing genome size and nuclear volume*. 

## Layout

This repository is split into two main folders: Data and Scripts and local functions. To look at and/or download raw data (nuclear export simulation results and statistics; period & amplitude of expression across all protein half-life and total delay combinations; critical total delay values across increasing protein half-life), navigate to the Data folder. To download MATLAB scripts used for the simulations described in the paper (which produce the raw data given in the Data folder), navigate to the Scripts_localFunctions folder. 

### Data

#### Nuclear export simulation results and statistics

#### Period of oscillation

#### Amplitude of oscillation

#### Impact of increased gene stability on period of oscilltion

#### Analytical results: Tcrit across increasing protein stability 

### Scripts and local functions

#### Scripts 

DDE_solns.m 
Gives the code and general XX used to solve the DDE systems that different parameter combinations yield; we include the local functions 
called to solved DDE systems and assess oscillatory behavior (period and amplitude of gene expression).
We also use this general method of solving to show that for the same DDE system, the period of mRNA and protein expression are very nearly
equivalent (Supplement 1); we show this for a DDE system corresponding to X. laevis and A. mexicanum Brownian Motion models. 

amplitude.m
Script used to generate amplitude plots. Resulting amplitude data matrices are given for every model shown in Figures S6 and S7 in Amplitude of expression.xlsx 

analytical_results.m
Script that uses steady state solution to Lewis' DDE system to generate model specific plots of critical delay vale Tcrit across increasing 
protein half-life. Can be used to generate plots in Figures S1 and S2 (when model-specific parameters are plugged in).

inc_genestability.m
Script used to test if additional increases in mRNA and protein stability can act to recapitulate the 155 minute period of gene expression/
rate of somite segmentation seen in A. mexicanum. Can be used to generate plots seen in Figure S6. Corresponding period data matrices are given in Amex_BM_hmTexp_incGeneProductStability.xlsx

nuclear_export_simulations.m
Script used for nuclear export simulations. Used to generate species- and diffusion-specific estimates for export delay (Texp), and used to 
generate plots in Figures 3, S3 and S4. Data matrices for each model type (Brownian Motion with initial position at the origin; fractional Brownian 
Motion with initial position at the origin; Brownian Motion with initial position drawn from a uniform distribution; and fractional Brownian 
Motion with initial position drawn from a uniform distribution) are given in nuclear export simulations.xlsx. Matrices are 10,000 by 26 where rows correspondto the 10,000 iterations and columns correspond to increasing nuclear radii (0.5 to 13 µm with increments of 0.5 µm); this data is used to generate export time disributions, and the mean export time for each nuclear radius. Statistics (mean, standard deviation, and coefficient of variation for each radius length for each model type) are also given in nuclear export simulations.xlsx

period.m 
Script used to generate period plots in Figures 1 and 2. Resulting period data matrices for each model are given in Period of expression.xlsx

sensitivity analysis.m
Script used to carry out sensitivty analyses shown in Tables 3 and S2. 

#### Local functions

Amp.m - gives average amplitude of a DDE solution

P.m - gives average periof of a DDE solution 

ddefun_nested.m - uses ddesd to solve DDE systems corresponding to the different models and parameter combinations

history.m - sets intial mRNA and protein molecule values to 100 

osc_behavior.m - stores local extrema and time stamps in vectors, cuts first 5 cycles of oscillation, and uses local functions Amp.m and P.m
                  to retrieve period and amplitude of gene expression for each solution
                  
## How to use



