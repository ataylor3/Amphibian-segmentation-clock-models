# Amphibian-segmentation-clock-models

## Overview

This repository includes simluation code and raw data used to generate figures and results in the manuscript *Amphibian segmentation clock models suggest mechanisms of slowed development across increasing genome size and nuclear volume*. 

## Layout

This repository is split into two main folders: Data and Scripts and local functions. To look at and/or download raw data (nuclear export simulation results and statistics; period & amplitude of expression across all protein half-life and total delay combinations; critical total delay values across increasing protein half-life), navigate to the Data folder. To download MATLAB scripts used for the simulations described in the paper (which produce the raw data given in the Data folder), navigate to the Scripts_localFunctions folder. 

### Data

__nuclear export simulations.xlsx__

Excel book that includes raw data of first exit times (in steps or seconds) for particles simulated in sphere of corresponding nuclear radii ranging from 0.5 to 13 µm with intervals of 0.5 µm. First exit time (step, seconds) for particles is given for 10,000 iterations for each nuclear radius size. Different sheets correspond to different model types: FPT_BM_origin gives the raw data for Brownian Motion simulations for which particles start at the origin; FPT_fBM_origin gives the raw data for fractional Brownian Motion simulations for which particles start at the origin; FPT_BM_uniformDist gives the raw data for Brownian Motion simulations for which initial particle position is drawn from a uniform distribution; FPT_fBM_uniformDist gives the raw data for fractional Brownian Motion simulations for which initial particle position is drawn from a uniform distribution. Each model data sheet is followed by a correpsonding statistics sheet for which raw data has been converted into mean exit/export time (minutes), standard deviation (minutes), and the coefficient of variation for each radius size: BM_origin_statistics, fBM_origin_statistics, BM_uniformDist_statistics, and fBM_uniformDist_statistics. 

__Period of expression.xlsx__

Excel book that gives periods of gene expression for each protein half-life/total delay time combination for every model shown in figures 1 and 2:

Xlaevis_BM corresponds to figure 1A; Xlaevis_fBM corresponds to figure 1B; Amexicanum_BM_hm3 corresponds to figure 1C; and Amexicanum_fBM_hm3 corresponds to figure 1D.

Amexicanum_BM_hmTexp corresponds to figure 2A; Amexicanum_fBM_hmTexp corresponds to figure 2B; Amexicanum_BM_hm1 2Texp corresponds to figure 2C; Amexicanum_fBM_hm1 2Texp corresponds to figure 2D; Amexicanum_BM_hm1 4Texp corresponds to figure 2E; and Amexicanum_fBM_hm1 4Texp corresponds to figure 2F.  

__Amplitude of expression.xlsx__

Excel book that gives amplitudes of gene expression for each protein half-life/total delay time combination for every model shown in figures S6 and S7:

Xlaevis_BM corresponds to figure S6A; Xlaevis_fBM corresponds to figure S6B; Amexicanum_BM_hm3 corresponds to figure S6C; and Amexicanum_fBM_hm3 corresponds to figure S6D.

Amexicanum_BM_hmTexp corresponds to figure S7A; Amexicanum_fBM_hmTexp corresponds to figure S7B; Amexicanum_BM_hm1 2Texp corresponds to figure S7C; Amexicanum_fBM_hm1 2Texp corresponds to figure S7D; Amexicanum_BM_hm1 4Texp corresponds to figure S7E; and Amexicanum_fBM_hm1 4Texp corresponds to figure S7F. 

__AmexBM_hmTexp_incGeneProductStability__

Excel book that gives periods of gene expression for each protein half-life/total delay time combination for every model shown in figure S5:

Amexicanum_BM_hmTexp corresponds to figure S5A; Amexicanum_BM_25%inc_hmTexp corresponds to figure S5B; and Amexicanum_BM_50%inc_hmTexp corresponds to figure S5C. 

__Analytical results: Tcrit across increasing protein stability__

ADD

### Scripts and local functions

#### Scripts 

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

#### Local functions

__Amp.m__ Gives average amplitude of a DDE solution

__P.m__ Gives average periof of a DDE solution 

__ddefun_nested.m__ Uses ddesd to solve DDE systems corresponding to the different models and parameter combinations

__history.m__ Sets intial mRNA and protein molecule values to 100 

__osc_behavior.m__ Stores local extrema and time stamps in vectors, cuts first 5 cycles of oscillation, and uses local functions Amp.m and P.m
to retrieve period and amplitude of gene expression for each solution
                  
## How to use scripts and local functions

### Nuclear export simulations
How to use nuclear_export_simulations.m :
- choose radius range (r)

- scale trajectories to match empirical observations (i.e. nuclear export time of a transcript of interest within a nucleus of known/well-estimated
radius size)


### Solving DDE system and assessing period & amplitude 
How to use: DDE_solns.m ; amplitude.m ; inc_genestability ; period.m


### Sensitivity analysis
How to use: sensitivity_analysis.m 

### Analytical results
How to use: analytical_results.m




