# Amphibian-segmentation-clock-models

## Overview

This repository includes simluation code and raw data used to generate figures and results in the manuscript *Amphibian segmentation clock models suggest mechanisms of slowed development across increasing genome size and nuclear volume*. 

## Layout

This repository is split into two main folders: Data and Scripts and local functions __(Data and Scripts_localFunctions folders)__. To look at and/or download raw data (nuclear export simulation results and statistics; period & amplitude of expression across all protein half-life and total delay combinations; critical total delay values across increasing protein half-life), navigate to the Data folder. To download MATLAB scripts used for the simulations described in the paper (which produce the raw data given in the Data folder), navigate to the Scripts_localFunctions folder. 

                  
## Workflow: How to generate a species-specific clock model and associated results

Start with running nuclear export simulations. You may adjust the radius range r and intervals to match your own range of interest; you may also adjust the parameter used to scale x, y, and z particle positions (i.e. the overall trajectory) to match empirical observations (i.e. nuclear export time of a somitogenesis transcript within a nucleus of known/well-estimated radius size).  To incorporate simulation results into a species-specific segmentation clock model, measure or estimate the nuclear radius of PSM cells in your species of interest, and find the corresponding mean export time for that radius. Within nuclear_export_simulations.m, there is a module that generates statistics based on simulation results. The resulting mean export time, for your radius, can be plugged into other scripts as the Texp value - the average transcript export time for that species. The nuclear_export_simulations.m script gives results for Brownian Motion (normal diffusion) and fractional Brownian motion (obstructed diffusion) models with initial particle positions at the origin and drawn from a uniform distribution exlcuding the nuclear periphery. You may run individual modules corresponding to these models based on empirical observations or assumptions about your species of interest, or you may run all modules/models and compare results.

After getting Texp, the following scripts can be used to solve DDE systems and assess the period and/or amplitude of expression: DDE_solns.m, amplitude.m, inc_genestability, and period.m. DDE systems are generated by plugging in species-specific model parameters. Texp should come from nuclear export simulations (keeping in mind which model you would like to use/assume, and what your estimated PSM nuclear radius is). Other model parameters in these scripts include: pcrit, a, k, HL_m, HL_p, Ttx, Tin, and Tp or Ttl (protein synthesis delay = translation delay)  - these should be based on species-speicifc estimates or empirical results and should be adjusted accordingly before running scripts.

Species-specific parameters are also needed to run the senstivity analysis (which takes all model parameters) and analtyical results (which takes mRNA half-life, a range of protein half-life and a crititcal protein threshold value) scripts. 

In depth explanations of what each script does can be found in the README file in the scripts and local functions folder. 





