# GrayWhaleCalfProduction

This repository contains R code to extract and analyze gray whale calf production data that are collected at Piedras Blancas, CA, annually. 

Necessary functions are stored in GrayWhaleCalfProduction_fcns.R. Some of the functions in the file may be obsolete and unused. The analytical part was conducted using JAGS (Just Another Gibbs Sampler, which can be downloaded from https://sourceforge.net/projects/mcmc-jags/files/). The model was developed by Josh Stewart. In 2026, an updated analytical approach was introduced due to the previous model's failure when counts declined and were spread over several weeks. The new approach uses the Negative Binomial model instead of the Poisson-Binomial and is coded in Stan (https://mc-stan.org/). The updated model was developed by Tomo Eguchi.

It is assumed that this repository is downloaded to your local computer and a project is created in RStudio. Please create three empty folders within the project folder and name them "data," "figures," and "RData". They are case sensitive. 

In 2023, new data extraction code was developed to address inconsistencies in the previous version (*Collating and Formatting Passdown.R*). The comparison of the two extraction methods is provided in an RStudio notebook file named *compare_data_extraction.Rmd*. The new extraction method was first used for the 2023 season. For the 2023 report, even though the new extraction method was used, the table was replaced with the results from the old extraction version. These reports can be found in *Report_calf_production_2022.Rmd* and *Report_calf_production_2023.Rmd*. The new extraction method was reported in Lang et al. (2024). 

The order of the data processing and analysis should be: 

1. Raw data should be saved in the 'data/All data' directory. Folder names are 'C_C YYYY', where YYYY indicates a four-digit year. There is a space between C_C and YYYY.
2. Edit list.sheet.names.inshore, list.col.types.inshore, and list.col.names.inshore in *GrayWhaleCalfProduction_fcns.R* (if you pull the most recent version, this should be done) 
3. Run *Extract Excel Data_v2.Rmd*. This will create two subfolders within the data folder: 'Processed data v3' and 'Formatted Annual Data v3'. 'v3' refers to the extraction method version, and the user should ignore it. 
4. To run the Jags model (*GWCalfCount_v1.jags*), use *calf_production_jags.R*.
5. Run *Combine inshore and offshore sightings.R* This file combines sightings and creates one file per year. Output files from this script is used to create the input object for the next script.
6. Run *calf_production_hierarchical_stan.R*. This script can run all models below. See Eguchi et al. (2026) for details. This step requires the following Stan code in the 'models' subdirectory:
   1. *GWCalfCount_nb_fixedyear.stan*
   2. *GWCalfCount_nb_shiftingpeak.stan*
   3. *GWCalfCount_nb_shiftingpeak_centered.stan*
   4. *GWCalfCount_nb_singleyear.stan*
7. Edit *gray whale calf production tech memo 2026 v4.qmd* to create a new report. This is a Quarto document. 
8. In the Quarto document, *reference.bib* is used to create citations and the literature cited section. Because it includes my reference library as of the time of creating the 2026 report, this file is too large for GitHub (approx. 7MB) and has not been uploaded. Your reference manager (e.g., Zotero, Mendeley, etc.) can create a .bib reference file. To facilitate formatting the literature cited section, two .csl files for Marine Mammal Science and Marine Ecology Progress Series are included. 
