# 
Datasets and codes for the manuscript: 

Linking mode of seed dispersal and climatic niche evolution in flowering plants

Thais Vasconcelos, James Boyko, and Jeremy Beaulieu

# files 

Rosaceae_niche_partial.csv - removed Cotoneaster_morrisonensis and Sorbus_randaiensis as outliers in aridity index

21.01.11-runCorHMM.R - runs the corhmm models 

21.01.12-analyzeCorHMM.R - preliminary analysis of corhmm models 

21.05.04-CompleteRun.R - a script which has all the necessary code to complete the OU part of the analyses. generates simmaps proportional to the model averages of corHMM modeling. then runs the OU models on 1000 generated simmaps for each dataset (prec and temp with and without SE). model averages the tips so that we can values associated with abiotic and biotic parameters. produces preliminary figures.

21.02.03-reconCorHMM.R - model averages the corhmm fits and then does an ancestral state reconstruction. this is the same model averaged corhmm rate matrix used to generate simmaps.

21.02.04-cleaningGBIF.R - filters GBIF points and organize climatic data for OUwie

21.02.04-treePlots.R - code for some descriptive plots

21.02.10-functionsClimate.R - functions used to filter GBIF data and to extract summary statistics from climate layers

summ_results.R - code to summarize table of parameters

21.03.04-plotAICwts.R - code to make plot with aicw for OUwie results

plot_corhmm.R - code to plot aicw for corHMM

SummaryTable_corHMM.xslx - a one stop shop to navigate through the results saved in table_corhmm/ 

# folders

trees/ - time-calibrated phylogenies in newick format

trait_data/ - fruit types and organized climate niche datasets for OUwie

climate_data/ - raw climatic data from filtered points

res_corhmm/ - save files from the corHMM generated by 21.01.11-runCorHMM.R 

res_ouwie/ - save files from the OUwie generated by 21.05.04-CompleteRun.R 

table_corhmm/ - csv files of the rate matrices, weighted ASR, and model summary tables from the corhmm analysis. generated by 21.01.12-analyzeCorHMM.R 

recon_corhmm/ - model averaged ancestral state reconstructions save files generated by 21.02.03-reconCorHMM

tip_avg_tables/ - csv files which contain the tip model averaged parameters. generated at the end of 21.05.04-CompleteRun.R and can be easily run locally.

# Citations for GBIF datasets:
Apocynaceae: GBIF.org (15 February 2021) GBIF Occurrence Download https://doi.org/10.15468/dl.sbde8b 

Ericaceae: GBIF.org (14 December 2020) GBIF Occurrence Download https://doi.org/10.15468/dl.a56dvd 

Melastomataceae: GBIF.org (04 March 2021) GBIF Occurrence Download https://doi.org/10.15468/dl.agmsjd

Rosaceae: GBIF.org (14 December 2020) GBIF Occurrence Download https://doi.org/10.15468/dl.67u798 

Solanaceae: GBIF.org (14 December 2020) GBIF Occurrence Download https://doi.org/10.15468/dl.hphjbs 
 

