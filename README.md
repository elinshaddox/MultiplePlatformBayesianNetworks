# MultiplePlatformBayesianNetworks
Associated files for Shaddox, E., Peterson, C.B., Stingo, F.C., Hanania, N., Cruickshank-Quinn, C., Kechris, K., Bowler, R. and Vannucci, M. (2018) Bayesian Inference of Networks Across Multiple Sample Groups and Data Types. Biostatistics, invited revision.


Author: Elin Shaddox
Contact: elin@rice.edu

The provided Matlab, R, and supplementary files for Bayesian inference of multiple graphical models are associated with the following publication:

Shaddox, E., Peterson, C.B., Stingo, F.C., Hanania, N., Cruickshank-Quinn, C., Kechris, K., Bowler, R. and Vannucci, M. (2018) Bayesian Inference of Networks Across Multiple Sample Groups and Data Types. Biostatistics, invited revision.

These scripts rely on functions from the Matlab code for G-wishart sampling provided by Hao Wang at https://msu.edu/~haowang/ and are associated with the following publications
Associated publications:
H. Wang, Scaling It Up: Stochastic Search Structure Learning in Graphical Models Bayesian Analysis 10 (2015): 351-377

Wang, H. and Li, S. (2012). Efficient Gaussian graphical model determination
under G-Wishart prior distributions. Electronic Journal of Statistics.
6: 168—198.

Shaddox, E., Stingo, F., Peterson, C.B., Jacobson, S., Cruickshank-Quinn, C., Kechris, K., Bowler, R. and Vannucci, M. (2016). A Bayesian Approach for Learning Gene Networks Underlying Disease Severity in COPD. Statistics in Biosciences, in press.
Please cite all publications if you use this code. Thanks!

OVERVIEW OF FILES 
—————————————————

Example_multiple_graphs_SSVS_with_platforms.m
Basic example of running the MCMC sampler and generating results summaries on a simple setting with:
Platform 1 - 3 groups with identical dependence structure
Platform 2 - 3 groups with differing dependence structure
—————————————————

MCMC_multiple_graphs_SSVS_with_platforms.m
Code for running the MCMC sampler
—————————————————

calc_mrf_C.m
Helper function for calculating the normalizing constant for the MRF prior
—————————————————

Data Generation for Simulations
Scripts to generate matrices similar to those from simulation Settings:
- ScaleFreeSimTruths_p40.R provides input (true precision matrices) for generating data for p=40 simulations in ScaleFreeSimDataGeneration.m
- count_shared_edges.m is a matlab function for counting edges between two adjacency matrix inputs
- AR2_SimulationTruthsAndDataGeneration.m provides set up for true precision matrices and data generation in AR(2) p=80 simulation settings
—————————————————

Metabolite Selection Example
Script, plots, and generated data to provide an example of correlated metabolite data and the iterative metabolite selection process from the publication listed above:
- MetaboliteExampleSelection.R script for data generation and iterative process
- corData_MetabSelectionExample.csv generated data based on correlation of a metabolite subset
- BEFORESELECTIONgeneratedCorrDataPlots.png plot of correlated data before selection
- AFTERSELECTIONgeneratedCorrDataPlots.png plot of less correlated data subset after selection
