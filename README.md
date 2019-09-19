# biomarker_surrogacy
Evaluating the surrogacy of multiple vaccine-induced immune response biomarkers in HIV vaccine trials

The Folder contains the following files:

1. fun_general.R - Contains all the major functions to implement our method, including the main function BI_SS which runs the bootstrap imputed stability selection procedure.
2. simu_set_gauss_contw_mainlow.R - Simulation setting I from our paper (Table 1).
3. simu_set_gauss_contw_mainmedium.R - Simulation setting II from our paper (Table 2).
4. simu_set_gauss_contw_mainhigh.R - Simulation setting III from our paper (Table S1).
5. simu_set_gauss_disw.R - Simulation setting IV from our paper (Table 3).
6. simu_set_nongauss_disw.R - Simulation setting V from our paper (Table S2).
7. selection_estimation.R - The main code to produce results given in Tables 1-3, S1 and S2. Need to be run after generating data using one of the previous 5 simulation setting codes.
8. sampledat_analysis.R - Sample code to run the BI_SS procedure with one of the several imputation methods (parametric, nonparametric or MICE) on the synthetic data (sampledat.csv) provided in this folder.
9. sampledat.csv - Synthetic data on which our methods can be applied (please consult sampledat_analysis.R to run the BI_SS procedure on this dataset).
