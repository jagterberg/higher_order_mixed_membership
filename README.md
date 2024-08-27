# Code for "Estimating Higher-Order Mixed Memberships via the \ell_{2,\infty} Tensor Perturbation Bound" by Joshua Agterberg and Anru Zhang

This repository contains all the code needed to analyze the data and reproduce the figures in the paper referenced above.  

## Key code:

The following code contains the main methodology of our paper
- HOOI.R (runs the higher-order orthogonal iteration algorithm with diagonal-deletion initialization)
- membership_estimation.R (runs the successive projection algorithm to obtain approximate vertices)

In addition, the file estimate_rank.R contains the get_Elbows function used to estimate the rank of each mode of the tensor.

The simulations are run using sim_two_infty.R, which calls tensor_two_infty.R for a single simulation, and then the results are analyzed in analyze_results.R. To run the simulations, all that is needed is sim_two_infty.R, which runs all the simulations using tensor_two_infty.R.  The outputs of sim_two_infty.R are used to produce figure 2. 

The file trade_data.R is used to analyze the global trade dataset in the main paper and is used to produce Figures 1, 6, and 7.  It uses the files fao_trade_layers.txt, fao_trade_multiplex.edges, and fao_trade_nodes.txt in the data folder.

The file global_flights.R is used to analyze the global flight network in the supplementary material.  This uses only the data flight_route.RData.

Finally, the file USA_flights.R is  used to analyze the time series of USA flights in the supplementary material.  This uses the data US_airport_networks-only48states.RData and the file uscities.csv, which contains helpful information on latitude and longitude for plotting.  This file also uses the code pltfunctions.R for additional plotting.

Finally, the file misc.R contains additional helper functions to calculate ell_{2,\infty} distances between orthonormal matrices.  