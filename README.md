# Code for "Estimating Higher-Order Mixed Memberships via the $\ell_{2,\infty}$ Tensor Perturbation Bound" by Joshua Agterberg and Anru Zhang

This repository contains all the code required to analyze the data and reproduce the figures in “Estimating Higher-Order Mixed Memberships via the $\ell_{2,\infty}$ Tensor Perturbation Bound.”

## Key code:
- *HOOI.R*: Implements the higher-order orthogonal iteration algorithm with diagonal-deletion initialization.
- *membership_estimation.R*: Runs the successive projection algorithm to obtain approximate vertices.
- *estimate_rank.R*: Contains the get_Elbows function used to estimate the rank of each mode of the tensor.

## Simulations:
- *sim_two_infty.R*: Runs all the simulations by calling tensor_two_infty.R for individual simulations. The results are then analyzed in analyze_results.R. The outputs from sim_two_infty.R are used to produce Figure 2 in the paper.

## Data Analysis:
- *trade_data.R*: Analyzes the global trade dataset used in the main paper, producing Figures 1, 6, and 7. It uses the following files in the “data” folder:
	- fao_trade_layers.txt
	- fao_trade_multiplex.edges
	- fao_trade_nodes.txt
- *global_flights.R*: Analyzes the global flight network, as discussed in the supplementary material, using the flight_route.RData file.
- *USA_flights.R*: Analyzes the USA flights data presented in the supplementary material. This script uses the US_airport_networks-only48states.RData and uscities.csv, which provides latitude and longitude information for plotting. Additional plotting functions are collected in pltfunctions.R.

## Additional Utilities:
- *misc.R*: Contains additional helper functions for calculating  $\ell_{2,\infty}$ distances between orthonormal matrices.