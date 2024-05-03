# composite-dyadic-models
Julia and R Code for Schwob, Hooten, Narasimhan, "Composite Dyadic Models for Spatio-Temporal Data," in review.

## aDNA
Contains the code used to analyze the Bronze Age Europe human genomes.
- changePops.R : changes the "population identifier" for each genome to an individual identifier to compute genetic dissimilarites between all individuals
- filter.R : takes the subset of data that we used for analysis in our paper
- CRPS.jl : computes the continuous rank probability scores (CRPSs) for the competing models
- mcmc_both/mcmc_ds/mcmc_dt/mcmc_none.jl : contains the MCMC algorithm for our Bayesian hierarchical model fully specified in Web Appendix B with the corresponding composite weighting scheme
- meta.jl : the master script that runs all other scripts
- plotting.jl : plots the potential surface and MCMC diagnostics
- script.jl : structures the data for use in the MCMC algorithm

The aDNA data is made publicly available at the [Allen Ancient DNA Resource](https://reich.hms.harvard.edu/allen-ancient-dna-resource-aadr-downloadable-genotypes-present-day-and-ancient-dna-data).

## ar1
Contains the code used to simulate the mechanistic population dynamic model in Web Appendix A.
- ar1.r : simulates data from the model, analyzes the simulated data via MCMC, and (when using composite weights) adequately recovers the true parametric values

## rupicapra
Contains the code used to analyze the *R. rupicapra* landscape genomic data in Web Appendix C.
- CRPS.jl : computes the continuous rank probability scores (CRPSs) for the competing models
- mcmc_ds/mcmc_none.jl : contains the MCMC algorithm for our Bayesian hierarchical model fully specified in Web Appendix B with the corresponding composite weighting scheme
- plotting.jl : plots the potential surface and MCMC diagnostics
- script.jl : structures the data for use in the MCMC algorithm
