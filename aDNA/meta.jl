################
####
#### Meta Script
####
################

using RCall

nMCMC = 500000
startDate = 6000
endDate = 4500

stringDate = "$startDate-$endDate"
mkdir("$stringDate")
mkdir("$stringDate/none")
mkdir("$stringDate/dt")
mkdir("$stringDate/ds")
mkdir("$stringDate/both")

@rput startDate endDate stringDate
R"""
print(Sys.time())

## Population Identifiers; turns all population identifiers the same as the individual ID
source("changePops.R")

## Filter data ; creates anno_R.csv
source("filter.R")
"""

## Script to Prepare Data & Run MCMC Algorithm
include("script.jl")

## Script to Print MCMC Output ; trace plots, potential surface, and eta surface
include("plotting.jl")

## Script for model selection using CRPS
include("CRPS.jl")

print(Sys.time())