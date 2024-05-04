#######################################
####
#### Script for Alpine chamois LG study
####
#######################################

println("\n\n%%%%%%%%%%%%%%%%%%%%%%%%%%% Data Preparation %%%%%%%%%%%%%%%%%%%%%%%%%%%")

###
### Set-up
###

using LinearAlgebra, Random, Statistics, Distributions, SparseArrays, StatsBase # stats
using RCall # R
using CSV, DataFrames, Query # data
using Dates, ProgressBars # progress
using JLD2 # environments
using Distances # for euclidean distance
using MCMCChains, MCMCDiagnosticTools, StatsPlots, Plots # for chains
using GMT # to obtain elevation data
using BlockDiagonals # to quickly construct block-diagonal matrices
Random.seed!(702)

println("Time started: ", Dates.Time(Dates.now()), ".")
nMCMC = 100000

###
### Prepare data in R
###

R"""
## load data
library(adegenet)
data(rupica)

## count observations
n = 335
N = n*(n-1)/2

## genetic dissimilarities
y = dist.genpop(as.genpop(rupica$tab), diag = TRUE, upper = TRUE)

## geographic covariates
library(raster)
elev = raster(rupica$other$mnt)
rupicaLocs = rupica$other$xy
locIdx = cellFromXY(elev, rupicaLocs)

## get elevation for each cell 
nodeElev = extract(elev, rupicaLocs)

## grid to predict on 
gridPred = rasterToPoints(elev)
"""
@rget y locIdx rupicaLocs n N nodeElev gridPred

println("Number of observations: ", n)

## unique locations
obs_locs = rupicaLocs
unique_locs = unique(obs_locs, dims = 1)
s = size(unique_locs)[1]
N = Int(N)
n = Int(n)

## get indices for unique locations
uniqueIdx = zeros(n)
for i in 1:n 
    for j in 1:s 
        if (rupicaLocs[i, 1] == unique_locs[j, 1]) & (rupicaLocs[i, 2] == unique_locs[j, 2])
            uniqueIdx[i] = j 
        end
    end
end
uniqueIdx = Int.(uniqueIdx)

###
### Covariate Matrix 
###

Xnode = [nodeElev nodeElev.^2]

@rput Xnode
R"""
XnodeSC = Xnode
for(i in 1:dim(XnodeSC)[2]){
    XnodeSC[, i] = scale(Xnode[, i])
}
e1Scale = unlist(attributes(scale(Xnode[, 1]))[3])
e1Center = unlist(attributes(scale(Xnode[, 1]))[2])
e2Scale = unlist(attributes(scale(Xnode[, 2]))[3])
e2Center = unlist(attributes(scale(Xnode[, 2]))[2])
"""
@rget XnodeSC e1Scale e1Center e2Scale e2Center

###
### Scale Prediction Grid to Same Scale
###

gridCovs = zeros(size(gridPred)[1], 2)
gridCovs[:, 1] = (gridPred[:, 3] .- e1Center)./e1Scale
gridCovs[:, 2] = (gridPred[:, 3].^2 .- e2Center)./e2Scale

###
### Get Regression Variables 
###

p = 3 # intercept, elev, elev^2
H = zeros(N, n) # mapping matrix for node-level
K = zeros(N, s) # mapping matrix for unique locations
locMat = zeros(N, s) # for variogram analysis
X = ones(N, p) # design matrix
dist = zeros(N) # vector of euclidean distances

idxi = 1 # ith index
idxj = 2 # jth index
for i in 1:N

    ## for phi list
    dist[i] = euclidean(rupicaLocs[idxi, :], rupicaLocs[idxj, :])

    ## X matrix
    X[i, 2] = XnodeSC[idxi, 1] - XnodeSC[idxj, 1] # difference in elevation
    X[i, 3] = XnodeSC[idxi, 2] - XnodeSC[idxj, 2] # difference in elevation^2

    ## mapping matrix H
    H[i, idxi] = 1
    H[i, idxj] = 1 

    ## mapping matrix for variogram analysis
    if uniqueIdx[idxi] == uniqueIdx[idxj] # if the nodes are co-located
        locMat[i, uniqueIdx[idxi]] = 2 # we take two of the same etas
    else # if the nodes are in different locations
        locMat[i, uniqueIdx[idxi]] = 1 # we take one eta...
        locMat[i, uniqueIdx[idxj]] = 1 # ... for each location
    end

    # new mapping matrix K
    if dist[i] == 0 # if the nodes are co-located
        K[i, uniqueIdx[idxi]] = 0 # the etas cancel out and we get nothing
    else # if the nodes are in different locations
        K[i, uniqueIdx[idxi]] = 1
        K[i, uniqueIdx[idxj]] = -1
    end
    
    ## update individual indices
    if idxj == n
        global idxi += 1
        global idxj = idxi + 1
    else
        global idxj += 1
    end
end

###
### Pre-compute ϕ and R
###

## pre-compute phi
phi = 0:5:(maximum(dist)/4)

println("ϕ range: ", phi)

## pre-compute R
R_list = ones(s, s, length(phi))
println("Precomputing R(ϕ):")
for k in ProgressBar(1:length(phi))
    for i in 1:(s-1) # for each unique location
        for j in (i+1):s # for each unique location
            R_list[i, j, k] = exp(-euclidean(unique_locs[i, 1:2], unique_locs[j, 1:2])/phi[k])
            R_list[j, i, k] = R_list[i, j, k]
        end
    end
    R_list[:, :, k] = R_list[:, :, k] + 0.0001*I
end

###
### For covariance weights 
###

@rput dist 
R"""
library(scales)
dtij = rep(0, length(dist))
dsij = rescale(dist, to = c(0.00001, 1))
"""
@rget dtij dsij

###
### Get a variogram
###

@rput y X dist locMat unique_locs s
R"""
euclidean <- function(a, b) sqrt(sum((a - b)^2))
library(dplyr)
png("vario.png")

nSeq = 25
dSeq = seq(0, max(dist)/4, length.out = nSeq)
midPoints = dSeq[-length(dSeq)] + diff(dSeq)/2

fit = lm(y ~ X[, -1] + locMat)
eta = fit$coefficients[-(1:3)]
eta[length(eta)] = fit$coefficients[1]

pairNo = rep(0, nSeq-1)
vario1 = rep(0, nSeq-1)
for(i in 1:(nSeq-1)){
    tmpSum = 0
    counter = 0
    for(j in 1:(s-1)){
        for(k in j:s){
            if(between(euclidean(unique_locs[j, 1:2], unique_locs[k, 1:2]), dSeq[i], dSeq[i+1])){
                counter = counter + 1
                tmpSum = tmpSum + (eta[j] - eta[k])^2
            }
        }
    }
    vario1[i] = tmpSum/counter 
    pairNo[i] = counter
}
woT = plot(midPoints, vario1, xlab = "distance", ylab = "Variogram", type = "l", main = "Variogram")
dev.off()

png("resids.png")
plot(dist, fit$residuals^2, main = "Residuals v. Distance")
dev.off()

## print number of pairs in each bin 
cat("Number of pairs in bins: ", pairNo)
"""

###
### Save Julia Output
###

@save "juliaVars.jld2"

###
### Run MCMC Algorithm
###

println("\n%%%%%%%%%%%%%%%%%%%%%%%%%%% Running MCMC Algorithms %%%%%%%%%%%%%%%%%%%%%%%%%%%")
startTime = Dates.Time(Dates.now())

include("mcmc_none.jl")
include("mcmc_ds.jl")

println("MCMC Algorithm (Neither):")
out = mcmc_none(y, X, H, K, phi, R_list, nMCMC, dtij, dsij)

println("MCMC Algorithm (ds):")
out3 = mcmc_ds(y, X, H, K, phi, R_list, nMCMC, dtij, dsij)

println("Time took: ", Dates.Time(Dates.now()) - startTime, ".")

###
### Plotting 
###

include("plotting.jl")

###
### Model Comparisons w/ CRPS
###

include("CRPS.jl")