######################################
####
#### Script for Dyadic Regression MCMC
####
######################################

println("\n\n%%%%%%%%%%%%%%%%%%%%%%%%%%% Preparing Data %%%%%%%%%%%%%%%%%%%%%%%%%%%")

###
### Set-up
###

using LinearAlgebra, Random, Statistics, Distributions, SparseArrays, StatsBase # standard stats libraries
using RCall # R
using CSV, DataFrames, Query # data
using Dates, ProgressBars # progress
using JLD2 # environments
using Distances # for distance functions
using MCMCChains, MCMCDiagnosticTools, StatsPlots, Plots # for chains
using GMT # to obtain elevation data
using BlockDiagonals # to quickly construct block-diagonal matrices
Random.seed!(123)

println("Time started: ", Dates.Time(Dates.now()), ".")
include("mcmc_none.jl")
include("mcmc_dt.jl")
include("mcmc_ds.jl")
include("mcmc_both.jl")

###
### Load Data
###

## Get population information & PCA values
R"""
library(geodist)
anno = read.csv("Data/anno_R.csv")
pca = read.delim("Data/pca.txt")
"""
@rget anno pca

println("Number of observations: ", size(anno)[1])

###
### Filter Data
###

## find which individuals are in both anno and pca
in_pca = in(pca[:, 1])
keepIdx = in_pca.(anno.Genetic_ID)
individuals = anno.Genetic_ID[keepIdx]
noInd = length(individuals)

## filter 
@rput individuals noInd keepIdx
R"""
library(tidyverse)
anno = anno %>% filter(keepIdx)
"""
@rget anno

###
### Compute GDM using PCA
### Order is (ind1/ind1, ind1/ind2, ind1/ind3, ..., ind1/ind1920;
###           ind2/ind1, ind2/ind2, ind2/ind3, ..., ind2/ind1920; etc.)
###

weights = vec([10.678 4.259 3.139 2.666 2.397 2.319 2.286 2.242 2.235 2.205]) # from Vagheesh

## idx for each individual in pca data set
pcaIdx = zeros(noInd) # will hold the row index for each individual in the PCA data
println("\nGet indices:")
for i in ProgressBar(1:noInd)
    for j in 1:size(pca)[1]
        if(pca[j, 1] == individuals[i]) # if the row contains the individual
            pcaIdx[i] = j
        end
    end
end

## compute GDM
GDM = zeros(noInd, noInd)
println("\nCompute GDM:")
for i in ProgressBar(1:(noInd-1))
    for j in (i+1):noInd
        GDM[i, j] = wsqeuclidean(pca[Int(pcaIdx[i]), 2:11], pca[Int(pcaIdx[j]), 2:11], weights)
        GDM[j, i] = GDM[i, j]
    end
end

###
### Location of each sample
###

## locations & times of each individual
lon = zeros(noInd)
lat = zeros(noInd)
times = zeros(noInd)
for i in 1:noInd
    lon[i] = anno |> @filter(_.Genetic_ID == individuals[i]) |> @select(:Long_) |> DataFrame |> Matrix |> mean
    lat[i] = anno |> @filter(_.Genetic_ID == individuals[i]) |> @select(:Lat_) |> DataFrame |> Matrix |> mean
    times[i] = anno |> @filter(_.Genetic_ID == individuals[i]) |> @select(:Date_BP) |> DataFrame |> Matrix |> mean
end

## matrix of observed and unique locations 
obs_locs = [lon lat]
unique_locs = unique(obs_locs, dims = 1)
s = size(unique_locs)[1]

## get unique locations per individual 
loc_idx = ones(noInd)
for i in 1:noInd 
    loc_idx[i] = only((findall(Bool[[obs_locs[i, :]] == [unique_locs[j, :]] for j=1:size(unique_locs, 1)])))
end
loc_idx = Int.(loc_idx)

###
### Locations that we'll predict on
###

lonPts = collect(range(minimum(lon), maximum(lon), length = 150))
latPts = collect(range(minimum(lat), maximum(lat), length = 150))
gridVec = vec(collect(Base.product(lonPts, latPts)))
gridLength = size(gridVec)[1]
gridVecDF = zeros(gridLength, 2)
for i in 1:gridLength
    gridVecDF[i, 1] = gridVec[i][1] # longitude
    gridVecDF[i, 2] = gridVec[i][2] # latitude
end

###
### Elevation 
###

lonmin = floor(minimum(lon)-1)
lonmax = ceil(maximum(lon)+1)
latmin = floor(minimum(lat)-1)
latmax = ceil(maximum(lat)+1)
G = grdcut("@earth_relief_30s", limits = (lonmin, lonmax, latmin, latmax))
elev = zeros(noInd)
println("\nObtain elevation data for observed locations.")
for i in ProgressBar(1:noInd) 
    elev[i] = grdtrack([lon[i] lat[i]], G)[3]
end

elevPts = zeros(gridLength)
println("\nObtain elevation data for gridded locations.")
for i in ProgressBar(1:gridLength)
    elevPts[i] = grdtrack([gridVecDF[i, 1] gridVecDF[i, 2]], G)[3]
end

## keep gridded points that are not absurdly deep (likely in the ocean)
aboveIdx = zeros(0)
for i in 1:gridLength 
    if elevPts[i] >= -30 # Caspian Sea is -29m
        append!(aboveIdx, i)
    end
end
aboveIdx = Int.(aboveIdx)
resetIdx = findall(x -> x <= -30, elev)
elev[resetIdx] .= 0

###
### Get Shortest Distance to Modern Lakes & Reservoirs
###

## get coordinates of naturally-occuring lake's centroids
R"""
library(foreign)
EL_Labels = read.dbf("Data/ne_10m_lakes_europe/ne_10m_lakes_europe.dbf")[9]
library(rgdal)
EL = readOGR("Data/ne_10m_lakes_europe/ne_10m_lakes_europe.shp")
EL_coords = coordinates(EL)[which(EL_Labels != "Reservoir"), ] # remove reservoirs
"""
@rget EL_coords

## get dist2water for samples
dist2water = zeros(noInd)
for i in 1:noInd 
    dist2water[i] = 999999999
    for j in 1:size(EL_coords)[1]
        currDist = euclidean([lon[i] lat[i]], EL_coords[j, :])
        if(currDist < dist2water[i])
            dist2water[i] = copy(currDist)
        end
    end
end

## get dist2water for grid
dist2waterPts = zeros(gridLength)
for i in 1:gridLength 
    dist2waterPts[i] = 999999999
    for j in 1:size(EL_coords)[1]
        currDist = euclidean([gridVecDF[i, 1] gridVecDF[i, 2]], EL_coords[j, :])
        if(currDist < dist2waterPts[i])
            dist2waterPts[i] = copy(currDist)
        end
    end
end

## construct design matrices
xy = [lon lat elev dist2water]
xyg = [gridVecDF[aboveIdx, 1] gridVecDF[aboveIdx, 2] elevPts[aboveIdx] dist2waterPts[aboveIdx]] # only keep points above -30m in altitude
xyc = [xy; xyg]

###
### Scale covariates together
###

@rput xyc
R"""
xycsc = xyc
for(i in 1:dim(xyc)[2]){
    xycsc[, i] = scale(xyc[, i])
}
lonScale = unlist(attributes(scale(xyc[, 1]))[3])
lonCenter = unlist(attributes(scale(xyc[, 1]))[2])
latScale = unlist(attributes(scale(xyc[, 2]))[3])
latCenter = unlist(attributes(scale(xyc[, 2]))[2])
elevScale = unlist(attributes(scale(xyc[, 3]))[3])
elevCenter = unlist(attributes(scale(xyc[, 3]))[2])
dist2waterScale = unlist(attributes(scale(xyc[, 4]))[3])
dist2waterCenter = unlist(attributes(scale(xyc[, 4]))[2])
"""
@rget xycsc lonScale lonCenter latScale latCenter elevScale elevCenter dist2waterScale dist2waterCenter

## extract scaled design matrices
xysc = xycsc[1:noInd, :] # scaled design matrix for sampled locations
xygsc = xycsc[(noInd + 1):end, :] # scaled design matrix for simulated locations; for use in plotting.jl

###
### Native Geodesic Distance (Haversine) for Julia
###

function geodesic_dist(lon1, lat1, lon2, lat2)
    R = 6371.0 # Earth's radius in kilometers

    ## convert degrees to radians
    φ1, λ1 = deg2rad(lat1), deg2rad(lon1)
    φ2, λ2 = deg2rad(lat2), deg2rad(lon2)

    ## Haversine formula
    Δφ = φ2 - φ1
    Δλ = λ2 - λ1
    a = sin(Δφ / 2)^2 + cos(φ1) * cos(φ2) * sin(Δλ / 2)^2
    c = 2 * atan(sqrt(a), sqrt(1 - a))

    distance = R * c # Distance in kilometers
    return distance
end

###
### Get Regression Variables 
###

n = noInd # number of individuals
N = Int(n*(n-1)/2) # number of connections
p = 8 # number of parameters + intercept

H = zeros(N, n) # theta mapping matrix
K = zeros(N, s) # eta mapping matrix
X = ones(N, p + 1) # design matrix
Ygd = zeros(N) # vector of genetic distances
dist = zeros(N) # vector of geodesic distances
dt = zeros(N) # vector of distances in time
idxiVec = zeros(N) # vector of ith node for each connection
idxjVec = zeros(N) # vector of jth node for each connection

idxi = 1 # ith index
idxj = 2 # jth index
println("Structing Data for Dyadic Analysis: ")
for i in ProgressBar(1:N)

    ## y vector
    Ygd[i] = copy(GDM[idxi, idxj])
    
    ## node vectors
    idxiVec[i] = idxi
    idxjVec[i] = idxj

    ## for phi list
    dist[i] = geodesic_dist(xy[idxi, 1], xy[idxi, 2], xy[idxj, 1], xy[idxj, 2])

    ## for weights
    dt[i] = abs(times[idxi] - times[idxj])

    ## X matrix
    if times[idxj] > times[idxi] # the greater the years BP, the further back in time; if i is most recent
        X[i, 2] = xysc[idxi, 1] - xysc[idxj, 1] # difference in longitude
        X[i, 3] = xysc[idxi, 2] - xysc[idxj, 2] # difference in latitude
        X[i, 4] = xysc[idxi, 3] - xysc[idxj, 3] # difference in elevation
        X[i, 5] = xysc[idxi, 4] - xysc[idxj, 4] # difference in proximity to water
        X[i, 6] = xysc[idxi, 1]^2 - xysc[idxj, 1]^2 # difference in longitude squared
        X[i, 7] = xysc[idxi, 2]^2 - xysc[idxj, 2]^2 # difference in latitude squared
        X[i, 8] = xysc[idxi, 1]*xysc[idxi, 2] - xysc[idxj, 1]*xysc[idxj, 2] # difference in longitude*latitude
        X[i, 9] = times[idxi] - times[idxj]
    else # if j is most recent
        X[i, 2] = xysc[idxj, 1] - xysc[idxi, 1] # difference in longitude
        X[i, 3] = xysc[idxj, 2] - xysc[idxi, 2] # difference in latitude
        X[i, 4] = xysc[idxj, 3] - xysc[idxi, 3] # difference in elevation
        X[i, 5] = xysc[idxj, 4] - xysc[idxi, 4] # difference in proximity to water
        X[i, 6] = xysc[idxj, 1]^2 - xysc[idxi, 1]^2 # difference in longitude squared
        X[i, 7] = xysc[idxj, 2]^2 - xysc[idxi, 2]^2 # difference in latitude squared
        X[i, 8] = xysc[idxj, 1]*xysc[idxj, 2] - xysc[idxi, 1]*xysc[idxi, 2] # difference in longitude*latitude
        X[i, 9] = times[idxj] - times[idxi]
    end

    ## mapping matrix H ; added regardless of order to induce dyadic dependence
    if times[idxj] > times[idxi] # if i is most recent
        H[i, idxi] = 1
        H[i, idxj] = 1
    else 
        H[i, idxi] = 1
        H[i, idxj] = 1
    end

    # new mapping matrix K
    if loc_idx[idxi] == loc_idx[idxj] # if the nodes are co-located
        K[i, loc_idx[idxi]] = 0 # the etas cancel out and we get nothing
    else # if the nodes are in different locations
        if times[idxj] > times[idxi] # if node i is most recent
            K[i, loc_idx[idxi]] = 1
            K[i, loc_idx[idxj]] = -1
        else # if node j is the most recent
            K[i, loc_idx[idxi]] = -1
            K[i, loc_idx[idxj]] = 1
        end
    end
    
    ## update individual indices
    if idxj == noInd
        global idxi += 1
        global idxj = idxi + 1
    else
        global idxj += 1
    end
end
X = X[:, 2:end]

###
### Pre-compute ϕ and R
###

## pre-compute phi
phi = 0:10:(maximum(dist)/3)

println("ϕ range: ", phi)
println("Length of phi: ", length(phi))
println("Maximum distance: ", maximum(dist))

## pre-compute R
println("Computing R(ϕ) for all ϕ: ")
R_list = ones(s, s, length(phi))
for k in ProgressBar(1:length(phi))
    for i in 1:(s-1) # for each unique location
        for j in (i+1):s # for each unique location
            tmpDist = geodesic_dist(unique_locs[i, 1], unique_locs[i, 2], unique_locs[j, 1], unique_locs[j, 2])
            R_list[i, j, k] = exp(-tmpDist/phi[k])
            R_list[j, i, k] = R_list[i, j, k]
        end
    end
    R_list[:, :, k] = R_list[:, :, k] + 0.0001*I
end

###
### For covariance weights 
###

@rput dt dist 
R"""
library(scales)
dtij = rescale(dt, to = c(0.00001, 1))
dsij = rescale(dist, to = c(0.00001, 1))
"""
@rget dtij dsij

###
### Set Directory for Output 
###

cd("$stringDate")

###
### Save environment
###

@save "juliaVars.jld2"

###
### Run MCMC Algorithm
###

println("\n%%%%%%%%%%%%%%%%%%%%%%%%%%% Running MCMC %%%%%%%%%%%%%%%%%%%%%%%%%%%")
startTime = Dates.Time(Dates.now())

println("MCMC Algorithm (Both):")
out4 = mcmc_both(Ygd, X, H, K, phi, R_list, nMCMC, dtij, dsij)

println("MCMC Algorithm (Neither):")
out = mcmc_none(Ygd, X, H, K, phi, R_list, nMCMC, dtij, dsij)

println("MCMC Algorithm (dt):")
out2 = mcmc_dt(Ygd, X, H, K, phi, R_list, nMCMC, dtij, dsij)

println("MCMC Algorithm (ds):")
out3 = mcmc_ds(Ygd, X, H, K, phi, R_list, nMCMC, dtij, dsij)

println("Time took: ", Dates.Time(Dates.now()) - startTime, ".")