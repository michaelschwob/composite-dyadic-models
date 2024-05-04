#################################################
####
#### Script for Plotting Dyadic Regression Output
####
#################################################

println("\n%%%%%%%%%%%%%%%%%%%%%%%%%%% Plotting w/ MCMC Output %%%%%%%%%%%%%%%%%%%%%%%%%%%")

###
### Set-up
###

using LinearAlgebra, Random, Statistics, Distributions, SparseArrays # stats
using RCall # R
using CSV, DataFrames, Query # data
using Dates, ProgressBars # progress
using JLD2 # environments
using Distances
using MCMCChains, MCMCDiagnosticTools, StatsPlots, Plots # for chains
import PlotlyJS # for contours
using Mads # for kriging
Random.seed!(702)

println("Time started: ", Dates.Time(Dates.now()), ".")
parentPath = pwd()

###
### Load MCMC Output
###

println("Reading Environments.")
@load "juliaVars.jld2"

###
### For each model
###

for m in 1:2

    cd(parentPath)

    ###
    ### Load MCMC Output 
    ###

    if m == 1
        println("\n(Neither)")
        @load "none/MCMCoutput.jld2"
        cd(parentPath*"/none")
    else 
        println("\n(ds)")
        @load "ds/MCMCoutput.jld2"
        cd(parentPath*"/ds")
    end

    ###
    ### Construct Dictionary
    ###

    outDict = Dict("beta" => betaSave, "theta" => thetaSave, "eta" => etaSave, "s2_y" => s2_ySave, "s2_eta" => s2_etaSave, "s2_theta" => s2_thetaSave, "gamma" => gammaSave, "phi" => phiSave)

    ###
    ### Compute Posterior Means
    ###

    etaMeans = mean(outDict["eta"], dims = 2)
    thetaMeans = mean(outDict["theta"], dims = 2)
    betaMeans = mean(outDict["beta"], dims = 2)
    s2_yMean = mean(outDict["s2_y"])
    s2_etaMean = mean(outDict["s2_eta"])
    s2_thetaMean = mean(outDict["s2_theta"])
    gammaMean = mean(outDict["gamma"], dims = 2)
    phiMean = mean(outDict["phi"])

    ###
    ### Construct Chains
    ###

    println("Plotting trace plots.")

    nBurn = Int(0.2*k)

    etaChs = Chains(permutedims(outDict["eta"], [2, 1]))
    betaChs = Chains(permutedims(outDict["beta"], [2, 1]))
    thetaChs = Chains(permutedims(outDict["theta"], [2, 1]))
    s2_yCh = Chains(outDict["s2_y"])
    s2_etaCh = Chains(outDict["s2_eta"])
    s2_thetaCh = Chains(outDict["s2_theta"])
    gammaChs = Chains(permutedims(outDict["gamma"], [2, 1]))
    phiCh = Chains(outDict["phi"])

    ## specifically for visualization
    etaChsPlot = Chains(permutedims(outDict["eta"][1:10, nBurn:k], [2, 1]), ["η"*string(i) for i in 1:10])
    thetaChsPlot = Chains(permutedims(outDict["theta"][1:10, nBurn:k], [2, 1]), ["θ"*string(i) for i in 1:10])
    betaChsPlot = Chains(permutedims(outDict["beta"][:, nBurn:k], [2, 1]), ["β"*string(i) for i in (1-1):(size(outDict["beta"])[1]-1)])
    s2_yChPlot = Chains(outDict["s2_y"][nBurn:k], ["σ²(y)"])
    s2_etaChPlot = Chains(outDict["s2_eta"][nBurn:k], ["σ²(η)"])
    s2_thetaChPlot = Chains(outDict["s2_theta"][nBurn:k], ["σ²(θ)"])
    gammaChsPlot = Chains(permutedims(outDict["gamma"][:, nBurn:k], [2, 1]), ["γ"*string(i) for i in 1:2])
    phiChPlot = Chains(outDict["phi"][nBurn:k], ["ϕ"])

    ###
    ### Plot Traces
    ###

    savefig(Plots.plot(etaChsPlot), "etaChs.png")
    savefig(Plots.plot(thetaChsPlot), "thetaChs.png")
    savefig(Plots.plot(betaChsPlot), "betaChs.png")
    savefig(Plots.plot(s2_yChPlot), "s2_yCh.png")
    savefig(Plots.plot(s2_etaChPlot), "s2_etaCh.png")
    savefig(Plots.plot(s2_thetaChPlot), "s2_thetaCh.png")
    savefig(Plots.plot(gammaChsPlot), "gammaChs.png")
    savefig(Plots.plot(phiChPlot), "phiCh.png")

    ###
    ### Plot Potential Surface
    ###

    println("Plotting potential surface.")

    ## construct grid data frame to plot
    gridDf = zeros(size(gridPred)[1], 3)
    for i in 1:size(gridPred)[1]
        gridDf[i, 1] = gridPred[i, 1] # longitude
        gridDf[i, 2] = gridPred[i, 2] # latitude
        gridDf[i, 3] = gridCovs[i, 1]*betaMeans[2] + gridCovs[i, 2]*betaMeans[3] # potential function effect
    end

    ## plot potential surface
    @rput gridDf rupicaLocs
    R"""
    library(ggplot2)
    potPlot = ggplot(as.data.frame(gridDf), aes(gridDf[, 1], gridDf[, 2], z = gridDf[, 3])) + geom_contour_filled() + geom_point(as.data.frame(rupicaLocs), mapping = aes(x = rupicaLocs[, 1], y = rupicaLocs[, 2], z = "black")) + labs(x = "Longitude", y = "Latitude", title = "Xβ Surface", fill = "Xβ") + theme_classic()
    ggsave("XbPlot.png", potPlot)
    """

    ###
    ### Write Covariance Function
    ###

    noNew = size(gridPred)[1]
    etaKrige = zeros(noNew, length((nMCMC-1000):nMCMC))

    println("Kriging.")
    for k in ProgressBar((nMCMC-1000):nMCMC)
        function covFunc(dist)
            return s2_etaSave[k]*exp(-dist/phiSave[k])
        end
        etaKrige[:, k - (nMCMC-1000) + 1] = Mads.krige(gridPred[:, 1:2]', unique_locs', etaSave[:, k], covFunc)
    end
    etaKrigePost = mean(etaKrige, dims = 2)

    etaDf = DataFrame(x = gridDf[:, 1], y = gridDf[:, 2], z = vec(etaKrigePost))

    ## plot in R 
    @rput etaDf
    R"""
    etaPlot = ggplot(etaDf, aes(x, y, z = z)) + geom_contour_filled() + geom_point(as.data.frame(rupicaLocs), mapping = aes(x = rupicaLocs[, 1], y = rupicaLocs[, 2], z = "black")) + labs(x = "Longitude", y = "Latitude", title = "η Surface", fill = "η") + theme_classic()
    ggsave("etaPlot.png", etaPlot)
    """

    ###
    ### Plot η and potential together
    ###

    togDF = DataFrame(x = gridDf[:, 1], y = gridDf[:, 2], z = vec(etaKrigePost) + gridDf[:, 3])

    ## plot in R 
    @rput togDF
    R"""
    togPlot = ggplot(togDF, aes(x, y, z = z)) + geom_contour_filled() + geom_point(as.data.frame(rupicaLocs), mapping = aes(x = rupicaLocs[, 1], y = rupicaLocs[, 2], z = "black")) + labs(x = "Longitude", y = "Latitude", title = "Potential Surface", fill = "ρ") + theme_classic()
    ggsave("potPlot.png", togPlot)
    """

    ###
    ### Plot elevation by potential
    ###

    elev = sort(gridPred[:, 3])
    minElev = elev[1]
    maxElev = elev[end]
    elevSeq = minElev:1:maxElev
    elevSeqSC = (elevSeq .- e1Center)./e1Scale
    elev2SeqSC = ((elevSeq.^2) .- e2Center)./e2Scale
    potSeq = elevSeqSC.*betaMeans[2] + elev2SeqSC.*betaMeans[3]
    @rput elevSeq potSeq
    R"""
    png("elevPot.png")
    plot(elevSeq, potSeq, main = "Posterior Potential per Elevation", type = "l", xlab = "Elevation (m)", ylab = expression(rho))
    dev.off()
    """

end

## switch back to parent directory
cd(parentPath)