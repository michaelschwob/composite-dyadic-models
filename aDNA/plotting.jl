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
using Distances # compute distances
using MCMCChains, MCMCDiagnosticTools, StatsPlots, Plots # for MCMC chain management and plotting
import PlotlyJS # for contours
using Mads # for kriging
Random.seed!(512)

println("Time started: ", Dates.Time(Dates.now()), ".")
parentPath = pwd()

###
### Load MCMC Output
###

println("Reading Environments.")
@load "juliaVars.jld2"

###
### R trace plot function
###

R"""
library(ggplot2)
mrs.trace <- function(chain, par_name, thin = 1, burn = floor(length(chain)*0.1)){
  new.chain <- chain[seq(burn, length(chain), by=thin)]

  if(missing(par_name)){
    title <- "Trace Plot"
    par <- "Parameter"
  }else{
    title <- "Trace Plot"
    par <- parse(text=par_name)
  }

  dfmrs.trace <- data.frame(cbind(new.chain = new.chain, iter = 1:length(new.chain)))

  plot <- ggplot(dfmrs.trace, aes(x = iter, y = new.chain)) + theme_classic() + geom_line() + xlab("Iterations") + ylab(par) + ggtitle(title) + theme(plot.title = element_text(hjust = 0.5))
  return(plot)
}
"""

###
### For each model
###

for m in 1:4

    cd(parentPath)
    @rput parentPath
    R"""
    setwd(parentPath)
    """

    ###
    ### Load MCMC Output 
    ###

    if m == 1
        println("\n(Neither)")
        @load "none/MCMCoutput.jld2"
        cd(parentPath*"/none")
        R"""
        setwd(paste0(parentPath, "/none"))
        """
    elseif m == 2
        println("\n(dt)")
        @load "dt/MCMCoutput.jld2"
        cd(parentPath*"/dt")
        R"""
        setwd(paste0(parentPath, "/dt"))
        """
    elseif m == 3
        println("\n(ds)")
        @load "ds/MCMCoutput.jld2"
        cd(parentPath*"/ds")
        R"""
        setwd(paste0(parentPath, "/ds"))
        """
    else 
        println("\n(Both)")
        @load "both/MCMCoutput.jld2"
        cd(parentPath*"/both")
        R"""
        setwd(paste0(parentPath, "/both"))
        """
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
    gammaChs = Chains(permutedims(outDict["gamma"], [2, 1]))
    s2_yCh = Chains(outDict["s2_y"])
    s2_etaCh = Chains(outDict["s2_eta"])
    s2_thetaCh = Chains(outDict["s2_theta"])
    phiCh = Chains(outDict["phi"])

    ## specifically for visualization
    etaChsPlot = Chains(permutedims(outDict["eta"][1:10, nBurn:k], [2, 1]), ["η"*string(i) for i in 1:10])
    thetaChsPlot = Chains(permutedims(outDict["theta"][1:10, nBurn:k], [2, 1]), ["θ"*string(i) for i in 1:10])
    betaChsPlot = Chains(permutedims(outDict["beta"][:, nBurn:k], [2, 1]), ["β"*string(i) for i in (1-1):(size(outDict["beta"])[1]-1)])
    gammaChsPlot = Chains(permutedims(outDict["gamma"][:, nBurn:k], [2, 1]), ["γ"*string(i) for i in 1:2])
    s2_yChPlot = Chains(outDict["s2_y"][nBurn:k], ["σ²(y)"])
    s2_etaChPlot = Chains(outDict["s2_eta"][nBurn:k], ["σ²(η)"])
    s2_thetaChPlot = Chains(outDict["s2_theta"][nBurn:k], ["σ²(θ)"])
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

    ## ϕ trace plot (in R as an example)
    phiCh = outDict["phi"]
    @rput phiCh
    R"""
    png("phiTrace.png")
    phi = mrs.trace(phiCh, "phi")
    print(phi)
    dev.off()
    """

    ###
    ### Plot Potential Surface
    ###

    println("Plotting potential surface.")

    ## construct grid data frame to plot
    gridDf = zeros(size(xygsc)[1], 3)
    for i in 1:size(xygsc)[1]
        gridDf[i, 1] = xygsc[i, 1] # longitude
        gridDf[i, 2] = xygsc[i, 2] # latitude
        gridDf[i, 3] = xygsc[i, 1]*betaMeans[2] + xygsc[i, 2]*betaMeans[3] + xygsc[i, 3]*betaMeans[4] + xygsc[i, 4]*betaMeans[5] + (xygsc[i, 1]^2)*betaMeans[6] + (xygsc[i, 2]^2)*betaMeans[7] + xygsc[i, 1]*xygsc[i, 2]*betaMeans[8]
    end

    ## scale lake coordinates
    lakesDF = copy(EL_coords)
    lakesDF[:, 1] = (EL_coords[:, 1] .- lonCenter)./lonScale
    lakesDF[:, 2] = (EL_coords[:, 2] .- latCenter)./latScale

    ## keep lake centroids within our study area
    inArea = zeros(0)
    for i in 1:size(lakesDF)[1]
        if minimum(xysc[:, 1]) <= lakesDF[i, 1] <= maximum(xysc[:, 1]) && minimum(xysc[:, 2]) <= lakesDF[i, 2] <= maximum(xysc[:, 2])
            append!(inArea, i)
        end
    end 
    lakesDF = lakesDF[Int.(inArea), :]

    ## plot potential surface
    @rput gridDf xysc lakesDF
    R"""
    library(ggplot2)
    potPlot = ggplot(as.data.frame(gridDf), aes(gridDf[, 1], gridDf[, 2], z = gridDf[, 3])) + geom_contour_filled() + geom_point(as.data.frame(xysc), mapping = aes(x = xysc[, 1], y = xysc[, 2], z = "black")) + labs(x = "Longitude", y = "Latitude", title = "Xβ Surface", fill = "Xβ") + theme_classic()
    ggsave("XbPlot.png", potPlot, width = 10, height = 8)
    potPlotLakes = potPlot + geom_point(data = as.data.frame(lakesDF), mapping = aes(x = lakesDF[, 1], y = lakesDF[, 2], z = "blue"), col = "blue")
    ggsave("XbPlotLakes.png", potPlotLakes, width = 10, height = 8)
    """

    ###
    ### Get Predictive Covariance matrices
    ###

    println("Kriging the η surface.")
    noNew = size(xygsc)[1]
    etaKrige = zeros(noNew, length((nMCMC-1000):nMCMC))
    println("Kriging.")
    for k in ProgressBar((nMCMC-1000):nMCMC)
        function covFunc(dist)
            return s2_etaSave[k]*exp(-dist/phiSave[k])
        end
        etaKrige[:, k - (nMCMC-1000) + 1] = Mads.krige(xygsc[:, 1:2]', unique_locs', etaSave[:, k], covFunc)
    end
    etaKrigePost = mean(etaKrige, dims = 2)
    etaDf = DataFrame(x = gridDf[:, 1], y = gridDf[:, 2], z = vec(etaKrigePost))

    ## plot in R 
    @rput etaDf
    R"""
    etaPlot = ggplot(etaDf, aes(x, y, z = z)) + geom_contour_filled() + geom_point(as.data.frame(xysc), mapping = aes(x = xysc[, 1], y = xysc[, 2], z = "black")) + labs(x = "Longitude", y = "Latitude", title = "η Surface", fill = "η") + theme_classic()
    ggsave("etaPlot.png", etaPlot, width = 10, height = 8)
    """

    ###
    ### Plot η and potential together
    ###

    Epi = vec(etaKrigePost) + gridDf[:, 3]
    togDF = DataFrame(x = gridDf[:, 1], y = gridDf[:, 2], z = Epi)

    ## plot in R 
    @rput togDF
    R"""
    togPlot = ggplot(togDF, aes(x, y, z = z)) + geom_contour_filled() + geom_point(as.data.frame(xysc), mapping = aes(x = xysc[, 1], y = xysc[, 2], z = "black")) + labs(x = "", y = "", title = "Potential Surface", fill = "ρ") + theme_classic()
    ggsave("potPlot.png", togPlot, width = 10, height = 8)
    """

    ###
    ### Get Standard Deviation Across Potential Surface
    ###

    ## construct grid data frame to plot
    sds = zeros(size(xygsc)[1])
    println("Standard Deviation Computations:")
    for i in ProgressBar(1:size(xygsc)[1])
        xstar = [xygsc[i, 1] xygsc[i, 2] xygsc[i, 3] xygsc[i, 4] (xygsc[i, 1]^2) (xygsc[i, 2]^2) xygsc[i, 1]*xygsc[i, 2]]
        tmpsum = 0
        for k in (nMCMC-1000):nMCMC
            tmpsum += ((xstar*outDict["beta"][2:8, k])[1] + etaKrige[i, k - (nMCMC-1000) + 1] - Epi[i])^2
        end
        sds[i] = sqrt(tmpsum/1000)
    end
    sdsDF = DataFrame(x = gridDf[:, 1], y = gridDf[:, 2], z = sds)

    ## plot in R 
    @rput sdsDF
    R"""
    sdsPlot = ggplot(sdsDF, aes(x, y, z = z)) + geom_contour_filled() + geom_point(as.data.frame(xysc), mapping = aes(x = xysc[, 1], y = xysc[, 2], z = "black")) + labs(x = "", y = "", title = "Potential Surface", fill = "Standard Deviation") + theme_classic()
    ggsave("sdsPlot.png", sdsPlot, width = 10, height = 8)
    """

end

## re-set directory
cd(parentPath)