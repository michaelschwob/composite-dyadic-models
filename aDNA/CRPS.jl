##############################
####
#### Script for Computing CRPS
####
##############################

println("\n%%%%%%%%%%%%%%%%%%%%%%%%%%% Plotting & Computing CRPS %%%%%%%%%%%%%%%%%%%%%%%%%%%")

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
using PlotlyJS, DataFrames # for boxplots

println("Time started: ", Dates.Time(Dates.now()), ".")

###
### Load Global Environment
###

@load "juliaVars.jld2"
Hsp = sparse(H)
Ksp = sparse(K)

###
### CRPS Function 
###

function myCRPS(y, mu, sd)
    ω = (y-mu)/sd
    snCDF = cdf(Normal(0, 1), ω)
    snPDF = pdf(Normal(0, 1), ω)
    sd*(ω*(2*snCDF - 1) + 2*snPDF - π^(-1/2))
end

###
### Initialize Weights Matrix and Grid
###

weights = [0 0]

R"""
seqList = seq(0.0001, 1, length.out = 1000)
gridLocs = expand.grid(seqList, seqList)
"""
@rget gridLocs
wGrid = zeros(size(gridLocs)[1], 4) # will contain the estimated downweight

###
### For Each Model 
###

saved = zeros(4)

for i in 1:4

    ###
    ### Load MCMC Output 
    ###

    if i == 1
        println("\n(Neither)")
        @load "none/MCMCoutput.jld2"
    elseif i == 2
        println("(dt)")
        @load "dt/MCMCoutput.jld2"
    elseif i == 3
        println("(ds)")
        @load "ds/MCMCoutput.jld2"
    else 
        println("(Both)")
        @load "both/MCMCoutput.jld2"
    end

    ###
    ### Compute Posterior Means
    ###

    etaMeans = mean(etaSave, dims = 2)
    thetaMeans = mean(thetaSave, dims = 2)
    betaMeans = mean(betaSave, dims = 2)
    s2_yMean = mean(s2_ySave)
    gammaMean = mean(gammaSave, dims = 2)

    ###
    ### Calculate CRPS
    ###

    println("\nCalculating CRPS.")

    if i == 1 # without weights
        mu = X*betaMeans + K*etaMeans + H*thetaMeans

        post_weights = ones(N)
        post_sd = sqrt.(s2_yMean./post_weights)

        testCRPS = zeros(N)
        [testCRPS[i] = myCRPS(Ygd[i], mu[i], post_sd[i]) for i in 1:N]
        crps = mean(testCRPS)
    else # with weights
        mu = X*betaMeans + K*etaMeans + H*thetaMeans

        post_weights = exp.(-gammaMean[1].*dtij - gammaMean[2].*dsij)
        post_sd = sqrt.(s2_yMean./post_weights)

        testCRPS = zeros(N)
        [testCRPS[i] = myCRPS(Ygd[i], mu[i], post_sd[i]) for i in 1:N]
        crps = mean(testCRPS)
    end

    roundedCRPS = round(crps, digits = 6)
    saved[i] = copy(roundedCRPS)

    ###
    ### Compute Weights
    ###

    post_weights = exp.(-gammaMean[1].*dtij - gammaMean[2].*dsij)
    if i == 1
        tmpID = ["No Weights <br>(CRPS = $roundedCRPS)" for i in 1:length(post_weights)]
    elseif i == 2
        tmpID = ["Weights: dt <br>(CRPS = $roundedCRPS)" for i in 1:length(post_weights)]
    elseif i == 3
        tmpID = ["Weights: ds <br>(CRPS = $roundedCRPS)" for i in 1:length(post_weights)]
    else 
        tmpID = ["Weights: dt & ds <br>(CRPS = $roundedCRPS)" for i in 1:length(post_weights)]
    end
    global weights = vcat(weights, [post_weights tmpID])

    ###
    ### Gridded Weights
    ###

    wGrid[:, i] = exp.(-gammaMean[1].*gridLocs[:, 1] - gammaMean[2].*gridLocs[:, 2])

end

###
### Structure Weights Data Frame for Box plots
###

weights = weights[2:end, :] # get rid of initial row
weightsDF = DataFrame(weights, ["Weights", "Model"])

###
### Grouped Boxplots
###

BP = PlotlyJS.plot(weightsDF, x = :Model, y = :Weights, kind = "box", Layout(title = attr(text = "Model Comparison", x = 0.5, y = 0.95), yaxis_title = "Distribution of Weights", template = "simple_white"))
PlotlyJS.savefig(BP, "weights.png")
BP_wop = PlotlyJS.plot(weightsDF, x = :Model, y = :Weights, kind = "box", Layout(title = attr(text = "Model Comparison", x = 0.5, y = 0.95), yaxis_title = "Distribution of Weights", template = "simple_white"), boxpoints = false) # exclude outliers
PlotlyJS.savefig(BP_wop, "weights_wop.png")