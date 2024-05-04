########################################################
####
#### Script for AR(1) Example w/ Centered-Mean Structure
####
########################################################

###
### Set-up
###

set.seed(702)
library(mvtnorm)
library(LaplacesDemon)
library(ggplot2)
library(progress)
library(dplyr)

###
### MCMC Algorithm
###

## standard AR(1) implementation
mcmc = function(y, X, n.mcmc){
  
  ###
  ### Set-up
  ###
  
  T = length(y)
  p = dim(X)[2]
  
  Sig_beta_inv = solve(10^6*diag(p))
  s2.alpha = 0.01
  s2.beta = 0.01
  
  beta = rep(0, p)
  s2 = 1
  
  betaSave = matrix(0, n.mcmc, p)
  betaSave[1, ] = beta 
  s2Save = rep(0, n.mcmc)
  s2Save[1] = s2
  
  ###
  ### Gibbs Sampler
  ###
  
  pb1 = progress_bar$new(format = paste0("Correct AR(1) MCMC", " [:bar] :percent eta: :eta"), total = n.mcmc, clear = FALSE)
  for(k in 2:n.mcmc){
    
    pb1$tick()
    
    ###
    ### beta update
    ###
    
    sum1 = X[1, ]%*%t(X[1, ])
    sum2 = y[1]*X[1, ]
    
    for(t in 2:T){
      sum1 = sum1 + (X[t, ] - X[t-1, ])%*%t(X[t, ] - X[t-1, ])
      sum2 = sum2 + (y[t] - y[t-1])*(X[t, ] - X[t-1, ])
    }
    
    A_inv = solve(sum1/s2 + Sig_beta_inv)
    b = sum2/s2
    beta = rmvnorm(1, A_inv%*%b, A_inv)
    
    ###
    ### s2 update 
    ###
    
    sum3 = (y[1] - X[1, ]%*%t(beta))^2
    for(t in 2:T){
      sum3 = sum3 + (y[t] - y[t-1] - (X[t, ] - X[t-1, ])%*%t(beta))^2
    }  
    
    alpha.star = s2.alpha + T/2
    beta.star = s2.beta + sum3/2
    s2 = rinvgamma(1, shape = alpha.star, scale = 1/beta.star)
    
    ###
    ### Save updates
    ###
    
    betaSave[k, ] = beta
    s2Save[k] = s2
    
  }
  
  return(list(betaSave = betaSave, s2Save = s2Save))
}

## (non-weighted) dyadic implementation
mcmcDyad = function(yij, Xij, n.mcmc){
  
  ###
  ### Set-up
  ###
  
  T = length(yij)
  p = dim(Xij)[2]
  
  Sig_beta_inv = solve(10^6*diag(p))
  s2.alpha = 0.01
  s2.beta = 0.01
  
  beta = rep(0, p)
  s2 = 1
  
  betaSave = matrix(0, n.mcmc, p)
  betaSave[1, ] = beta 
  s2Save = rep(0, n.mcmc)
  s2Save[1] = s2
  
  ###
  ### Gibbs Sampler
  ###
  
  pb2 = progress_bar$new(format = paste0("MCMC w/o Composite Weights", " [:bar] :percent eta: :eta"), total = n.mcmc, clear = FALSE)
  for(k in 2:n.mcmc){
    
    pb2$tick()
    
    ###
    ### beta update
    ###
    
    sum1 = Xij[1, ]%*%t(Xij[1, ])
    sum2 = yij[1]*Xij[1, ]
    
    for(t in 2:T){
      sum1 = sum1 + Xij[t, ]%*%t(Xij[t, ])
      sum2 = sum2 + yij[t]*Xij[t, ]
    }
    
    A_inv = solve(sum1/s2 + Sig_beta_inv)
    b = sum2/s2
    beta = rmvnorm(1, A_inv%*%b, A_inv)
    
    ###
    ### s2 update 
    ###
    
    sum3 = (yij[1] - Xij[1, ]%*%t(beta))^2
    for(t in 2:T){
      sum3 = sum3 + (yij[t] - Xij[t, ]%*%t(beta))^2
    }  
    
    alpha.star = s2.alpha + T/2
    beta.star = s2.beta + sum3/2
    s2 = rinvgamma(1, shape = alpha.star, scale = 1/beta.star)
    
    ###
    ### Save updates
    ###
    
    betaSave[k, ] = beta
    s2Save[k] = s2
    
  }
  
  return(list(betaSave = betaSave, s2Save = s2Save))

}

## (weighted) dyadic implementation
mcmcWeighted = function(yij, Xij, nudt, n.mcmc){
  
  ###
  ### Set-up
  ###
  
  T = length(yij)
  p = dim(Xij)[2]
  p.g = dim(nudt)[2]
  
  Sig_beta_inv = solve(10^6*diag(p))
  s2.alpha = 0.01
  s2.beta = 0.01
  gamma.alpha = 10
  gamma.beta = 2
  
  beta = rep(0, p)
  s2 = 1
  gamma = rep(1, p.g)
  gij = exp(nudt%*%gamma)
  
  betaSave = matrix(0, n.mcmc, p)
  betaSave[1, ] = beta 
  s2Save = rep(0, n.mcmc)
  s2Save[1] = s2
  gammaSave = matrix(0, n.mcmc, p.g) 
  gammaSave[1, ] = gamma
  
  ls_gamma = rep(0, p.g) # from Roberts & Rosenthal (2009)
  M = 5 #  global maximal parameter value for ls
  accGamma = rep(0, p.g) # the number of accepted gammas in the last 50 iterations
  nBatch = 0 # the number of 50-iteration batches that we've processed
  
  ###
  ### Gibbs Sampler
  ###
  
  pb3 = progress_bar$new(format = paste0("MCMC w/ Composite Weights", " [:bar] :percent eta: :eta"), total = n.mcmc, clear = FALSE)
  for(k in 2:n.mcmc){
    
    pb3$tick()
    
    ###
    ### beta update
    ###
    
    sum1 = Xij[1, ]%*%t(Xij[1, ])/gij[1]
    sum2 = yij[1]*Xij[1, ]/gij[1]
    
    for(t in 2:T){
      sum1 = sum1 + Xij[t, ]%*%t(Xij[t, ])/gij[t]
      sum2 = sum2 + yij[t]*Xij[t, ]/gij[t]
    }
    
    A_inv = solve(sum1/s2 + Sig_beta_inv)
    b = sum2/s2
    beta = rmvnorm(1, A_inv%*%b, A_inv)
    
    ###
    ### s2 update 
    ###
    
    sum3 = (yij[1] - Xij[1, ]%*%t(beta))^2/gij[1]
    for(t in 2:T){
      sum3 = sum3 + (yij[t] - Xij[t, ]%*%t(beta))^2/gij[t]
    }  
    
    alpha.star = s2.alpha + T/2
    beta.star = s2.beta + sum3/2
    s2 = rinvgamma(1, shape = alpha.star, scale = 1/beta.star)
    
    ###
    ### gamma update
    ###
    
    ## update number of batches
    if(k%%50 == 0){ # every 50 iterations
      nBatch = nBatch + 1 # count up number of 50-iteration batches
      delta_nBatch = min(c(0.01, nBatch^(-0.5))) # the modification to ls
    }
    
    for(i in 1:p.g){
      ## adaptive tuning from Roberts & Rosenthal (2009)
      if(k%%50 == 0){ # every 50 iterations
        if(accGamma[i]/50 > 0.44){ # if we are accepting too many proposals...
          if((ls_gamma[i] + delta_nBatch) < M){ # ...and if the potential modification is within bounds,
            ls_gamma[i] = ls_gamma[i] + delta_nBatch # then increase ls
          }
        }else{ # if we are accepting too few proposals...
          if((ls_gamma[i] - delta_nBatch) > (-M)){ # ...and if the potential modification is within bounds,
            ls_gamma[i] = ls_gamma[i] - delta_nBatch # then decrease ls
          }
        }
        accGamma[i] = 0 # reset the count of accepted gammas
      }
      
      ## obtain valid proposal for gamma
      gammaStar_i = rnorm(1, gamma[i], exp(ls_gamma[i]))
      while(gammaStar_i < 0){ # restrict to proper support
        gammaStar_i = rnorm(1, gamma[i], exp(ls_gamma[i]))
      }
      
      ## compute variance scaler
      gammaStar = gamma
      gammaStar[i] = gammaStar_i
      gijStar = exp(nudt%*%gammaStar)
      
      ## compute M-H ratio
      mh1 = sum(dnorm(yij, Xij%*%t(beta), sqrt(s2*gijStar), log = TRUE)) + dgamma(gammaStar_i, gamma.alpha, gamma.beta, log = TRUE)
      mh2 = sum(dnorm(yij, Xij%*%t(beta), sqrt(s2*gij), log = TRUE)) + dgamma(gamma[i], gamma.alpha, gamma.beta, log = TRUE)
      mh = exp(mh1 - mh2)
      
      if(mh > runif(1)){
        accGamma[i] = accGamma[i] + 1
        gamma[i] = gammaStar_i
        gij = gijStar
      }
    }
    
    ###
    ### Save updates
    ###
    
    betaSave[k, ] = beta
    s2Save[k] = s2
    gammaSave[k, ] = gamma
    
  }
  
  return(list(betaSave = betaSave, s2Save = s2Save, gammaSave = gammaSave))
  
}

###
### Simulate Data
###

T = 15
p = 2
beta_true = c(1.3, 0.8)
s2_true = 0.01
X = matrix(rnorm(T*p, 0, 3), T, p) # random initialization
for(t in 2:T){
  X[t, 1] = rnorm(1, X[t-1, 1], sd = 0.5)
  X[t, 2] = rnorm(1, X[t-1, 2], sd = 0.5)
}
plot(1:T, X[, 1], xlab = "Time", ylab = "X[,1]", type = "l")
plot(1:T, X[, 2], xlab = "Time", ylab = "X[,2]", type = "l")

y = rnorm(1, mean = X[1, ]%*%beta_true, sd = sqrt(s2_true))
for(t in 2:T){
  y = c(y, rnorm(1, mean = y[t-1] + (X[t, ] - X[t-1, ])%*%beta_true, sd = sqrt(s2_true)))
}
plot(1:T, y, type = "l")

###
### Convert to (Correct) Dyadic Data 
###

N_t = T-1
yij_t = c(y[1], y[2:T] - y[1:(T-1)])
Xij_t = rbind(X[1, ], X[2:T, ] - X[1:(T-1), ])

###
### Convert to (Incorrect) Fully-connected Dyadic Data
###

N = T*(T-1)/2
yij = rep(0, N)
Xij = matrix(0, N, p)
dtij = rep(0, N)

idxi = 1 # ith index
idxj = 2 # jth index
for(i in 1:N){
  yij[i] = y[idxj] - y[idxi] # dyadic outcomes
  Xij[i, ] = X[idxj, ] - X[idxi, ] # pairwise difference in X
  dtij[i] = abs(idxj - idxi) # pairwise distance in times
  
  ## update individual indices
  if(idxj == T){
    idxi = idxi + 1
    idxj = idxi + 1
  }else{
    idxj = idxj + 1
  }
}

###
### Basis Functions for Composite Weights
###

## standard basis functions
basis_functions <- list(
  b6 = function(x) x^2,
  b5 = function(x) x^(5/3),
  b4 = function(x) x^(4/3),
  b3 = function(x) x,
  b2 = function(x) x^(2/3),
  b1 = function(x) x^(1/3)
)
nudt <- matrix(nrow = length(dtij), ncol = length(basis_functions))
for(i in seq_along(basis_functions)){
  func = basis_functions[[i]]
  nudt[, i] = func(dtij)
}
nudt.grid <- matrix(nrow = length(1:19), ncol = length(basis_functions))
for(i in seq_along(basis_functions)){
  func = basis_functions[[i]]
  nudt.grid[, i] = func(1:19)
}

###
### Run MCMC Algorithms
###

n.mcmc = 100000
out = mcmc(y, X, n.mcmc) # correctly assuming the AR(1) specification
out_dyad_full = mcmcDyad(yij, Xij, n.mcmc) # using all edges without composite weights
out_dyad_weighted = mcmcWeighted(yij, Xij, nudt, n.mcmc) # using all edges with composite weights

###
### Trace Plots
###

tracePlot = function(chain, par_name, truth, thin = 1, burn = floor(length(chain)*0.1)){
  new.chain <- chain[seq(burn, length(chain), by = thin)]
  title <- "Trace Plot"
  par <- parse(text = par_name)
  dftracePlot <- data.frame(cbind(new.chain = new.chain, iter = 1:length(new.chain)))
  plot <- ggplot(dftracePlot, aes(x = iter, y = new.chain)) + 
    theme_classic() + geom_line() + xlab("Iterations") + 
    ylab(par) + ggtitle(title) + theme(plot.title = element_text(hjust = 0.5))
  if(!missing(truth)){
    plot <- plot + geom_hline(yintercept = truth, color = "red")
  }
  return(plot)
}

tracePlot(out$betaSave[,1], "beta[1]", beta_true[1])
tracePlot(out$betaSave[,2], "beta[2]", beta_true[2])
tracePlot(out$s2Save, "s^2", s2_true)

tracePlot(out_dyad_full$betaSave[,1], "beta[1]", beta_true[1])
tracePlot(out_dyad_full$betaSave[,2], "beta[2]", beta_true[2])

tracePlot(out_dyad_weighted$betaSave[,1], "beta[1]", beta_true[1])
tracePlot(out_dyad_weighted$betaSave[,2], "beta[2]", beta_true[2])
tracePlot(out_dyad_weighted$gammaSave[,1], "gamma[1]")
tracePlot(out_dyad_weighted$gammaSave[,2], "gamma[2]")
tracePlot(out_dyad_weighted$gammaSave[,3], "gamma[3]")
tracePlot(out_dyad_weighted$gammaSave[,4], "gamma[4]")
tracePlot(out_dyad_weighted$gammaSave[,5], "gamma[5]")
tracePlot(out_dyad_weighted$gammaSave[,6], "gamma[6]")

###
### Potential Surface Plotting
###

## computations for plotting
betaMean = apply(out$betaSave, 2, mean)
betaLower = apply(out$betaSave, 2, function(column) quantile(column, probs = 0.025))
betaUpper = apply(out$betaSave, 2, function(column) quantile(column, probs = 0.975))
betaMean_dyad_full = apply(out_dyad_full$betaSave, 2, mean)
betaLower_dyad_full = apply(out_dyad_full$betaSave, 2, function(column) quantile(column, probs = 0.025))
betaUpper_dyad_full = apply(out_dyad_full$betaSave, 2, function(column) quantile(column, probs = 0.975))
betaMean_dyad_weighted = apply(out_dyad_weighted$betaSave, 2, mean)
betaLower_dyad_weighted = apply(out_dyad_weighted$betaSave, 2, function(column) quantile(column, probs = 0.025))
betaUpper_dyad_weighted = apply(out_dyad_weighted$betaSave, 2, function(column) quantile(column, probs = 0.975))

## data preparation (w/ AR(1) specification)
mn = 4
Model = rep(c("Truth", "AR(1)", "Full", "Weighted"), each = T)
Time = rep(1:T, times = mn)
PosteriorMean <- c(X%*%beta_true, X%*%betaMean, X%*%betaMean_dyad_full, X%*%betaMean_dyad_weighted)
LowerCI <- c(rep(NA, T), X%*%betaLower, X%*%betaLower_dyad_full, X%*%betaLower_dyad_weighted)
UpperCI <- c(rep(NA, T), X%*%betaUpper, X%*%betaUpper_dyad_full, X%*%betaUpper_dyad_weighted)
models_data <- data.frame(Model, Time, PosteriorMean, LowerCI, UpperCI)
truth_data <- dplyr::filter(models_data, Model == "Truth")
non_truth_data <- dplyr::filter(models_data, Model != "Truth")

## plotting (w/ AR(1) specification)
ggplot(as.data.frame(models_data), aes(x = Time, y = PosteriorMean, group = Model, color = Model)) +
  geom_line(data = truth_data, linetype = "dashed") +
  geom_line(data = non_truth_data, linetype = 3) + 
  geom_ribbon(aes(ymin = LowerCI, ymax = UpperCI, fill = Model), alpha = 0.2) +
  labs(y = "Posterior Mean", x = "Time", title = "Potential Over Time") +
  theme_minimal() +
  scale_color_manual(values = c("Truth" = "black", "AR(1)" = "green", "Full" = "red", "Weighted" = "blue")) +
  scale_fill_manual(values = c("Truth" = "black", "AR(1)" = "green", "Full" = "red", "Weighted" = "blue"))

## data preparation (w/o AR(1) specification)
mn = 3
Model = rep(c("Truth", "Scenario 1", "Scenario 2"), each = T)
Time = rep(1:T, times = mn)
PosteriorMean <- c(X%*%beta_true, X%*%betaMean_dyad_full, X%*%betaMean_dyad_weighted)
LowerCI <- c(rep(NA, T), X%*%betaLower_dyad_full, X%*%betaLower_dyad_weighted)
UpperCI <- c(rep(NA, T), X%*%betaUpper_dyad_full, X%*%betaUpper_dyad_weighted)
models_data <- as.data.frame(data.frame(Model, Time, PosteriorMean, LowerCI, UpperCI))
truth_data <- dplyr::filter(models_data, Model == "Truth")
non_truth_data <- dplyr::filter(models_data, Model != "Truth")

png("potPlot_ar1.png", width = 800, height = 600)
## plotting (w/o AR(1) specification)
ggplot(data = as.data.frame(models_data), aes(x = Time, y = PosteriorMean, group = Model, color = Model)) +
  geom_line(data = truth_data, linetype = "dashed") +
  geom_line(data = non_truth_data, linetype = 3) + 
  geom_ribbon(aes(ymin = LowerCI, ymax = UpperCI, fill = Model), alpha = 0.2) +
  labs(y = expression(rho), x = "Time", title = "Potential Over Time") +
  theme_minimal() +
  scale_color_manual(values = c("Truth" = "black", "Scenario 1" = "red", "Scenario 2" = "blue"),
                     breaks = c("Truth", "Scenario 1", "Scenario 2")) +
  scale_fill_manual(values = c("Truth" = "black", "Scenario 1" = "red", "Scenario 2" = "blue"),
                    breaks = c("Truth", "Scenario 1", "Scenario 2"))
dev.off()

###
### Image Plot
###

postGamma = apply(out_dyad_weighted$gammaSave, 2, mean)
## set-up image data
library(lattice)
times = 1:T
imageData = expand.grid(X = times, Y = times)

## compute weights given distance in time 
imageDataZ = rep(0, dim(imageData)[1])
for(i in 1:length(imageDataZ)){
  tmpdt = abs(imageData$X[i] - imageData$Y[i])
  imageDataZ[i] = ifelse(tmpdt == 0, 0, 1/exp(nudt.grid[tmpdt, ]%*%postGamma))
}
imageData$Z = imageDataZ

## get nice colors
library(viridisLite)
colfunc <- colorRampPalette(c("#a6cfdd", "#4E598C"))
coul <- colfunc(100)
png("imagePlot.png", width = 600, height = 600)
print(levelplot(Z ~ X*Y, data = imageData, xlab = "Times", ylab = "Times", col.regions = coul))
dev.off()