############################################################
####
#### MCMC Function to Implement Dyadic Regression Model
#### (Spatial and Temporal Lags in Composite Weights)
####
############################################################

function mcmc_both(y, X, H, K, phi_list, R, nMCMC, dt, ds)

    ## y is the vector of (dyadic) genetic dissimilarities
    ## X is the design matrix of pairwise differences
    ## H is the mapping matrix for theta
    ## K is the mapping matrix for location-level spatial random effects
    ## phi_list is a list of possible phi's to sample from
    ## R is a list of covariance functions for various specifications of phi for eta's covariance matrix
    ## nMCMC is the number of iterations in the algorithm
    ## dt is a vector of differences in time for dyadic responses
    ## ds is a vector of pairwise Euclidean distances
    
    ###
    ### Set values
    ###

    N = length(y)
    n = size(H)[2]
    p = size(X)[2]
    s = size(K)[2]

    ###
    ### Hyperparameters
    ###

    mu_beta = zeros(p)
    Sig_beta = 10^6*I

    mu_eta = zeros(s)

    alpha_y = 0.01
    beta_y = 0.01

    alpha_eta = 0.01
    beta_eta = 0.01

    alpha_s2_theta = 0.01
    beta_s2_theta = 0.01
    
    ## for gamma M-H update
    alpha_gamma = 2
    beta_gamma = 2
    ls_gamma = 0 # from Roberts & Rosenthal (2009)
    M = 5 # global maximal parameter value for ls

    ###
    ### Initial values
    ###

    beta = rand(p)
    theta = rand(n)
    eta = rand(s)
    s2_y = 1
    s2_eta = 1
    s2_theta = 100
    gamma = rand(2)
    phi_idx = Int(round(length(phi_list)/2))
    phi = phi_list[phi_idx]
    Rinv = inv(R[:, :, phi_idx])

    weights = exp.(-gamma[1].*dt - gamma[2].*ds)
    Sig_y = Diagonal(ones(N)./weights)
    Sig_y_i = inv(Sig_y)
    
    accGamma = 0 # the number of accepted gammas in the last 50 iterations
    nBatch = 0 # the number of 50-iteration batches that we've processed

    ###
    ### Save matrices
    ###

    betaSave = zeros(p, nMCMC)
    betaSave[:, 1] = copy(beta)

    thetaSave = zeros(n, nMCMC) 
    thetaSave[:, 1] = copy(theta)

    etaSave = zeros(s, nMCMC)
    etaSave[:, 1] = copy(eta)

    s2_ySave = zeros(nMCMC)
    s2_ySave[1] = copy(s2_y)

    s2_etaSave = zeros(nMCMC)
    s2_etaSave[1] = copy(s2_eta)

    s2_thetaSave = zeros(nMCMC)
    s2_thetaSave[1] = copy(s2_theta)

    gammaSave = zeros(2, nMCMC)
    gammaSave[:, 1] = copy(gamma)

    phiSave = zeros(nMCMC)
    phiSave[1] = copy(phi)

    phiIdxSave = zeros(nMCMC)
    phiIdxSave[1] = copy(phi_idx) 

    ###
    ### Repeat Calculations
    ###

    ## RSR terms
    H0 = [H K]
    H0X = H0'*X

    ## eta constraints
    constraint_eta = inv(K'*K)*K'*X
    constraint_theta = inv(H'*H)*H'*ones(N)

    ## other repeated computations
    H = sparse(H)
    K = sparse(K)
    Sig_beta_inv = inv(Sig_beta)
    alphaStar_y = alpha_y + N/2
    alphaStar_eta = alpha_eta + s/2
    alphaStar_s2_theta = alpha_s2_theta + n/2

    ## RSR terms
    XtXiXtH = inv(X'*X)*X'*H0
    HtXXtXi = H0'*X*inv(X'*X)

    ###
    ### Gibbs Sampler
    ###

    for k in ProgressBar(2:nMCMC)

        ###
        ### Update β
        ###

        beta_Ainv = inv(X'*Sig_y_i*X/s2_y + Sig_beta_inv)
        b = ((y' - eta'*K' - theta'*H')*Sig_y_i*X/s2_y)'
        beta = vec(rand(MvNormal(vec(beta_Ainv*b), Hermitian(beta_Ainv)), 1))

        ###
        ### Update θ
        ###

        theta_Ainv = inv(Matrix(H'*Sig_y_i*H/s2_y + I/s2_theta))
        b = ((y' - beta'*X' - eta'*K')*Sig_y_i*H/s2_y)'
        theta = vec(rand(MvNormal(vec(theta_Ainv*b), Hermitian(theta_Ainv)), 1))

        ###
        ### Update η
        ###

        eta_Ainv = inv(K'*Sig_y_i*K/s2_y + Rinv/s2_eta)
        b = ((y' - beta'*X' - theta'*H')*Sig_y_i*K/s2_y)'
        eta = vec(rand(MvNormal(vec(eta_Ainv*b), Hermitian(eta_Ainv)), 1))
        eta = eta - eta_Ainv*constraint_eta*inv(constraint_eta'*eta_Ainv*constraint_eta)*constraint_eta'*eta

        ###
        ### Update σ²_y
        ###

        alphaStar_y = alpha_y + N/2
        r = y - X*beta - K*eta - H*theta
        betaStar_y = beta_y + 1/2*r'*Sig_y_i*r
        s2_y = rand(InverseGamma(alphaStar_y, betaStar_y), 1)[1]

        ###
        ### Update σ²_η
        ###

        betaStar_eta = beta_eta + 1/2*eta'*Rinv*eta
        s2_eta = rand(InverseGamma(alphaStar_eta, betaStar_eta), 1)[1]

        ###
        ### Update σ²_θ
        ###

        betaStar_s2_theta = beta_s2_theta + 1/2*theta'*theta
        s2_theta = rand(InverseGamma(alphaStar_s2_theta, betaStar_s2_theta), 1)[1]
    
        ###
        ### Update γ
        ###

        ## proposed value via lower truncated random-walk
        gammaStar = zeros(2)
        for indgamma in 1:2
            gammaStar[indgamma] = rand(Normal(gamma[indgamma], 0.001), 1)[1]
            while gammaStar[indgamma] < 0
                gammaStar[indgamma] = rand(Normal(gamma[indgamma], 0.001), 1)[1]
            end
        end

        ## compute starred weights
        weightsStar = exp.(-gammaStar[1].*dt - gammaStar[2].*ds)
        Sig_y_star = Diagonal(ones(N)./weightsStar)

        ## M-H ratio
        mh1 = logpdf(MvNormal(vec(X*beta + K*eta + H*theta), Hermitian(s2_y*Sig_y_star)), y) + logpdf(Gamma(alpha_gamma, beta_gamma), gammaStar[1]) + logpdf(Gamma(alpha_gamma, beta_gamma), gammaStar[2])
        mh2 = logpdf(MvNormal(vec(X*beta + K*eta + H*theta), Hermitian(s2_y*Sig_y)), y) + logpdf(Gamma(alpha_gamma, beta_gamma), gamma[1]) + logpdf(Gamma(alpha_gamma, beta_gamma), gamma[2])
        mh = exp(mh1 - mh2)

        if mh > rand(1)[1]
            accGamma += 1 # we accept a gamma
            gamma = copy(gammaStar)
            weights = copy(weightsStar)
            Sig_y = copy(Sig_y_star)
            Sig_y_i = inv(Sig_y)
        end

        ###
        ### Update ϕ
        ###

        alpha_phi = 400
        beta_phi = 250

        phiStar_idx = sample((phi_idx - 150):(phi_idx + 150), 1)[1]
        while phiStar_idx < 1 || phiStar_idx > length(phi_list) # make sure it's in the proper support
            phiStar_idx = sample((phi_idx - 150):(phi_idx + 150), 1)[1]
        end

        mh1 = logpdf(MvNormal(zeros(s), Hermitian(s2_eta*R[:, :, phiStar_idx])), eta) + logpdf(Gamma(alpha_phi, beta_phi), phi_list[phiStar_idx])
        mh2 = logpdf(MvNormal(zeros(s), Hermitian(s2_eta*R[:, :, phi_idx])), eta) + logpdf(Gamma(alpha_phi, beta_phi), phi_list[phi_idx])
        mh = exp(mh1 - mh2)

        if mh > rand(1)[1]
            phi_idx = copy(phiStar_idx)
            phi = copy(phi_list[phiStar_idx])
            Rinv = inv(R[:, :, phiStar_idx])
        end

        ###
        ### Save values
        ###

        betaSave[:, k] = copy(beta)
        thetaSave[:, k] = copy(theta)
        etaSave[:, k] = copy(eta)
        s2_ySave[k] = copy(s2_y)
        s2_etaSave[k] = copy(s2_eta)
        s2_thetaSave[k] = copy(s2_theta)
        gammaSave[:, k] = copy(gamma)
        phiSave[k] = copy(phi)
        phiIdxSave[k] = copy(phi_idx)

        if(k % 100 == 0) # every 100 iterations, save the information
            @save "both/MCMCoutput.jld2" k betaSave thetaSave etaSave s2_ySave s2_etaSave s2_thetaSave gammaSave phiSave phiIdxSave
        end

    end

    ###
    ### Write output
    ###

    return [betaSave, thetaSave, etaSave, s2_ySave, s2_etaSave, s2_thetaSave, gammaSave, phiSave, phiIdxSave]

end