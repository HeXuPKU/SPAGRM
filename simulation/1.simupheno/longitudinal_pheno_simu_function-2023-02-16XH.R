
library(dplyr)
#### function to generate longitudinal phenotype.
pheno_simu = function(nSub,                           # number of independent individuals.
                      nFam = 0,                       # number of families.
                      FamMode = "4-members",          # Structure of families, default is 4.
                      mi = 6:15,                      # Vector of observations for each individual.
                      betaG = 0,                      # Genetic effect for BS variance, a number.
                      tauG = 0,                       # Genetic effect for WS variance, a number.
                      Geno = NULL,                    # Genotype vector, default equal NULL.
                      cov = rbind(c(2, 0, 0.2),
                                  c(0, 1.2, 0.1),
                                  c(0.2, 0.1, 1)),    # covariance matrix of random effect gamma and omega.
                      bvectoBS = NULL,                # random effect passes to BS variance, match the family structure.
                      bvectoWS = NULL)                # random effect passes to WS variance, match the family structure.

{
  if(nFam == 0){
    N = nSub
    SubjID = UnrelatSubjID = paste0("U-", 1:nSub)
  }else{
    if(FamMode == "4-members"){
      N = nSub + 4 * nFam
      SubjID = c(paste0("F4-", rep(1:nFam, each=4), "-", 1:4), paste0("U-", 1:nSub))
      UnrelatSubjID = c(paste0("F4-", rep(1:nFam, each=2), "-", 1:2), paste0("U-", 1:nSub))
    }else if(FamMode == "10-members"){
      N = nSub + 10 * nFam
      SubjID = c(paste0("F10-", rep(1:nFam, each=10), "-", 1:10), paste0("U-", 1:nSub))
      UnrelatSubjID = c(paste0("F10-", rep(1:nFam, each=4), "-", 1:4), paste0("U-", 1:nSub))
    }else if(FamMode == "20-members"){
      N = nSub + 20 * nFam
      SubjID = c(paste0("F20-", rep(1:nFam, each=20), "-", 1:20), paste0("U-", 1:nSub))
      UnrelatSubjID = c(paste0("F20-", rep(1:nFam, each=8), "-", 1:8), paste0("U-", 1:nSub))
    }else{
      stop("Please check if the FamMode is in '4-members', '10-members' or '20-members'!")
    }
  }
  
  if(is.null(Geno)){
    Geno = rep(0, N)
  }else if(length(Geno) != N){
    stop("Length of geno should correspond to number of individuals.")
  }
  
  if(is.null(bvectoBS)){
    bvectoBS = rep(0, N)
  }else if(length(bvectoBS) != N){
    stop("Length of random vector should correspond to number of individuals.")
  }
  if(is.null(bvectoWS)){
    bvectoWS = rep(0, N)
  }else if(length(bvectoWS) != N){
    stop("Length of random vector should correspond to number of individuals.")
  }
  
  p = 4; # ncol of X.
  l = 4; # ncol of W.
  q = 2; # ncol of Z.
  
  # generate covariates and corresponding coefficients.
  xone = rbinom(N, 1, 0.5) # binary variable (0 or 1).
  xtwo = rnorm(N, 0, 1) # time-invariant standard normal variable.
  # xthree will be generated in the loop.
  beta = c(1, 0.5, 0.5, -0.3) # corresponding coefficients.
  tau = c(0.25, 0.3, -0.15, 0.1) # corresponding coefficients.
  
  # generate random effect gamma and omega
  Lgamma_omega = t(chol(cov))
  lgamma_omega = Lgamma_omega[q+1, 1:q]
  gamma_omega = MASS::mvrnorm(n = N, mu = rep(0, nrow(cov)), Sigma = cov)

  # generate mi_i.
  m = sample(mi, N, replace = TRUE)
  covdata = data.frame()
  
  for(i in 1:N){
    covdata_i = indivi_pheno_simu(mi_i = m[i],
                                  Gi = Geno[i],
                                  xonei = xone[i],
                                  xtwoi = xtwo[i],
                                  bvectoBSi = bvectoBS[i],
                                  bvectoWSi = bvectoWS[i],
                                  beta = beta, tau = tau,
                                  betaGi = betaG, tauGi = tauG,
                                  gamma_omegai = gamma_omega[i,],
                                  lgamma_omega = lgamma_omega)
    
    covdata = dplyr::bind_rows(covdata,
                               covdata_i)
    
    if(i %% 1000 ==0) cat("Have completed", i,"individuals.\n")
  }
  covdata = covdata %>% mutate(SubjID = rep(SubjID, times = m)) %>%
    mutate(UNRELATED = SubjID %in% UnrelatSubjID)
  
  return(covdata)
}

indivi_pheno_simu = function(mi_i, Gi, xonei, xtwoi,
                             bvectoBSi,  bvectoWSi,
                             beta, tau, betaGi, tauGi,
                             gamma_omegai, lgamma_omega)
{
  xthreei = rnorm(mi_i, 0, 1) # a vector
  zonei = rnorm(mi_i, 0, 1) # a vector
  
  meanyi = beta[1] + xonei*beta[2] + xtwoi*beta[3] + xthreei*beta[4] +
    gamma_omegai[1] + zonei*gamma_omegai[2] + betaGi*Gi + bvectoBSi
  
  sdyi = exp(0.5*(tau[1] + xonei*tau[2] + xtwoi*tau[3] + xthreei*tau[4] +
                    c(lgamma_omega %*% gamma_omegai[1:2]) + gamma_omegai[3] + tauGi*Gi + bvectoWSi))
  
  yi = meanyi + sdyi * rnorm(mi_i)
  
  covdata_i = data.table::data.table(xone = xonei, xtwo = xtwoi, xthree = xthreei,
                                     zone = zonei, pheno = yi)
  
  return(covdata_i)
}
