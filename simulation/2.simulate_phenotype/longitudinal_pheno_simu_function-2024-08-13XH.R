
library(dplyr)
#### function to generate longitudinal phenotype.
pheno_simu = function(nSub,                           # number of independent individuals.
                      nFam = 0,                       # number of families.
                      FamMode = "4-members",          # Structure of families, e.g., 4, 10, 20, default is 4.
                      mi = 6:15,                      # Vector of observations for each individual.
                      betaG = 0,                      # Genetic effect for mean level trajectories, a number.
                      Geno = NULL,                    # Genotype vector, default equal NULL.
                      bvectoBS = NULL,                # random effect passes to BS variance, match the family structure.
                      corstr = "exc",                 # correlation structures for repeated measurements, e.g., exc, ar1.
                      rho)                            # parameter used to generate correlation structure.
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
  
  if(corstr != "exc" & corstr != "ar1")
  {
    stop("corstr should be 'exc' or 'ar1'.")
  }
  
  p = 4; # ncol of X.
  
  # generate covariates and corresponding coefficients.
  xone = rbinom(N, 1, 0.5) # binary variable (0 or 1).
  xtwo = rnorm(N, 0, 1) # time-invariant standard normal variable.
  # xthree will be generated in the loop.
  beta = c(1, 0.5, 0.5, -0.3) # corresponding coefficients.
  
  # generate mi_i.
  m = sample(mi, N, replace = TRUE)
  covdata = data.frame()
  
  for(i in 1:N){
    covdata_i = indivi_pheno_simu(mi_i = m[i],
                                  Gi = Geno[i],
                                  xonei = xone[i],
                                  xtwoi = xtwo[i],
                                  bvectoBSi = bvectoBS[i],
                                  beta = beta, betaGi = betaG,
                                  corstr = corstr, rho = rho)
    
    covdata = dplyr::bind_rows(covdata,
                               covdata_i)
    
    if(i %% 1000 ==0) cat("Have completed", i,"individuals.\n")
  }
  covdata = covdata %>% mutate(SubjID = rep(SubjID, times = m)) %>%
    mutate(UNRELATED = SubjID %in% UnrelatSubjID)
  
  return(covdata)
}

indivi_pheno_simu = function(mi_i, Gi, xonei, xtwoi,
                             bvectoBSi, beta, betaGi, corstr, rho)
{
  if(corstr == "exc")
  {
    working.correlation = matrix(rho, nrow = mi_i, ncol = mi_i); 
    diag(working.correlation) = 1
  }else
  {
    working.correlation = toeplitz(1*rho^(0:(mi_i-1)))
  }
  xthreei = rnorm(mi_i, 0, 1) # a vector
  
  meanyi = beta[1] + xonei*beta[2] + xtwoi*beta[3] + xthreei*beta[4] + betaGi*Gi + bvectoBSi
  
  sdyi = MASS::mvrnorm(n = 1, mu = rep(0, mi_i), Sigma = 2 * working.correlation)
  
  yi = meanyi + sdyi
  
  covdata_i = data.table::data.table(xone = xonei, xtwo = xtwoi, xthree = xthreei, pheno = yi)
  
  return(covdata_i)
}
