
### Here we use pedigree data consisting of 6,250 4-member families (25,000 individuals) to show how to use Chow-Liu algorithm 
### to approximate the joint distribution of genotypes and their corresponding moment generating function (MGF), and compare it
### with that using a normal distribution. For how to use Chow-Liu algorithm to approximate genotype distribution approximation 
### in actual analysis, please refer to 'chou.liu.tree()' in "https://github.com/HeXuPKU/SPAGRM/blob/main/old_version/subfunc-2023-03-22XH.R".

library(dplyr)
library(tidyr)
library(igraph)
library(ggplot2)

# function to calculate mutual information entropy
Mutual_Information = function(a, b, c,
                              pa, pb, pc)
{
  result = c()
  for(i in 1:length(a)){
    pi = a[i]*pa + b[i]*pb + c[i]*pc
    temp = pi*log(pi/pa)
    result = c(result, sum(temp, na.rm = T))
  }
  return(result)
}

# function to approximate the discrete union probability.
get.arr.prob = function(n,              # the number of family members
                        CLT,            # a data.frame containing ID information and a, b, c
                        p0, pa, pb, pc) # Marginal probability and joint probability
{
  arr.dims = rep(3, n)
  arr.prob = array(1, arr.dims)
  for(i in 1:n)
    dimnames(arr.prob)[[i]] = paste0("ID",i,":",0:2) 
  
  m1 = nrow(CLT); m2 = m1 - 1
  vec = as.numeric(c(CLT$ID1, CLT$ID2))
  vec = vec[duplicated(vec)]
  
  for(k in 1:m1)
  {
    matrix.prob = matrix(CLT$a[k]*pa + CLT$b[k]*pb + CLT$c[k]*pc, 3, 3)
    matrix.index = as.numeric(c(CLT$ID1[k], CLT$ID2[k]))
    for(i in 1:3){
      for(j in 1:3){
        indexString = rep("", n)
        indexString[matrix.index[1]] = i
        indexString[matrix.index[2]] = j
        indexString = paste0(indexString, collapse = ",")
        cmd = paste0("arr.prob[",indexString,"] = arr.prob[", indexString, "] * matrix.prob[",i,",",j,"]")
        # "arr.prob[1,1,] = arr.prob[1,1,] * matrix.prob[1,1]"
        eval(parse(text = cmd))
      }
    }
  }
  
  for(k in 1:m2)
  {
    vector.prob = p0
    vector.index = vec[k]
    for(i in 1:3){
      indexString = rep("", n)
      indexString[vector.index] = i
      indexString = paste0(indexString, collapse = ",")
      cmd = paste0("arr.prob[",indexString,"] = arr.prob[", indexString, "] / vector.prob[",i,"]")
      # arr.prob[,,1] = arr.prob[,,1] / vector.prob[1]"
      eval(parse(text = cmd))
    }
  }
  return(arr.prob)
}


N = 25000
n = 4
nFam = 6250

### calculate GRM and IBD
# since pedigree information is known, we use the theoretical GRM and IBD.

GRM = rbind(c(1, 0, 0.5, 0.5),
            c(0, 1, 0.5, 0.5),
            c(0.5, 0.5, 1, 0.5),
            c(0.5, 0.5, 0.5, 1))

IBD = rbind(data.frame(ID1 = 1, ID2 = 2, a = 1, b = 0, c = 0),
            data.frame(ID1 = 1, ID2 = 3, a = 0, b = 1, c = 0),
            data.frame(ID1 = 1, ID2 = 4, a = 0, b = 1, c = 0),
            data.frame(ID1 = 2, ID2 = 3, a = 0, b = 1, c = 0),
            data.frame(ID1 = 2, ID2 = 4, a = 0, b = 1, c = 0),
            data.frame(ID1 = 3, ID2 = 4, a = 0.25, b = 0.5, c = 0.25))

# later used for construct S.est
arr.index = c()
for(i in 1:n)
{
  indexString = rep("c(1, 1, 1)", n)
  indexString[i] = "0:2"
  indexString = paste0(indexString, collapse = "%o%")
  cmd = paste0("arr.index = c(arr.index, list(arr.index", i, "=", indexString, "))")
  eval(parse(text = cmd))
}


### Looping over different scenarios; calculate log of MGF, which is also called CGF.

results = c()

for(para in c("beta", "tau"))
{
  cat(para, "\n")
  
  if(para == "beta")
  {
    # ResidMat = data.frame(SubjID = 1:25000, Resid = rnorm(25000)) # use normal distribution to simulate residuals for testing beta = 0.
    
    ResidMat = data.table::fread("/gdata01/user/xuhe/SPA-GRM/simulation-2023-09-19/scenario2/residuals/ResidB-1.beta.txt")
  }else
  {
    # ResidMat = data.frame(SubjID = 1:25000, Resid = invgamma::rinvgamma(n = 25000, shape = 5, rate = 10)) # use inverse-gamma distribution to simulate residuals for testing tau = 0.
    
    ResidMat = data.table::fread("/gdata01/user/xuhe/SPA-GRM/simulation-2023-09-19/scenario2/residuals/ResidB-1.tau.txt")
  }
  
  ResidMat$Resid = scale(ResidMat$Resid)
  
  for(mu in c(0.3, 0.01, 0.001))
  {
    cat(mu, "\n")
    
    # Marginal probability distribution
    p0 = c((1-mu)^2, 2*mu*(1-mu), mu^2)
    
    # p = c(G00, G10, G20, G01, G11, G21, G02, G12, G22)
    # share 0 allele.
    pa = c((1-mu)^4, 2*mu*(1-mu)^3, mu^2*(1-mu)^2, 2*mu*(1-mu)^3, 4*mu^2*(1-mu)^2, 
           2*mu^3*(1-mu), mu^2*(1-mu)^2, 2*mu^3*(1-mu), mu^4)
    # share 1 allele.
    pb = c((1-mu)^3, mu*(1-mu)^2, 0, mu*(1-mu)^2, mu*(1-mu), mu^2*(1-mu), 0, mu^2*(1-mu), mu^3)
    # share 2 alleles.
    pc = c((1-mu)^2, 0, 0, 0, 2*mu*(1-mu), 0, 0, 0, mu^2)
    
    ### construct CLT
    Weight = Mutual_Information(IBD$a, IBD$b, IBD$c,
                                pa, pb, pc)
    
    graph = graph_from_data_frame(IBD, directed = T)
    chow.liu.tree = mst(graph, weights = - Weight, algorithm = "prim")
    chow.liu.tree = chow.liu.tree %>% get.edgelist() %>% as.data.frame()
    colnames(chow.liu.tree) = c("ID1", "ID2");
    chow.liu.tree$ID1 = as.numeric(chow.liu.tree$ID1); 
    chow.liu.tree$ID2 = as.numeric(chow.liu.tree$ID2); 
    
    CLT = merge(chow.liu.tree, IBD, all.x = T)
    
    arr.prob = get.arr.prob(n, CLT, p0, pa, pb, pc); print(sum(arr.prob))
    
    
    # user later for normal approximation
    Mean = 2 * mu * sum(ResidMat$Resid)
    Var = 2*mu*(1 - mu) * as.numeric(t(ResidMat$Resid) %*% Matrix::bdiag(rep(list(GRM), nFam)) %*% ResidMat$Resid)
    
    interval1 = seq(0.1, 3, 0.05); interval2 = seq(-0.1, 0.1, 0.05); interval3 = -1 * interval1
    interval = sort(c(interval1, interval2, interval3))
    
    for(t in interval)
    {
      Resid = split(ResidMat$Resid, cut(seq_along(ResidMat$Resid), nFam, labels = FALSE))
      
      MGF.norm = Mean*t + 0.5*Var*t^2
      
      results = rbind(results,
                      c(MGF.norm, t, Methods = "MGF.norm", Para = para, MAF = mu))
      
      MGF.est = c(); MGF.true = c()
      
      for(i in 1:nFam)
      {
        R = Resid[[i]]
        
        S.est = array(rowSums(mapply(function(x, y) x*y, arr.index, R)), dim = rep(3, n))
        
        temp_MGF.est = sum(exp(t*S.est)*arr.prob)
        
        temp_MGF.true = 0.25 * sum(prod(1 - mu + mu * exp(t * c(sum(R[c(1,3,4)]), sum(R[c(2,3,4)]), R[1], R[2]))),
                                 prod(1 - mu + mu * exp(t * c(sum(R[c(1,3,4)]), sum(R[c(2,3)]), sum(R[c(2,4)]), R[1]))),
                                 prod(1 - mu + mu * exp(t * c(sum(R[c(2,3,4)]), sum(R[c(1,3)]), sum(R[c(1,4)]), R[2]))),
                                 prod(1 - mu + mu * exp(t * c(sum(R[c(1,3)]), sum(R[c(1,4)]), sum(R[c(2,3)]), sum(R[c(2,4)])))))
        
        MGF.est = c(MGF.est, log(temp_MGF.est))
        MGF.true = c(MGF.true, log(temp_MGF.true))
      }
      
      MGF.est = sum(MGF.est); MGF.true = sum(MGF.true)
      
      results = rbind(results,
                      c(MGF.est, t, Methods = "MGF.est", Para = para, MAF = mu),
                      c(MGF.true, t, Methods = "MGF.true", Para = para, MAF = mu))
    }
  }
}


colnames(results) = c("logMGF", "t", "Methods", "Para", "MAF")
results = as.data.frame(results)

results$logMGF = as.numeric(results$logMGF)
results$t = as.numeric(results$t)
results$MAF = as.numeric(results$MAF)

results = results %>% drop_na %>% filter(logMGF != Inf)

results$Para = factor(results$Para, levels = c("beta", "tau"),
                      labels = c(expression(paste("Residuals for testing ", italic(beta[g]) == 0)), 
                                 expression(paste("Residuals for testing ", italic(tau[g]) == 0))))
results$MAF = factor(results$MAF, levels = c(0.3, 0.01, 0.001),
                     labels = c(expression(paste("MAF = 0.3")), expression(paste("MAF = 0.01")), expression(paste("MAF = 0.001"))))

results$Methods = factor(results$Methods, levels = c("MGF.norm", "MGF.est", "MGF.true"))

p = ggplot(results, aes(x = t, y = logMGF, color = Methods, shape = Methods)) +
  geom_point() + geom_line() + geom_abline(slope = 0, intercept = 0) + 
  ggh4x::facet_grid2(Para ~ MAF, scales = "free", labeller = labeller(MAF = label_parsed, Para = label_parsed), axes = "all") + 
  theme_bw() + scale_x_continuous(breaks = c(-3:3)) +
  xlab(expression(paste("t"))) + ylab(expression(paste("Natural logarithm of MGF"))) +
  scale_color_discrete(breaks = c("MGF.norm", "MGF.est", "MGF.true"),
                       labels = c(expression(paste("Normal distribution")), expression(paste("Saddlepoint approximation based on Chow-Liu algorithm")), expression(paste("Theoretical distribution")))) + 
  scale_shape_discrete(breaks = c("MGF.norm", "MGF.est", "MGF.true"),
                       labels = c(expression(paste("Normal distribution")), expression(paste("Saddlepoint approximation based on Chow-Liu algorithm")), expression(paste("Theoretical distribution")))) + 
  theme(legend.position = "bottom", legend.title = element_blank(), 
        plot.title = element_text(hjust = 0, size = 15, vjust = -0.05),
        axis.title = element_text(size = 10), axis.text = element_text(size = 11),
        strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12))

ggsave("/gdata01/user/xuhe/SPA-GRM/simulation-2024-08-13/evaluate_CLT.png", 
       p, width = 10.5, height = 8.5)

