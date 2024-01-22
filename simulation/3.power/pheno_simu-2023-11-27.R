
# cd "/gdata01/user/xuhe/family_relatedness/simulation-2023-09-19/"
# sbatch -J phenoimu --mem=5G -t 1-0:0 --array=1-200 -o log/%A_%a.log --wrap='Rscript /gdata01/user/xuhe/family_relatedness/simulation-2023-09-19/code/power_pheno_simu-2023-11-27.R $SLURM_ARRAY_TASK_ID'
args = commandArgs(TRUE)
print(args)
print(sessionInfo())
n.cpu = as.numeric(args)
print(n.cpu)

source("/gdata01/user/xuhe/family_relatedness/simulation-2023-02-15/code/longitudinal_pheno_simu_function-2023-02-16XH.R")

# D : 4-member families.
nSub = 25e3
nFam = 6250
randMat = data.table::fread("/gdata02/master_data1/Related_Subjects/RandomEffect/randMat4Members.txt")
bvectoBS = c(as.numeric(t(randMat[sample(2.5e6, nFam),])), rnorm(nSub))
bvectoWS = c(as.numeric(t(randMat[sample(2.5e6, nFam),])), rnorm(nSub))

load("/gdata01/user/xuhe/SPA-GRM/simulation-2023-09-19/GenoMatD.RData")
MAF = colMeans(GenoMat)/2
GenoMatSD=t((t(GenoMat) - 2*MAF)/sqrt(2*MAF*(1-MAF)))
Geno = GenoMatSD %*% (c(rep(-0.1, 5), rep(-0.02, 5)) * log10(MAF)) %>% as.numeric()

covfileA3 = pheno_simu(nSub = nSub, nFam = nFam, FamMode = "4-members",
                       mi = 6:15, betaG = 0, tauG = 0.5, Geno = Geno,
                       cov = matrix(c(2, 0, 0.2, 0 ,1.2, 0.1, 0.2, 0.1, 1), 3, 3),
                       bvectoBS = bvectoBS, bvectoWS = bvectoWS)
write.csv(covfileA3, file = paste0("/gdata01/user/xuhe/SPA-GRM/simulation-2023-09-19/scenario3/pheno_new/phenoD3-", n.cpu, ".csv"))

# E : 10-member families.
nSub = 25e3
nFam = 2500
randMat = data.table::fread("/gdata02/master_data1/Related_Subjects/RandomEffect/randMat10Members.txt")
bvectoBS = c(as.numeric(t(randMat[sample(1e6, nFam),])), rnorm(nSub))
bvectoWS = c(as.numeric(t(randMat[sample(1e6, nFam),])), rnorm(nSub))

load("/gdata01/user/xuhe/SPA-GRM/simulation-2023-09-19/GenoMatE.RData")
MAF = colMeans(GenoMat)/2
GenoMatSD=t((t(GenoMat) - 2*MAF)/sqrt(2*MAF*(1-MAF)))
Geno = GenoMatSD %*% (c(rep(-0.1, 5), rep(-0.02, 5)) * log10(MAF)) %>% as.numeric()

covfileB3 = pheno_simu(nSub = nSub, nFam = nFam, FamMode = "10-members",
                       mi = 6:15, betaG = 0, tauG = 0.5, Geno = Geno,
                       cov = matrix(c(2, 0, 0.2, 0 ,1.2, 0.1, 0.2, 0.1, 1), 3, 3),
                       bvectoBS = bvectoBS, bvectoWS = bvectoWS)
write.csv(covfileB3, file = paste0("/gdata01/user/xuhe/SPA-GRM/simulation-2023-09-19/scenario3/pheno_new/phenoE3-", n.cpu, ".csv"))

# F : unrelated.
nSub = 50e3
nFam = 0
bvectoBS = rnorm(nSub)
bvectoWS = rnorm(nSub)

load("/gdata01/user/xuhe/SPA-GRM/simulation-2023-09-19/GenoMatF.RData")
MAF = colMeans(GenoMat)/2
GenoMatSD=t((t(GenoMat) - 2*MAF)/sqrt(2*MAF*(1-MAF)))
Geno = GenoMatSD %*% (c(rep(-0.1, 5), rep(-0.02, 5)) * log10(MAF)) %>% as.numeric()

covfileC3 = pheno_simu(nSub = nSub, nFam = nFam, FamMode = "4-members",
                       mi = 6:15, betaG = 0, tauG = 0.5, Geno = Geno,
                       cov = matrix(c(2, 0, 0.2, 0 ,1.2, 0.1, 0.2, 0.1, 1), 3, 3),
                       bvectoBS = bvectoBS, bvectoWS = bvectoWS)
write.csv(covfileC3, file = paste0("/gdata01/user/xuhe/SPA-GRM/simulation-2023-09-19/scenario3/pheno_new/phenoF3-", n.cpu, ".csv"))
