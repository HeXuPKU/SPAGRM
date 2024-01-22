
# cd "/gdata01/user/xuhe/family_relatedness/simulation-2023-02-15/"
# sbatch -J pheno_simu --mem=4G -t 1-0:0 --array=1-10000 -o log/%A_%a.log --wrap='Rscript /gdata01/user/xuhe/family_relatedness/simulation-2023-02-15/code/type1error_pheno_simu-2023-02-16.R $SLURM_ARRAY_TASK_ID'
args = commandArgs(TRUE)
print(args)
print(sessionInfo())
n.cpu = as.numeric(args)
print(n.cpu)

source("/gdata01/user/xuhe/family_relatedness/simulation-2023-02-15/code/longitudinal_pheno_simu_function-2023-02-16XH.R")

# A : 4-member families.
nSub = 25e3
nFam = 6250
randMat = data.table::fread("/gdata02/master_data1/Related_Subjects/RandomEffect/randMat4Members.txt")
bvectoBS = c(as.numeric(t(randMat[sample(2.5e6, nFam),])), rnorm(nSub))
bvectoWS = c(as.numeric(t(randMat[sample(2.5e6, nFam),])), rnorm(nSub))

covfileA = pheno_simu(nSub = nSub, nFam = nFam, FamMode = "4-members",
                      mi = 6:15, betaG = 0, tauG = 0, Geno = NULL,
                      cov = matrix(c(2, 0, 0.2, 0 ,1.2, 0.1, 0.2, 0.1, 1), 3, 3),
                      bvectoBS = bvectoBS, bvectoWS = bvectoWS)

write.csv(covfileA, file = paste0("/gdata01/user/xuhe/SPA-GRM/simulation-2023-02-15/scenario1/pheno/phenoA-", n.cpu, ".csv"))

# B : 10-member families.
nSub = 25e3
nFam = 25e2 # 10-members family.
randMat = data.table::fread("/gdata02/master_data1/Related_Subjects/RandomEffect/randMat10Members.txt")
bvectoBS = c(as.numeric(t(randMat[sample(1e6, nFam),])), rnorm(nSub))
bvectoWS = c(as.numeric(t(randMat[sample(1e6, nFam),])), rnorm(nSub))

covfileB = pheno_simu(nSub = nSub, nFam = nFam, FamMode = "10-members",
                      mi = 6:15, betaG = 0, tauG = 0, Geno = NULL,
                      cov = matrix(c(2, 0, 0.2, 0 ,1.2, 0.1, 0.2, 0.1, 1), 3, 3),
                      bvectoBS = bvectoBS, bvectoWS = bvectoWS)

write.csv(covfileB, file = paste0("/gdata01/user/xuhe/SPA-GRM/simulation-2023-02-15/scenario1/pheno/phenoB-", n.cpu, ".csv"))

# C : unrelated.
nSub = 50e3
nFam = 0
bvectoBS = rnorm(nSub)
bvectoWS = rnorm(nSub)

covfileC = pheno_simu(nSub = nSub, nFam = nFam, FamMode = "4-members",
                      mi = 6:15, betaG = 0, tauG = 0, Geno = NULL,
                      cov = matrix(c(2, 0, 0.2, 0 ,1.2, 0.1, 0.2, 0.1, 1), 3, 3),
                      bvectoBS = bvectoBS, bvectoWS = bvectoWS)

write.csv(covfileC, file = paste0("/gdata01/user/xuhe/SPA-GRM/simulation-2023-02-15/scenario1/pheno/phenoC-", n.cpu, ".csv"))
