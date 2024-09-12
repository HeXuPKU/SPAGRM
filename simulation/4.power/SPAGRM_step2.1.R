
# ###### part1 preconditioning ######
# cd "/gdata01/user/xuhe/family_relatedness/simulation-2023-09-19/"
# sbatch -J step1 --mem=3G -t 1-0:0 --array=1-1000 -o log/%A_%a.log --wrap='Rscript /gdata01/user/xuhe/family_relatedness/simulation-2023-09-19/code/SPAGRM_step1-2023-11-27.R $SLURM_ARRAY_TASK_ID'
args=commandArgs(TRUE)
print(args)
print(sessionInfo())
n.cpu=as.numeric(args)

table = data.frame(nrep = 1:200,
                   scr = rep(3, each = 200),
                   type = rep(c("A", "B", "C", "D", "E"), each = 200),
                   dir = "residuals_new")

nrep = table$nrep[n.cpu]; scr = table$scr[n.cpu]
type = table$type[n.cpu]; dir = table$dir[n.cpu]

cat("ncpu is", n.cpu, "; scenario is", scr, "; type is", type, "; rep is", nrep ,"; dir is", dir, ".\n")

library(GRAB)

ResidMatFile.beta = paste0("/gdata01/user/xuhe/SPA-GRM/simulation-2023-09-19/scenario", scr, "/", dir, "/Resid", type, "-", nrep, ".beta.txt")

ResidMatFile.tau = paste0("/gdata01/user/xuhe/SPA-GRM/simulation-2023-09-19/scenario", scr, "/", dir, "/Resid", type, "-", nrep, ".tau.txt")

SparseGRMFile = paste0("/gdata01/user/xuhe/SPA-GRM/simulation-2023-09-19/SparseGRM_", type, ".txt")

PairwiseIBDFile = paste0("/gdata01/user/xuhe/SPA-GRM/simulation-2023-09-19/PairwiseIBD_", type, ".txt")

if(file.exists(ResidMatFile.beta))
{
  SPAGRM.beta = SPAGRM.NullModel(ResidMatFile.beta, SparseGRMFile, PairwiseIBDFile)
  
  save(SPAGRM.beta,
       file = paste0("/gdata01/user/xuhe/SPA-GRM/simulation-2023-09-19/scenario", scr, "/working.data_new/SPAGRM.beta", type, "-", nrep, ".RData"))
}

if(file.exists(ResidMatFile.tau))
{
  SPAGRM.tau = SPAGRM.NullModel(ResidMatFile.tau, SparseGRMFile, PairwiseIBDFile)
  
  save(SPAGRM.tau,
       file = paste0("/gdata01/user/xuhe/SPA-GRM/simulation-2023-09-19/scenario", scr, "/working.data_new/SPAGRM.tau", type, "-", nrep, ".RData"))
}