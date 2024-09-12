
###### part2 analyzing ######
# cd "/gdata01/user/xuhe/family_relatedness/simulation-2023-09-19/"
# sbatch -J step2 --mem=2G -t 1-0:0 --array=1-1000 -o log/%A_%a.log --wrap='Rscript /gdata01/user/xuhe/family_relatedness/simulation-2023-09-19/code/SPAGRM_step2-2023-11-27.R $SLURM_ARRAY_TASK_ID'
args=commandArgs(TRUE)
print(args)
print(sessionInfo())
n.cpu=as.numeric(args)
print(n.cpu)

source("/gdata01/user/xuhe/family_relatedness/simulation-2023-03-22/code/subfunc-2023-03-22XH.R")

table = data.frame(nrep = 1:200,
                   scr = rep(3, each = 200),
                   type = rep(c("A", "B", "C", "D", "E"), each = 200),
                   dir = "working.data_new")

nrep = table$nrep[n.cpu]; scr = table$scr[n.cpu]
type = table$type[n.cpu]; dir = table$dir[n.cpu]

cat("ncpu is", n.cpu, "; scenario is", scr, "; type is", type, "; rep is", nrep ,"; dir is", dir, ".\n")

GenoMat = paste0("/gdata01/user/xuhe/SPA-GRM/simulation-2023-09-19/GenoMat", type, ".RData")
load(GenoMat)

SPAGRM.beta = paste0("/gdata01/user/xuhe/SPA-GRM/simulation-2023-09-19/scenario", scr, "/", dir, "/SPAGRM.beta", type, "-", nrep, ".RData")
SPAGRM.tau = paste0("/gdata01/user/xuhe/SPA-GRM/simulation-2023-09-19/scenario", scr, "/", dir, "/SPAGRM.tau", type, "-", nrep, ".RData")

if(file.exists(SPAGRM.beta))
{
  load(SPAGRM.beta)
  
  output.beta = SPA_G(SPAGRM.beta,
                      GenoMat,
                      Cutoff = 2,
                      impute.method = "fixed",
                      missing.cutoff = 0.15,
                      hwepval.cutoff = 0,
                      min.maf = 0.0001,
                      G.model = "Add",
                      zeta = 0,
                      tol = 1e-4)
  
  data.table::fwrite(output.beta,
                     file = paste0("/gdata01/user/xuhe/SPA-GRM/simulation-2023-09-19/scenario", scr, "/SPA.GRMpval_new/pvalue.beta", type, "-", nrep, ".txt"),
                     row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
}

if(file.exists(SPAGRM.tau))
{
  load(SPAGRM.tau)
  
  output.tau = SPA_G(SPAGRM.tau,
                     GenoMat,
                     Cutoff = 2,
                     impute.method = "fixed",
                     missing.cutoff = 0.15,
                     hwepval.cutoff = 0,
                     min.maf = 0.0001,
                     G.model = "Add",
                     zeta = 0,
                     tol = 1e-5)
  
  data.table::fwrite(output.tau,
                     file = paste0("/gdata01/user/xuhe/SPA-GRM/simulation-2023-09-19/scenario", scr, "/SPA.GRMpval_new/pvalue.tau", type, "-", nrep, ".txt"),
                     row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
}