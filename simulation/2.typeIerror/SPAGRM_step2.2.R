
###### part2 analyzing ######
# cd "/gdata01/user/xuhe/family_relatedness/simulation-2023-03-22/"
# sbatch -J step2 --mem=12G -t 1-0:0 --array=1-60000 -o log/%A_%a.log --wrap='Rscript /gdata01/user/xuhe/family_relatedness/simulation-2023-03-22/code/type1error_EmpSPAGRM_step2-2023-07-01.R $SLURM_ARRAY_TASK_ID'
args=commandArgs(TRUE)
print(args)
print(sessionInfo())
n.cpu=as.numeric(args)

table = data.frame(nrep = rep(1:10000, times = 3),
                   type = rep(c("A", "B", "C"), each = 10000),
                   genofile = rep(c("nSub_25000_nFam_6250", "nSub_25000_nFam_2500", "nSub_50000"), each = 9000),
                   MAF = rep(c("common", "rare"), each = 30000))
nrep = table[n.cpu, "nrep"]; type = table[n.cpu, "type"];
genofile = table[n.cpu, "genofile"]; MAF = table[n.cpu, "MAF"]
cat("ncpu is ", n.cpu, "; type is ", type, "; nrep is ", nrep, "; genofile is ", genofile, "; MAF is ", MAF,".\n")

library(GRAB)
source("/gdata01/user/xuhe/family_relatedness/simulation-2023-03-22/code/subfunc-2023-03-22XH.R")

obj.precond.beta = paste0("/gdata01/user/xuhe/SPA-GRM/simulation-2023-03-22/scenario1/working.data/obj.precond.beta", type, "-", nrep, ".RData")
obj.precond.tau = paste0("/gdata01/user/xuhe/SPA-GRM/simulation-2023-03-22/scenario1/working.data/obj.precond.tau", type, "-", nrep, ".RData")

if(!file.exists(obj.precond.beta) & !file.exists(obj.precond.tau))
  stop("No obj.precond.beta and obj.precond.tau data.")

if(file.exists(obj.precond.beta))
  load(obj.precond.beta);
if(file.exists(obj.precond.tau))
  load(obj.precond.tau);

pvalue.beta = c(); pvalue.tau = c()

bedFiledir = paste0("/gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/", genofile, MAF, "/")
bedFile = paste0(bedFiledir, genofile, ".bed")

for(i in 1:10)
{
  IDsToIncludeFile = paste0(bedFiledir, "SNPs_group", i, ".txt")

  GenoList = GRAB.ReadGeno(GenoFile = bedFile,
                           control = list(IDsToIncludeFile = IDsToIncludeFile,
                                          ImputeMethod = "mean"))
  
  if(class(obj.precond.beta) == "SPAGRM_NULL_Model")
  {
    output.beta = SPA_G(obj.precond.beta,
                        GenoList$GenoMat,
                        Cutoff = 2,
                        impute.method = "fixed",
                        hwepval.cutoff = 0,
                        missing.cutoff = 0.15,
                        min.maf = 0.0001,
                        G.model = "Add",
                        zeta = - 0.02,
                        tol = 1e-4)
    pvalue.beta = rbind(pvalue.beta,
                        output.beta)
  }

  if(class(obj.precond.tau) == "SPAGRM_NULL_Model")
  {
    output.tau = SPA_G(obj.precond.tau,
                       GenoList$GenoMat,
                       Cutoff = 2,
                       impute.method = "fixed",
                       missing.cutoff = 0.15,
                       min.maf = 0.0001,
                       G.model = "Add",
                       zeta = - 0.0002,
                       tol = 1e-5)
    pvalue.tau = rbind(pvalue.tau,
                       output.tau)
  }
}

data.table::fwrite(pvalue.beta,
                   file = paste0("/gdata01/user/xuhe/SPA-GRM/simulation-2023-03-22/scenario1/SPA.GRMpval/beta", MAF, "pval", type, "-", nrep, ".txt"),
                   row.names = F, col.names = T, quote = F, sep = "\t")

data.table::fwrite(pvalue.tau,
                   file = paste0("/gdata01/user/xuhe/SPA-GRM/simulation-2023-03-22/scenario1/SPA.GRMpval/tau", MAF, "pval", type, "-", nrep, ".txt"),
                   row.names = F, col.names = T, quote = F, sep = "\t")
