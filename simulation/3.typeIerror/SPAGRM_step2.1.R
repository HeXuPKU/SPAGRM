
###### part1 preconditioning ######
# cd "/gdata01/user/xuhe/family_relatedness/simulation-2023-03-22/"
# sbatch -J step1 --mem=5G -t 1-0:0 --array=1-30000 -o log/%A_%a.log --wrap='Rscript /gdata01/user/xuhe/family_relatedness/simulation-2023-03-22/code/type1error_EmpSPAGRM_step1-2023-07-01.R $SLURM_ARRAY_TASK_ID'
args=commandArgs(TRUE)
print(args)
print(sessionInfo())
n.cpu=as.numeric(args)

library(GRAB)

table = data.frame(nrep = rep(1:10000, times = 3),
                   type = rep(c("A", "B", "C"), each = 10000),
                   genofile = rep(c("nSub_25000_nFam_6250", "nSub_25000_nFam_2500", "nSub_50000"), each = 10000))

n.rep = table[n.cpu, "nrep"]; type = table[n.cpu, "type"]; genofile = table[n.cpu, "genofile"]
cat("ncpu is ", n.cpu, "; n.rep is ", n.rep, "; type is ", type, "; genofile is ", genofile, ".\n")

Resid = paste0("/gdata01/user/xuhe/SPA-GRM/simulation-2023-03-22/scenario1/residuals/Resid", type, "-", n.rep, ".txt")

SparseGRMFile = "/gdata02/master_data1/Related_Subjects/Merged_150K_Members/SparseGRM/SparseGRM_0.05.txt"

PairwiseIBDFile = paste0("/gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/PairwiseIBD_", type, ".txt")

if(file.exists(Resid))
{
  Resid = data.table::fread(Resid)
  
  ResidMat.beta = Resid %>% select(SubjID, R_beta) %>% rename("Resid" = "R_beta")
  ResidMat.tau = Resid %>% select(SubjID, R_tau) %>% rename("Resid" = "R_tau")
  
  ResidMatFile.beta = paste0("/gdata01/user/xuhe/SPA-GRM/simulation-2023-03-22/scenario1/residuals/Resid", type, ".beta-", n.rep, ".txt")
  ResidMatFile.tau = paste0("/gdata01/user/xuhe/SPA-GRM/simulation-2023-03-22/scenario1/residuals/Resid", type, ".tau-", n.rep, ".txt")
  
  data.table::fwrite(ResidMat.beta, file = ResidMatFile.beta,
                     row.names = F, col.names = T, quote = F, sep = "\t")
  
  data.table::fwrite(ResidMat.tau, file = ResidMatFile.tau,
                     row.names = F, col.names = T, quote = F, sep = "\t")
  
  obj.precond.beta = SPAGRM.NullModel(ResidMatFile.beta, SparseGRMFile, PairwiseIBDFile)
  
  save(obj.precond.beta,
       file = paste0("/gdata01/user/xuhe/SPA-GRM/simulation-2023-03-22/scenario1/working.data/obj.precond.beta", type, "-", n.rep, ".RData"))
  
  obj.precond.tau = SPAGRM.NullModel(ResidMatFile.tau, SparseGRMFile, PairwiseIBDFile)
  
  save(obj.precond.tau,
       file = paste0("/gdata01/user/xuhe/SPA-GRM/simulation-2023-03-22/scenario1/working.data/obj.precond.tau", type, "-", n.rep, ".RData"))
  
  file.remove(ResidMatFile.beta); file.remove(ResidMatFile.tau)
  
  library(GMMAT)
  ## Load GRM.
  SGRM = data.table::fread(SparseGRMFile)
  SubjID1 = paste0("F4-", rep(1:12500, each=4), "-", 1:4) # 4-members
  SubjID2 = paste0("F10-", rep(1:5000, each=10), "-", 1:10)# 10-members
  SubjID3 = paste0("U-",1:50e3) #unrelated subjects
  SubjID = c(SubjID1, SubjID2, SubjID3)
  SparseGRM = Matrix::sparseMatrix(i = match(SGRM$ID1,SubjID),
                                   j = match(SGRM$ID2,SubjID),
                                   x = SGRM$Value,
                                   symmetric = T)
  rownames(SparseGRM) = colnames(SparseGRM) = SubjID
  
  IDs = Resid$SubjID
  GRM = SparseGRM[IDs, IDs]
  
  model.beta = glmmkin(R_beta ~ NULL, data = Resid, id = "SubjID",
                       kins = GRM, family = gaussian(link = "identity"))
  
  model.tau = glmmkin(R_tau ~ NULL, data = Resid, id = "SubjID",
                      kins = GRM, family = gaussian(link = "identity"))
  
  adj_residuals = data.frame(R_beta = model.beta$residuals,
                             R_tau = model.tau$residuals,
                             SubjID = IDs)
  
  adj.ResidMat.beta = adj_residuals %>% select(SubjID, R_beta) %>% rename("Resid" = "R_beta")
  adj.ResidMat.tau = adj_residuals %>% select(SubjID, R_tau) %>% rename("Resid" = "R_tau")
  
  adj.ResidMatFile.beta = paste0("/gdata01/user/xuhe/SPA-GRM/simulation-2023-03-22/scenario1/adj_residuals/Resid", type, ".beta-", n.rep, ".txt")
  adj.ResidMatFile.tau = paste0("/gdata01/user/xuhe/SPA-GRM/simulation-2023-03-22/scenario1/adj_residuals/Resid", type, ".tau-", n.rep, ".txt")
  
  data.table::fwrite(adj.ResidMat.beta, file = adj.ResidMatFile.beta,
                     row.names = F, col.names = T, quote = F, sep = "\t")
  
  data.table::fwrite(adj.ResidMat.tau, file = adj.ResidMatFile.tau,
                     row.names = F, col.names = T, quote = F, sep = "\t")
  
  if(max(adj.ResidMat.beta$Resid) != min(adj.ResidMat.beta$Resid))
  {
    obj.precond.beta = SPAGRM.NullModel(adj.ResidMatFile.beta, SparseGRMFile, PairwiseIBDFile)
    
    save(obj.precond.beta,
         file = paste0("/gdata01/user/xuhe/SPA-GRM/simulation-2023-03-22/scenario1/adj_working.data/obj.precond.beta", type, "-", n.rep, ".RData"))
  }
  
  if(max(adj.ResidMat.tau$Resid) != min(adj.ResidMat.tau$Resid))
  {
    obj.precond.tau = SPAGRM.NullModel(adj.ResidMatFile.tau, SparseGRMFile, PairwiseIBDFile)
    
    save(obj.precond.tau,
         file = paste0("/gdata01/user/xuhe/SPA-GRM/simulation-2023-03-22/scenario1/adj_working.data/obj.precond.tau", type, "-", n.rep, ".RData"))
  }
  
  file.remove(adj.ResidMatFile.beta); file.remove(adj.ResidMatFile.tau)
}
