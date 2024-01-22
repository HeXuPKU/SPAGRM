
###### part1 preconditioning ######
# library(GRAB)
# 
# ResidMatFile.beta = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/null_fitter/Hct_beta_resid.txt"
# ResidMatFile.tau = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/null_fitter/Hct_tau_resid.txt"
# 
# SparseGRMFile = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/GRM_IBD/UKB410K_sparseGRM_0.05.txt"
# 
# PairwiseIBDFile = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/GRM_IBD/UKB410K_pairwiseIBD.txt"
# 
# SPAGRM.beta = SPAGRM.NullModel(ResidMatFile.beta,
#                                SparseGRMFile,
#                                PairwiseIBDFile)
# 
# save(SPAGRM.beta,
#      file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preconding/Hct.SPAGRM.beta.RData")
# 
# SPAGRM.tau = SPAGRM.NullModel(ResidMatFile.tau,
#                               SparseGRMFile,
#                               PairwiseIBDFile)
# 
# save(SPAGRM.tau,
#      file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preconding/Hct.SPAGRM.tau.RData")


###### part2 analyzing ######
# cd "/gdata01/user/xuhe/family_relatedness/real_data-2022-12-29/"
# sbatch -J Hctstep2 --mem=5G -t 3-0:0 --array=1-22 -o log/%A_%a.log --wrap='Rscript /gdata01/user/xuhe/family_relatedness/real_data-2022-12-29/4.EmpSPA/Hct_SPAGRM.R $SLURM_ARRAY_TASK_ID'
args=commandArgs(TRUE)
print(args)
print(sessionInfo())
n.cpu=as.numeric(args)
print(n.cpu)

cat("This is Hct.")

library(GRAB)

GenoFile = paste0("/gdata02/master_data1/UK_Biobank/ukb22828_imp/ukb22828_c", n.cpu, "_b0_v3.bgen")
bgenbgiFile = paste0("/gdata02/master_data1/UK_Biobank/ukb22828_imp/ukb_imp_chr", n.cpu, "_v3.bgen.bgi")
bgensampleFile = "/gdata02/master_data1/UK_Biobank/ukb22828_imp/ukb22828_c1_b0_v3_s487203.sample"

OutputDir = paste0("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/temp_pvalues/Hct/chr", n.cpu)

SNPprefilter = data.table::fread(paste0("/gdata02/master_data1/UK_Biobank/ukb22828_imp/ukb_mfi_chr", n.cpu, "_v3.txt"))
SNPprefilter = SNPprefilter %>% filter(V6 > 5e-5 & V8 >= 0.6)

IDsToIncludeFile = paste0(OutputDir, ".SNPprefilter.txt")

data.table::fwrite(data.table::data.table(SNPprefilter$V2), file = IDsToIncludeFile,
                   col.names = F, row.names = F, quote = F, sep = "\t")

precond.beta = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preconding/Hct.SPAGRM.beta.RData"
load(precond.beta)

GRAB.Marker(objNull = SPAGRM.beta,
            GenoFile = GenoFile,
            GenoFileIndex = c(bgenbgiFile, bgensampleFile),
            OutputFile = paste0(OutputDir, ".beta.txt"),
            control = list(IDsToIncludeFile = IDsToIncludeFile,
                           SPA_Cutoff = 2,
                           zeta = -0.02,
                           tol = 1e-4,
                           impute_method = "mean",
                           missing_cutoff = 0.05,
                           min_maf_marker = 0.0002,
                           min_mac_marker = 20,
                           nMarkersEachChunk = 2500))


precond.tau = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preconding/Hct.SPAGRM.tau.RData"
load(precond.tau)

GRAB.Marker(objNull = SPAGRM.tau,
            GenoFile = GenoFile,
            GenoFileIndex = c(bgenbgiFile, bgensampleFile),
            OutputFile = paste0(OutputDir, ".tau.txt"),
            control = list(IDsToIncludeFile = IDsToIncludeFile,
                           SPA_Cutoff = 2,
                           zeta = 0,
                           tol = 1e-5,
                           impute_method = "mean",
                           missing_cutoff = 0.05,
                           min_maf_marker = 0.0002,
                           min_mac_marker = 20,
                           nMarkersEachChunk = 2500))
