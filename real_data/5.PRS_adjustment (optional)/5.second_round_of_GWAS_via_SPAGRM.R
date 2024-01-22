
###### part1 preconditioning ######
# cd /gdata01/user/xuhe/family_relatedness/real_data-2023-08-29/
# sbatch -J PRS --mem=8G -t 5-0:0 --array=1-66 -o log/test_PRS_longitudinal_%a.log --wrap='Rscript /gdata01/user/xuhe/family_relatedness/real_data-2023-08-29/4.PRS_new/SPAGRM_PRS_longitudinal.R $SLURM_ARRAY_TASK_ID'
args=commandArgs(TRUE)
print(args)
print(sessionInfo())
n.cpu=as.numeric(args)
print(n.cpu)

Sys.time()

pheno_df = data.frame(nrep = 1:66,
                      pheno = rep(c("eGFR", "ferritin", "TSH"), each = 22),
                      chr = rep(1:22, times = 3))

pheno = pheno_df$pheno[n.cpu]; chr = pheno_df$chr[n.cpu]

library(GRAB)

ResidMatFile.beta = paste0("/gdata01/user/xuhe/SPA-GRM/real_data-2023-08-29/PRS_new/longitudinal/ResidMat/", pheno, ".", chr, ".betaresid.txt")

ResidMatFile.tau = paste0("/gdata01/user/xuhe/SPA-GRM/real_data-2023-08-29/PRS_new/longitudinal/ResidMat/", pheno, ".", chr, ".tauresid.txt")

SparseGRMFile = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/GRM_IBD/UKB410K_sparseGRM_0.05.txt"

PairwiseIBDFile = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/GRM_IBD/UKB410K_pairwiseIBD.txt"

SPAGRM.beta = SPAGRM.NullModel(ResidMatFile.beta,
                               SparseGRMFile,
                               PairwiseIBDFile)

save(SPAGRM.beta,
     file = paste0("/gdata01/user/xuhe/SPA-GRM/real_data-2023-08-29/PRS_new/longitudinal/ResidMat/", pheno, ".", chr, ".beta.RData"))

SPAGRM.tau = SPAGRM.NullModel(ResidMatFile.tau,
                              SparseGRMFile,
                              PairwiseIBDFile)

save(SPAGRM.tau,
     file = paste0("/gdata01/user/xuhe/SPA-GRM/real_data-2023-08-29/PRS_new/longitudinal/ResidMat/", pheno, ".", chr, ".tau.RData"))

Sys.time()

GenoFile = paste0("/gdata02/master_data1/UK_Biobank/ukb22828_imp/ukb22828_c", chr, "_b0_v3.bgen")
bgenbgiFile = paste0("/gdata02/master_data1/UK_Biobank/ukb22828_imp/ukb_imp_chr", chr, "_v3.bgen.bgi")
bgensampleFile = "/gdata02/master_data1/UK_Biobank/ukb22828_imp/ukb22828_c1_b0_v3_s487203.sample"

OutputDir = paste0("/gdata01/user/xuhe/SPA-GRM/real_data-2023-08-29/PRS_new/longitudinal/pvalues/", pheno, ".", chr)

SNPprefilter = data.table::fread(paste0("/gdata02/master_data1/UK_Biobank/ukb22828_imp/ukb_mfi_chr", chr, "_v3.txt"))
SNPprefilter = SNPprefilter %>% filter(V6 > 5e-3 & V8 >= 0.6)

IDsToIncludeFile = paste0(OutputDir, ".SNPprefilter.txt")

data.table::fwrite(data.table::data.table(SNPprefilter$V2), file = IDsToIncludeFile,
                   col.names = F, row.names = F, quote = F, sep = "\t")

GRAB.Marker(objNull = SPAGRM.beta,
            GenoFile = GenoFile,
            GenoFileIndex = c(bgenbgiFile, bgensampleFile),
            OutputFile = paste0(OutputDir, ".pvalue.beta.txt"),
            control = list(IDsToIncludeFile = IDsToIncludeFile,
                           SPA_Cutoff = 2,
                           zeta = 0,
                           tol = 1e-4,
                           impute_method = "mean",
                           missing_cutoff = 0.05,
                           min_maf_marker = 2e-4,
                           min_mac_marker = 20,
                           nMarkersEachChunk = 10000))

# GRAB.Marker(objNull = SPAGRM.tau,
#             GenoFile = GenoFile,
#             GenoFileIndex = c(bgenbgiFile, bgensampleFile),
#             OutputFile = paste0(OutputDir, ".pvalue.tau.txt"),
#             control = list(IDsToIncludeFile = IDsToIncludeFile,
#                            SPA_Cutoff = 2,
#                            zeta = 0,
#                            tol = 1e-5,
#                            impute_method = "mean",
#                            missing_cutoff = 0.05,
#                            min_maf_marker = 2e-4,
#                            min_mac_marker = 20,
#                            nMarkersEachChunk = 10000))

file.remove(IDsToIncludeFile)