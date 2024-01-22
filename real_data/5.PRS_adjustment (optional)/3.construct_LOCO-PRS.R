
# cd "/gdata01/user/xuhe/family_relatedness/real_data-2023-08-29/"
# sbatch -J PRS --mem=10G -t 1-0:0 --array=1-3 -o log/test_for_PRS_longitudinal_%a.log --wrap='Rscript /gdata01/user/xuhe/family_relatedness/real_data-2023-08-29/4.PRS_new/select_markers_longitudinal_new_step3.R $SLURM_ARRAY_TASK_ID'
args=commandArgs(TRUE)
print(args)
print(sessionInfo())
n.cpu=as.numeric(args)
print(n.cpu)

phenotypes = c("eGFR", "ferritin", "TSH")

pheno = phenotypes[n.cpu]

cat("This is", pheno, ".\n")

load(paste0("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preconding/", pheno, ".SPAGRM.beta.RData"))

library(dplyr)
library(tidyr)
library(stringr)
library(GRAB)

for(chr in 1:22)
{
  cat(chr, "\n")
  
  if(!file.exists(paste0("/gdata01/user/xuhe/SPA-GRM/real_data-2023-08-29/PRS_new/longitudinal/lead_SNPs/", pheno, ".", chr, ".effectsize.txt")))
    next;
  
  # 8 columns if using the bgen file as input and set disable_wsvar = true
  # V1: chr; V2: pos; V3: snpid; V4: ref; V5: alt; V6: betaeffect; V7: betastderr (maybe it is betase); V8: betapval
  effectsize = data.table::fread(paste0("/gdata01/user/xuhe/SPA-GRM/real_data-2023-08-29/PRS_new/longitudinal/lead_SNPs/", pheno, ".", chr, ".effectsize.txt"), header = FALSE)
  effectsize = effectsize %>% replace_na(list(V7 = 0))
  
  data.table::fwrite(data.frame(effectsize$V3),
                     file = paste0("/gdata01/user/xuhe/SPA-GRM/real_data-2023-08-29/PRS_new/longitudinal/lead_SNPs/", pheno, ".", chr, ".leadsnp.txt"),
                     row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
  
  temp_genes = GRAB.ReadGeno(GenoFile = paste0("/gdata02/master_data1/UK_Biobank/ukb22828_imp/ukb22828_c", chr, "_b0_v3.bgen"),
                             GenoFileIndex = c(paste0("/gdata02/master_data1/UK_Biobank/ukb22828_imp/ukb_imp_chr", chr, "_v3.bgen.bgi"),
                                               "/gdata02/master_data1/UK_Biobank/ukb22828_imp/ukb22828_c1_b0_v3_s487203.sample"),
                             SampleIDs = SPAGRM.beta$subjData,
                             control = list(IDsToIncludeFile = paste0("/gdata01/user/xuhe/SPA-GRM/real_data-2023-08-29/PRS_new/longitudinal/lead_SNPs/", pheno, ".", chr, ".leadsnp.txt")))
  
  if(nrow(effectsize) != nrow(temp_genes$markerInfo))
  {
    temp_effect = effectsize %>% mutate(Info = paste(V2, V5, V6, sep = ":"))
    temp_marker = temp_genes$markerInfo %>% mutate(Info = paste(POS, REF, ALT, sep = ":"))
    
    remain_pos = temp_marker$Info %in% temp_effect$Info
    
    temp_genes$markerInfo = temp_genes$markerInfo[remain_pos,]
    temp_genes$GenoMat = temp_genes$GenoMat[,remain_pos]
    gc()
  }
  
  PRS_single_chr = t(effectsize$V7 %*% t(temp_genes$GenoMat)) %>%
    as.data.frame() %>% tibble::rownames_to_column(var = "FID") %>% rename(PRS = V1)
  
  data.table::fwrite(PRS_single_chr, paste0("/gdata01/user/xuhe/SPA-GRM/real_data-2023-08-29/PRS_new/longitudinal/PRS/", pheno, chr, ".single.txt"),
                     row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
}


tot_PRS = c()

for(chr in 1:22)
{
  cat(chr, "\n")
  
  single_PRS_file = paste0("/gdata01/user/xuhe/SPA-GRM/real_data-2023-08-29/PRS_new/longitudinal/PRS/", pheno, chr, ".single.txt")
  
  if(file.exists(single_PRS_file))
  {
    temp_PRS = data.table::fread(single_PRS_file)
    
    tot_PRS = cbind(tot_PRS, temp_PRS[,2])
  }
}

tot_PRS = rowSums(tot_PRS)


for(chr in 1:22)
{
  cat(chr, "\n")
  
  single_PRS_file = paste0("/gdata01/user/xuhe/SPA-GRM/real_data-2023-08-29/PRS_new/longitudinal/PRS/", pheno, chr, ".single.txt")
  
  if(file.exists(single_PRS_file))
  {
    LOCO_PRS = data.table::fread(single_PRS_file)
    
    LOCO_PRS$PRS = tot_PRS - LOCO_PRS$PRS
  }else
  {
    LOCO_PRS = data.frame(FID = SPAGRM.beta$subjData, PRS = tot_PRS)
  }
  
  data.table::fwrite(LOCO_PRS, paste0("/gdata01/user/xuhe/SPA-GRM/real_data-2023-08-29/PRS_new/longitudinal/PRS/", pheno, chr, ".LOCO.txt"),
                     row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
}
