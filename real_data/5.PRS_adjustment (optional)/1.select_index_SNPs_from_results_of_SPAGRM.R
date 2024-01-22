
# cd "/gdata01/user/xuhe/family_relatedness/real_data-2023-08-29/"
# sbatch -J PRS --mem=10G -t 1-0:0 --array=1-3 -o log/test_for_PRS_longitudinal_%a.log --wrap='Rscript /gdata01/user/xuhe/family_relatedness/real_data-2023-08-29/4.PRS_new/select_markers_longitudinal_new_step1.R $SLURM_ARRAY_TASK_ID'
args=commandArgs(TRUE)
print(args)
print(sessionInfo())
n.cpu=as.numeric(args)
print(n.cpu)

phenotypes = c("eGFR", "ferritin", "TSH")

pheno = phenotypes[n.cpu]

cat("This is", pheno, "QQ plot combined.\n")

library(dplyr)
library(tidyr)
library(stringr)
library(GRAB)

setwd("/gdata01/user/xuhe/SPA-GRM/real_data-2023-08-29/PRS_new/longitudinal/lead_SNPs/")

load(paste0("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preconding/", pheno, ".SPAGRM.beta.RData"))

spagrm_noGRM = c()

for(chr in 1:22)
{
  cat(chr, "\n")
  
  temp = data.table::fread(paste0("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/temp_pvalues/", pheno, "/chr", chr, ".beta.txt"))
  
  spagrm_noGRM = rbind(spagrm_noGRM,
                       temp %>% drop_na)
}

spagrm_noGRM = spagrm_noGRM %>%
  filter(hwepvalVec > 1e-6) %>%
  filter(AltFreq > 0.01 & AltFreq < 0.99) %>%
  select(Marker, Info, AltFreq, zScore, Pvalue) %>%
  mutate(CHR = as.numeric(str_split(Info, ":", simplify = TRUE)[,1])) %>%
  mutate(POS = as.numeric(str_split(Info, ":", simplify = TRUE)[,2])) %>%
  arrange(Pvalue) %>%
  filter(Pvalue < 5e-8)

for(chr in 1:22)
{
  cat(chr, "\n")
  
  temp_spagrm_noGRM = spagrm_noGRM %>% 
    filter(CHR == chr) %>%
    distinct(Marker, .keep_all = TRUE)
  
  data.table::fwrite(data.frame(temp_spagrm_noGRM$Marker), 
                     file = paste0("/gdata01/user/xuhe/SPA-GRM/real_data-2023-08-29/PRS_new/longitudinal/lead_SNPs/", pheno, ".", chr, "markers.txt"),
                     row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
  
  if(nrow(temp_spagrm_noGRM) == 0)
    next;
  
  data.table::fwrite(temp_spagrm_noGRM,
                     paste0("/gdata01/user/xuhe/SPA-GRM/real_data-2023-08-29/PRS_new/longitudinal/lead_SNPs/", pheno, chr, "pvalues.txt"),
                     row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
  
  temp_genes = GRAB.ReadGeno(GenoFile = paste0("/gdata02/master_data1/UK_Biobank/ukb22828_imp/ukb22828_c", chr, "_b0_v3.bgen"),
                             GenoFileIndex = c(paste0("/gdata02/master_data1/UK_Biobank/ukb22828_imp/ukb_imp_chr", chr, "_v3.bgen.bgi"),
                                               "/gdata02/master_data1/UK_Biobank/ukb22828_imp/ukb22828_c1_b0_v3_s487203.sample"),
                             SampleIDs = SPAGRM.beta$subjData,
                             control = list(IDsToIncludeFile = paste0("/gdata01/user/xuhe/SPA-GRM/real_data-2023-08-29/PRS_new/longitudinal/lead_SNPs/", pheno, ".", chr, "markers.txt")))
  
  temp_genes$markerInfo = temp_genes$markerInfo %>%
    mutate(Info = paste(CHROM, POS, REF, ALT, sep = ":"))
  
  remain_pos = temp_genes$markerInfo$Info %in% temp_spagrm_noGRM$Info
  
  temp_genes$markerInfo = temp_genes$markerInfo[remain_pos,]
  temp_genes$GenoMat = temp_genes$GenoMat[,remain_pos]
  gc()
  
  if(nrow(temp_genes$markerInfo) == 1)
  {
    temp_genes$GenoMat = as.matrix(temp_genes$GenoMat)
    colnames(temp_genes$GenoMat) = temp_genes$markerInfo$ID
  }
  
  bgen_to_plink = temp_genes$GenoMat
  if(ncol(bgen_to_plink) != 1)
  {
    for(i in 1:ncol(bgen_to_plink))
    {
      # cat(i, "\n")
      pos_na = which(bgen_to_plink[,i] != 0 & bgen_to_plink[,i] != 1 & bgen_to_plink[,i] != 2)
      bgen_to_plink[pos_na,i] = -9
    }
  }else
  {
    bgen_to_plink[which(bgen_to_plink != 0 & bgen_to_plink != 1 & bgen_to_plink != 2)] = -9
  }
  gc()
  
  if(ncol(bgen_to_plink) > 1000)
  {
    nparts = ncol(bgen_to_plink) %/% 500
    
    cycle_num = 0; part_num = ncol(bgen_to_plink) %/% nparts
    for(i in 1:(nparts - 1))
    {
      GRAB.makePlink(GenoMat = bgen_to_plink[, (cycle_num + 1):(cycle_num + part_num)],
                     OutputPrefix = paste0("/gdata01/user/xuhe/SPA-GRM/real_data-2023-08-29/PRS_new/longitudinal/lead_SNPs/", pheno, chr, "part", i),
                     # A1 = temp_genes$markerInfo$ALT,
                     # A2 = temp_genes$markerInfo$REF,
                     CHR = rep(chr, ncol(bgen_to_plink[, (cycle_num + 1):(cycle_num + part_num)])),
                     BP = temp_genes$markerInfo$POS[(cycle_num + 1):(cycle_num + part_num)])
      cycle_num = cycle_num + part_num
      gc()
    }
    
    i = i + 1;
    GRAB.makePlink(GenoMat = bgen_to_plink[, (cycle_num + 1):ncol(bgen_to_plink)],
                   OutputPrefix = paste0("/gdata01/user/xuhe/SPA-GRM/real_data-2023-08-29/PRS_new/longitudinal/lead_SNPs/", pheno, chr, "part", i),
                   # A1 = temp_genes$markerInfo$ALT,
                   # A2 = temp_genes$markerInfo$REF,
                   CHR = rep(chr, ncol(bgen_to_plink[, (cycle_num + 1):ncol(bgen_to_plink)])),
                   BP = temp_genes$markerInfo$POS[(cycle_num + 1):ncol(bgen_to_plink)])
    gc()
    
    merge_frame = data.frame(V1 = paste0("/gdata01/user/xuhe/SPA-GRM/real_data-2023-08-29/PRS_new/longitudinal/lead_SNPs/", pheno, chr, "part", 1:nparts, ".ped"),
                             V2 = paste0("/gdata01/user/xuhe/SPA-GRM/real_data-2023-08-29/PRS_new/longitudinal/lead_SNPs/", pheno, chr, "part", 1:nparts, ".map"))
    
    data.table::fwrite(merge_frame, paste0("/gdata01/user/xuhe/SPA-GRM/real_data-2023-08-29/PRS_new/longitudinal/lead_SNPs/", pheno, chr, "merge.txt"),
                       row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
    
    cmd0 = paste0("plink --merge-list ", pheno, chr, "merge.txt --make-bed --out ", pheno, chr)
    system(cmd0)
    
    rm(bgen_to_plink); gc()
    
    file.remove(paste0("/gdata01/user/xuhe/SPA-GRM/real_data-2023-08-29/PRS_new/longitudinal/lead_SNPs/", pheno, chr,"part", 1:nparts, ".ped"))
    file.remove(paste0("/gdata01/user/xuhe/SPA-GRM/real_data-2023-08-29/PRS_new/longitudinal/lead_SNPs/", pheno, chr,"part", 1:nparts, ".map"))
  }else
  {
    GRAB.makePlink(GenoMat = bgen_to_plink,
                   OutputPrefix = paste0("/gdata01/user/xuhe/SPA-GRM/real_data-2023-08-29/PRS_new/longitudinal/lead_SNPs/", pheno, chr),
                   # A1 = temp_genes$markerInfo$ALT,
                   # A2 = temp_genes$markerInfo$REF,
                   CHR = rep(chr, ncol(bgen_to_plink)),
                   BP = temp_genes$markerInfo$POS)
    
    rm(bgen_to_plink); gc()
    
    cmd1 = paste0("plink --file ", pheno, chr, " --make-bed --out ", pheno, chr)
    system(cmd1)
  }
  
  cmd2 = paste0("plink --bfile ", pheno, chr, " --clump-p1 5e-8 --clump-r2 0.01 --clump-kb 2000 --clump ", pheno, chr, "pvalues.txt --clump-snp-field Marker --clump-field Pvalue --out ", pheno, chr)
  system(cmd2)
  
  clumped_SNPs = data.table::fread(paste0("/gdata01/user/xuhe/SPA-GRM/real_data-2023-08-29/PRS_new/longitudinal/lead_SNPs/", pheno, chr,".clumped"), select = 1:5)
  
  remain_pos = temp_genes$markerInfo$ID %in% clumped_SNPs$SNP
  
  temp_genes$markerInfo = temp_genes$markerInfo[remain_pos,]
  temp_genes$GenoMat = temp_genes$GenoMat[,remain_pos]
  gc()
  
  if(nrow(temp_genes$markerInfo) == 1)
  {
    temp_genes$GenoMat = as.matrix(temp_genes$GenoMat)
    colnames(temp_genes$GenoMat) = temp_genes$markerInfo$ID
  }
  
  if(sum(remain_pos) != nrow(temp_genes$markerInfo))
    stop("sum(remain_pos) != nrow(temp_genes$markerInfo)")
  
  mfi_file = data.table::fread(paste0("/gdata02/master_data1/UK_Biobank/ukb22828_imp/ukb_mfi_chr", chr, "_v3.txt"))
  
  mfi_file = mfi_file %>% mutate(CHR = ifelse(chr < 10, paste0(0, chr), chr)) %>%
    mutate(Info = paste(CHR, V3, V4, V5, sep = ":"))
  
  snpmask = mfi_file$Info %in% temp_genes$markerInfo$Info
  
  data.table::fwrite(as.data.frame(snpmask), 
                     file = paste0("/gdata01/user/xuhe/SPA-GRM/real_data-2023-08-29/PRS_new/longitudinal/lead_SNPs/", pheno, ".", chr, ".snpmask.txt"),
                     row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
  
  file.remove(paste0("/gdata01/user/xuhe/SPA-GRM/real_data-2023-08-29/PRS_new/longitudinal/lead_SNPs/", pheno, chr, ".ped"))
  file.remove(paste0("/gdata01/user/xuhe/SPA-GRM/real_data-2023-08-29/PRS_new/longitudinal/lead_SNPs/", pheno, chr, ".map"))
  file.remove(paste0("/gdata01/user/xuhe/SPA-GRM/real_data-2023-08-29/PRS_new/longitudinal/lead_SNPs/", pheno, chr, ".bed"))
  file.remove(paste0("/gdata01/user/xuhe/SPA-GRM/real_data-2023-08-29/PRS_new/longitudinal/lead_SNPs/", pheno, chr, ".bim"))
  file.remove(paste0("/gdata01/user/xuhe/SPA-GRM/real_data-2023-08-29/PRS_new/longitudinal/lead_SNPs/", pheno, chr, ".fam"))
  gc()
}
