
# This code is similar to /gdata02/master_data1/UK_Biobank/ukb23155_b0_v1_s200604/SimulationUsingWES/makeData.R
#  and /gdata02/master_data1/UK_Biobank/ukb22418_b0_v2_s488180/FilesForSimulation/makeData1.R.
#
# Three sample settings:
#   (A). 25,000 unrelated subjects and 25,000 (6,250×4) related subjects;
#   (B). 25,000 unrelated subjects and 25,000 (2,500×10) related subjects;
#   (C). 50,000 unrelated subjects.
# 
# Two MAF settings:
#   common variants: 100K with MAFs range from 0.05 to 0.5;
#   rare variants: 100K with MAFs range from 0.0002 to 0.05.


############################################################### common variants.
# this part of code is similar to makeData1.R.

### randomly select 120K SNPs for common variants.
PlinkPrefix = "/gdata02/master_data1/UK_Biobank/ukb22418_b0_v2_s488180/PLINK_for_450K_WES/ukb_allchr_v2_newID_passedQC_white.British_geno0.05_poly_500_50_0.2.pruned_woWithdrawl.chr.1_22_maf_0.05_WES450K"
bimFile = paste0(PlinkPrefix, ".bim")
bimData = data.table::fread(bimFile)

set.seed(123)
posInBIM = sort(sample(nrow(bimData), 12e4))

SNPs = bimData$V2[posInBIM]

for(group in 1:100)
{
  TempSNPs = SNPs[1:1200+(group-1)*1200]
  FileSplit = paste0("/gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/tempsplitdFiles/commonSNPs_group_",group,".txt")
  data.table::fwrite(data.table::data.table(TempSNPs), FileSplit,
                     row.names = F, col.names = F, quote = F, sep = "\t")
}

### (A). 25,000 unrelated subjects and 25,000 (6,250×4) related subjects
# first for common variants
library(tidyr)
library(dplyr)
library(GRAB)

PlinkPrefix = "/gdata02/master_data1/UK_Biobank/ukb22418_b0_v2_s488180/PLINK_for_450K_WES/ukb_allchr_v2_newID_passedQC_white.British_geno0.05_poly_500_50_0.2.pruned_woWithdrawl.chr.1_22_maf_0.05_WES450K"
DirSplit = "/gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/tempsplitdFiles/"

bimFile = paste0(PlinkPrefix, ".bim")
famFile = paste0(PlinkPrefix, ".fam")
bedFile = paste0(PlinkPrefix, ".bed")

bimData = data.table::fread(bimFile)
famData = data.table::fread(famFile)

SubjData = data.table::fread("/gdata02/master_data1/UK_Biobank/ukb23155_b0_v1_s200604/SimulationUsingWES/Random_100e3_Subjects.sample")  
SubjData = SubjData %>% filter(V1 %in% famData$V2)

nSub = 25e3
nFam = 6250

# please run the code below on the slurm.
# cd "/gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/"
# sbatch -J simugeno --mem=5G -t 1-0:0 --array=1-100 -o log/%A_%a.log --wrap='Rscript /gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/SimulatedGenotypeUsingWES-2023-02-16XH.R $SLURM_ARRAY_TASK_ID'
args = commandArgs(TRUE)
print(args)
print(sessionInfo())
group = as.numeric(args)

IDsToIncludeFile = paste0(DirSplit, "commonSNPs_group_", group, ".txt");
GenoList = GRAB.SimuGMatFromGenoFile(nFam = nFam,
                                     nSub = nSub,
                                     FamMode = "4-members",
                                     GenoFile = bedFile,
                                     SampleIDs = SubjData$V1,
                                     control = list(IDsToIncludeFile = IDsToIncludeFile))

GenoList$GenoMat[is.na(GenoList$GenoMat)] = -9

# update FID and IID
rownames(GenoList$GenoMat) = gsub("Subj-", "U-", rownames(GenoList$GenoMat))
rownames(GenoList$GenoMat) = gsub("_", "-", rownames(GenoList$GenoMat))
rownames(GenoList$GenoMat) = gsub("f", "F4-", rownames(GenoList$GenoMat))

OutputPrefix = paste0(DirSplit, "nSub_25000_nFam_6250/", "SNPs_group_", group)

GRAB.makePlink(GenoList$GenoMat, OutputPrefix)

cmd1 = paste0("plink --file ", OutputPrefix, " --make-bed --out ", OutputPrefix);
cmd2 = paste0("rm ", OutputPrefix, ".map ", OutputPrefix, ".ped");

system(cmd1)
system(cmd2)

bimFileGroup = paste0(OutputPrefix, ".bim")
markersInGroup = paste0(DirSplit, "commonSNPs_group_", group, ".txt")
markersInGroup = data.table::fread(markersInGroup, header = F)

pos = match(markersInGroup$V1, bimData$V2)
bimDataInGroup = bimData[pos,]
data.table::fwrite(bimDataInGroup, bimFileGroup, 
                   sep = "\t", row.names = F, col.names = F, quote = F)

# plink command.
if(FALSE)
{
  cd "/gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/tempsplitdFiles/nSub_25000_nFam_6250/"
  plink --merge-list ../group.mergelist --out nSub_25000_nFam_6250 --make-bed
  plink --bfile nSub_25000_nFam_6250 --freq --out nSub_25000_nFam_6250
  plink --bfile nSub_25000_nFam_6250 --missing --out nSub_25000_nFam_6250
}

# QC for missing rates <= 0.05; MAF >= 0.05.
library(dplyr)
DirSplit = "/gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/tempsplitdFiles/"
missing = data.table::fread(paste0(DirSplit, "nSub_25000_nFam_6250/nSub_25000_nFam_6250.lmiss"))
freq = data.table::fread(paste0(DirSplit, "nSub_25000_nFam_6250/nSub_25000_nFam_6250.frq"))

allinfo = inner_join(missing, freq, by = c("CHR", "SNP"))

filterinfo = allinfo %>% filter(MAF >= 0.05 & F_MISS < 0.05) # remaining 115760 SNPs.

SNPIDs = sample(filterinfo$SNP, 1e5)
data.table::fwrite(data.table::data.table(SNPIDs),
                   "/gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/nSub_25000_nFam_6250common/SNPIDs.txt",
                   row.names = F, col.names = F, quote = F, sep = "\t")

for(i in 1:10){
  SNPIDi = SNPIDs[1:10000+10000*(i-1)]
  data.table::fwrite(data.table::data.table(SNPIDi),
                     paste0("/gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/nSub_25000_nFam_6250common/SNPs_group", i, ".txt"),
                     row.names = F, col.names = F, quote = F, sep = "\t")
}

# plink command.
if(FALSE)
{
  plink --bfile /gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/tempsplitdFiles/nSub_25000_nFam_6250/nSub_25000_nFam_6250 \
  --extract /gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/nSub_25000_nFam_6250common/SNPIDs.txt \
  --allow-no-vars --make-bed --out /gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/nSub_25000_nFam_6250common/nSub_25000_nFam_6250
}


### (B). 25,000 unrelated subjects and 25,000 (2,500×4) related subjects
# first for common variants
library(tidyr)
library(dplyr)
library(GRAB)

PlinkPrefix = "/gdata02/master_data1/UK_Biobank/ukb22418_b0_v2_s488180/PLINK_for_450K_WES/ukb_allchr_v2_newID_passedQC_white.British_geno0.05_poly_500_50_0.2.pruned_woWithdrawl.chr.1_22_maf_0.05_WES450K"
DirSplit = "/gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/tempsplitdFiles/"

bimFile = paste0(PlinkPrefix, ".bim")
famFile = paste0(PlinkPrefix, ".fam")
bedFile = paste0(PlinkPrefix, ".bed")

bimData = data.table::fread(bimFile)
famData = data.table::fread(famFile)

SubjData = data.table::fread("/gdata02/master_data1/UK_Biobank/ukb23155_b0_v1_s200604/SimulationUsingWES/Random_100e3_Subjects.sample")  
SubjData = SubjData %>% filter(V1 %in% famData$V2)

nSub = 25e3
nFam = 2500

# please run the code below on the slurm.
# cd "/gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/"
# sbatch -J simugeno --mem=5G -t 1-0:0 --array=1-100 -o log/%A_%a.log --wrap='Rscript /gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/SimulatedGenotypeUsingWES-2023-02-16XH.R $SLURM_ARRAY_TASK_ID'
args = commandArgs(TRUE)
print(args)
print(sessionInfo())
group = as.numeric(args)

IDsToIncludeFile = paste0(DirSplit, "commonSNPs_group_", group, ".txt");
GenoList = GRAB.SimuGMatFromGenoFile(nFam = nFam,
                                     nSub = nSub,
                                     FamMode = "10-members",
                                     GenoFile = bedFile,
                                     SampleIDs = SubjData$V1,
                                     control = list(IDsToIncludeFile = IDsToIncludeFile))

GenoList$GenoMat[is.na(GenoList$GenoMat)] = -9

# update FID and IID
rownames(GenoList$GenoMat) = gsub("Subj-", "U-", rownames(GenoList$GenoMat))
rownames(GenoList$GenoMat) = gsub("_", "-", rownames(GenoList$GenoMat))
rownames(GenoList$GenoMat) = gsub("f", "F10-", rownames(GenoList$GenoMat))

OutputPrefix = paste0(DirSplit, "nSub_25000_nFam_2500/", "SNPs_group_", group)

GRAB.makePlink(GenoList$GenoMat, OutputPrefix)

cmd1 = paste0("plink --file ", OutputPrefix, " --make-bed --out ", OutputPrefix);
cmd2 = paste0("rm ", OutputPrefix, ".map ", OutputPrefix, ".ped");

system(cmd1)
system(cmd2)

bimFileGroup = paste0(OutputPrefix, ".bim")
markersInGroup = paste0(DirSplit, "commonSNPs_group_", group, ".txt")
markersInGroup = data.table::fread(markersInGroup, header = F)

pos = match(markersInGroup$V1, bimData$V2)
bimDataInGroup = bimData[pos,]
data.table::fwrite(bimDataInGroup, bimFileGroup, 
                   sep = "\t", row.names = F, col.names = F, quote = F)

# plink command.
if(FALSE)
{
  cd "/gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/tempsplitdFiles/nSub_25000_nFam_2500/"
  plink --merge-list ../group.mergelist --out nSub_25000_nFam_2500 --make-bed
  plink --bfile nSub_25000_nFam_2500 --freq --out nSub_25000_nFam_2500
  plink --bfile nSub_25000_nFam_2500 --missing --out nSub_25000_nFam_2500
}

# QC for missing rates <= 0.05; MAF >= 0.05.
library(dplyr)
DirSplit = "/gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/tempsplitdFiles/"
missing = data.table::fread(paste0(DirSplit, "nSub_25000_nFam_2500/nSub_25000_nFam_2500.lmiss"))
freq = data.table::fread(paste0(DirSplit, "nSub_25000_nFam_2500/nSub_25000_nFam_2500.frq"))

allinfo = inner_join(missing, freq, by = c("CHR", "SNP"))

filterinfo = allinfo %>% filter(MAF >= 0.05 & F_MISS < 0.05) # remaining 115760 SNPs.

SNPIDs = sample(filterinfo$SNP, 1e5)
data.table::fwrite(data.table::data.table(SNPIDs),
                   "/gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/nSub_25000_nFam_2500common/SNPIDs.txt",
                   row.names = F, col.names = F, quote = F, sep = "\t")

for(i in 1:10){
  SNPIDi = SNPIDs[1:10000+10000*(i-1)]
  data.table::fwrite(data.table::data.table(SNPIDi),
                     paste0("/gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/nSub_25000_nFam_2500common/SNPs_group", i, ".txt"),
                     row.names = F, col.names = F, quote = F, sep = "\t")
}

# plink command.
if(FALSE)
{
  plink --bfile /gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/tempsplitdFiles/nSub_25000_nFam_2500/nSub_25000_nFam_2500 \
  --extract /gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/nSub_25000_nFam_2500common/SNPIDs.txt \
  --allow-no-vars --make-bed --out /gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/nSub_25000_nFam_2500common/nSub_25000_nFam_2500
}


### (C). 50,000 unrelated subjects
# first for common variants
library(tidyr)
library(dplyr)
library(GRAB)

PlinkPrefix = "/gdata02/master_data1/UK_Biobank/ukb22418_b0_v2_s488180/PLINK_for_450K_WES/ukb_allchr_v2_newID_passedQC_white.British_geno0.05_poly_500_50_0.2.pruned_woWithdrawl.chr.1_22_maf_0.05_WES450K"
DirSplit = "/gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/tempsplitdFiles/"

bimFile = paste0(PlinkPrefix, ".bim")
famFile = paste0(PlinkPrefix, ".fam")
bedFile = paste0(PlinkPrefix, ".bed")

bimData = data.table::fread(bimFile)
famData = data.table::fread(famFile)

SubjData = data.table::fread("/gdata02/master_data1/UK_Biobank/ukb23155_b0_v1_s200604/SimulationUsingWES/Random_100e3_Subjects.sample")  
SubjData = SubjData %>% filter(V1 %in% famData$V2)

nSub = 50e3
nFam = 0

# please run the code below on the slurm.
# cd "/gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/"
# sbatch -J simugeno --mem=5G -t 1-0:0 --array=1-100 -o log/%A_%a.log --wrap='Rscript /gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/SimulatedGenotypeUsingWES-2023-02-16XH.R $SLURM_ARRAY_TASK_ID'
args = commandArgs(TRUE)
print(args)
print(sessionInfo())
group = as.numeric(args)

IDsToIncludeFile = paste0(DirSplit, "commonSNPs_group_", group, ".txt");
GenoList = GRAB.SimuGMatFromGenoFile(nFam = nFam,
                                     nSub = nSub,
                                     FamMode = "4-members",
                                     GenoFile = bedFile,
                                     SampleIDs = SubjData$V1,
                                     control = list(IDsToIncludeFile = IDsToIncludeFile))

GenoList$GenoMat[is.na(GenoList$GenoMat)] = -9

# update FID and IID
rownames(GenoList$GenoMat) = gsub("Subj-", "U-", rownames(GenoList$GenoMat))
rownames(GenoList$GenoMat) = gsub("_", "-", rownames(GenoList$GenoMat))

OutputPrefix = paste0(DirSplit, "nSub_50000/", "SNPs_group_", group)

GRAB.makePlink(GenoList$GenoMat, OutputPrefix)

cmd1 = paste0("plink --file ", OutputPrefix, " --make-bed --out ", OutputPrefix);
cmd2 = paste0("rm ", OutputPrefix, ".map ", OutputPrefix, ".ped");

system(cmd1)
system(cmd2)

bimFileGroup = paste0(OutputPrefix, ".bim")
markersInGroup = paste0(DirSplit, "commonSNPs_group_", group, ".txt")
markersInGroup = data.table::fread(markersInGroup, header = F)

pos = match(markersInGroup$V1, bimData$V2)
bimDataInGroup = bimData[pos,]
data.table::fwrite(bimDataInGroup, bimFileGroup, 
                   sep = "\t", row.names = F, col.names = F, quote = F)

# plink command.
if(FALSE)
{
  cd "/gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/tempsplitdFiles/nSub_50000/"
  plink --merge-list ../group.mergelist --out nSub_50000 --make-bed
  plink --bfile nSub_50000 --freq --out nSub_50000
  plink --bfile nSub_50000 --missing --out nSub_50000
}

# QC for missing rates <= 0.05; MAF >= 0.05.
library(dplyr)
DirSplit = "/gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/tempsplitdFiles/"
missing = data.table::fread(paste0(DirSplit, "nSub_50000/nSub_50000.lmiss"))
freq = data.table::fread(paste0(DirSplit, "nSub_50000/nSub_50000.frq"))

allinfo = inner_join(missing, freq, by = c("CHR", "SNP"))

filterinfo = allinfo %>% filter(MAF >= 0.05 & F_MISS < 0.05) # remaining 115760 SNPs.

SNPIDs = sample(filterinfo$SNP, 1e5)
data.table::fwrite(data.table::data.table(SNPIDs),
                   "/gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/nSub_50000common/SNPIDs.txt",
                   row.names = F, col.names = F, quote = F, sep = "\t")

for(i in 1:10){
  SNPIDi = SNPIDs[1:10000+10000*(i-1)]
  data.table::fwrite(data.table::data.table(SNPIDi),
                     paste0("/gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/nSub_50000common/SNPs_group", i, ".txt"),
                     row.names = F, col.names = F, quote = F, sep = "\t")
}

# plink command.
if(FALSE)
{
  plink --bfile /gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/tempsplitdFiles/nSub_50000/nSub_50000 \
  --extract /gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/nSub_50000common/SNPIDs.txt \
  --allow-no-vars --make-bed --out /gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/nSub_50000common/nSub_50000
}


################################################################# rare variants.
# this part of code is similar to makeData.R.

# Randomly select 1e4 genes from chromosomes 1-22
library(tidyr)
library(dplyr)

GeneInfo = data.table::fread("/gdata02/master_data1/UK_Biobank/ukb23155_b0_v1_s200604/SimulationUsingWES/SimplifiedMarkerAnnoInfo/gene_count_simplified.csv")

GeneInfo1 = GeneInfo %>% filter(n > 100) %>% 
  filter(!grepl(";", Gene.refGene))  

set.seed(23333)
n = nrow(GeneInfo1)
pos = sample(n, 10000)
GeneInfo2 = GeneInfo1[pos,]

GeneSelected = GeneInfo2$Gene.refGene

summary(GeneInfo2$n);sum(GeneInfo2$n)

data.table::fwrite(GeneInfo2,
                   "/gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/tempsplitdFiles/raregene_selected.csv",
                   row.names = F, col.names = T, sep = ",", quote = F)

###
library(tidyr)
library(dplyr)

GeneInfo2 = data.table::fread("/gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/tempsplitdFiles/raregene_selected.csv")

for(c in 1:22)
{
  cat("chr:\t", c, "\n")
  
  tempGeneInfo2 = GeneInfo2 %>% filter(chr == c)
  tempGenes = tempGeneInfo2$Gene.refGene
  
  SimpleFile = paste0("/gdata02/master_data1/UK_Biobank/ukb23155_b0_v1_s200604/SimulationUsingWES/SimplifiedMarkerAnnoInfo/UKBexomeOQFE_chr",
                      c,".avoutput.hg38_multianno_simplied.csv")
  SimpleInfo = data.table::fread(SimpleFile)
  
  for(g in tempGenes){
    # cat("gene:\t", g, "\n")
    tempSimpleInfo = SimpleInfo %>% filter(Gene.refGene == g)
    data.table::fwrite(tempSimpleInfo, 
                       paste0("/gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/tempsplitdFiles/GeneLevelGenotype/", g, ".csv"))
  }
}

###
library(dplyr)
library(tidyr)

genes = list.files("/gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/tempsplitdFiles/GeneLevelGenotype",
                   full.names = T)

markerInfo = data.table::data.table()
for(g in genes)
{
  cat("Gene:\t",g,"\n")
  data0 = data.table::fread(g)
  markerInfo = rbind(markerInfo, data0)
}

for(chr in 1:22)
{
  cat("chr:\t",chr,"\n")
  markerInfoTemp = markerInfo %>% filter(Chr == chr)
  markerInfoTemp %>% 
    select(Otherinfo1) %>% 
    data.table::fwrite(paste0("/gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/tempsplitdFiles/GeneLevelPLINK/chr",chr,".marker"))
}

Random_100e3_Subjects = data.table::fread("/gdata02/master_data1/UK_Biobank/ukb23155_b0_v1_s200604/SimulationUsingWES/Random_100e3_Subjects.sample")
Random_50e3_Subjects = Random_100e3_Subjects[sample(100e3, 50e3),]
data.table::fwrite(Random_50e3_Subjects, 
                   "/gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/tempsplitdFiles/Random_50e3_Subjects.sample",
                   col.names = F, row.names = F, quote = F, sep = "\t")

# plink command
if(FALSE){
  cd /gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/tempsplitdFiles/GeneLevelPLINK/
    for CHR in {1..22}
  do
  echo $CHR
  plink --bfile /gdata02/master_data1/UK_Biobank/ukb23155_b0_v1_s200604/ukb23155_c${CHR}_b0_v1 \
  --extract chr${CHR}.marker \
  --keep /gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/tempsplitdFiles/Random_50e3_Subjects.sample \
  --maf 0.0002 \
  --make-bed --out chr${CHR}_10K_Genes
  done
}

###
mergeMat = matrix(paste0("chr", 1:22,"_10K_Genes",
                         rep(c(".bed",".bim",".fam"), each=22)),ncol=3)

write.table(mergeMat, 
            "/gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/tempsplitdFiles/GeneLevelPLINK/merge_10K_Genes_50K_Subj.txt", 
            row.names = F, col.names = F, quote = F)

# plink command
if(FALSE){
  cd /gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/tempsplitdFiles/GeneLevelPLINK/
    plink --merge-list merge_10K_Genes_50K_Subj.txt --make-bed --out merge_10K_Genes_50K_Subj
}

### remaining 157505 SNPs for rare variants.
PlinkPrefix = "/gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/tempsplitdFiles/GeneLevelPLINK/merge_10K_Genes_50K_Subj"
bimFile = paste0(PlinkPrefix, ".bim")
bimData = data.table::fread(bimFile)

set.seed(123)
posInBIM = sort(sample(nrow(bimData), 157500))

SNPs = bimData$V2[posInBIM]

for(group in 1:100)
{
  TempSNPs = SNPs[1:1575+(group-1)*1575]
  FileSplit = paste0("/gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/tempsplitdFiles/rareSNPs_group_",group,".txt")
  data.table::fwrite(data.table::data.table(TempSNPs), FileSplit,
                     row.names = F, col.names = F, quote = F, sep = "\t")
}

### (A). 25,000 unrelated subjects and 25,000 (6,250×4) related subjects
# next for rare variants
library(tidyr)
library(dplyr)
library(GRAB)

PlinkPrefix = "/gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/tempsplitdFiles/GeneLevelPLINK/merge_10K_Genes_50K_Subj"
DirSplit = "/gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/tempsplitdFiles/"

bimFile = paste0(PlinkPrefix, ".bim")
famFile = paste0(PlinkPrefix, ".fam")
bedFile = paste0(PlinkPrefix, ".bed")

bimData = data.table::fread(bimFile)
famData = data.table::fread(famFile)

nSub = 25e3
nFam = 6250

# please run the code below on the slurm.
# cd "/gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/"
# sbatch -J simugeno --mem=5G -t 1-0:0 --array=1-100 -o log/%A_%a.log --wrap='Rscript /gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/SimulatedGenotypeUsingWES-2023-02-16XH.R $SLURM_ARRAY_TASK_ID'
args = commandArgs(TRUE)
print(args)
print(sessionInfo())
group = as.numeric(args)

IDsToIncludeFile = paste0(DirSplit, "rareSNPs_group_", group, ".txt");
GenoList = GRAB.SimuGMatFromGenoFile(nFam = nFam,
                                     nSub = nSub,
                                     FamMode = "4-members",
                                     GenoFile = bedFile,
                                     control = list(IDsToIncludeFile = IDsToIncludeFile))

GenoList$GenoMat[is.na(GenoList$GenoMat)] = -9

# update FID and IID
rownames(GenoList$GenoMat) = gsub("Subj-", "U-", rownames(GenoList$GenoMat))
rownames(GenoList$GenoMat) = gsub("_", "-", rownames(GenoList$GenoMat))
rownames(GenoList$GenoMat) = gsub("f", "F4-", rownames(GenoList$GenoMat))

OutputPrefix = paste0(DirSplit, "nSub_25000_nFam_6250r/", "SNPs_group_", group)

GRAB.makePlink(GenoList$GenoMat, OutputPrefix)

cmd1 = paste0("plink --file ", OutputPrefix, " --make-bed --out ", OutputPrefix);
cmd2 = paste0("rm ", OutputPrefix, ".map ", OutputPrefix, ".ped");

system(cmd1)
system(cmd2)

bimFileGroup = paste0(OutputPrefix, ".bim")
markersInGroup = paste0(DirSplit, "rareSNPs_group_", group, ".txt")
markersInGroup = data.table::fread(markersInGroup, header = F)

pos = match(markersInGroup$V1, bimData$V2)
bimDataInGroup = bimData[pos,]
data.table::fwrite(bimDataInGroup, bimFileGroup, 
                   sep = "\t", row.names = F, col.names = F, quote = F)

# plink command.
if(FALSE)
{
  cd "/gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/tempsplitdFiles/nSub_25000_nFam_6250r/"
  plink --merge-list ../group.mergelist --out nSub_25000_nFam_6250 --make-bed
  plink --bfile nSub_25000_nFam_6250 --freq --out nSub_25000_nFam_6250
  plink --bfile nSub_25000_nFam_6250 --missing --out nSub_25000_nFam_6250
}

# QC for missing rates <= 0.05; MAF >= 0.0002 $ MAF < 0.05.
library(dplyr)
DirSplit = "/gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/tempsplitdFiles/"
missing = data.table::fread(paste0(DirSplit, "nSub_25000_nFam_6250r/nSub_25000_nFam_6250.lmiss"))
freq = data.table::fread(paste0(DirSplit, "nSub_25000_nFam_6250r/nSub_25000_nFam_6250.frq"))

allinfo = inner_join(missing, freq, by = c("CHR", "SNP"))

filterinfo = allinfo %>% filter(MAF >= 0.0002 & MAF < 0.05 & F_MISS < 0.05) # ramaining 104758 SNPs.

SNPIDs = sample(filterinfo$SNP, 1e5)
data.table::fwrite(data.table::data.table(SNPIDs),
                   "/gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/nSub_25000_nFam_6250rare/SNPIDs.txt",
                   row.names = F, col.names = F, quote = F, sep = "\t")

for(i in 1:10){
  SNPIDi = SNPIDs[1:10000+10000*(i-1)]
  data.table::fwrite(data.table::data.table(SNPIDi),
                     paste0("/gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/nSub_25000_nFam_6250rare/SNPs_group", i, ".txt"),
                     row.names = F, col.names = F, quote = F, sep = "\t")
}

# plink command.
if(FALSE)
{
  plink --bfile /gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/tempsplitdFiles/nSub_25000_nFam_6250r/nSub_25000_nFam_6250 \
  --extract /gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/nSub_25000_nFam_6250rare/SNPIDs.txt \
  --make-bed --out /gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/nSub_25000_nFam_6250rare/nSub_25000_nFam_6250
}


### (B). 25,000 unrelated subjects and 25,000 (2,500×4) related subjects
# first for common variants
library(tidyr)
library(dplyr)
library(GRAB)

PlinkPrefix = "/gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/tempsplitdFiles/GeneLevelPLINK/merge_10K_Genes_50K_Subj"
DirSplit = "/gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/tempsplitdFiles/"

bimFile = paste0(PlinkPrefix, ".bim")
famFile = paste0(PlinkPrefix, ".fam")
bedFile = paste0(PlinkPrefix, ".bed")

bimData = data.table::fread(bimFile)
famData = data.table::fread(famFile)

nSub = 25e3
nFam = 2500

# please run the code below on the slurm.
# cd "/gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/"
# sbatch -J simugeno --mem=5G -t 1-0:0 --array=1-100 -o log/%A_%a.log --wrap='Rscript /gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/SimulatedGenotypeUsingWES-2023-02-16XH.R $SLURM_ARRAY_TASK_ID'
args = commandArgs(TRUE)
print(args)
print(sessionInfo())
group = as.numeric(args)

IDsToIncludeFile = paste0(DirSplit, "rareSNPs_group_", group, ".txt");
GenoList = GRAB.SimuGMatFromGenoFile(nFam = nFam,
                                     nSub = nSub,
                                     FamMode = "10-members",
                                     GenoFile = bedFile,
                                     control = list(IDsToIncludeFile = IDsToIncludeFile))

GenoList$GenoMat[is.na(GenoList$GenoMat)] = -9

# update FID and IID
rownames(GenoList$GenoMat) = gsub("Subj-", "U-", rownames(GenoList$GenoMat))
rownames(GenoList$GenoMat) = gsub("_", "-", rownames(GenoList$GenoMat))
rownames(GenoList$GenoMat) = gsub("f", "F10-", rownames(GenoList$GenoMat))

OutputPrefix = paste0(DirSplit, "nSub_25000_nFam_2500r/", "SNPs_group_", group)

GRAB.makePlink(GenoList$GenoMat, OutputPrefix)

cmd1 = paste0("plink --file ", OutputPrefix, " --make-bed --out ", OutputPrefix);
cmd2 = paste0("rm ", OutputPrefix, ".map ", OutputPrefix, ".ped");

system(cmd1)
system(cmd2)

bimFileGroup = paste0(OutputPrefix, ".bim")
markersInGroup = paste0(DirSplit, "rareSNPs_group_", group, ".txt")
markersInGroup = data.table::fread(markersInGroup, header = F)

pos = match(markersInGroup$V1, bimData$V2)
bimDataInGroup = bimData[pos,]
data.table::fwrite(bimDataInGroup, bimFileGroup, 
                   sep = "\t", row.names = F, col.names = F, quote = F)

# plink command.
if(FALSE)
{
  cd "/gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/tempsplitdFiles/nSub_25000_nFam_2500r/"
  plink --merge-list ../group.mergelist --out nSub_25000_nFam_2500 --make-bed
  plink --bfile nSub_25000_nFam_2500 --freq --out nSub_25000_nFam_2500
  plink --bfile nSub_25000_nFam_2500 --missing --out nSub_25000_nFam_2500
}

# QC for missing rates <= 0.05; MAF >= 0.0002 $ MAF < 0.05.
library(dplyr)
DirSplit = "/gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/tempsplitdFiles/"
missing = data.table::fread(paste0(DirSplit, "nSub_25000_nFam_2500r/nSub_25000_nFam_2500.lmiss"))
freq = data.table::fread(paste0(DirSplit, "nSub_25000_nFam_2500r/nSub_25000_nFam_2500.frq"))

allinfo = inner_join(missing, freq, by = c("CHR", "SNP"))

filterinfo = allinfo %>% filter(MAF >= 0.0002 & MAF < 0.05 & F_MISS < 0.05) # remaining 102878 SNPs.

SNPIDs = sample(filterinfo$SNP, 1e5)
data.table::fwrite(data.table::data.table(SNPIDs),
                   "/gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/nSub_25000_nFam_2500rare/SNPIDs.txt",
                   row.names = F, col.names = F, quote = F, sep = "\t")

for(i in 1:10){
  SNPIDi = SNPIDs[1:10000+10000*(i-1)]
  data.table::fwrite(data.table::data.table(SNPIDi),
                     paste0("/gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/nSub_25000_nFam_2500rare/SNPs_group", i, ".txt"),
                     row.names = F, col.names = F, quote = F, sep = "\t")
}

# plink command.
if(FALSE)
{
  plink --bfile /gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/tempsplitdFiles/nSub_25000_nFam_2500r/nSub_25000_nFam_2500 \
  --extract /gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/nSub_25000_nFam_2500rare/SNPIDs.txt \
  --allow-no-vars --make-bed --out /gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/nSub_25000_nFam_2500rare/nSub_25000_nFam_2500
}


### (C). 50,000 unrelated subjects
# first for common variants
library(tidyr)
library(dplyr)
library(GRAB)

PlinkPrefix = "/gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/tempsplitdFiles/GeneLevelPLINK/merge_10K_Genes_50K_Subj"
DirSplit = "/gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/tempsplitdFiles/"

bimFile = paste0(PlinkPrefix, ".bim")
famFile = paste0(PlinkPrefix, ".fam")
bedFile = paste0(PlinkPrefix, ".bed")

bimData = data.table::fread(bimFile)
famData = data.table::fread(famFile)

nSub = 50e3
nFam = 0

# please run the code below on the slurm.
# cd "/gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/"
# sbatch -J simugeno --mem=5G -t 1-0:0 --array=1-100 -o log/%A_%a.log --wrap='Rscript /gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/SimulatedGenotypeUsingWES-2023-02-16XH.R $SLURM_ARRAY_TASK_ID'
args = commandArgs(TRUE)
print(args)
print(sessionInfo())
group = as.numeric(args)

IDsToIncludeFile = paste0(DirSplit, "rareSNPs_group_", group, ".txt");
GenoList = GRAB.SimuGMatFromGenoFile(nFam = nFam,
                                     nSub = nSub,
                                     FamMode = "4-members",
                                     GenoFile = bedFile,
                                     control = list(IDsToIncludeFile = IDsToIncludeFile))

GenoList$GenoMat[is.na(GenoList$GenoMat)] = -9

# update FID and IID
rownames(GenoList$GenoMat) = gsub("Subj-", "U-", rownames(GenoList$GenoMat))
rownames(GenoList$GenoMat) = gsub("_", "-", rownames(GenoList$GenoMat))

OutputPrefix = paste0(DirSplit, "nSub_50000r/", "SNPs_group_", group)

GRAB.makePlink(GenoList$GenoMat, OutputPrefix)

cmd1 = paste0("plink --file ", OutputPrefix, " --make-bed --out ", OutputPrefix);
cmd2 = paste0("rm ", OutputPrefix, ".map ", OutputPrefix, ".ped");

system(cmd1)
system(cmd2)

bimFileGroup = paste0(OutputPrefix, ".bim")
markersInGroup = paste0(DirSplit, "rareSNPs_group_", group, ".txt")
markersInGroup = data.table::fread(markersInGroup, header = F)

pos = match(markersInGroup$V1, bimData$V2)
bimDataInGroup = bimData[pos,]
data.table::fwrite(bimDataInGroup, bimFileGroup, 
                   sep = "\t", row.names = F, col.names = F, quote = F)

# plink command.
if(FALSE)
{
  cd "/gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/tempsplitdFiles/nSub_50000r/"
  plink --merge-list ../group.mergelist --out nSub_50000 --make-bed
  plink --bfile nSub_50000 --freq --out nSub_50000
  plink --bfile nSub_50000 --missing --out nSub_50000
}

# QC for missing rates <= 0.05; MAF >= 0.0002 $ MAF < 0.05.
library(dplyr)
DirSplit = "/gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/tempsplitdFiles/"
missing = data.table::fread(paste0(DirSplit, "nSub_50000r/nSub_50000.lmiss"))
freq = data.table::fread(paste0(DirSplit, "nSub_50000r/nSub_50000.frq"))

allinfo = inner_join(missing, freq, by = c("CHR", "SNP"))

filterinfo = allinfo %>% filter(MAF >= 0.0002 & MAF < 0.05 & F_MISS < 0.05) # remaining 109655 SNPs.

SNPIDs = sample(filterinfo$SNP, 1e5)
data.table::fwrite(data.table::data.table(SNPIDs),
                   "/gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/nSub_50000rare/SNPIDs.txt",
                   row.names = F, col.names = F, quote = F, sep = "\t")

for(i in 1:10){
  SNPIDi = SNPIDs[1:10000+10000*(i-1)]
  data.table::fwrite(data.table::data.table(SNPIDi),
                     paste0("/gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/nSub_50000rare/SNPs_group", i, ".txt"),
                     row.names = F, col.names = F, quote = F, sep = "\t")
}

# plink command.
if(FALSE)
{
  plink --bfile /gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/tempsplitdFiles/nSub_50000r/nSub_50000 \
  --extract /gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/nSub_50000rare/SNPIDs.txt \
  --allow-no-vars --make-bed --out /gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/nSub_50000rare/nSub_50000
}