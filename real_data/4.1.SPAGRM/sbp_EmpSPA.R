
###### part1 preconditioning ######
# cd "/gdata01/user/xuhe/family_relatedness/real_data-2022-12-29/"
# sbatch -J sbpstep1 --mem=100G -t 1-0:0 --array=1 -o log/%A_%a.log --wrap='Rscript /gdata01/user/xuhe/family_relatedness/real_data-2022-12-29/4.EmpSPA/sbp_EmpSPA.R $SLURM_ARRAY_TASK_ID'
# library(GRAB)
# source("/gdata01/user/xuhe/family_relatedness/simulation-2022-12-04/code/subfunc-2022-12-05XH.R")
# 
# Resid = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/null_fitter/sbp_residuals.txt")
# ResidMat.beta = Resid %>% select(SubjID, R_beta) %>% rename("Resid" = "R_beta")
# ResidMat.tau = Resid %>% select(SubjID, R_tau) %>% rename("Resid" = "R_tau")
# 
# GenoMat = c()
# for(chr in 1:22)
# {
#   bgenFile = paste0("/gdata02/master_data1/UK_Biobank/ukb22828_imp/ukb22828_c", chr, "_b0_v3.bgen")
# 
#   bgenbgiFile = paste0("/gdata02/master_data1/UK_Biobank/ukb22828_imp/ukb_imp_chr", chr, "_v3.bgen.bgi")
#   bgensampleFile = "/gdata02/master_data1/UK_Biobank/ukb22828_imp/ukb22828_c1_b0_v3_s487203.sample"
# 
#   SampleIDs = ResidMat.beta$SubjID; # SampleIDs = ResidMat.tau$SubjID
# 
#   IDsToIncludeFile = paste0("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/split_geno/chr", chr, ".txt")
# 
#   GenoList = GRAB.ReadGeno(GenoFile = bgenFile,
#                            GenoFileIndex = c(bgenbgiFile, bgensampleFile),
#                            SampleIDs = SampleIDs,
#                            control = list(IDsToIncludeFile = IDsToIncludeFile,
#                                           ImputeMethod = "mean"))
#   GenoMat = cbind(GenoMat, GenoList$GenoMat[, !duplicated(colnames(GenoList$GenoMat))])
#   GenoList = NULL
# }
# 
# GenoMat = GenoMat[, sample(1:ncol(GenoMat), 20e3)]
# 
# GRM.file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/split_geno/UKB500K_sparseGRM_0.05.txt"
# 
# obj.precond.beta = precond(GRM.file,
#                            GenoMat,
#                            ResidMat = ResidMat.beta,
#                            MaxQuantile = 0.75,
#                            MinQuantile = 0.25,
#                            OutlierRatio = 1.5,
#                            output.noGRM  = TRUE)
# save(obj.precond.beta,
#      file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preconding/obj.precond.beta.sbp.RData")
# 
# obj.precond.tau = precond(GRM.file,
#                           GenoMat,
#                           ResidMat = ResidMat.tau,
#                           MaxQuantile = 0.75,
#                           MinQuantile = 0.25,
#                           OutlierRatio = 1.5,
#                           output.noGRM  = TRUE)
# save(obj.precond.tau,
#      file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preconding/obj.precond.tau.sbp.RData")


###### part2 analyzing ######
# cd "/gdata01/user/xuhe/family_relatedness/real_data-2022-12-29/"
# sbatch -J sbpstep2 --mem=5G -t 1-0:0 --array=1-2336 -o log/%A_%a.log --wrap='Rscript /gdata01/user/xuhe/family_relatedness/real_data-2022-12-29/4.EmpSPA/sbp_EmpSPA.R $SLURM_ARRAY_TASK_ID'
# args=commandArgs(TRUE)
# print(args)
# print(sessionInfo())
# n.cpu=as.numeric(args)
# print(n.cpu)
# 
# 
# library(GRAB)
# source("/gdata01/user/xuhe/family_relatedness/simulation-2022-12-04/code/subfunc-2022-12-05XH.R")
# 
# obj.precond.beta = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preconding/obj.precond.beta.sbp.RData"
# obj.precond.tau = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preconding/obj.precond.tau.sbp.RData"
# load(obj.precond.beta); load(obj.precond.tau)
# 
# load("/gdata01/user/home/wenjianb/LA_PD/SPAGE_2022-06-23/version-1/before-step2-table.RData")
# table = table %>% filter(group == n.cpu)
# 
# output.beta = c(); output.tau = c()
# for(i in 1:nrow(table))
# {
#   filename = stringr::str_split(table$file[i], pattern = "/", simplify = TRUE)[6]
#   chr = stringr::str_split(filename, pattern = "_", n = 3, simplify = TRUE)[2]
#   
#   bgenFile = paste0("/gdata02/master_data1/UK_Biobank/ukb22828_imp/ukb22828_c", chr, "_b0_v3.bgen")
#   bgenbgiFile = paste0("/gdata02/master_data1/UK_Biobank/ukb22828_imp/ukb_imp_chr", chr, "_v3.bgen.bgi")
#   bgensampleFile = "/gdata02/master_data1/UK_Biobank/ukb22828_imp/ukb22828_c1_b0_v3_s487203.sample"
#   
#   Resid.infomation = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/null_fitter/sbp_residuals.txt")
#   SampleIDs = Resid.infomation$SubjID
#   
#   IDsToIncludeFile = paste0("/gdata02/master_data1/UK_Biobank/ukb22828_imp/split_R2_0.6_MAF_5e-04/", filename)
#   
#   GenoList = GRAB.ReadGeno(GenoFile = bgenFile,
#                            GenoFileIndex = c(bgenbgiFile, bgensampleFile),
#                            SampleIDs = SampleIDs,
#                            control = list(IDsToIncludeFile = IDsToIncludeFile,
#                                           ImputeMethod = "mean"))
#   temp.output.beta = SPA_G(obj.precond.beta,
#                            GenoList$GenoMat,
#                            Cutoff = 2,
#                            impute.method = "fixed",
#                            missing.cutoff = 0.15,
#                            min.maf = 0.002,
#                            G.model = "Add",
#                            output.noGRM = F)
#   output.beta = rbind(output.beta,
#                       temp.output.beta)
#   
#   temp.output.tau = SPA_G(obj.precond.tau,
#                           GenoList$GenoMat,
#                           Cutoff = 2,
#                           impute.method = "fixed",
#                           missing.cutoff = 0.15,
#                           min.maf = 0.002,
#                           G.model = "Add",
#                           output.noGRM = F)
#   output.tau = rbind(output.tau,
#                      temp.output.tau)
# }
# 
# data.table::fwrite(output.beta %>% as.data.frame() %>% drop_na(Score),
#                    file = paste0("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/temp_pvalues/sbp.pvalue.beta-", n.cpu, ".txt"),
#                    row.names = F, col.names = T, quote = F, sep = "\t")
# 
# data.table::fwrite(output.tau %>% as.data.frame() %>% drop_na(Score),
#                    file = paste0("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/temp_pvalues/sbp.pvalue.tau-", n.cpu, ".txt"),
#                    row.names = F, col.names = T, quote = F, sep = "\t")


###### part3 merging and ploting ######
library(dplyr)
library(tidyr)

empspa.grm = c(); adj.empspa.grm = c()
for(n.rep in 1:2336)
{
  if(n.rep %% 100 == 0) cat(paste0("Have read ", n.rep, " files.\n"))
  
  betapval = data.table::fread(paste0("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/temp_pvalues/sbp.pvalue.beta-", n.rep, ".txt"))
  taupval = data.table::fread(paste0("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/temp_pvalues/sbp.pvalue.tau-", n.rep, ".txt"))
  
  if(all(betapval$SNPID == taupval$SNPID)){
    empspa.grm = rbind(empspa.grm,
                       data.frame(SNPID = betapval$SNPID,
                                  MAF = betapval$MAF,
                                  betascore = betapval$Score,
                                  tauscore = taupval$Score,
                                  betaspapval = betapval$p.value.spa.G.GRM,
                                  betanormpval = betapval$p.value.norm.GRM,
                                  tauspapval = taupval$p.value.spa.G.GRM,
                                  taunormpval = taupval$p.value.norm.GRM))
  }else{
    cat(paste0(n.rep, "th betapval and taupval file's SNPID mismatch.\n"))
  }
  
  adjbetapval = data.table::fread(paste0("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/temp_pvalues/adj.sbp.pvalue.beta-", n.rep, ".txt"))
  adjtaupval = data.table::fread(paste0("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/temp_pvalues/adj.sbp.pvalue.tau-", n.rep, ".txt"))
  
  if(all(adjbetapval$SNPID == adjtaupval$SNPID)){
    adj.empspa.grm = rbind(adj.empspa.grm,
                           data.frame(SNPID = adjbetapval$SNPID,
                                      MAF = adjbetapval$MAF,
                                      betascore = adjbetapval$Score,
                                      tauscore = adjtaupval$Score,
                                      betaspapval = adjbetapval$p.value.spa.G.GRM,
                                      betanormpval = adjbetapval$p.value.norm.GRM,
                                      tauspapval = adjtaupval$p.value.spa.G.GRM,
                                      taunormpval = adjtaupval$p.value.norm.GRM))
  }else{
    cat(paste0(n.rep, "th adjbetapval and adjtaupval file's SNPID mismatch.\n"))
  }
}

data.table::fwrite(empspa.grm,
                   "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/output/sbp_SPAGRM.txt",
                   row.names = F, col.names = T, quote = F, sep = "\t")
data.table::fwrite(adj.empspa.grm %>% filter(MAF >= 0.002),
                   "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/output/sbp_adjSPAGRM.txt",
                   row.names = F, col.names = T, quote = F, sep = "\t")

trajgwaspval = c()
for(chr in 1:22)
{
  temppval = data.table::fread(paste0("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/temp_pvalues/trajgwas_sbp_chr", chr, ".txt"))
  
  trajgwaspval = rbind(trajgwaspval,
                       temppval)
}

data.table::fwrite(trajgwaspval,
                   "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/output/sbp_trajgwas.txt",
                   row.names = F, col.names = T, quote = F, sep = "\t")

###
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)

trajgwas1 = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/output/sbp_trajgwas.txt")
trajgwas2 = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/TrajGWAS/ukbbgen_sbp_clean.csv")

all(trajgwas1$snpid == trajgwas2$snpid)

empspagrm = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/output/sbp_SPAGRM.txt")
adj.empspagrm = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/output/sbp_adjSPAGRM.txt")

all(empspagrm$SNPID == adj.empspagrm$SNPID)

snpintrajgwas = trajgwas1$snpid[!do::duplicated_all(trajgwas1$snpid)]
snpinempspagrm = empspagrm$SNPID[!do::duplicated_all(empspagrm$SNPID)]
snpinall = intersect(snpintrajgwas, snpinempspagrm)

test.trajgwas1 = trajgwas1 %>% filter(snpid %in% snpinall); test.trajgwas1 = test.trajgwas1[match(snpinall, test.trajgwas1$snpid)]
test.trajgwas2 = trajgwas2 %>% filter(snpid %in% snpinall); test.trajgwas2 = test.trajgwas2[match(snpinall, test.trajgwas2$snpid)]
test.empspagrm = empspagrm %>% filter(SNPID %in% snpinall); test.empspagrm = test.empspagrm[match(snpinall, test.empspagrm$SNPID)] 
test.adj.empspagrm = adj.empspagrm %>% filter(SNPID %in% snpinall); test.adj.empspagrm = test.adj.empspagrm[match(snpinall, test.adj.empspagrm$SNPID)]

all(test.trajgwas1$snpid == test.empspagrm$SNPID)

outputpval1 = rbind(data.frame(snpid = test.trajgwas2$snpid, 
                               pubbetapval = test.trajgwas2$betapval,
                               pubtaupval = test.trajgwas2$taupval,
                               betapval = test.trajgwas1$betapval,
                               taupval = test.trajgwas1$taupval,
                               method = "reprotraj"),
                    data.frame(snpid = test.empspagrm$SNPID, 
                               pubbetapval = test.trajgwas2$betapval,
                               pubtaupval = test.trajgwas2$taupval,
                               betapval = test.empspagrm$betaspapval,
                               taupval = test.empspagrm$tauspapval,
                               method = "empspa"),
                    data.frame(snpid = test.adj.empspagrm$SNPID, 
                               pubbetapval = test.trajgwas2$betapval,
                               pubtaupval = test.trajgwas2$taupval,
                               betapval = test.adj.empspagrm$betaspapval,
                               taupval = test.adj.empspagrm$tauspapval,
                               method = "adjspa"))

p.beta = ggplot(outputpval1, aes(x = -log10(pubbetapval), y = -log10(betapval), color = method, shape = method, group = method)) +
  geom_point(size = 1, alpha=0.8) + geom_abline(slope = 1, intercept = 0) + theme_bw() +
  labs(xlab = "-log10(Public TrajGWAS Pvalues for beta)", ylab = "-log10(Experimental Pvalues for beta)") +
  theme(legend.position = "bottom") + theme(legend.title=element_blank()) + theme(plot.title=element_text(hjust = 0.5))
p.tau = ggplot(outputpval1, aes(x = -log10(pubtaupval), y = -log10(taupval), color = method, shape = method, group = method)) +
  geom_point(size = 1, alpha=0.8) + geom_abline(slope = 1, intercept = 0) + theme_bw() +
  labs(xlab = "-log10(Public TrajGWAS Pvalues for tau)", ylab = "-log10(Experimental Pvalues for tau)") +
  theme(legend.position = "bottom") + theme(legend.title=element_blank()) + theme(plot.title=element_text(hjust = 0.5))

p1 = plot_grid(p.beta, p.tau, ncol = 2, align = "h")

ggsave("QQplot1_sbp-2023-02-09.png",
       plot = p1, width = 10.5, height = 6.5,
       path = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/output/")

outputpval2 = rbind(data.frame(snpid = test.trajgwas1$snpid, 
                               betapval = test.trajgwas1$betapval,
                               taupval = test.trajgwas1$taupval,
                               method = "reprotraj"),
                    data.frame(snpid = test.trajgwas2$snpid, 
                               betapval = test.trajgwas2$betapval,
                               taupval = test.trajgwas2$taupval,
                               method = "pubtraj"),
                    data.frame(snpid = test.empspagrm$SNPID, 
                               betapval = test.empspagrm$betaspapval,
                               taupval = test.empspagrm$tauspapval,
                               method = "empspa"),
                    data.frame(snpid = test.adj.empspagrm$SNPID, 
                               betapval = test.adj.empspagrm$betaspapval,
                               taupval = test.adj.empspagrm$tauspapval,
                               method = "adjspa"),
                    data.frame(snpid = test.empspagrm$SNPID, 
                               betapval = test.empspagrm$betanormpval,
                               taupval = test.empspagrm$taunormpval,
                               method = "empnorm"),
                    data.frame(snpid = test.adj.empspagrm$SNPID, 
                               betapval = test.adj.empspagrm$betanormpval,
                               taupval = test.adj.empspagrm$taunormpval,
                               method = "adjnorm"))

outputpval2 = outputpval2 %>% 
  group_by(method) %>% 
  mutate(estbetapval = rank(betapval)/(length(betapval) + 1)) %>%
  mutate(esttaupval = rank(taupval)/(length(taupval) + 1))

p2 = ggplot(outputpval2, aes(x = -log10(estbetapval), y = -log10(betapval), color = method, shape = method, group = method)) + 
  geom_point(alpha=0.8) + geom_abline(slope = 1, intercept = 0) + theme_bw() +
  xlab("-log10(Estimated pvalues for beta)") + ylab("-log10(Empirical pvalues for beta)") +
  theme(legend.position = "bottom") + theme(legend.title=element_blank()) + theme(plot.title=element_text(hjust = 0.5))

p3 = ggplot(outputpval2, aes(x = -log10(esttaupval), y = -log10(taupval), color = method, shape = method, group = method)) + 
  geom_point(alpha=0.8) + geom_abline(slope = 1, intercept = 0) + theme_bw() +
  xlab("-log10(Estimated pvalues for tau)") + ylab("-log10(Empirical pvalues for tau)") +
  theme(legend.position = "bottom") + theme(legend.title=element_blank()) + theme(plot.title=element_text(hjust = 0.5))

p4 = plot_grid(p2, p3, ncol = 2, align = "h")

ggsave("QQplot2_sbp-2023-02-09.png",
       plot = p4, width = 10.5, height = 6.5,
       path = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/output/")
