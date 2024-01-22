
###### pre-clean blood protein phenotype gp_clinical table.
library(dplyr)
library(tidyr)

### read in the Hct, Hb, MCH, MCHC, MCV, RDW, ESR phecode.
source("/gdata01/user/xuhe/family_relatedness/real_data-2022-12-29/1.extract_pheno/phecode-2023-05-10XH.R")

### load the cleaned gp_clinical table.
gp_clinical = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/cleaned_gp_clinical.txt")

### extract Haematocrit from gp_clinical data.
Hct = gp_clinical %>%
  filter(grepl(Hct_codes, code)) %>% 
  mutate(value1_numeric = as.numeric(value1), value2_numeric = as.numeric(value2), value3_numeric = as.numeric(value3)) %>%
  mutate(hct = coalesce(value1_numeric, value2_numeric, value3_numeric)) %>%
  filter((hct > 20 & hct < 70) | (hct > 0.2 & hct < 0.7)) %>%
  mutate(hct = ifelse(hct < 1, 100 * hct, hct)) %>%
  mutate(value3 = toupper(value3)) %>%
  group_by(eid) %>%
  mutate(mi = n()) %>%
  mutate(mean = mean(hct)) %>%
  mutate(sd = ifelse(mi == 1, 0, 
                     ifelse(sd(hct) > 4, 4, sd(hct)))) %>%
  mutate(lower_tail = ifelse(mi == 1, mean, mean - 4 * sd),
         higher_tail = ifelse(mi == 1, mean, mean + 4 * sd)) %>%
  filter(hct >= lower_tail & hct <= higher_tail) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  mutate(hct_new = ifelse(rep == 1, hct, 
                          ifelse(max(hct) - min(hct) <= 1, 
                                 mean(hct), hct[which.min(abs(hct - mean))]))) %>%
  mutate(hct_new = round(hct_new, 0)) %>%
  select(eid, data_provider, event_dt, hct_new) %>% 
  dplyr::rename(hct = hct_new) %>%
  distinct(eid, event_dt, .keep_all = TRUE)

### write the pre-cleaned Haematocrit into file.
data.table::fwrite(Hct %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/Hct.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### extract Haemoglobin from gp_clinical data.
Hb = gp_clinical %>%
  filter(grepl(Hb_codes, code)) %>% 
  mutate(value1_numeric = as.numeric(value1), value2_numeric = as.numeric(value2), value3_numeric = as.numeric(value3)) %>%
  mutate(hb = coalesce(value1_numeric, value2_numeric, value3_numeric)) %>%
  filter((hb > 5 & hb < 25) | (hb > 50 & hb < 250)) %>%
  mutate(hb = ifelse(hb > 50, 0.1 * hb, hb)) %>%
  mutate(value3 = toupper(value3)) %>%
  filter(value3 != "MEA001") %>%
  group_by(eid) %>%
  mutate(mi = n()) %>%
  mutate(mean = mean(hb)) %>%
  mutate(sd = ifelse(mi == 1, 0, 
                     ifelse(sd(hb) > 1.5, 1.5, sd(hb)))) %>%
  mutate(lower_tail = ifelse(mi == 1, mean, mean - 4 * sd),
         higher_tail = ifelse(mi == 1, mean, mean + 4 * sd)) %>%
  filter(hb >= lower_tail & hb <= higher_tail) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  mutate(hb_new = ifelse(rep == 1, hb, 
                         ifelse(max(hb) - min(hb) <= 1, 
                                mean(hb), hb[which.min(abs(hb - mean))]))) %>%
  mutate(hb_new = round(hb_new, 1)) %>%
  select(eid, data_provider, event_dt, hb_new) %>% 
  dplyr::rename(hb = hb_new) %>%
  distinct(eid, event_dt, .keep_all = TRUE)

### write the pre-cleaned Haematocrit into file.
data.table::fwrite(Hb %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/Hb.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### extract mean cell haemoglobin from gp_clinical data.
MCH = gp_clinical %>%
  filter(grepl(MCH_codes, code)) %>% 
  mutate(value1_numeric = as.numeric(value1), value2_numeric = as.numeric(value2), value3_numeric = as.numeric(value3)) %>%
  mutate(mch = coalesce(value1_numeric, value2_numeric, value3_numeric)) %>%
  filter(mch > 15 & mch < 45) %>%
  mutate(value3 = toupper(value3)) %>%
  filter(value3 %in% c("", "MEA000", "MEA114", "MEA116", "PG")) %>%
  group_by(eid) %>%
  mutate(mi = n()) %>%
  mutate(mean = mean(mch)) %>%
  mutate(sd = ifelse(mi == 1, 0, 
                     ifelse(sd(mch) > 2, 2, sd(mch)))) %>%
  mutate(lower_tail = ifelse(mi == 1, mean, mean - 4 * sd),
         higher_tail = ifelse(mi == 1, mean, mean + 4 * sd)) %>%
  filter(mch >= lower_tail & mch <= higher_tail) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  mutate(mch_new = ifelse(rep == 1, mch, 
                         ifelse(max(mch) - min(mch) <= 1, 
                                mean(mch), mch[which.min(abs(mch - mean))]))) %>%
  mutate(mch_new = round(mch_new, 1)) %>%
  select(eid, data_provider, event_dt, mch_new) %>% 
  dplyr::rename(mch = mch_new) %>%
  distinct(eid, event_dt, .keep_all = TRUE)

### write the pre-cleaned mean cell haemoglobin into file.
data.table::fwrite(MCH %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/MCH.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### extract Mean cell haemoglobin concentration from gp_clinical data.
MCHC = gp_clinical %>%
  filter(grepl(MCHC_codes, code)) %>% 
  mutate(value1_numeric = as.numeric(value1), value2_numeric = as.numeric(value2), value3_numeric = as.numeric(value3)) %>%
  mutate(mchc = coalesce(value1_numeric, value2_numeric, value3_numeric)) %>%
  filter((mchc > 2.5 & mchc < 4) | (mchc > 25 & mchc < 40) | (mchc > 250 & mchc < 400)) %>%
  mutate(mchc = ifelse(mchc < 4, 10 * mchc, mchc)) %>%
  mutate(mchc = ifelse(mchc > 250, 0.1 * mchc, mchc)) %>%
  mutate(value3 = toupper(value3)) %>%
  filter(value3 != "MEA096") %>%
  group_by(eid) %>%
  mutate(mi = n()) %>%
  mutate(mean = mean(mchc)) %>%
  mutate(sd = ifelse(mi == 1, 0, 
                     ifelse(sd(mchc) > 1.25, 1.25, sd(mchc)))) %>%
  mutate(lower_tail = ifelse(mi == 1, mean, mean - 4 * sd),
         higher_tail = ifelse(mi == 1, mean, mean + 4 * sd)) %>%
  filter(mchc >= lower_tail & mchc <= higher_tail) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  mutate(mchc_new = ifelse(rep == 1, mchc, 
                           ifelse(max(mchc) - min(mchc) <= 1, 
                                  mean(mchc), mchc[which.min(abs(mchc - mean))]))) %>%
  mutate(mchc_new = round(mchc_new, 1)) %>%
  select(eid, data_provider, event_dt, mchc_new) %>% 
  dplyr::rename(mchc = mchc_new) %>%
  distinct(eid, event_dt, .keep_all = TRUE)

### write the pre-cleaned Mean cell haemoglobin concentration into file.
data.table::fwrite(MCHC %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/MCHC.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### extract Mean corpuscular volume from gp_clinical data.
MCV = gp_clinical %>%
  filter(grepl(MCV_codes, code)) %>% 
  mutate(value1_numeric = as.numeric(value1), value2_numeric = as.numeric(value2), value3_numeric = as.numeric(value3)) %>%
  mutate(mcv = coalesce(value1_numeric, value2_numeric, value3_numeric)) %>%
  filter(mcv > 50 & mcv < 130) %>%
  mutate(value3 = toupper(value3)) %>%
  filter(value3 != "MICRONS") %>%
  group_by(eid) %>%
  mutate(mi = n()) %>%
  mutate(mean = mean(mcv)) %>%
  mutate(sd = ifelse(mi == 1, 0, 
                     ifelse(sd(mcv) > 4.5, 4.5, sd(mcv)))) %>%
  mutate(lower_tail = ifelse(mi == 1, mean, mean - 4 * sd),
         higher_tail = ifelse(mi == 1, mean, mean + 4 * sd)) %>%
  filter(mcv >= lower_tail & mcv <= higher_tail) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  mutate(mcv_new = ifelse(rep == 1, mcv, 
                           ifelse(max(mcv) - min(mcv) <= 1, 
                                  mean(mcv), mcv[which.min(abs(mcv - mean))]))) %>%
  mutate(mcv_new = round(mcv_new, 0)) %>%
  select(eid, data_provider, event_dt, mcv_new) %>% 
  dplyr::rename(mcv = mcv_new) %>%
  distinct(eid, event_dt, .keep_all = TRUE)

### write the pre-cleaned Mean corpuscular volume into file.
data.table::fwrite(MCV %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/MCV.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### extract Red blood cell distribut width from gp_clinical data.
RDW = gp_clinical %>%
  filter(grepl(RDW_codes, code)) %>% 
  mutate(value1_numeric = as.numeric(value1), value2_numeric = as.numeric(value2), value3_numeric = as.numeric(value3)) %>%
  mutate(rdw = coalesce(value1_numeric, value2_numeric, value3_numeric)) %>%
  filter(rdw > 8 & rdw < 30) %>%
  mutate(value3 = toupper(value3)) %>%
  filter(value3 %in% c("", "MEA000", "MEA001")) %>%
  group_by(eid) %>%
  mutate(mi = n()) %>%
  mutate(mean = mean(rdw)) %>%
  mutate(sd = ifelse(mi == 1, 0, 
                     ifelse(sd(rdw) > 1.5, 1.5, sd(rdw)))) %>%
  mutate(lower_tail = ifelse(mi == 1, mean, mean - 4 * sd),
         higher_tail = ifelse(mi == 1, mean, mean + 4 * sd)) %>%
  filter(rdw >= lower_tail & rdw <= higher_tail) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  mutate(rdw_new = ifelse(rep == 1, rdw, 
                          ifelse(max(rdw) - min(rdw) <= 1, 
                                 mean(rdw), rdw[which.min(abs(rdw - mean))]))) %>%
  mutate(rdw_new = round(rdw_new, 1)) %>%
  select(eid, data_provider, event_dt, rdw_new) %>% 
  dplyr::rename(rdw = rdw_new) %>%
  distinct(eid, event_dt, .keep_all = TRUE)

### write the pre-cleaned Red blood cell distribut width into file.
data.table::fwrite(RDW %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/RDW.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### extract Erythrocyte sedimentation rate from gp_clinical data.
ESR = gp_clinical %>%
  filter(grepl(ESR_codes, code)) %>% 
  mutate(value1_numeric = as.numeric(value1), value2_numeric = as.numeric(value2), value3_numeric = as.numeric(value3)) %>%
  mutate(esr = coalesce(value1_numeric, value2_numeric, value3_numeric)) %>%
  filter(esr >= 1 & esr < 100) %>%
  mutate(value3 = toupper(value3)) %>%
  filter(value3 != "MEA194") %>%
  group_by(eid) %>%
  mutate(mi = n()) %>%
  mutate(mean = mean(esr)) %>%
  mutate(sd = ifelse(mi == 1, 0, 
                     ifelse(sd(esr) > 15, 15, sd(esr)))) %>%
  mutate(lower_tail = ifelse(mi == 1, mean, mean - 4 * sd),
         higher_tail = ifelse(mi == 1, mean, mean + 4 * sd)) %>%
  filter(esr >= lower_tail & esr <= higher_tail) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  mutate(esr_new = ifelse(rep == 1, esr, 
                          ifelse(max(esr) - min(esr) <= 5, 
                                 mean(esr), esr[which.min(abs(esr - mean))]))) %>%
  mutate(esr_new = round(esr_new, 1)) %>%
  select(eid, data_provider, event_dt, esr_new) %>% 
  dplyr::rename(esr = esr_new) %>%
  distinct(eid, event_dt, .keep_all = TRUE)

### write the pre-cleaned Erythrocyte sedimentation rate into file.
data.table::fwrite(ESR %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/ESR.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)


###### merge the gp_clinical table and covariates.
library(dplyr)
library(tidyr)
library(lubridate)

### read in the gp_clinical table and covariates.
Hct_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/Hct.txt")

Hb_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/Hb.txt")

MCH_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/MCH.txt")

MCHC_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/MCHC.txt")

MCV_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/MCV.txt")

RDW_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/RDW.txt")

ESR_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/ESR.txt")

cova_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/assessment_center/covariate.txt")

bmi_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/assessment_center/bmi_from_UKB.txt")

### combine the gp_clinical and covariates data.
Hct_joined = Hct_info %>%
  inner_join(cova_info, by = c('eid' = 'Wenjian')) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18) %>%
  left_join(bmi_info %>% select(feid, f21001_0_0) %>% rename(eid = feid, BMI = f21001_0_0))

### write the final-cleaned Haematocrit into file.
data.table::fwrite(Hct_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/Hct.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### combine the gp_clinical and covariates data.
Hb_joined = Hb_info %>%
  inner_join(cova_info, by = c('eid' = 'Wenjian')) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18) %>%
  left_join(bmi_info %>% select(feid, f21001_0_0) %>% rename(eid = feid, BMI = f21001_0_0))

### write the final-cleaned Haemoglobin into file.
data.table::fwrite(Hb_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/Hb.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### combine the gp_clinical and covariates data.
MCH_joined = MCH_info %>%
  inner_join(cova_info, by = c('eid' = 'Wenjian')) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18) %>%
  left_join(bmi_info %>% select(feid, f21001_0_0) %>% rename(eid = feid, BMI = f21001_0_0))

### write the final-cleaned mean cell haemoglobin into file.
data.table::fwrite(MCH_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/MCH.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### combine the gp_clinical and covariates data.
MCHC_joined = MCHC_info %>%
  inner_join(cova_info, by = c('eid' = 'Wenjian')) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18) %>%
  left_join(bmi_info %>% select(feid, f21001_0_0) %>% rename(eid = feid, BMI = f21001_0_0))

### write the final-cleaned Mean cell haemoglobin concentration into file.
data.table::fwrite(MCHC_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/MCHC.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### combine the gp_clinical and covariates data.
MCV_joined = MCV_info %>%
  inner_join(cova_info, by = c('eid' = 'Wenjian')) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18) %>%
  left_join(bmi_info %>% select(feid, f21001_0_0) %>% rename(eid = feid, BMI = f21001_0_0))

### write the final-cleaned Mean corpuscular volume into file.
data.table::fwrite(MCV_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/MCV.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### combine the gp_clinical and covariates data.
RDW_joined = RDW_info %>%
  inner_join(cova_info, by = c('eid' = 'Wenjian')) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18) %>%
  left_join(bmi_info %>% select(feid, f21001_0_0) %>% rename(eid = feid, BMI = f21001_0_0))

### write the final-cleaned Red blood cell distribut width into file.
data.table::fwrite(RDW_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/RDW.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### combine the gp_clinical and covariates data.
ESR_joined = ESR_info %>%
  inner_join(cova_info, by = c('eid' = 'Wenjian')) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18) %>%
  left_join(bmi_info %>% select(feid, f21001_0_0) %>% rename(eid = feid, BMI = f21001_0_0))

### write the final-cleaned Erythrocyte sedimentation rate into file.
data.table::fwrite(ESR_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/ESR.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)
