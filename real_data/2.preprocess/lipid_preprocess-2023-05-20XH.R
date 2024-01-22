
library(dplyr)
library(tidyr)

### read in the HDL, LDL, NHDL, total cholesterol and Triglyceride phecode.
source("/gdata01/user/xuhe/family_relatedness/real_data-2022-12-29/1.extract_pheno/phecode-2023-05-10XH.R")

### load the cleaned gp_clinical table
gp_clinical = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/cleaned_gp_clinical.txt")

### extract HDL from gp_clinical data.
HDL = gp_clinical %>%
  filter(grepl(HDL_codes, code)) %>% 
  mutate(value1_numeric = as.numeric(value1), value2_numeric = as.numeric(value2), value3_numeric = as.numeric(value3)) %>%
  mutate(HDL = coalesce(value1_numeric, value2_numeric, value3_numeric)) %>%
  filter(HDL > 0.2 & HDL < 5) %>%
  mutate(value3 = toupper(value3)) %>%
  filter(value3 %in% c("", "MEA000", "MEA096", "MMOL/L", "UNKNOWN")) %>%
  group_by(eid) %>%
  mutate(mi = n()) %>%
  mutate(mean = mean(HDL)) %>%
  mutate(sd = ifelse(mi == 1, 0, 
                     ifelse(sd(HDL) > 0.4, 0.4, sd(HDL)))) %>%
  mutate(lower_tail = ifelse(mi == 1, mean, mean - 4 * sd),
         higher_tail = ifelse(mi == 1, mean, mean + 4 * sd)) %>%
  filter(HDL >= lower_tail & HDL <= higher_tail) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  mutate(HDL_new = ifelse(rep == 1, HDL, 
                            ifelse(max(HDL) - min(HDL) <= 0.5, 
                                   mean(HDL), HDL[which.min(abs(HDL - mean))]))) %>%
  mutate(HDL_new = round(HDL_new, 2)) %>%
  select(eid, data_provider, event_dt, HDL_new) %>% 
  dplyr::rename(HDL = HDL_new) %>%
  distinct(eid, event_dt, .keep_all = TRUE)

### write the pre-cleaned HDL into file.
data.table::fwrite(HDL %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/HDL.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### extract LDL from gp_clinical data.
LDL = gp_clinical %>%
  filter(grepl(LDL_codes, code)) %>% 
  mutate(value1_numeric = as.numeric(value1), value2_numeric = as.numeric(value2), value3_numeric = as.numeric(value3)) %>%
  mutate(LDL = coalesce(value1_numeric, value2_numeric, value3_numeric)) %>%
  filter(LDL > 0.2 & LDL < 10) %>%
  mutate(value3 = toupper(value3)) %>%
  filter(value3 %in% c("", "MEA000", "MEA096", "MMOL/L", "UNKNOWN")) %>%
  group_by(eid) %>%
  mutate(mi = n()) %>%
  mutate(mean = mean(LDL)) %>%
  mutate(sd = ifelse(mi == 1, 0, 
                     ifelse(sd(LDL) > 0.8, 0.8, sd(LDL)))) %>%
  mutate(lower_tail = ifelse(mi == 1, mean, mean - 4 * sd),
         higher_tail = ifelse(mi == 1, mean, mean + 4 * sd)) %>%
  filter(LDL >= lower_tail & LDL <= higher_tail) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  mutate(LDL_new = ifelse(rep == 1, LDL, 
                          ifelse(max(LDL) - min(LDL) <= 0.5, 
                                 mean(LDL), LDL[which.min(abs(LDL - mean))]))) %>%
  mutate(LDL_new = round(LDL_new, 2)) %>%
  select(eid, data_provider, event_dt, LDL_new) %>% 
  dplyr::rename(LDL = LDL_new) %>%
  distinct(eid, event_dt, .keep_all = TRUE)

### write the pre-cleaned LDL into file.
data.table::fwrite(LDL %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/LDL.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### extract NHDL from gp_clinical data.
NHDL = gp_clinical %>%
  filter(grepl(NHDL_codes, code)) %>% 
  mutate(value1_numeric = as.numeric(value1), value2_numeric = as.numeric(value2), value3_numeric = as.numeric(value3)) %>%
  mutate(NHDL = coalesce(value1_numeric, value2_numeric, value3_numeric)) %>%
  filter(NHDL > 0.2 & NHDL < 10) %>%
  mutate(value3 = toupper(value3)) %>%
  group_by(eid) %>%
  mutate(mi = n()) %>%
  mutate(mean = mean(NHDL)) %>%
  mutate(sd = ifelse(mi == 1, 0, 
                     ifelse(sd(NHDL) > 0.8, 0.8, sd(NHDL)))) %>%
  mutate(lower_tail = ifelse(mi == 1, mean, mean - 4 * sd),
         higher_tail = ifelse(mi == 1, mean, mean + 4 * sd)) %>%
  filter(NHDL >= lower_tail & NHDL <= higher_tail) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  mutate(NHDL_new = ifelse(rep == 1, NHDL, 
                          ifelse(max(NHDL) - min(NHDL) <= 0.5, 
                                 mean(NHDL), NHDL[which.min(abs(NHDL - mean))]))) %>%
  mutate(NHDL_new = round(NHDL_new, 2)) %>%
  select(eid, data_provider, event_dt, NHDL_new) %>% 
  dplyr::rename(NHDL = NHDL_new) %>%
  distinct(eid, event_dt, .keep_all = TRUE)

### write the pre-cleaned NHDL into file.
data.table::fwrite(NHDL %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/NHDL.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### extract total cholesterol from gp_clinical data.
TC = gp_clinical %>%
  filter(grepl(totchol_codes, code)) %>% 
  mutate(value1_numeric = as.numeric(value1), value2_numeric = as.numeric(value2), value3_numeric = as.numeric(value3)) %>%
  mutate(TC = coalesce(value1_numeric, value2_numeric, value3_numeric)) %>%
  filter(TC > 0.6 & TC < 16) %>%
  mutate(value3 = toupper(value3)) %>%
  filter(value3 %in% c("", "MEA000", "mea096", "MMOL/L", "UNKNOWN")) %>%
  group_by(eid) %>%
  mutate(mi = n()) %>%
  mutate(mean = mean(TC)) %>%
  mutate(sd = ifelse(mi == 1, 0, 
                     ifelse(sd(TC) > 1, 1, sd(TC)))) %>%
  mutate(lower_tail = ifelse(mi == 1, mean, mean - 4 * sd),
         higher_tail = ifelse(mi == 1, mean, mean + 4 * sd)) %>%
  filter(TC >= lower_tail & TC <= higher_tail) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  mutate(TC_new = ifelse(rep == 1, TC, 
                           ifelse(max(TC) - min(TC) <= 0.5, 
                                  mean(TC), TC[which.min(abs(TC - mean))]))) %>%
  mutate(TC_new = round(TC_new, 2)) %>%
  select(eid, data_provider, event_dt, TC_new) %>% 
  dplyr::rename(TC = TC_new) %>%
  distinct(eid, event_dt, .keep_all = TRUE)

### write the pre-cleaned TC into file.
data.table::fwrite(TC %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/totchol.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### extract triglyceride from gp_clinical data.
triglyc = gp_clinical %>%
  filter(grepl(triglyc_codes, code)) %>% 
  mutate(value1_numeric = as.numeric(value1), value2_numeric = as.numeric(value2), value3_numeric = as.numeric(value3)) %>%
  mutate(tri = coalesce(value1_numeric, value2_numeric, value3_numeric)) %>%
  filter(tri > 0.2 & tri < 12) %>%
  mutate(value3 = toupper(value3)) %>%
  filter(value3 %in% c("", "MEA000", "MEA096", "MMOL/L", "UNKNOWN")) %>%
  group_by(eid) %>%
  mutate(mi = n()) %>%
  mutate(mean = mean(tri)) %>%
  mutate(sd = ifelse(mi == 1, 0, 
                     ifelse(sd(tri) > 1, 1, sd(tri)))) %>%
  mutate(lower_tail = ifelse(mi == 1, mean, mean - 4 * sd),
         higher_tail = ifelse(mi == 1, mean, mean + 4 * sd)) %>%
  filter(tri >= lower_tail & tri <= higher_tail) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  mutate(tri_new = ifelse(rep == 1, tri, 
                         ifelse(max(tri) - min(tri) <= 0.5, 
                                mean(tri), tri[which.min(abs(tri - mean))]))) %>%
  mutate(tri_new = round(tri_new, 2)) %>%
  select(eid, data_provider, event_dt, tri_new) %>% 
  dplyr::rename(tri = tri_new) %>%
  distinct(eid, event_dt, .keep_all = TRUE)

### write the pre-cleaned triglyceride into file.
data.table::fwrite(triglyc %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/triglyceride.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)


###### merge the gp_clinical table and covariates.
library(dplyr)
library(tidyr)
library(lubridate)

### read in the gp_clinical table and covariates.
HDL_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/HDL.txt")

LDL_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/LDL.txt")

NHDL_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/NHDL.txt")

TC_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/totchol.txt")

triglyc_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/triglyceride.txt")

cova_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/assessment_center/covariate.txt")

bmi_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/assessment_center/bmi_from_UKB.txt")

drug_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/assessment_center/covariate_drugs.txt")

### combine the gp_clinical and covariates data.
HDL_joined = HDL_info %>%
  inner_join(cova_info, by = c('eid' = 'Wenjian')) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18) %>%
  left_join(drug_info %>% select(feid, self_cholesteroldrugs), by = c('eid' = 'feid')) %>%
  left_join(bmi_info %>% select(feid, f21001_0_0) %>% rename(eid = feid, BMI = f21001_0_0)) %>%
  mutate(HDL_shifted = ifelse(self_cholesteroldrugs == 1, HDL - 0.06, HDL))

### write the final-cleaned HDL into file.
data.table::fwrite(HDL_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/HDL.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### combine the gp_clinical and covariates data.
LDL_joined = LDL_info %>%
  inner_join(cova_info, by = c('eid' = 'Wenjian')) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18) %>%
  left_join(drug_info %>% select(feid, self_cholesteroldrugs), by = c('eid' = 'feid')) %>%
  left_join(bmi_info %>% select(feid, f21001_0_0) %>% rename(eid = feid, BMI = f21001_0_0)) %>%
  mutate(LDL_shifted = ifelse(self_cholesteroldrugs == 1, LDL + 1.29, LDL))

### write the final-cleaned LDL into file.
data.table::fwrite(LDL_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/LDL.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### combine the gp_clinical and covariates data.
NHDL_joined = NHDL_info %>%
  inner_join(cova_info, by = c('eid' = 'Wenjian')) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18) %>%
  left_join(drug_info %>% select(feid, self_cholesteroldrugs), by = c('eid' = 'feid')) %>%
  left_join(bmi_info %>% select(feid, f21001_0_0) %>% rename(eid = feid, BMI = f21001_0_0))

### write the final-cleaned NHDL into file.
data.table::fwrite(NHDL_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/NHDL.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### combine the gp_clinical and covariates data.
TC_joined = TC_info %>%
  inner_join(cova_info, by = c('eid' = 'Wenjian')) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18) %>%
  left_join(drug_info %>% select(feid, self_cholesteroldrugs), by = c('eid' = 'feid')) %>%
  left_join(bmi_info %>% select(feid, f21001_0_0) %>% rename(eid = feid, BMI = f21001_0_0)) %>%
  mutate(TC_shifted = ifelse(self_cholesteroldrugs == 1, TC + 1.347, TC))

### write the final-cleaned total cholesterol into file.
data.table::fwrite(TC_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/totchol.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### combine the gp_clinical and covariates data.
triglyc_joined = triglyc_info %>%
  inner_join(cova_info, by = c('eid' = 'Wenjian')) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18) %>%
  left_join(drug_info %>% select(feid, self_cholesteroldrugs), by = c('eid' = 'feid')) %>%
  left_join(bmi_info %>% select(feid, f21001_0_0) %>% rename(eid = feid, BMI = f21001_0_0)) %>%
  mutate(tri_shifted = ifelse(self_cholesteroldrugs == 1, tri + 0.208, tri))

### write the final-cleaned triglyceride into file.
data.table::fwrite(triglyc_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/triglyceride.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)
