
###### pre-clean the lung function phenotype gp_clinical table.
library(dplyr)
library(tidyr)

### read in the spo2, FEV1, FVC and PEFR phecode.
source("/gdata01/user/xuhe/family_relatedness/real_data-2022-12-29/1.extract_pheno/phecode-2023-05-10XH.R")

### load the cleaned gp_clinical table.
gp_clinical = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/cleaned_gp_clinical.txt")

### extract Oxygen saturation at periphery (spo2) from gp_clinical data.
spo2 = gp_clinical %>%
  filter(grepl(spo2_codes, code)) %>% 
  mutate(value1_numeric = as.numeric(value1), value2_numeric = as.numeric(value2), value3_numeric = as.numeric(value3)) %>%
  mutate(spo2 = coalesce(value1_numeric, value2_numeric, value3_numeric)) %>%
  filter(spo2 > 80 & spo2 <= 100) %>%
  mutate(value3 = toupper(value3)) %>%
  filter(value3 == "") %>%
  group_by(eid) %>%
  mutate(mi = n()) %>%
  mutate(mean = mean(spo2)) %>%
  mutate(sd = ifelse(mi == 1, 0, 
                     ifelse(sd(spo2) > 2.5, 2.5, sd(spo2)))) %>%
  mutate(lower_tail = ifelse(mi == 1, mean, mean - 4 * sd),
         higher_tail = ifelse(mi == 1, mean, mean + 4 * sd)) %>%
  filter(spo2 >= lower_tail & spo2 <= higher_tail) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  mutate(spo2_new = ifelse(rep == 1, spo2, 
                           ifelse(max(spo2) - min(spo2) <= 5, 
                                  mean(spo2), spo2[which.min(abs(spo2 - mean))]))) %>%
  mutate(spo2_new = round(spo2_new, 0)) %>%
  select(eid, data_provider, event_dt, spo2_new) %>% 
  dplyr::rename(spo2 = spo2_new) %>%
  distinct(eid, event_dt, .keep_all = TRUE)

### write the pre-cleaned Oxygen saturation at periphery (spo2) into file.
data.table::fwrite(spo2 %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/spo2.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### extract forced expiratory volume in one second (FEV1) from gp_clinical data.
FEV1 = gp_clinical %>%
  filter(grepl(FEV1_codes, code)) %>% 
  mutate(value1_numeric = as.numeric(value1), value2_numeric = as.numeric(value2), value3_numeric = as.numeric(value3)) %>%
  mutate(fev1 = coalesce(value1_numeric, value2_numeric, value3_numeric)) %>%
  filter(fev1 > 0.1 & fev1 < 6) %>%
  mutate(value3 = toupper(value3)) %>%
  filter(value3 %in% c("", "L", "LITRE", "LITRES", "MEA000", "MEA154", "UNKNOWN")) %>%
  group_by(eid) %>%
  mutate(mi = n()) %>%
  mutate(mean = mean(fev1)) %>%
  mutate(sd = ifelse(mi == 1, 0, 
                     ifelse(sd(fev1) > 0.5, 0.5, sd(fev1)))) %>%
  mutate(lower_tail = ifelse(mi == 1, mean, mean - 4 * sd),
         higher_tail = ifelse(mi == 1, mean, mean + 4 * sd)) %>%
  filter(fev1 >= lower_tail & fev1 <= higher_tail) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  mutate(fev1_new = ifelse(rep == 1, fev1, 
                           ifelse(max(fev1) - min(fev1) <= 1, 
                                  mean(fev1), fev1[which.min(abs(fev1 - mean))]))) %>%
  mutate(fev1_new = round(fev1_new, 2)) %>%
  select(eid, data_provider, event_dt, fev1_new) %>% 
  dplyr::rename(fev1 = fev1_new) %>%
  distinct(eid, event_dt, .keep_all = TRUE)

### write the pre-cleaned forced expiratory volume in one second (FEV1) into file.
data.table::fwrite(FEV1 %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/FEV1.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### extract forced vital capacity (FVC) from gp_clinical data.
FVC = gp_clinical %>%
  filter(grepl(FVC_codes, code)) %>% 
  mutate(value1_numeric = as.numeric(value1), value2_numeric = as.numeric(value2), value3_numeric = as.numeric(value3)) %>%
  mutate(fvc = coalesce(value1_numeric, value2_numeric, value3_numeric)) %>%
  filter(fvc > 0.1 & fvc < 8) %>%
  mutate(value3 = toupper(value3)) %>%
  filter(value3 %in% c("", "L", "LITRE", "LITRES", "MEA154")) %>%
  group_by(eid) %>%
  mutate(mi = n()) %>%
  mutate(mean = mean(fvc)) %>%
  mutate(sd = ifelse(mi == 1, 0, 
                     ifelse(sd(fvc) > 0.5, 0.5, sd(fvc)))) %>%
  mutate(lower_tail = ifelse(mi == 1, mean, mean - 4 * sd),
         higher_tail = ifelse(mi == 1, mean, mean + 4 * sd)) %>%
  filter(fvc >= lower_tail & fvc <= higher_tail) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  mutate(fvc_new = ifelse(rep == 1, fvc, 
                           ifelse(max(fvc) - min(fvc) <= 1, 
                                  mean(fvc), fvc[which.min(abs(fvc - mean))]))) %>%
  mutate(fvc_new = round(fvc_new, 2)) %>%
  select(eid, data_provider, event_dt, fvc_new) %>% 
  dplyr::rename(fvc = fvc_new) %>%
  distinct(eid, event_dt, .keep_all = TRUE)

### write the pre-cleaned forced vital capacity (FVC) into file.
data.table::fwrite(FVC %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/FVC.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### extract FEV1/FVC ratio from gp_clinical data.
FEV1_FVC = gp_clinical %>%
  filter(grepl(FEV1_FVC_codes, code)) %>%
  mutate(value1_numeric = as.numeric(value1), value2_numeric = as.numeric(value2), value3_numeric = as.numeric(value3)) %>%
  mutate(ratio_report = coalesce(value1_numeric, value2_numeric, value3_numeric)) %>%
  filter(ratio_report > 0) %>%
  mutate(ratio_report = ifelse(ratio_report < 1, 100*ratio_report, ratio_report)) %>%
  filter(ratio_report > 20 & ratio_report < 100) %>%
  full_join(inner_join(FEV1, FVC) %>% 
              mutate(ratio_calcu = 100 * fev1 / fvc) %>%
              filter(ratio_calcu > 20 & ratio_calcu < 100)) %>%
  mutate(ratio_both = coalesce(ratio_calcu, ratio_report)) %>%
  filter(!(!is.na(ratio_report) & abs(ratio_both - ratio_report) > 10)) %>%
  group_by(eid) %>%
  mutate(mi = n()) %>%
  mutate(mean = mean(ratio_both)) %>%
  mutate(sd = ifelse(mi == 1, 0, 
                     ifelse(sd(ratio_both) > 7.5, 7.5, sd(ratio_both)))) %>%
  mutate(lower_tail = ifelse(mi == 1, mean, mean - 4 * sd),
         higher_tail = ifelse(mi == 1, mean, mean + 4 * sd)) %>%
  filter(ratio_both >= lower_tail & ratio_both <= higher_tail) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  mutate(ratio_new = ifelse(rep == 1, ratio_both, 
                            ifelse(max(ratio_both) - min(ratio_both) <= 10, 
                                   mean(ratio_both), ratio_both[which.min(abs(ratio_both - mean))]))) %>%
  mutate(ratio_new = round(ratio_new, 0)) %>%
  select(eid, data_provider, event_dt, ratio_new) %>% 
  dplyr::rename(ratio = ratio_new) %>%
  distinct(eid, event_dt, .keep_all = TRUE)
  
### write the pre-cleaned FEV1/FVC ratio into file.
data.table::fwrite(FEV1_FVC %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/FEV1_FVC.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### extract peak expiratory flow (PEF) from gp_clinical data.
PEF = gp_clinical %>%
  filter(grepl(PEF_codes, code)) %>% 
  mutate(value1_numeric = as.numeric(value1), value2_numeric = as.numeric(value2), value3_numeric = as.numeric(value3)) %>%
  mutate(pef = coalesce(value1_numeric, value2_numeric, value3_numeric)) %>%
  filter(pef > 10 & pef < 900) %>%
  mutate(value3 = toupper(value3)) %>%
  filter(value3 %in% c("", "L/MIN", "L/MINUTE", "LITRES/MIN", "LITRES/MINUTE", 
                       "MEA000", "MEA071", "PEFR", "UNKNOWN")) %>%
  group_by(eid) %>%
  mutate(mi = n()) %>%
  mutate(mean = mean(pef)) %>%
  mutate(sd = ifelse(mi == 1, 0, 
                     ifelse(sd(pef) > 40, 40, sd(pef)))) %>%
  mutate(lower_tail = ifelse(mi == 1, mean, mean - 4 * sd),
         higher_tail = ifelse(mi == 1, mean, mean + 4 * sd)) %>%
  filter(pef >= lower_tail & pef <= higher_tail) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  mutate(pef_new = ifelse(rep == 1, pef, 
                          ifelse(max(pef) - min(pef) <= 100, 
                                 mean(pef), pef[which.min(abs(pef - mean))]))) %>%
  mutate(pef_new = round(pef_new, 0)) %>%
  select(eid, data_provider, event_dt, pef_new) %>% 
  dplyr::rename(pef = pef_new) %>%
  distinct(eid, event_dt, .keep_all = TRUE)

### write the pre-cleaned peak expiratory flow (PEF) into file.
data.table::fwrite(PEF %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/PEF.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)


###### merge the gp_clinical table and covariates.
library(dplyr)
library(tidyr)
library(lubridate)

### read in the gp_clinical table and covariates.
spo2_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/spo2.txt")

FEV1_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/FEV1.txt")

FVC_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/FVC.txt")

FEV1_FVC_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/FEV1_FVC.txt")

PEF_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/PEF.txt")

cova_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/assessment_center/covariate.txt")

bmi_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/assessment_center/bmi_from_UKB.txt")

smoke_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/assessment_center/covariate_smoke.txt")

### combine the gp_clinical and covariates data.
spo2_joined = spo2_info %>%
  inner_join(cova_info, by = c('eid' = 'Wenjian')) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18) %>%
  left_join(bmi_info %>% select(feid, f21001_0_0) %>% rename(eid = feid, BMI = f21001_0_0)) %>%
  left_join(smoke_info %>% select(feid, f20160_0_0) %>% rename(eid = feid, smoke = f20160_0_0))

### write the final-cleaned spo2 into file.
data.table::fwrite(spo2_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/spo2.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### combine the gp_clinical and covariates data.
FEV1_joined = FEV1_info %>%
  inner_join(cova_info, by = c('eid' = 'Wenjian')) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18) %>%
  left_join(bmi_info %>% select(feid, f21001_0_0) %>% rename(eid = feid, BMI = f21001_0_0)) %>%
  left_join(smoke_info %>% select(feid, f20160_0_0) %>% rename(eid = feid, smoke = f20160_0_0))

### write the final-cleaned FEV1 into file.
data.table::fwrite(FEV1_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/FEV1.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### combine the gp_clinical and covariates data.
FVC_joined = FVC_info %>%
  inner_join(cova_info, by = c('eid' = 'Wenjian')) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18) %>%
  left_join(bmi_info %>% select(feid, f21001_0_0) %>% rename(eid = feid, BMI = f21001_0_0)) %>%
  left_join(smoke_info %>% select(feid, f20160_0_0) %>% rename(eid = feid, smoke = f20160_0_0))

### write the final-cleaned FVC into file.
data.table::fwrite(FVC_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/FVC.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### combine the gp_clinical and covariates data.
FEV1_FVC_joined = FEV1_FVC_info %>%
  inner_join(cova_info, by = c('eid' = 'Wenjian')) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18) %>%
  left_join(bmi_info %>% select(feid, f21001_0_0) %>% rename(eid = feid, BMI = f21001_0_0)) %>%
  left_join(smoke_info %>% select(feid, f20160_0_0) %>% rename(eid = feid, smoke = f20160_0_0))

### write the final-cleaned FEV1/FVC ratio into file.
data.table::fwrite(FEV1_FVC_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/FEV1_FVC.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### combine the gp_clinical and covariates data.
PEF_joined = PEF_info %>%
  inner_join(cova_info, by = c('eid' = 'Wenjian')) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18) %>%
  left_join(bmi_info %>% select(feid, f21001_0_0) %>% rename(eid = feid, BMI = f21001_0_0)) %>%
  left_join(smoke_info %>% select(feid, f20160_0_0) %>% rename(eid = feid, smoke = f20160_0_0))

### write the final-cleaned PEF into file.
data.table::fwrite(PEF_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/PEF.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)
