
###### pre-clean blood immunity phenotype gp_clinical table.
library(dplyr)
library(tidyr)

### read in the C-reactive protein (CRP), Prostate-specific antigen (PSA), creatine kinase (CK), amylase, Carbohydrate antigen 125 phecode.
source("/gdata01/user/xuhe/family_relatedness/real_data-2022-12-29/1.extract_pheno/phecode-2023-05-10XH.R")

### load the cleaned gp_clinical table.
gp_clinical = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/cleaned_gp_clinical.txt")

### extract C-reactive protein from gp_clinical data.
CRP = gp_clinical %>%
  filter(grepl(CRP_codes, code)) %>% 
  mutate(value1_numeric = as.numeric(value1), value2_numeric = as.numeric(value2), value3_numeric = as.numeric(value3)) %>%
  mutate(crp = coalesce(value1_numeric, value2_numeric, value3_numeric)) %>%
  filter(crp >= 0.1 & crp < 80) %>%
  mutate(value3 = toupper(value3)) %>%
  filter(value3 %in% c("", "MG/L", "MEA000", "MEA083")) %>%
  group_by(eid) %>%
  mutate(mi = n()) %>%
  mutate(mean = mean(crp)) %>%
  mutate(sd = ifelse(mi == 1, 0, 
                     ifelse(sd(crp) > 10, 10, sd(crp)))) %>%
  mutate(lower_tail = ifelse(mi == 1, mean, mean - 4 * sd),
         higher_tail = ifelse(mi == 1, mean, mean + 4 * sd)) %>%
  filter(crp >= lower_tail & crp <= higher_tail) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  mutate(crp_new = ifelse(rep == 1, crp, 
                          ifelse(max(crp) - min(crp) <= 1, 
                                 mean(crp), crp[which.min(abs(crp - mean))]))) %>%
  mutate(crp_new = round(crp_new, 1)) %>%
  select(eid, data_provider, event_dt, crp_new) %>% 
  dplyr::rename(crp = crp_new) %>%
  distinct(eid, event_dt, .keep_all = TRUE)

### write the pre-cleaned C-reactive protein into file.
data.table::fwrite(CRP %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/CRP.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### extract Prostate-specific antigen (PSA) from gp_clinical data.
PSA = gp_clinical %>%
  filter(grepl(PSA_codes, code)) %>% 
  filter(value3 != '\xb5g/L') %>%
  mutate(value1_numeric = as.numeric(value1), value2_numeric = as.numeric(value2), value3_numeric = as.numeric(value3)) %>%
  mutate(psa = coalesce(value1_numeric, value2_numeric, value3_numeric)) %>%
  filter(psa > 0.01 & psa < 10) %>% # analyze low-PSA, see https://www.nature.com/articles/s41591-023-02277-9#Sec10
  mutate(value3 = toupper(value3)) %>%
  filter(value3 %in% c("", "UG/L", "NANOGRAMS/ML", "MEA000", "MEA083", "MEA106", "MEA133")) %>%
  group_by(eid) %>%
  mutate(mi = n()) %>%
  mutate(mean = mean(psa)) %>%
  mutate(sd = ifelse(mi == 1, 0, 
                     ifelse(sd(psa) > 2, 2, sd(psa)))) %>%
  mutate(lower_tail = ifelse(mi == 1, mean, mean - 4 * sd),
         higher_tail = ifelse(mi == 1, mean, mean + 4 * sd)) %>%
  filter(psa >= lower_tail & psa <= higher_tail) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  mutate(psa_new = ifelse(rep == 1, psa, 
                          ifelse(max(psa) - min(psa) <= 1, 
                                 mean(psa), psa[which.min(abs(psa - mean))]))) %>%
  mutate(psa_new = round(psa_new, 2)) %>%
  select(eid, data_provider, event_dt, psa_new) %>% 
  dplyr::rename(psa = psa_new) %>%
  distinct(eid, event_dt, .keep_all = TRUE)

### write the pre-cleaned PSA into file.
data.table::fwrite(PSA %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/PSA.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### extract blood Carbohydrate antigen 125 from gp_clinical data.
CA125 = gp_clinical %>%
  filter(grepl(CA125_codes, code)) %>% 
  mutate(value1_numeric = as.numeric(value1), value2_numeric = as.numeric(value2), value3_numeric = as.numeric(value3)) %>%
  mutate(ca = coalesce(value1_numeric, value2_numeric, value3_numeric)) %>%
  filter(ca > 0 & ca < 100) %>% 
  group_by(eid) %>%
  mutate(mi = n()) %>%
  mutate(mean = mean(ca)) %>%
  mutate(sd = ifelse(mi == 1, 0, 
                     ifelse(sd(ca) > 8, 8, sd(ca)))) %>%
  mutate(lower_tail = ifelse(mi == 1, mean, mean - 4 * sd),
         higher_tail = ifelse(mi == 1, mean, mean + 4 * sd)) %>%
  filter(ca >= lower_tail & ca <= higher_tail) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  mutate(ca_new = ifelse(rep == 1, ca, 
                         ifelse(max(ca) - min(ca) <= 1, 
                                mean(ca), ca[which.min(abs(ca - mean))]))) %>%
  mutate(ca_new = round(ca_new, 0)) %>%
  select(eid, data_provider, event_dt, ca_new) %>% 
  dplyr::rename(ca = ca_new) %>%
  distinct(eid, event_dt, .keep_all = TRUE)

### write the pre-cleaned blood Carbohydrate antigen 125 into file.
data.table::fwrite(CA125 %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/CA125.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### extract creatine kinase (CK) from gp_clinical data.
CK = gp_clinical %>%
  filter(grepl(CK_codes, code)) %>% 
  mutate(value1_numeric = as.numeric(value1), value2_numeric = as.numeric(value2), value3_numeric = as.numeric(value3)) %>%
  mutate(ck = coalesce(value1_numeric, value2_numeric, value3_numeric)) %>%
  filter(ck > 10 & ck < 1000) %>% 
  mutate(value3 = toupper(value3)) %>%
  filter(value3 %in% c("", "U/L", "IU/L", "MEA000", "MEA061", "MEA127")) %>%
  group_by(eid) %>%
  mutate(mi = n()) %>%
  mutate(mean = mean(ck)) %>%
  mutate(sd = ifelse(mi == 1, 0, 
                     ifelse(sd(ck) > 80, 80, sd(ck)))) %>%
  mutate(lower_tail = ifelse(mi == 1, mean, mean - 4 * sd),
         higher_tail = ifelse(mi == 1, mean, mean + 4 * sd)) %>%
  filter(ck >= lower_tail & ck <= higher_tail) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  mutate(ck_new = ifelse(rep == 1, ck, 
                         ifelse(max(ck) - min(ck) <= 1, 
                                mean(ck), ck[which.min(abs(ck - mean))]))) %>%
  mutate(ck_new = round(ck_new, 0)) %>%
  select(eid, data_provider, event_dt, ck_new) %>% 
  dplyr::rename(ck = ck_new) %>%
  distinct(eid, event_dt, .keep_all = TRUE)

### write the pre-cleaned CK into file.
data.table::fwrite(CK %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/CK.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### extract blood amylase from gp_clinical data.
amylase = gp_clinical %>%
  filter(grepl(amylase_codes, code)) %>% 
  mutate(value1_numeric = as.numeric(value1), value2_numeric = as.numeric(value2), value3_numeric = as.numeric(value3)) %>%
  mutate(amy = coalesce(value1_numeric, value2_numeric, value3_numeric)) %>%
  filter(amy > 10 & amy < 150) %>% 
  mutate(value3 = toupper(value3)) %>%
  filter(value3 %in% c("", "U/L", "IU/L", "MEA000", "MEA061", "MEA127")) %>%
  group_by(eid) %>%
  mutate(mi = n()) %>%
  mutate(mean = mean(amy)) %>%
  mutate(sd = ifelse(mi == 1, 0, 
                     ifelse(sd(amy) > 8, 8, sd(amy)))) %>%
  mutate(lower_tail = ifelse(mi == 1, mean, mean - 4 * sd),
         higher_tail = ifelse(mi == 1, mean, mean + 4 * sd)) %>%
  filter(amy >= lower_tail & amy <= higher_tail) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  mutate(amy_new = ifelse(rep == 1, amy, 
                          ifelse(max(amy) - min(amy) <= 1, 
                                 mean(amy), amy[which.min(abs(amy - mean))]))) %>%
  mutate(amy_new = round(amy_new, 0)) %>%
  select(eid, data_provider, event_dt, amy_new) %>% 
  dplyr::rename(amy = amy_new) %>%
  distinct(eid, event_dt, .keep_all = TRUE)

### write the pre-cleaned blood amylase into file.
data.table::fwrite(amylase %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/amylase.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)


###### merge the gp_clinical table and covariates.
library(dplyr)
library(tidyr)
library(lubridate)

### read in the gp_clinical table and covariates.
CRP_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/CRP.txt")

PSA_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/PSA.txt")

CA125_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/CA125.txt")

CK_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/CK.txt")

amylase_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/amylase.txt")

cova_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/assessment_center/covariate.txt")

cova_PSA_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/assessment_center/covariate_PSA_operations.txt")

bmi_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/assessment_center/bmi_from_UKB.txt")

### combine the gp_clinical and covariates data.
CRP_joined = CRP_info %>%
  inner_join(cova_info, by = c('eid' = 'Wenjian')) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18) %>%
  left_join(bmi_info %>% select(feid, f21001_0_0) %>% rename(eid = feid, BMI = f21001_0_0))

### write the final-cleaned CRP into file.
data.table::fwrite(CRP_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/CRP.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### combine the gp_clinical and covariates data.
PSA_opera_id = cova_PSA_info%>%filter(opera==1)

PSA_joined = PSA_info %>%
  filter(!(eid %in% PSA_opera_id$feid)) %>%
  inner_join(cova_info, by = c('eid' = 'Wenjian')) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18) %>%
  filter(sex_genetic == 1) %>%
  left_join(bmi_info %>% select(feid, f21001_0_0) %>% rename(eid = feid, BMI = f21001_0_0))

### write the final-cleaned PSA into file.
data.table::fwrite(PSA_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/PSA.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### combine the gp_clinical and covariates data.
CA125_joined = CA125_info %>%
  inner_join(cova_info, by = c('eid' = 'Wenjian')) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18) %>%
  filter(sex_genetic == 0) %>%
  left_join(bmi_info %>% select(feid, f21001_0_0) %>% rename(eid = feid, BMI = f21001_0_0))

### write the final-cleaned CA125 into file.
data.table::fwrite(CA125_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/CA125.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### combine the gp_clinical and covariates data.
CK_joined = CK_info %>%
  inner_join(cova_info, by = c('eid' = 'Wenjian')) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18) %>%
  left_join(bmi_info %>% select(feid, f21001_0_0) %>% rename(eid = feid, BMI = f21001_0_0))

### write the final-cleaned CK into file.
data.table::fwrite(CK_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/CK.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### combine the gp_clinical and covariates data.
amylase_joined = amylase_info %>%
  inner_join(cova_info, by = c('eid' = 'Wenjian')) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18) %>%
  left_join(bmi_info %>% select(feid, f21001_0_0) %>% rename(eid = feid, BMI = f21001_0_0))

### write the final-cleaned amylase into file.
data.table::fwrite(amylase_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/amylase.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)
