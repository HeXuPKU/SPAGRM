
###### pre-clean blood protein phenotype gp_clinical table.
library(dplyr)
library(tidyr)

### read in the total protein, blood globulin, blood albumin phecode.
source("/gdata01/user/xuhe/family_relatedness/real_data-2022-12-29/1.extract_pheno/phecode-2023-05-10XH.R")

### load the cleaned gp_clinical table.
gp_clinical = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/cleaned_gp_clinical.txt")

### extract total protein from gp_clinical data.
total_protein = gp_clinical %>%
  filter(grepl(total_protein_codes, code)) %>% 
  mutate(value1_numeric = as.numeric(value1), value2_numeric = as.numeric(value2), value3_numeric = as.numeric(value3)) %>%
  mutate(pro = coalesce(value1_numeric, value2_numeric, value3_numeric)) %>%
  filter(pro > 30 & pro < 120) %>%
  mutate(value3 = toupper(value3)) %>%
  filter(value3 %in% c("", "G/L", "MEA000", "MEA057")) %>%
  group_by(eid) %>%
  mutate(mi = n()) %>%
  mutate(mean = mean(pro)) %>%
  mutate(sd = ifelse(mi == 1, 0, 
                     ifelse(sd(pro) > 4, 4, sd(pro)))) %>%
  mutate(lower_tail = ifelse(mi == 1, mean, mean - 4 * sd),
         higher_tail = ifelse(mi == 1, mean, mean + 4 * sd)) %>%
  filter(pro >= lower_tail & pro <= higher_tail) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  mutate(pro_new = ifelse(rep == 1, pro, 
                          ifelse(max(pro) - min(pro) <= 4, 
                                 mean(pro), pro[which.min(abs(pro - mean))]))) %>%
  mutate(pro_new = round(pro_new, 0)) %>%
  select(eid, data_provider, event_dt, pro_new) %>% 
  dplyr::rename(pro = pro_new) %>%
  distinct(eid, event_dt, .keep_all = TRUE)

### write the pre-cleaned total protein into file.
data.table::fwrite(total_protein %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/total_protein.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### extract blood globulin from gp_clinical data.
blood_globulin = gp_clinical %>%
  filter(grepl(blood_globulin_codes, code)) %>% 
  mutate(value1_numeric = as.numeric(value1), value2_numeric = as.numeric(value2), value3_numeric = as.numeric(value3)) %>%
  mutate(glo = coalesce(value1_numeric, value2_numeric, value3_numeric)) %>%
  filter(glo > 10 & glo < 60) %>%
  mutate(value3 = toupper(value3)) %>%
  filter(value3 %in% c("", "G/L", "MEA000", "MEA057")) %>%
  group_by(eid) %>%
  mutate(mi = n()) %>%
  mutate(mean = mean(glo)) %>%
  mutate(sd = ifelse(mi == 1, 0, 
                     ifelse(sd(glo) > 4, 4, sd(glo)))) %>%
  mutate(lower_tail = ifelse(mi == 1, mean, mean - 4 * sd),
         higher_tail = ifelse(mi == 1, mean, mean + 4 * sd)) %>%
  filter(glo >= lower_tail & glo <= higher_tail) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  mutate(glo_new = ifelse(rep == 1, glo, 
                          ifelse(max(glo) - min(glo) <= 4, 
                                 mean(glo), glo[which.min(abs(glo - mean))]))) %>%
  mutate(glo_new = round(glo_new, 0)) %>%
  select(eid, data_provider, event_dt, glo_new) %>% 
  dplyr::rename(glo = glo_new) %>%
  distinct(eid, event_dt, .keep_all = TRUE)

### write the pre-cleaned blood globulin into file.
data.table::fwrite(blood_globulin %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/blood_globulin.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### extract blood albumin from gp_clinical data.
blood_albumin = gp_clinical %>%
  filter(grepl(blood_albumin_codes, code)) %>% 
  mutate(value1_numeric = as.numeric(value1), value2_numeric = as.numeric(value2), value3_numeric = as.numeric(value3)) %>%
  mutate(alb = coalesce(value1_numeric, value2_numeric, value3_numeric)) %>%
  filter(alb > 10 & alb < 60) %>%
  mutate(value3 = toupper(value3)) %>%
  filter(value3 %in% c("", "G/L", "MEA000", "MEA057")) %>%
  group_by(eid) %>%
  mutate(mi = n()) %>%
  mutate(mean = mean(alb)) %>%
  mutate(sd = ifelse(mi == 1, 0, 
                     ifelse(sd(alb) > 4, 4, sd(alb)))) %>%
  mutate(lower_tail = ifelse(mi == 1, mean, mean - 4 * sd),
         higher_tail = ifelse(mi == 1, mean, mean + 4 * sd)) %>%
  filter(alb >= lower_tail & alb <= higher_tail) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  mutate(alb_new = ifelse(rep == 1, alb, 
                          ifelse(max(alb) - min(alb) <= 4, 
                                 mean(alb), alb[which.min(abs(alb - mean))]))) %>%
  mutate(alb_new = round(alb_new, 0)) %>%
  select(eid, data_provider, event_dt, alb_new) %>% 
  dplyr::rename(alb = alb_new) %>%
  distinct(eid, event_dt, .keep_all = TRUE)

### write the pre-cleaned blood albumin into file.
data.table::fwrite(blood_albumin %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/blood_albumin.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)


###### merge the gp_clinical table and covariates.
library(dplyr)
library(tidyr)
library(lubridate)

### read in the gp_clinical table and covariates.
total_protein_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/total_protein.txt")

blood_globulin_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/blood_globulin.txt")

blood_albumin_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/blood_albumin.txt")

cova_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/assessment_center/covariate.txt")

bmi_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/assessment_center/bmi_from_UKB.txt")

### combine the gp_clinical and covariates data.
total_protein_joined = total_protein_info %>%
  inner_join(cova_info, by = c('eid' = 'Wenjian')) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18) %>%
  left_join(bmi_info %>% select(feid, f21001_0_0) %>% rename(eid = feid, BMI = f21001_0_0))

### write the final-cleaned blood total protein into file.
data.table::fwrite(total_protein_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/total_protein.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### combine the gp_clinical and covariates data.
blood_globulin_joined = blood_globulin_info %>%
  inner_join(cova_info, by = c('eid' = 'Wenjian')) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18) %>%
  left_join(bmi_info %>% select(feid, f21001_0_0) %>% rename(eid = feid, BMI = f21001_0_0))

### write the final-cleaned blood globulin into file.
data.table::fwrite(blood_globulin_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/blood_globulin.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### combine the gp_clinical and covariates data.
blood_albumin_joined = blood_albumin_info %>%
  inner_join(cova_info, by = c('eid' = 'Wenjian')) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18) %>%
  left_join(bmi_info %>% select(feid, f21001_0_0) %>% rename(eid = feid, BMI = f21001_0_0))

### write the final-cleaned blood globulin into file.
data.table::fwrite(blood_albumin_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/blood_albumin.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)
