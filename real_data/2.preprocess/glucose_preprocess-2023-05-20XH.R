
###### pre-clean the glucose gp_clinical table.
library(dplyr)
library(tidyr)

### read in the fast glucose, random glucose and HbA1c phecode.
source("/gdata01/user/xuhe/family_relatedness/real_data-2022-12-29/1.extract_pheno/phecode-2023-05-10XH.R")

### load the cleaned gp_clinical table.
gp_clinical = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/cleaned_gp_clinical.txt")

### extract fastglu from gp_clinical data.
fastglu = gp_clinical %>%
  filter(grepl(fastgluc_codes, code)) %>% 
  mutate(value1_numeric = as.numeric(value1), value2_numeric = as.numeric(value2), value3_numeric = as.numeric(value3)) %>%
  mutate(fglu = coalesce(value1_numeric, value2_numeric, value3_numeric)) %>%
  filter(fglu > 1 & fglu < 40) %>%
  mutate(value3 = toupper(value3)) %>%
  filter(value3 %in% c("", "MEA000", "MEA096", "MMOL/L", "UNKNOWN")) %>%
  group_by(eid) %>%
  mutate(mi = n()) %>%
  mutate(mean = mean(fglu)) %>%
  mutate(sd = ifelse(mi == 1, 0, 
                     ifelse(sd(fglu) > 1.5, 1.5, sd(fglu)))) %>%
  mutate(lower_tail = ifelse(mi == 1, mean, mean - 4 * sd),
         higher_tail = ifelse(mi == 1, mean, mean + 4 * sd)) %>%
  filter(fglu >= lower_tail & fglu <= higher_tail) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  mutate(fglu_new = ifelse(rep == 1, fglu, 
                           ifelse(max(fglu) - min(fglu) <= 1.5, 
                                  mean(fglu), fglu[which.min(abs(fglu - mean))]))) %>%
  mutate(fglu_new = round(fglu_new, 1)) %>%
  select(eid, data_provider, event_dt, fglu_new) %>% 
  dplyr::rename(fglu = fglu_new) %>%
  distinct(eid, event_dt, .keep_all = TRUE)

### write the pre-cleaned fast glucose into file.
data.table::fwrite(fastglu %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/fastglu.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### extract randomglu from gp_clinical data.
randomglu = gp_clinical %>%
  filter(grepl(randgluc_codes, code)) %>% 
  mutate(value1_numeric = as.numeric(value1), value2_numeric = as.numeric(value2), value3_numeric = as.numeric(value3)) %>%
  mutate(rglu = coalesce(value1_numeric, value2_numeric, value3_numeric)) %>%
  filter(rglu > 1 & rglu < 40) %>%
  mutate(value3 = toupper(value3)) %>%
  filter(value3 %in% c("", "MEA000", "MEA096", "MMOL/L", "UNKNOWN")) %>%
  group_by(eid) %>%
  mutate(mi = n()) %>%
  mutate(mean = mean(rglu)) %>%
  mutate(sd = ifelse(mi == 1, 0, 
                     ifelse(sd(rglu) > 2.5, 2.5, sd(rglu)))) %>%
  mutate(lower_tail = ifelse(mi == 1, mean, mean - 4 * sd),
         higher_tail = ifelse(mi == 1, mean, mean + 4 * sd)) %>%
  filter(rglu >= lower_tail & rglu <= higher_tail) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  mutate(rglu_new = ifelse(rep == 1, rglu, 
                           ifelse(max(rglu) - min(rglu) <= 2.5, 
                                  mean(rglu), rglu[which.min(abs(rglu - mean))]))) %>%
  mutate(rglu_new = round(rglu_new, 1)) %>%
  select(eid, data_provider, event_dt, rglu_new) %>% 
  dplyr::rename(rglu = rglu_new) %>%
  distinct(eid, event_dt, .keep_all = TRUE)

### write the pre-cleaned random glucose into file.
data.table::fwrite(randomglu %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/randomglu.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### extract hba1c from gp_clinical data.
hba1c = gp_clinical %>%
  filter(grepl(a1c_codes, code)) %>%
  mutate(value1_numeric = as.numeric(value1), value2_numeric = as.numeric(value2), value3_numeric = as.numeric(value3)) %>%
  mutate(hba1c = coalesce(value1_numeric, value2_numeric, value3_numeric)) %>%
  filter(hba1c > 0) %>%
  mutate(value3 = toupper(value3)) %>%
  mutate(value3 = ifelse(value3 %in% c("MEA000", "MEA097", "UNKNOWN", "MEA001", "%", "HBA1C", 
                                       "%TOTAL HB", "% TOTAL HB", "MEA215", "MMOL/MOL HB", "PER CENT", "%TOTAL"), "", value3)) %>%
  mutate(units = ifelse(value3 != "", value3, 
                        ifelse(code %in% c("XaPbt", "42W5."), "MMOL/MOL", "%"))) %>%
  filter(units %in% c("%", "MMOL/MOL")) %>%
  mutate(hba1c_percent = ifelse(units == "%", round(hba1c, 1), round(hba1c/10.929 + 2.15, 1))) %>%
  mutate(hba1c_mmol_mol = ifelse(units =="%", round(10.929 * (hba1c - 2.15), 1), round(hba1c, 1))) %>%
  filter(hba1c_percent > 3 & hba1c_percent < 18) %>%
  group_by(eid) %>%
  mutate(mi = n()) %>%
  mutate(mean = mean(hba1c_percent)) %>%
  mutate(sd = ifelse(mi == 1, 0, 
                     ifelse(sd(hba1c_percent) > 1.25, 1.25, sd(hba1c_percent)))) %>%
  mutate(lower_tail = ifelse(mi == 1, mean, mean - 4 * sd),
         higher_tail = ifelse(mi == 1, mean, mean + 4 * sd)) %>%
  filter(hba1c_percent >= lower_tail & hba1c_percent <= higher_tail) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  mutate(hba1c_percent_new = ifelse(rep == 1, hba1c_percent, 
                                    ifelse(max(hba1c_percent) - min(hba1c_percent) <= 0.5, 
                                           mean(hba1c_percent), hba1c_percent[which.min(abs(hba1c_percent - mean))]))) %>%
  mutate(hba1c_percent_new = round(hba1c_percent_new, 1)) %>%
  select(eid, data_provider, event_dt, hba1c_percent_new) %>%
  dplyr::rename(hba1c = hba1c_percent_new) %>%
  distinct(eid, event_dt, .keep_all = TRUE)

### write the pre-cleaned hba1c into file.
data.table::fwrite(hba1c %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/hba1c.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)


###### merge the gp_clinical table and covariates.
library(dplyr)
library(tidyr)
library(lubridate)

### read in the gp_clinical table and covariates.
fastglu_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/fastglu.txt")

randomglu_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/randomglu.txt")

hba1c_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/hba1c.txt")

cova_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/assessment_center/covariate.txt")

bmi_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/assessment_center/bmi_from_UKB.txt")

drug_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/assessment_center/covariate_drugs.txt")

### combine the gp_clinical and covariates data.
fastglu_joined = fastglu_info %>%
  inner_join(cova_info, by = c('eid' = 'Wenjian')) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18) %>%
  left_join(drug_info %>% select(feid, self_insulin), by = c('eid' = 'feid')) %>%
  left_join(bmi_info %>% select(feid, f21001_0_0) %>% rename(eid = feid, BMI = f21001_0_0))

### write the final-cleaned fasting glucose into file.
data.table::fwrite(fastglu_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/fastglu.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### combine the gp_clinical and covariates data.
randomglu_joined = randomglu_info %>%
  inner_join(cova_info, by = c('eid' = 'Wenjian')) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18) %>%
  left_join(drug_info %>% select(feid, self_insulin), by = c('eid' = 'feid')) %>%
  left_join(bmi_info %>% select(feid, f21001_0_0) %>% rename(eid = feid, BMI = f21001_0_0))

### write the final-cleaned random glucose into file.
data.table::fwrite(randomglu_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/randomglu.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### combine the gp_clinical and covariates data.
hba1c_joined = hba1c_info %>%
  inner_join(cova_info, by = c('eid' = 'Wenjian')) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18) %>%
  left_join(drug_info %>% select(feid, self_insulin), by = c('eid' = 'feid')) %>%
  left_join(bmi_info %>% select(feid, f21001_0_0) %>% rename(eid = feid, BMI = f21001_0_0))

### write the final-cleaned hba1c into file.
data.table::fwrite(hba1c_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/hba1c.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)
