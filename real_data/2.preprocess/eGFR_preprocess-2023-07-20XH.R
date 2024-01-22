
###### pre-clean blood protein phenotype gp_clinical table.
library(dplyr)
library(tidyr)

### read in the blood creatinine, eGFR phecode.
source("/gdata01/user/xuhe/family_relatedness/real_data-2022-12-29/1.extract_pheno/phecode-2023-05-10XH.R")

### load the cleaned gp_clinical table.
gp_clinical = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/cleaned_gp_clinical.txt")

### extract blood creatinine from gp_clinical data.
blood_creatinine = gp_clinical %>%
  filter(grepl(blood_creatinine_codes, code)) %>% 
  filter(value3 != "\xb5mol/L") %>%
  mutate(value1_numeric = as.numeric(value1), value2_numeric = as.numeric(value2), value3_numeric = as.numeric(value3)) %>%
  mutate(crea = coalesce(value1_numeric, value2_numeric, value3_numeric)) %>%
  filter(crea > 20 & crea < 300) %>%
  mutate(value3 = toupper(value3)) %>%
  filter(value3 %in% c("", "MEA000", "MEA142", "UMOL/L", "UNKNOWN")) %>%
  group_by(eid) %>%
  mutate(mi = n()) %>%
  mutate(mean = mean(crea)) %>%
  mutate(sd = ifelse(mi == 1, 0, 
                     ifelse(sd(crea) > 20, 20, sd(crea)))) %>%
  mutate(lower_tail = ifelse(mi == 1, mean, mean - 4 * sd),
         higher_tail = ifelse(mi == 1, mean, mean + 4 * sd)) %>%
  filter(crea >= lower_tail & crea <= higher_tail) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  mutate(crea_new = ifelse(rep == 1, crea, 
                          ifelse(max(crea) - min(crea) <= 5, mean(crea), 0))) %>%
  filter(crea_new > 0) %>%
  mutate(crea_new = round(crea_new, 0)) %>%
  select(eid, data_provider, event_dt, crea_new) %>% 
  dplyr::rename(crea = crea_new) %>%
  distinct(eid, event_dt, .keep_all = TRUE)

### write the pre-cleaned blood creatinine into file.
data.table::fwrite(blood_creatinine %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/blood_creatinine.txt", 
                   row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

### since eGFR were calculated by different formulas, we derived the eGFR using blood creatinine.
library(dplyr)
library(tidyr)
library(lubridate)

blood_creatinine_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/blood_creatinine.txt")

cova_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/assessment_center/covariate.txt")

bmi_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/assessment_center/bmi_from_UKB.txt")

blood_creatinine_joined = blood_creatinine_info %>%
  inner_join(cova_info, by = c('eid' = 'Wenjian')) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18) %>%
  rowwise() %>%
  mutate(egfr = ifelse(sex_genetic == 1, 141*min(crea/88.4/0.9, 1)^(-0.329)*max(crea/88.4/0.9, 1)^(-1.209)*0.993^age,
                       141*min(crea/88.4/0.7, 1)^(-0.411)*max(crea/88.4/0.7, 1)^(-1.209)*0.993^age*1.018)) %>%
  mutate(egfr = round(egfr, 1)) %>%
  filter(egfr > 15 & egfr < 200) %>% # https://www.nature.com/articles/s41467-021-24491-0
  left_join(bmi_info %>% select(feid, f21001_0_0) %>% rename(eid = feid, BMI = f21001_0_0))

### write the final-cleaned blood_creatinine into file.
data.table::fwrite(blood_creatinine_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/eGFR.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)
