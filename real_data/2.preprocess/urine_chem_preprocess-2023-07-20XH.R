
###### pre-clean blood protein phenotype gp_clinical table.
library(dplyr)
library(tidyr)

### read in the urine albumin, urine creatinine, UACR phecode.
source("/gdata01/user/xuhe/family_relatedness/real_data-2022-12-29/1.extract_pheno/phecode-2023-05-10XH.R")

### load the cleaned gp_clinical table.
gp_clinical = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/cleaned_gp_clinical.txt")

extracted_records = gp_clinical %>%
  filter(grepl(paste(urine_albumin_codes, urine_creatinine_codes, UACR_codes, sep = "|"), code)) %>%
  mutate(value1_numeric = as.numeric(value1), value2_numeric = as.numeric(value2), value3_numeric = as.numeric(value3)) %>%
  mutate(value = coalesce(value1_numeric, value2_numeric, value3_numeric)) %>%
  mutate(isualb = grepl(urine_albumin_codes, code)) %>%
  mutate(isucr = grepl(urine_creatinine_codes, code)) %>%
  mutate(isuacr = grepl(UACR_codes, code)) %>%
  filter(value > 0) %>% 
  distinct(eid, event_dt, value, isualb, isucr, isuacr, .keep_all = TRUE) %>%
  group_by(eid, event_dt) %>%
  mutate(n = n())

urine_albumin = extracted_records %>% 
  filter(isualb == TRUE) %>%
  filter(n > 1) %>%
  filter(value < 5000) %>%
  mutate(value3 = toupper(value3)) %>%
  filter(value3 %in% c("", "MEA000", "MEA083", "MG/L")) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  mutate(value_new = ifelse(rep == 1, value, 
                            ifelse(max(value) - min(value) <= 1, mean(value), 0))) %>%
  filter(value_new > 0) %>%
  mutate(value_new = round(value_new, 3)) %>%
  select(eid, data_provider, event_dt, value_new) %>% 
  dplyr::rename(value = value_new) %>%
  distinct(eid, event_dt, .keep_all = TRUE)

urine_creatinine = extracted_records %>% 
  filter(isucr == TRUE) %>%
  filter(n > 1) %>%
  filter(value < 100) %>%
  mutate(value3 = toupper(value3)) %>%
  filter(value3 %in% c("", "MEA000", "MEA096", "MMOL/L")) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  mutate(value_new = ifelse(rep == 1, value, 
                            ifelse(max(value) - min(value) <= 1, mean(value), 0))) %>%
  filter(value_new > 0) %>%
  mutate(value_new = round(value_new, 3)) %>%
  select(eid, data_provider, event_dt, value_new) %>% 
  dplyr::rename(value = value_new) %>%
  distinct(eid, event_dt, .keep_all = TRUE)

UACR = extracted_records %>% 
  filter(isuacr == TRUE) %>%
  filter(n > 1) %>%
  mutate(value3 = toupper(value3)) %>%
  filter(!(value3 %in% c("MEA083", "MEA096", "MEA156", "MEA241", "MG/L", 
                         "M1/MIN", "ML/MIN", "MMOL/L", "NG/ML", "UMOL/L", "UNKNOWN"))) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  mutate(value_new = ifelse(rep == 1, value, 
                            ifelse(max(value) - min(value) <= 1, mean(value), 0))) %>%
  filter(value_new > 0) %>%
  mutate(value_new = round(value_new, 3)) %>%
  select(eid, data_provider, event_dt, value_new) %>% 
  dplyr::rename(value = value_new) %>%
  distinct(eid, event_dt, .keep_all = TRUE)

combined_records = full_join(urine_albumin %>% rename(ualb_report = value),
                            urine_creatinine %>% rename(ucr_report = value)) %>%
  full_join(UACR %>% rename(uacr_report = value)) %>%
  rowwise %>%
  mutate(n = sum(!is.na(ualb_report), !is.na(ucr_report), !is.na(uacr_report))) %>%
  filter(n > 1)
  
combined_records_3  = combined_records %>%
  filter(n == 3) %>%
  mutate(uacr_calcu = round(ualb_report/ucr_report, 1)) %>%
  filter(abs(uacr_report - uacr_calcu) < 1) %>%
  select(eid, data_provider, event_dt, ualb_report, ucr_report, uacr_calcu) %>%
  rename(ualb = ualb_report, ucr = ucr_report, uacr = uacr_calcu)

combined_records_2  = combined_records %>%
  filter(n == 2) %>%
  mutate(uacr_calcu = ifelse(is.na(uacr_report), round(ualb_report / ucr_report, 1), uacr_report)) %>%
  mutate(ualb_calcu = ifelse(is.na(ualb_report), round(ucr_report * uacr_report, 1), ualb_report)) %>%
  mutate(ucr_calcu = ifelse(is.na(ucr_report), round(ualb_report/uacr_report, 1), ucr_report)) %>%
  select(eid, data_provider, event_dt, ualb_calcu, ucr_calcu, uacr_calcu) %>%
  rename(ualb = ualb_calcu, ucr = ucr_calcu, uacr = uacr_calcu)

merged_records = rbind(combined_records_3, combined_records_2) %>%
  ungroup %>%
  arrange(eid, event_dt) %>%
  filter(uacr < 30) %>%
  filter(ualb > 0.1 & ualb < 300) %>%
  filter(ucr > 0.8 & ucr < 30) %>%
  group_by(eid) %>%
  mutate(mi = n()) %>%
  mutate(mean = mean(uacr)) %>%
  mutate(sd = ifelse(mi == 1, 0, 
                     ifelse(sd(uacr) > 4, 4, sd(uacr)))) %>%
  mutate(lower_tail = ifelse(mi == 1, mean, mean - 4 * sd),
         higher_tail = ifelse(mi == 1, mean, mean + 4 * sd)) %>%
  filter(uacr >= lower_tail & uacr <= higher_tail)

### write the pre-cleaned UACR into file.
precleaned_UACR = merged_records %>%
  select(eid, data_provider, event_dt, uacr)

data.table::fwrite(precleaned_UACR %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/UACR.txt", 
                   row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

### write the pre-cleaned urine albumin into file.
precleaned_urine_albumin = merged_records %>%
  select(eid, data_provider, event_dt, ualb)

data.table::fwrite(precleaned_urine_albumin %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/urine_albumin.txt", 
                   row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

### write the pre-cleaned urine creatinine into file.
precleaned_urine_creatinine = merged_records %>%
  select(eid, data_provider, event_dt, ucr)

data.table::fwrite(precleaned_urine_creatinine %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/urine_creatinine.txt", 
                   row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)


###### merge the gp_clinical table and covariates.
library(dplyr)
library(tidyr)
library(lubridate)

### read in the gp_clinical table and covariates.
urine_albumin_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/urine_albumin.txt")

urine_creatinine_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/urine_creatinine.txt")

UACR_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/UACR.txt")

cova_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/assessment_center/covariate.txt")

bmi_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/assessment_center/bmi_from_UKB.txt")

### combine the gp_clinical and covariates data.
urine_albumin_joined = urine_albumin_info %>%
  inner_join(cova_info, by = c('eid' = 'Wenjian')) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18) %>%
  left_join(bmi_info %>% select(feid, f21001_0_0) %>% rename(eid = feid, BMI = f21001_0_0))

### write the final-cleaned urine_albumin into file.
data.table::fwrite(urine_albumin_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/urine_albumin.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### combine the gp_clinical and covariates data.
urine_creatinine_joined = urine_creatinine_info %>%
  inner_join(cova_info, by = c('eid' = 'Wenjian')) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18) %>%
  left_join(bmi_info %>% select(feid, f21001_0_0) %>% rename(eid = feid, BMI = f21001_0_0))

### write the final-cleaned urine_creatinine into file.
data.table::fwrite(urine_creatinine_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/urine_creatinine.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### combine the gp_clinical and covariates data.
UACR_joined = UACR_info %>%
  inner_join(cova_info, by = c('eid' = 'Wenjian')) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18) %>%
  left_join(bmi_info %>% select(feid, f21001_0_0) %>% rename(eid = feid, BMI = f21001_0_0))

### write the final-cleaned UACR into file.
data.table::fwrite(UACR_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/UACR.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)
