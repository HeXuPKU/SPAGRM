
###### pre-clean blood protein phenotype gp_clinical table.
library(dplyr)
library(tidyr)

### read in the blood TSH, FT3, FT4 phecode.
source("/gdata01/user/xuhe/family_relatedness/real_data-2022-12-29/1.extract_pheno/phecode-2023-05-10XH.R")

### load the cleaned gp_clinical table.
gp_clinical = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/cleaned_gp_clinical.txt")

### extract TSH from gp_clinical data.
TSH = gp_clinical %>%
  filter(grepl(TSH_codes, code)) %>% 
  mutate(value1_numeric = as.numeric(value1), value2_numeric = as.numeric(value2), value3_numeric = as.numeric(value3)) %>%
  mutate(tsh = coalesce(value1_numeric, value2_numeric, value3_numeric)) %>%
  filter(tsh >= 0.01 & tsh < 10) %>% # slightly different from https://www.medrxiv.org/content/10.1101/2022.12.22.22283779v2.full-text
  mutate(value3 = toupper(value3)) %>%
  filter(value3 %in% c("", "MEA000", "MEA149", "MEA160", "MEA170", "MICROU/L", "MIU/L", "MU/L")) %>%
  group_by(eid) %>%
  mutate(mi = n()) %>%
  mutate(mean = mean(tsh)) %>%
  mutate(sd = ifelse(mi == 1, 0, 
                     ifelse(sd(tsh) > 2, 2, sd(tsh)))) %>%
  mutate(lower_tail = ifelse(mi == 1, mean, mean - 4 * sd),
         higher_tail = ifelse(mi == 1, mean, mean + 4 * sd)) %>%
  filter(tsh >= lower_tail & tsh <= higher_tail) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  mutate(tsh_new = ifelse(rep == 1, tsh, 
                           ifelse(max(tsh) - min(tsh) <= 0.2, mean(tsh), 0))) %>%
  filter(tsh_new > 0) %>%
  mutate(tsh_new = round(tsh_new, 2)) %>%
  select(eid, data_provider, event_dt, tsh_new) %>% 
  dplyr::rename(tsh = tsh_new) %>%
  distinct(eid, event_dt, .keep_all = TRUE)

### write the pre-cleaned TSH into file.
data.table::fwrite(TSH %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/TSH.txt", 
                   row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

### extract T3 from gp_clinical data.
FT3 = gp_clinical %>%
  filter(grepl(FT3_codes, code)) %>% 
  mutate(value1_numeric = as.numeric(value1), value2_numeric = as.numeric(value2), value3_numeric = as.numeric(value3)) %>%
  mutate(ft3 = coalesce(value1_numeric, value2_numeric, value3_numeric)) %>%
  filter(ft3 > 1 & ft3 < 10) %>%
  mutate(value3 = toupper(value3)) %>%
  filter(value3 %in% c("", "MEA000", "MEA120", "PMOL/L")) %>%
  group_by(eid) %>%
  mutate(mi = n()) %>%
  mutate(mean = mean(ft3)) %>%
  mutate(sd = ifelse(mi == 1, 0, 
                     ifelse(sd(ft3) > 1.2, 1.2, sd(ft3)))) %>%
  mutate(lower_tail = ifelse(mi == 1, mean, mean - 4 * sd),
         higher_tail = ifelse(mi == 1, mean, mean + 4 * sd)) %>%
  filter(ft3 >= lower_tail & ft3 <= higher_tail) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  mutate(ft3_new = ifelse(rep == 1, ft3, 
                          ifelse(max(ft3) - min(ft3) <= 1, mean(ft3), 0))) %>%
  filter(ft3_new > 0) %>%
  mutate(ft3_new = round(ft3_new, 1)) %>%
  select(eid, data_provider, event_dt, ft3_new) %>% 
  dplyr::rename(ft3 = ft3_new) %>%
  distinct(eid, event_dt, .keep_all = TRUE)

### write the pre-cleaned FT3 into file.
data.table::fwrite(FT3 %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/FT3.txt", 
                   row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

### extract FT4 from gp_clinical data.
FT4 = gp_clinical %>%
  filter(grepl(FT4_codes, code)) %>% 
  mutate(value1_numeric = as.numeric(value1), value2_numeric = as.numeric(value2), value3_numeric = as.numeric(value3)) %>%
  mutate(ft4 = coalesce(value1_numeric, value2_numeric, value3_numeric)) %>%
  filter(ft4 > 2 & ft4 < 40) %>%
  mutate(value3 = toupper(value3)) %>%
  filter(value3 %in% c("", "MEA000", "MEA120", "PMOL/L")) %>%
  group_by(eid) %>%
  mutate(mi = n()) %>%
  mutate(mean = mean(ft4)) %>%
  mutate(sd = ifelse(mi == 1, 0, 
                     ifelse(sd(ft4) > 4, 4, sd(ft4)))) %>%
  mutate(lower_tail = ifelse(mi == 1, mean, mean - 4 * sd),
         higher_tail = ifelse(mi == 1, mean, mean + 4 * sd)) %>%
  filter(ft4 >= lower_tail & ft4 <= higher_tail) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  mutate(ft4_new = ifelse(rep == 1, ft4, 
                          ifelse(max(ft4) - min(ft4) <= 2, mean(ft4), 0))) %>%
  filter(ft4_new > 0) %>%
  mutate(ft4_new = round(ft4_new, 1)) %>%
  select(eid, data_provider, event_dt, ft4_new) %>% 
  dplyr::rename(ft4 = ft4_new) %>%
  distinct(eid, event_dt, .keep_all = TRUE)

### write the pre-cleaned FT4 into file.
data.table::fwrite(FT4 %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/FT4.txt", 
                   row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)


###### merge the gp_clinical table and covariates.
library(dplyr)
library(tidyr)
library(lubridate)

### read in the gp_clinical table and covariates.
TSH_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/TSH.txt")

FT3_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/FT3.txt")

FT4_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/FT4.txt")

cova_TSH_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/assessment_center/covariate_TSH_operations.txt")

cova_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/assessment_center/covariate.txt")

bmi_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/assessment_center/bmi_from_UKB.txt")

cova_TSH_id = cova_TSH_info %>% filter(opera == 1)

### combine the gp_clinical and covariates data.
TSH_joined = TSH_info %>%
  filter(!(eid %in% cova_TSH_id$feid)) %>%
  inner_join(cova_info, by = c('eid' = 'Wenjian')) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18) %>%
  left_join(bmi_info %>% select(feid, f21001_0_0) %>% rename(eid = feid, BMI = f21001_0_0))

### write the final-cleaned TSH into file.
data.table::fwrite(TSH_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/TSH.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### combine the gp_clinical and covariates data.
FT3_joined = FT3_info %>%
  filter(!(eid %in% cova_TSH_id$feid)) %>%
  inner_join(cova_info, by = c('eid' = 'Wenjian')) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18) %>%
  left_join(bmi_info %>% select(feid, f21001_0_0) %>% rename(eid = feid, BMI = f21001_0_0))

### write the final-cleaned FT3 into file.
data.table::fwrite(FT3_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/FT3.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### combine the gp_clinical and covariates data.
FT4_joined = FT4_info %>%
  filter(!(eid %in% cova_TSH_id$feid)) %>%
  inner_join(cova_info, by = c('eid' = 'Wenjian')) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18) %>%
  left_join(bmi_info %>% select(feid, f21001_0_0) %>% rename(eid = feid, BMI = f21001_0_0))

### write the final-cleaned FT4 into file.
data.table::fwrite(FT4_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/FT4.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)
