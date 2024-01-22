
library(dplyr)
library(tidyr)

### read in the pulse rate phecode.
source("/gdata01/user/xuhe/family_relatedness/real_data-2022-12-29/1.extract_pheno/phecode-2023-05-10XH.R")

### load the cleaned gp_clinical table
gp_clinical = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/cleaned_gp_clinical.txt")

### extract pulse rate.
pr = gp_clinical %>%
  filter(grepl(pulse_rate_codes, code)) %>%
  filter(value1 != "" | value2 != "" | value3 != "")  %>%
  mutate(value1_num = as.numeric(value1), value2_num = as.numeric(value2), value3_num = as.numeric(value3)) %>%
  mutate(value = coalesce(value1_num, value2_num, value3_num)) %>%
  filter(value > 30 & value < 200) %>%
  filter(!(value3 %in% c("%", "mmHG"))) %>%
  group_by(eid) %>%
  mutate(mi = n()) %>%
  mutate(mean = mean(value)) %>%
  mutate(sd = ifelse(mi == 1, 0, 
                     ifelse(sd(value) > 10, 10, sd(value)))) %>%
  mutate(lower_tail = ifelse(mi == 1, mean, mean - 4 * sd),
         higher_tail = ifelse(mi == 1, mean, mean + 4 * sd)) %>%
  filter(value >= lower_tail & value <= higher_tail) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  mutate(value_new = ifelse(rep == 1, value, 
                            ifelse(max(value) - min(value) <= 20, 
                                   mean(value), value[which.min(abs(value - mean))]))) %>%
  mutate(value_new = round(value_new, 0)) %>%
  select(eid, event_dt, data_provider, value_new) %>%
  dplyr::rename(pr = value_new) %>%
  distinct(eid, event_dt, pr, .keep_all = TRUE)

### write the pre-cleaned pulse rate into file.
data.table::fwrite(pr %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/pulse_rate.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)


###### merge the gp_clinical table and covariates.
library(dplyr)
library(tidyr)
library(lubridate)

### read in the gp_clinical table and covariates.
pr_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/pulse_rate.txt")

cova_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/assessment_center/covariate.txt")

bmi_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/assessment_center/bmi_from_UKB.txt")

### combine the gp_clinical and covariates data.
pr_joined = inner_join(pr_info, cova_info, by = c('eid' = 'Wenjian')) %>%
  left_join(bmi_info %>% select(feid, f21001_0_0) %>% drop_na %>% rename(eid = feid, BMI = f21001_0_0)) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18)

### write the final-cleaned pulse rate into file.
data.table::fwrite(pr_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/pulse_rate.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

