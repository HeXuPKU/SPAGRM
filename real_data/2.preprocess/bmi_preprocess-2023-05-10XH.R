
###### pre-clean the bmi gp_clinical table.
library(dplyr)
library(tidyr)

### read in the bmi phecode.
source("/gdata01/user/xuhe/family_relatedness/real_data-2022-12-29/1.extract_pheno/phecode-2023-05-10XH.R")

### load the cleaned gp_clinical table.
gp_clinical = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/cleaned_gp_clinical.txt")

### extract height, weight and BMI.
hwbmi = gp_clinical %>%
  filter(grepl(height_weight_BMI_codes, code)) %>%
  mutate(value = coalesce(as.numeric(value1), as.numeric(value2), as.numeric(value3))) %>%
  filter(!is.na(value) & value > 0) %>%
  mutate(trait = ifelse(grepl("BMI|Body Mass Index", term_description, ignore.case=T), "BMI",
                        ifelse(grepl("Height", term_description, ignore.case=T), "Height",
                               "Weight"))) %>%
  mutate(value = ifelse(trait == "Height" & value > 2.1, value/100, value)) %>% # cm to meters
  filter((trait == "Weight" & value < 200 & value > 30) |
           (trait == "Height" & value < 2.1 & value > 1.25 )|
           (trait == "BMI" & value < 75 & value > 12))

### display the data.
hwbmi %>% group_by(code, term_description, trait) %>% 
  summarize(n = n(), mean=round(mean(value), 1)) %>% arrange(trait, desc(n))

### separate the phenotypes.
height = hwbmi %>% 
  filter(trait == "Height") %>% 
  arrange(eid, event_dt)

weight = hwbmi %>% 
  filter(trait == "Weight") %>% 
  arrange(eid, event_dt)

bmi = hwbmi %>% 
  filter(trait == "BMI")

bmi_extra = weight %>% 
  filter(data_provider == 2) %>% 
  mutate(bmi_extra = as.numeric(value3)) %>% 
  filter(!is.na(bmi_extra)) %>%
  filter(bmi_extra > 12 & bmi_extra < 75) %>%
  select(-value) %>%
  dplyr::rename(value = bmi_extra)

bmi_combine = full_join(bmi, bmi_extra) %>%
  arrange(eid, event_dt)

### pre-clean the height.
height_pre = height %>%
  group_by(eid) %>%
  mutate(mi = n()) %>%
  mutate(mean = mean(value)) %>%
  mutate(sd = ifelse(mi == 1, 0, 
                     ifelse(sd(value) > 0.025, 0.025, sd(value)))) %>%
  mutate(lower_tail = ifelse(mi == 1, mean, mean - 4 * sd),
         higher_tail = ifelse(mi == 1, mean, mean + 4 * sd)) %>%
  filter(value >= lower_tail & value <= higher_tail) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  mutate(value_new = ifelse(rep == 1, value, 
                            ifelse(max(value) - min(value) <= 0.05, 
                                   mean(value), value[which.min(abs(value - mean))]))) %>%
  mutate(value_new = round(value_new, 2)) %>%
  select(eid, event_dt, value_new, code, term_description) %>%
  dplyr::rename(height = value_new, height_code = code, height_term_description = term_description) %>%
  distinct(eid, event_dt, height, .keep_all = TRUE)

### pre-clean the weight.
weight_pre = weight %>% 
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
                            ifelse(max(value) - min(value) <= 1, 
                                   mean(value), value[which.min(abs(value - mean))]))) %>%
  mutate(value_new = round(value_new, 2)) %>%
  select(eid, event_dt, value_new, code, term_description) %>%
  dplyr::rename(weight = value_new, weight_code = code, weight_term_description = term_description) %>%
  distinct(eid, event_dt, weight, .keep_all = TRUE)

### pre-clean the BMI.
bmi_pre = bmi_combine %>% 
  group_by(eid) %>%
  mutate(mi = n()) %>%
  mutate(mean = mean(value)) %>%
  mutate(sd = ifelse(mi == 1, 0, 
                     ifelse(sd(value) > 4, 4, sd(value)))) %>%
  mutate(lower_tail = ifelse(mi == 1, mean, mean - 4 * sd),
         higher_tail = ifelse(mi == 1, mean, mean + 4 * sd)) %>%
  filter(value >= lower_tail & value <= higher_tail) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  mutate(value_new = ifelse(rep == 1, value, 
                            ifelse(max(value) - min(value) <= 2, 
                                   mean(value), value[which.min(abs(value - mean))]))) %>%
  mutate(value_new = round(value_new, 2)) %>%
  select(eid, event_dt, value_new, code, term_description) %>%
  dplyr::rename(bmi = value_new, bmi_code = code, bmi_term_description = term_description) %>%
  distinct(eid, event_dt, bmi, .keep_all = TRUE)

### join the phenotypes.
joined_bmi = full_join(height_pre, weight_pre) %>%
  full_join(bmi_pre)

### clean the combined bmi.
cleaned_bmi = joined_bmi %>%
  arrange(eid, event_dt) %>%
  group_by(eid) %>%
  fill(height, .direction = "downup") %>%
  ungroup() %>%
  mutate(bmi_observed = round(bmi, 1),
         bmi_calculated = round(weight/(height^2), 1),
         bmi_diff = abs(bmi_observed - bmi_calculated)) %>%
  filter(!(!is.na(bmi_diff) & bmi_diff > 2)) %>%
  mutate(bmi_final = coalesce(bmi_calculated, bmi_observed)) %>%
  mutate(height_carried = ifelse(!is.na(height), height, sqrt(weight/bmi_final))) %>%
  filter(height_carried < 2.1 & height_carried > 1.25) %>%
  group_by(eid) %>%
  mutate(mi = n()) %>%
  mutate(mean = mean(height_carried)) %>%
  mutate(sd = ifelse(mi == 1, 0, 
                     ifelse(sd(height_carried) > 0.025, 0.025, sd(height_carried)))) %>%
  mutate(lower_tail = ifelse(mi == 1, mean, mean - 4 * sd),
         higher_tail = ifelse(mi == 1, mean, mean + 4 * sd)) %>%
  ungroup() %>%
  filter(!is.na(height) | 
           (is.na(height) & height_carried >= lower_tail & height_carried <= higher_tail)) %>%
  filter(bmi_final < 75 & bmi_final > 12) %>%
  select(eid, event_dt, height_carried, weight, bmi_final) %>%
  dplyr::rename(height = height_carried, bmi = bmi_final) %>%
  distinct(eid, event_dt, bmi, .keep_all = TRUE)

### write the pre-cleaned bmi into file.
data.table::fwrite(cleaned_bmi %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/bmi_height_weight.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)


###### merge the gp_clinical table and covariates.
library(dplyr)
library(tidyr)
library(lubridate)

### read in the gp_clinical table and covariates.
bmi_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/bmi_height_weight.txt")

cova_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/assessment_center/covariate.txt")

### combine the gp_clinical and covariates data.
bmi_joined = inner_join(bmi_info, cova_info, by = c('eid' = 'Wenjian')) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18)

### write the final-cleaned bmi into file.
data.table::fwrite(bmi_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/bmi.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)
