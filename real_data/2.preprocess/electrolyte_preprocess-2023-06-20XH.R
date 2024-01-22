
library(dplyr)
library(tidyr)

### read in the blood sodium, potassium, calcium, chloride, phosphate, bicarbonate, ferritin, iron, lithium phecode.
source("/gdata01/user/xuhe/family_relatedness/real_data-2022-12-29/1.extract_pheno/phecode-2023-05-10XH.R")

### load the cleaned gp_clinical table
gp_clinical = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/cleaned_gp_clinical.txt")

### extract blood sodium.
sodium = gp_clinical %>%
  filter(grepl(sodium_codes, code)) %>%
  filter(value1 != "" | value2 != "" | value3 != "")  %>%
  mutate(value1_num = as.numeric(value1), value2_num = as.numeric(value2), value3_num = as.numeric(value3)) %>%
  mutate(value = coalesce(value1_num, value2_num, value3_num)) %>%
  filter(value > 100 & value < 200) %>%
  mutate(value3 = toupper(value3)) %>%
  filter(value3 %in% c("", "MEA000", "MEA096", "MMOL/L", "UNKNOWN")) %>%
  group_by(eid) %>%
  mutate(mi = n()) %>%
  mutate(mean = mean(value)) %>%
  mutate(sd = ifelse(mi == 1, 0, 
                     ifelse(sd(value) > 5, 5, sd(value)))) %>%
  mutate(lower_tail = ifelse(mi == 1, mean, mean - 4 * sd),
         higher_tail = ifelse(mi == 1, mean, mean + 4 * sd)) %>%
  filter(value >= lower_tail & value <= higher_tail) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  mutate(value_new = ifelse(rep == 1, value, 
                            ifelse(max(value) - min(value) <= 5, 
                                   mean(value), value[which.min(abs(value - mean))]))) %>%
  mutate(value_new = round(value_new, 0)) %>%
  select(eid, event_dt, data_provider, value_new) %>%
  dplyr::rename(na = value_new) %>%
  distinct(eid, event_dt, na, .keep_all = TRUE)

### write the pre-cleaned blood sodium into file.
data.table::fwrite(sodium %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/sodium.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### extract blood potassium.
potassium = gp_clinical %>%
  filter(grepl(potassium_codes, code)) %>%
  filter(value1 != "" | value2 != "" | value3 != "")  %>%
  mutate(value1_num = as.numeric(value1), value2_num = as.numeric(value2), value3_num = as.numeric(value3)) %>%
  mutate(value = coalesce(value1_num, value2_num, value3_num)) %>%
  filter(value > 1.5 & value < 15) %>%
  mutate(value3 = toupper(value3)) %>%
  filter(value3 %in% c("", "MEA000", "MEA096", "MMOL/L", "UNKNOWN")) %>%
  group_by(eid) %>%
  mutate(mi = n()) %>%
  mutate(mean = mean(value)) %>%
  mutate(sd = ifelse(mi == 1, 0, 
                     ifelse(sd(value) > 1, 1, sd(value)))) %>%
  mutate(lower_tail = ifelse(mi == 1, mean, mean - 4 * sd),
         higher_tail = ifelse(mi == 1, mean, mean + 4 * sd)) %>%
  filter(value >= lower_tail & value <= higher_tail) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  mutate(value_new = ifelse(rep == 1, value, 
                            ifelse(max(value) - min(value) <= 1, 
                                   mean(value), value[which.min(abs(value - mean))]))) %>%
  mutate(value_new = round(value_new, 1)) %>%
  select(eid, event_dt, data_provider, value_new) %>%
  dplyr::rename(k = value_new) %>%
  distinct(eid, event_dt, k, .keep_all = TRUE)

### write the pre-cleaned blood potassium into file.
data.table::fwrite(potassium %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/potassium.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### extract blood calcium.
calcium = gp_clinical %>%
  filter(grepl(calcium_codes, code)) %>%
  mutate(adjcalcium = ifelse(grepl("adjusted|corrected", term_description, ignore.case = TRUE), TRUE, FALSE)) %>%
  filter(value1 != "" | value2 != "" | value3 != "")  %>%
  mutate(value1_num = as.numeric(value1), value2_num = as.numeric(value2), value3_num = as.numeric(value3)) %>%
  mutate(value = coalesce(value1_num, value2_num, value3_num)) %>%
  filter(value > 1 & value < 5) %>%
  distinct(eid, event_dt, adjcalcium, value, .keep_all = TRUE) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n(), adjnum = sum(adjcalcium)) %>%
  filter(rep == 1 | (rep == 2 & adjnum == 1))

calcium_total = calcium %>% 
  filter(adjcalcium == FALSE) %>%
  select(eid, event_dt, data_provider, value) %>%
  rename(ca_tol = value) %>%
  distinct(eid, event_dt, ca_tol, .keep_all = TRUE)

### write the pre-cleaned blood total calcium into file.
data.table::fwrite(calcium_total %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/calcium_total.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

calcium_adj = calcium %>%
  filter(adjcalcium == TRUE) %>%
  select(eid, event_dt, data_provider, value) %>%
  rename(ca_adj = value) %>%
  distinct(eid, event_dt, ca_adj, .keep_all = TRUE)

### write the pre-cleaned blood adjusted calcium into file.
data.table::fwrite(calcium_adj %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/calcium_adj.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### extract blood chloride.
chloride = gp_clinical %>%
  filter(grepl(chloride_codes, code)) %>%
  filter(value1 != "" | value2 != "" | value3 != "") %>%
  mutate(value1_num = as.numeric(value1), value2_num = as.numeric(value2), value3_num = as.numeric(value3)) %>%
  mutate(value = coalesce(value1_num, value2_num, value3_num)) %>%
  filter(value > 80 & value < 120) %>%
  mutate(value3 = toupper(value3)) %>%
  filter(value3 != "MMOL") %>%
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
                            ifelse(max(value) - min(value) <= 4, 
                                   mean(value), value[which.min(abs(value - mean))]))) %>%
  mutate(value_new = round(value_new, 0)) %>%
  select(eid, event_dt, data_provider, value_new) %>%
  dplyr::rename(cl = value_new) %>%
  distinct(eid, event_dt, cl, .keep_all = TRUE)

### write the pre-cleaned blood chloride into file.
data.table::fwrite(chloride %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/chloride.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### extract blood inorganic phosphate.
phosphate = gp_clinical %>%
  filter(grepl(phosphate_codes, code)) %>%
  filter(value1 != "" | value2 != "" | value3 != "")  %>%
  mutate(value1_num = as.numeric(value1), value2_num = as.numeric(value2), value3_num = as.numeric(value3)) %>%
  mutate(value = coalesce(value1_num, value2_num, value3_num)) %>%
  filter(value > 0.3 & value < 3) %>%
  group_by(eid) %>%
  mutate(mi = n()) %>%
  mutate(mean = mean(value)) %>%
  mutate(sd = ifelse(mi == 1, 0, 
                     ifelse(sd(value) > 0.2, 0.2, sd(value)))) %>%
  mutate(lower_tail = ifelse(mi == 1, mean, mean - 4 * sd),
         higher_tail = ifelse(mi == 1, mean, mean + 4 * sd)) %>%
  filter(value >= lower_tail & value <= higher_tail) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  mutate(value_new = ifelse(rep == 1, value, 
                            ifelse(max(value) - min(value) <= 0.2, 
                                   mean(value), value[which.min(abs(value - mean))]))) %>%
  mutate(value_new = round(value_new, 2)) %>%
  select(eid, event_dt, data_provider, value_new) %>%
  dplyr::rename(pi = value_new) %>%
  distinct(eid, event_dt, pi, .keep_all = TRUE)

### write the pre-cleaned blood phosphate into file.
data.table::fwrite(phosphate %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/phosphate.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### extract blood bicarbonate.
bicarbonate = gp_clinical %>%
  filter(grepl(bicarbonate_codes, code)) %>%
  filter(value1 != "" | value2 != "" | value3 != "")  %>%
  mutate(value1_num = as.numeric(value1), value2_num = as.numeric(value2), value3_num = as.numeric(value3)) %>%
  mutate(value = coalesce(value1_num, value2_num, value3_num)) %>%
  filter(value > 10 & value < 40) %>%
  mutate(value3 = toupper(value3)) %>%
  filter(value3 != "MEA026" & value3 != "MEA079") %>%
  group_by(eid) %>%
  mutate(mi = n()) %>%
  mutate(mean = mean(value)) %>%
  mutate(sd = ifelse(mi == 1, 0, 
                     ifelse(sd(value) > 3, 3, sd(value)))) %>%
  mutate(lower_tail = ifelse(mi == 1, mean, mean - 4 * sd),
         higher_tail = ifelse(mi == 1, mean, mean + 4 * sd)) %>%
  filter(value >= lower_tail & value <= higher_tail) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  mutate(value_new = ifelse(rep == 1, value, 
                            ifelse(max(value) - min(value) <= 3, 
                                   mean(value), value[which.min(abs(value - mean))]))) %>%
  mutate(value_new = round(value_new, 0)) %>%
  select(eid, event_dt, data_provider, value_new) %>%
  dplyr::rename(bi = value_new) %>%
  distinct(eid, event_dt, bi, .keep_all = TRUE)

### write the pre-cleaned blood bicarbonate into file.
data.table::fwrite(bicarbonate %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/bicarbonate.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### extract blood ferritin.
ferritin = gp_clinical %>%
  filter(grepl(ferritin_codes, code)) %>%
  filter(value1 != "" | value2 != "" | value3 != "")  %>%
  mutate(value1_num = as.numeric(value1), value2_num = as.numeric(value2), value3_num = as.numeric(value3)) %>%
  mutate(value = coalesce(value1_num, value2_num, value3_num)) %>%
  filter(value > 1 & value < 3000) %>%
  mutate(value3 = toupper(value3)) %>%
  filter(value3 %in% c("", "MEA000", "MEA106", "MEA133", "NG/ML", "UG/L", "UNKNOWN")) %>%
  group_by(eid) %>%
  mutate(mi = n()) %>%
  mutate(mean = mean(value)) %>%
  mutate(sd = ifelse(mi == 1, 0, 
                     ifelse(sd(value) > 100, 100, sd(value)))) %>%
  mutate(lower_tail = ifelse(mi == 1, mean, mean - 4 * sd),
         higher_tail = ifelse(mi == 1, mean, mean + 4 * sd)) %>%
  filter(value >= lower_tail & value <= higher_tail) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  filter(rep <= 2) %>%
  mutate(value_new = ifelse(rep == 1, value, 
                            ifelse(max(value) - min(value) <= 5, 
                                   mean(value), value[which.min(abs(value - mean))]))) %>%
  mutate(value_new = round(value_new, 1)) %>%
  select(eid, event_dt, data_provider, value_new) %>%
  dplyr::rename(sf = value_new) %>%
  distinct(eid, event_dt, sf, .keep_all = TRUE)

### write the pre-cleaned blood ferritin into file.
data.table::fwrite(ferritin %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/ferritin.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### extract blood iron.
iron = gp_clinical %>%
  filter(grepl(iron_codes, code)) %>%
  filter(value1 != "" | value2 != "" | value3 != "")  %>%
  mutate(value1_num = as.numeric(value1), value2_num = as.numeric(value2), value3_num = as.numeric(value3)) %>%
  mutate(value = coalesce(value1_num, value2_num, value3_num)) %>%
  filter(value > 1 & value < 50) %>%
  mutate(value3 = toupper(value3)) %>%
  filter(value3 %in% c("", "MEA000", "MEA142", "UMOL/L")) %>%
  group_by(eid) %>%
  mutate(mi = n()) %>%
  mutate(mean = mean(value)) %>%
  mutate(sd = ifelse(mi == 1, 0, 
                     ifelse(sd(value) > 5, 5, sd(value)))) %>%
  mutate(lower_tail = ifelse(mi == 1, mean, mean - 4 * sd),
         higher_tail = ifelse(mi == 1, mean, mean + 4 * sd)) %>%
  filter(value >= lower_tail & value <= higher_tail) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  mutate(value_new = ifelse(rep == 1, value, 
                            ifelse(max(value) - min(value) <= 5, 
                                   mean(value), value[which.min(abs(value - mean))]))) %>%
  mutate(value_new = round(value_new, 1)) %>%
  select(eid, event_dt, data_provider, value_new) %>%
  dplyr::rename(fe = value_new) %>%
  distinct(eid, event_dt, fe, .keep_all = TRUE)

### write the pre-cleaned blood iron into file.
data.table::fwrite(iron %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/iron.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)


###### merge the gp_clinical table and covariates.
library(dplyr)
library(tidyr)
library(lubridate)

### read in the gp_clinical table and covariates.
sodium_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/sodium.txt")

potassium_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/potassium.txt")

calcium_total_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/calcium_total.txt")

calcium_adj_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/calcium_adj.txt")

chloride_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/chloride.txt")

phosphate_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/phosphate.txt")

bicarbonate_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/bicarbonate.txt")

ferritin_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/ferritin.txt")

iron_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/iron.txt")

cova_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/assessment_center/covariate.txt")

bmi_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/assessment_center/bmi_from_UKB.txt")

### combine the gp_clinical and covariates data.
sodium_joined = sodium_info %>%
  inner_join(cova_info, by = c('eid' = 'Wenjian')) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18) %>%
  left_join(bmi_info %>% select(feid, f21001_0_0) %>% rename(eid = feid, BMI = f21001_0_0))

### write the final-cleaned sodium into file.
data.table::fwrite(sodium_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/sodium.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### combine the gp_clinical and covariates data.
potassium_joined = potassium_info %>%
  inner_join(cova_info, by = c('eid' = 'Wenjian')) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18) %>%
  left_join(bmi_info %>% select(feid, f21001_0_0) %>% rename(eid = feid, BMI = f21001_0_0))

### write the final-cleaned potassium into file.
data.table::fwrite(potassium_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/potassium.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### combine the gp_clinical and covariates data.
calcium_total_joined = calcium_total_info %>%
  inner_join(cova_info, by = c('eid' = 'Wenjian')) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18) %>%
  left_join(bmi_info %>% select(feid, f21001_0_0) %>% rename(eid = feid, BMI = f21001_0_0))

### write the final-cleaned calcium_total into file.
data.table::fwrite(calcium_total_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/calcium_total.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### combine the gp_clinical and covariates data.
calcium_adj_joined = calcium_adj_info %>%
  inner_join(cova_info, by = c('eid' = 'Wenjian')) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18) %>%
  left_join(bmi_info %>% select(feid, f21001_0_0) %>% rename(eid = feid, BMI = f21001_0_0))

### write the final-cleaned calcium_adj into file.
data.table::fwrite(calcium_adj_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/calcium_adj.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### combine the gp_clinical and covariates data.
chloride_joined = chloride_info %>%
  inner_join(cova_info, by = c('eid' = 'Wenjian')) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18) %>%
  left_join(bmi_info %>% select(feid, f21001_0_0) %>% rename(eid = feid, BMI = f21001_0_0))

### write the final-cleaned chloride into file.
data.table::fwrite(chloride_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/chloride.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### combine the gp_clinical and covariates data.
phosphate_joined = phosphate_info %>%
  inner_join(cova_info, by = c('eid' = 'Wenjian')) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18) %>%
  left_join(bmi_info %>% select(feid, f21001_0_0) %>% rename(eid = feid, BMI = f21001_0_0))

### write the final-cleaned phosphate into file.
data.table::fwrite(phosphate_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/phosphate.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### combine the gp_clinical and covariates data.
bicarbonate_joined = bicarbonate_info %>%
  inner_join(cova_info, by = c('eid' = 'Wenjian')) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18) %>%
  left_join(bmi_info %>% select(feid, f21001_0_0) %>% rename(eid = feid, BMI = f21001_0_0))

### write the final-cleaned bicarbonate into file.
data.table::fwrite(bicarbonate_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/bicarbonate.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### combine the gp_clinical and covariates data.
ferritin_joined = ferritin_info %>%
  inner_join(cova_info, by = c('eid' = 'Wenjian')) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18) %>%
  left_join(bmi_info %>% select(feid, f21001_0_0) %>% rename(eid = feid, BMI = f21001_0_0))

### write the final-cleaned ferritin into file.
data.table::fwrite(ferritin_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/ferritin.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### combine the gp_clinical and covariates data.
iron_joined = iron_info %>%
  inner_join(cova_info, by = c('eid' = 'Wenjian')) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18) %>%
  left_join(bmi_info %>% select(feid, f21001_0_0) %>% rename(eid = feid, BMI = f21001_0_0))

### write the final-cleaned iron into file.
data.table::fwrite(iron_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/iron.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)
