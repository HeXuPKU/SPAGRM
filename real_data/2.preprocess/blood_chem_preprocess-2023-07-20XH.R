
###### pre-clean blood protein phenotype gp_clinical table.
library(dplyr)
library(tidyr)

### read in the blood bilirubin, urea, urate, folate, vitamin B12, VD phecode.
source("/gdata01/user/xuhe/family_relatedness/real_data-2022-12-29/1.extract_pheno/phecode-2023-05-10XH.R")

### load the cleaned gp_clinical table.
gp_clinical = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/cleaned_gp_clinical.txt")

### extract bilirubin from gp_clinical data.
bilirubin = gp_clinical %>%
  filter(grepl(bilirubin_codes, code)) %>% 
  filter(value3 != "\xb5mol/L") %>%
  mutate(value1_numeric = as.numeric(value1), value2_numeric = as.numeric(value2), value3_numeric = as.numeric(value3)) %>%
  mutate(bili = coalesce(value1_numeric, value2_numeric, value3_numeric)) %>%
  filter(bili > 1 & bili < 100) %>%
  mutate(value3 = toupper(value3)) %>%
  filter(value3 %in% c("", "MEA000", "MEA142", "UMOL/L")) %>%
  group_by(eid) %>%
  mutate(mi = n()) %>%
  mutate(mean = mean(bili)) %>%
  mutate(sd = ifelse(mi == 1, 0, 
                     ifelse(sd(bili) > 4.5, 4.5, sd(bili)))) %>%
  mutate(lower_tail = ifelse(mi == 1, mean, mean - 4 * sd),
         higher_tail = ifelse(mi == 1, mean, mean + 4 * sd)) %>%
  filter(bili >= lower_tail & bili <= higher_tail) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  mutate(bili_new = ifelse(rep == 1, bili, 
                           ifelse(max(bili) - min(bili) <= 2, 
                                  mean(bili), bili[which.min(abs(bili - mean))]))) %>%
  mutate(bili_new = round(bili_new, 0)) %>%
  select(eid, data_provider, event_dt, bili_new) %>% 
  dplyr::rename(bili = bili_new) %>%
  distinct(eid, event_dt, .keep_all = TRUE)

### write the pre-cleaned bilirubin into file.
data.table::fwrite(bilirubin %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/bilirubin.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### extract urea from gp_clinical data.
urea = gp_clinical %>%
  filter(grepl(urea_codes, code)) %>% 
  mutate(value1_numeric = as.numeric(value1), value2_numeric = as.numeric(value2), value3_numeric = as.numeric(value3)) %>%
  mutate(urea = coalesce(value1_numeric, value2_numeric, value3_numeric)) %>%
  filter(urea > 0.8 & urea < 40) %>%
  mutate(value3 = toupper(value3)) %>%
  filter(value3 %in% c("", "MEA000", "MEA096", "MMOL/L")) %>%
  group_by(eid) %>%
  mutate(mi = n()) %>%
  mutate(mean = mean(urea)) %>%
  mutate(sd = ifelse(mi == 1, 0, 
                     ifelse(sd(urea) > 4, 4, sd(urea)))) %>%
  mutate(lower_tail = ifelse(mi == 1, mean, mean - 4 * sd),
         higher_tail = ifelse(mi == 1, mean, mean + 4 * sd)) %>%
  filter(urea >= lower_tail & urea <= higher_tail) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  mutate(urea_new = ifelse(rep == 1, urea, 
                           ifelse(max(urea) - min(urea) <= 1, 
                                  mean(urea), urea[which.min(abs(urea - mean))]))) %>%
  mutate(urea_new = round(urea_new, 1)) %>%
  select(eid, data_provider, event_dt, urea_new) %>% 
  dplyr::rename(urea = urea_new) %>%
  distinct(eid, event_dt, .keep_all = TRUE)

### write the pre-cleaned urea into file.
data.table::fwrite(urea %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/urea.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### extract urate from gp_clinical data.
urate = gp_clinical %>%
  filter(grepl(urate_codes, code)) %>% 
  mutate(value1_numeric = as.numeric(value1), value2_numeric = as.numeric(value2), value3_numeric = as.numeric(value3)) %>%
  mutate(urate = coalesce(value1_numeric, value2_numeric, value3_numeric)) %>%
  filter((urate > 80 & urate < 1000) | (urate > 0.08 & urate < 1)) %>%
  mutate(value3 = toupper(value3)) %>%
  filter(!(value3 %in% c("MEA061", "MEA095", "MEA099", "MEA156", "MG/MMOL"))) %>%
  mutate(urate = ifelse(urate < 1, 1000 * urate, urate)) %>%
  group_by(eid) %>%
  mutate(mi = n()) %>%
  mutate(mean = mean(urate)) %>%
  mutate(sd = ifelse(mi == 1, 0, 
                     ifelse(sd(urate) > 90, 90, sd(urate)))) %>%
  mutate(lower_tail = ifelse(mi == 1, mean, mean - 4 * sd),
         higher_tail = ifelse(mi == 1, mean, mean + 4 * sd)) %>%
  filter(urate >= lower_tail & urate <= higher_tail) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  mutate(urate_new = ifelse(rep == 1, urate, 
                           ifelse(max(urate) - min(urate) <= 10, 
                                  mean(urate), urate[which.min(abs(urate - mean))]))) %>%
  mutate(urate_new = round(urate_new, 0)) %>%
  select(eid, data_provider, event_dt, urate_new) %>% 
  dplyr::rename(urate = urate_new) %>%
  distinct(eid, event_dt, .keep_all = TRUE)

### write the pre-cleaned urate into file.
data.table::fwrite(urate %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/urate.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### extract folate from gp_clinical data.
folate = gp_clinical %>%
  filter(grepl(folate_codes, code)) %>% 
  mutate(value1_numeric = as.numeric(value1), value2_numeric = as.numeric(value2), value3_numeric = as.numeric(value3)) %>%
  mutate(fo = coalesce(value1_numeric, value2_numeric, value3_numeric)) %>%
  filter(fo > 0.8 & fo < 30) %>%
  mutate(value3 = toupper(value3)) %>%
  filter(value3 %in% c("", "MEA000", "MEA106", "MEA133", "UG/L")) %>%
  group_by(eid) %>%
  mutate(mi = n()) %>%
  mutate(mean = mean(fo)) %>%
  mutate(sd = ifelse(mi == 1, 0, 
                     ifelse(sd(fo) > 4, 4, sd(fo)))) %>%
  mutate(lower_tail = ifelse(mi == 1, mean, mean - 4 * sd),
         higher_tail = ifelse(mi == 1, mean, mean + 4 * sd)) %>%
  filter(fo >= lower_tail & fo <= higher_tail) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  mutate(fo_new = ifelse(rep == 1, fo, 
                         ifelse(max(fo) - min(fo) <= 1, mean(fo), 0))) %>%
  filter(fo_new > 0) %>%
  mutate(fo_new = round(fo_new, 1)) %>%
  select(eid, data_provider, event_dt, fo_new) %>% 
  dplyr::rename(fo = fo_new) %>%
  distinct(eid, event_dt, .keep_all = TRUE)

### write the pre-cleaned folate into file.
data.table::fwrite(folate %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/folate.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### extract vitamin B12 from gp_clinical data.
B12 = gp_clinical %>%
  filter(grepl(B12_codes, code)) %>% 
  mutate(value1_numeric = as.numeric(value1), value2_numeric = as.numeric(value2), value3_numeric = as.numeric(value3)) %>%
  mutate(b12 = coalesce(value1_numeric, value2_numeric, value3_numeric)) %>%
  filter(b12 > 50 & b12 < 1500) %>%
  mutate(value3 = toupper(value3)) %>%
  filter(value3 %in% c("", "MEA000", "MEA105", "MEA116", "NG/L", "PG/ML")) %>%
  group_by(eid) %>%
  mutate(mi = n()) %>%
  mutate(mean = mean(b12)) %>%
  mutate(sd = ifelse(mi == 1, 0, 
                     ifelse(sd(b12) > 200, 200, sd(b12)))) %>%
  mutate(lower_tail = ifelse(mi == 1, mean, mean - 4 * sd),
         higher_tail = ifelse(mi == 1, mean, mean + 4 * sd)) %>%
  filter(b12 >= lower_tail & b12 <= higher_tail) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  mutate(b12_new = ifelse(rep == 1, b12, 
                         ifelse(max(b12) - min(b12) <= 50, mean(b12), 0))) %>%
  filter(b12_new > 0) %>%
  mutate(b12_new = round(b12_new, 0)) %>%
  select(eid, data_provider, event_dt, b12_new) %>% 
  dplyr::rename(b12 = b12_new) %>%
  distinct(eid, event_dt, .keep_all = TRUE)

### write the pre-cleaned vitamin B12 into file.
data.table::fwrite(B12 %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/B12.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### extract vitamin B12 from gp_clinical data.
VD = gp_clinical %>%
  filter(grepl(VD_codes, code)) %>% 
  mutate(value1_numeric = as.numeric(value1), value2_numeric = as.numeric(value2), value3_numeric = as.numeric(value3)) %>%
  mutate(vd = coalesce(value1_numeric, value2_numeric, value3_numeric)) %>%
  filter(vd > 10 & vd < 200) %>%
  mutate(value3 = toupper(value3)) %>%
  filter(value3 %in% c("", "MEA110", "NMOL/L")) %>%
  group_by(eid) %>%
  mutate(mi = n()) %>%
  mutate(mean = mean(vd)) %>%
  mutate(sd = ifelse(mi == 1, 0, 
                     ifelse(sd(vd) > 20, 20, sd(vd)))) %>%
  mutate(lower_tail = ifelse(mi == 1, mean, mean - 4 * sd),
         higher_tail = ifelse(mi == 1, mean, mean + 4 * sd)) %>%
  filter(vd >= lower_tail & vd <= higher_tail) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  mutate(vd_new = ifelse(rep == 1, vd, 
                          ifelse(max(vd) - min(vd) <= 2, mean(vd), 0))) %>%
  filter(vd_new > 0) %>%
  mutate(vd_new = round(vd_new, 0)) %>%
  select(eid, data_provider, event_dt, vd_new) %>% 
  dplyr::rename(vd = vd_new) %>%
  distinct(eid, event_dt, .keep_all = TRUE)

### write the pre-cleaned vitamin D into file.
data.table::fwrite(VD %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/VD.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)


###### merge the gp_clinical table and covariates.
library(dplyr)
library(tidyr)
library(lubridate)

### read in the gp_clinical table and covariates.
bilirubin_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/bilirubin.txt")

urea_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/urea.txt")

urate_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/urate.txt")

folate_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/folate.txt")

B12_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/B12.txt")

VD_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/VD.txt")

cova_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/assessment_center/covariate.txt")

bmi_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/assessment_center/bmi_from_UKB.txt")

### combine the gp_clinical and covariates data.
bilirubin_joined = bilirubin_info %>%
  inner_join(cova_info, by = c('eid' = 'Wenjian')) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18) %>%
  left_join(bmi_info %>% select(feid, f21001_0_0) %>% rename(eid = feid, BMI = f21001_0_0))

### write the final-cleaned bilirubin into file.
data.table::fwrite(bilirubin_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/bilirubin.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### combine the gp_clinical and covariates data.
urea_joined = urea_info %>%
  inner_join(cova_info, by = c('eid' = 'Wenjian')) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18) %>%
  left_join(bmi_info %>% select(feid, f21001_0_0) %>% rename(eid = feid, BMI = f21001_0_0))

### write the final-cleaned urea into file.
data.table::fwrite(urea_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/urea.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### combine the gp_clinical and covariates data.
urate_joined = urate_info %>%
  inner_join(cova_info, by = c('eid' = 'Wenjian')) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18) %>%
  left_join(bmi_info %>% select(feid, f21001_0_0) %>% rename(eid = feid, BMI = f21001_0_0))

### write the final-cleaned urate into file.
data.table::fwrite(urate_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/urate.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### combine the gp_clinical and covariates data.
folate_joined = folate_info %>%
  inner_join(cova_info, by = c('eid' = 'Wenjian')) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18) %>%
  left_join(bmi_info %>% select(feid, f21001_0_0) %>% rename(eid = feid, BMI = f21001_0_0))

### write the final-cleaned folate into file.
data.table::fwrite(folate_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/folate.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### combine the gp_clinical and covariates data.
B12_joined = B12_info %>%
  inner_join(cova_info, by = c('eid' = 'Wenjian')) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18) %>%
  left_join(bmi_info %>% select(feid, f21001_0_0) %>% rename(eid = feid, BMI = f21001_0_0))

### write the final-cleaned VB12 into file.
data.table::fwrite(B12_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/B12.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### combine the gp_clinical and covariates data.
VD_joined = VD_info %>%
  inner_join(cova_info, by = c('eid' = 'Wenjian')) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18) %>%
  left_join(bmi_info %>% select(feid, f21001_0_0) %>% rename(eid = feid, BMI = f21001_0_0))

### write the final-cleaned VD into file.
data.table::fwrite(VD_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/VD.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)
