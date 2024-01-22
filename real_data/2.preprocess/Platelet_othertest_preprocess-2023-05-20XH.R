
###### pre-clean blood protein phenotype gp_clinical table.
library(dplyr)
library(tidyr)

### read in the MPV, PV, PT, INR phecode.
source("/gdata01/user/xuhe/family_relatedness/real_data-2022-12-29/1.extract_pheno/phecode-2023-05-10XH.R")

### load the cleaned gp_clinical table.
gp_clinical = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/cleaned_gp_clinical.txt")

### extract Mean platelet volume from gp_clinical data.
MPV = gp_clinical %>%
  filter(grepl(MPV_codes, code)) %>% 
  mutate(value1_numeric = as.numeric(value1), value2_numeric = as.numeric(value2), value3_numeric = as.numeric(value3)) %>%
  mutate(mpv = coalesce(value1_numeric, value2_numeric, value3_numeric)) %>%
  filter(mpv > 5 & mpv < 15) %>%
  group_by(eid) %>%
  mutate(mi = n()) %>%
  mutate(mean = mean(mpv)) %>%
  mutate(sd = ifelse(mi == 1, 0, 
                     ifelse(sd(mpv) > 1.2, 1.2, sd(mpv)))) %>%
  mutate(lower_tail = ifelse(mi == 1, mean, mean - 4 * sd),
         higher_tail = ifelse(mi == 1, mean, mean + 4 * sd)) %>%
  filter(mpv >= lower_tail & mpv <= higher_tail) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  mutate(mpv_new = ifelse(rep == 1, mpv, 
                          ifelse(max(mpv) - min(mpv) <= 1.2, 
                                 mean(mpv), mpv[which.min(abs(mpv - mean))]))) %>%
  mutate(mpv_new = round(mpv_new, 1)) %>%
  select(eid, data_provider, event_dt, mpv_new) %>% 
  dplyr::rename(mpv = mpv_new) %>%
  distinct(eid, event_dt, .keep_all = TRUE)

### write the pre-cleaned Mean platelet volume into file.
data.table::fwrite(MPV %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/MPV.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### extract Plasma viscosity from gp_clinical data.
PV = gp_clinical %>%
  filter(grepl(PV_codes, code)) %>% 
  mutate(value1_numeric = as.numeric(value1), value2_numeric = as.numeric(value2), value3_numeric = as.numeric(value3)) %>%
  mutate(pv = coalesce(value1_numeric, value2_numeric, value3_numeric)) %>%
  filter(pv > 1.2 & pv < 2.5) %>%
  mutate(value3 = toupper(value3)) %>%
  filter(value3 %in% c("", "MEA158", "MEA172")) %>%
  group_by(eid) %>%
  mutate(mi = n()) %>%
  mutate(mean = mean(pv)) %>%
  mutate(sd = ifelse(mi == 1, 0, 
                     ifelse(sd(pv) > 0.12, 0.12, sd(pv)))) %>%
  mutate(lower_tail = ifelse(mi == 1, mean, mean - 4 * sd),
         higher_tail = ifelse(mi == 1, mean, mean + 4 * sd)) %>%
  filter(pv >= lower_tail & pv <= higher_tail) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  mutate(pv_new = ifelse(rep == 1, pv, 
                     ifelse(max(pv) - min(pv) <= 0.1, 
                            mean(pv), pv[which.min(abs(pv - mean))]))) %>%
  mutate(pv_new = round(pv_new, 2)) %>%
  select(eid, data_provider, event_dt, pv_new) %>% 
  dplyr::rename(pv = pv_new) %>%
  distinct(eid, event_dt, .keep_all = TRUE)

### write the pre-cleaned Plasma viscosity into file.
data.table::fwrite(PV %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/PV.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### extract International normalised ratio from gp_clinical data.
INR = gp_clinical %>%
  filter(grepl(INR_codes, code)) %>% 
  mutate(value1_numeric = as.numeric(value1), value2_numeric = as.numeric(value2), value3_numeric = as.numeric(value3)) %>%
  mutate(inr = coalesce(value1_numeric, value2_numeric, value3_numeric)) %>%
  filter(inr > 0.8 & inr < 6) %>%
  mutate(value3 = toupper(value3)) %>%
  filter(value3 %in% c("", "MEA151", "MEA161", "RATIO")) %>%
  group_by(eid) %>%
  mutate(mi = n()) %>%
  mutate(mean = mean(inr)) %>%
  mutate(sd = ifelse(mi == 1, 0, 
                     ifelse(sd(inr) > 1, 1, sd(inr)))) %>%
  mutate(lower_tail = ifelse(mi == 1, mean, mean - 4 * sd),
         higher_tail = ifelse(mi == 1, mean, mean + 4 * sd)) %>%
  filter(inr >= lower_tail & inr <= higher_tail) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  mutate(inr_new = ifelse(rep == 1, inr, 
                          ifelse(max(inr) - min(inr) <= 1, 
                                 mean(inr), inr[which.min(abs(inr - mean))]))) %>%
  mutate(inr_new = round(inr_new, 2)) %>%
  select(eid, data_provider, event_dt, inr_new) %>% 
  dplyr::rename(inr = inr_new) %>%
  distinct(eid, event_dt, .keep_all = TRUE)

### write the pre-cleaned International normalised ratio into file.
data.table::fwrite(INR %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/INR.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### extract Activated partial thromboplastin time from gp_clinical data.
APTT = gp_clinical %>%
  filter(grepl(APTT_codes, code)) %>% 
  mutate(value1_numeric = as.numeric(value1), value2_numeric = as.numeric(value2), value3_numeric = as.numeric(value3)) %>%
  mutate(aptt = coalesce(value1_numeric, value2_numeric, value3_numeric)) %>%
  filter(aptt > 15 & aptt < 60) %>%
  mutate(value3 = toupper(value3)) %>%
  group_by(eid) %>%
  mutate(mi = n()) %>%
  mutate(mean = mean(aptt)) %>%
  mutate(sd = ifelse(mi == 1, 0, 
                     ifelse(sd(aptt) > 2.5, 2.5, sd(aptt)))) %>%
  mutate(lower_tail = ifelse(mi == 1, mean, mean - 4 * sd),
         higher_tail = ifelse(mi == 1, mean, mean + 4 * sd)) %>%
  filter(aptt >= lower_tail & aptt <= higher_tail) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  mutate(aptt_new = ifelse(rep == 1, aptt, 
                          ifelse(max(aptt) - min(aptt) <= 2.5, mean(aptt), 0))) %>%
  filter(aptt_new > 0) %>%
  mutate(aptt_new = round(aptt_new, 1)) %>%
  select(eid, data_provider, event_dt, aptt_new) %>% 
  dplyr::rename(aptt = aptt_new) %>%
  distinct(eid, event_dt, .keep_all = TRUE)

### write the pre-cleaned Activated partial thromboplastin time into file.
data.table::fwrite(APTT %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/APTT.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)


###### merge the gp_clinical table and covariates.
library(dplyr)
library(tidyr)
library(lubridate)

### read in the gp_clinical table and covariates.
MPV_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/MPV.txt")

PV_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/PV.txt")

INR_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/INR.txt")

APTT_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/APTT.txt")

cova_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/assessment_center/covariate.txt")

bmi_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/assessment_center/bmi_from_UKB.txt")

### combine the gp_clinical and covariates data.
MPV_joined = MPV_info %>%
  inner_join(cova_info, by = c('eid' = 'Wenjian')) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18) %>%
  left_join(bmi_info %>% select(feid, f21001_0_0) %>% rename(eid = feid, BMI = f21001_0_0))

### write the final-cleaned Mean platelet volume into file.
data.table::fwrite(MPV_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/MPV.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### combine the gp_clinical and covariates data.
PV_joined = PV_info %>%
  inner_join(cova_info, by = c('eid' = 'Wenjian')) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18) %>%
  left_join(bmi_info %>% select(feid, f21001_0_0) %>% rename(eid = feid, BMI = f21001_0_0))

### write the final-cleaned Plasma viscosity into file.
data.table::fwrite(PV_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/PV.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### combine the gp_clinical and covariates data.
INR_joined = INR_info %>%
  inner_join(cova_info, by = c('eid' = 'Wenjian')) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18) %>%
  left_join(bmi_info %>% select(feid, f21001_0_0) %>% rename(eid = feid, BMI = f21001_0_0))

### write the final-cleaned International normalised ratio into file.
data.table::fwrite(INR_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/INR.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### combine the gp_clinical and covariates data.
APTT_joined = APTT_info %>%
  inner_join(cova_info, by = c('eid' = 'Wenjian')) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18) %>%
  left_join(bmi_info %>% select(feid, f21001_0_0) %>% rename(eid = feid, BMI = f21001_0_0))

### write the final-cleaned International normalised ratio into file.
data.table::fwrite(APTT_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/APTT.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

