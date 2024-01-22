
###### pre-clean the lung function phenotype gp_clinical table.
library(dplyr)
library(tidyr)

### read in the RBC, WBC and Platelet phecode.
source("/gdata01/user/xuhe/family_relatedness/real_data-2022-12-29/1.extract_pheno/phecode-2023-05-10XH.R")

### load the cleaned gp_clinical table.
gp_clinical = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/cleaned_gp_clinical.txt")

### extract red blood cell (RBC) from gp_clinical data.
RBC = gp_clinical %>%
  filter(grepl(RBC_codes, code)) %>% 
  mutate(value1_numeric = as.numeric(value1), value2_numeric = as.numeric(value2), value3_numeric = as.numeric(value3)) %>%
  mutate(rbc = coalesce(value1_numeric, value2_numeric, value3_numeric)) %>%
  filter(rbc > 1 & rbc < 8) %>%
  mutate(value3 = toupper(value3)) %>%
  filter(!(value3 %in% c("10*-2", "MEA035", "MEA154", "MEA184"))) %>%
  group_by(eid) %>%
  mutate(mi = n()) %>%
  mutate(mean = mean(rbc)) %>%
  mutate(sd = ifelse(mi == 1, 0, 
                     ifelse(sd(rbc) > 0.4, 0.4, sd(rbc)))) %>%
  mutate(lower_tail = ifelse(mi == 1, mean, mean - 4 * sd),
         higher_tail = ifelse(mi == 1, mean, mean + 4 * sd)) %>%
  filter(rbc >= lower_tail & rbc <= higher_tail) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  mutate(rbc_new = ifelse(rep == 1, rbc, 
                           ifelse(max(rbc) - min(rbc) <= 0.4, 
                                  mean(rbc), rbc[which.min(abs(rbc - mean))]))) %>%
  mutate(rbc_new = round(rbc_new, 2)) %>%
  select(eid, data_provider, event_dt, rbc_new) %>% 
  dplyr::rename(rbc = rbc_new) %>%
  distinct(eid, event_dt, .keep_all = TRUE)

### write the pre-cleaned red blood cell (RBC) into file.
data.table::fwrite(RBC %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/RBC.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### extract white blood cell (WBC) from gp_clinical data.
WBC = gp_clinical %>%
  filter(grepl(WBC_codes, code)) %>% 
  mutate(value1_numeric = as.numeric(value1), value2_numeric = as.numeric(value2), value3_numeric = as.numeric(value3)) %>%
  mutate(wbc = coalesce(value1_numeric, value2_numeric, value3_numeric)) %>%
  filter(wbc > 1 & wbc < 16) %>%
  mutate(value3 = toupper(value3)) %>%
  filter(value3 != "/UL") %>%
  group_by(eid) %>%
  mutate(mi = n()) %>%
  mutate(mean = mean(wbc)) %>%
  mutate(sd = ifelse(mi == 1, 0, 
                     ifelse(sd(wbc) > 2, 2, sd(wbc)))) %>%
  mutate(lower_tail = ifelse(mi == 1, mean, mean - 4 * sd),
         higher_tail = ifelse(mi == 1, mean, mean + 4 * sd)) %>%
  filter(wbc >= lower_tail & wbc <= higher_tail) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  mutate(wbc_new = ifelse(rep == 1, wbc, 
                          ifelse(max(wbc) - min(wbc) <= 1, 
                                 mean(wbc), wbc[which.min(abs(wbc - mean))]))) %>%
  mutate(wbc_new = round(wbc_new, 2)) %>%
  select(eid, data_provider, event_dt, wbc_new) %>% 
  dplyr::rename(wbc = wbc_new) %>%
  distinct(eid, event_dt, .keep_all = TRUE)

### write the pre-cleaned white blood cell (WBC) into file.
data.table::fwrite(WBC %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/WBC.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### extract platelet from gp_clinical data.
Platelet = gp_clinical %>%
  filter(grepl(Platelet_codes, code)) %>% 
  mutate(value1_numeric = as.numeric(value1), value2_numeric = as.numeric(value2), value3_numeric = as.numeric(value3)) %>%
  mutate(pla = coalesce(value1_numeric, value2_numeric, value3_numeric)) %>%
  filter(pla > 10 & pla < 800) %>%
  mutate(value3 = toupper(value3)) %>%
  filter(value3 != "MEA038") %>%
  group_by(eid) %>%
  mutate(mi = n()) %>%
  mutate(mean = mean(pla)) %>%
  mutate(sd = ifelse(mi == 1, 0, 
                     ifelse(sd(pla) > 75, 75, sd(pla)))) %>%
  mutate(lower_tail = ifelse(mi == 1, mean, mean - 4 * sd),
         higher_tail = ifelse(mi == 1, mean, mean + 4 * sd)) %>%
  filter(pla >= lower_tail & pla <= higher_tail) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  mutate(pla_new = ifelse(rep == 1, pla, 
                          ifelse(max(pla) - min(pla) <= 10, 
                                 mean(pla), pla[which.min(abs(pla - mean))]))) %>%
  mutate(pla_new = round(pla_new, 0)) %>%
  select(eid, data_provider, event_dt, pla_new) %>% 
  dplyr::rename(pla = pla_new) %>%
  distinct(eid, event_dt, .keep_all = TRUE)

### write the pre-cleaned platelet into file.
data.table::fwrite(Platelet %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/Platelet.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)


###### merge the gp_clinical table and covariates.
library(dplyr)
library(tidyr)
library(lubridate)

### read in the gp_clinical table and covariates.
RBC_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/RBC.txt")

WBC_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/WBC.txt")

Platelet_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/Platelet.txt")

cova_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/assessment_center/covariate.txt")

bmi_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/assessment_center/bmi_from_UKB.txt")

### combine the gp_clinical and covariates data.
RBC_joined = RBC_info %>%
  inner_join(cova_info, by = c('eid' = 'Wenjian')) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18) %>%
  left_join(bmi_info %>% select(feid, f21001_0_0) %>% rename(eid = feid, BMI = f21001_0_0))

### write the final-cleaned spo2 into file.
data.table::fwrite(RBC_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/RBC.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### combine the gp_clinical and covariates data.
WBC_joined = WBC_info %>%
  inner_join(cova_info, by = c('eid' = 'Wenjian')) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18) %>%
  left_join(bmi_info %>% select(feid, f21001_0_0) %>% rename(eid = feid, BMI = f21001_0_0))

### write the final-cleaned spo2 into file.
data.table::fwrite(WBC_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/WBC.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### combine the gp_clinical and covariates data.
Platelet_joined = Platelet_info %>%
  inner_join(cova_info, by = c('eid' = 'Wenjian')) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18) %>%
  left_join(bmi_info %>% select(feid, f21001_0_0) %>% rename(eid = feid, BMI = f21001_0_0))

### write the final-cleaned spo2 into file.
data.table::fwrite(Platelet_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/Platelet.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)
