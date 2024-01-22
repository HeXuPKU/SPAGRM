
###### pre-clean the white blood cell classification phenotype gp_clinical table.
library(dplyr)
library(tidyr)

### read in the neutrophil, eosinophil, basophil, lymphocyte, monocyte phecode.
source("/gdata01/user/xuhe/family_relatedness/real_data-2022-12-29/1.extract_pheno/phecode-2023-05-10XH.R")

### load the cleaned gp_clinical table.
gp_clinical = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/cleaned_gp_clinical.txt")

### extract NEUT, EOS and BASO from gp_clinical data.
GRAN = gp_clinical %>%
  filter(grepl(paste(NEUT_codes, EOS_codes, BASO_codes, sep = "|"), code)) %>% 
  mutate(value1_numeric = as.numeric(value1), value2_numeric = as.numeric(value2), value3_numeric = as.numeric(value3)) %>%
  mutate(value = coalesce(value1_numeric, value2_numeric, value3_numeric)) %>%
  mutate(isNEUT = grepl("NEUT", term_description, ignore.case = T)) %>%
  mutate(isEOS = grepl("EOS", term_description, ignore.case = T)) %>%
  mutate(isBASO = grepl("BASO", term_description, ignore.case = T)) %>%
  group_by(eid, event_dt) %>%
  mutate(hasNEUT = any(isNEUT), hasEOS = any(isEOS), hasBASO = any(isBASO)) 

NEUT = GRAN %>%
  filter(isNEUT == TRUE) %>%
  filter(value > 0.5 & value < 15) %>%
  group_by(eid) %>%
  mutate(mi = n()) %>%
  mutate(mean = mean(value)) %>%
  mutate(sd = ifelse(mi == 1, 0, 
                     ifelse(sd(value) > 2, 2, sd(value)))) %>%
  mutate(lower_tail = ifelse(mi == 1, mean, mean - 3 * sd),
         higher_tail = ifelse(mi == 1, mean, mean + 4 * sd)) %>%
  filter(value >= lower_tail & value <= higher_tail) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  mutate(value_new = ifelse(rep == 1, value, 
                            ifelse(max(value) - min(value) <= 1, 
                                   mean(value), value[which.min(abs(value - mean))]))) %>%
  mutate(value_new = round(value_new, 2)) %>%
  select(eid, data_provider, event_dt, value_new) %>% 
  dplyr::rename(neut = value_new) %>%
  distinct(eid, event_dt, .keep_all = TRUE)

### write the pre-cleaned NEUT into file.
data.table::fwrite(NEUT %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/NEUT.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

EOS = GRAN %>%
  filter(hasNEUT == T & hasEOS == T & hasBASO == T) %>%
  filter(!all(value %in% c(NA, 0))) %>%
  filter(isEOS == T) %>%
  filter(value <= 0.6) %>%
  group_by(eid) %>%
  mutate(mi = n()) %>%
  mutate(mean = mean(value)) %>%
  mutate(sd = ifelse(mi == 1, 0, 
                     ifelse(sd(value) > 0.125, 0.125, sd(value)))) %>%
  mutate(lower_tail = ifelse(mi == 1, mean, mean - 4 * sd),
         higher_tail = ifelse(mi == 1, mean, mean + 4 * sd)) %>%
  filter(value >= lower_tail & value <= higher_tail) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  mutate(value_new = ifelse(rep == 1, value, 
                            ifelse(max(value) - min(value) <= 0.1, 
                                   mean(value), value[which.min(abs(value - mean))]))) %>%
  mutate(value_new = round(value_new, 2)) %>%
  select(eid, data_provider, event_dt, value_new) %>% 
  dplyr::rename(eos = value_new) %>%
  distinct(eid, event_dt, .keep_all = TRUE)

### write the pre-cleaned EOS into file.
data.table::fwrite(EOS %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/EOS.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

BASO = GRAN %>%
  filter(hasNEUT == T & hasEOS == T & hasBASO == T) %>%
  filter(!all(value %in% c(NA, 0))) %>%
  filter(isBASO == T) %>%
  filter(value <= 0.2) %>%
  group_by(eid) %>%
  mutate(mi = n()) %>%
  mutate(mean = mean(value)) %>%
  mutate(sd = ifelse(mi == 1, 0, 
                     ifelse(sd(value) > 0.024, 0.024, sd(value)))) %>%
  mutate(lower_tail = ifelse(mi == 1, mean, mean - 4 * sd),
         higher_tail = ifelse(mi == 1, mean, mean + 4 * sd)) %>%
  filter(value >= lower_tail & value <= higher_tail) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  mutate(value_new = ifelse(rep == 1, value, 
                            ifelse(max(value) - min(value) <= 0.024, 
                                   mean(value), value[which.min(abs(value - mean))]))) %>%
  mutate(value_new = round(value_new, 2)) %>%
  select(eid, data_provider, event_dt, value_new) %>% 
  dplyr::rename(baso = value_new) %>%
  distinct(eid, event_dt, .keep_all = TRUE)

### write the pre-cleaned BASO into file.
data.table::fwrite(BASO %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/BASO.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### extract LYMPH from gp_clinical data.
LYMPH = gp_clinical %>%
  filter(grepl(LYMPH_codes, code)) %>% 
  mutate(value1_numeric = as.numeric(value1), value2_numeric = as.numeric(value2), value3_numeric = as.numeric(value3)) %>%
  mutate(lymph = coalesce(value1_numeric, value2_numeric, value3_numeric)) %>%
  filter(lymph > 0.01 & lymph < 15) %>%
  group_by(eid) %>%
  mutate(mi = n()) %>%
  mutate(mean = mean(lymph)) %>%
  mutate(sd = ifelse(mi == 1, 0, 
                     ifelse(sd(lymph) > 1, 1, sd(lymph)))) %>%
  mutate(lower_tail = ifelse(mi == 1, mean, mean - 4 * sd),
         higher_tail = ifelse(mi == 1, mean, mean + 4 * sd)) %>%
  filter(lymph >= lower_tail & lymph <= higher_tail) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  mutate(lymph_new = ifelse(rep == 1, lymph, 
                            ifelse(max(lymph) - min(lymph) <= 0.5, 
                                   mean(lymph), lymph[which.min(abs(lymph - mean))]))) %>%
  mutate(lymph_new = round(lymph_new, 2)) %>%
  select(eid, data_provider, event_dt, lymph_new) %>% 
  dplyr::rename(lymph = lymph_new) %>%
  distinct(eid, event_dt, .keep_all = TRUE)

### write the pre-cleaned LYMPH into file.
data.table::fwrite(LYMPH %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/LYMPH.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### extract MONO from gp_clinical data.
MONO = gp_clinical %>%
  filter(grepl(MONO_codes, code)) %>% 
  mutate(value1_numeric = as.numeric(value1), value2_numeric = as.numeric(value2), value3_numeric = as.numeric(value3)) %>%
  mutate(mono = coalesce(value1_numeric, value2_numeric, value3_numeric)) %>%
  filter(mono > 0.01 & mono < 1.5) %>%
  group_by(eid) %>%
  mutate(mi = n()) %>%
  mutate(mean = mean(mono)) %>%
  mutate(sd = ifelse(mi == 1, 0, 
                     ifelse(sd(mono) > 0.15, 0.15, sd(mono)))) %>%
  mutate(lower_tail = ifelse(mi == 1, mean, mean - 4 * sd),
         higher_tail = ifelse(mi == 1, mean, mean + 4 * sd)) %>%
  filter(mono >= lower_tail & mono <= higher_tail) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  mutate(mono_new = ifelse(rep == 1, mono, 
                            ifelse(max(mono) - min(mono) <= 0.15, 
                                   mean(mono), mono[which.min(abs(mono - mean))]))) %>%
  mutate(mono_new = round(mono_new, 2)) %>%
  select(eid, data_provider, event_dt, mono_new) %>% 
  dplyr::rename(mono = mono_new) %>%
  distinct(eid, event_dt, .keep_all = TRUE)

### write the pre-cleaned MONO into file.
data.table::fwrite(MONO %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/MONO.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)


### extract NEUT_perc, EOS_perc and BASO_perc from gp_clinical data.
GRAN_perc = gp_clinical %>%
  filter(grepl(paste(perc_NEUT_codes, perc_EOS_codes, perc_BASO_codes, sep = "|"), code)) %>% 
  mutate(value1_numeric = as.numeric(value1), value2_numeric = as.numeric(value2), value3_numeric = as.numeric(value3)) %>%
  mutate(value = coalesce(value1_numeric, value2_numeric, value3_numeric)) %>%
  mutate(isNEUT = grepl("NEUT", term_description, ignore.case = T)) %>%
  mutate(isEOS = grepl("EOS", term_description, ignore.case = T)) %>%
  mutate(isBASO = grepl("BASO", term_description, ignore.case = T)) %>%
  group_by(eid, event_dt) %>%
  mutate(hasNEUT = any(isNEUT), hasEOS = any(isEOS), hasBASO = any(isBASO)) 

NEUT_perc = GRAN_perc %>%
  filter(isNEUT == TRUE) %>%
  filter(value > 20) %>%
  group_by(eid) %>%
  mutate(mi = n()) %>%
  mutate(mean = mean(value)) %>%
  mutate(sd = ifelse(mi == 1, 0, 
                     ifelse(sd(value) > 8, 8, sd(value)))) %>%
  mutate(lower_tail = ifelse(mi == 1, mean, mean - 4 * sd),
         higher_tail = ifelse(mi == 1, mean, mean + 4 * sd)) %>%
  filter(value >= lower_tail & value <= higher_tail) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  mutate(value_new = ifelse(rep == 1, value, 
                            ifelse(max(value) - min(value) <= 8, 
                                   mean(value), value[which.min(abs(value - mean))]))) %>%
  mutate(value_new = round(value_new, 2)) %>%
  select(eid, data_provider, event_dt, value_new) %>% 
  dplyr::rename(neut_perc = value_new) %>%
  distinct(eid, event_dt, .keep_all = TRUE)

### write the pre-cleaned NEUT into file.
data.table::fwrite(NEUT_perc %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/NEUT_perc.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

EOS_perc = GRAN_perc %>%
  filter(hasNEUT == T & hasEOS == T & hasBASO == T) %>%
  filter(!all(value %in% c(NA, 0))) %>%
  filter(isEOS == T) %>%
  filter(value >= 0.3 & value < 15) %>%
  group_by(eid) %>%
  mutate(mi = n()) %>%
  mutate(mean = mean(value)) %>%
  mutate(sd = ifelse(mi == 1, 0, 
                     ifelse(sd(value) > 1.5, 1.5, sd(value)))) %>%
  mutate(lower_tail = ifelse(mi == 1, mean, mean - 4 * sd),
         higher_tail = ifelse(mi == 1, mean, mean + 4 * sd)) %>%
  filter(value >= lower_tail & value <= higher_tail) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  mutate(value_new = ifelse(rep == 1, value, 
                            ifelse(max(value) - min(value) <= 1.5, 
                                   mean(value), value[which.min(abs(value - mean))]))) %>%
  mutate(value_new = round(value_new, 2)) %>%
  select(eid, data_provider, event_dt, value_new) %>% 
  dplyr::rename(eos_perc = value_new) %>%
  distinct(eid, event_dt, .keep_all = TRUE)

### write the pre-cleaned NEUT into file.
data.table::fwrite(EOS_perc %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/EOS_perc.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

BASO_perc = GRAN_perc %>%
  filter(hasNEUT == T & hasEOS == T & hasBASO == T) %>%
  filter(!all(value %in% c(NA, 0))) %>%
  filter(isBASO == T) %>%
  filter(value >= 0.1 & value < 3) %>%
  group_by(eid) %>%
  mutate(mi = n()) %>%
  mutate(mean = mean(value)) %>%
  mutate(sd = ifelse(mi == 1, 0, 
                     ifelse(sd(value) > 0.4, 0.4, sd(value)))) %>%
  mutate(lower_tail = ifelse(mi == 1, mean, mean - 4 * sd),
         higher_tail = ifelse(mi == 1, mean, mean + 4 * sd)) %>%
  filter(value >= lower_tail & value <= higher_tail) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  mutate(value_new = ifelse(rep == 1, value, 
                            ifelse(max(value) - min(value) <= 0.2 & min(value) > 0, 
                                   mean(value), value[which.min(abs(value - mean))]))) %>%
  mutate(value_new = round(value_new, 2)) %>%
  select(eid, data_provider, event_dt, value_new) %>% 
  dplyr::rename(baso_perc = value_new) %>%
  distinct(eid, event_dt, .keep_all = TRUE)

### write the pre-cleaned NEUT into file.
data.table::fwrite(BASO_perc %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/BASO_perc.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### extract LYMPH_perc from gp_clinical data.
LYMPH_perc = gp_clinical %>%
  filter(grepl(perc_LYMPH_codes, code)) %>% 
  mutate(value1_numeric = as.numeric(value1), value2_numeric = as.numeric(value2), value3_numeric = as.numeric(value3)) %>%
  mutate(lymph_perc = coalesce(value1_numeric, value2_numeric, value3_numeric)) %>%
  filter(lymph_perc > 5 & lymph_perc < 100) %>%
  group_by(eid) %>%
  mutate(mi = n()) %>%
  mutate(mean = mean(lymph_perc)) %>%
  mutate(sd = ifelse(mi == 1, 0, 
                     ifelse(sd(lymph_perc) > 8, 8, sd(lymph_perc)))) %>%
  mutate(lower_tail = ifelse(mi == 1, mean, mean - 4 * sd),
         higher_tail = ifelse(mi == 1, mean, mean + 4 * sd)) %>%
  filter(lymph_perc >= lower_tail & lymph_perc <= higher_tail) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  mutate(lymph_perc_new = ifelse(rep == 1, lymph_perc, 
                                 ifelse(max(lymph_perc) - min(lymph_perc) <= 8, 
                                        mean(lymph_perc), lymph_perc[which.min(abs(lymph_perc - mean))]))) %>%
  mutate(lymph_perc_new = round(lymph_perc_new, 1)) %>%
  select(eid, data_provider, event_dt, lymph_perc_new) %>% 
  dplyr::rename(lymph_perc = lymph_perc_new) %>%
  distinct(eid, event_dt, .keep_all = TRUE)

### write the pre-cleaned LYMPH_perc into file.
data.table::fwrite(LYMPH_perc %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/LYMPH_perc.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### extract MONO_perc from gp_clinical data.
MONO_perc = gp_clinical %>%
  filter(grepl(perc_MONO_codes, code)) %>% 
  mutate(value1_numeric = as.numeric(value1), value2_numeric = as.numeric(value2), value3_numeric = as.numeric(value3)) %>%
  mutate(mono_perc = coalesce(value1_numeric, value2_numeric, value3_numeric)) %>%
  filter(mono_perc > 1.5 & mono_perc < 25) %>%
  group_by(eid) %>%
  mutate(mi = n()) %>%
  mutate(mean = mean(mono_perc)) %>%
  mutate(sd = ifelse(mi == 1, 0, 
                     ifelse(sd(mono_perc) > 2, 2, sd(mono_perc)))) %>%
  mutate(lower_tail = ifelse(mi == 1, mean, mean - 4 * sd),
         higher_tail = ifelse(mi == 1, mean, mean + 4 * sd)) %>%
  filter(mono_perc >= lower_tail & mono_perc <= higher_tail) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  mutate(mono_perc_new = ifelse(rep == 1, mono_perc, 
                                 ifelse(max(mono_perc) - min(mono_perc) <= 2, 
                                        mean(mono_perc), mono_perc[which.min(abs(mono_perc - mean))]))) %>%
  mutate(mono_perc_new = round(mono_perc_new, 1)) %>%
  select(eid, data_provider, event_dt, mono_perc_new) %>% 
  dplyr::rename(mono_perc = mono_perc_new) %>%
  distinct(eid, event_dt, .keep_all = TRUE)

### write the pre-cleaned MONO_perc into file.
data.table::fwrite(MONO_perc %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/MONO_perc.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)


###### merge the gp_clinical table and covariates.
library(dplyr)
library(tidyr)
library(lubridate)

### read in the gp_clinical table and covariates.
BASO_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/BASO.txt")

EOS_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/EOS.txt")

LYMPH_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/LYMPH.txt")

MONO_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/MONO.txt")

NEUT_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/NEUT.txt")

cova_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/assessment_center/covariate.txt")

bmi_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/assessment_center/bmi_from_UKB.txt")

### combine the gp_clinical and covariates data.
BASO_joined = BASO_info %>%
  inner_join(cova_info, by = c('eid' = 'Wenjian')) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18) %>%
  left_join(bmi_info %>% select(feid, f21001_0_0) %>% rename(eid = feid, BMI = f21001_0_0))

### write the final-cleaned BASO into file.
data.table::fwrite(BASO_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/BASO.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### combine the gp_clinical and covariates data.
EOS_joined = EOS_info %>%
  inner_join(cova_info, by = c('eid' = 'Wenjian')) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18) %>%
  left_join(bmi_info %>% select(feid, f21001_0_0) %>% rename(eid = feid, BMI = f21001_0_0))

### write the final-cleaned EOS into file.
data.table::fwrite(EOS_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/EOS.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### combine the gp_clinical and covariates data.
LYMPH_joined = LYMPH_info %>%
  inner_join(cova_info, by = c('eid' = 'Wenjian')) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18) %>%
  left_join(bmi_info %>% select(feid, f21001_0_0) %>% rename(eid = feid, BMI = f21001_0_0))

### write the final-cleaned LYMPH into file.
data.table::fwrite(LYMPH_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/LYMPH.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### combine the gp_clinical and covariates data.
MONO_joined = MONO_info %>%
  inner_join(cova_info, by = c('eid' = 'Wenjian')) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18) %>%
  left_join(bmi_info %>% select(feid, f21001_0_0) %>% rename(eid = feid, BMI = f21001_0_0))

### write the final-cleaned MONO into file.
data.table::fwrite(MONO_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/MONO.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### combine the gp_clinical and covariates data.
NEUT_joined = NEUT_info %>%
  inner_join(cova_info, by = c('eid' = 'Wenjian')) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18) %>%
  left_join(bmi_info %>% select(feid, f21001_0_0) %>% rename(eid = feid, BMI = f21001_0_0))

### write the final-cleaned NEUT into file.
data.table::fwrite(NEUT_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/NEUT.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)


BASO_perc_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/BASO_perc.txt")

EOS_perc_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/EOS_perc.txt")

LYMPH_perc_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/LYMPH_perc.txt")

MONO_perc_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/MONO_perc.txt")

NEUT_perc_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/NEUT_perc.txt")

### combine the gp_clinical and covariates data.
BASO_perc_joined = BASO_perc_info %>%
  inner_join(cova_info, by = c('eid' = 'Wenjian')) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18) %>%
  left_join(bmi_info %>% select(feid, f21001_0_0) %>% rename(eid = feid, BMI = f21001_0_0))

### write the final-cleaned BASO_perc into file.
data.table::fwrite(BASO_perc_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/BASO_perc.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### combine the gp_clinical and covariates data.
EOS_perc_joined = EOS_perc_info %>%
  inner_join(cova_info, by = c('eid' = 'Wenjian')) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18) %>%
  left_join(bmi_info %>% select(feid, f21001_0_0) %>% rename(eid = feid, BMI = f21001_0_0))

### write the final-cleaned EOS_perc into file.
data.table::fwrite(EOS_perc_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/EOS_perc.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### combine the gp_clinical and covariates data.
LYMPH_perc_joined = LYMPH_perc_info %>%
  inner_join(cova_info, by = c('eid' = 'Wenjian')) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18) %>%
  left_join(bmi_info %>% select(feid, f21001_0_0) %>% rename(eid = feid, BMI = f21001_0_0))

### write the final-cleaned LYMPH_perc into file.
data.table::fwrite(LYMPH_perc_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/LYMPH_perc.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### combine the gp_clinical and covariates data.
MONO_perc_joined = MONO_perc_info %>%
  inner_join(cova_info, by = c('eid' = 'Wenjian')) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18) %>%
  left_join(bmi_info %>% select(feid, f21001_0_0) %>% rename(eid = feid, BMI = f21001_0_0))

### write the final-cleaned MONO_perc into file.
data.table::fwrite(MONO_perc_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/MONO_perc.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### combine the gp_clinical and covariates data.
NEUT_perc_joined = NEUT_perc_info %>%
  inner_join(cova_info, by = c('eid' = 'Wenjian')) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18) %>%
  left_join(bmi_info %>% select(feid, f21001_0_0) %>% rename(eid = feid, BMI = f21001_0_0))

### write the final-cleaned NEUT_perc into file.
data.table::fwrite(NEUT_perc_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/NEUT_perc.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)
