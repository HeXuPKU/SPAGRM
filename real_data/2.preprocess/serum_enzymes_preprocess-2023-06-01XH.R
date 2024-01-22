
###### pre-clean the white blood cell classification phenotype gp_clinical table.
library(dplyr)
library(tidyr)

### read in the alkaline phosphatase (ALP), alanine aminotransferase (ALT), aspartate aminotransferase (AST), gamma-glutamyl transferase (GGT) phecode.
source("/gdata01/user/xuhe/family_relatedness/real_data-2022-12-29/1.extract_pheno/phecode-2023-05-10XH.R")

### load the cleaned gp_clinical table.
gp_clinical = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/cleaned_gp_clinical.txt")

### extract ALP from gp_clinical data.
ALP = gp_clinical %>%
  filter(grepl(ALP_codes, code)) %>% 
  mutate(value1_numeric = as.numeric(value1), value2_numeric = as.numeric(value2), value3_numeric = as.numeric(value3)) %>%
  mutate(alp = coalesce(value1_numeric, value2_numeric, value3_numeric)) %>%
  filter(alp > 10 & alp < 500) %>%
  mutate(value3 = toupper(value3)) %>%
  filter(value3 %in% c("", "IU/L", "MEA000", "MEA061", "U/L")) %>%
  group_by(eid) %>%
  mutate(mi = n()) %>%
  mutate(mean = mean(alp)) %>%
  mutate(sd = ifelse(mi == 1, 0, 
                     ifelse(sd(alp) > 40, 40, sd(alp)))) %>%
  mutate(lower_tail = ifelse(mi == 1, mean, mean - 4 * sd),
         higher_tail = ifelse(mi == 1, mean, mean + 4 * sd)) %>%
  filter(alp >= lower_tail & alp <= higher_tail) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  mutate(alp_new = ifelse(rep == 1, alp, 
                          ifelse(max(alp) - min(alp) <= 5, 
                                 mean(alp), alp[which.min(abs(alp - mean))]))) %>%
  mutate(alp_new = round(alp_new, 0)) %>%
  select(eid, data_provider, event_dt, alp_new) %>% 
  dplyr::rename(alp = alp_new) %>%
  distinct(eid, event_dt, .keep_all = TRUE)

### write the pre-cleaned ALP into file.
data.table::fwrite(ALP %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/ALP.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### extract ALT from gp_clinical data.
ALT = gp_clinical %>%
  filter(grepl(ALT_codes, code)) %>% 
  mutate(value1_numeric = as.numeric(value1), value2_numeric = as.numeric(value2), value3_numeric = as.numeric(value3)) %>%
  mutate(alt = coalesce(value1_numeric, value2_numeric, value3_numeric)) %>%
  filter(alt > 3 & alt < 200) %>%
  mutate(value3 = toupper(value3)) %>%
  filter(value3 %in% c("", "IU/L", "MEA000", "MEA061", "MEA127", "U/L")) %>%
  group_by(eid) %>%
  mutate(mi = n()) %>%
  mutate(mean = mean(alt)) %>%
  mutate(sd = ifelse(mi == 1, 0, 
                     ifelse(sd(alt) > 15, 15, sd(alt)))) %>%
  mutate(lower_tail = ifelse(mi == 1, mean, mean - 4 * sd),
         higher_tail = ifelse(mi == 1, mean, mean + 4 * sd)) %>%
  filter(alt >= lower_tail & alt <= higher_tail) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  mutate(alt_new = ifelse(rep == 1, alt, 
                          ifelse(max(alt) - min(alt) <= 5, 
                                 mean(alt), alt[which.min(abs(alt - mean))]))) %>%
  mutate(alt_new = round(alt_new, 0)) %>%
  select(eid, data_provider, event_dt, alt_new) %>% 
  dplyr::rename(alt = alt_new) %>%
  distinct(eid, event_dt, .keep_all = TRUE)

### write the pre-cleaned ALT into file.
data.table::fwrite(ALT %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/ALT.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### extract AST from gp_clinical data.
AST = gp_clinical %>%
  filter(grepl(AST_codes, code)) %>% 
  mutate(value1_numeric = as.numeric(value1), value2_numeric = as.numeric(value2), value3_numeric = as.numeric(value3)) %>%
  mutate(ast = coalesce(value1_numeric, value2_numeric, value3_numeric)) %>%
  filter(ast > 3 & ast < 150) %>%
  mutate(value3 = toupper(value3)) %>%
  group_by(eid) %>%
  mutate(mi = n()) %>%
  mutate(mean = mean(ast)) %>%
  mutate(sd = ifelse(mi == 1, 0, 
                     ifelse(sd(ast) > 10, 10, sd(ast)))) %>%
  mutate(lower_tail = ifelse(mi == 1, mean, mean - 4 * sd),
         higher_tail = ifelse(mi == 1, mean, mean + 4 * sd)) %>%
  filter(ast >= lower_tail & ast <= higher_tail) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  mutate(ast_new = ifelse(rep == 1, ast, 
                          ifelse(max(ast) - min(ast) <= 5, 
                                 mean(ast), ast[which.min(abs(ast - mean))]))) %>%
  mutate(ast_new = round(ast_new, 0)) %>%
  select(eid, data_provider, event_dt, ast_new) %>% 
  dplyr::rename(ast = ast_new) %>%
  distinct(eid, event_dt, .keep_all = TRUE)

### write the pre-cleaned AST into file.
data.table::fwrite(AST %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/AST.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### extract GGT from gp_clinical data.
GGT = gp_clinical %>%
  filter(grepl(GGT_codes, code)) %>% 
  mutate(value1_numeric = as.numeric(value1), value2_numeric = as.numeric(value2), value3_numeric = as.numeric(value3)) %>%
  mutate(ggt = coalesce(value1_numeric, value2_numeric, value3_numeric)) %>%
  filter(ggt > 3 & ggt < 500) %>%
  mutate(value3 = toupper(value3)) %>%
  filter(value3 %in% c("", "I.U./L", "IU/L", "MEA000", "MEA061", "MEA127", "U/L")) %>%
  group_by(eid) %>%
  mutate(mi = n()) %>%
  mutate(mean = mean(ggt)) %>%
  mutate(sd = ifelse(mi == 1, 0, 
                     ifelse(sd(ggt) > 40, 40, sd(ggt)))) %>%
  mutate(lower_tail = ifelse(mi == 1, mean, mean - 4 * sd),
         higher_tail = ifelse(mi == 1, mean, mean + 4 * sd)) %>%
  filter(ggt >= lower_tail & ggt <= higher_tail) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  mutate(ggt_new = ifelse(rep == 1, ggt, 
                          ifelse(max(ggt) - min(ggt) <= 5, 
                                 mean(ggt), ggt[which.min(abs(ggt - mean))]))) %>%
  mutate(ggt_new = round(ggt_new, 0)) %>%
  select(eid, data_provider, event_dt, ggt_new) %>% 
  dplyr::rename(ggt = ggt_new) %>%
  distinct(eid, event_dt, .keep_all = TRUE)

### write the pre-cleaned GGT into file.
data.table::fwrite(GGT %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/GGT.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)


###### merge the gp_clinical table and covariates.
library(dplyr)
library(tidyr)
library(lubridate)

### read in the gp_clinical table and covariates.
ALP_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/ALP.txt")

ALT_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/ALT.txt")

AST_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/AST.txt")

GGT_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/GGT.txt")

cova_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/assessment_center/covariate.txt")

bmi_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/assessment_center/bmi_from_UKB.txt")

### combine the gp_clinical and covariates data.
ALP_joined = ALP_info %>%
  inner_join(cova_info, by = c('eid' = 'Wenjian')) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18) %>%
  left_join(bmi_info %>% select(feid, f21001_0_0) %>% rename(eid = feid, BMI = f21001_0_0))

### write the final-cleaned ALP into file.
data.table::fwrite(ALP_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/ALP.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### combine the gp_clinical and covariates data.
ALT_joined = ALT_info %>%
  inner_join(cova_info, by = c('eid' = 'Wenjian')) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18) %>%
  left_join(bmi_info %>% select(feid, f21001_0_0) %>% rename(eid = feid, BMI = f21001_0_0))

### write the final-cleaned ALT into file.
data.table::fwrite(ALT_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/ALT.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### combine the gp_clinical and covariates data.
AST_joined = AST_info %>%
  inner_join(cova_info, by = c('eid' = 'Wenjian')) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18) %>%
  left_join(bmi_info %>% select(feid, f21001_0_0) %>% rename(eid = feid, BMI = f21001_0_0))

### write the final-cleaned AST into file.
data.table::fwrite(AST_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/AST.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)

### combine the gp_clinical and covariates data.
GGT_joined = GGT_info %>%
  inner_join(cova_info, by = c('eid' = 'Wenjian')) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18) %>%
  left_join(bmi_info %>% select(feid, f21001_0_0) %>% rename(eid = feid, BMI = f21001_0_0))

### write the final-cleaned GGT into file.
data.table::fwrite(GGT_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/GGT.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)
