
library(dplyr)
library(tidyr)

### read in the bp phecode.
source("/gdata01/user/xuhe/family_relatedness/real_data-2022-12-29/1.extract_pheno/phecode-2023-05-10XH.R")

### load the cleaned gp_clinical table
gp_clinical = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/cleaned_gp_clinical.txt")

### extract blood pressure.
bp = gp_clinical %>%
  filter(grepl(BP_codes, code)) %>%
  filter(value1 != "" | value2 != "" | value3 != "")  %>%
  mutate(value1 = as.numeric(value1), value2 = as.numeric(value2), value3 = as.numeric(value3)) # this warning doesn't make sense.

### display the data.
bp %>% group_by(code, term_description) %>% summarize(n = n()) %>% arrange(desc(n))

### 172 rows of value3 contains meaningless values.
bp$value3[bp$value3 < 10] = NA

### multiple values per record. Take the larger value to be systolic and the smaller value to be diastolic. Filter out any records where either of these values are 0.
bp_mult = bp %>% 
  rowwise() %>% 
  filter(sum(!is.na(value1), !is.na(value2), !is.na(value3)) == 2) %>% 
  ungroup() %>% 
  mutate(Systolic_bp_pc = pmax(value1, value2, value3, na.rm=T)) %>% 
  mutate(Diastolic_bp_pc = pmin(value1, value2, value3, na.rm=T)) %>% 
  filter(Systolic_bp_pc != 0 & Diastolic_bp_pc != 0)

### pre-QC according to blood pressure distribution in UKB.
bp_mult_less = bp_mult %>% 
  filter(Systolic_bp_pc - Diastolic_bp_pc >= 10) %>%
  filter(Systolic_bp_pc > 50 & Diastolic_bp_pc > 30 & Systolic_bp_pc < 270) %>% # pre-filter unreasonable values.
  select(eid, event_dt, data_provider, terminology, Systolic_bp_pc, Diastolic_bp_pc, code, term_description) %>%
  dplyr::rename(term_description_both = term_description, code_both = code)

### one value per record. create one column with the value, where previously there were three columns that may contain the value. 
bp_single = bp %>% 
  rowwise() %>%
  filter(sum(!is.na(value1), !is.na(value2), !is.na(value3)) == 1) %>%
  ungroup() %>%
  mutate(value = coalesce(value1, value2, value3)) %>%
  filter(value >= 30 & value < 270) %>%
  mutate(isSBP = grepl("systolic", term_description, ignore.case=T)) %>%
  mutate(isDBP = grepl("diastolic", term_description, ignore.case=T)) %>%
  group_by(eid, event_dt) %>%
  mutate(bothBP = any(isSBP) & any(isDBP))

### pre-QC according to blood pressure distribution in UKB.
bp_single_less = bp_single %>%
  filter(bothBP == TRUE) %>%
  filter(isSBP == TRUE | isDBP == TRUE) %>%
  mutate(numSBP = sum(isSBP), numDBP = sum(isDBP)) %>%
  filter(abs(numSBP - numDBP) <= 1) %>%
  group_by(eid, event_dt, isSBP) %>%
  mutate(meanSBP = ifelse(isSBP == TRUE, mean(value), NA),
         meanDBP = ifelse(isSBP == FALSE, mean(value), NA)) %>%
  group_by(eid, event_dt) %>%
  fill(meanSBP, meanDBP, .direction = "downup") %>%
  filter(meanSBP - meanDBP >= 10) %>%
  filter((value - meanDBP >= 10 & isSBP ==TRUE) | (meanSBP - value >= 10 & isDBP ==TRUE)) %>%
  distinct(eid, data_provider, event_dt, value, isSBP, isDBP)
  
### look at the 'Unknown' codes to see if we can figure out whether they are systolic or diastolic.  
# For many of these, the same code is given twice, each with a different value. 
# Sometimes an 'unknown' code is a duplicate of a systolic or diastolic measurement. 
# If there are two *unique* values given per ID/date, then we can assume they are systolic (higher) and diastolic (lower). 
# Otherwise, we discard that set of values. 
bp_unknown = bp_single %>%
  filter(bothBP == FALSE) %>%
  rowwise() %>%
  mutate(isunknown = !any(isSBP, isDBP)) %>%
  group_by(eid, event_dt) %>%
  filter(any(isunknown)) %>%
  mutate(n = length(unique(value))) %>%
  filter(n == 2) %>% 
  mutate(SBP = max(value), DBP = min(value)) %>%
  filter(SBP - DBP > 15) %>%
  filter(SBP > 50 & DBP > 30 & SBP < 270) %>%
  distinct(eid, data_provider, event_dt, SBP, DBP) %>%
  pivot_longer(cols = c("SBP", "DBP"), names_to = "bp_type", values_to = "value") %>%
  mutate(isSBP = grepl("SBP", bp_type), isDBP = grepl("DBP", bp_type))

### combine the bp_single_less and bp_unknown.
bp_single_more = full_join(bp_single_less, bp_unknown) %>%
  arrange(eid,event_dt)

### noting that ~37w rows of bp_single data has appeared in bp_mult, most are repeat records.
repeat_eid_event = inner_join(bp_mult_less %>% group_by(eid, event_dt) %>% summarize(n_mult=n()),
                              bp_single_more %>% group_by(eid, event_dt) %>% summarize(n_single=n()))

### thus, we only extract bp_single data with unique eid and unique event_dt.
bp_single_unique = bp_single_more %>%
  left_join(repeat_eid_event) %>%
  filter(is.na(n_single)) %>%
  select(- c(n_mult, n_single))

### extract SBP from bp_mult and bp_single for final QC.
systolic_single = full_join(bp_mult_less %>% 
                              select(eid, event_dt, data_provider, Systolic_bp_pc) %>%
                              dplyr::rename(SBP = Systolic_bp_pc),
                            bp_single_unique %>%
                              filter(isSBP == TRUE) %>%
                              select(eid, event_dt, data_provider, value) %>%
                              dplyr::rename(SBP = value)) %>%
  arrange(eid, event_dt) %>%
  group_by(eid) %>%
  mutate(mi = n()) %>%
  mutate(mean = mean(SBP)) %>%
  mutate(sd = ifelse(mi == 1, 0, 
                     ifelse(sd(SBP) > 17.5, 17.5, sd(SBP)))) %>%
  mutate(lower_tail = ifelse(mi == 1, mean, mean - 4 * sd),
         higher_tail = ifelse(mi == 1, mean, mean + 4 * sd)) %>%
  filter(SBP >= lower_tail & SBP <= higher_tail) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  mutate(SBP_new = ifelse(rep == 1, SBP, mean(SBP))) %>%
  mutate(SBP_new = round(SBP_new, 0)) %>%
  select(eid, event_dt, data_provider, SBP_new) %>%
  dplyr::rename(SBP = SBP_new) %>%
  distinct(eid, event_dt, SBP, .keep_all = TRUE)

### extract DBP from bp_mult and bp_single for final QC.
diastolic_single = full_join(bp_mult_less %>% 
                              select(eid, event_dt, data_provider, Diastolic_bp_pc) %>%
                              dplyr::rename(DBP = Diastolic_bp_pc),
                            bp_single_unique %>%
                              filter(isDBP == TRUE) %>%
                              select(eid, event_dt, data_provider, value) %>%
                              dplyr::rename(DBP = value)) %>%
  arrange(eid, event_dt) %>%
  group_by(eid) %>%
  mutate(mi = n()) %>%
  mutate(mean = mean(DBP)) %>%
  mutate(sd = ifelse(mi == 1, 0, 
                     ifelse(sd(DBP) > 12.5, 12.5, sd(DBP)))) %>%
  mutate(lower_tail = ifelse(mi == 1, mean, mean - 4 * sd),
         higher_tail = ifelse(mi == 1, mean, mean + 4 * sd)) %>%
  filter(DBP >= lower_tail & DBP <= higher_tail) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  mutate(DBP_new = ifelse(rep == 1, DBP, mean(DBP))) %>%
  mutate(DBP_new = round(DBP_new, 0)) %>%
  select(eid, event_dt, data_provider, DBP_new) %>%
  dplyr::rename(DBP = DBP_new) %>%
  distinct(eid, event_dt, DBP, .keep_all = TRUE)

### combine cleaned SBP and DBP.
bp_clean = full_join(systolic_single, diastolic_single) %>%
  drop_na(SBP, DBP) %>%
  arrange(eid, event_dt)

### write the pre-cleaned bp into file.
data.table::fwrite(bp_clean, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/bp.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)


###### merge the gp_clinical table and covariates.
library(dplyr)
library(tidyr)
library(lubridate)

### read in the gp_clinical table and covariates.
bp_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/bp.txt")

cova_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/assessment_center/covariate.txt")

drug_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/assessment_center/covariate_drugs.txt")

bmi_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/assessment_center/bmi_from_UKB.txt")

### combine the gp_clinical and covariates data.
bp_joined = bp_info %>%
  mutate(PP = SBP - DBP) %>%
  inner_join(cova_info, by = c('eid' = 'Wenjian')) %>%
  arrange(eid, event_dt) %>%
  mutate(age = time_length(interval(birth_dt, event_dt), 'year')) %>%
  filter(event_dt >= "1980-01-01") %>%
  filter(age >= 18)

### add the bmi and medication data.
bp_joined = bp_joined %>%
  left_join(drug_info %>% select(feid, self_bpdrugs), by = c('eid' = 'feid')) %>%
  left_join(bmi_info %>% select(feid, f21001_0_0) %>% rename(eid = feid, BMI = f21001_0_0)) %>%
  mutate(SBP_shifted = ifelse(self_bpdrugs == 1, SBP + 15, SBP)) %>%
  mutate(DBP_shifted = ifelse(self_bpdrugs == 1, DBP + 10, DBP)) %>%
  mutate(PP_shifted = ifelse(self_bpdrugs == 1, PP + 5, PP))

# samples with missing BMI: 585.
length(unique(bp_joined$eid))
length(unique(bp_joined %>% filter(!is.na(BMI)) %>% select(eid) %>% unlist))

# samples with missing drug information: 1008.
length(unique(bp_joined$eid))
length(unique(bp_joined %>% filter(self_bpdrugs != -1) %>% select(eid) %>% unlist))

# remaining samples: 182833 -> 181245.
bp_left = bp_joined %>%
  filter(!is.na(BMI)) %>%
  filter(self_bpdrugs != -1)

length(unique(bp_joined$eid))
length(unique(bp_left$eid))

### write the final-cleaned bp into file.
data.table::fwrite(bp_joined, 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/bp.txt", 
                   row.names = F, col.names = T, sep = "\t", quote = T)
