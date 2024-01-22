
###### pre-clean the gp_clinical table.
library(dplyr)
library(tidyr)
library(lubridate)

# read in the dictonary of read2 and read3 codes and descriptions.
raw_readv2_dictinary = openxlsx::read.xlsx("/gdata02/master_data1/UK_Biobank/ukb669197/all_lkps_maps_v3.xlsx",
                                           sheet = "read_v2_lkp")

raw_readctv3_dictinary = openxlsx::read.xlsx("/gdata02/master_data1/UK_Biobank/ukb669197/all_lkps_maps_v3.xlsx",
                                             sheet = "read_ctv3_lkp")

# delete the last row that contains meaningless information.
raw_readv2_dictinary = raw_readv2_dictinary[-nrow(raw_readv2_dictinary),]
raw_readctv3_dictinary = raw_readctv3_dictinary[-nrow(raw_readctv3_dictinary),]

# check whether Read v2 and Read v3 codes overlap.
all(raw_readv2_dictinary$read_code %in% raw_readctv3_dictinary$read_code)

# import TPP local codes which are also used in primary care data, 
# download from https://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=8708
raw_tpp_local = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/tpp_local_code.tsv", quote="")
tpp_local = raw_tpp_local[-c(1:2),]

# TPP codes do not appear in Read v2 or Read v3 dictionaries.
tpp_local %>% filter(coding %in% raw_readv2_dictinary$read_code)
tpp_local %>% filter(coding %in% raw_readctv3_dictinary$read_code)

# merge all of the dictionaries.
full_dict =
  full_join(raw_readv2_dictinary %>% rename(code = read_code) %>% mutate(terminology="read2"),
            raw_readctv3_dictinary %>% rename(code = read_code) %>% mutate(terminology="read3")) %>%
  full_join(tpp_local %>% rename(code = coding, term_description = meaning) %>% mutate(terminology="read3")) %>%
  distinct(code, term_description, terminology) %>%
  group_by(code) %>% 
  summarize(term_description = paste(term_description, collapse = " | "))

# read in the raw gp_clinical table.
raw_gp_clinical = data.table::fread("/gdata02/master_data1/UK_Biobank/ukb669197/gp_clinical.txt")

# special_date is downloaded from https://biobank.ndph.ox.ac.uk/showcase/coding.cgi?tk=68EMKkCB8skqob929k0lJTa0eLbPoYnY157881&id=819
special_date = c("1900-01-01", "1901-01-01", "1902-02-02", "1903-03-03", "2037-07-07")

# fix the event date and create code and terminology fields. (this might take minutes)
gp_clinical = raw_gp_clinical %>% 
  mutate(event_dt = dmy(event_dt)) %>% # Fix the date format
  filter(!is.na(event_dt)) %>% # filter out missing event_dt
  filter(!(event_dt %in% ymd(special_date))) %>% # filter out special data
  mutate(code = ifelse(read_2 != "", read_2, read_3)) %>% # merge the code
  mutate(terminology = ifelse(read_2 != "", "read2", "read3")) %>% # add the terminology
  select(-read_2, -read_3) %>%
  distinct()

# add the term descriptions to the gp_clinical table.
gp_clinical = gp_clinical %>%
  left_join(full_dict)

# write the pre-cleaned gp_clinical data into files. (after extracting all phenotypes, this file should be removed, ~12GB)
data.table::fwrite(gp_clinical,
                   "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/cleaned_gp_clinical.txt",
                   row.names = F, col.names = T, sep = "\t")

# there are 230056 subjects in gp_clinical, but 230036 in UKBB website (https://biobank.ndph.ox.ac.uk/ukb/field.cgi?id=42040)
SubjID = unique(gp_clinical$eid)
length(SubjID) # 229912 subjects with proper event_dt remained.

data.table::fwrite(data.table::data.table(SubjID),
                   "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/GRM_IBD/SubjIDinGP_clinical.txt",
                   row.names = F, col.names = F, quote = F, sep = "\t")


### pre-clean the covariate data.
library(dplyr)
library(tidyr)
library(lubridate)
source("/gdata02/master_data1/UK_Biobank/Data_Tools/extractPheno.R")

# Read and filter the covariate data.
raw_covariate = data.table::fread("/gdata02/master_data1/UK_Biobank/Data_Tools/dataMergeUmich_2022-05-09.csv")

covariate = raw_covariate %>% filter(Wenjian > 0 & WhiteBritish == 1) %>%
  select(Wenjian, sex_genetic, birthYear, birthMonth, "PC1.wb-Umich":"PC10.wb-Umich", UNRELATED_WB, genoArray) %>%
  dplyr::rename(PC = "PC1.wb-Umich":"PC10.wb-Umich") %>%
  drop_na(sex_genetic, birthYear, birthMonth, PC1, genoArray) %>%
  rowwise() %>% mutate(birth_dt = paste(c(birthYear, birthMonth, 1), collapse = "-")) %>%
  ungroup() %>% mutate(birth_dt = ymd(birth_dt)) %>%
  select(-c(birthYear, birthMonth)) %>%
  select(Wenjian, sex_genetic, birth_dt, PC1:PC10, genoArray, UNRELATED_WB)

covariate$genoArray = factor(covariate$genoArray, 
                             levels = c("BiLEVEAX", "Axiom"), labels = c(1, 0))

data.table::fwrite(covariate,
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/assessment_center/covariate.txt",
                   row.names = F, col.names = T, quote = F, sep = "\t")


### extract bmi data from UKB assessment center.
library(dplyr)
library(tidyr)
source("/gdata02/master_data1/UK_Biobank/Data_Tools/extractPheno.R")

raw_data_center_height = extractPhenoFromUKBB(fieldID = 50)

raw_data_center_weight = extractPhenoFromUKBB(fieldID = 21002)

raw_data_center_BMI = extractPhenoFromUKBB(fieldID = 21001)

data.table::fwrite(raw_data_center_height$fieldData,
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/assessment_center/height_from_UKB.txt",
                   row.names = F, col.names = T, quote = F, sep = "\t")

data.table::fwrite(raw_data_center_weight$fieldData,
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/assessment_center/weight_from_UKB.txt",
                   row.names = F, col.names = T, quote = F, sep = "\t")

data.table::fwrite(raw_data_center_BMI$fieldData,
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/assessment_center/bmi_from_UKB.txt",
                   row.names = F, col.names = T, quote = F, sep = "\t")

raw_data_center_WC = extractPhenoFromUKBB(fieldID = 48)

data.table::fwrite(raw_data_center_WC$fieldData,
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/assessment_center/waistcircu_from_UKB.txt",
                   row.names = F, col.names = T, quote = F, sep = "\t")

### extract drug data from UKB assessment center.
library(dplyr)
library(tidyr)
source("/gdata02/master_data1/UK_Biobank/Data_Tools/extractPheno.R")

raw_data_center_drugs = extractPhenoFromUKBB(fieldID = c(6177, 6153)) # 6177 for male and 6153 for female.

data.table::fwrite(raw_data_center_drugs$fieldData,
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/assessment_center/drugs_from_UKB.txt",
                   row.names = F, col.names = T, quote = F, sep = "\t")

### extract five medication data from raw_data_center_drugs.
drugs_info = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/assessment_center/drugs_from_UKB.txt")

data_center_drugs = drugs_info %>%
  select(feid, f6153_0_0:f6153_0_3, f6177_0_0:f6177_0_2) %>%
  mutate(across(f6153_0_0:f6177_0_2, as.numeric)) %>%
  rowwise() %>%
  mutate(self_cholesteroldrugs = case_when(
    any(c_across(f6153_0_0:f6177_0_2) == 1) ~ 1,
    all(c_across(f6153_0_0:f6177_0_2) %in% c(NA, -1, -3)) ~ -1)) %>%
  mutate(self_bpdrugs = case_when(
    any(c_across(f6153_0_0:f6177_0_2) == 2) ~ 1,
    all(c_across(f6153_0_0:f6177_0_2) %in% c(NA, -1, -3)) ~ -1)) %>%
  mutate(self_insulin = case_when(
    any(c_across(f6153_0_0:f6177_0_2) == 3) ~ 1,
    all(c_across(f6153_0_0:f6177_0_2) %in% c(NA, -1, -3)) ~ -1)) %>%
  mutate(self_HRT = case_when(
    any(c_across(f6153_0_0:f6177_0_2) == 4) ~ 1,
    all(c_across(f6153_0_0:f6177_0_2) %in% c(NA, -1, -3)) ~ -1)) %>%
  mutate(self_contraceptivedrugs = case_when(
    any(c_across(f6153_0_0:f6177_0_2) == 5) ~ 1,
    all(c_across(f6153_0_0:f6177_0_2) %in% c(NA, -1, -3)) ~ -1)) %>%
  ungroup() %>% 
  mutate(across(self_cholesteroldrugs:self_contraceptivedrugs, replace_na, 0)) %>%
  dplyr::select(feid, self_cholesteroldrugs, self_bpdrugs, self_insulin, self_HRT, self_contraceptivedrugs)

data.table::fwrite(data_center_drugs,
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/assessment_center/covariate_drugs.txt",
                   row.names = F, col.names = T, quote = F, sep = "\t")

### extract smoking status from UKB assessment center.
library(dplyr)
library(tidyr)
source("/gdata02/master_data1/UK_Biobank/Data_Tools/extractPheno.R")

raw_data_center_smoke = extractPhenoFromUKBB(fieldID = 20160)

data.table::fwrite(raw_data_center_smoke$fieldData,
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/assessment_center/covariate_smoke.txt",
                   row.names = F, col.names = T, quote = F, sep = "\t")

### extract prostate surgery and thyroidectomy status from UKB assessment center.
library(dplyr)
library(tidyr)
source("/gdata02/master_data1/UK_Biobank/Data_Tools/extractPheno.R")

raw_data_center_thyroidectomy = extractPhenoFromUKBB(fieldID = 41272)

PSA_operations = raw_data_center_thyroidectomy$fieldData %>%
  rowwise() %>%
  mutate(opera = ifelse(any(c_across(f41272_0_0:f41272_0_123) %in% 
                              c("M61", "M611", "M612", "M613", "M614", "M618", "M619",
                                "M62", "M621", "M622", "M623", "M624", "M628", "M629")), 1, 0)) %>%
  select(feid, opera)

data.table::fwrite(PSA_operations,
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/assessment_center/covariate_PSA_operations.txt",
                   row.names = F, col.names = T, quote = F, sep = "\t")


TSH_operations = raw_data_center_thyroidectomy$fieldData %>%
  rowwise() %>%
  mutate(opera = ifelse(any(c_across(f41272_0_0:f41272_0_123) %in% 
                              c("B01", "B011", "B012", "B013", "B014", "B018", "B019",
                                "B02", "B021", "B022", "B023", "B028", "B029",
                                "B08", "B081", "B082", "B083", "B084", "B085", "B086", "B088", "B089",
                                "B09", "B091", "B092", "B098", "B099")), 1, 0)) %>%
  select(feid, opera)

data.table::fwrite(TSH_operations,
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/assessment_center/covariate_TSH_operations.txt",
                   row.names = F, col.names = T, quote = F, sep = "\t")

