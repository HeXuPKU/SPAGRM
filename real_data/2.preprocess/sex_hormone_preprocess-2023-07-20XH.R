
###### pre-clean blood protein phenotype gp_clinical table.
library(dplyr)
library(tidyr)

### read in the blood FSH, LH, E2, T phecode.
source("/gdata01/user/xuhe/family_relatedness/real_data-2022-12-29/1.extract_pheno/phecode-2023-05-10XH.R")

### load the cleaned gp_clinical table.
gp_clinical = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/cleaned_gp_clinical.txt")

# follow https://www.nature.com/articles/ejhg2015102

### extract TSH from gp_clinical data.
FSH = gp_clinical %>%
  filter(grepl(FSH_codes, code)) %>% 
  mutate(value1_numeric = as.numeric(value1), value2_numeric = as.numeric(value2), value3_numeric = as.numeric(value3)) %>%
  mutate(fsh = coalesce(value1_numeric, value2_numeric, value3_numeric)) %>%
  filter(fsh > 1 & fsh < 200) %>% 
  mutate(value3 = toupper(value3)) %>%
  filter(value3 %in% c("", "MEA000", "MEA149", "MEA160", "MEA170", "MICROU/L", "MIU/L", "MU/L")) %>%
  group_by(eid) %>%
  mutate(mi = n()) %>%
  mutate(mean = mean(fsh)) %>%
  mutate(sd = ifelse(mi == 1, 0, 
                     ifelse(sd(fsh) > 2, 2, sd(fsh)))) %>%
  mutate(lower_tail = ifelse(mi == 1, mean, mean - 4 * sd),
         higher_tail = ifelse(mi == 1, mean, mean + 4 * sd)) %>%
  filter(fsh >= lower_tail & fsh <= higher_tail) %>%
  group_by(eid, event_dt) %>%
  mutate(rep = n()) %>%
  mutate(fsh_new = ifelse(rep == 1, fsh, 
                          ifelse(max(fsh) - min(fsh) <= 0.2, mean(fsh), 0))) %>%
  filter(fsh_new > 0) %>%
  mutate(fsh_new = round(fsh_new, 2)) %>%
  select(eid, data_provider, event_dt, fsh_new) %>% 
  dplyr::rename(fsh = fsh_new) %>%
  distinct(eid, event_dt, .keep_all = TRUE)

### write the pre-cleaned FSH into file.
data.table::fwrite(FSH %>% arrange(eid, event_dt), 
                   file = "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/preclean_pheno/FSH.txt", 
                   row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
