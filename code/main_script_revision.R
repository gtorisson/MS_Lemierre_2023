library(tidyverse) # data wrangling
library(readxl) # read in .xlsx files
library(stringi) # string manipulations
library(ggsci) # for plots in ggplot2
library(here) # cleaner file paths
library(UniprotR) # to retrieve protein names from Uniprot IDs

# setting the accepted fraction of missing per sample
missing_accepted = 0.3

# sourcing local functions
source(here("code", "functions.R"))

# file import of raw data
raw_data <- read_excel(here("data", "20220629_Karin_Holm_linkopingLib_Hong_Report_rawdata_nygruppering_avidentifierat.xlsx"), 2) %>%
  select(starts_with("PG") | starts_with("LS") | starts_with( "Sepsis"))



# data is in wide format in the original data
df_wide <-  raw_data %>%
  # removing double names (separated by semicolons) in protein name
  mutate(PG.ProteinGroups = stri_extract_first_regex(PG.ProteinGroups, "^[^;]+")) %>%
  dplyr::rename(protein = PG.ProteinGroups) %>%
  # setting NA if filtered (low signal)
  mutate_all(~na_if(., "Filtered")) %>%
  # making readings from MS numeric (imported as character)
  mutate(across(2:24,~as.numeric(.))) %>%
  # setting NA if NaN 
  mutate_all(~ifelse(is.nan(.), NA, .))

# making a long data frame 
df_long <- pivot_longer(df_wide, cols = c(2:24), names_to = c("group", "id"),
                        names_sep = -2, values_to = "intensity") %>% 
  # tidying the names from ... created by pivoting
  mutate_at(vars(2,3), ~str_remove_all(.,"[.]")) %>%
  # renaming for clarity
  mutate(group = ifelse(group == "NonLSthrombosis", "Other_thrombosis",group),
         # factorising
         group = factor(group),
         # naming individual samples
         sample = factor(str_c(group, id, sep = "_"))) %>%
  # selecting only relevant variables
  select(protein, group,sample, intensity)


# LOG2 TRANSFORMATION - ADDING A LOGGED VARIABLE
df_long <- df_long %>% mutate(log2_intensity = log2(intensity))

# FILTERING 
## filtering away those with > 30% NA
df_long_filt <- df_long %>% group_by(protein) %>% 
  filter(sum(is.na(intensity)) < missing_accepted*length(levels(df_long$sample)))

# NORMALISATION 
# normalisation using median subtraction per sample. Adding a normalised variable
df_long_norm <- df_long_filt %>% group_by (sample) %>% 
  mutate(norm_log2_intensity = log2_intensity - median(log2_intensity, na.rm = T))

# IMPUTATION
# the imputation - using perseus defaults
df_long_imp <- impute_as_perseus(df_long_norm, sample = "sample", value = "norm_log2_intensity")


## DIFF EXPRESSION
#  LS vs Sepsis
table1 <- diff_expression(df_long_imp, protein_names = "protein", test_value = "norm_log2_intensity_imputed", 
                          group_var = "group", grp1 = "LS", grp2 = "Sepsis")










