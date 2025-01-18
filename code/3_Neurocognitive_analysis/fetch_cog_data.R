################################################################################
# Written by Clément Guichet, PhD Student
# LPNC - CNRS UMR 5105
# 2024

################################################################################

library(tidyverse)
library(readxl)
library(data.table)

rm(list = ls())
setwd("E:/Research_Projects/MEG_CamCAN/TDE_HMM")
################################################################################
# CREATE N x B matrix where N is #subjects and B is #behavioral/neuropsy scores
################################################################################

# Load all subjects'ID
participants <- read_excel("./meta_data/participant_data_T1.xlsx")[, c(1:4, 6:7)] %>%
  mutate(gender_code = case_when(gender_code == 1 ~ -0.5, gender_code == 2 ~ 0.5)) %>% #-0.5 is MALE, 0.5 is FEMALE
  mutate_at(vars(tiv_cubicmm), funs(as.numeric(.))) %>% 
  dplyr::rename(Subj_ID = Subject) %>%
  replace("Subj_ID", seq_len(628))

# Replace missing TIV value by the mean-gender-matched value
participants %>% group_by(gender_code) %>% rstatix::get_summary_stats(tiv_cubicmm, type = "mean_sd")

participants[84,'tiv_cubicmm'] <- 1367718
participants[535,'tiv_cubicmm'] <- 1367718

# Load CAMCAN Cognitive data
CAMCAN_cognitive_data <- read_excel("./meta_data/CognitiveData_CamCAN_Apr2022.xlsx") %>%
  filter(Observations %in% participants$Observations) %>%
  dplyr::rename(Age_Cog = Age) %>%
  dplyr::select(c(
    Observations,
    MMSE,
    Age_Cog,
    Cattell, # Cattell Fluid intelligence
    Proverbs_Summary__Score, # Proverb comprehension (abstraction & EF)
    Picture__Primming_Summary_ACC_baseline_all, # Picture-picture priming (word production)
    TOT_Summary_ToT_ratio # Tip-of-the-tongue
  )) %>%
  dplyr::rename(Proverb = Proverbs_Summary__Score) %>%
  dplyr::rename(Naming = Picture__Primming_Summary_ACC_baseline_all) %>%
  dplyr::rename(ToT_Ratio = TOT_Summary_ToT_ratio) %>%
  mutate_at(vars(MMSE, Age_Cog), funs(as.numeric(.)))


# Load CAMCAN Supplementary Cognitive data
CAMCAN_cognitive_data_supp <- read_excel("./meta_data/CognitiveData_CamCAN_Supplement.xlsx") %>%
  filter(Observations %in% participants$Observations) %>%
  dplyr::select(c(
    Observations,
    Hotel_Task, # EF
    Sentence_Comprehension_c, # Semantic
    Story_Recall, # Memory
    Verbal_Fluency # Language in interaction
  )) %>%
  mutate_at(vars(Hotel_Task, Sentence_Comprehension_c, Story_Recall, Verbal_Fluency), funs(as.numeric(.)))


CAMCAN_CogData_FULL <- merge(CAMCAN_cognitive_data, CAMCAN_cognitive_data_supp, by = "Observations") %>%
  merge(., participants, by = "Observations")

################################################################################
# DATA IMPUTATION

# Exclude subjects with > 3 missing scores
# Impute the median age-decile score for the remaining
# Excluding CC510015 from middle-aged sample
################################################################################

# Exclusion
CAMCAN_CogData_FULL_exclusion <- CAMCAN_CogData_FULL %>%
  mutate(
    E0 = ifelse(is.na(MMSE) | is.nan(MMSE), 1, 0),
    E1 = ifelse(is.na(Cattell) | is.nan(Cattell), 1, 0),
    E2 = ifelse(is.na(Proverb) | is.nan(Cattell), 1, 0),
    E3 = ifelse(is.na(Naming) | is.nan(Cattell), 1, 0),
    E4 = ifelse(is.na(ToT_Ratio) | is.nan(ToT_Ratio), 1, 0),
    E5 = ifelse(is.na(Hotel_Task) | is.nan(Hotel_Task), 1, 0),
    E6 = ifelse(is.na(Sentence_Comprehension_c) | is.nan(Sentence_Comprehension_c), 1, 0),
    E7 = ifelse(is.na(Story_Recall) | is.nan(Story_Recall), 1, 0),
    E8 = ifelse(is.na(Verbal_Fluency) | is.nan(Verbal_Fluency), 1, 0),
    E_tot = E0 + E1 + E2 + E3 + E4 + E5 + E6 + E7 + E8
  ) %>% # Total number of missing scores
  filter(E_tot <= 3) # Discard subjects with more 3 or more



# Age-decile binning
binned <- CAMCAN_CogData_FULL_exclusion %>%
  mutate(Age_Cog_decile = ifelse(Age_Cog <= 29, 25,
                                 ifelse(Age_Cog <= 39, 35,
                                        ifelse(Age_Cog <= 49, 45,
                                               ifelse(Age_Cog <= 59, 55,
                                                      ifelse(Age_Cog <= 69, 65,
                                                             ifelse(Age_Cog <= 79, 75, 85)
                                                      )
                                               )
                                        )
                                 )
  )) %>%
  group_by(Age_Cog_decile, .all = TRUE) %>%
  group_split()

# For each age decile
imputed_list <- list()
for (i in 1:length(binned)) {
  tmp <- rbindlist(lapply(binned[i], as.data.table)) %>% as.data.frame()
  tmp_select <- tmp[, c(2, 4:11)] # Select the cog measures
  # getting median of each column using apply()
  all_column_median <- apply(tmp_select, 2, median, na.rm = TRUE)
  # imputing median value
  for (j in colnames(tmp_select)) {
    tmp_select[, j][is.na(tmp_select[, j])] <- all_column_median[j]
  }
  tmp_imputed <- cbind(
    tmp[, c(12, # Subj_index
            1, # Subj_ID
            3, # Age at time of assessment
            13:16) # Demographics of interest
        ],
    tmp_select
  )
  imputed_list[[i]] <- tmp_imputed
}

CAMCAN_cognitive_data_imputed <- rbindlist(imputed_list)

################################################################################
################################################################################
# Load subjects used for HMM inference
participants_hmm <- read.table("./sflip/subj.txt", quote="\"", comment.char="")

participants_meg <- rio::import("./meta_data/participants_MEG.tsv") %>% 
  dplyr::filter(participant_id %in% participants_hmm$V1) %>% # Keep participants used for training
  dplyr::mutate(dplyr::across(c('participant_id'), substr, 5, nchar(participant_id))) # Removing "sub-"

write.csv(participants_meg, "./output/3_Neurocognitive_analysis/participants_hmm_demographics.csv")

participants_meg %>% rstatix::get_summary_stats(age)
participants_meg %>% group_by(sex) %>% summarise(n = n())

# Adding in Education levels
participants_education <- rio::import("./meta_data/YA_Education_CamCAN.csv")[,-1] %>% 
  dplyr::rename(Observations = SUBJ_ID) %>% 
  na.omit()
# Not keeping physical information because it results in keeping only 243 participants
# participants_physical <- rio::import("./meta_data/epaq.csv")[-1] %>% 
#   na.omit()

cog_data <- CAMCAN_cognitive_data_imputed %>% 
  dplyr::filter(Observations %in% participants_meg$participant_id) %>% # Keep participants with complete cognitive data
  dplyr::select(-c(Subj_ID, age, hand)) %>% # Remove duplicate variables before merging
  dplyr::mutate(ToT_Ratio = 1-ToT_Ratio, # Staying consistent with higher is better
                Hotel_Task = 1-log(Hotel_Task)) # Staying consistent with higher is better / log transformation because reaction times

write.csv(cog_data, "./output/3_Neurocognitive_analysis/cog_data.csv")

cog_data_education <- cog_data %>%
  merge(., participants_education, by = "Observations")  # Keep participants with education information in their younger adulthood

cog_data_education %>% rstatix::get_summary_stats(Age_Cog)
cog_data_education %>% group_by(as.factor(gender_code)) %>% summarise(n = n())
write.csv(cog_data_education, "./output/cog_data_education.csv")


# For Raphael: cog_data_education for 621 MEG particiapnts with NA values for missing education values
# raphael <- cog_data %>%
#   left_join(., participants_education, by = "Observations")  %>% # Keep participants with education information in their younger adulthood
#   left_join((participants_meg %>% rename(Observations = participant_id)) %>% select(Observations), ., by = "Observations")
# 
# write.csv(raphael, "E:/Research_Projects/MEG-MHA/cog_data_education_for_raphael.csv")

