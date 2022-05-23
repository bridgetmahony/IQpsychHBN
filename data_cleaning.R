##################################################################
##################################################################
####### IMPORTING AND CLEANING HBN BEHAVIORAL INSTRUMENTS ########
##################################################################
##################################################################

################
#### SET UP ####
################

library(magrittr)
library(dplyr)
library(ggplot2)
library(gplots)
library(forcats)
library(tidyverse)
library(visdat)
library(naniar)


## READ IN DATA ##
setwd('/Users/mahonybw/Desktop/HBN/iqanalyses_rawdata') # path to HBN behavioral data

basic.info <- read.csv('./9994_Basic_Demos_20210310.csv')       # basic info about the visit/subject

cbcl1 <- read.csv('./9994_CBCL_20210310.csv')                   # Child Behavior Checklist item-level and syndrome-scale scores, ages 6-18 
cbcl.desc1 <- read.csv('../modified_data_dictionaries/CBCL_data_dictionary.csv', skip = 1)  # CBCL item descriptions

ysr1 <- read.csv('./9994_YSR_20210310.csv')                     # Youth Self-Report item and syndrome scores, ages 11-17
ysr.desc1 <- read.csv('../modified_data_dictionaries/YSR_data_dictionary.csv', skip = 1)    # YSR descriptions

scared.p1 <- read.csv('./9994_SCARED_P_20210310.csv')           # Screen for Anxiety-Related Disorders item and syndrome scores, parent report, ages 8-18
scared.p.desc1 <- read.csv('../modified_data_dictionaries/SCARED_P_data_dictionary.csv', skip = 1)   # SCARED parent report descriptions, includes clinical cutoffs for syndrome scores
scared.sr1 <- read.csv('./9994_SCARED_SR_20210310.csv')         # SCARED Self-Report, ages 8-18
scared.sr.desc1 <- read.csv('../modified_data_dictionaries/SCARED_SR_data_dictionary.csv', skip = 1)  # SCARED self-report descriptions

conners.sr1 <- read.csv('./9994_C3SR_20210310.csv')             # Conners 3 Self-Report for ADHD items and syndrome scores, ages 8-18
conners.sr.desc1 <- read.csv('../modified_data_dictionaries/C3SR_data_dictionary.csv', skip = 1)  # Conners 3 descriptions

rbs1 <- read.csv('./9994_RBS_20210310.csv')                     # Repetitive Behavior Scale items and syndrome questions, ages 5-21
rbs.desc1 <- read.csv('../modified_data_dictionaries/RBS_Data_dictionary.csv', skip = 1)    # RBS descriptions

srs1 <- read.csv('./9994_SRS_20210310.csv')                     # Social Responsiveness Scale items and syndrome scores, ages 5+
srs.desc1 <- read.csv('../modified_data_dictionaries/SRS_data_dictionary.csv', skip = 1)    # SRS descriptions

sdq1 <- read.csv('./9994_SDQ_20210310.csv')                      # Strengths and Difficulties Questionnaire items and syndrome scores, ages 5+
sdq.desc1 <- read.csv('../modified_data_dictionaries/SDQ_data_dictionary.csv', skip = 1)     # SDQ descriptions

icu.p1 <- read.csv('./9994_ICU_P_20210310.csv')                  # Inventory of Callous-Unemotional Traits items and syndrome scores, parent report, ages 5+
icu.p.desc1 <- read.csv('../modified_data_dictionaries/ICU_P_data_dictionary.csv', skip = 1) # ICU parent-report descriptions
icu.sr1 <- read.csv('./9994_ICU_SR_20210310.csv')                # ICU Self-Report items and syndrome scores, ages 5+
icu.sr.desc1 <- read.csv('../modified_data_dictionaries/ICU_SR_data_dictionary.csv', skip = 1)      # ICU self-report descriptions

# IQ data
wisc1 <- read.csv('./9994_WISC_20210310.csv')                    # most participants, N = 1991, took the WISC
wisc.desc1 <- read.csv('../modified_data_dictionaries/WISC_data_dictionary.csv', skip = 1)   # WISC descriptions


## EDITING/TIDYING DATA ##

# BASIC INFO
basic.info <- basic.info[-1, ] %>%
  mutate(Sex = as.factor(Sex), 
         Age = as.numeric(Age), 
         Participant_Status = as.factor(Participant_Status),
         Enrollment.Year = as.factor(Enrollment.Year), 
         Study.Site = as.factor(Study.Site), 
         Release.Number = as.factor(Release.Number))

# FUNCTION TO TIDY EACH DATA SET #
TIDYING_FXN <- function(x) {
  # assign interim variable name
  intermed <- x
  # coerce columns containing "score" info to numeric
  intermed[12:ncol(x)] <- apply(x[12:ncol(x)], 2, as.numeric)   # score info begins in col 12 for all data sets
  # remove the first row, which contains unnecessary descriptions
  intermed1 <- intermed[-1, ]
  # merge with the age and sex columns from basic.info
  output <- merge(intermed1, basic.info[, c(1, 10, 11)], by.x = "Anonymized.ID") %>%
    # reorder the columns
    select(1, (ncol(x) + 1), ncol(x) + 2, everything())
  # return the output
  return(output)
}

# CBCL
cbcl <- TIDYING_FXN(cbcl1)
# remove rows that do not directly correspond to a column in cbcl, only keep the item text and variable name cols
cbcl.desc <- cbcl.desc1[-c(56, 121, 125), c(1, 2)] %>%
  dplyr::rename(Item.or.Scale = Question) 

# YSR
ysr <-TIDYING_FXN(ysr1)
ysr.desc <- ysr.desc1[-120, c(1, 2)] %>%
  dplyr::rename(Item.or.Scale = Question)

# SCARED - PARENT REPORT
scared.p <- TIDYING_FXN(scared.p1)
scared.p.desc <- scared.p.desc1[-42, ] %>%
  dplyr::rename(Item.or.Scale = Question)
# SCARED - SELF-REPORT
scared.sr <- TIDYING_FXN(scared.sr1) 
scared.sr.desc <- scared.sr.desc1[-42, ] %>%
  dplyr::rename(Item.or.Scale = Question)

# CONNERS - SELF-REPORT
conners.sr <- TIDYING_FXN(conners.sr1)
# need to replace "C3SR" with "CSR" so there are no numbers
colnames(conners.sr)  <- str_replace(colnames(conners.sr), "C3SR", "CSR")
conners.sr.desc <- conners.sr.desc1[-40, c(1, 2)] %>%
  dplyr::rename(Item.or.Scale = Question)
# also replace variable names here
conners.sr.desc$Variable.Name <-  str_replace(conners.sr.desc$Variable.Name, "C3SR", "CSR")

# RBS
rbs <- TIDYING_FXN(rbs1) %>%
  # need to rename some of the columns so they contain the instrument name
  dplyr::rename(RBS_ST = Score_01, RBS_SI = Score_02, RBS_COM = Score_03, 
                RBS_RIT = Score_04, RBS_RI = Score_05, RBS_Tot = Score_Total)
rbs.desc <- rbs.desc1[-45, c(1, 2)] %>%
  dplyr::rename(Item.or.Scale = Question)
# also rename the variables in the description dataset
rbs.desc$Variable.Name[45:50] <- c("RBS_ST", "RBS_SI", "RBS_COMP", "RBS_RIT", "RBS_RI", "RBS_Tot")
rbs.desc$Item.or.Scale[50] <- "Total Score"

# SRS
srs <- srs1
srs <- srs[-1, ]
srs[12:ncol(srs)] <- apply(srs[12:ncol(srs)], 2, as.numeric) 
srs <- merge(basic.info[, c(1, 10, 11)], srs, by.x = "Anonymized.ID") 
srs.desc <- srs.desc1[-66, c(1, 2)] %>%
  dplyr::rename(Item.or.Scale = Question)

# SDQ 
sdq <- TIDYING_FXN(sdq1) %>%
  dplyr::rename(SDQ_CP = SDQ_Conduct_Problems, SDQ_Diff_tot = SDQ_Difficulties_Total, 
                SDQ_EP = SDQ_Emotional_Problems, SDQ_Ext = SDQ_Externalizing, 
                SDQ_GI = SDQ_Generating_Impact, 
                SDQ_HY = SDQ_Hyperactivity, SDQ_Int = SDQ_Internalizing, 
                SDQ_PP = SDQ_Peer_Problems, SDQ_PS_rev = SDQ_Prosocial) 
sdq.desc <- sdq.desc1[-34, c(1, 2)] %>%
  dplyr::rename(Item.or.Scale = Question)
# reverse-score SDQ_Prosocial so that higher score means greater impairment
sdq$SDQ_PS_rev <- 10- sdq$SDQ_PS_rev

# ICU - PARENT REPORT
icu.p <- TIDYING_FXN(icu.p1) %>%
  # shorten some var names
  dplyr::rename(ICU_P_CA = ICU_P_Callous, 
                ICU_P_UC = ICU_P_Uncaring, 
                ICU_P_UE = ICU_P_Unemotional, 
                ICU_P_Tot = ICU_P_Total)
icu.p.desc <- icu.p.desc1[-25, c(1, 2)] %>%
  dplyr::rename(Item.or.Scale = Question) 
# ICU - SELF-REPORT
icu.sr <- TIDYING_FXN(icu.sr1) %>%
  # shorten some var names
  dplyr::rename(ICU_SR_CA = ICU_SR_Callous, 
                ICU_SR_UC = ICU_SR_Uncaring, 
                ICU_SR_UE = ICU_SR_Unemotional, 
                ICU_SR_Tot = ICU_SR_Total)
icu.sr.desc <- icu.sr.desc1[-25, c(1, 2)] %>%
  dplyr::rename(Item.or.Scale = Question)

# IQ VARIABLES
wisc <- TIDYING_FXN(wisc1)
wisc.desc <- wisc.desc1[-3, c(1, 2)] %>%
  dplyr::rename(Item.or.Scale = Question)

## REMOVE ORIGINAL DATA FILES ##
rm(cbcl1, cbcl.desc1, conners.sr.desc1, conners.sr1, icu.p.desc1, icu.p1, icu.sr.desc1, icu.sr1, rbs1, rbs.desc1, scared.p.desc1, 
   scared.p1, scared.sr.desc1, scared.sr1, sdq.desc1, sdq1, srs.desc1, srs1, wisc.desc1, wisc1, ysr.desc1, ysr1)


## COMPILE INTO MASTER DATASET ##
# ALL BEHAVIORAL AND IQ DATA, INCLUDING ITEMS
all_data <- full_join(basic.info, cbcl[, c(1, 14:ncol(cbcl))], by = "Anonymized.ID") %>%
  full_join(., ysr[, c(1, 14:ncol(ysr))], by = "Anonymized.ID") %>%
  full_join(., conners.sr[, c(1, 14:ncol(conners.sr))], by = "Anonymized.ID") %>%
  full_join(., srs[, c(1, 14:ncol(srs))], by = "Anonymized.ID") %>%
  full_join(., scared.p[, c(1, 14:ncol(scared.p))], by = "Anonymized.ID") %>%
  full_join(., scared.sr[, c(1, 14:ncol(scared.sr))], by = "Anonymized.ID") %>%
  full_join(., rbs[, c(1, 14:ncol(rbs))], by = "Anonymized.ID") %>%
  full_join(., sdq[, c(1, 14:ncol(sdq))], by = "Anonymized.ID") %>%
  full_join(., icu.p[, c(1, 14:ncol(icu.p))], by = "Anonymized.ID") %>%
  full_join(., icu.sr[, c(1, 14:ncol(icu.sr))], by = "Anonymized.ID") %>%
  full_join(., wisc[, c(1, 14:ncol(wisc))], by = "Anonymized.ID")
# export
write.csv(all_data, file = "../modified_data/ninth_release_behavioral_data2.csv")

# ALL BEHAVIORAL DATA, SUBSCALES ONLY, T SCORES INSTEAD OF RAW WHERE POSSIBLE,
# AND ONLY KEEPING THOSE PARTICIPANTS FOR WHOM WISC IS COMPLETE
subscales_only <- all_data %>%
  # filter based on the wisc being complete
  filter(WISC_complete == 1) %>%
  select(1:15, !matches("[0123456789]")) %>% # this line removes items
  select(1:15, starts_with(c("CBCL_", "YSR_", "CSR", "SRS_")) & ends_with("_T") | # these lines removes raw scores for subscales that also report t scores
           starts_with(c("SCARED_", "RBS_", "SDQ_", "ICU_")) | 
           WISC_VSI, WISC_VCI, WISC_FRI, WISC_WMI, WISC_PSI, WISC_FSIQ) # this line keeps the wisc subscales of interest
  

# EXPORT THIS DATA FILE
write.csv(subscales_only, file = "../modified_data/ninth_release_iq_analyses2.csv")

# ALL BEHAVIORAL DATA, RAW SUBSCALES ONLY, 
# AND ONLY KEEPING THOSE PARTICIPANTS FOR WHOM WISC IS COMPLETE
raw_subscales_only <- all_data %>%
  # filter based on the wisc being complete
  filter(WISC_complete == 1) %>%
  select(1:15, !matches("[0123456789]")) %>% # removes items
  # remove normed subscales and unwanted wisc subscales
  select(!ends_with("_T") & !starts_with("WISC"), 
         WISC_VSI, WISC_VCI, WISC_FRI, WISC_WMI, WISC_PSI, WISC_FSIQ) %>%
  # remove vars that don't correspond to an actual subscale
  select(!c(CBCL_OP, CBCL_C, YSR_OP, YSR_C, CSR_NI, CSR_PI))

# save
write.csv(raw_subscales_only, file = "../modified_data/ninth_release_rawscores2.csv")








