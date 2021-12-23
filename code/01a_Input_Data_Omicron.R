##########################################################
# Name of file: 01_Input_Data_Omicron.R
# Data release (if applicable):
# Original author(s): Chris Robertson chrisobertson@nhs.net
# Original date: 03 December 2021
# Latest update author (if not using version control) - Chris Robertson chrisobertson@nhs.net
# Latest update date (if not using version control) - 
# Latest update description (if not using version control)
# Type of script: Descriptive stats
# Written/run onL: R Studio SERVER
# Version of R that the script was most recently run on: R 3.6.1
# Description of content: reads in the cohort and merges in the risk groups 
#                         
# Approximate run time: Unknown
##########################################################

# 01 Setup ####
#Libraries
library(plyr)
library(tidyverse)
library(survival)
library(lubridate)
#Load data

Location <- "/conf/"  # Server
#Location <- "//isdsf00d03/"  # Desktop
project_path <- paste0(Location,"EAVE/GPanalysis/progs/CR/SGene_Omicron")
project_path_vaccine <- paste0(Location,"EAVE/GPanalysis/progs/CR/Vaccine")

#just use the demographics from here
EAVE_cohort <- readRDS(paste0(Location,"EAVE/GPanalysis/outputs/temp/Cohort_Demog_Endpoints_Times2021-07-28.rds"))
EAVE_cohort <- filter(EAVE_cohort, !duplicated(EAVE_LINKNO)) %>% dplyr::select(EAVE_LINKNO:ur6_2016_name) %>% 
  mutate(ageYear=ageYear+1)  # make age as at March 2021

a_begin <- as.Date("2021-11-01")  #start date for the analysis
#read in all deaths and then omit any who have died prior to a_begin
all_deaths  <- readRDS(paste0(Location,"EAVE/GPanalysis/data/all_deaths.rds"))
summary(all_deaths)
EAVE_cohort <- EAVE_cohort %>% left_join(dplyr::select(all_deaths, EAVE_LINKNO, NRS.Date.Death), by="EAVE_LINKNO")
EAVE_cohort <- EAVE_cohort %>% filter(is.na(NRS.Date.Death) | !is.na(NRS.Date.Death) & NRS.Date.Death <= a_begin)

EAVE_Weights <- readRDS(paste0(Location,"EAVE/GPanalysis/outputs/temp/CR_Cohort_Weights.rds"))
EAVE_cohort  <- EAVE_cohort %>% left_join(EAVE_Weights, by="EAVE_LINKNO")
EAVE_cohort$eave_weight[is.na(EAVE_cohort$eave_weight)] <- mean(EAVE_cohort$eave_weight, na.rm=T)

#read in the Previous Tests data
cdw_full  <- readRDS(paste0(Location,"EAVE/GPanalysis/data/CDW_full.rds"))
cdw_full <- cdw_full %>% mutate(date_ecoss_specimen = as_date(date_ecoss_specimen))
cdw_full <- cdw_full %>%  
  arrange(EAVE_LINKNO, date_ecoss_specimen, desc(test_result)) %>%
  filter(!duplicated(paste(EAVE_LINKNO, date_ecoss_specimen)))  #get one test per person per day - preferentially positive test 
cdw_full <- filter(cdw_full, date_ecoss_specimen <= Sys.Date())  
summary(cdw_full)

Positive_Tests <- cdw_full %>% filter(test_result=="POSITIVE") %>% 
  dplyr::select(EAVE_LINKNO, test_id, date_ecoss_specimen) 
summary(Positive_Tests)
length(unique(paste(Positive_Tests$EAVE_LINKNO, Positive_Tests$date_ecoss_specimen))) # no duplicates on the same date

#read in the GP vaccination data
source("../00_Read_GP_Vaccinations.R")



#get covid death certificate deaths
z <- all_deaths %>%  
  mutate(across(UNDERLYING_CAUSE_OF_DEATH:CAUSE_OF_DEATH_CODE_9, ~if_else(. %in% c("U071","U072"), 1,0)))
z <- z %>% rowwise() %>% mutate(rowsum = sum(c_across(UNDERLYING_CAUSE_OF_DEATH:CAUSE_OF_DEATH_CODE_9))) %>% 
  mutate(covid_death_cert = if_else(rowsum>=1,1,0)) 
all_deaths <- z %>% dplyr::select(-rowsum)
summary(dplyr::select(all_deaths, NRS.Date.Death, NRS.Reg.Date, covid_death_cert))


all_hospitalisations  <- readRDS(paste0(Location,"EAVE/GPanalysis/data/automated_any_hospitalisation_post_01022020.rds"))
summary(all_hospitalisations)

#cohort + risk groups
rg <- readRDS(paste0(project_path_vaccine,"/output/temp/Qcovid_all.rds"))
rg <- rg %>% dplyr::select(-(Sex:ur6_2016_name), -Q_BMI)
rg <- filter(rg, !duplicated(EAVE_LINKNO))

z <- readRDS(paste0(Location,"EAVE/GPanalysis/outputs/temp/CR_Cohort_RG_EAVE_BP_Smoke.rds"))
z <- filter(z, !duplicated(EAVE_LINKNO))
z <- z %>% dplyr::select(EAVE_LINKNO, EAVE_Smoking_Status_Worst, EAVE_BP) %>% 
  dplyr::rename(EAVE_Smoke = EAVE_Smoking_Status_Worst)

rg <- rg %>% left_join(z, by="EAVE_LINKNO")
rg <- rg %>% mutate(EAVE_Smoke = if_else(!is.na(EAVE_Smoke), as.character(EAVE_Smoke), "Unknown"),
                    EAVE_BP = if_else(!is.na(EAVE_BP), as.character(EAVE_BP), "No Investigation"))

#update weights
#those with a pis records over the last 12 months before March 2020
bnf <- readRDS(paste0(Location,"EAVE/GPanalysis/data/BNF_paragraphs.rds"))
z_ids <- c(Vaccinations$EAVE_LINKNO, all_deaths$EAVE_LINKNO,  bnf$EAVE_LINKNO,
           cdw_full$EAVE_LINKNO, all_hospitalisations$EAVE_LINKNO) %>% unique()
#summary(filter(EAVE_cohort, !(EAVE_LINKNO %in% z_ids))$eave_weight)
z_N <- round(sum(EAVE_cohort$eave_weight) )
z_k <- sum(EAVE_cohort$EAVE_LINKNO %in% z_ids)
z_m <- round(sum(filter(EAVE_cohort, (EAVE_LINKNO %in% z_ids))$eave_weight))
z <- EAVE_cohort %>% mutate(ew = if_else(EAVE_LINKNO %in% z_ids, 1, eave_weight*(z_N - z_k)/(z_N - z_m)) )
EAVE_cohort <- z %>% dplyr::select(-eave_weight) %>% dplyr::rename(eave_weight=ew)

z <- read_csv(paste0(Location,"/EAVE/GPanalysis/data/restored/map_files/Datazone2011Lookup.csv")) %>% 
  dplyr::select(DataZone, InterZone, Council, HB)
EAVE_cohort <- EAVE_cohort %>% left_join(z, by="DataZone") %>% 
  mutate(HB = if_else(is.na(HB),"Unknown", HB),
         InterZone = if_else(is.na(InterZone),"Unknown", InterZone),
         Council = if_else(is.na(Council),"Unknown", Council))


wgs <- readRDS(paste0(Location,"EAVE/GPanalysis/data/WGS_latest.rds")) %>% 
  mutate_at(c("Collection_Date","Sequencing_Date","Alignment_Date"), ~ as.Date(. , format="%d/%m/%Y")) %>% 
  filter(Collection_Date >= a_begin)

#z_id <- filter(wgs, duplicated(EAVE_LINKNO)) %>% pull(EAVE_LINKNO)
#z <- filter(wgs, EAVE_LINKNO %in% z_id) %>% arrange(EAVE_LINKNO)
#duplicates all have the same lineage
wgs <- wgs %>% arrange(EAVE_LINKNO, Collection_Date) %>% filter(!duplicated(EAVE_LINKNO))
a_end_wgs <- max(wgs$Collection_Date) - 2  #table(wgs$Collection_Date)  #check each time

sgene <- readRDS(paste0(Location,"/EAVE/GPanalysis/data/omicron_ctvals.rds"))
summary(sgene)  

shielding <- readRDS(paste0(Location,"EAVE/GPanalysis/data/Shielding_list.rds"))
immuno <- readRDS(paste0(Location,"EAVE/GPanalysis/data/cleaned_data/Imm_supp_cohort_Nov2021.rds"))
