source("01a_Input_Data_Omicron.R")
#reads in all the hospitalisations, deaths testing and cohort data
#need to check the format of the vaccination data in case new codes are added
#this file also accesses 
#source("00_Read_GP_Vaccinations.R") to read in vaccination data

#main output data frames are 
#Vaccinations
head(Vaccinations,1)
#  EAVE_LINKNO vacc_type vacc_type_2 date_vacc_1 date_vacc_2 flag_incon
#1 EAVE0000001        AZ          AZ  2021-04-15  2021-06-11          0
#flag incon=1 for inconsistent vacciantion records
#one row per vaccinated person

#cdw_full  - all testing data
#Positive_Tests - a subset of cdw_full
# date_ecoss_specimen is the same as SpecimenDate in other data frames
# multiple positive tests are here.
#all_deaths - date and cause of death
#all_hospitalisations - date of admission to hospital - no cause
# multiple admissions
#wgs - genomic Sequencing data
#sgene - S gene status information and Ct values


source("01b_get_all_casest.R")
# define the cohort of individuals who tested positive in the study period
# output is
# z_df - endpoints, demographic and clinical information for the tested positive cohort - one row per person
 z_df is saved 

source("02a_analysis.R")
# tabulations, graphs, fitting cox model and sensitivity analyses
# all outputs are saved

#02b_analysis.Rmd
#report file for the analysis in 02a_analysis.R

#################################################################################################

#Below for the Vaccine Effect Analysis using a TND design of individuals symptomatic at the time of test
#Same setup from the 01a_Input_Data_Omicron.R file

source("04a_Get_TND.R")
# define the individuals in the study with the covariates and positive test split by S_Pos and S_Neg as the outcome
# df - demographic and clinical information  - one row per person
# df are saved with the endpoint name
# there is some glm and gam fitting code at the bottom of the file

#04a_Fit_TND_Models.Rmd
#Fit the models, print out odds ratios and VE estimates

#################################################################################################

#Below for the whole population description
#Same setup from the 01a_Input_Data_Omicron.R file

source("05a_Get_Full_Cohort_3.R")
# define the individuals in the whole population with the covariates and positive test split by S_Pos and S_Neg as the outcome
# df_all - demographic and clinical information  - one row per person
# df_all is saved with the endpoint name

source("05a_Full_Cohort_Description.R")
#tables and charts

