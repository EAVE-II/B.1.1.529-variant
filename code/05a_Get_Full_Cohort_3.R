##########################################################
# Name of file: 05a_Get_Full_Cohort.R
# Data release (if applicable):
# Original author(s): Chris Robertson chris.robertson@phs.scot
# Original date: 12 Dec 2021
# Latest update author (if not using version control) - Chris Robertson chris.robertson@phs.scot
# Latest update date (if not using version control) - 
# Latest update description (if not using version control)
# Type of script: Descriptive stats
# Written/run onL: R Studio SERVER
# Version of R that the script was most recently run on: R 3.6.1
# Description of content: reads in the cohort and merges in vaccination data 
# Approximate run time: Unknown
#
#   run 01a_Data_Omicron###############################

library(survival)

output_list <- list()

a_begin <- as.Date("2021-11-15")  #beginning of omicron period
output_list$a_begin <- a_begin
#output_list$dose_2_periods <- 15

z_df <- Vaccinations 
z_df <- z_df %>% right_join(dplyr::select(EAVE_cohort, EAVE_LINKNO), by="EAVE_LINKNO")  #EAVE_cohort is the whole population

#link in deaths and remove any who died before vaccination started
z_df <- z_df %>% left_join(dplyr::select(all_deaths, EAVE_LINKNO, NRS.Date.Death), by="EAVE_LINKNO")
z_df <- z_df %>% filter(NRS.Date.Death > a_begin | is.na(NRS.Date.Death))

#covid_hospitalisations is unique

output_list$endpoint <- "covid_pos"  # "covid_death", "covid_hosp", "covid_hosp_death"

if (output_list$endpoint == "covid_hosp") {z_event <- covid_hospitalisations
  z_event <- z_event %>% dplyr::rename(event_date=admission_date) 
}
if (output_list$endpoint == "covid_hosp_death") {z_event <- covid_hosp_death
  z_event <- z_event %>% dplyr::rename(event_date=admission_date) 
}
if (output_list$endpoint == "covid_pos") {
  z_event <- Positive_Tests %>% filter(date_ecoss_specimen >= a_begin) %>% 
    arrange(EAVE_LINKNO, date_ecoss_specimen) %>% 
    filter(!duplicated(EAVE_LINKNO)) %>% 
    dplyr::rename(event_date=date_ecoss_specimen)
}


print(nrow(z_event))
print(table(z_event$EAVE_LINKNO %in% EAVE_cohort$EAVE_LINKNO))

a_end <- max(z_event$event_date)
output_list$a_end <- a_end

z <- z_df %>% left_join(z_event, by="EAVE_LINKNO")
if (output_list$endpoint %in% c("covid_death", "covid_hosp", "covid_hosp_death")) {
  z <- z %>% filter(is.na(event_date) | event_date > a_begin) #omit those with event before vaccination
  }  # don't omit those with a previous positive test
z <- z %>% mutate(event = if_else(is.na(event_date), 0L,1L)) %>% 
  mutate(event_date = if_else(is.na(event_date), a_end ,event_date))
print(sum(z$event))

#change event date to NRS.Date.Death for any who died before the study end and who did not have an event
z1 <- z %>% mutate(event_date = if_else(!is.na(NRS.Date.Death) & event==0 & NRS.Date.Death < event_date, NRS.Date.Death, event_date))

df_all <- z1 #keep df_all - this has all subjects
#add in the covariates to df_all for descriptives of the cohort
df_all <- df_all %>% left_join(dplyr::select(EAVE_cohort, - NRS.Date.Death), by="EAVE_LINKNO")
df_all <- df_all %>%  filter(!is.na(ageYear))  # drop those who do not link into EAVE
#add in the covariates to df_all for descriptives of the cohort
df_all <- df_all %>% left_join(rg, by="EAVE_LINKNO")
#df_all <- df_all %>%  mutate(days_between_vacc = as.numeric(date_vacc_2-date_vacc_1)) %>% 
#  mutate(vacc_gap = cut(days_between_vacc, breaks=c(1,48, 62, 76, 90, max(days_between_vacc)),
#                                     labels=c("<7wk","7-8 wk","9-10 wk","11-12 wk", "13+ wk")))
#get positive test before a_begin variable
z <- df_all %>% dplyr::select(EAVE_LINKNO)
z <- z %>% left_join(Positive_Tests, by="EAVE_LINKNO") %>% filter(!is.na(date_ecoss_specimen)) %>% 
   filter(date_ecoss_specimen < a_begin) #keep the samples before a_begin
z <- z %>% arrange(EAVE_LINKNO, desc(date_ecoss_specimen)) %>% 
  filter(!duplicated(EAVE_LINKNO)) %>% dplyr::select(-test_id) %>% 
  rename(date_ecoss_specimen_prior = date_ecoss_specimen)
#now link back
df_all <- df_all %>% left_join(z, by="EAVE_LINKNO") 
#group days from previous positive test to a_begin
z1 <- df_all %>% 
  mutate(days_from_previous_pos = as.numeric(a_begin - date_ecoss_specimen_prior))  # all previous tests before a_begin
z1 <- z1 %>% mutate(pos_before_start = case_when(is.na(days_from_previous_pos) ~ "not_prev_pos",
                                                days_from_previous_pos >= 1 & days_from_previous_pos <= 28 ~ "pos_1:28",
                                                days_from_previous_pos >= 29 & days_from_previous_pos <= 90 ~ "pos_29:90",
                                                TRUE ~ "pos_91+")) 
df_all <- z1
  

#get number of tests before vaccination variable
z <- df_all %>% dplyr::select(EAVE_LINKNO)
z <- z %>% left_join(dplyr::select(cdw_full, EAVE_LINKNO, date_ecoss_specimen), by="EAVE_LINKNO") %>%
  filter(!is.na(date_ecoss_specimen)) %>% #omit those never tested
  filter(date_ecoss_specimen < a_begin) #keep those sample before a_begin
z <- z %>% group_by(EAVE_LINKNO) %>% dplyr::summarise(N=n()) 
#now link back
df_all <- df_all %>% left_join(z, by="EAVE_LINKNO") %>% 
  mutate(N=if_else(is.na(N), 0L,N))
df_all <- df_all %>% mutate(n_tests_gp = cut(N, breaks=c(-1,0,1,2,3,4,9,19, max(N)), labels=c("0","1","2","3","4","5-9","10-19","20+")))
df_all$N <- NULL 

#bmi  - DONT USE FOR CHILDREN
df_all <- df_all %>% mutate(bmi.gp = cut(bmi_impute, breaks=c(-1, 20,25,30,35,40,51), labels=FALSE))
df_all <- df_all %>% mutate(bmi_gp = case_when(ageYear <= 17 ~ "NA_too_young",
                                           ageYear >= 18 & !is.na(bmi.gp) ~ as.character(bmi.gp),
                                           TRUE ~ "Unknown")) %>% 
  mutate(bmi_gp = factor(bmi_gp, levels=c("1","2","3","4","5","6","NA_too_young"),
                         labels=c("<20","20-24","25-29","30-34","35-39","40+","NA_too_young")))
df_all$bmi.gp <- NULL


df_all <- df_all %>% mutate(Q_DIAG_DIABETES = if_else(Q_DIAG_DIABETES_1==1,Q_DIAG_DIABETES_1,Q_DIAG_DIABETES_2 ),
                            Q_DIAG_CKD = if_else(Q_DIAG_CKD_LEVEL >=1, 1, Q_DIAG_CKD_LEVEL))

#df_all$event_date <- as.Date(df_all$event_date, origin=as.Date("1970-01-01"))
print(sum(df_all$event))

df_all <- df_all %>% mutate(age_gp = cut(ageYear, breaks=c(seq(-1,84, by=5),max(ageYear)), 
                                     labels=c(paste0(seq(0,80, by=5),"-",seq(4,84, by=5)),"85+")))



#vaccination status at a_begin
z1 <- df_all %>% 
  filter(is.na(flag_incon) | flag_incon==0) %>% #omit those vaccinated with inconsistent records
  #correct missing date_vacc_2 - need to clean up Vaccinations - chaked and the vacc 3 looks OK
  mutate(date_vacc_2 = if_else(!is.na(date_vacc_2) & !is.na(date_vacc_3) & date_vacc_2==date_vacc_3, date_vacc_1+77, date_vacc_2)) %>%
  filter(!(!is.na(date_vacc_2) & !is.na(date_vacc_3) & date_vacc_3 <= date_vacc_2 + 28)) %>% 
  mutate(day_1 = as.numeric(a_begin - date_vacc_1),
         day_2 = as.numeric(a_begin - date_vacc_2),
         day_3 = as.numeric(a_begin - date_vacc_3)) 
z1 <- z1 %>% mutate(day_1 = if_else(!is.na(day_1) & day_1 <= 0, NA_real_, day_1),
                    day_2 = if_else(!is.na(day_2) & day_2 <= 0, NA_real_, day_2),
                    day_3 = if_else(!is.na(day_3) & day_3 <= 0, NA_real_, day_3))
z1 <- z1 %>%   mutate(vacc_1_gp = cut(day_1, breaks= c( 0, 27, max(day_1, na.rm=T)),labels=FALSE),
                      vacc_2_gp = cut(day_2, breaks= c(0, 13, 69, 104, 139, 174, max(day_2, na.rm=T)),labels=FALSE),
                      vacc_3_gp = cut(day_3, breaks= c( 0, 7, 13, max(day_3, na.rm=T)),labels=FALSE))
z1 <- z1 %>% mutate(vacc_1_gp = case_when(is.na(vacc_1_gp) ~ "v1_uv",
                                          vacc_1_gp==1 ~ "v1_0:3",
                                          vacc_1_gp==2 ~ "v1_4+") )
z1 <- z1 %>% mutate(vacc_2_gp = case_when(is.na(vacc_2_gp) ~ "v2_uv",
                                          vacc_2_gp==1 ~ "v2_0:1",
                                          vacc_2_gp==2 ~ "v2_2:9",
                                          vacc_2_gp==3 ~ "v2_10:14",
                                          vacc_2_gp==4 ~ "v2_15:19",
                                          vacc_2_gp==5 ~ "v2_20:24",
                                          vacc_2_gp==6 ~ "v2_25+") )
z1 <- z1 %>% mutate(vacc_3_gp = case_when(is.na(vacc_3_gp) ~ "v3_uv",
                                          vacc_3_gp==1 ~ "v3_0",
                                          vacc_3_gp==2 ~ "v3_1",
                                          vacc_3_gp==3 ~ "v3_2+") )

z1 <- z1 %>% mutate(vs=case_when(vacc_3_gp != "v3_uv" ~ vacc_3_gp,
                                 vacc_2_gp != "v2_uv" ~ vacc_2_gp,
                                 TRUE ~ vacc_1_gp)) %>% 
  mutate(vs = if_else(vs=="v1_uv", "uv",vs))
#z_labs <- sort(unique(z1$vs))
z_labs <- c("uv" ,"v1_0:3", "v1_4+", "v2_0:1", "v2_2:9", "v2_10:14", "v2_15:19", "v2_20:24", "v2_25+", "v3_0", "v3_1", "v3_2+" )
z1 <- z1 %>% mutate(vs=factor(vs, levels=z_labs))

z1 %>% group_by(vs) %>% dplyr::summarise(N=round(sum(eave_weight)))

df_all <- z1

z_df <- df_all %>% 
  mutate(vacc_type = if_else(is.na(vacc_type) , "uv", as.character(vacc_type)) ) %>% 
  filter(vacc_type != "UNK") %>% 
  mutate(vacc_type = factor(vacc_type, levels=c("uv","AZ","Mo","PB") ) ) %>% 
  mutate(vt = fct_cross(vs,vacc_type, sep="_")) %>% 
  mutate(vt = fct_recode(vt, "uv" ="uv_uv", "uv" = "uv_AZ", "uv" = "uv_Mo","uv" = "uv_PB"))
df_all <- z_df

#add in s-gene status
#read in the all positive cases data
#don't use the date filter for sgene as that removes those who were post a_being but also tested positiv in period after Nov 01
#correct those with an S_gene value but negative after
z_sg <- readRDS("./output/temp/Sgene_all_positive.rds") %>% dplyr::select(EAVE_LINKNO, s_gene, date_ecoss_specimen)
z_sg <- z_sg %>% #filter(date_ecoss_specimen >= a_begin) %>%   #EAVE_LINKNO is missing for about 10% of sgene - there are a few cases where the date ecoss specimen is not the same (273/59357)
  dplyr::select(-date_ecoss_specimen)
z <- df_all %>%  left_join(z_sg, by="EAVE_LINKNO")
z <- z %>%  mutate(s_gene = fct_explicit_na(s_gene, "Not_Done"))
z_levels <- levels(z$s_gene)
z <- z %>%  mutate(s_gene = case_when(event==0 ~ factor("Not_Done",levels=z_levels),
                            TRUE ~ s_gene))
table(z$event, z$s_gene, exclude=NULL)
df_all <- z

df_all <- df_all %>% mutate(shielding = if_else(EAVE_LINKNO %in% shielding$EAVE_LINKNO, 1, 0))
df_all <- df_all %>% mutate(immuno = if_else(EAVE_LINKNO %in% immuno$EAVE_LINKNO, 1, 0))

#get the endpoints
df_all <- df_all %>%  mutate(event_s_neg = if_else(event==1 & !is.na(s_gene) & s_gene=="S_Neg", 1,0),
                             event_s_pos = if_else(event==1 & !is.na(s_gene) & s_gene=="S_Pos", 1,0))

print(table(df_all$event, df_all$s_gene, exclude=NULL))
print(table(df_all$event_s_pos, df_all$s_gene, exclude=NULL))
print(table(df_all$event_s_neg, df_all$s_gene, exclude=NULL))

saveRDS(df_all, paste0("./output/temp/df_all_full_",output_list$endpoint,".RDS"))


#sample cohort - z_n_controls_per_event non events per event
event_ids <- filter(df_all, event==1 & event_date > a_begin)%>% pull(EAVE_LINKNO) %>% unique()
output_list$seed <- 21021954
set.seed(output_list$seed)
if (output_list$endpoint == "covid_pos") output_list$controls_per_event <- 5 else output_list$controls_per_event <- 10 

no_event_ids <- filter(df_all, event==0 & event_date > a_begin) %>% 
  filter(if (output_list$endpoint == "covid_pos") pos_before_start != "pos_1:28" else TRUE) %>% 
  slice_sample(n=output_list$controls_per_event*length(event_ids)) %>%
  pull(EAVE_LINKNO) %>% unique()
z_df <- df_all %>% dplyr::select(EAVE_LINKNO, date_vacc_1, date_vacc_2, date_vacc_3, event_date, event) %>% #keep minimal columns
  filter(EAVE_LINKNO %in% c(event_ids, no_event_ids))
print(sum(z_df$event))

#z_df <- slice_sample(z_df, n=10000)

#add in deaths and change event_date
#have checked no event_dates after date death
#this bit just checks and does no changes to the z_df data frame
z <- z_df %>% left_join(all_deaths, by="EAVE_LINKNO") %>% 
  mutate(event_date = if_else(!is.na(NRS.Date.Death) & (NRS.Date.Death < event_date), NRS.Date.Death, event_date))
print(sum(z$event))
print(summary(z$event_date))
print(summary(z_df$event_date))
print(summary(z))
table(z$event_date > z$NRS.Date.Death, exclude=NULL)

#set up the survival times
z <- tmerge(z_df,z_df,id=EAVE_LINKNO, endpt = event(event_date, event), tstart=a_begin, tstop=event_date)

#get the post vacc 1 periods
z_df <- z_df %>% mutate(period_end_date = date_vacc_1)
z <- tmerge(z,z_df, id=EAVE_LINKNO, per1=tdc(period_end_date))
names(z)[names(z)=="per1"] <- "pv_uv"

#get the post vacc 1 periods
z_df <- z_df %>% mutate(period_end_date = date_vacc_1 + 27)
z <- tmerge(z,z_df, id=EAVE_LINKNO, per1=tdc(period_end_date))
names(z)[names(z)=="per1"] <- "pv_v1_1"

z_df <- z_df %>% mutate(period_end_date = date_vacc_2)
z <- tmerge(z,z_df, id=EAVE_LINKNO, per1=tdc(period_end_date))
names(z)[names(z)=="per1"] <- "pv_v1_2"

#get the post vaccination periods - they do not have to be the same length
z_df <- z_df %>% mutate(period_end_date = date_vacc_2 + 13)
z <- tmerge(z,z_df, id=EAVE_LINKNO, per1=tdc(period_end_date))
names(z)[names(z)=="per1"] <- "pv_v2_1"

z_df <- z_df %>% mutate(period_end_date = date_vacc_2 + 69)
z <- tmerge(z,z_df, id=EAVE_LINKNO, per1=tdc(period_end_date))
names(z)[names(z)=="per1"] <- "pv_v2_2"

z_df <- z_df %>% mutate(period_end_date = date_vacc_2 + 104)
z <- tmerge(z,z_df, id=EAVE_LINKNO, per1=tdc(period_end_date))
names(z)[names(z)=="per1"] <- "pv_v2_3"

z_df <- z_df %>% mutate(period_end_date = date_vacc_2 + 139)
z <- tmerge(z,z_df, id=EAVE_LINKNO, per1=tdc(period_end_date))
names(z)[names(z)=="per1"] <- "pv_v2_4"

z_df <- z_df %>% mutate(period_end_date = date_vacc_2 + 174)
z <- tmerge(z,z_df, id=EAVE_LINKNO, per1=tdc(period_end_date))
names(z)[names(z)=="per1"] <- "pv_v2_5"

z_df <- z_df %>% mutate(period_end_date = date_vacc_3)
z <- tmerge(z,z_df, id=EAVE_LINKNO, per1=tdc(period_end_date))
names(z)[names(z)=="per1"] <- "pv_v2_6"

#get the post vacc 3 periods
z_df <- z_df %>% mutate(period_end_date = date_vacc_3 + 7)
z <- tmerge(z,z_df, id=EAVE_LINKNO, per1=tdc(period_end_date))
names(z)[names(z)=="per1"] <- "pv_v3_1"

z_df <- z_df %>% mutate(period_end_date = date_vacc_3 + 13)
z <- tmerge(z,z_df, id=EAVE_LINKNO, per1=tdc(period_end_date))
names(z)[names(z)=="per1"] <- "pv_v3_2"


#get the temporal periods
z_period_length <- 7
z_n_periods <- trunc((as.numeric(a_end-a_begin)/z_period_length))
for (i in 1:z_n_periods){
  #i <- 1
  print(i)
  z_df <- z_df %>% mutate(period_end_date = a_begin + z_period_length*i)
  z <- tmerge(z,z_df, id=EAVE_LINKNO, per1=tdc(period_end_date))
  names(z)[names(z)=="per1"] <- paste0("p_",i)
}

z_names <- names(z)[grepl("p_", names(z))]
z1 <- z %>% mutate(period = apply(z[,z_names], 1, sum))
z1 <- z1 %>% dplyr::select(-all_of(z_names))
z_names <- names(z)[grepl("pv_", names(z))]
#z_names <- z_names[!grepl("pv_v3",z_names)] # omit dose 3
z1 <- z1 %>% mutate(pv_period = apply(z[,z_names], 1, sum)) 
z1 <- z1 %>% dplyr::select(-all_of(z_names))
#z_names <- names(z)[grepl("pv_v3", names(z))]
#z1 <- z1 %>% mutate(pv_v3_period = apply(z[,z_names], 1, sum)) 
#z1 <- z1 %>% dplyr::select(-all_of(z_names))
z <- z1

#having got the time dependent periods move the data back to z_df
#z_df in long format
#endpt is the response variable; the origianl one event is duplicated when the periods and pv_periods are added
#don't use event.
#tstart and tstop are the start and stop dates of the intervals
z_df <- z
print(sum(z_df$endpt))

#checking
#z <- z_df %>% filter(!duplicated(EAVE_LINKNO)) %>% left_join(dplyr::select(df_all, EAVE_LINKNO, vs), by="EAVE_LINKNO")
#z1 <- filter(z, vs=="v1_0:3")

#z_df <- z_df %>% mutate(tstop=as.Date(tstop, origin=as.Date("1970-01-01")))
#calculate the interval lengths 
z_df <- z_df %>% mutate(pyears = as.numeric(tstop-tstart))

#get the dates for the beginning of the periods
z <- z_df %>% group_by(period) %>% 
  dplyr::summarise(date_period_begin = min(tstart), n_events =sum(endpt), pyears=sum(pyears))

#merge in the covariates and calculate other covariates
df <- z_df %>% 
  left_join(dplyr::select(df_all, -date_vacc_1, -date_vacc_2, -vacc_type_2, -event_date, -event), by="EAVE_LINKNO") %>%
  filter(!is.na(ageYear)) # %>%  #drop those who do not link into EAVE

z_levs <- sort(unique(df$pv_period))
z_labs <- levels(df$vs)
#make period and pv_period factors
z_df <- df %>%
  mutate(period_f = factor(period, levels=0:max(period), labels=z$date_period_begin, ordered=FALSE),
         pv_period_f = factor(pv_period, levels=z_levs,labels=z_labs))

df <- z_df  #df is the working data set for the analysis

#calculate the weights for the sampling - event is correct here
df <- df %>% mutate(weight = if_else(event==1,1,nrow(df_all)/(output_list$controls_per_event*length(event_ids)) ) )  
#merge the sampling weights and eave_weights
df <- df %>% mutate(ew = eave_weight*weight)


variables_hosp <- c("age_gp" , "Sex", "simd2020_sc_quintile", "ur6_2016_name", "pos_before_start", "n_tests_gp", "bmi_gp",
                    "Q_DIAG_AF", "Q_DIAG_ASTHMA", "Q_DIAG_CCF", "Q_DIAG_CHD", "Q_DIAG_COPD" , "Q_DIAG_DIABETES", "Q_DIAG_EPILEPSY", 
                    "Q_DIAG_FRACTURE", "Q_DIAG_PVD", "Q_DIAG_RA_SLE", "Q_DIAG_SEV_MENT_ILL", "Q_DIAG_STROKE", "Q_DIAG_VTE", "Q_DIAG_CKD",
                    "EAVE_Smoke", "EAVE_BP", "shielding","immuno", "n_oth_risk_gps")

z_vars <- names(df_all)[grepl("Q_DIAG", names(df_all))]
z_vars <- z_vars[!(z_vars %in% variables_hosp)]
z_vars <- z_vars[!(z_vars %in% c("Q_DIAG_DIABETES_1"  ,  "Q_DIAG_DIABETES_2" , "Q_DIAG_CKD_LEVEL"))]
z_vars <- c(z_vars, "Q_HOME_CAT","Q_LEARN_CAT")

z <- df_all %>% dplyr::select_at(c("EAVE_LINKNO", z_vars)) %>% 
  mutate(Q_HOME_CAT = if_else(Q_HOME_CAT >=1, 1, Q_HOME_CAT), 
         Q_LEARN_CAT = if_else(Q_LEARN_CAT >=1, 1, Q_LEARN_CAT))
z <- z %>% mutate(N = apply(z[,z_vars], 1, sum))
z <- z %>% mutate(n_oth_risk_gps = cut(N, breaks=c(-1,0,1,max(N)), labels=c("0","1","2+")) )
z <- z %>% dplyr::select(EAVE_LINKNO, n_oth_risk_gps)

df_all <- df_all %>% left_join(z, by="EAVE_LINKNO")
df <- df %>% left_join(z, by="EAVE_LINKNO")

#endpt is the response, pyears is the offset


#add flag for symptomatic and lh - not done
z <- cdw_full  %>% filter(test_result=="POSITIVE") %>% 
  dplyr::select(EAVE_LINKNO, date_ecoss_specimen, test_result_record_source, date_onset_of_symptoms, flag_covid_symptomatic) %>% 
  arrange(EAVE_LINKNO, date_ecoss_specimen) %>%
  filter(!duplicated(paste(EAVE_LINKNO, date_ecoss_specimen))) %>% 
  mutate(lab = if_else(test_result_record_source == "ECOSS", "nhs","lh") ) %>% 
  mutate(date_onset_of_symptoms = as.Date(date_onset_of_symptoms )) %>% 
  mutate(flag_covid_symptomatic = if_else(!is.na(flag_covid_symptomatic) & flag_covid_symptomatic=="true", 1L, 0L))

z_df <- df_all %>% left_join(dplyr::select(z, EAVE_LINKNO, date_ecoss_specimen, lab, flag_covid_symptomatic),
                             by=c("EAVE_LINKNO","date_ecoss_specimen"))
df_all <- z_df
z_df <- df %>% left_join(dplyr::select(z, EAVE_LINKNO, date_ecoss_specimen, lab, flag_covid_symptomatic),
                         by=c("EAVE_LINKNO","SpecimenDate"="date_ecoss_specimen"))
df <- z_df

saveRDS(df_all, paste0("./output/temp/df_all_full_",output_list$endpoint,".RDS"))
saveRDS(df, paste0("./output/temp/df_full_",output_list$endpoint,".RDS"))

remove(list=ls(pa="^z"))

#get the df for  hosp from the df for hosp_death
output_list$endpoint <- "covid_hosp"
#find the hosp_deaths that don't have a hosp
z <- df_all %>% filter(event==1) %>% dplyr::select(EAVE_LINKNO, event_date) %>% 
  left_join(df_endpoints, by="EAVE_LINKNO") %>% 
  filter(is.na(covid_hosp_date)) %>% 
  mutate(death_only=1) %>% 
  dplyr::select(EAVE_LINKNO, death_only)

df <- df %>% left_join(z, by="EAVE_LINKNO") %>% #filter(!is.na(death_only))
  mutate(endpt = if_else(!is.na(death_only) & death_only==1 & endpt==1, 0,endpt)) %>% 
  dplyr::select(-death_only)
saveRDS(df, paste0(project_path,"/output/temp/df_full_",output_list$endpoint,".RDS"))

