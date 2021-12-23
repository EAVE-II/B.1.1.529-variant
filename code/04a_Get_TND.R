##########################################################
# Name of file: 03a_Get_Matched_TND.R
# Data release (if applicable):
# Original author(s): Chris Robertson chris.robertson@phs.scot
# Original date: 07 Dec 2021
# Latest update author (if not using version control) - Chris Robertson chris.robertson@phs.scot
# Latest update date (if not using version control) - 
# Latest update description (if not using version control)
# Type of script: Descriptive stats
# Written/run onL: R Studio SERVER
# Version of R that the script was most recently run on: R 3.6.1
# Description of content: creates a data set for a matched TND study
#                         matching on date of test/symp
#                         run 01a_Input_Data_Omicron.R to get all the data sets
# Approximate run time: Unknown
##########################################################

output_list <- list()
output_list$endpoint <- "pos_test"
a_begin <- as.Date("2021-11-15")  #beginning of omicron - 21 Nov
output_list$a_begin <- a_begin
output_list$event_date <- "symptom"

#only use symptomatic tests - make symptom onset 5 days before specimen date if it is missing or after specimen date (not many after)
z <- cdw_full %>%  filter(date_ecoss_specimen >= a_begin) %>% 
  filter(test_result_record_source == "NHS DIGITAL") %>% #omit non community cases
  filter(flag_covid_symptomatic == "true") %>% 
  dplyr::select(EAVE_LINKNO, test_id, subject_sex, age, date_ecoss_specimen, test_result, date_onset_of_symptoms) %>% 
  dplyr::rename(age_year=age) %>% 
  mutate(date_onset_of_symptoms = as_date(date_onset_of_symptoms)) %>% 
  mutate(date_onset_of_symptoms = if_else(is.na(date_onset_of_symptoms), date_ecoss_specimen -5, date_onset_of_symptoms)) %>% 
  mutate(date_onset_of_symptoms = if_else(!is.na(date_onset_of_symptoms) & date_onset_of_symptoms > date_ecoss_specimen, date_ecoss_specimen -5, date_onset_of_symptoms)) %>% 
  filter(date_onset_of_symptoms >= date_ecoss_specimen-7) # select those with symptom onset within previous 7 days of test
z$event_date <- if (output_list$event_date == "symptom") z$date_onset_of_symptoms else z$date_ecoss_specimen
z <- z %>% filter(event_date >= a_begin)

#link in sgene status
z_sg <- sgene %>% dplyr::select(test_id, sgene_classification)  #EAVE_LINKNO is missing for about 10% of sgene - there are a few cases where the date ecoss specimen is not the same (273/59357)
z <- z %>%  left_join(z_sg, by="test_id")

#now select those who are in EAVE cohort 
z <- z %>% filter(EAVE_LINKNO %in% EAVE_cohort$EAVE_LINKNO)
z <- z %>% left_join(EAVE_cohort, by="EAVE_LINKNO")
z <- z %>%  filter(is.na(NRS.Date.Death) | NRS.Date.Death >= event_date)

#get the first positive test in the period
z_pos <- z %>% filter(test_result=="POSITIVE") %>% 
  arrange(EAVE_LINKNO, event_date) %>% 
  filter(!duplicated(EAVE_LINKNO))

#negative tests - omit anyone who is positive
#keep random one test per person
z_neg <- z %>% filter(test_result=="NEGATIVE") %>% 
  filter(!(EAVE_LINKNO %in% z_pos$EAVE_LINKNO)) 
z_neg <- z_neg %>%   mutate(random_id = runif(nrow(z_neg))) %>% 
  arrange(EAVE_LINKNO, random_id) %>% 
  filter(!duplicated(EAVE_LINKNO))
  
z1 <- bind_rows(z_pos, dplyr::select(z_neg, -random_id)) 
print(table(z1$test_result, z1$sgene_classification, exclude=NULL))

#using the time of event
#z_n_pv_pers <- 7
df <- z1 %>% left_join(Vaccinations, by="EAVE_LINKNO") %>%
  filter(is.na(flag_incon) | flag_incon==0) %>% #omit those vaccinated with inconsistent records
  #correct missing date_vacc_2 - need to clean up Vaccinations - chaked and the vacc 3 looks OK
  mutate(date_vacc_2 = if_else(!is.na(date_vacc_2) & !is.na(date_vacc_3) & date_vacc_2==date_vacc_3, date_vacc_1+77, date_vacc_2)) %>%
  filter(!(!is.na(date_vacc_2) & !is.na(date_vacc_3) & date_vacc_3 <= date_vacc_2 + 28)) %>% 
  mutate(day_1 = as.numeric(event_date - date_vacc_1),
         day_2 = as.numeric(event_date - date_vacc_2),
         day_3 = as.numeric(event_date - date_vacc_3)) 
df <- df %>% mutate(day_1 = if_else(!is.na(day_1) & day_1 <= 0, NA_real_, day_1),
                          day_2 = if_else(!is.na(day_2) & day_2 <= 0, NA_real_, day_2),
                          day_3 = if_else(!is.na(day_3) & day_3 <= 0, NA_real_, day_3))
df <- df %>%   mutate(vacc_1_gp = cut(day_1, breaks= c( 0, 27, max(day_1, na.rm=T)),labels=FALSE),
                            vacc_2_gp = cut(day_2, breaks= c(0, 13, 69, 104, 139, 174, max(day_2, na.rm=T)),labels=FALSE),
                            vacc_3_gp = cut(day_3, breaks= c( 0, 7, 13, max(day_3, na.rm=T)),labels=FALSE))
df <- df %>% mutate(vacc_1_gp = case_when(is.na(vacc_1_gp) ~ "v1_uv",
                                                vacc_1_gp==1 ~ "v1_0:3",
                                                vacc_1_gp==2 ~ "v1_4+") )
df <- df %>% mutate(vacc_2_gp = case_when(is.na(vacc_2_gp) ~ "v2_uv",
                                                vacc_2_gp==1 ~ "v2_0:1",
                                                vacc_2_gp==2 ~ "v2_2:9",
                                                vacc_2_gp==3 ~ "v2_10:14",
                                                vacc_2_gp==4 ~ "v2_15:19",
                                                vacc_2_gp==5 ~ "v2_20:24",
                                                vacc_2_gp==6 ~ "v2_25+") )
df <- df %>% mutate(vacc_3_gp = case_when(is.na(vacc_3_gp) ~ "v3_uv",
                                                vacc_3_gp==1 ~ "v3_0",
                                                vacc_3_gp==2 ~ "v3_1",
                                                vacc_3_gp==3 ~ "v3_2+") )

df <- df %>% mutate(vs=case_when(vacc_3_gp != "v3_uv" ~ vacc_3_gp,
                                       vacc_2_gp != "v2_uv" ~ vacc_2_gp,
                                       TRUE ~ vacc_1_gp)) %>% 
  mutate(vs = if_else(vs=="v1_uv", "uv",vs))
#z_labs <- sort(unique(df_cc$vs))
z_labs <- c("uv" ,"v1_0:3", "v1_4+", "v2_0:1", "v2_2:9", "v2_10:14", "v2_15:19", "v2_20:24", "v2_25+", "v3_0", "v3_1", "v3_2+" )
df <- df %>% mutate(vs=factor(vs, levels=z_labs))


print(table(df$test_result, df$sgene_classification, exclude=NULL))


#add in days from event date to latest positive test before a_begin
z <- Positive_Tests %>% filter(date_ecoss_specimen < a_begin) %>% arrange(EAVE_LINKNO, desc(date_ecoss_specimen) ) %>% 
  filter(!duplicated(EAVE_LINKNO)) %>% 
  dplyr::select(-test_id)
z1 <- df %>% left_join(z, by="EAVE_LINKNO", suffix=c("","_prior")) %>% 
  mutate(days_from_previous_pos = as.numeric(event_date - date_ecoss_specimen_prior)) %>% 
  mutate(days_from_previous_pos = if_else(!is.na(days_from_previous_pos) & days_from_previous_pos <=0, NA_real_, days_from_previous_pos))
z1 <- z1 %>% mutate(pos_before_test = case_when(is.na(days_from_previous_pos) ~ "not_prev_pos",
                                                days_from_previous_pos >= 1 & days_from_previous_pos <= 28 ~ "pos_1:28",
                                                days_from_previous_pos >= 29 & days_from_previous_pos <= 90 ~ "pos_29:90",
                                                TRUE ~ "pos_91+")) %>% 
  dplyr::select(-date_ecoss_specimen_prior)
table(z1$pos_before_test, z1$test_result,exclude=NULL)
df <- z1

#add in number of prior tests
#get number of tests before vaccination variable
z <- df %>% dplyr::select(EAVE_LINKNO, event_date)
z <- z %>% left_join(dplyr::select(cdw_full, EAVE_LINKNO, date_ecoss_specimen), by="EAVE_LINKNO") %>%
  filter(!is.na(date_ecoss_specimen)) %>% #omit those never tested
  filter(date_ecoss_specimen < event_date) #keep those sample before event date
z <- z %>% group_by(EAVE_LINKNO) %>% dplyr::summarise(N=n()) 
#now link back
df <- df %>% left_join(z, by=c("EAVE_LINKNO")) %>% 
  mutate(N=if_else(is.na(N), 0L,N))
df <- df %>% mutate(n_tests_gp = cut(N, breaks=c(-1,0,1,2,3,4,9,19, max(N)), labels=c("0","1","2","3","4","5-9","10-19","20+")))
df$N <- NULL 


z <- df %>% left_join(rg, by="EAVE_LINKNO")
df <- z

df <- df %>% mutate(shielding = if_else(EAVE_LINKNO %in% shielding$EAVE_LINKNO, 1, 0))
df <- df %>% mutate(immuno = if_else(EAVE_LINKNO %in% immuno$EAVE_LINKNO, 1, 0))


df <- df %>% 
  mutate(vacc_type = if_else(is.na(vacc_type) , "uv", as.character(vacc_type)) ) %>% 
  filter(vacc_type != "UNK") %>% 
  mutate(vacc_type = factor(vacc_type, levels=c("uv","AZ","Mo","PB") ) ) %>% 
  mutate(vt = fct_cross(vs,vacc_type, sep="_")) %>% 
  mutate(vt = fct_recode(vt, "uv" ="uv_uv", "uv" = "uv_AZ", "uv" = "uv_Mo","uv" = "uv_PB"))

#get the endpoints
df <- df%>%  mutate(event= if_else(test_result=="POSITIVE",1,0)) %>% 
  mutate(event_s_neg = if_else(event==1 & !is.na(sgene_classification) & sgene_classification=="True S Gene Dropout", 1,0),
        event_s_pos = if_else(event==1 & !is.na(sgene_classification) & sgene_classification=="Positive S Gene", 1,0))

print(table(df$event, df$sgene_classification, exclude=NULL))
print(table(df$event_s_pos, df$sgene_classification, exclude=NULL))
print(table(df$event_s_neg, df$sgene_classification, exclude=NULL))

df <- df %>% mutate(days = as.numeric(date_ecoss_specimen - min(date_ecoss_specimen)))

saveRDS(df, paste0("./output/temp/df_tnd_",output_list$endpoint,".RDS"))

z_df <- df 

z_ids <- df %>% filter(event==1 & pos_before_test %in% c("pos_1:28")) %>% pull(EAVE_LINKNO)

z_df <- filter(df, event_s_neg==1 | event==0 )  #s negs

#z_df <- filter(df, event_s_pos==1 | event==0 )  #s pos
z_levs <- levels(z_df$vs)
z_ref <- "v2_25+"
z_df <- z_df %>% filter(!(EAVE_LINKNO %in% z_ids)) #%>% filter(age_year >= 16 & age_year < 50) 
z_df <- z_df %>% mutate(vs = fct_relevel(vs, z_ref))
table(z_df$event, z_df$sgene_classification, exclude=NULL)
table(z_df$vs, z_df$event, exclude=NULL)


library(mgcv)
z <- gam(event ~ s(days) + s(age_year) + subject_sex + simd2020_sc_quintile + n_risk_gps + vs + pos_before_test + HB + n_tests_gp + shielding + immuno,
         data=z_df, family="binomial")

z_est <- fun_ve_glm(z)
z1 <- z_est %>% filter(grepl("vs", names)) %>% dplyr::rename(vs=names) %>% 
     mutate(vs = gsub("^vs","", vs)) %>% 
     mutate(across(c("est","lcl","ucl"), ~exp(.)))
z1 <- bind_rows(data.frame(est=1, lcl=1,ucl=1, vs=z_ref), z1)

z2 <- z_df %>% group_by(vs) %>% dplyr::summarise(N=n(), R=sum(event)) %>% as.data.frame()
z3 <- z2 %>% left_join(z1, by="vs") %>% mutate(vs=factor(vs, levels=z_levs)) %>% arrange(vs) #%>% 
   #mutate(ve=100*(1-est), x=100*(1-ucl), ucl=100*(1-lcl), lcl=x)
 #z3 <- z3 %>% dplyr::select( vs, N, R, ve, lcl, ucl)
z3



###########################################################################
#multinomial
library(nnet)

z_ids <- df %>% filter(event==1 & pos_before_test %in% c("pos_1:28")) %>% pull(EAVE_LINKNO)
z_df <- df %>% filter(!(EAVE_LINKNO %in% z_ids)) %>% filter(age_year >= 16) 
#z_df <- z_df %>% filter(vs != "uv") %>%  mutate(vs = fct_relevel(vs, "v2_2:9"))
z_df <- z_df %>% 
  mutate(response=case_when(event==0 ~ "neg",
                            event_s_pos==1 ~ "s_pos",
                            event_s_neg==1 ~ "s_neg",
                            TRUE ~ "other")) %>% 
  filter(response!="other") %>% 
  mutate(response=factor(response, levels=c("neg","s_pos","s_neg"))) %>% 
  mutate(response=as.numeric(response) -1) %>% 
  mutate(age_gp=cut(ageYear, breaks = c(-1,11,17,29,39,49,59,69,120), labels=c("0-11","12-17","18-29","30-39","40-49","50-59","60-69","70+"))) %>% 
  mutate(days_gp = cut(days, breaks=c(-1,6,13,20,27)))

table(z_df$response, z_df$sgene_classification, exclude=NULL)

z <-  gam(list(response ~ s(days) + s(age_year) + subject_sex  + simd2020_sc_quintile + n_risk_gps + pos_before_test + HB + n_tests_gp + shielding + immuno + vs,
               ~ s(days) + s(age_year) + subject_sex  + simd2020_sc_quintile + n_risk_gps + pos_before_test + HB + n_tests_gp + shielding + immuno + vs),
              data=z_df, family=multinom(K=2))

