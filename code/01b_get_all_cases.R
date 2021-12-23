##########################################################
# Name of file: 07_01d_s_gene_positive_all_cases.R
# Data release (if applicable):
# Original author(s): Chris Robertson chrisobertson@nhs.net
# Original date: 06 August 2020
# Latest update author (if not using version control) - Chris Robertson chrisobertson@nhs.net
# Latest update date (if not using version control) - 
# Latest update description (if not using version control)
# Type of script: Descriptive stats
# Written/run onL: R Studio SERVER
# Version of R that the script was most recently run on: R 3.6.1
# Description of content: sets up analysis for all positive tests
#                         run 07_01a_s_gene_positive_description.R
#                         to read in the data - line 291
#                         then get first positive test post 01 April
# Approximate run time: Unknown
##########################################################

#Need to get Positive Tests going right back
a_begin <- as_date("2021-11-01")  #start of omicron
z_df <- cdw_full %>% filter(test_result == "POSITIVE") %>% filter(date_ecoss_specimen >= a_begin)
#get the first positive test in the period
z_df <- z_df %>% arrange(EAVE_LINKNO, date_ecoss_specimen) %>% 
  filter(!duplicated(EAVE_LINKNO))
#Ecoss is the nhs labs
z_df <- z_df %>% mutate(lab = if_else(test_result_record_source == "ECOSS", "nhs","lh") )
#  dplyr::select(EAVE_LINKNO, ecossid, sex, age_year, specimen_date, lab)
z_df <- z_df %>%  
  dplyr::select(EAVE_LINKNO, test_id, subject_sex, age, date_ecoss_specimen, lab, test_result_record_source, flag_covid_symptomatic) %>% 
  dplyr::rename(age_year=age)

#get the last positive test for anyone before a_begin
z <- Positive_Tests %>% filter(date_ecoss_specimen < a_begin) %>% 
  dplyr::select(EAVE_LINKNO, date_ecoss_specimen) %>% 
  arrange(EAVE_LINKNO, desc(date_ecoss_specimen)) %>% 
  filter(!duplicated(EAVE_LINKNO)) 

z_df <- z_df %>% left_join(z, by="EAVE_LINKNO", suffix=c("","_first"))
z_df <- z_df %>% mutate(days_from_previous_pos = as.numeric(date_ecoss_specimen - date_ecoss_specimen_first)) %>% 
  dplyr::select(-date_ecoss_specimen_first)

#for those with a test_id - lighthouse - match on this 
z <- sgene %>% dplyr::select(test_id, sgene_classification) %>% filter(!is.na(test_id))
z_lh <- z_df %>% left_join(z, by="test_id")
z_oth <- filter(z_lh, is.na(sgene_classification))  # those who did not get a match using test_id
z_lh <- filter(z_lh, !is.na(sgene_classification))  #those who did get a match using test_id

#for those who do not link on test_id link on eave_linkno and date ecoss specimen
#first of all deduplicate to get the first positive test in the study period
z <- sgene %>% dplyr::select(EAVE_LINKNO, date_ecoss_specimen, sgene_classification) %>% 
  arrange(EAVE_LINKNO, date_ecoss_specimen) %>% 
  filter(!duplicated(EAVE_LINKNO))
z_oth <- z_oth %>% dplyr::select(-sgene_classification) %>% left_join(z, by=c("EAVE_LINKNO","date_ecoss_specimen"))

print(nrow(z_lh) + nrow(z_oth) == nrow(z_df))
table(z_lh$EAVE_LINKNO %in% z_oth$EAVE_LINKNO)
z_df <- bind_rows(z_lh,z_oth)

z_df <- z_df %>%  mutate(sgene_classification = if_else(is.na(sgene_classification), "unknown", sgene_classification))
z_df <- z_df %>% mutate(s_gene = factor(sgene_classification, 
                  levels=c("Positive S Gene","True S Gene Dropout","Weak Positive","other", "unknown"),
                  labels=c("S_Pos","S_Neg","Weak_S_Pos", "Other", "Unknown")))
z_df$sgene_classification <- NULL

#link in hospitalisations - use all as some will be in hospital at time of test
#keep all admissions post a_begin and those in hospital at a_begin
z_h <- all_hospitalisations %>% dplyr::select(-validchi) %>% 
  filter(!(!is.na(discharge_date) & discharge_date <= a_begin))
#summary(filter(z_h, admission_date < a_begin))  # all admitted before are discharged after a_begin or no discharge date

# Steven: I have taken the end date of hospitalisation records before we remove any records
# for people who have had multiple hospitalisations.
a_end_hosp <- max(c(max(z_h$admission_date, na.rm=T), max(z_h$discharge_date, na.rm=T)))


z <- z_df %>%
  left_join(z_h, by="EAVE_LINKNO") %>%
  # Create variable for if they were in hospital at time of test
  # I think Chris's original code can miss if they were in hospital at time of test for people
  # with multiple admissions, since it ony keeps the first admission post test.
  mutate(In_Hosp_At_Test = case_when(date_ecoss_specimen > a_end_hosp ~ 'unknown',
                                     date_ecoss_specimen >= admission_date & date_ecoss_specimen <= discharge_date ~ 'yes',
                                     TRUE ~ 'no')) %>%
  group_by(EAVE_LINKNO) %>%
  mutate(In_Hosp_At_Test = case_when( any(In_Hosp_At_Test == 'yes') ~ 'yes',
                                      TRUE ~ In_Hosp_At_Test)) %>%
  ungroup()


z <- z %>% mutate(ad=admission_date, dd=discharge_date, em = emergency) %>% 
  mutate(ad = as.Date(ifelse(!is.na(dd) & dd < date_ecoss_specimen, NA, ad), origin="1970-01-01"),
         em = ifelse(!is.na(dd) & dd < date_ecoss_specimen, NA, em), 
         dd = as.Date(ifelse(!is.na(dd) & dd < date_ecoss_specimen, NA, dd), origin="1970-01-01") )
#multiple admissions - take the first admission in the study period for those with multiple
#this will be the first after specimen date or if before will have no discharge
z_id <- z %>% filter(duplicated(EAVE_LINKNO)) %>% pull(EAVE_LINKNO) %>% unique()
z_01 <- z %>% filter(!(EAVE_LINKNO %in% z_id)) # 0 or 1 admission
z_m <- z %>% filter((EAVE_LINKNO %in% z_id))  #multiple admissions
z_m <- z_m %>%  arrange(EAVE_LINKNO, ad) %>% 
  filter(!duplicated(EAVE_LINKNO))  # NA on ad go to the end, pick the first admission post specimen

z <- bind_rows(z_01, z_m) %>% 
  dplyr::select(-admission_date, -discharge_date, -emergency) %>% 
  dplyr::rename(admission_date=ad, discharge_date=dd, emergency=em)

z_df <- z

#now modify the admission dates for those admitted a long time ago but with no discharge and
#evidence of a recent test in the community 
#impute discharge dates for people admitted a long time ago with no discharge
#so different rules for those who are tested in teh community as opposed to thos tests in the nhs labs (ECOSS)
#z <- z_df %>% filter(!is.na(admission_date) & admission_date < date_ecoss_specimen - 7 & is.na(discharge_date) & test_result_record_source != "ECOSS")
#if there is an admission date more than 7 days before the specimen date and the test was a lighhouse
#or more than 30 days before and the test was in NHS
#- asssume the person was discharged and set admission date to NA so it won't count as a covid admission (in hosp at time of test)
z_df <- z_df %>% mutate(discharge_date_orig = discharge_date, admission_date_orig = admission_date)
z_df <- z_df %>%
  mutate(admission_date = case_when(
    !is.na(admission_date_orig) & admission_date_orig < date_ecoss_specimen - 7 & is.na(discharge_date_orig) & test_result_record_source != "ECOSS" ~ NA_Date_,
    !is.na(admission_date_orig) & admission_date_orig < date_ecoss_specimen - 30 & is.na(discharge_date_orig) & test_result_record_source == "ECOSS" ~ NA_Date_,
    TRUE ~ admission_date_orig))
#z_df <- z_df %>% 
#  mutate(discharge_date = case_when(
#    !is.na(admission_date_orig) & admission_date_orig < a_begin - 7 & is.na(discharge_date_orig) & test_result_record_source != "ECOSS" ~ admission_date_orig+6,
#    !is.na(admission_date_orig) & admission_date_orig < a_begin - 30 & is.na(discharge_date_orig) & test_result_record_source == "ECOSS" ~ admission_date_orig+14,
#    TRUE ~ discharge_date_orig))
#change the admission and discharge dates to NA if the imputed discharge date is before the speciment date
#use discharge_date as that has the imputed discharges - use admission_data_orig as that is unchanged
#z_df <- z_df %>%  mutate(admission_date = if_else(
#    !is.na(admission_date_orig) & !is.na(discharge_date) & discharge_date < date_ecoss_specimen,  NA_Date_, admission_date_orig)) %>% 
#  mutate(discharge_date = if_else(
#     !is.na(discharge_date) & discharge_date < date_ecoss_specimen,  NA_Date_, discharge_date))

#a_end_hosp <- max(c(max(z_df$admission_date, na.rm=T), max(z_df$discharge_date, na.rm=T)))

z_df <- z_df %>% mutate(Time.To.Hosp = case_when(
  !is.na(admission_date) & !is.na(discharge_date) & discharge_date <= date_ecoss_specimen ~ as.numeric(a_end_hosp - date_ecoss_specimen),
  !is.na(admission_date) & is.na(discharge_date)  ~  as.numeric(admission_date - date_ecoss_specimen),
  !is.na(admission_date) & !is.na(discharge_date) & discharge_date > date_ecoss_specimen ~ as.numeric(admission_date - date_ecoss_specimen),
  is.na(admission_date) ~ as.numeric(a_end_hosp - date_ecoss_specimen),
  TRUE ~ NA_real_) ) %>% 
 # mutate(In_Hosp_At_Test = ifelse(Time.To.Hosp <= -3  ~ "yes", In_Hosp_At_Test) %>% 
  mutate(hosp_covid = if_else(!is.na(admission_date) & Time.To.Hosp <= 14,  1L, 0L ),
         Time.To.Hosp = case_when(Time.To.Hosp < 0 ~ 0,
                                  Time.To.Hosp >= 15 ~ 15,
                                  TRUE ~ Time.To.Hosp))
z_df <- z_df %>% mutate(hosp_covid_emerg = if_else(!is.na(emergency) & emergency ,hosp_covid, 0L ) )
z_df <- z_df %>% mutate(days = as.numeric(date_ecoss_specimen- min(date_ecoss_specimen)) )

#add in the vaccinations and risk groups
z_df <- z_df %>%  
  left_join(Vaccinations, by="EAVE_LINKNO" )
z_df <- z_df %>%  mutate(vs1 = case_when(is.na(date_vacc_1) | date_vacc_1 > date_ecoss_specimen ~ "uv",
                                         date_vacc_1 <= date_ecoss_specimen &  date_vacc_1 > date_ecoss_specimen - 28 ~ "v1_0:3",
                                         TRUE ~ "v1_4+"),
                         vs2 = case_when(is.na(date_vacc_2) | date_vacc_2 > date_ecoss_specimen ~ "uv",
                                         date_vacc_2 <= date_ecoss_specimen &  date_vacc_2 > date_ecoss_specimen - 14 ~ "v2_0:1",
                                         date_vacc_2 <= date_ecoss_specimen-14 &  date_vacc_2 > date_ecoss_specimen - 42 ~ "v2_2-5",
                                         date_vacc_2 <= date_ecoss_specimen-42 &  date_vacc_2 > date_ecoss_specimen - 69 ~ "v2_6-9",
                                         TRUE ~ "v2_10+"),
                         vs3 = case_when(is.na(date_vacc_3) | date_vacc_3 > date_ecoss_specimen ~ "uv",
                                         date_vacc_3 <= date_ecoss_specimen &  date_vacc_3 > date_ecoss_specimen - 14 ~ "v3_0:1",
                                         TRUE ~ "v3_2+") ) %>% 
  mutate(vs = case_when(vs3 !="uv" ~ vs3,
                        vs2 !="uv" ~ vs2,
                        TRUE ~ vs1))
z_df <- z_df %>% mutate(vs=factor(vs, levels=c("uv","v1_0:3","v1_4+","v2_0:1","v2_2-5","v2_6-9","v2_10+","v3_0:1","v3_2+")))


z_df <- z_df %>% 
  mutate(vacc_type = if_else(is.na(vacc_type) , "uv", as.character(vacc_type)) ) %>% 
  filter(vacc_type != "UNK") %>% 
  mutate(vacc_type = factor(vacc_type, levels=c("uv","AZ","Mo","PB") ) ) %>% 
  mutate(vt = fct_cross(vs,vacc_type, sep="_")) %>% 
  mutate(vt = fct_recode(vt, "uv" ="uv_uv", "uv" = "uv_AZ", "uv" = "uv_Mo","uv" = "uv_PB"))

#add in demographics
z_df <- z_df %>%  
  left_join(dplyr::select(EAVE_cohort, -Sex, -ageYear, -NRS.Date.Death), by="EAVE_LINKNO" )
z_df <- z_df %>% mutate(in_eave = if_else(is.na(eave_weight), 0L,1L))


z_df <- z_df %>% mutate(age_gp = cut(age_year, breaks = c(-1, 11, 19, 39, 59, 74, 120),
                                     labels=c("0-11", "12-19","20-39","40-59","60-74","75+")))
z_df <- z_df %>% dplyr::rename(sex=subject_sex)
z_df <- z_df %>% filter(!is.na(age_year))

#add in risk groups as at Dec 2020
z_df <- z_df %>% left_join(rg, by="EAVE_LINKNO")
z_df <- z_df %>% mutate(n_risk_gps = fct_explicit_na(n_risk_gps, na_level = "Unknown"))

#add in icu/death  - icu dates depend upon the endpoints linkage

#covid death derived from all deaths
a_end_death <- max(all_deaths$NRS.Date.Death)
z <- z_df %>% left_join(dplyr::select(all_deaths, EAVE_LINKNO, NRS.Date.Death, covid_death_cert ), by="EAVE_LINKNO") %>% 
  mutate(death = if_else(is.na(NRS.Date.Death), 0L, 1L)) %>% 
  mutate(Time.To.Death = if_else(is.na(NRS.Date.Death), as.numeric(a_end_death - date_ecoss_specimen),
                                       as.numeric(NRS.Date.Death - date_ecoss_specimen))) %>% 
  mutate(Time.To.Death = if_else(Time.To.Death < 0, 0, Time.To.Death)) %>% 
  mutate(covid_death = if_else(is.na(covid_death_cert), 0, covid_death_cert)) %>% 
  mutate(covid_death = if_else(!is.na(NRS.Date.Death) & (NRS.Date.Death - date_ecoss_specimen <= 28) , 1, covid_death))
z_df <- z

z_df <- z_df %>% mutate(prev_pos = case_when(is.na(days_from_previous_pos) ~ "not_prev_pos",
                                             days_from_previous_pos >= 1 & days_from_previous_pos <= 28 ~ "pos_1:28",
                                             days_from_previous_pos >= 29 & days_from_previous_pos <= 90 ~ "pos_29:90",
                                             TRUE ~ "pos_91+"))

z_df <- z_df %>% mutate(sgx = if_else(s_gene == "Other", "Unknown", as.character(s_gene))) %>% 
  mutate(sg_lab = case_when(sgx != "Unknown" ~ as.character(sgx),
                                        TRUE ~ paste(sgx, lab, sep="_"))) %>% 
  mutate(sg_lab=factor(sg_lab, levels =c("S_Neg","S_Pos","Weak_S_Pos","Unknown_lh","Unknown_nhs")))
z_df$sgx <- NULL

z.df <- z_df
z.df <- mutate(z.df, Time.To.Hosp = if_else(Time.To.Hosp==0,0.1,Time.To.Hosp))  #add 0.1 to 0 days so survival models work
z.df <- mutate(z.df, Time.To.Death = if_else(Time.To.Death==0,0.1,Time.To.Death))

z.df <- z.df %>% left_join(dplyr::select(wgs, EAVE_LINKNO, VariantofInterest), by="EAVE_LINKNO")   #collection date is date ecoss specimen
z.df <- z.df %>% mutate(variant = case_when(is.na(VariantofInterest) ~ "not_sequenced",
                                            VariantofInterest=="VOC-20DEC-01" ~ "alpha",
                                            VariantofInterest=="VOC-21APR-02" ~ "delta",
                                            VariantofInterest=="VOC-21NOV-01" ~ "omicron",
                                            TRUE ~ "other")) %>% 
  dplyr::select(-VariantofInterest)


z.df <- z.df %>% mutate(bmi.gp = cut(bmi_impute, breaks=c(-1, 20,25,30,35,40,51), labels=FALSE))
z.df <- z.df %>% mutate(bmi_gp = case_when(age_year <= 17 ~ "NA_too_young",
                                           age_year >= 18 & !is.na(bmi.gp) ~ as.character(bmi.gp),
                                           TRUE ~ "Unknown")) %>% 
  mutate(bmi_gp = factor(bmi_gp, levels=c("1","2","3","4","5","6","NA_too_young","Unknown"),
                         labels=c("<20","20-24","25-29","30-34","35-39","40+","NA_too_young","Unknown")))

z_df <- z.df
z_df$eave_weight[is.na(z_df$eave_weight)] <- 0
z_df$wght_all <- 1

z_df <- z_df %>% mutate(shielding = if_else(EAVE_LINKNO %in% shielding$EAVE_LINKNO, 1, 0))
z_df <- z_df %>% mutate(immuno = if_else(EAVE_LINKNO %in% immuno$EAVE_LINKNO, 1, 0))


rgs <- colnames(z_df)[startsWith(colnames(z_df), "Q")]
z_df <- z_df %>% mutate(across(all_of(rgs), ~ as.factor(.)))

saveRDS(z_df, "output/temp/Sgene_all_positive.rds") 


