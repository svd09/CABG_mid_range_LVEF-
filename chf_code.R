##################################################################
##                         CHF analysis                         ##
##################################################################



library(easypackages)
libraries(c('tidyverse','rms','naniar','Hmisc',
            'MASS', 'tableone','haven',"lubridate", "survival", "mstate",
            "naniar", "survminer", "ggsci", "reda","reReg", "naniar", "splines"))

# get my main dataset

dt = 
  read_csv('P:/ORD_Perez_201602128D/Deo/CABG_HF/JTCVS_paper/data/dt.csv')

# get the readmits information from SAS 

re = read_sas("P:\\ORD_Perez_201602128D\\Deo\\CABG_HF\\data\\deo_hf_readmits.sas7bdat")

# now am going to first limit this to the patients that are part of my group
# to do that I need to use the crosswalk to join the scrssn and patientsid 

cw = read_sas("P:\\ORD_Perez_201602128D\\Deo\\CABG_HF\\data\\deo_cwalk.sas7bdat")

# now from the cw , limit it to the patients that are in my group

mypat <- dt$scrssn

cw$mine <- with(cw, ifelse(ScrSSN %in% mypat, 1, 0))

cw %>% count(mine)

# keep those patients that are in my main group

cw2 = cw %>% filter(mine == 1)

# now this is the list of patientsid are part of my data.

mypat2 <- cw2$PatientSID 

re$mine <- with(re, ifelse(PatientSID %in% mypat2, 1, 0))

re %>% count(mine)

# now to limit the readmission data to those patients that are in my group

re2 = re %>% filter(mine == 1)

# see unique scrssn in this 

re3 = left_join(re2, cw2, by = "PatientSID")

# keep only the required col

re4 = re3[, c("PatientSID", "AdmitDatetime", "PrincipalDiagnosisICD9SID", 
              "PrincipalDiagnosisICD10SID",  "ScrSSN", "PatientICN")]

# now that I have scrssn, see the actual # of patients that had readmissions.

length(unique(re4$ScrSSN))



df = dt %>% dplyr::select(scrssn, ACT_LAST_DT, surgdate)

df$ACT_LAST_DT = as_date(df$ACT_LAST_DT)

df$surgdate = as_date(df$surgdate)

names(re4) = tolower(names(re4))


re5 = left_join(re4, df, by = "scrssn")


re5$admitdatetime = as_date(re5$admitdatetime)


re5$within = with(re5, ifelse((admitdatetime > surgdate & admitdatetime <= ACT_LAST_DT), 1, 0
))

re5 %>% count(within)

re6 = re5 %>% filter(within == 1)

# now re6 has those readmissions for the patients within my cohort and also within the dates.
# now we need to flag those that were for CHF readmissions.

# now to get the codes for HF admissions 

icd9_hf = read_sas("P:\\ORD_Perez_201602128D\\Deo\\CABG_HF\\data\\deo_icd9_hf.sas7bdat")

icd10_hf = read_sas("P:\\ORD_Perez_201602128D\\Deo\\CABG_HF\\data\\deo_icd10_hf.sas7bdat")


hf_icd9 = icd9_hf$icd9sid

re6$hf_icd9 = with(re6, ifelse((principaldiagnosisicd9sid %in% hf_icd9), 1, 0))

re6 %>% count(hf_icd9)

hf_icd10 = icd10_hf$icd10sid

re6$hf_icd10 = with(re6, ifelse((principaldiagnosisicd10sid %in% hf_icd10), 1, 0))

re6 %>% count(hf_icd10)

re6$hf_readmit = with(re6, ifelse((hf_icd9 == 1|hf_icd10 == 1), 1, 0))

re6 %>% count(hf_readmit)

glimpse(re6)

# now am going to only limit this dataset to those with hf 
# can then convert it to wide format and then tmerge it with the base data
# have already limited it to the dates between surgery and last date of followup 

# keep on hf readmitted patients 

hf_r = re6 %>% filter(hf_readmit == 1) # THIS IS WHERE TO CONTINUE FOR AG MODEL/TIME GAP MODEL.

length(unique(hf_r$scrssn))

# 673 patients from 6533 readmitted for HF during the follow-up
# obtain counts for each patients 

total = hf_r %>% group_by(scrssn) %>% 
  summarise(total = n())

describe(total$total)

# now to get the total back into the hf_r data 

hf_r2 = left_join(hf_r, total, by = "scrssn")

# now to get the hf_r2 data in a format to join to the main data. 
# need to obtain 1st CHF event rate as a CIF model.

hf_r3 = hf_r2 %>% dplyr::select(scrssn, admitdatetime, hf_readmit)

hf_r3 = hf_r3 %>% arrange(scrssn, admitdatetime)

# keep only the first event 

hf_r4 = hf_r3[!duplicated(hf_r3$scrssn), ] # keeping only the first row

# now to get that into the main dt dataset.

dt2 = left_join(dt, hf_r4, by = "scrssn")

dt3 = dt2 %>% dplyr::select(scrssn, admitdatetime, hf_readmit, 
                            ACT_LAST_DT, surgdate, died, cohort)

# now dt3 only contains information for CIF calculation 

glimpse(dt3)

dt3$surgdate = as_date(dt3$surgdate)

dt3$ACT_LAST_DT = as_date(dt3$ACT_LAST_DT)

dt3$admitdatetime = as_date(dt3$admitdatetime)

# calculate survival time 

dt3$survdays = (dt3$surgdate %--% dt3$ACT_LAST_DT)/ddays(1)

dt3$hftime = (dt3$surgdate %--% dt3$admitdatetime)/ddays(1)

summary(dt3$survdays)

summary(dt3$hftime)

# create time column with time from both 
# convert NA in hf_readmit == 0

dt3$hf_readmit[is.na(dt3$hf_readmit)]<- 0

table(dt3$hf_readmit)

# first the time 

dt3$timevar = with(dt3, ifelse(hf_readmit == 1, hftime, survdays))

summary(dt3$timevar)

dt3$time_y = (dt3$timevar + 1)/365.25

summary(dt3$time_y)

# now for event 

dt3$event = with(dt3, ifelse(hf_readmit == 1, 1, 
                             ifelse(died == 1, 2, 0)))

dt3 %>% count(event)

# now to obtain CIF for the whole group and then by cohort


dt3$cohort_n[dt3$cohort == "NORMAL"]<- 0
dt3$cohort_n[dt3$cohort == "MID"]<- 1
dt3$cohort_n[dt3$cohort == "LOW"]<- 2

dt3$cohort_n <- factor(dt3$cohort_n)

summary(dt3$cohort_n)

chf_g = survfit(Surv(time_y, event, type = "mstate") ~ 1,
                data = dt3)

plot(chf_g)

summary(chf_g, times = c(1,5,10))

# by group 


chf_cohort = survfit(Surv(time_y, event, type = "mstate") ~ cohort_n,
                data = dt3)

plot(chf_cohort, col = c("green","blue","red"))

summary(chf_cohort, times = c(1,5,10))


# would need to create a pub graph for the paper.


dt3$event = factor(dt3$event)

chf_p = npsurv(Surv(time_y, event) ~ cohort_n, 
              data = dt3, conf.int = 0.68)


mycols = pal_jama("default")(4)

mycols2 = pal_jama("default", alpha = 0.2)(4)


tiff(filename = 
       "P:/ORD_Perez_201602128D/Deo/CABG_HF/JTCVS_paper/fig_chf.tiff",
     height = 5, width = 7, units = "in", res = 1200)



survplot(chf_p, state = "1",
         xlim = c(0,10),
         ylim = c(0,0.4),
         col = mycols,
         col.fill = mycols2,
         lwd = 2,
         lty = c(1,1,1),
         ylab = "Cumulative Estimate for CHF (First Event)",
         xlab = "Time:Years", 
         n.risk = T, label.curves = F,
         adj.n.risk = 0.5,
         time.inc = 5)

dev.off()


# AM GOING TO MODEL CHF AS RECURRENT EVENT ANALYSIS # 

# given the high rate of CHF in LOw patients, will do  segmented Cox model 
# am going to create the dataset using tmerge

df.event = dt3 %>% dplyr::select(scrssn, hftime, hf_readmit)

df.event = df.event[df.event$hf_readmit == 1, ]

dt$survdays = dt$survdays + 1

df.chf = tmerge(data1 = dt, data2 = dt, id = scrssn, tstop = survdays)

df.chf2 = tmerge(data1 = df.chf, data2 = df.event, id = scrssn, chf = event(hftime))

df.chf3 = df.chf2 %>% dplyr::select(scrssn, tstart, tstop , chf, cohort, died)


df.chf4 = df.chf3 %>% group_by(scrssn) %>% 
  mutate(rowid = row_number())

df.chf4 = data.frame(df.chf4)

df.chf4$rowid = factor(df.chf4$rowid)

df.chf4$cohort_n[df.chf4$cohort == "NORMAL"]<- 0
df.chf4$cohort_n[df.chf4$cohort == "MID"]<- 1
df.chf4$cohort_n[df.chf4$cohort == "LOW"]<- 2

df.chf4$cohort_n <- factor(df.chf4$cohort_n)

summary(df.chf4$cohort_n)


# obtain HR before and after first CHF event 


chf  = coxph(Surv(tstart, tstop, died)  ~ cohort_n, data = df.chf4[df.chf4$rowid == 1, ],
          x = TRUE, method = "breslow")
summary(chf)



chf2 = coxph(Surv(tstart, tstop, died)  ~ cohort_n, data = df.chf4[df.chf4$rowid == 2, ],
           x= TRUE, method = "breslow")

summary(chf2)

summary(b2)


#############################################################

# getting the data with tmerge to do the Andersen Gill model.

#############################################################

# OBTAIN DATA FROM HF_R TO CONTINUE FOR TIME GAP MODEL.

library(reda);library(reReg)

glimpse(hf_r)

# found that if patients transferred then they have very close admit dates; so to ensure 
# that we are not capturing that, am going to limit to having at least 3 days gap between admissions.

# first limit the col to those that we need.

dr = hf_r %>% dplyr::select(scrssn, admitdatetime, surgdate, hf_readmit)

# now need to convert this data into wide format.

dr$readmit_time = (dr$surgdate %--% dr$admitdatetime)/ddays(1)

glimpse(dr)

dr2 = dr %>% dplyr::select(scrssn, readmit_time, hf_readmit)

# first create an row_id for each patient 

dr2 = dr2 %>% arrange(scrssn, readmit_time)

dr3 = dr2 %>% group_by(scrssn) %>%
  mutate(rowid = paste0("readmit", row_number(), sep = ""))

dr_w = dr3 %>% 
  pivot_wider(id = scrssn, values_from = readmit_time, names_from = rowid)
# 
# # now the dr_w contains these columns --- data is in the wide format and 
# contains 1 row per patients, it contains readmit# col and that contains the time
# for readmission from the surgery date for each patient.
# if the time between values is < 3, then am going to consider that as a transfer and hence convert to NA
# 


# start with the first column and then go from there.

dr_w = as_tibble(dr_w)

dr_w = data.frame(dr_w)

dr_w$re2.n = with(dr_w, ifelse((readmit2 - readmit1) < 3, NA, readmit2))

dr_w$re3.n = with(dr_w, ifelse((readmit3 - readmit2) < 3, NA, readmit3))

dr_w$re4.n = with(dr_w, ifelse((readmit4 - readmit3) < 3, NA, readmit4))

dr_w$re5.n = with(dr_w, ifelse((readmit5 - readmit4) < 3, NA, readmit5))

dr_w$re6.n = with(dr_w, ifelse((readmit6 - readmit5) < 3, NA, readmit6))

dr_w$re7.n = with(dr_w, ifelse((readmit7 - readmit6) < 3, NA, readmit7))

dr_w$re8.n = with(dr_w, ifelse((readmit8 - readmit7) < 3, NA, readmit8))


dr_w$re9.n = with(dr_w, ifelse((readmit10 - readmit9) < 3, NA, readmit9))


dr_w$re10.n = with(dr_w, ifelse((readmit11 - readmit10) < 3, NA, readmit10))


dr_w$re11.n = with(dr_w, ifelse((readmit12 - readmit11) < 3, NA, readmit11))


dr_w$re12.n = with(dr_w, ifelse((readmit13 - readmit12) < 3, NA, readmit12))


dr_w$re13.n = with(dr_w, ifelse((readmit14 - readmit13) < 3, NA, readmit13))


dr_w$re14.n = with(dr_w, ifelse((readmit15 - readmit14) < 3, NA, readmit14))

dr_w$re15.n = dr_w$readmit15

dr_w$re1.n = dr_w$readmit1

dr_w.n = dr_w %>% dplyr::select(scrssn, contains (".n"))

# now the dr_w.n can be tmerged with the main dataset to do the time gap analysis
# to get the recurrent cumulative mean function,we need to format the data first.
# get the main dataaset and then start tmerge with this wide dataset.

dt = 
  read_csv('P:/ORD_Perez_201602128D/Deo/CABG_HF/JTCVS_paper/data/dt.csv')

# dt now contains 6533 patients, am going to now tmerge the data...

dt$survdays2 = dt$survdays + 1

dt1 = tmerge(data1 = dt, data2 = dt, id = scrssn, tstop = survdays2)

dt2 = tmerge(data1 = dt1, data2 = dr_w.n, id = scrssn, chf = event(re1.n))

dt3 = tmerge(data1 = dt2, data2 = dr_w.n, id = scrssn, chf = event(re2.n))

dt4 = tmerge(data1 = dt3, data2 = dr_w.n, id = scrssn, chf = event(re3.n))

dt5 = tmerge(data1 = dt4, data2 = dr_w.n, id = scrssn, chf = event(re4.n))

dt6 = tmerge(data1 = dt5, data2 = dr_w.n, id = scrssn, chf = event(re5.n))

dt7 = tmerge(data1 = dt6, data2 = dr_w.n, id = scrssn, chf = event(re6.n))

dt8 = tmerge(data1 = dt7, data2 = dr_w.n, id = scrssn, chf = event(re7.n))

dt9 = tmerge(data1 = dt8, data2 = dr_w.n, id = scrssn, chf = event(re8.n))

dt10 = tmerge(data1 = dt9, data2 = dr_w.n, id = scrssn, chf = event(re9.n))

dt11 = tmerge(data1 = dt10, data2 = dr_w.n, id = scrssn, chf = event(re10.n))

dt12 = tmerge(data1 = dt11, data2 = dr_w.n, id = scrssn, chf = event(re11.n))

dt13 = tmerge(data1 = dt12, data2 = dr_w.n, id = scrssn, chf = event(re12.n))

dt14 = tmerge(data1 = dt13, data2 = dr_w.n, id = scrssn, chf = event(re13.n))

dt15 = tmerge(data1 = dt14, data2 = dr_w.n, id = scrssn, chf = event(re14.n))

dt16 = tmerge(data1 = dt15, data2 = dr_w.n, id = scrssn, chf = event(re15.n))

# dt16 now contains the data in long format
# am going to save dt16 and then get the other variables for calculating the MCF and other analyses.

write_csv(dt16, 
          'P:/ORD_Perez_201602128D/Deo/CABG_HF/JTCVS_paper/data/dt_chf_core.csv')

glimpse(dt16)

# so as I tmerged with the main df dt, I have all the variables needed already.
# will have to change some variables and will limit the dataset for MCF analysis to only the variables needed.

# get the dataset again and then do MCF

dt16 = read_csv('P:/ORD_Perez_201602128D/Deo/CABG_HF/JTCVS_paper/data/dt_chf_core.csv')


dt_mcf = dt16 %>% dplyr::select(id, tstart, tstop, chf, survdays, died, cohort)


dt_mcf$cohort_n[dt_mcf$cohort == "NORMAL"]<- 0
dt_mcf$cohort_n[dt_mcf$cohort == "MID"]<- 1
dt_mcf$cohort_n[dt_mcf$cohort == "LOW"]<- 2

dt_mcf$cohort_n <- factor(dt_mcf$cohort_n)

#- see if we can convert the results to years 

dt_mcf$tstart.y = dt_mcf$tstart/365.25
dt_mcf$tstop.y = dt_mcf$tstop/365.25

attach(dt_mcf)

g <- Recur(tstart.y %to% tstop.y, id = id, event = chf)

plot(g2)


# obtain overall rates for CHF readmissions

mcf_overall = mcf(g ~ 1, data = dt_mcf,
                  variance = "bootstrap", level = 0.68,
                  )

str(mcf_overall)

# am going to extract the information from MCF and then create a plot

res = mcf_overall@MCF

# res is actually a dataframe, so am going to save this for further use

write_csv(res, 
          'P:/ORD_Perez_201602128D/Deo/CABG_HF/JTCVS_paper/overall_mcf.csv')

# now this can be used to plot the graphs & also get estimates.

# can use this to plot as base R or ggplot2
# excellent plot - am going to provide the estimates from the data and then
# paste them into the plot later.


tiff('P:/ORD_Perez_201602128D/Deo/CABG_HF/JTCVS_paper/mcf_overall.tiff',
     height = 7, width = 5, units = "in", res = 1200)

plot(x = res$time, y = 100*res$MCF, type = "s",
     xlab = "Followtime:Years",
     ylab = "Event Rate/100 Patient-Years Followup",
     xlim = c(0,10))
polygon(c(res$time, rev(res$time)), c(100*res$upper, rev(100*res$lower)),
        col = mycol, border = NA)

dev.off()


# plot with shaded polygon for CI


# now to get the plot and results for each group

mcf_group = mcf(g ~ cohort_n, data = dt_mcf,
                level = 0.68, variance = "bootstrap")

# now going to again save the results so that it can be plotted and presented.

res_group = mcf_group@MCF
res_group0 = mcf_group@MCF %>% filter(cohort_n == 0)
res_group1 = mcf_group@MCF %>% filter(cohort_n == 1)
res_group2 = mcf_group@MCF %>% filter(cohort_n == 2)

write_csv(res_group, 
'P:/ORD_Perez_201602128D/Deo/CABG_HF/JTCVS_paper/group_mcf.csv')


# plot for MCF for each group.


tiff('P:/ORD_Perez_201602128D/Deo/CABG_HF/JTCVS_paper/mcf_group.tiff',
     height = 7, width = 5, units = "in", res = 1200)



plot(x = res_group0$time, y = res_group0$MCF*100, 
     type = "s", col = "black", ylim = c(0,120),
     xlab = "Followup Time:Years",
     ylab = "Event Rate/Per 100 Patient-Years Followup")

polygon(c(res_group0$time, rev(res_group0$time)), c(res_group0$lower*100, rev(res_group0$upper*100)),
        col = t_col("black"), border = NA)

lines(x = res_group1$time, y = res_group1$MCF*100, col = "blue", lty = 1)


polygon(c(res_group1$time, rev(res_group1$time)), 
        c(res_group1$lower*100, rev(res_group1$upper*100)),
        col = t_col("blue"), border = NA)

lines(x = res_group2$time, y = res_group2$MCF*100, col = "red", lty = 1)

polygon(c(res_group2$time, rev(res_group2$time)), 
        c(res_group2$lower*100, rev(res_group2$upper*100)),
        col = t_col("red"), border = NA)

dev.off()

#######################################


detach(dt_mcf)


# now going to use the dt16 data to do the model for CHF events
# going to use the Andersen Gill model.
# first make a simple df of only the var needed for the  model.

dt16$cohort_n[dt16$cohort == "NORMAL"]<- 0
dt16$cohort_n[dt16$cohort == "MID"]<- 1
dt16$cohort_n[dt16$cohort == "LOW"]<- 2

dt16$cohort_n <- factor(dt16$cohort_n)

# now get the variables 

ag = dt16 %>% dplyr::select(scrssn, id, tstart, tstop,chf,
                            obese, anemia, prior_mi, priorpci, cr, alb, age, 
                            ltm, sex, diabetes, cohort_n, pvd)

# see missing infomation 

describe(ag)

miss_var_summary(ag)

# variable n_miss pct_miss
# <chr>     <int>    <dbl>
#   1 alb         877  11.2   
# 2 ltm         549   7.02  
# 3 cr           10   0.128 
# 4 anemia        9   0.115 
# 5 pvd           4   0.0512
# 6 scrssn        0   0     
# 7 id            0   0     
# 8 tstart        0   0     
# 9 tstop         0   0     
# 10 obese         0   0     
# 11 prior_mi      0   0     
# 12 priorpci      0   0     
# 13 age           0   0     
# 14 sex           0   0     
# 15 diabetes      0   0     
# 16 cohort_n      0   0     

# am going to do simple imputation for albumin and other variables 

mean(ag$alb, na.rm = T)

ag$alb[is.na(ag$alb)]<- 3.81


ag$anemia[is.na(ag$anemia)]<- 0
ag$pvd[is.na(ag$pvd)]<- 0
ag$lmcad[is.na(ag$lmcad)]<- 0

# use creatinine rather than CKD

describe(ag$cr)

ag$cr <- with(ag, ifelse(cr > 5, 5, cr))
ag$cr[is.na(ag$cr)]<- 1.17

ag$ltm[is.na(ag$ltm)]<- 0

ag$diabetes = factor(ag$diabetes)

# now to fit the AG model for all variables together. 
# fit age, cr, albumin as splines in the model with 3 df each.

# ag model:

andersen = coxph(Surv(tstart, tstop, chf) ~ obese + anemia +  prior_mi + 
      priorpci +  ns(cr,3) + ns(alb,3)  + ns(age,3) + 
           ltm +  sex +  diabetes +  cohort_n +  pvd,
      method = "breslow",
      data = ag)


summary(andersen)


# marginal means model:


marg = coxph(Surv(tstart, tstop, chf) ~ pspline(cr,df = 3) + pspline(alb,df = 3)  + 
               pspline(age,df = 3) + obese + anemia +  prior_mi + 
             priorpci +   + 
             ltm +  sex +  diabetes +  cohort_n +  pvd + cluster(id),
           data = ag)

marg

summary(marg)

# all the cont vars have important spline terms, but for the HF group, albumin would be 
# very important 


t = termplot(marg, term = 2, se = TRUE
,plot = F)

t_alb = t$alb %>% tbl_df()



t_alb$HR = exp(t_alb$y)
t_alb$lower = exp(t_alb$y - t_alb$se)
t_alb$upper = exp(t_alb$y + t_alb$se)



tiff('P:/ORD_Perez_201602128D/Deo/CABG_HF/JTCVS_paper/alb_spline.tiff',
     height = 5, width = 5, units = "in", res = 1200)



plot(x = t_alb$x, y = t_alb$HR, type = "l",
     col = "blue", size = 2, xlim = c(2,5),
     ylim = c(0.8, 2), lwd = 2,
    ylab = "Relative Hazard Ratio",
    xlab = "Serum Albumin (mg/dl)")
polygon(c(t_alb$x, rev(t_alb$x)),
        c(t_alb$lower, rev(t_alb$upper)),
        col = t_col("blue"), border = NA)
abline(h = 1, lty = 3, col = "black")

dev.off()

###_ using plotHR

library(Greg)

plotHR(marg, term = 2, se = T)

# result for the marginal means model:

# coef      se(coef) se2       Chisq            DF                    p      
# pspline(cr, df = 3), line -0.141577 0.061160 0.0688214   5.36 1.00 2.1e-02
# pspline(cr, df = 3), nonl                               48.02 2.05 4.1e-11
# pspline(alb, df = 3), lin -0.279663 0.093770 0.0737044   8.89 1.00 2.9e-03
# pspline(alb, df = 3), non                               18.53 2.04 1.0e-04
# pspline(age, df = 3), lin  0.017069 0.005752 0.0087254   8.81 1.00 3.0e-03
# pspline(age, df = 3), non                               13.58 2.03 1.2e-03
# obese                      0.226418 0.096558 0.0587105   5.50 1.00 1.9e-02
# anemia                     0.290666 0.109520 0.0606472   7.04 1.00 8.0e-03
# prior_mi                   0.082962 0.100056 0.0585979   0.69 1.00 4.1e-01
# priorpci                   0.448361 0.280291 0.1284637   2.56 1.00 1.1e-01
# ltm                       -0.001325 0.001560 0.0009188   0.72 1.00 4.0e-01
# sex                       -0.336970 0.303987 0.3044707   1.23 1.00 2.7e-01
# diabetes1                  0.290296 0.136491 0.0797100   4.52 1.00 3.3e-02
# diabetes2                  0.562292 0.109513 0.0670463  26.36 1.00 2.8e-07
# cohort_n1                  1.399028 0.108307 0.0711074 166.85 1.00 3.6e-38
# cohort_n2                  1.998718 0.127766 0.0805560 244.72 1.00 3.7e-55
# pvd                        0.326793 0.101402 0.0597638  10.39 1.00 1.3e-03


# obese        1.2541     0.7974    1.0379     1.515
# anemia       1.3373     0.7478    1.0790     1.658
# prior_mi     1.0865     0.9204    0.8930     1.322
# priorpci     1.5657     0.6387    0.9039     2.712
# ltm          0.9987     1.0013    0.9956     1.002
# sex          0.7139     1.4007    0.3935     1.295
# diabetes1    1.3368     0.7480    1.0230     1.747
# diabetes2    1.7547     0.5699    1.4157     2.175
# cohort_n1    4.0513     0.2468    3.2764     5.009
# cohort_n2    7.3796     0.1355    5.7448     9.480
# pvd          1.3865     0.7212    1.1366     1.691
