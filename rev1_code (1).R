##################################################################
##               Outcome of CABG in HF with HFmEF               ##
##################################################################

# JTCVS
# revision 1
# 2020-12-14


# 1. redo figures for
# all-cause mortality, MI, CHF as CIF in the competing risk model.

# all-cause mortality.




library(easypackages)
libraries(c('tidyverse','rms','naniar','Hmisc',"rstpm2","etm","riskRegression",
            'MASS', 'tableone','haven',"lubridate", "survival", "mstate",
            "naniar", "survminer", "ggsci", "mice", "miceadds","psfmi"))

library(reda);library(reReg)

dt = 
  read_csv('P:/ORD_Perez_201602128D/Deo/CABG_HF/JTCVS_paper/data/dt.csv')

# first look at survival 

s = survfit(Surv(survyears, died) ~ cohort, data = dt)



tiff('P:\\ORD_Perez_201602128D\\Deo\\CABG_HF\\JTCVS_paper\\rev1\\figures_rev\\allcm.tiff',
     height = 7, width = 7, units = "in", res = 1200)



f1 <- ggsurvplot(s, conf.int = T, 
                 censor.size = 0,
                 risk.table = T,
                 surv.scale = "percent",
                 palette = "jama",
                 legend.labs = c("HFrEF", "HFmEF", "Control"),
                 ylab = "Survival",
                 xlab = "Followup Time:Years")



print(f1)

dev.off()

# MI - CIF plot.


# T_COL 
# t_col function for polygon code

t_col <- function(color, percent = 80, name = NULL) {
  #      color = color name
  #    percent = % transparency
  #       name = an optional name for the color
  
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  
  ## Save the color
  invisible(t.col)
}


mypal <- pal_jama("default", alpha = 0.8)(3)

mycol1 <- t_col(mypal[1])

mycol2 <- t_col(mypal[2])

mycol3 <- t_col(mypal[3])



# ETM TO DF

etm_to_df <- function(object, ci.fun = "cloglog", level = 0.95, ...) {
  l.X <- ncol(object$X)
  l.trans <- nrow(object[[1]]$trans)
  res <- list()
  for (i in seq_len(l.X)) {
    temp <- summary(object[[i]], ci.fun = ci.fun, level = level)
    res[[i]] <- data.table::rbindlist(
      temp[object$failcode + 1], idcol = "CIF"
    )[, CIF := paste0("CIF", CIF, "; ", names(object)[i])]
  }
  do.call(rbind, res)
}


df3 = read_csv("P:\\ORD_Perez_201602128D\\Deo\\CABG_HF\\JTCVS_paper\\mi_cif.csv")


df3 %>% count(cohort_n)
# 2 - hfref, 1  - hfref, 0 - normal.

describe(df3$event_y)

df3$event_y = df3$event_y + 0.01

# base R plot for the MI CIF plot.
# plot with 3 groups
# MI with 95% CI for each group.

library(etm)

mi_res = etmCIF(Surv(event_y, event != 0) ~ cohort_n, 
              failcode = 1,
              etype = event,
              data = df3)


mi_res2 = etm_to_df(mi_res)

glimpse(mi_res2)

mi_res2 %>% count(CIF)

mi_normal = mi_res2 %>% filter(CIF == "CIF0 1; cohort_n=0")

mi_hfmef = mi_res2 %>% filter(CIF == "CIF0 1; cohort_n=1")

mi_hfref = mi_res2 %>% filter(CIF == "CIF0 1; cohort_n=2")


# convert the mi_res to df to use for figures.


tiff('P:\\ORD_Perez_201602128D\\Deo\\CABG_HF\\JTCVS_paper\\rev1\\figures_rev\\mi_cif.tiff',
     height = 5, width = 7, units = "in", res = 1200)

plot(x = mi_normal$time, y = 100*mi_normal$P, type = "s",
     xlim = c(0,10), ylim = c(0,15),
     xlab = "Followup Time:Years",
     ylab = "Cumulative Incidence: Myocardial Infarction",
     frame.plot = F,
     col = mypal[1], yaxt = "n")

axis(2, at = c(0, 5, 10, 15),
     labels = c("0%","5%","10%","15%"),
     las = 0)

polygon(c(mi_normal$time, rev(mi_normal$time)), 
        c(mi_normal$lower*100, rev(mi_normal$upper*100)),
        col = mycol1, border = NA)

lines(mi_hfmef$time, 100*mi_hfmef$P, col = mypal[2], lwd = 1.5)


polygon(c(mi_hfmef$time, rev(mi_hfmef$time)), 
        c(mi_hfmef$lower*100, rev(mi_hfmef$upper*100)),
        col = mycol2, border = NA)


lines(mi_hfref$time, 100*mi_hfref$P, col = mypal[3], lwd = 1.5)

polygon(c(mi_hfref$time, rev(mi_hfref$time)), 
        c(mi_hfref$lower*100, rev(mi_hfref$upper*100)),
        col = mycol3, border = NA)

dev.off()


# now to that I have MI information am going to use this to 
# obtain unadjusted subHR for this model.
# using crr from the cmprsk package.

df3$cohort_n = factor(df3$cohort_n)

df3 %>% count(cohort_n)

mi_subhr <- crr(fstatus = df3$event,
                ftime = df3$event_y,
                failcode = 1,
                cencode = 0,
                cov1 = df3$cohort_n)


# obtain subHR for MI: 

df3 %>% count(event)

df3$event = factor(df3$event)

fgdata_mi <- finegray(Surv(event_y, event) ~ ., data = df3)

mi_fg <- coxph(Surv(fgstart, fgstop, fgstatus) ~ cohort_n, data = fgdata_mi)

summary(mi_fg)

# subHR compared to the Control group.

# Call:
#   coxph(formula = Surv(fgstart, fgstop, fgstatus) ~ cohort_n, data = fgdata_mi)
# 
# n= 210924, number of events= 317 
# 
# coef exp(coef) se(coef)      z Pr(>|z|)  
# cohort_n1  0.2465    1.2796   0.1213  2.033   0.0421 *
#   cohort_n2 -0.5565    0.5732   0.2600 -2.140   0.0323 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# exp(coef) exp(-coef) lower .95 upper .95
# cohort_n1    1.2796     0.7815    1.0089    1.6229
# cohort_n2    0.5732     1.7445    0.3443    0.9542
# 
# Concordance= 0.544  (se = 0.014 )
# Likelihood ratio test= 11.47  on 2 df,   p=0.003
# Wald test            = 10.39  on 2 df,   p=0.006
# Score (logrank) test = 10.68  on 2 df,   p=0.005

#----------------------------------------------------------------#



# to get HFmEF vs HFrEF 

df3$cohort_r[df3$cohort_n == 1]<- 0
df3$cohort_r[df3$cohort_n == 2]<- 1
df3$cohort_r[df3$cohort_n == 0]<- 2

df3$cohort_r = factor(df3$cohort_r)

# FG HR for HFmEF vs HFrEF

df3$event = factor(df3$event)

fgdata_mi <- finegray(Surv(event_y, event) ~ ., data = df3)

mi_fg2 <- coxph(Surv(fgstart, fgstop, fgstatus) ~ cohort_r, data = fgdata_mi)0.44

summary(mi_fg2)

# 
# Call:
#   coxph(formula = Surv(fgstart, fgstop, fgstatus) ~ cohort_r, data = fgdata_mi)
# 
# n= 210924, number of events= 317 
# 
# coef exp(coef) se(coef)      z Pr(>|z|)   
# cohort_r1 -0.8030    0.4480   0.2686 -2.990  0.00279 ** --- HFmEF vs HFrEF (HFmEF index)
#   cohort_r2 -0.2465    0.7815   0.1213 -2.033  0.04206 * 
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
#            exp(coef) exp(-coef) lower .95 upper .95
# cohort_r1    0.4480      2.232    0.2646    0.7583
# cohort_r2    0.7815      1.280    0.6162    0.9912
# 
# Concordance= 0.544  (se = 0.014 )
# Likelihood ratio test= 11.47  on 2 df,   p=0.003
# Wald test            = 10.39  on 2 df,   p=0.006
# Score (logrank) test = 10.68  on 2 df,   p=0.005


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

# am going to save this dataset for chf as first event.

write_csv(dt3,
          'P:/ORD_Perez_201602128D/Deo/CABG_HF/JTCVS_paper/rev1/chf_cif.csv')


dt <- read_csv("P:/ORD_Perez_201602128D/Deo/CABG_HF/JTCVS_paper/rev1/chf_cif.csv")


# now to open this dataset and then cif for CHF.


chf_res = etmCIF(Surv(time_y, event != 0) ~ cohort_n, 
                failcode = 1,
                etype = event,
                data = dt)


chf_res2 = etm_to_df(chf_res)

glimpse(chf_res2)

chf_res2 %>% count(CIF)

chf_normal = chf_res2 %>% filter(CIF == "CIF0 1; cohort_n=0")

chf_hfmef = chf_res2 %>% filter(CIF == "CIF0 1; cohort_n=1")

chf_hfref = chf_res2 %>% filter(CIF == "CIF0 1; cohort_n=2")


# convert the mi_res to df to use for figures.


tiff('P:\\ORD_Perez_201602128D\\Deo\\CABG_HF\\JTCVS_paper\\rev1\\figures_rev\\chf_cif.tiff',
     height = 5, width = 7, units = "in", res = 1200)

plot(x = chf_normal$time, y = 100*chf_normal$P, type = "s",
     xlim = c(0,10), ylim = c(0,40),
     xlab = "Followup Time:Years",
     ylab = "Cumulative Incidence: CHF (1st Readmission)",
     frame.plot = F,
     col = mypal[1], yaxt = "n")

axis(2, at = c(0, 10, 20, 30, 40),
     labels = c("0%","10%","20%","30%", "40%"),
     las = 0)

polygon(c(chf_normal$time, rev(chf_normal$time)), 
        c(chf_normal$lower*100, rev(chf_normal$upper*100)),
        col = mycol1, border = NA)

lines(chf_hfmef$time, 100*chf_hfmef$P, col = mypal[2], lwd = 1.5)


polygon(c(chf_hfmef$time, rev(chf_hfmef$time)), 
        c(chf_hfmef$lower*100, rev(chf_hfmef$upper*100)),
        col = mycol2, border = NA)


lines(chf_hfref$time, 100*chf_hfref$P, col = mypal[3], lwd = 1.5)

polygon(c(chf_hfref$time, rev(chf_hfref$time)), 
        c(chf_hfref$lower*100, rev(chf_hfref$upper*100)),
        col = mycol3, border = NA)

dev.off()


#################################################################
##                        CPH Model: with multiple imputation  ##
#################################################################

# Need to impute using mice and then coxph model.
# Q regarding variable selection.
# all variables included, no issue with overfit as plenty of events
# adequate sample size.


dt = 
  read_csv('P:/ORD_Perez_201602128D/Deo/CABG_HF/JTCVS_paper/data/dt.csv')


impvars <- c('renfail', "age", "csmok", "diabetes", "prior_mi", "sex", "mitreg", 
             "pvd", "priorhs", "priorstroke", "priorpci", "cr", "cohort", "bmi", "obese",
             "alb","hgb", "anemia", "race.mod")



# look at missing for these variables 

describe(dt[, c(impvars)]) # see missing data # very minimal here

library(naniar)

miss_var_summary(dt[,c(impvars)]) # good code to list missing % in order.

# need to identify CKD using eGFR


dt$race_n <- with(dt, ifelse(race.mod == 'white', 1, 0
))


gfr <- function(age, scr,sex, race){
  male <- age^(-0.203)*scr^(-1.154)*175
  female <- age^(-0.203)*scr^(-1.154)*175*0.742
  
  a <- ifelse(sex == 1, female , male)
  b <- ifelse(race == 1, a, a*1.212)
  return(b)
}

dt$gfr = with(dt, gfr(age = age, scr = cr,
                      sex = sex, race = race_n))

# using gfr to divide into CKD groups

dt$ckd = with(dt, ifelse(gfr > 60, 1, 0))

# add ckd to table 1

ckd_table = tableone::CreateCatTable(vars = c('ckd'), strata = c("cohort"),
                                     data = dt)

ckd_table

# need to make some changes before the regression 

dt$smoking = with(dt, ifelse(csmok == 0, 0 , 1))

dt %>% count(smoking)

dt$diabetes_f = with(dt, ifelse(diabetes == 0, 0, 1))

dt %>% count(diabetes_f)

describe(dt$ltm)

dt$lmcad = with(dt, ifelse(ltm > 50, 1, 0))

dt %>% count(lmcad)

describe(dt$ckd)


# now to see the variables and then MI for them.
# this was the prior model that was presented.
# 
# mw2 <- coxph(Surv(survyears, died) ~ pspline(cr,df=3) + cohort_n + age + diabetes + prior_mi + 
#                sex + pvd + priorstroke + priorpci +  obese + 
#                anemia + lmcad + alb+ race.mod , data = dt)
# 


impvars <- c('renfail', "age", "csmok", "diabetes", "prior_mi", "sex", "mitreg", 
             "pvd", "priorhs", "priorstroke", "priorpci", "cr", "cohort", "bmi", "obese",
             "alb","hgb", "anemia", "ckd", "lmcad", "diabetes_f", "race.mod",
             )



# look at missing for these variables 

describe(dt[, c(impvars)]) # see missing data # very minimal here

library(naniar)

miss_var_summary(dt[,c(impvars)]) # good code to list missing % in order.

# using mice to impute and create 10 datasets.
# mice, miceadds used to create the imputed dataasets.
# using the MI datasets, psfmi used to do the Cox model and then pool the results 
# contains 1 cat variable, other factors and then continuous variables too.

# create a separate dataset with only the variables needed for the model.
# need to convert race.mod to race_i

dt$race_i[dt$race.mod == "black"]<- 0
dt$race_i[dt$race.mod == "others"]<- 1
dt$race_i[dt$race.mod == "white"]<- 2

keep <- c("survyears", 'renfail', "age", "csmok", "prior_mi", "sex", "mitreg", 
             "pvd", "priorhs", "priorstroke", "priorpci", "cr",  "obese",
             "alb", "anemia", "ckd", "lmcad", "diabetes_f", "cohort", "died",
          "race_i")


dti <- dt[, c(keep)]

glimpse(dti)

dti$cohort_n[dti$cohort == "NORMAL"]<- 0
dti$cohort_n[dti$cohort == "MID"]<- 1
dti$cohort_n[dti$cohort == "LOW"]<- 2

dti2 = dti %>% dplyr::select(-cohort)

glimpse(dti2)

miss_var_summary(dti2)

# now to use mice to create 10 imputed datasets.
# before that need to identify type of variables.

dti2$ckd <- factor(dti2$ckd)
dti2$pvd <- factor(dti2$pvd)
dti2$anemia <- factor(dti2$anemia)
dti2$lmcad <- factor(dti2$lmcad)
dti2$race_i <- factor(dti2$race_i)
dti2$prior_mi <- factor(dti2$prior_mi)
dti2$priorpci <- factor(dti2$priorpci)
dti2$priorhs <- factor(dti2$priorhs)
dti2$priorstroke <- factor(dti2$priorstroke)
dti2$diabetes_f <- factor(dti2$diabetes_f)
dti2$renfail <- factor(dti2$renfail)


dti2$cohort_n = factor(dti2$cohort_n)


dti2$cr <- with(dti2, ifelse(cr < 0.8, 0.8, 
                             ifelse(cr > 5, 5, cr)))


describe(dti2$cr)

describe(dti2$alb)


# mitreg is categorical 

dti2$mitreg <- factor(dti2$mitreg)

dti2$survyears <- as.numeric(dti2$survyears)


dti2 <- data.frame(dti2)


dti2$nelson <- nelsonaalen(data = dti2,
                           timevar = survyears,
                           statusvar = died)


dti2$cohort_r[dti2$cohort_n == 1]<- 0
dti2$cohort_r[dti2$cohort_n == 2]<- 1
dti2$cohort_r[dti2$cohort_n == 0]<- 2

dti2$cohort_r = factor(dti2$cohort_r)


mi10 <- mice(data = dti2, m = 10, seed = 1974)

mi10

mi10$imp$cr # see creatinine imputed values.
mi10$imp$alb # see albumin imputed values.

# plot for creatinine

tiff('P:\\ORD_Perez_201602128D\\Deo\\CABG_HF\\JTCVS_paper\\rev1\\figures_rev\\mi_cr.tiff',
     height = 5, width = 7, units = "in", res = 1200)

stripplot(mi10, cr, pch = 19,
          xlab = "Imputed Dataset",
          ylab = "Serum Creatinine (mg/dl)") # plot for creatinine

dev.off()


# plot for albumin.

tiff('P:\\ORD_Perez_201602128D\\Deo\\CABG_HF\\JTCVS_paper\\rev1\\figures_rev\\mi_alb.tiff',
     height = 5, width = 7, units = "in", res = 1200)


stripplot(mi10, alb, pch = 19,
         xlab = "Imputed Dataset",
         ylab = "Serum Albumin (mg/dl)") # plot for albumin

dev.off()


# looking categorical variable like mitreg.

stripplot(mi10, mitreg)

# now am going to use the psfmi package to model the data.
# will try with routine mice package first.

model <- with(mi10, coxph(Surv(survyears, died) ~ age + sex + cohort_n + 
        race_i + obese + cr + alb + mitreg + anemia + lmcad + 
        ckd + diabetes_f + prior_mi + priorhs + priorpci +  renfail + 
          csmok + priorpci + priorstroke))

summary(pool(model), conf.int = T, exponentiate = T)

# > summary(pool(model), conf.int = T, exponentiate = T)
#            term  estimate   std.error   statistic        df      p.value     2.5 %    97.5 %
# 1           age 1.0403543 0.003295465 12.00479342 1735.5088 0.000000e+00 1.0336517 1.0471005
# 2           sex 0.9919996 0.245555127 -0.03271183 1744.2726 9.739081e-01 0.6128433 1.6057338
# 3     cohort_n1 1.3109001 0.055729908  4.85760767 1737.7984 1.294850e-06 1.1751658 1.4623120
# 4     cohort_n2 1.4608228 0.081031034  4.67721844 1674.9678 3.141290e-06 1.2461603 1.7124629
# 5       race_i1 0.9037174 0.103360394 -0.97947182 1745.8467 3.274826e-01 0.7378891 1.1068129
# 6       race_i2 1.0987170 0.085987058  1.09485202 1739.4583 2.737330e-01 0.9282021 1.3005562
# 7         obese 0.8465776 0.051009991 -3.26511304 1692.5400 1.116058e-03 0.7659773 0.9356591
# 8            cr 1.2527282 0.037125846  6.06918746 1721.6656 1.577654e-09 1.1647509 1.3473507
# 9           alb 0.6575985 0.056140043 -7.46634017  240.9692 1.502798e-12 0.5887530 0.7344945
# 10      mitreg1 1.1941870 0.059475763  2.98383094  548.7756 2.973327e-03 1.0625138 1.3421780
# 11      mitreg2 1.4743590 0.095069683  4.08356583  450.4889 5.248480e-05 1.2230995 1.7772344
# 12      mitreg3 2.0676076 0.142480590  5.09818351  211.0323 7.611325e-07 1.5613091 2.7380877
# 13      anemia1 1.2770728 0.053118924  4.60420880 1699.7254 4.447884e-06 1.1507170 1.4173033
# 14       lmcad1 1.0385762 0.058215255  0.65018551  386.8514 5.159583e-01 0.9262536 1.1645196
# 15         ckd1 0.9781641 0.066067557 -0.33417077 1730.7606 7.382912e-01 0.8592817 1.1134940
# 16  diabetes_f1 1.3041012 0.050451463  5.26276289 1670.0690 1.602691e-07 1.1812336 1.4397491
# 17    prior_mi1 1.1830826 0.050271075  3.34433719 1715.5296 8.425655e-04 1.0719982 1.3056781
# 18     priorhs1 1.0318394 0.148190670  0.21150470 1748.2325 8.325181e-01 0.7715864 1.3798746
# 19    priorpci1 0.8018262 0.132716809 -1.66417067 1747.8386 9.625758e-02 0.6180636 1.0402250
# 20      renfail 4.2652213 0.134346643 10.79665291 1625.5674 0.000000e+00 3.2771800 5.5511486
# 21        csmok 0.9967359 0.019705650 -0.16591542 1748.5593 8.682427e-01 0.9589478 1.0360130
# 22 priorstroke1 0.8684803 0.164956528 -0.85483343 1750.8647 3.927603e-01 0.6284223 1.2002407


# select only 1 dataset to run model and test for PH.


long <- complete(mi10, 'long', inc = TRUE)

df1 <- complete(mi10, 1)

model.df1 <- coxph(Surv(survyears, died) ~ age + sex + cohort_n + 
        race_i + obese + cr + alb + mitreg + anemia + lmcad + 
        ckd + diabetes_f + prior_mi + priorhs + priorpci +  renfail + 
        csmok + priorpci + priorstroke, data = df1)


cox.zph(model.df1)

# > cox.zph(model.df1)
# chisq df       p
# age          1.0939  1  0.2956
# sex          0.0574  1  0.8107
# cohort_n     3.3187  2  0.1903
# race_i       3.1654  2  0.2054
# obese        0.1557  1  0.6931
# cr           0.0159  1  0.8996
# alb          4.5820  1  0.0323
# mitreg       1.6139  3  0.6562
# anemia       2.8310  1  0.0925
# lmcad        6.0166  1  0.0142
# ckd          0.0136  1  0.9073
# diabetes_f   4.3971  1  0.0360
# prior_mi     0.4201  1  0.5169
# priorhs      6.9619  1  0.0083
# priorpci     1.9652  1  0.1610
# renfail     43.4609  1 4.3e-11
# csmok        0.7366  1  0.3907
# priorstroke  1.0047  1  0.3162
# GLOBAL      74.9103 22 1.1e-07



# now to obtain pairwise comparison for the main model here.
# to obtain the pairwise comparison here, will need to now make HFmEF as the control group.


model_r <- with(mi10, coxph(Surv(survyears, died) ~ age + sex + cohort_r + 
                            race_i + obese + cr + alb + mitreg + anemia + lmcad + 
                            ckd + diabetes_f + prior_mi + priorhs + priorpci +  renfail + 
                            csmok + priorpci + priorstroke))

summary(pool(model_r), conf.int = T, exponentiate = T)

# HFmEF vs HFrEF

# cohort_r1 1.1074424 (0.94 - 1.29); p = 0.2

# now to see the # of grafts/ # radial used / # complete revasc information.

glimpse(dt)

dt %>% count(grafts)

dt$grafts[is.na(dt$grafts)]<- 3

prop.table(table(dt$cohort, dt$grafts), 1)

# > prop.table(table(dt$cohort, dt$grafts), 1)
# 
# 1          2          3 more than 3
# LOW    0.07773852 0.13780919 0.66254417  0.12190813
# MID    0.07346939 0.14169096 0.65947522  0.12536443
# NORMAL 0.07455315 0.15945437 0.65286924  0.11312324

CreateCatTable(vars = 'grafts', data = dt, strata = c("cohort"))


# > CreateCatTable(vars = 'grafts', data = dt, strata = c("cohort"))
# Stratified by cohort
#              LOW         MID          NORMAL       p      test
# n              566         1715         4252                    
# grafts (%)                                            0.498     
# 1            44 ( 7.8)   126 ( 7.3)   317 ( 7.5)             
# 2            78 (13.8)   243 (14.2)   678 (15.9)             
# 3           375 (66.3)  1131 (65.9)  2776 (65.3)             
# more than 3  69 (12.2)   215 (12.5)   481 (11.3)             



art <- c("33534","33535","33536")

dt$art <- with(dt, ifelse((cpt01 %in% art|cpt02 %in% art|cpt03 %in% art), 1, 0))

dt %>% count(art)

prop.table(table(dt$art))

CreateCatTable(vars = 'art', data = dt, strata = c("cohort"))

# > CreateCatTable(vars = 'art', data = dt, strata = c("cohort"))
# Stratified by cohort
#             LOW       MID        NORMAL     p      test
# n           566       1715       4252                  
# art = 1 (%) 33 (5.8)  102 (5.9)  275 (6.5)   0.680

# calculate the LOS overall and per group.

glimpse(dt)

# no direct method to calculate LOS, so going to use surgdate and disd.

dt$surgdate <- as_date(dt$surgdate)

dt$disd <- as_date(dt$disd)

dt$po_days = (dt$surgdate %--% dt$disd)/ddays(1)

describe(dt$po_days)

quantile(dt$po_days, 0.99, na.rm = T)

dt$po_days = with(dt, ifelse(po_days > 45, 45, po_days))

los <- CreateContTable(vars = 'po_days', data = dt, strata = c("cohort"))

print(los, nonnormal = 'po_days')

# > print(los, nonnormal = 'po_days')
# Stratified by cohort
#                        LOW                MID                NORMAL           
# n                      566                1715               4252             
# po_days (median [IQR]) 8.00 [6.00, 14.00] 7.00 [6.00, 11.00] 7.00 [5.00, 9.00]
# Stratified by cohort
# p      test   
# n                                    
# po_days (median [IQR]) <0.001 nonnorm


# shapiro wilk test for all the cont vars from table 1.

# age 

qqnorm(dt$age)
qqline(dt$age)

# serum creat.

qqnorm(dt$cr)
qqline(dt$cr)

# bmi 

hist(dt$bmi, breaks = 100)

qqnorm(dt$bmi)
qqline(dt$bmi)

# serum albumin 

hist(dt$alb, breaks = 100)

# hb

hist(dt$hgb, breaks = 100)

qqnorm(dt$hgb)
qqline(dt$hgb)

# none of the cont variables are normally distributed; hence use median + IQR.

table(dt$curdiur)

prop.table(table(dt$cohort, dt$curdiur), 1)


CreateCatTable(vars = 'curdiur', data = dt, strata = c("cohort"))


# 
# Stratified by cohort
#                 LOW         MID          NORMAL   p      test
# n               566         1715         4252                
# curdiur = 1 (%) 467 (82.5)  1228 (71.6)  0 (0.0)  <0.001     
# Warning messages:


# C index for the cox model.
# use cph on the dataset 1 from the MI data.

m_h <- cph(Surv(survyears, died) ~ age + sex + cohort_n + 
             race_i + obese + cr + alb + mitreg + anemia + lmcad + 
             ckd + diabetes_f + prior_mi + priorhs + priorpci +  renfail + 
             csmok + priorpci + priorstroke, data = df1, x = T, y = T)

m_h


# Model Tests    Discrimination    
# Indexes    
# Obs     6533    LR chi2     888.69    R2       0.129    
# Events  1780    d.f.            22    Dxy      0.397    
# Center 1.705    Pr(> chi2)  0.0000    g        0.730    
# Score chi2 1195.64    gr       2.074    
# Pr(> chi2)  0.0000  


t <- cox.zph(m_h,transform = "km")

# > t
# chisq df       p
# age         1.09e+00  1  0.2975
# sex         3.42e-02  1  0.8532
# cohort_n    3.42e+00  2  0.1805
# race_i      3.10e+00  2  0.2126
# obese       7.89e-02  1  0.7788
# cr          6.56e-04  1  0.9796
# alb         3.46e+00  1  0.0629
# mitreg      2.79e+00  3  0.4244
# anemia      2.72e+00  1  0.0993
# lmcad       6.04e+00  1  0.0140
# ckd         1.64e-03  1  0.9677
# diabetes_f  4.45e+00  1  0.0348
# prior_mi    3.33e-01  1  0.5641
# priorhs     6.72e+00  1  0.0095
# priorpci    1.75e+00  1  0.1862
# renfail     4.40e+01  1 3.3e-11
# csmok       7.35e-01  1  0.3912
# priorstroke 1.03e+00  1  0.3106
# GLOBAL      7.35e+01 22 1.9e-07
 

### medications :-
  
# get each medication and then identify
# will need cohort crosswalk too to get patientsid.

cw <- read_sas("P:\\ORD_Perez_201602128D\\Deo\\CABG_HF\\data\\deo_cwalk.sas7bdat")

# now to limit to only my patients.

cw$mine <- with(cw, ifelse(ScrSSN %in% dt$scrssn, 1, 0))

cw2 = cw %>% filter(mine == 1)

psid = cw2$PatientSID

# now to get ap data.

ap = read_sas('P:\\ORD_Perez_201602128D\\Deo\\CABG_HF\\JTCVS_paper\\rev1\\data\\aptherapy.sas7bdat')

# limit to only my patients.

ap$mine = with(ap, ifelse(PatientSID %in% psid, 1, 0))

ap %>% count(mine)

ap2 = ap %>% filter(mine == 1)

ap2$ActionDateTime = as_date(ap2$ActionDateTime)

dt$surgdate = as_date(dt$surgdate)

dt$disd = as_date(dt$disd)

names(cw2) = tolower(names(cw2))

names(ap2) = tolower(names(ap2))

ap3 = left_join(ap2, cw2, by = "patientsid")

date = dt %>% dplyr::select(scrssn, disd, surgdate)

ap4 = left_join(ap3, date, by = "scrssn")

ap4$after = with(ap4, ifelse(actiondatetime > surgdate, 1, 0))

ap5 = ap4 %>% filter(after == 1)

ap5$before = with(ap5, ifelse(actiondatetime <= disd, 1, 0))

ap6 = ap5 %>% filter(before == 1)

ap_scrssn = ap6$scrssn

# add anti-platelet therapy to the main dataset.

dt$antiplatelet = with(dt, ifelse(scrssn %in% ap_scrssn, 1, 0))

dt %>% count(antiplatelet)

# now to go with statin therapy in the same manner...


st = read_sas('P:\\ORD_Perez_201602128D\\Deo\\CABG_HF\\JTCVS_paper\\rev1\\data\\statintherapy.sas7bdat')


st$mine = with(st, ifelse(PatientSID %in% psid, 1, 0))

st %>% count(mine)

st2 = st %>% filter(mine == 1)

st2$ActionDateTime = as_date(st2$ActionDateTime)

dt$surgdate = as_date(dt$surgdate)

dt$disd = as_date(dt$disd)

names(cw2) = tolower(names(cw2))

names(st2) = tolower(names(st2))

st3 = left_join(st2, cw2, by = "patientsid")

date = dt %>% dplyr::select(scrssn, disd, surgdate)

st4 = left_join(st3, date, by = "scrssn")

st4$after = with(st4, ifelse(actiondatetime > surgdate, 1, 0))

st5 = st4 %>% filter(after == 1)

st5$before = with(st5, ifelse(actiondatetime <= disd, 1, 0))

st6 = st5 %>% filter(before == 1)

st_scrssn = st6$scrssn

# add anti-lipid therapy to the main dataset.

dt$antilipid = with(dt, ifelse(scrssn %in% st_scrssn, 1, 0))

dt %>% count(antilipid)


# beta blocker therapy


bb = 
read_sas('P:\\ORD_Perez_201602128D\\Deo\\CABG_HF\\JTCVS_paper\\rev1\\data\\bbtherapy.sas7bdat')


bb$mine = with(bb, ifelse(PatientSID %in% psid, 1, 0))

bb %>% count(mine)

bb2 = bb %>% filter(mine == 1)

bb2$ActionDateTime = as_date(bb2$ActionDateTime)

dt$surgdate = as_date(dt$surgdate)

dt$disd = as_date(dt$disd)

names(cw2) = tolower(names(cw2))

names(bb2) = tolower(names(bb2))

bb3 = left_join(bb2, cw2, by = "patientsid")

date = dt %>% dplyr::select(scrssn, disd, surgdate)

bb4 = left_join(bb3, date, by = "scrssn")

bb4$after = with(bb4, ifelse(actiondatetime > surgdate, 1, 0))

bb5 = bb4 %>% filter(after == 1)

bb5$before = with(bb5, ifelse(actiondatetime <= disd, 1, 0))

bb6 = bb5 %>% filter(before == 1)

bb_scrssn = bb6$scrssn

# add betablocker therapy to the main dataset.

dt$betablocker = with(dt, ifelse(scrssn %in% bb_scrssn, 1, 0))

dt %>% count(betablocker)

# ACE/ARB therapy 



ac = 
  read_sas('P:\\ORD_Perez_201602128D\\Deo\\CABG_HF\\JTCVS_paper\\rev1\\data\\acetherapy.sas7bdat')


ac$mine = with(ac, ifelse(PatientSID %in% psid, 1, 0))

ac %>% count(mine)

ac2 = ac %>% filter(mine == 1)

ac2$ActionDateTime = as_date(ac2$ActionDateTime)

dt$surgdate = as_date(dt$surgdate)

dt$disd = as_date(dt$disd)

names(cw2) = tolower(names(cw2))

names(ac2) = tolower(names(ac2))

ac3 = left_join(ac2, cw2, by = "patientsid")

date = dt %>% dplyr::select(scrssn, disd, surgdate)

ac4 = left_join(ac3, date, by = "scrssn")

ac4$after = with(ac4, ifelse(actiondatetime > surgdate, 1, 0))

ac5 = ac4 %>% filter(after == 1)

ac5$before = with(ac5, ifelse(actiondatetime <= disd, 1, 0))

ac6 = ac5 %>% filter(before == 1)

ac_scrssn = ac6$scrssn

# add betablocker therapy to the main dataset.

dt$ace = with(dt, ifelse(scrssn %in% ac_scrssn, 1, 0))

dt %>% count(ace)

# spironolactone



sp = 
  read_sas('P:\\ORD_Perez_201602128D\\Deo\\CABG_HF\\JTCVS_paper\\rev1\\data\\spirotherapy.sas7bdat')


sp$mine = with(sp, ifelse(PatientSID %in% psid, 1, 0))

sp %>% count(mine)

sp2 = sp %>% filter(mine == 1)

sp2$ActionDateTime = as_date(sp2$ActionDateTime)

dt$surgdate = as_date(dt$surgdate)

dt$disd = as_date(dt$disd)

names(cw2) = tolower(names(cw2))

names(sp2) = tolower(names(sp2))

sp3 = left_join(sp2, cw2, by = "patientsid")

date = dt %>% dplyr::select(scrssn, disd, surgdate)

sp4 = left_join(sp3, date, by = "scrssn")

sp4$after = with(sp4, ifelse(actiondatetime > surgdate, 1, 0))

sp5 = sp4 %>% filter(after == 1)

sp5$before = with(sp5, ifelse(actiondatetime <= disd, 1, 0))

sp6 = sp5 %>% filter(before == 1)

sp_scrssn = sp6$scrssn

# add betablocker therapy to the main dataset.

dt$spiro = with(dt, ifelse(scrssn %in% sp_scrssn, 1, 0))

dt %>% count(spiro)

# now to save this dataset in the rev1 folder.

write_csv(dt,
  'P:\\ORD_Perez_201602128D\\Deo\\CABG_HF\\JTCVS_paper\\rev1\\data\\dt_meds.csv')

# now am going to present the medications in the table 1.

meds = CreateCatTable(vars = c('spiro',"ace","betablocker","antilipid","antiplatelet"),
               data = dt,
               strata = c("cohort"))


# add the medications to the table 1.


# > meds
# Stratified by cohort
#                      LOW         MID          NORMAL       p      test
# n                    566         1715         4252                    
# spiro = 1 (%)         82 (14.5)   100 ( 5.8)    90 ( 2.1)  <0.001     
# ace = 1 (%)          197 (34.8)   456 (26.6)   946 (22.2)  <0.001     
# betablocker = 1 (%)  524 (92.6)  1629 (95.0)  4099 (96.4)  <0.001     
# antilipid = 1 (%)    500 (88.3)  1529 (89.2)  3788 (89.1)   0.853     
# antiplatelet = 1 (%) 548 (96.8)  1672 (97.5)  4210 (99.0)  <0.001


# pairwise comparison --- 
# run model for Cox with HFmEF as control now.
# pairwise comparisons done and results pasted at the the correct locations in the rscript.

# now for the HF analysis.

##################################################################
##              CHF marginal means model: Pairwise              ##
##################################################################


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

plot(g)


# obtain overall rates for CHF readmissions

mcf_overall = mcf(g ~ 1, data = dt_mcf,
                  variance = "bootstrap", level = 0.95)




str(mcf_overall)

# am going to extract the information from MCF and then create a plot

res = mcf_overall@MCF

# res is actually a dataframe, so am going to save this for further use

write_csv(res, 
          'P:/ORD_Perez_201602128D/Deo/CABG_HF/JTCVS_paper/rev1/overall_mcf.csv')

# now this can be used to plot the graphs & also get estimates.

# can use this to plot as base R or ggplot2
# excellent plot - am going to provide the estimates from the data and then
# paste them into the plot later.


tiff('P:/ORD_Perez_201602128D/Deo/CABG_HF/JTCVS_paper/rev1/figures_rev/mcf_overall.tiff',
     height = 7, width = 5, units = "in", res = 1200)

plot(x = res$time, y = 100*res$MCF, type = "s",
     xlab = "Follow-up Time:Years",
     ylab = "Event Rate/100 Patient-Years Followup",
     xlim = c(0,10),
     ylim = c(0,40))

polygon(c(res$time, rev(res$time)), c(100*res$upper, rev(100*res$lower)),
        col = t_col("gray"), border = NA)

dev.off()


# plot with shaded polygon for CI


# now to get the plot and results for each group

mcf_group = mcf(g ~ cohort_n, data = dt_mcf,
                level = 0.95, variance = "bootstrap")

# now going to again save the results so that it can be plotted and presented.

res_group = mcf_group@MCF
res_group0 = mcf_group@MCF %>% filter(cohort_n == 0)
res_group1 = mcf_group@MCF %>% filter(cohort_n == 1)
res_group2 = mcf_group@MCF %>% filter(cohort_n == 2)

write_csv(res_group, 
          'P:/ORD_Perez_201602128D/Deo/CABG_HF/JTCVS_paper/group_mcf.csv')


# plot for MCF for each group.


tiff('P:/ORD_Perez_201602128D/Deo/CABG_HF/JTCVS_paper/rev1/mcf_group.tiff',
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
                            ltm, sex, diabetes, cohort_n, pvd,                            
                            priorstroke, priorpci, mitreg,csmok, race.mod,
                            priorhs, anemia, csmok)


# need to make ckd variable again.
# make egfr first and then ckd.



# need to identify CKD using eGFR


ag$race_n <- with(ag, ifelse(race.mod == 'white', 1, 0
))


gfr <- function(age, scr,sex, race){
  male <- age^(-0.203)*scr^(-1.154)*175
  female <- age^(-0.203)*scr^(-1.154)*175*0.742
  
  a <- ifelse(sex == 1, female , male)
  b <- ifelse(race == 1, a, a*1.212)
  return(b)
}

ag$gfr = with(ag, gfr(age = age, scr = cr,
                      sex = sex, race = race_n))

# using gfr to divide into CKD groups

ag$ckd = with(ag, ifelse(gfr > 60, 1, 0))

ag$lmcad = with(ag, ifelse(ltm > 50, 1, 0))

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


# plan to do mice imputation here too ...

ag$mitreg = factor(ag$mitreg)
ag$ltm = factor(ag$ltm)
ag$ckd = factor(ag$ckd)
ag$anemia = factor(ag$anemia)
ag$pvd = factor(ag$pvd)
ag$lmcad = factor(ag$lmcad)

describe(ag$alb)
describe(ag$cr)

ag$cr = with(ag, ifelse(cr > 5, 5,
                        ifelse(cr < 0.8, 0.8, cr)))


# now to do the mice for the data.


ag2 = ag %>% dplyr::select(-scrssn, -ltm)


ag_imp = mice(m = 10, data = ag2, seed = 1974)


# MM model with MI:



marg_i = with(ag_imp, coxph(Surv(tstart, tstop, chf) ~ cr + alb  + 
               age + obese + anemia +  prior_mi + csmok + ckd + 
               priorpci + prior_mi  + priorstroke + priorpci + race.mod + 
               lmcad +  sex +  diabetes +  cohort_n +  pvd + cluster(id)))




marg_i

summary(pool(marg_i), conf.int = T, exponentiate = T)

# results of the CHF model with imputed datasets.

# > summary(pool(marg_i), conf.int = T, exponentiate = T)
# term  estimate   std.error  statistic        df      p.value
# 1              cr 0.8362034 0.053243226 -3.3597406 1253.3027 8.036218e-04
# 2             alb 0.7252198 0.070492767 -4.5576376   88.7081 1.646019e-05
# 3             age 1.0202064 0.003960011  5.0517362 1253.6317 5.026309e-07
# 4           obese 1.2722972 0.058502691  4.1164623 1251.5291 4.099288e-05
# 5         anemia1 1.3047652 0.061960277  4.2934458 1118.6975 1.911707e-05
# 6        prior_mi 1.1098567 0.058863545  1.7707210 1238.3513 7.685308e-02
# 7           csmok 0.9971445 0.023780457 -0.1202499 1250.6776 9.043045e-01
# 8            ckd1 0.6541124 0.074691443 -5.6830614 1241.5397 1.647218e-08
# 9        priorpci 1.5469206 0.128499973  3.3950689 1260.6639 7.074501e-04
# 10    priorstroke 0.9185884 0.203624845 -0.4170272 1245.9057 6.767303e-01
# 11 race.modothers 0.5282727 0.113701993 -5.6124145 1260.1997 2.451153e-08
# 12  race.modwhite 0.6243007 0.082760892 -5.6925811 1259.0294 1.555445e-08
# 13         lmcad1 1.0437412 0.066066175  0.6480100  685.9384 5.171953e-01
# 14            sex 0.6941334 0.303969404 -1.2010783 1261.2109 2.299463e-01
# 15       diabetes 1.3333827 0.033485554  8.5923341 1254.2917 0.000000e+00
# 16      cohort_n1 4.1242463 0.071087915 19.9314226 1259.8955 0.000000e+00
# 17      cohort_n2 7.2940025 0.081077287 24.5081267 1250.6511 0.000000e+00
# 18           pvd1 1.3419678 0.059816915  4.9172886 1256.4680 9.935955e-07
# 2.5 %    97.5 %
#   1  0.7532643 0.9282745
# 2  0.6304292 0.8342631
# 3  1.0123111 1.0281632
# 4  1.1343388 1.4270342
# 5  1.1554058 1.4734322
# 6  0.9888107 1.2457207
# 7  0.9516922 1.0447675
# 8  0.5649537 0.7573418
# 9  1.2022189 1.9904557
# 10 0.6160652 1.3696679
# 11 0.4226509 0.6602897
# 12 0.5307376 0.7343580
# 13 0.9167645 1.1883048
# 14 0.3823448 1.2601746
# 15 1.2486028 1.4239190
# 16 3.5873707 4.7414690
# 17 6.2213687 8.5515704
# 18 1.1933745 1.5090632
# > summary(pool(marg_i), conf.int = T, exponentiate = T)
#              term  estimate   std.error  statistic        df      p.value     2.5 %    97.5 %
# 1              cr 0.8362034 0.053243226 -3.3597406 1253.3027 8.036218e-04 0.7532643 0.9282745
# 2             alb 0.7252198 0.070492767 -4.5576376   88.7081 1.646019e-05 0.6304292 0.8342631
# 3             age 1.0202064 0.003960011  5.0517362 1253.6317 5.026309e-07 1.0123111 1.0281632
# 4           obese 1.2722972 0.058502691  4.1164623 1251.5291 4.099288e-05 1.1343388 1.4270342
# 5         anemia1 1.3047652 0.061960277  4.2934458 1118.6975 1.911707e-05 1.1554058 1.4734322
# 6        prior_mi 1.1098567 0.058863545  1.7707210 1238.3513 7.685308e-02 0.9888107 1.2457207
# 7           csmok 0.9971445 0.023780457 -0.1202499 1250.6776 9.043045e-01 0.9516922 1.0447675
# 8            ckd1 0.6541124 0.074691443 -5.6830614 1241.5397 1.647218e-08 0.5649537 0.7573418
# 9        priorpci 1.5469206 0.128499973  3.3950689 1260.6639 7.074501e-04 1.2022189 1.9904557
# 10    priorstroke 0.9185884 0.203624845 -0.4170272 1245.9057 6.767303e-01 0.6160652 1.3696679
# 11 race.modothers 0.5282727 0.113701993 -5.6124145 1260.1997 2.451153e-08 0.4226509 0.6602897
# 12  race.modwhite 0.6243007 0.082760892 -5.6925811 1259.0294 1.555445e-08 0.5307376 0.7343580
# 13         lmcad1 1.0437412 0.066066175  0.6480100  685.9384 5.171953e-01 0.9167645 1.1883048
# 14            sex 0.6941334 0.303969404 -1.2010783 1261.2109 2.299463e-01 0.3823448 1.2601746
# 15       diabetes 1.3333827 0.033485554  8.5923341 1254.2917 0.000000e+00 1.2486028 1.4239190
# 16      cohort_n1 4.1242463 0.071087915 19.9314226 1259.8955 0.000000e+00 3.5873707 4.7414690
# 17      cohort_n2 7.2940025 0.081077287 24.5081267 1250.6511 0.000000e+00 6.2213687 8.5515704
# 18           pvd1 1.3419678 0.059816915  4.9172886 1256.4680 9.935955e-07 1.1933745 1.5090632


# would need to again do the imputation prior to fitting this model.


ag$cohort_r[ag$cohort_n == 1]<- 0
ag$cohort_r[ag$cohort_n == 2]<- 1
ag$cohort_r[ag$cohort_n == 0]<- 2

ag$cohort_r = factor(ag$cohort_r)

# now to again do MI for this dataset.



ag2_pw = ag %>% dplyr::select(-scrssn, -ltm)


ag2_pw_imp = mice(m = 10, data = ag2_pw, seed = 1974)


# now to do the model for pairwise comparisons.



marg_pw_i = with(ag2_pw_imp, coxph(Surv(tstart, tstop, chf) ~ cr + alb  + 
age + obese + anemia +  prior_mi + csmok + ckd + 
priorpci + prior_mi  + priorstroke + priorpci + race.mod + 
lmcad +  sex +  diabetes +  cohort_n +  pvd + cluster(id)))




marg_pw_i

summary(pool(marg_pw_i), conf.int = T, exponentiate = T)

# results for the pairwise comparison.

# coxph(formula = Surv(tstart, tstop, chf) ~ cr + alb + age + obese + 
#         anemia + prior_mi + csmok + ckd + priorpci + priorstroke + 
#         race.mod + lmcad + sex + diabetes + cohort_n + pvd, cluster = id)
# 
# coef exp(coef)  se(coef) robust se      z        p
# cr             -0.180932  0.834492  0.053239  0.070660 -2.561 0.010449
# alb            -0.311791  0.732135  0.059042  0.096508 -3.231 0.001235
# age             0.020197  1.020402  0.003950  0.007434  2.717 0.006594
# obese           0.240625  1.272044  0.058402  0.096308  2.498 0.012472
# anemia1         0.273453  1.314495  0.061048  0.111750  2.447 0.014405
# prior_mi        0.108342  1.114428  0.058503  0.098770  1.097 0.272683
# csmok          -0.002950  0.997054  0.023699  0.037656 -0.078 0.937558
# ckd1           -0.421689  0.655938  0.074481  0.125052 -3.372 0.000746
# priorpci        0.430571  1.538135  0.128385  0.274794  1.567 0.117142
# priorstroke    -0.068156  0.934115  0.202946  0.382820 -0.178 0.858694
# race.modothers -0.630292  0.532436  0.113630  0.190644 -3.306 0.000946
# race.modwhite  -0.465322  0.627933  0.082679  0.161679 -2.878 0.004001
# lmcad1         -0.001321  0.998680  0.064843  0.100520 -0.013 0.989513
# sex            -0.346920  0.706862  0.303815  0.306541 -1.132 0.257750
# diabetes        0.286542  1.331815  0.033393  0.054250  5.282 1.28e-07
# cohort_n1       1.415227  4.117421  0.071091  0.106091 13.340  < 2e-16
# cohort_n2       1.990755  7.321058  0.080881  0.128256 15.522  < 2e-16
# pvd1            0.295432  1.343707  0.059780  0.102527  2.882 0.003958
# 
# Likelihood ratio test=1430  on 18 df, p=< 2.2e-16
# n= 7815, number of events= 1282 
# 
# 
# > summary(pool(marg_pw_i), conf.int = T, exponentiate = T)
# term  estimate   std.error  statistic        df      p.value
# 1              cr 0.8374220 0.053275784 -3.3303536 1254.5827 8.925740e-04
# 2             alb 0.7306507 0.062753241 -5.0008532  355.4036 8.985592e-07
# 3             age 1.0205467 0.003957088  5.1397656 1257.5492 3.186624e-07
# 4           obese 1.2719393 0.058434023  4.1164848 1259.0349 4.097357e-05
# 5         anemia1 1.3131672 0.061281464  4.4457479 1243.6466 9.539967e-06
# 6        prior_mi 1.1128804 0.058701334  1.8219619 1249.9174 6.869966e-02
# 7           csmok 0.9970172 0.023723813 -0.1259189 1258.2741 8.998162e-01
# 8            ckd1 0.6573037 0.074604492 -5.6244483 1254.4244 2.292329e-08
# 9        priorpci 1.5443939 0.128582468  3.3801774 1258.4063 7.465888e-04
# 10    priorstroke 0.9268030 0.203198921 -0.3740880 1257.3937 7.084019e-01
# 11 race.modothers 0.5288781 0.113699473 -5.6024645 1259.9413 2.592620e-08
# 12  race.modwhite 0.6248000 0.082779286 -5.6816599 1257.5458 1.655886e-08
# 13         lmcad1 0.9879270 0.069028354 -0.1759628  445.4337 8.604031e-01
# 14            sex 0.6975245 0.303962444 -1.1850730 1261.2358 2.362118e-01
# 15       diabetes 1.3310644 0.033479078  8.5420180 1253.2106 0.000000e+00
# 16      cohort_n1 4.1387159 0.072359576 19.6295450 1032.0487 0.000000e+00
# 17      cohort_n2 7.3510340 0.082587516 24.1542677  960.9672 0.000000e+00
# 18           pvd1 1.3504338 0.059750002  5.0280479 1259.0144 5.671328e-07
# 2.5 %    97.5 %
#   1  0.7543139 0.9296866
# 2  0.6458202 0.8266240
# 3  1.0126547 1.0285003
# 4  1.1341733 1.4264396
# 5  1.1644111 1.4809273
# 6  0.9918212 1.2487157
# 7  0.9516767 1.0445178
# 8  0.5678077 0.7609058
# 9  1.2000604 1.9875271
# 10 0.6220962 1.3807572
# 11 0.4231374 0.6610432
# 12 0.5311427 0.7349719
# 13 0.8625966 1.1314673
# 14 0.3842179 1.2663136
# 15 1.2464477 1.4214253
# 16 3.5908785 4.7701333
# 17 6.2511678 8.6444169
# 18 1.2010610 1.5183837
# 


