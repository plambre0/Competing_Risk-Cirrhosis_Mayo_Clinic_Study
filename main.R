library(dplyr)
library(VIM)
library(naniar)
library(ggplot2)
library(mice)
library(survival)
library(cmprsk)

set.seed(15634)

cirrhosis <- read.csv("C:/Users/Paolo/Desktop/Personal Research/Cirrhosis/cirrhosis.csv")

#--data coding
cirrhosis$Status[cirrhosis$Status == "CL"] <- 2
cirrhosis$Status[cirrhosis$Status == "D"] <- 1
cirrhosis$Status[cirrhosis$Status == "C"] <- 0
#convert to factors
cirrhosis[, c("ID", 
               "Status", 
               "Drug", 
               "Sex", 
               "Ascites", 
               "Hepatomegaly", 
               "Spiders", 
               "Edema",
               "Stage")] <- lapply(cirrhosis[, c("ID", 
                                                 "Status", 
                                                 "Drug", 
                                                 "Sex", 
                                                 "Ascites", 
                                                 "Hepatomegaly", 
                                                 "Spiders", 
                                                 "Edema",
                                                 "Stage")], as.factor)
#code vars for censoring and indicator for transplant 
cirrhosis$Stage <- factor(cirrhosis$Stage, levels = c(1, 2, 3, 4))

#--exploratory analysis
#for data exploration, using complete case analysis
cirrhosis_exc <- na.exclude(cirrhosis)
#plot comp case curves
plot(survfit(data = cirrhosis_exc, 
             Surv(N_Days, Status) ~ 1), 
     xlab = "Days", ylab = "Survival Rate", 
     main = "Survival in Cirrhosis Patients by Day (Complete Cases)", 
     cex.main = .95)
plot(survfit(data = cirrhosis_exc[cirrhosis_exc$Drug == "Placebo",], 
             Surv(N_Days, Status) ~ 1), 
     xlab = "Days", ylab = "Survival Rate", 
     main = "Survival in Cirrhosis Patients by Day on Placebo (Complete Cases)", 
     cex.main = .95)
plot(survfit(data = cirrhosis_exc[cirrhosis_exc$Drug == "D-penicillamine",], 
             Surv(N_Days, Status) ~ 1), 
     xlab = "Days", ylab = "Survival Rate", 
     main = "Survival in Cirrhosis Patients by Day on D-Penicillamine (Complete Cases)", 
     cex.main = .95)
#correlation analyses
cirrhosis_exc_cor <- polycor::hetcor(cirrhosis_exc)

#--missingness analysis
summary(cirrhosis)
naniar::miss_var_summary(cirrhosis)
naniar::vis_miss(cirrhosis) + theme(axis.text.x =  element_text(angle = 90))
naniar::gg_miss_upset(cirrhosis, nsets = 11)
naniar::mcar_test(cirrhosis)
VIM::parcoordMiss(cirrhosis)
VIM::marginplot(data.frame(cirrhosis$Stage, cirrhosis$Age))
VIM::marginplot(data.frame(cirrhosis$Stage, cirrhosis$Sex))
VIM::marginplot(data.frame(cirrhosis$Stage, cirrhosis$Status))
VIM::marginplot(data.frame(cirrhosis$N_Days, cirrhosis$Age))
VIM::marginplot(data.frame(cirrhosis$N_Days, cirrhosis$Sex))
VIM::marginplot(data.frame(cirrhosis$N_Days, cirrhosis$Status))
mice::fluxplot(cirrhosis)

cirrhosis_quickpred <- mice::quickpred(cirrhosis, minpuc = .3, mincor = .3)
cirrhosis_mice <- mice::mice(cirrhosis, m = 5, maxit = 50, predictorMatrix = cirrhosis_quickpred)

cirrhosis_mice_cox <- with(
  cirrhosis_mice,
  coxph(
    Surv(N_Days, Status) ~ 
      Drug + Age + Sex + Stage + Ascites + Hepatomegaly + Spiders +
      Edema + Bilirubin + Cholesterol + Albumin + Copper + Alk_Phos +
      SGOT + Tryglicerides + Platelets + Prothrombin +
      cluster(ID),
    id = ID
  )
)
cirrhosis_mice_cox_pooled <- pool(cirrhosis_mice_cox)
summary(cirrhosis_mice_cox_pooled)

cirrhosis_mice_fg <- with(
  cirrhosis_mice,
  crr(ftime = N_Days,
      fstatus = Status,
      cov1 = model.matrix(~ Drug + Age + Sex + Stage + Ascites + Hepatomegaly +
                            Spiders + Edema + Bilirubin + Cholesterol + Albumin +
                            Copper + Alk_Phos + SGOT + Tryglicerides + Platelets +
                            Prothrombin)[,-1]))
fg_pooled <- pool(cirrhosis_mice_fg)
summary(fg_pooled)
