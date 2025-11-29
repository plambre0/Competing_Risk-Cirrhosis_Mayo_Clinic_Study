library(VIM)
library(naniar)
library(ggplot2)
library(mice)
library(survival)

set.seed(15634)

cirrhosis <- read.csv("C:/Users/Paolo/Desktop/Personal Research/Cirrhosis/cirrhosis.csv")

#--data coding
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
cirrhosis$ID <- NULL
cirrhosis$Transplant <- ifelse(cirrhosis$Status == "CL", 0, 1)
cirrhosis$Status <- ifelse(cirrhosis$Status == "D", 1, 0)
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
VIM::scattmatrixMiss(data.frame(cirrhosis$Stage, cirrhosis$Age))
VIM::scattmatrixMiss(data.frame(cirrhosis$Stage, cirrhosis$Sex))
VIM::scattmatrixMiss(data.frame(cirrhosis$Stage, cirrhosis$Transplant))
VIM::scattmatrixMiss(data.frame(cirrhosis$N_Days, cirrhosis$Status))
mice::fluxplot(cirrhosis)

cirrhosis_mice <- mice::mice(cirrhosis, m = 5, maxit = 50)
cirrhosis_mice_km <- with(cirrhosis_mice, coxph(Surv(N_Days, Status) ~ Drug + Age + Sex + Stage + Ascites + Hepatomegaly + Spiders + Edema + Bilirubin + Cholesterol + Albumin + Copper + Alk_Phos + SGOT + Tryglicerides + Platelets + Prothrombin))
cirrhosis_mice_km_pooled <- pool(cirrhosis_mice_km)
summary(cirrhosis_mice_km_pooled)

lapply(cirrhosis_mice_km, survfit)