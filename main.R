library(dplyr)
library(VIM)
library(naniar)
library(ggplot2)
library(umap)
library(mice)
library(survival)
library(cmprsk)
library(gtsummary)

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
summary(cirrhosis)
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
#umap visualization
cirrhosis_exc_umap <- 
  umap::umap(cirrhosis_exc[,c("Bilirubin", "Cholesterol", "Albumin",
                              "Copper", "Alk_Phos", "SGOT", "Tryglicerides",
                              "Platelets", "Prothrombin", "Age")])
plot(cirrhosis_exc_umap$layout, col = as.factor(cirrhosis_exc$Stage))
plot(cirrhosis_exc_umap$layout, col = as.factor(cirrhosis_exc$Status))
plot(cirrhosis_exc_umap$layout, col = cirrhosis_exc$N_Days)
plot(survfit(Surv(N_Days, Status==1) ~ Stage, data = cirrhosis_exc))
plot(survfit(Surv(N_Days, Status==1) ~ Sex, data = cirrhosis_exc))
plot(survfit(Surv(N_Days, Status==1) ~ Drug, data = cirrhosis_exc))
plot(survfit(Surv(N_Days, Status==2) ~ Bilirubin > 1.2, data = cirrhosis_exc))
plot(survfit(Surv(N_Days, Status==2) ~ Cholesterol > 200, data = cirrhosis_exc))
plot(survfit(Surv(N_Days, Status==2) ~ Albumin < 3.5, data = cirrhosis_exc))

#distribution plots
plot(density(cirrhosis_exc$Bilirubin))
plot(density(cirrhosis_exc$Cholesterol))
plot(density(cirrhosis_exc$Albumin))
plot(density(cirrhosis_exc$Copper))
plot(density(cirrhosis_exc$Alk_Phos))
plot(density(cirrhosis_exc$SGOT))
plot(density(cirrhosis_exc$Tryglicerides))
plot(density(cirrhosis_exc$Platelets))
plot(density(cirrhosis_exc$Prothrombin))
#transformed distribution plots
plot(density(log(cirrhosis_exc$Bilirubin)))
plot(density(log(cirrhosis_exc$Cholesterol)))
plot(density(log(cirrhosis_exc$Albumin)))
plot(density(log(cirrhosis_exc$Copper)))
plot(density(log(cirrhosis_exc$Alk_Phos)))
plot(density(log(cirrhosis_exc$SGOT)))
plot(density(log(cirrhosis_exc$Tryglicerides)))
plot(density(log(cirrhosis_exc$Platelets)))
plot(density(log(cirrhosis_exc$Prothrombin)))

#transform vars
cirrhosis_transf <- cirrhosis
cirrhosis_transf[,c("Bilirubin", "Cholesterol", "Albumin",
                        "Copper", "Alk_Phos", "SGOT", "Tryglicerides",
                        "Platelets", "Prothrombin")] <-
  lapply(cirrhosis_transf[,c("Bilirubin", "Cholesterol", "Albumin",
                                 "Copper", "Alk_Phos", "SGOT", "Tryglicerides",
                                 "Platelets", "Prothrombin")], log)
  

#correlation analyses
cirrhosis_exc_cor <- polycor::hetcor(cirrhosis_transf[, !names(cirrhosis_transf) %in% "ID"])
heatmap(cirrhosis_exc_cor$correlations)

#--missingness analysis
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

cirrhosis_quickpred <- mice::quickpred(cirrhosis_transf, minpuc = .3, mincor = .3)
cirrhosis_mice <- mice::mice(cirrhosis_transf, m = 5, maxit = 50, predictorMatrix = cirrhosis_quickpred)
plot(cirrhosis_mice)

cirrhosis_mice_cox_td <- with(
  cirrhosis_mice,
  coxph(
    Surv(N_Days, Status == 1) ~
      Drug + Age + Sex + Stage + Ascites + Hepatomegaly + Spiders +
      Edema + Bilirubin + Cholesterol + Albumin + Copper + Alk_Phos +
      SGOT + Tryglicerides + Platelets + Prothrombin +
      
      # time interactions
      tt(Age) + tt(Bilirubin) + tt(Cholesterol) + tt(Albumin) + tt(Copper) +
      tt(Alk_Phos) + tt(SGOT) + tt(Tryglicerides) + tt(Platelets) + 
      tt(Prothrombin),
    
    tt = function(x, t, ...) x * log(t)
  )
)
ph_pooled <- pool(cirrhosis_mice_cox_td)
summary(ph_pooled)

#multi-state cox
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
summary(cirrhosis_mice_cox_pooled)

#fine-grey model
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
