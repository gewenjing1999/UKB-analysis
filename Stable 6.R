#Stable 6
#CPRS与全肿瘤发生风险

#加载数据
load("imputed_final_20240120_UKB_database_for_analysis_202801_without_prevalent_cancer_male.rdata")
load("imputed_final_20240120_UKB_database_for_analysis_239598_without_prevalent_cancer_female.rdata")


#加载R包
library("survival")
library("survminer")
library("rstpm2")
library("flexsurv")
library("rstpm2")
library("rms")


#male  CPRS #

Q20=quantile(male$CPRS,seq(0.05,1,0.05))[[4]]
Q80=quantile(male$CPRS,seq(0.05,1,0.05))[[16]]

####对原始数据的PRS三分位，<=Q20，1；>=Q80，3
male$cprs_3_Q20=0
male[male$CPRS<=Q20,]$cprs_3_Q20=1
male[male$CPRS>Q20&male$CPRS<Q80,]$cprs_3_Q20=2
male[male$CPRS>=Q80,]$cprs_3_Q20=3
table(male$cprs_3_Q20)
#1      2      3 
#40561 121679  40561 


#人年
fit2 <- pyears(Surv(survival_time,cancer_total_20)~1,data=male[male$cprs_3_Q20==1,],scale = 1)
fit2$pyears#
fit2$event#
#2647/429265

fit2 <- pyears(Surv(survival_time,cancer_total_20)~1,data=male[male$cprs_3_Q20==2,],scale = 1)
fit2$pyears#
fit2$event#
#11216/1273673

fit2 <- pyears(Surv(survival_time,cancer_total_20)~1,data=male[male$cprs_3_Q20==3,],scale = 1)
fit2$pyears#
fit2$event#
#5805/414689

#HR
model11 <- coxph(Surv(survival_time,cancer_total_20)~as.factor(cprs_3_Q20)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)
res <- summary(model11)
print(res$coefficients,digit=2)
print(res$conf.int,digit=2)
#coef exp(coef) se(coef)      z Pr(>|z|)
#as.factor(cprs_3_Q20)2               0.35556      1.43  0.02168  16.40  2.0e-60
#as.factor(cprs_3_Q20)3               0.84078      2.32  0.02358  35.65 2.0e-278

#                                        exp(coef) exp(-coef) lower .95 upper .95
#as.factor(cprs_3_Q20)2                   1.43       0.70      1.37      1.49
#as.factor(cprs_3_Q20)3                   2.32       0.43      2.21      2.43

#Concordance= 0.7  (se = 0.002 )
#Likelihood ratio test= 9813  on 20 df,   p=<2e-16
#Wald test            = 8529  on 20 df,   p=<2e-16
#Score (logrank) test = 9151  on 20 df,   p=<2e-16

#p for trend
model12 <- coxph(Surv(survival_time,cancer_total_20)~as.numeric(cprs_3_Q20)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)
res=summary(model12)
print(res$coefficients,digit=2)
print(res$conf.int,digit=2)
#                                        coef exp(coef) se(coef)       z Pr(>|z|)
#as.numeric(cprs_3_Q20)               0.43569      1.55  0.01148  37.954 3.3e-315
#                                       exp(coef) exp(-coef) lower .95 upper .95
#as.numeric(cprs_3_Q20)                   1.55       0.65      1.51      1.58




#female  CPRS #
Q20=quantile(female$CPRS,seq(0.05,1,0.05))[[4]]
Q80=quantile(female$CPRS,seq(0.05,1,0.05))[[16]]

####对原始数据的PRS三分位，<=Q20，1；>=Q80，3
female$cprs_3_Q20=0
female[female$CPRS<=Q20,]$cprs_3_Q20=1
female[female$CPRS>Q20&female$CPRS<Q80,]$cprs_3_Q20=2
female[female$CPRS>=Q80,]$cprs_3_Q20=3
table(female$cprs_3_Q20)
#1      2      3 
#47920 143758  47920


#人年
fit2 <- pyears(Surv(survival_time,cancer_total_20)~1,data=female[female$cprs_3_Q20==1,],scale = 1)
fit2$pyears# 
fit2$event#
#2825/511488

fit2 <- pyears(Surv(survival_time,cancer_total_20)~1,data=female[female$cprs_3_Q20==2,],scale = 1)
fit2$pyears#
fit2$event#
#10652/1532817

fit2 <- pyears(Surv(survival_time,cancer_total_20)~1,data=female[female$cprs_3_Q20==3,],scale = 1)
fit2$pyears#
fit2$event#
#4546/505678


#HR
model11 <- coxph(Surv(survival_time,cancer_total_20)~as.factor(cprs_3_Q20)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)
summary(model11)
res <- summary(model11)
print(res$coefficients,digit=2)
print(res$conf.int,digit=2)
#                                        coef exp(coef) se(coef)        z Pr(>|z|)
#as.factor(cprs_3_Q20)2               2.1e-01      1.23  0.02195   9.5011  2.1e-21
#as.factor(cprs_3_Q20)3               4.7e-01      1.60  0.02484  18.8540  2.7e-79

#                                      exp(coef) exp(-coef) lower .95 upper .95
#as.factor(cprs_3_Q20)2                   1.23       0.81      1.18      1.29
#as.factor(cprs_3_Q20)3                   1.60       0.63      1.52      1.68

#Concordance= 0.611  (se = 0.002 )
#Likelihood ratio test= 2776  on 20 df,   p=<2e-16
#Wald test            = 2676  on 20 df,   p=<2e-16
#Score (logrank) test = 2719  on 20 df,   p=<2e-16

#p for trend
model12 <- coxph(Surv(survival_time,cancer_total_20)~as.numeric(cprs_3_Q20)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)
res=summary(model12)
print(res$coefficients,digit=2)
print(res$conf.int,digit=2)
#                                        coef exp(coef) se(coef)        z Pr(>|z|)
#as.numeric(cprs_3_Q20)               2.4e-01      1.27  0.01222  19.5309  6.0e-85
#                                     exp(coef) exp(-coef) lower .95 upper .95
#as.numeric(cprs_3_Q20)                   1.27       0.79      1.24      1.30



