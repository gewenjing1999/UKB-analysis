#Table 2

#多肿瘤家族史与全肿瘤发生风险   四分类0，1，2，≥3


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


#多肿瘤家族史   FH_total_4

table(male$FH_total_4)
#0      1      2      3 
#146524  47087   8026   1164 

table(female$FH_total_4)
#0      1      2      3 
#167534  59088  11281   1695 


summary(male$survival_time)
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#0.002738 10.045175 10.899384 10.441894 11.646817 13.946612 

summary(female$survival_time)
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#0.002738 10.209446 10.962355 10.642759 11.682409 13.965777 



#男性多肿瘤家族史  男性全肿瘤发生风险

#case/人年#
fit2 <- pyears(Surv(survival_time,cancer_total_20)~1,data=male[male$FH_total_4==0,],scale = 1)
fit2$pyears#1534418
fit2$event#13255   0.008638454

fit2 <- pyears(Surv(survival_time,cancer_total_20)~1,data=male[male$FH_total_4==1,],scale = 1)
fit2$pyears#488890.4
fit2$event#5220    0.01067724

fit2 <- pyears(Surv(survival_time,cancer_total_20)~1,data=male[male$FH_total_4==2,],scale = 1)
fit2$pyears#82430.32
fit2$event#1033    0.0125318

fit2 <- pyears(Surv(survival_time,cancer_total_20)~1,data=male[male$FH_total_4==3,],scale = 1)
fit2$pyears#11887.58
fit2$event#160      0.01345943

#HR
model11 <- coxph(Surv(survival_time,cancer_total_20)~as.factor(FH_total_4)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)
summary(model11)
                                          #coef  exp(coef)   se(coef)      z Pr(>|z|)    
#as.factor(FH_total_4)1               0.1211888  1.1288380  0.0163834   7.397 1.39e-13 ***
#as.factor(FH_total_4)2               0.2144536  1.2391846  0.0323675   6.626 3.46e-11 ***
#as.factor(FH_total_4)3               0.2577050  1.2939571  0.0795965   3.238 0.001205 **  

#                                     exp(coef) exp(-coef) lower .95 upper .95
#as.factor(FH_total_4)1                 1.1288     0.8859    1.0932    1.1657
#as.factor(FH_total_4)2                 1.2392     0.8070    1.1630    1.3203
#as.factor(FH_total_4)3                 1.2940     0.7728    1.1070    1.5124

#p for trend
model12 <- coxph(Surv(survival_time,cancer_total_20)~as.numeric(FH_total_4)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)
summary(model12,digit=20)
#                                          coef  exp(coef)   se(coef)      z Pr(>|z|)    
#as.numeric(FH_total_4)               0.1093078  1.1155057  0.0113623   9.620  < 2e-16 ***




#女性多肿瘤家族史  女性全肿瘤发生风险  

#case/人年#
fit2 <- pyears(Surv(survival_time,cancer_total_20)~1,data=female[female$FH_total_4==0,],scale = 1)
fit2$pyears#1786561
fit2$event#11737    0.006569605

fit2 <- pyears(Surv(survival_time,cancer_total_20)~1,data=female[female$FH_total_4==1,],scale = 1)
fit2$pyears#626908.5
fit2$event#4975   0.007935767

fit2 <- pyears(Surv(survival_time,cancer_total_20)~1,data=female[female$FH_total_4==2,],scale = 1)
fit2$pyears#118877.2
fit2$event#1114       0.009371015

fit2 <- pyears(Surv(survival_time,cancer_total_20)~1,data=female[female$FH_total_4==3,],scale = 1)
fit2$pyears#17637.15
fit2$event#197     0.01116961

#HR
model11 <- coxph(Surv(survival_time,cancer_total_20)~as.factor(FH_total_4)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)
summary(model11)
#                                          coef  exp(coef)   se(coef)      z Pr(>|z|)    
#as.factor(FH_total_4)1               1.292e-01  1.138e+00  1.699e-02   7.603 2.90e-14 ***
#as.factor(FH_total_4)2               2.524e-01  1.287e+00  3.146e-02   8.022 1.04e-15 ***
#as.factor(FH_total_4)3               3.976e-01  1.488e+00  7.194e-02   5.527 3.25e-08 ***

#                                     exp(coef) exp(-coef) lower .95 upper .95
#as.factor(FH_total_4)1                 1.1379     0.8788    1.1006    1.1764
#as.factor(FH_total_4)2                 1.2870     0.7770    1.2101    1.3689
#as.factor(FH_total_4)3                 1.4883     0.6719    1.2926    1.7137

#p for trend
model12 <- coxph(Surv(survival_time,cancer_total_20)~as.numeric(FH_total_4)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)
summary(model12,digit=20)
#                                           coef  exp(coef)   se(coef)      z Pr(>|z|)    
#as.numeric(FH_total_4)               1.285e-01  1.137e+00  1.138e-02  11.287  < 2e-16 ***




