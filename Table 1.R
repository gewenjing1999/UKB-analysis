#Table 1
##单个肿瘤家族史与对应肿瘤发生风险   家族史三分类 0，1，≥2


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

#肠癌家族史三分类 FH_Bowel_3

#肠癌家族史与肠癌发生风险的关联

#male#
#人年
fit2 <- pyears(Surv(CRC_difftime_new,CRC_total.y)~1,data=male[male$FH_Bowel_3==0,],scale = 1)
fit2$pyears#1939718
fit2$event#2384

fit2 <- pyears(Surv(CRC_difftime_new,CRC_total.y)~1,data=male[male$FH_Bowel_3==1,],scale = 1)
fit2$pyears#229699.3
fit2$event#421

fit2 <- pyears(Surv(CRC_difftime_new,CRC_total.y)~1,data=male[male$FH_Bowel_3==2,],scale = 1)
fit2$pyears#14669.24
fit2$event#33

#HR
model11 <- coxph(Surv(CRC_difftime_new,CRC_total.y)~as.factor(FH_Bowel_3)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)
summary(model11)
#                                         coef  exp(coef)   se(coef)      z Pr(>|z|)    
#as.factor(FH_Bowel_3)1               0.2843083  1.3288425  0.0529574  5.369 7.93e-08 ***
#as.factor(FH_Bowel_3)2               0.4038155  1.4975277  0.1753739  2.303 0.021302 *    
#                                      exp(coef) exp(-coef) lower .95 upper .95
#as.factor(FH_Bowel_3)1                 1.3288     0.7525    1.1978    1.4742
#as.factor(FH_Bowel_3)2                 1.4975     0.6678    1.0619    2.1118

#p for trend
model12 <- coxph(Surv(CRC_difftime_new,CRC_total.y)~as.numeric(FH_Bowel_3)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)
summary(model12)
#coef  exp(coef)   se(coef)      z Pr(>|z|)    
#as.numeric(FH_Bowel_3)               0.2619647  1.2994807  0.0458345  5.715 1.09e-08 ***





#female#
#人年
fit2 <- pyears(Surv(CRC_difftime_new,CRC_total.y)~1,data=female[female$FH_Bowel_3==0,],scale = 1)
fit2$pyears#2339690
fit2$event#1853

fit2 <- pyears(Surv(CRC_difftime_new,CRC_total.y)~1,data=female[female$FH_Bowel_3==1,],scale = 1)
fit2$pyears#267055.9
fit2$event#290

fit2 <- pyears(Surv(CRC_difftime_new,CRC_total.y)~1,data=female[female$FH_Bowel_3==2,],scale = 1)
fit2$pyears#14943.2
fit2$event#24

#HR
model11 <- coxph(Surv(CRC_difftime_new,CRC_total.y)~as.factor(FH_Bowel_3)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)
summary(model11)
#                                          coef  exp(coef)   se(coef)      z Pr(>|z|)    
#as.factor(FH_Bowel_3)1               2.178e-01  1.243e+00  6.328e-02  3.441 0.000579 ***
#as.factor(FH_Bowel_3)2               5.308e-01  1.700e+00  2.056e-01  2.582 0.009835 **
#                                      exp(coef) exp(-coef) lower .95 upper .95
#as.factor(FH_Bowel_3)1                 1.2433     0.8043    1.0983    1.4075
#as.factor(FH_Bowel_3)2                 1.7002     0.5882    1.1363    2.5440

#p for trend
model12 <- coxph(Surv(CRC_difftime_new,CRC_total.y)~as.numeric(FH_Bowel_3)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)
summary(model12)
#coef  exp(coef)   se(coef)      z Pr(>|z|)    
#as.numeric(FH_Bowel_3)               2.300e-01  1.259e+00  5.520e-02  4.168 3.08e-05 ***





#肺癌单个家族史三分类   FH_Lung_3

#肠癌家族史  肠癌发生风险

#male#
#人年
fit2 <- pyears(Surv(LC_difftime_new,LC_total.y)~1,data=male[male$FH_Lung_3==0,],scale = 1)
fit2$pyears#1926980
fit2$event#1426

fit2 <- pyears(Surv(LC_difftime_new,LC_total.y)~1,data=male[male$FH_Lung_3==1,],scale = 1)
fit2$pyears#247071.1
fit2$event#344

fit2 <- pyears(Surv(LC_difftime_new,LC_total.y)~1,data=male[male$FH_Lung_3==2,],scale = 1)
fit2$pyears#18710.49
fit2$event#37

#HR
model11 <- coxph(Surv(LC_difftime_new,LC_total.y)~as.factor(FH_Lung_3)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)
summary(model11)
#                                           coef  exp(coef)   se(coef)      z Pr(>|z|)    
#as.factor(FH_Lung_3)1                0.4596130  1.5834610  0.0603385   7.617 2.59e-14 ***
#as.factor(FH_Lung_3)2                0.6277984  1.8734815  0.1668663   3.762 0.000168 *** 
#                                      exp(coef) exp(-coef) lower .95 upper .95
#as.factor(FH_Lung_3)1                  1.5835     0.6315    1.4068    1.7822
#as.factor(FH_Lung_3)2                  1.8735     0.5338    1.3509    2.5983

#p for trend
model12 <- coxph(Surv(LC_difftime_new,LC_total.y)~as.numeric(FH_Lung_3)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)
summary(model12)
#                                           coef  exp(coef)   se(coef)      z Pr(>|z|)    
#as.numeric(FH_Lung_3)                0.4073501  1.5028302  0.0496243   8.209 2.24e-16 ***


#female#
#人年
fit2 <- pyears(Surv(LC_difftime_new,LC_total.y)~1,data=female[female$FH_Lung_3==0,],scale = 1)
fit2$pyears#2294471
fit2$event#1300

fit2 <- pyears(Surv(LC_difftime_new,LC_total.y)~1,data=female[female$FH_Lung_3==1,],scale = 1)
fit2$pyears#311702.6
fit2$event#332

fit2 <- pyears(Surv(LC_difftime_new,LC_total.y)~1,data=female[female$FH_Lung_3==2,],scale = 1)
fit2$pyears#21205.51
fit2$event#50

#HR
model11 <- coxph(Surv(LC_difftime_new,LC_total.y)~as.factor(FH_Lung_3)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)
summary(model11)
#                                           coef  exp(coef)   se(coef)      z Pr(>|z|)    
#as.factor(FH_Lung_3)1                0.4238939  1.5278995  0.0617590   6.864 6.71e-12 ***
#as.factor(FH_Lung_3)2                1.0449494  2.8432545  0.1446585   7.224 5.06e-13 ***
#                                     exp(coef) exp(-coef) lower .95 upper .95
#as.factor(FH_Lung_3)1                  1.5279     0.6545    1.3537    1.7245
#as.factor(FH_Lung_3)2                  2.8433     0.3517    2.1413    3.7753

#p for trend
model12 <- coxph(Surv(LC_difftime_new,LC_total.y)~as.numeric(FH_Lung_3)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)
summary(model12)
#                                         coef  exp(coef)   se(coef)      z Pr(>|z|)    
#as.numeric(FH_Lung_3)                0.4628140  1.5885378  0.0497624   9.300  < 2e-16 ***




#前列腺癌单个家族史三分类  FH_Prostate_3

#male#

#人年
fit2 <- pyears(Surv(PRC_difftime_new,PRC_total.y)~1,data=male[male$FH_Prostate_3==0,],scale = 1)
fit2$pyears#1985878
fit2$event#8347

fit2 <- pyears(Surv(PRC_difftime_new,PRC_total.y)~1,data=male[male$FH_Prostate_3==1,],scale = 1)
fit2$pyears#160896.7
fit2$event# 1172

fit2 <- pyears(Surv(PRC_difftime_new,PRC_total.y)~1,data=male[male$FH_Prostate_3==2,],scale = 1)
fit2$pyears#3855.556
fit2$event#66

#HR
model11 <- coxph(Surv(PRC_difftime_new,PRC_total.y)~as.factor(FH_Prostate_3)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)
summary(model11)
#                                          coef  exp(coef)   se(coef)      z Pr(>|z|)    
#as.factor(FH_Prostate_3)1            5.264e-01  1.693e+00  3.122e-02 16.857  < 2e-16 ***
#as.factor(FH_Prostate_3)2            1.099e+00  3.000e+00  1.237e-01  8.885  < 2e-16 ***
#                                     exp(coef) exp(-coef) lower .95 upper .95
#as.factor(FH_Prostate_3)1              1.6927     0.5908    1.5923    1.7996
#as.factor(FH_Prostate_3)2              3.0003     0.3333    2.3546    3.8232

#p for trend
model12 <- coxph(Surv(PRC_difftime_new,PRC_total.y)~as.numeric(FH_Prostate_3)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)
summary(model12)
#                                          coef  exp(coef)   se(coef)      z Pr(>|z|)    
#as.numeric(FH_Prostate_3)            5.308e-01  1.700e+00  2.827e-02 18.775  < 2e-16 ***



#乳腺癌单个家族史三分类   FH_Breast_3

#female#

#人年
fit2 <- pyears(Surv(BC_difftime_new,BC_total)~1,data=female[female$FH_Breast_3==0,],scale = 1)
fit2$pyears#2306354
fit2$event#6867

fit2 <- pyears(Surv(BC_difftime_new,BC_total)~1,data=female[female$FH_Breast_3==1,],scale = 1)
fit2$pyears#268292.9
fit2$event#1232

fit2 <- pyears(Surv(BC_difftime_new,BC_total)~1,data=female[female$FH_Breast_3==2,],scale = 1)
fit2$pyears#12243.26
fit2$event#97

#HR
model11 <- coxph(Surv(BC_difftime_new,BC_total)~as.factor(FH_Breast_3)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)
summary(model11)
#                                          coef  exp(coef)   se(coef)      z Pr(>|z|)    
#as.factor(FH_Breast_3)1              0.4104795  1.5075405  0.0309762 13.251  < 2e-16 ***
#as.factor(FH_Breast_3)2              0.9254082  2.5228980  0.1023180  9.044  < 2e-16 ***
#                                       exp(coef) exp(-coef) lower .95 upper .95
#as.factor(FH_Breast_3)1                1.5075     0.6633    1.4187    1.6019
#as.factor(FH_Breast_3)2                2.5229     0.3964    2.0645    3.0831

#p for trend
model12 <- coxph(Surv(BC_difftime_new,BC_total)~as.numeric(FH_Breast_3)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)
summary(model12,digit=10)

#coef  exp(coef)   se(coef)      z Pr(>|z|)    
#as.numeric(FH_Breast_3)              0.4235111  1.5273147  0.0271871 15.578  < 2e-16 ***
