#Table 3
#家族史与PRS的关联 glm(FH)~PRS

#加载数据
load("imputed_final_20240120_UKB_database_for_analysis_202801_without_prevalent_cancer_male.rdata")
load("imputed_final_20240120_UKB_database_for_analysis_239598_without_prevalent_cancer_female.rdata")



library("survival")
library("survminer")
library("rstpm2")
library("flexsurv")
library("rstpm2")
library("rms")
library("epiDisplay")
library("performance")

##########################结直肠癌 COL_PRS##########

##male
control <- male[which(male$FH_Bowel_2==0),]
summary(control$COL_PRS)
sd(control$COL_PRS)  ##0.390638

male$COL_PRS_sd<-male$COL_PRS/0.390638

####单个肿瘤家族史与对应的PRS评分的关联####
model1=glm(FH_Bowel_2~COL_PRS_sd+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male,family=binomial())
logistic.display(model1)
summary(model1)$coef
#####OR 1.14 (1.12,1.15) P=4.353657e-73



##female#
control <- female[which(female$FH_Bowel_2==0),]
summary(control$COL_PRS)
sd(control$COL_PRS)  ##0.3913549

female$COL_PRS_sd<-female$COL_PRS/0.3913549
model1=glm(FH_Bowel_2~COL_PRS_sd+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female,family=binomial())
logistic.display(model1)
summary(model1)$coef
#####OR 1.16 (1.14,1.17) 3.150457e-112





##########################肺癌 LC_PRS##########

##male#
control <- male[which(male$FH_Lung_2==0),]
summary(control$LC_PRS)
sd(control$LC_PRS)  ##0.1093194

male$LC_PRS_sd<-male$LC_PRS/0.1093194
####单个肿瘤家族史与对应的PRS评分的关联####
model1=glm(FH_Lung_2~LC_PRS_sd+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male,family=binomial())
logistic.display(model1)
summary(model1)$coef
#####OR 1.09 (1.07,1.1) P=4.798650e-37


##female#
control <- female[which(female$FH_Lung_2==0),]
summary(control$LC_PRS)
sd(control$LC_PRS)  ##0.109029

female$LC_PRS_sd<-female$LC_PRS/0.109029

model1=glm(FH_Lung_2~LC_PRS_sd+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female,family=binomial())
logistic.display(model1)
summary(model1)$coef
#####OR 1.08 (1.06,1.09) P=4.672878e-34





##########################前列腺癌 PRO_PRS##########

##male#
control <- male[which(male$FH_Prostate_2==0),]
summary(control$PRO_PRS)
sd(control$PRO_PRS)  ##0.6953242

male$PRO_PRS_sd<-male$PRO_PRS/0.6953242

####单个肿瘤家族史与对应的PRS评分的关联####
model1=glm(FH_Prostate_2~PRO_PRS_sd+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male,family=binomial())
logistic.display(model1)
summary(model1)$coef
#####OR 1.21 (1.19,1.23) P=4.822965e-119


##########################乳腺腺癌 BC_PRS##########

##female#
control <- female[which(female$FH_Breast_2==0),]
summary(control$BC_PRS)
sd(control$BC_PRS)  ##0.1271848

female$BC_PRS_sd<-female$BC_PRS/0.1271848
####单个肿瘤家族史与对应的PRS评分的关联####
model1=glm(FH_Breast_2~BC_PRS_sd+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female,family=binomial())
logistic.display(model1)
summary(model1)$coef
#####OR 1.17 (1.15,1.18) P=6.231665e-118



##########################多肿瘤家族史 CPRS##########

##male#
control <- male[which(male$FH_total_2==0),]
summary(control$CPRS)
sd(control$CPRS)  ##0.1988989

male$CPRS_sd<-male$CPRS/0.1988989

####单个肿瘤家族史与对应的PRS评分的关联####
model1=glm(FH_total_2~CPRS_sd+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male,family=binomial())
logistic.display(model1)
summary(model1)$coef
#####OR 1.09 (1.08,1.1)  P=3.138729e-61


##female#
control <- female[which(female$FH_total_2==0),]
summary(control$CPRS)
sd(control$CPRS)  ##0.0782589

female$CPRS_sd<-female$CPRS/0.0782589
####单个肿瘤家族史与对应的PRS评分的关联####
model1=glm(FH_total_2~CPRS_sd+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female,family=binomial())
logistic.display(model1)
summary(model1)$coef
#####OR 1.09 (1.08,1.1) P=1.574716e-68


