#Stable 2
#家族史与肿瘤发生风险校正PRS

#加载数据
load("imputed_final_20240120_UKB_database_for_analysis_202801_without_prevalent_cancer_male.rdata")
load("imputed_final_20240120_UKB_database_for_analysis_239598_without_prevalent_cancer_female.rdata")

#加载R包
library("survival")
library("survminer")
library("rstpm2")
library("flexsurv")
library("rms")
library("lmtest")
library("pROC") 
library("meta")
library("magrittr")


#定义总结函数结果的代码
re=data.frame(matrix(NA,1,8))
names(re)=c('Crp_concentration','N_controls','N_cases','HR(95%CI)','Pvalue','HR(95%CI)','Pvalue','PY')

output1=function(x){  	 x <- summary(x)
HR <-round(x$coef[1,2],2);
HR.confint.lower <- round(x$conf.int[1,"lower .95"],2)
HR.confint.upper <- round(x$conf.int[1,"upper .95"],2) 
p.value<-signif(x$coefficients[1,"Pr(>|z|)"],2)
res<-c(paste0(HR,'(',HR.confint.lower,'-', HR.confint.upper,')'),p.value)
return(res)
}


#定义总结函数的结果进行meta的代码
output = function(x) {
  x <- summary(x)
  Beta = round(x$coef[1, 1],2)
  se = round(x$coef[1, 3],2)
  HR = round(x$coef[1, 2],2)
  HR.confint.lower = round(x$conf.int[1, "lower .95"], 2)
  HR.confint.upper = round(x$conf.int[1, "upper .95"], 2)
  p.value = signif(x$coefficients[1, "Pr(>|z|)"], 3)
  res <-
    c(Beta, se, HR, HR.confint.lower, HR.confint.upper, p.value)
  return(res)
}

result = data.frame(matrix(NA, 1, 6))
names(result) = c(
  'Beta',
  'se',
  'HR',
  'lower',
  'upper',
  'Pvalue'
)


########肠癌 CRC model1：家族史（1，0有无）~肿瘤发生风险   model2：家族史+PRS（SD）~肿瘤发生风险


#male#
control <- male[which(male$FH_Bowel_2==0),]
summary(control$COL_PRS)
sd(control$COL_PRS)  ##0.390638
male$COL_PRS_sd<-male$COL_PRS/0.390638

table(male$FH_Bowel_2)
#0      1 
#180048  22753 


model11 <- coxph(Surv(CRC_difftime_new,CRC_total.y)~as.factor(FH_Bowel_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)
summary(model11)
hr <- output1(model11)
hr
#"1.34(1.21-1.48)" "1.2e-08" 

male$prob=predict(model11,newdata=male)
roc2=roc(CRC_total.y~prob,data=male,type="expected")
auc(roc2)#Area under the curve: 0.6804

round(auc(roc2),3)#0.68

ci.auc(roc2,
       conf.level=0.95,
       partial.auc=c(1, .8),  partial.auc.focus="se", partial.auc.correct=TRUE,
       boot.n=10000, stratified=FALSE)
#95% CI: 0.6713-0.6895 (DeLong)



model22 <- coxph(Surv(CRC_difftime_new,CRC_total.y)~as.factor(FH_Bowel_2)+COL_PRS_sd+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health
                 +BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)
summary(model22)
hr <- output1(model22)
hr
#"1.28(1.15-1.41)" "2.1e-06"
#                                         coef  exp(coef)   se(coef)      z Pr(>|z|)    
#as.factor(FH_Bowel_2)1               0.2435461  1.2757652  0.0513575  4.742 2.11e-06 ***
#COL_PRS_sd                           0.3932755  1.4818265  0.0181758 21.637  < 2e-16 ***


male$prob=predict(model22,newdata=male)
roc3=roc(CRC_total.y~prob,data=male,type="expected")
auc(roc3)#Area under the curve: 0.7142

round(auc(roc3),3)#0.712

ci.auc(roc3,
       conf.level=0.95,
       partial.auc=c(1, .8),  partial.auc.focus="se", partial.auc.correct=TRUE,
       boot.n=10000, stratified=FALSE)
#95% CI: 0.7028-0.7203 (DeLong)

roc.test(roc3,roc2) #拿roc3-roc2

#DeLong检验是一种用于比较两个相关的ROC曲线之间差异的统计检验方法。
#该检验基于两个ROC曲线之间的协方差矩阵，通过计算Z统计量和对应的p值来评估两个曲线之间的显著性差异

#DeLong's test for two correlated ROC curves

#data:  roc3 and roc2
#Z = 11.666, p-value < 2.2e-16 
#alternative hypothesis: true difference in AUC is not equal to 0  
#95 percent confidence interval: 
# 0.02590204 0.03636280 
#sample estimates:
#AUC of roc1 AUC of roc2 
#  0.7115448   0.6804123




####异质性####
##要比较model1/2结果的异质性##
res = output(model11)
result[1, ] = c(res[1:6])
res = output(model22)
result[2, ] = c(res[1:6])
result=as.data.frame(lapply(result[,1:6],as.numeric))
metaresult <-
  metagen(
    result$Beta,
    result$se,
    lower = result$lower,
    upper = result$upper,
    sm = "β",
    comb.fixed = FALSE,
    data = result
  )

summary(metaresult)
#Test of heterogeneity:  
#  Q d.f. p-value
#0.50    1  0.4795




#female#
control <- female[which(female$FH_Bowel_2==0),]
summary(control$COL_PRS)
sd(control$COL_PRS)  ##0.3913549
female$COL_PRS_sd<-female$COL_PRS/0.3913549

table(female$FH_Bowel_2)
#0      1 
#213800  25798 


model11 <- coxph(Surv(CRC_difftime_new,CRC_total.y)~as.factor(FH_Bowel_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)
summary(model11)
hr <- output1(model11)
hr
#"1.27(1.13-1.43)" "9.7e-05"

female$prob=predict(model11,newdata=female)
roc2=roc(CRC_total.y~prob,data=female,type="expected")
auc(roc2)#Area under the curve: 0.651

round(auc(roc2),3)#0.651

ci.auc(roc2,
       conf.level=0.95,
       partial.auc=c(1, .8),  partial.auc.focus="se", partial.auc.correct=TRUE,
       boot.n=10000, stratified=FALSE)
#95% CI: 0.6399-0.6621 (DeLong)

model22 <- coxph(Surv(CRC_difftime_new,CRC_total.y)~as.factor(FH_Bowel_2)+COL_PRS_sd+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)
summary(model22)
hr <- output1(model22)
hr
#"1.2(1.06-1.35)" "0.0033"
#exp(coef) exp(-coef) lower .95 upper .95
#as.factor(FH_Bowel_2)1                 1.1971     0.8354    1.0616    1.3498
#COL_PRS_sd                             1.4648     0.6827    1.4063    1.5258


female$prob=predict(model22,newdata=female)
roc3=roc(CRC_total.y~prob,data=female,type="expected")
auc(roc3)#Area under the curve: 0.6842

round(auc(roc3),3)#0.684

ci.auc(roc3,
       conf.level=0.95,
       partial.auc=c(1, .8),  partial.auc.focus="se", partial.auc.correct=TRUE,
       boot.n=10000, stratified=FALSE)
#95% CI: 0.6759-0.6973 (DeLong)

roc.test(roc3,roc2) #拿roc3-roc2

#DeLong's test for two correlated ROC curves

#data:  roc3 and roc2
#Z = 9.5104, p-value < 2.2e-16
#alternative hypothesis: true difference in AUC is not equal to 0
#95 percent confidence interval:
# 0.02637793 0.04007239
##sample estimates:
#AUC of roc1 AUC of roc2 
# 0.6842348   0.6510096


##要比较model1/2结果的异质性##
res = output(model11)
result[1, ] = c(res[1:6])
res = output(model22)
result[2, ] = c(res[1:6])
result=as.data.frame(lapply(result[,1:6],as.numeric))
metaresult <-
  metagen(
    result$Beta,
    result$se,
    lower = result$lower,
    upper = result$upper,
    sm = "β",
    comb.fixed = FALSE,
    data = result
  )

summary(metaresult)
#Test of heterogeneity:
#  Q d.f. p-value
#0.50    1  0.4795



########肺癌 LC model1：家族史（1，0有无）~肿瘤发生风险   model2：家族史+PRS（SD）~肿瘤发生风险

#male#
control <- male[which(male$FH_Lung_2==0),]
summary(control$LC_PRS)
sd(control$LC_PRS)  ##0.1093194
male$LC_PRS_sd<-male$LC_PRS/0.1093194

table(male$FH_Lung_2)
#0      1 
#178146  24655 


model11 <- coxph(Surv(LC_difftime_new,LC_total.y)~as.factor(FH_Lung_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)
summary(model11)
hr <- output1(model11)
hr
#"1.61(1.43-1.8)" "2.7e-16"


male$prob=predict(model11,newdata=male)
roc2=roc(LC_total.y~prob,data=male,type="expected")
auc(roc2)#Area under the curve: 0.8145

round(auc(roc2),3)#0.814

ci.auc(roc2,
       conf.level=0.95,
       partial.auc=c(1, .8),  partial.auc.focus="se", partial.auc.correct=TRUE,
       boot.n=10000, stratified=FALSE)
#95% CI: 0.8059-0.8231 (DeLong)


model22 <- coxph(Surv(LC_difftime_new,LC_total.y)~as.factor(FH_Lung_2)+LC_PRS_sd+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)
summary(model22)
hr <- output1(model22)
hr
#"1.58(1.41-1.77)" "3.2e-15"

#coef  exp(coef)   se(coef)       z Pr(>|z|)    
#as.factor(FH_Lung_2)1                4.574e-01  1.580e+00  5.802e-02   7.883 3.19e-15 ***
#LC_PRS_sd                            1.966e-01  1.217e+00  2.082e-02   9.443  < 2e-16 ***

male$prob=predict(model22,newdata=male)
roc3=roc(LC_total.y~prob,data=male,type="expected")
auc(roc3)#Area under the curve: 0.8185

round(auc(roc3),3)#0.819

ci.auc(roc3,
       conf.level=0.95,
       partial.auc=c(1, .8),  partial.auc.focus="se", partial.auc.correct=TRUE,
       boot.n=10000, stratified=FALSE)
#95% CI: 0.8099-0.8271 (DeLong)

roc.test(roc3,roc2) #拿roc3-roc2
#DeLong's test for two correlated ROC curves

#data:  roc3 and roc2
#Z = 5.0625, p-value = 4.137e-07
#alternative hypothesis: true difference in AUC is not equal to 0
#95 percent confidence interval:
# 0.002472787 0.005597007
#sample estimates:
#AUC of roc1 AUC of roc2 
#0.8185167   0.8144818 




####异质性####
##要比较model1/2结果的异质性##
res = output(model11)
result[1, ] = c(res[1:6])
res = output(model22)
result[2, ] = c(res[1:6])
result=as.data.frame(lapply(result[,1:6],as.numeric))
metaresult <-
  metagen(
    result$Beta,
    result$se,
    lower = result$lower,
    upper = result$upper,
    sm = "β",
    comb.fixed = FALSE,
    data = result
  )

summary(metaresult)
#Test of heterogeneity:
#  Q d.f. p-value
#0.01    1  0.9062



#female#
control <- female[which(female$FH_Lung_2==0),]
summary(control$LC_PRS)
sd(control$LC_PRS)  ##0.109029
female$LC_PRS_sd<-female$LC_PRS/0.109029

table(female$FH_Lung_2)
#0      1 
#209249  30349 


model11 <- coxph(Surv(LC_difftime_new,LC_total.y)~as.factor(FH_Lung_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)
summary(model11)
hr <- output1(model11)
hr
#"1.63(1.45-1.82)" "1.1e-16"


female$prob=predict(model11,newdata=female)
roc2=roc(LC_total.y~prob,data=female,type="expected")
auc(roc2)#Area under the curve: 0.8039

round(auc(roc2),3)#0.804

ci.auc(roc2,
       conf.level=0.95,
       partial.auc=c(1, .8),  partial.auc.focus="se", partial.auc.correct=TRUE,
       boot.n=10000, stratified=FALSE)
#95% CI: 0.7938-0.8141 (DeLong)


model22 <- coxph(Surv(LC_difftime_new,LC_total.y)~as.factor(FH_Lung_2)+LC_PRS_sd+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)
summary(model22)
hr <- output1(model22)
hr
#""1.6(1.43-1.8)" "8.4e-16"
#coef  exp(coef)   se(coef)       z Pr(>|z|)    
#as.factor(FH_Lung_2)1                0.4712261  1.6019572  0.0585509   8.048 8.41e-16 ***
#LC_PRS_sd                            0.1991367  1.2203487  0.0217008   9.176  < 2e-16 ***


female$prob=predict(model22,newdata=female)
roc3=roc(LC_total.y~prob,data=female,type="expected")
auc(roc3)#Area under the curve: 0.8074

round(auc(roc3),3)#0.807

ci.auc(roc3,
       conf.level=0.95,
       partial.auc=c(1, .8),  partial.auc.focus="se", partial.auc.correct=TRUE,
       boot.n=10000, stratified=FALSE)
#95% CI: 0.7974-0.8175 (DeLong)


roc.test(roc3,roc2) #拿roc3-roc2

#DeLong's test for two correlated ROC curves

#data:  roc3 and roc2
#Z = 4.0087, p-value = 6.104e-05
#alternative hypothesis: true difference in AUC is not equal to 0
#95 percent confidence interval:
# 0.001787147 0.005206489
#sample estimates:
#  AUC of roc1 AUC of roc2 
#0.8074423   0.8039454



##要比较model1/2结果的异质性##
res = output(model11)
result[1, ] = c(res[1:6])
res = output(model22)
result[2, ] = c(res[1:6])
result=as.data.frame(lapply(result[,1:6],as.numeric))
metaresult <-
  metagen(
    result$Beta,
    result$se,
    lower = result$lower,
    upper = result$upper,
    sm = "β",
    comb.fixed = FALSE,
    data = result
  )

summary(metaresult)
#Test of heterogeneity:
#  Q d.f. p-value
#0.06    1  0.8137



########前列腺癌 PRC model1：家族史（1，0有无）~肿瘤发生风险   model2：家族史+PRS（SD）~肿瘤发生风险

#male#
control <- male[which(male$FH_Prostate_2==0),]
summary(control$PRO_PRS)
sd(control$PRO_PRS)  ##0.6953242
male$PRO_PRS_sd<-male$PRO_PRS/0.6953242

table(male$FH_Prostate_2)
#0      1 
#187028  15773 


model11 <- coxph(Surv(PRC_difftime_new,PRC_total.y)~as.factor(FH_Prostate_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)
summary(model11)
hr <- output1(model11)
hr
#"1.73(1.63-1.84)" "1e-72" 

male$prob=predict(model11,newdata=male)
roc2=roc(PRC_total.y~prob,data=male,type="expected")
auc(roc2)#Area under the curve: 0.7031

round(auc(roc2),3)#0.703

ci.auc(roc2,
       conf.level=0.95,
       partial.auc=c(1, .8),  partial.auc.focus="se", partial.auc.correct=TRUE,
       boot.n=10000, stratified=FALSE)
#95% CI: 0.6985-0.7078 (DeLong)



model22 <- coxph(Surv(PRC_difftime_new,PRC_total.y)~as.factor(FH_Prostate_2)+PRO_PRS_sd+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)
summary(model22)
hr <- output1(model22)
hr
#"1.57(1.47-1.66)" "9.9e-49"  

#coef  exp(coef)   se(coef)      z Pr(>|z|)    
#as.factor(FH_Prostate_2)1            0.4480488  1.5652550  0.0305398 14.671  < 2e-16 ***
#PRO_PRS_sd                           0.5985804  1.8195340  0.0100327 59.663  < 2e-16 ***

male$prob=predict(model22,newdata=male)
roc3=roc(PRC_total.y~prob,data=male,type="expected")
auc(roc3)#Area under the curve: 0.763

round(auc(roc3),3)#0.763

ci.auc(roc3,
       conf.level=0.95,
       partial.auc=c(1, .8),  partial.auc.focus="se", partial.auc.correct=TRUE,
       boot.n=10000, stratified=FALSE)
#95% CI: 0.7587-0.7674 (DeLong)

roc.test(roc3,roc2) #拿roc3-roc2
#DeLong's test for two correlated ROC curves

#data:  roc3 and roc2
#Z = 33.384, p-value < 2.2e-16
#alternative hypothesis: true difference in AUC is not equal to 0
#95 percent confidence interval:
# 0.05636767 0.06339916
#sample estimates:
#AUC of roc1 AUC of roc2 
#0.7630314   0.7031479 



####异质性####
##要比较model1/2结果的异质性##
res = output(model11)
result[1, ] = c(res[1:6])
res = output(model22)
result[2, ] = c(res[1:6])
result=as.data.frame(lapply(result[,1:6],as.numeric))
metaresult <-
  metagen(
    result$Beta,
    result$se,
    lower = result$lower,
    upper = result$upper,
    sm = "β",
    comb.fixed = FALSE,
    data = result
  )

summary(metaresult)
#Test of heterogeneity:
#  Q d.f. p-value
#5.56    1  0.0184




########乳腺癌 BRC model1：家族史（1，0有无）~肿瘤发生风险   model2：家族史+PRS（SD）~肿瘤发生风险
#female#
control <- female[which(female$FH_Breast_2==0),]
summary(control$BC_PRS)
sd(control$BC_PRS)  ##0.1271848
female$BC_PRS_sd<-female$BC_PRS/0.1271848

table(female$FH_Breast_2)
#0      1 
#213373  26225 


model11 <- coxph(Surv(BC_difftime_new,BC_total)~as.factor(FH_Breast_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)
summary(model11)
hr <- output1(model11)
hr
#"1.55(1.46-1.65)" "1e-48"

female$prob=predict(model11,newdata=female)
roc2=roc(BC_total~prob,data=female,type="expected")
auc(roc2)#Area under the curve: 0.5733

round(auc(roc2),3)#0.573

ci.auc(roc2,
       conf.level=0.95,
       partial.auc=c(1, .8),  partial.auc.focus="se", partial.auc.correct=TRUE,
       boot.n=10000, stratified=FALSE)
#95% CI: 0.567-0.5796 (DeLong)


model22 <- coxph(Surv(BC_difftime_new,BC_total)~as.factor(FH_Breast_2)+BC_PRS_sd+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)
summary(model22)
hr <- output1(model22)
hr
# "1.48(1.39-1.57)" "8.7e-39"  
#                                         coef  exp(coef)   se(coef)      z Pr(>|z|)    
#as.factor(FH_Breast_2)1              0.3914790  1.4791669  0.0300544 13.026  < 2e-16 ***
#BC_PRS_sd                            0.3256644  1.3849505  0.0111418 29.229  < 2e-16 ***

female$prob=predict(model22,newdata=female)
roc3=roc(BC_total~prob,data=female,type="expected")
auc(roc3)#Area under the curve: 0.6183

round(auc(roc3),3)#0.618

ci.auc(roc3,
       conf.level=0.95,
       partial.auc=c(1, .8),  partial.auc.focus="se", partial.auc.correct=TRUE,
       boot.n=10000, stratified=FALSE)
#95% CI: 0.6122-0.6244 (DeLong)

roc.test(roc3,roc2) #拿roc3-roc2

#DeLong's test for two correlated ROC curves

#data:  roc3 and roc2
#Z = 15.792, p-value < 2.2e-16
#alternative hypothesis: true difference in AUC is not equal to 0
#95 percent confidence interval:
# 0.03939560 0.05055984
#sample estimates:
#AUC of roc1 AUC of roc2 
#0.6182797   0.5733020
  



##要比较model1/2结果的异质性##
res = output(model11)
result[1, ] = c(res[1:6])
res = output(model22)
result[2, ] = c(res[1:6])
result=as.data.frame(lapply(result[,1:6],as.numeric))
metaresult <-
  metagen(
    result$Beta,
    result$se,
    lower = result$lower,
    upper = result$upper,
    sm = "β",
    comb.fixed = FALSE,
    data = result
  )

summary(metaresult)
#Test of heterogeneity:
#  Q d.f. p-value
#1.39    1  0.2386




#######多肿瘤  model1：家族史（1，0有无）~肿瘤发生风险   model2：家族史+PRS（SD）~肿瘤发生风险


#male#
control <- male[which(male$FH_total_2==0),]
summary(control$CPRS)
sd(control$CPRS)  ##0.1988989
male$CPRS_sd<-male$CPRS/0.1988989

table(male$FH_total_2)
#0      1 
#146524  56277 


model11 <- coxph(Surv(survival_time,cancer_total_20)~as.factor(FH_total_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)
summary(model11)
hr <- output1(model11)
hr
#"1.15(1.11-1.18)" "1e-19"


male$prob=predict(model11,newdata=male)
roc2=roc(cancer_total_20~prob,data=male,type="expected")
auc(roc2)#Area under the curve: 0.6862

round(auc(roc2),3)#0.686

ci.auc(roc2,
       conf.level=0.95,
       partial.auc=c(1, .8),  partial.auc.focus="se", partial.auc.correct=TRUE,
       boot.n=10000, stratified=FALSE)
#95% CI: 0.6827-0.6898 (DeLong)


model22 <- coxph(Surv(survival_time,cancer_total_20)~as.factor(FH_total_2)+CPRS_sd+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)
summary(model22)
hr <- output1(model22)
hr
#"1.12(1.09-1.16)" "2.5e-14" 

#coef  exp(coef)   se(coef)       z Pr(>|z|)    
#as.factor(FH_total_2)1               1.164e-01  1.123e+00  1.528e-02   7.621 2.51e-14 ***
#CPRS_sd                              3.116e-01  1.366e+00  7.147e-03  43.598  < 2e-16 ***


male$prob=predict(model22,newdata=male)
roc3=roc(cancer_total_20~prob,data=male,type="expected")
auc(roc3)#Area under the curve: 0.7075

round(auc(roc3),3)#0.707

ci.auc(roc3,
       conf.level=0.95,
       partial.auc=c(1, .8),  partial.auc.focus="se", partial.auc.correct=TRUE,
       boot.n=10000, stratified=FALSE)
#95% CI: 0.704-0.711 (DeLong)

roc.test(roc3,roc2) #拿roc3-roc2
#DeLong's test for two correlated ROC curves

#data:  roc3 and roc2
#Z = 23.567, p-value < 2.2e-16
#alternative hypothesis: true difference in AUC is not equal to 0
#95 percent confidence interval:
# 0.01948693 0.02302228
#sample estimates:
#AUC of roc1 AUC of roc2 
#0.7074741   0.6862195



####异质性####
##要比较model1/2结果的异质性##
res = output(model11)
result[1, ] = c(res[1:6])
res = output(model22)
result[2, ] = c(res[1:6])
result=as.data.frame(lapply(result[,1:6],as.numeric))
metaresult <-
  metagen(
    result$Beta,
    result$se,
    lower = result$lower,
    upper = result$upper,
    sm = "β",
    comb.fixed = FALSE,
    data = result
  )

summary(metaresult)
#Test of heterogeneity:
#  Q d.f. p-value
#0.50    1  0.4795


#female#
control <- female[which(female$FH_total_2==0),]
summary(control$CPRS)
sd(control$CPRS)  ##0.0782589
female$CPRS_sd<-female$CPRS/0.0782589

table(female$FH_total_2)
#0      1 
#167534  72064


model11 <- coxph(Surv(survival_time,cancer_total_20)~as.factor(FH_total_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)
summary(model11)
hr <- output1(model11)
hr
#"1.17(1.13-1.21)" "1.6e-23"

female$prob=predict(model11,newdata=female)
roc2=roc(cancer_total_20~prob,data=female,type="expected")
auc(roc2)#Area under the curve: 0.6078

round(auc(roc2),3)#0.608

ci.auc(roc2,
       conf.level=0.95,
       partial.auc=c(1, .8),  partial.auc.focus="se", partial.auc.correct=TRUE,
       boot.n=10000, stratified=FALSE)
#95% CI: 0.6036-0.612 (DeLong)


model22 <- coxph(Surv(survival_time,cancer_total_20)~as.factor(FH_total_2)+CPRS_sd+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)
summary(model22)
hr <- output1(model22)
hr
#"1.16(1.12-1.19)" "4.8e-20"
#coef  exp(coef)   se(coef)       z Pr(>|z|)    
#as.factor(FH_total_2)1               0.1443356  1.1552717  0.0157415   9.169  < 2e-16 ***
#CPRS_sd                              0.1736008  1.1895806  0.0079006  21.973  < 2e-16 ***


female$prob=predict(model22,newdata=female)
roc3=roc(cancer_total_20~prob,data=female,type="expected")
auc(roc3)#Area under the curve: 0.6173

round(auc(roc3),3)#0.617

ci.auc(roc3,
       conf.level=0.95,
       partial.auc=c(1, .8),  partial.auc.focus="se", partial.auc.correct=TRUE,
       boot.n=10000, stratified=FALSE)
#95% CI: 0.6131-0.6215 (DeLong)

roc.test(roc3,roc2) #拿roc3-roc2
#DeLong's test for two correlated ROC curves

#data:  roc3 and roc2
#Z = 10.515, p-value < 2.2e-16
#alternative hypothesis: true difference in AUC is not equal to 0
#95 percent confidence interval:
#  0.007720672 0.011258311
#sample estimates:
#  AUC of roc1 AUC of roc2 
#0.6172988   0.6078093




##要比较model1/2结果的异质性##s
res = output(model11)
result[1, ] = c(res[1:6])
res = output(model22)
result[2, ] = c(res[1:6])
result=as.data.frame(lapply(result[,1:6],as.numeric))
metaresult <-
  metagen(
    result$Beta,
    result$se,
    lower = result$lower,
    upper = result$upper,
    sm = "β",
    comb.fixed = FALSE,
    data = result
  )

summary(metaresult)
#Test of heterogeneity:
#  Q d.f. p-value
#0.50    1  0.4795
