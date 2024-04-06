#Stable 3
#PRS与肿瘤发生风险的关联校正家族史

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




########肠癌 CRC model1：PRS(sd)~肿瘤发生风险   model2：PRS(sd)+家族史(有无）~肿瘤发生风险

#male#
control <- male[which(male$FH_Bowel_2==0),]
summary(control$COL_PRS)
sd(control$COL_PRS)  ##0.390638
male$COL_PRS_sd<-male$COL_PRS/0.390638

table(male$FH_Bowel_2)
#0      1 
#180048  22753 


model11 <- coxph(Surv(CRC_difftime_new,CRC_total.y)~COL_PRS_sd+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)
summary(model11)
hr <- output1(model11)
hr
#"1.49(1.44-1.54)" "6.7e-106" 

male$prob=predict(model11,newdata=male)
roc2=roc(CRC_total.y~prob,data=male,type="expected")
auc(roc2)#Area under the curve: 0.7102

round(auc(roc2),3)#0.710

ci.auc(roc2,
       conf.level=0.95,
       partial.auc=c(1, .8),  partial.auc.focus="se", partial.auc.correct=TRUE,
       boot.n=10000, stratified=FALSE)
#95% CI: 0.7014-0.719 (DeLong)


model22 <- coxph(Surv(CRC_difftime_new,CRC_total.y)~COL_PRS_sd+as.factor(FH_Bowel_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)
summary(model22)
hr <- output1(model22)
hr
#"1.48(1.43-1.54)" "8e-104" 
#coef  exp(coef)   se(coef)      z Pr(>|z|)    
#COL_PRS_sd                           0.3932755  1.4818265  0.0181758 21.637  < 2e-16 ***
#as.factor(FH_Bowel_2)1               0.2435461  1.2757652  0.0513575  4.742 2.11e-06 ***

male$prob=predict(model22,newdata=male)
roc3=roc(CRC_total.y~prob,data=male,type="expected")
auc(roc3)#Area under the curve: 0.7115

round(auc(roc3),3)#0.712

ci.auc(roc3,
       conf.level=0.95,
       partial.auc=c(1, .8),  partial.auc.focus="se", partial.auc.correct=TRUE,
       boot.n=10000, stratified=FALSE)
#95% CI: 0.7028-0.7203 (DeLong)

roc.test(roc3,roc2) #拿roc3-roc2
#DeLong's test for two correlated ROC curves

#data:  roc3 and roc2
#Z = 2.4676, p-value = 0.0136
#alternative hypothesis: true difference in AUC is not equal to 0
#95 percent confidence interval:
# 0.0002833352 0.0024711984
#sample estimates:
#AUC of roc1 AUC of roc2 
#0.7115448   0.7101675




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
#0.13    1  0.7237



#female#
control <- female[which(female$FH_Bowel_2==0),]
summary(control$COL_PRS)
sd(control$COL_PRS)  ##0.3913549
female$COL_PRS_sd<-female$COL_PRS/0.3913549

table(female$FH_Bowel_2)
#0      1 
#213800  25798 


model11 <- coxph(Surv(CRC_difftime_new,CRC_total.y)~COL_PRS_sd+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)
summary(model11)
hr <- output1(model11)
hr
#"1.47(1.41-1.53)" "1.2e-76"  


female$prob=predict(model11,newdata=female)
roc2=roc(CRC_total.y~prob,data=female,type="expected")
auc(roc2)#Area under the curve: 0.6833

round(auc(roc2),3)#0.683

ci.auc(roc2,
       conf.level=0.95,
       partial.auc=c(1, .8),  partial.auc.focus="se", partial.auc.correct=TRUE,
       boot.n=10000, stratified=FALSE)
#95% CI: 0.6725-0.694 (DeLong)


model22 <- coxph(Surv(CRC_difftime_new,CRC_total.y)~COL_PRS_sd+as.factor(FH_Bowel_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)
summary(model22)
hr <- output1(model22)
hr
#"1.46(1.41-1.53)" "2.9e-75"
#coef  exp(coef)   se(coef)      z Pr(>|z|)    
#COL_PRS_sd                           3.817e-01  1.465e+00  2.079e-02 18.357  < 2e-16 ***
#as.factor(FH_Bowel_2)1               1.799e-01  1.197e+00  6.127e-02  2.936  0.00332 **  


female$prob=predict(model22,newdata=female)
roc3=roc(CRC_total.y~prob,data=female,type="expected")
auc(roc3)#Area under the curve: 0.6842

round(auc(roc3),3)#0.684

ci.auc(roc3,
       conf.level=0.95,
       partial.auc=c(1, .8),  partial.auc.focus="se", partial.auc.correct=TRUE,
       boot.n=10000, stratified=FALSE)
#95% CI: 0.6735-0.695 (DeLong)

roc.test(roc3,roc2) #拿roc3-roc2

#DeLong's test for two correlated ROC curves

#data:  roc3 and roc2
#Z = 1.7628, p-value = 0.07794
#alternative hypothesis: true difference in AUC is not equal to 0
#95 percent confidence interval:
# -0.0001063003  0.0020066722
#sample estimates:
#AUC of roc1 AUC of roc2 
#0.6842348   0.6832846




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
#0.00    1  1.0000





########肺癌 LC model1：PRS(sd)~肿瘤发生风险   model2：PRS(sd)+家族史(有无）~肿瘤发生风险

#male#
control <- male[which(male$FH_Lung_2==0),]
summary(control$LC_PRS)
sd(control$LC_PRS)  ##0.1093194
male$LC_PRS_sd<-male$LC_PRS/0.1093194

table(male$FH_Lung_2)
#0      1 
#178146  24655 


model11 <- coxph(Surv(LC_difftime_new,LC_total.y)~LC_PRS_sd+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)
summary(model11)
hr <- output1(model11)
hr
#"1.22(1.17-1.27)" "3.1e-22" 

male$prob=predict(model11,newdata=male)
roc2=roc(LC_total.y~prob,data=male,type="expected")
auc(roc2)#Area under the curve: 0.8164

round(auc(roc2),3)#0.816

ci.auc(roc2,
       conf.level=0.95,
       partial.auc=c(1, .8),  partial.auc.focus="se", partial.auc.correct=TRUE,
       boot.n=10000, stratified=FALSE)
#95% CI: 0.8078-0.825 (DeLong)


model22 <- coxph(Surv(LC_difftime_new,LC_total.y)~LC_PRS_sd+as.factor(FH_Lung_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)
summary(model22)
hr <- output1(model22)
hr
#"1.22(1.17-1.27)" "3.6e-21"
#coef  exp(coef)   se(coef)       z Pr(>|z|)    
#LC_PRS_sd                            1.966e-01  1.217e+00  2.082e-02   9.443  < 2e-16 ***
#as.factor(FH_Lung_2)1                4.574e-01  1.580e+00  5.802e-02   7.883 3.19e-15 ***

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
#Z = 3.0976, p-value = 0.001951
#alternative hypothesis: true difference in AUC is not equal to 0
#95 percent confidence interval:
#  0.0007858002 0.0034933034
#sample estimates:
#  AUC of roc1 AUC of roc2 
#0.8185167   0.8163772 




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
#0.00    1  1.0000


#female#
control <- female[which(female$FH_Lung_2==0),]
summary(control$LC_PRS)
sd(control$LC_PRS)  ##0.109029
female$LC_PRS_sd<-female$LC_PRS/0.109029

table(female$FH_Lung_2)
#0      1 
#209249  30349 


model11 <- coxph(Surv(LC_difftime_new,LC_total.y)~LC_PRS_sd+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)
summary(model11)
hr <- output1(model11)
hr
# "1.23(1.17-1.28)" "6e-21"  

female$prob=predict(model11,newdata=female)
roc2=roc(LC_total.y~prob,data=female,type="expected")
auc(roc2)#Area under the curve: 0.805

round(auc(roc2),3)#0.805

ci.auc(roc2,
       conf.level=0.95,
       partial.auc=c(1, .8),  partial.auc.focus="se", partial.auc.correct=TRUE,
       boot.n=10000, stratified=FALSE)
#95% CI: 95% CI: 0.795-0.8151 (DeLong)


model22 <- coxph(Surv(LC_difftime_new,LC_total.y)~LC_PRS_sd+as.factor(FH_Lung_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)
summary(model22)
hr <- output1(model22)
hr
#"1.22(1.17-1.27)" "4.5e-20"
#coef  exp(coef)   se(coef)       z Pr(>|z|)    
#LC_PRS_sd                            0.1991367  1.2203487  0.0217008   9.176  < 2e-16 ***
#as.factor(FH_Lung_2)1                0.4712261  1.6019572  0.0585509   8.048 8.41e-16 ***

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
#Z = 3.1482, p-value = 0.001643
#alternative hypothesis: true difference in AUC is not equal to 0
#95 percent confidence interval:
#0.0009112032 0.0039173178
#sample estimates:
#AUC of roc1 AUC of roc2 
#0.8074423   0.8050280 




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
#0.00    1  1.0000




########前列腺癌 PRC model1：PRS(sd)~肿瘤发生风险   model2：PRS(sd)+家族史(有无）~肿瘤发生风险

#male#
control <- male[which(male$FH_Prostate_2==0),]
summary(control$PRO_PRS)
sd(control$PRO_PRS)  ##0.6953242
male$PRO_PRS_sd<-male$PRO_PRS/0.6953242

table(male$FH_Prostate_2)
#0      1 
#187028  15773 


model11 <- coxph(Surv(PRC_difftime_new,PRC_total.y)~PRO_PRS_sd+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)
summary(model11)
hr <- output1(model11)
hr
#"1.83(1.8-1.87)" "0"

male$prob=predict(model11,newdata=male)
roc2=roc(PRC_total.y~prob,data=male,type="expected")
auc(roc2)#Area under the curve: 0.7602

round(auc(roc2),3)#0.76

ci.auc(roc2,
       conf.level=0.95,
       partial.auc=c(1, .8),  partial.auc.focus="se", partial.auc.correct=TRUE,
       boot.n=10000, stratified=FALSE)
#95% CI: 0.7558-0.7646 (DeLong)


model22 <- coxph(Surv(PRC_difftime_new,PRC_total.y)~PRO_PRS_sd+as.factor(FH_Prostate_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)
summary(model22)
hr <- output1(model22)
hr

#"1.82(1.78-1.86)" "0" 
#coef  exp(coef)   se(coef)      z Pr(>|z|)    
#PRO_PRS_sd                           0.5985804  1.8195340  0.0100327 59.663  < 2e-16 ***
#as.factor(FH_Prostate_2)1            0.4480488  1.5652550  0.0305398 14.671  < 2e-16 ***


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
#Z = 7.26, p-value = 3.871e-13
#alternative hypothesis: true difference in AUC is not equal to 0
#95 percent confidence interval:
#0.002081433 0.003620872
#sample estimates:
#AUC of roc1 AUC of roc2 
#0.7630314   0.7601802




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
control <- female[which(female$FH_Breast_2==0),]
summary(control$BC_PRS)
sd(control$BC_PRS)  ##0.1271848
female$BC_PRS_sd<-female$BC_PRS/0.1271848

table(female$FH_Breast_2)
#0      1 
#213373  26225 


model11 <- coxph(Surv(BC_difftime_new,BC_total)~BC_PRS_sd+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)
summary(model11)
hr <- output1(model11)
hr
#"1.39(1.36-1.43)" "2e-196"

female$prob=predict(model11,newdata=female)
roc2=roc(BC_total~prob,data=female,type="expected")
auc(roc2)#Area under the curve: 0.6117

round(auc(roc2),3)#0.612

ci.auc(roc2,
       conf.level=0.95,
       partial.auc=c(1, .8),  partial.auc.focus="se", partial.auc.correct=TRUE,
       boot.n=10000, stratified=FALSE)
#95% CI: 0.6056-0.6178 (DeLong)


model22 <- coxph(Surv(BC_difftime_new,BC_total)~BC_PRS_sd+as.factor(FH_Breast_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)
summary(model22)
hr <- output1(model22)
hr
#"1.38(1.36-1.42)" "8.3e-188" 

#coef  exp(coef)   se(coef)      z Pr(>|z|)    
#BC_PRS_sd                            0.3256644  1.3849505  0.0111418 29.229  < 2e-16 ***
#as.factor(FH_Breast_2)1              0.3914790  1.4791669  0.0300544 13.026  < 2e-16 ***


female$prob=predict(model22,newdata=female)
roc3=roc(BC_total~prob,data=female,type="expected")
auc(roc3)#Area under the curve:0.6183

round(auc(roc3),3)#0.618

ci.auc(roc3,
       conf.level=0.95,
       partial.auc=c(1, .8),  partial.auc.focus="se", partial.auc.correct=TRUE,
       boot.n=10000, stratified=FALSE)
#95% CI: 0.6122-0.6244 (DeLong)

roc.test(roc3,roc2) #拿roc3-roc2
#DeLong's test for two correlated ROC curves

#data:  roc3 and roc2
#Z = 6.1268, p-value = 8.969e-10
#alternative hypothesis: true difference in AUC is not equal to 0
#95 percent confidence interval:
#0.004477143 0.008689035
#sample estimates:
#AUC of roc1 AUC of roc2 
#0.6182797   0.6116966 




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
#0.00    1  1.0000



#######多肿瘤  CPRS(sd)~肿瘤发生风险   model2：CPRS(sd)+家族史(有无）~肿瘤发生风险

#male#
control <- male[which(male$FH_total_2==0),]
summary(control$CPRS)
sd(control$CPRS)  ##0.1988989
male$CPRS_sd<-male$CPRS/0.1988989

table(male$FH_total_2)
#0      1 
#146524  56277 


model11 <- coxph(Surv(survival_time,cancer_total_20)~CPRS_sd+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)
summary(model11)
hr <- output1(model11)
hr
#"1.37(1.35-1.39)" "0" 

male$prob=predict(model11,newdata=male)
roc2=roc(cancer_total_20~prob,data=male,type="expected")
auc(roc2)#Area under the curve: 0.7069

round(auc(roc2),3)#0.707

ci.auc(roc2,
       conf.level=0.95,
       partial.auc=c(1, .8),  partial.auc.focus="se", partial.auc.correct=TRUE,
       boot.n=10000, stratified=FALSE)
#95% CI: 0.7034-0.7104 (DeLong)


model22 <- coxph(Surv(survival_time,cancer_total_20)~CPRS_sd+as.factor(FH_total_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)
summary(model22)
hr <- output1(model22)
hr
#"1.37(1.35-1.38)" "0"  

#coef  exp(coef)   se(coef)       z Pr(>|z|)    
#CPRS_sd                              3.116e-01  1.366e+00  7.147e-03  43.598  < 2e-16 ***
#as.factor(FH_total_2)1               1.164e-01  1.123e+00  1.528e-02   7.621 2.51e-14 ****


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
#Z = 3.7121, p-value = 0.0002056
#alternative hypothesis: true difference in AUC is not equal to 0
#95 percent confidence interval:
#0.0002611921 0.0008455412
#sample estimates:
#AUC of roc1 AUC of roc2 
#0.7074741   0.7069207 




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
#0.00    1  1.0000




#female#
control <- female[which(female$FH_total_2==0),]
summary(control$CPRS)
sd(control$CPRS)  ##0.0782589
female$CPRS_sd<-female$CPRS/0.0782589

table(female$FH_total_2)
#0      1 
#167534  72064


model11 <- coxph(Surv(survival_time,cancer_total_20)~CPRS_sd+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)
summary(model11)
hr <- output1(model11)
hr
#"1.19(1.17-1.21)" "2.3e-110"


female$prob=predict(model11,newdata=female)
roc2=roc(cancer_total_20~prob,data=female,type="expected")
auc(roc2)#Area under the curve: 0.6157

round(auc(roc2),3)#0.616

ci.auc(roc2,
       conf.level=0.95,
       partial.auc=c(1, .8),  partial.auc.focus="se", partial.auc.correct=TRUE,
       boot.n=10000, stratified=FALSE)
#95% CI: 0.6115-0.6199 (DeLong)


model22 <- coxph(Surv(survival_time,cancer_total_20)~CPRS_sd+as.factor(FH_total_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)
summary(model22)
hr <- output1(model22)
hr
#"1.19(1.17-1.21)" "5.2e-107"  
#coef  exp(coef)   se(coef)       z Pr(>|z|)    
#CPRS_sd                              0.1736008  1.1895806  0.0079006  21.973  < 2e-16 ***
#as.factor(FH_total_2)1               0.1443356  1.1552717  0.0157415   9.169  < 2e-16 ***


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
#Z = 4.2902, p-value = 1.785e-05
#alternative hypothesis: true difference in AUC is not equal to 0
#95 percent confidence interval:
#0.0008571846 0.0022991374
#sample estimates:
#AUC of roc1 AUC of roc2 
#  0.6172988   0.6157206 




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
