#Figure 3


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
library("epiDisplay")
library("epiR")
library("forestplot")
library("eoffice") 



############多肿瘤家族史与CPRS联合效应############

#male#

control <- male[which(male$FH_total_2==0),]
summary(control$CPRS)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#2.878   3.603   3.735   3.738   3.870   4.801
sd(control$CPRS)  ##0.1988989
male$CPRS_sd<-male$CPRS/0.1988989


Q20=quantile(male$CPRS_sd,seq(0.05,1,0.05))[[4]]
Q80=quantile(male$CPRS_sd,seq(0.05,1,0.05))[[16]]

####对原始数据的PRS三分位，<=Q20，1；>=Q80，3
male$CPRS_sd_3_Q20=0
male[male$CPRS_sd<=Q20,]$CPRS_sd_3_Q20=1
male[male$CPRS_sd>Q20&male$CPRS_sd<Q80,]$CPRS_sd_3_Q20=2
male[male$CPRS_sd>=Q80,]$CPRS_sd_3_Q20=3

table(male$CPRS_sd_3_Q20)
#1      2      3 
#40561 121679  40561 

#生成新变量PRS_CRC_FH，组合PRS及家族史讲人群分为6层
male$CPRS_FH=0
male[male$CPRS_sd_3_Q20==1&male$FH_total_2==0,]$CPRS_FH=1  #对照（1，2，3，4）（1，2，5，6）
male[male$CPRS_sd_3_Q20==1&male$FH_total_2==1,]$CPRS_FH=2  #家族史效应
male[male$CPRS_sd_3_Q20==2&male$FH_total_2==0,]$CPRS_FH=3  #PRS中位效应
male[male$CPRS_sd_3_Q20==2&male$FH_total_2==1,]$CPRS_FH=4  #PRS中位联合家族史效应
male[male$CPRS_sd_3_Q20==3&male$FH_total_2==0,]$CPRS_FH=5  #PRS高位效应
male[male$CPRS_sd_3_Q20==3&male$FH_total_2==1,]$CPRS_FH=6  #PRS高位联合家族史效应
table(male$CPRS_FH)
#1      2      3      4      5      6 
#30414 10147 87988 33691 28122 12439 


#case/人年
fit2 <- pyears(Surv(survival_time,cancer_total_20)~1,data=male[male$CPRS_FH==1,],scale = 1)
fit2$pyears#
fit2$event#
#round(1837/322756.6*100000,2)=569.16

fit2 <- pyears(Surv(survival_time,cancer_total_20)~1,data=male[male$CPRS_FH==2,],scale = 1)
fit2$pyears#
fit2$event#  
#round(810/106508.8*100000,2)=760.50

fit2 <- pyears(Surv(survival_time,cancer_total_20)~1,data=male[male$CPRS_FH==3,],scale = 1)
fit2$pyears#
fit2$event#   
#round(7662/922929.1*100000,2)=830.18

fit2 <- pyears(Surv(survival_time,cancer_total_20)~1,data=male[male$CPRS_FH==4,],scale = 1)
fit2$pyears#
fit2$event#  
#round(3554/350743.4*100000,2)=1013.28

fit2 <- pyears(Surv(survival_time,cancer_total_20)~1,data=male[male$CPRS_FH==5,],scale = 1)
fit2$pyears#
fit2$event#   
#round(3756/288732.5*100000,2)=1300.86

fit2 <- pyears(Surv(survival_time,cancer_total_20)~1,data=male[male$CPRS_FH==6,],scale = 1)
fit2$pyears#
fit2$event#  
#round(2049/125956.1*100000,2)=1626.76

1837+810+7662+3554+3756+2049  #19668
table(male$cancer_total_20)
#0      1 
#183133  19668 

model1 <- coxph(Surv(survival_time,cancer_total_20)~as.factor(CPRS_FH)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)
res <- summary(model1)
print(res$coefficients,digit=2)
print(res$conf.int,digit=2)
#                                    coef  exp(coef)   se(coef)       z Pr(>|z|) 
#as.factor(CPRS_FH)2                  0.16050      1.17  0.04222   3.80  1.4e-04
#as.factor(CPRS_FH)3                  0.37115      1.45  0.02605  14.25  4.7e-46
#as.factor(CPRS_FH)4                  0.47352      1.61  0.02887  16.40  1.9e-60
#as.factor(CPRS_FH)5                  0.84136      2.32  0.02860  29.42 3.2e-190
#as.factor(CPRS_FH)6                  0.97961      2.66  0.03229  30.34 3.9e-202

#                                     exp(coef) exp(-coef) lower .95 upper .95
#as.factor(CPRS_FH)2                      1.17       0.85      1.08      1.28
#as.factor(CPRS_FH)3                      1.45       0.69      1.38      1.53
#as.factor(CPRS_FH)4                      1.61       0.62      1.52      1.70
#as.factor(CPRS_FH)5                      2.32       0.43      2.19      2.45
#as.factor(CPRS_FH)6                      2.66       0.38      2.50      2.84


model2 <- coxph(Surv(survival_time,cancer_total_20)~as.numeric(CPRS_FH)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)
res=summary(model2)
round(res$coefficients,2)
round(res$conf.int,2)
#                                   coef  exp(coef)   se(coef)    z Pr(>|z|)    
#as.numeric(CPRS_FH)                  0.20      1.23     0.01  38.43     0.00
#                                      exp(coef) exp(-coef) lower .95 upper .95
#as.numeric(CPRS_FH)                      1.23       0.82      1.21      1.24

#计算家族史（无/有）和遗传风险（低/中/高）之间的相加交互作用
#设HR11表示2 个危险因素同时存在的HR值，HR10,HR01分别表示存在1个因素时HR值。相加交互作用的评价指标：RERI=HR1-HR01-HR10+1,AP=RERI/HR11

##  相加交互作用检验
#coef的2，3，4，5，6分别对应分组时的2-6组，但是却是生成系数的第1-5行

RERI1=epi.interaction(model=model1,coef=c(1,2,3),param = c("dummy"),conf.level=0.95) 
RERI1  

#$reri
#est      lower      upper
#1 -0.01787749 -0.1298895 0.09413455

#$apab
#est       lower      upper
#1 -0.01113425 -0.08085994 0.05859144

#$s
#est     lower   upper
#1 0.9713276 0.8117819 1.16223

#$multiplicative
#est    lower    upper
#1 1.605631 1.517302 1.699102

RERI2=epi.interaction(model=model1,coef=c(1,4,5),param = c("dummy"),conf.level=0.95)
RERI2
#$reri
#est      lower    upper
#1 0.1697992 0.00612131 0.333477

#$apab
#est       lower    upper
#1 0.06375237 0.003448781 0.124056

#$s
#est    lower    upper
#1 1.113683 1.001158 1.238855

#$multiplicative
#est    lower    upper
#1 2.663417 2.500069 2.837438


setwd("D://2022.11.10   PRS & 家族史/20240120最新PRS/分析过程与结果/Figure 4/中间文件")


###森林图
   HR <- c(1.0,1.17,NA,1.45,1.61,NA,2.32,2.66)
lower <- c(1.0,1.08,NA,1.38,1.52,NA,2.19,2.50)
upper <- c(1.0,1.28,NA,1.53,1.70,NA,2.45,2.84)

foredata <- data.frame(HR,lower,upper)
foredata
class(foredata)
write.csv(foredata,file = "CPRS_FH_foredata_overall_male.csv",row.names = F)

pdata <- read.csv('CPRS_FH_foredata_overall_male.csv',header=T)
pdata
colnames(pdata)

labeltext <- as.matrix(pdata[,c(1:3)])#将数据框的XX列转换成矩阵
attach(pdata)#绑定数据框
names(pdata)
##森林图
forestplot <- forestplot(
  labeltext,
  new_page = TRUE,
  graph.pos = 4,
  mean = HR,
  lower = lower,
  upper = upper,
  zero = 1,
  xlog = FALSE,
  fn.ci_norm = fpDrawNormalCI,
  col = fpColors(box = "#00468BFF", line = "#1B1919FF", zero = "gray50"),
  boxsize = 0.4,
  lty.ci = 1,
  lwd.ci = 1,
  ci.vertices = FALSE,
  txt_gp = fpTxtGp(label=gpar(fontfamily="sans"),
                   ticks = gpar(cex = 1),
                   xlab = gpar(cex = 1),
                   cex = 1),
  lineheight = "auto",
  line.margin = .1,
  clip = c(0.95, 1.4),
  xticks = c(0.5,1.0,1.5,2.0,2.5,3.0),
)
forestplot
topptx(forestplot,"forestplot_overall_male_FH_CPRS.pptx",width = 5,height = 4)





###female###
control <- female[which(female$FH_total_2==0),]
summary(control$CPRS)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.070   1.386   1.438   1.437   1.490   1.825 

sd(control$CPRS)  ##0.0782589
female$CPRS_sd<-female$CPRS/0.0782589

Q20=quantile(female$CPRS_sd,seq(0.05,1,0.05))[[4]]
Q80=quantile(female$CPRS_sd,seq(0.05,1,0.05))[[16]]

####对原始数据的PRS三分位，<=Q20，1；>=Q80，3
female$CPRS_sd_3_Q20=0
female[female$CPRS_sd<=Q20,]$CPRS_sd_3_Q20=1
female[female$CPRS_sd>Q20&female$CPRS_sd<Q80,]$CPRS_sd_3_Q20=2
female[female$CPRS_sd>=Q80,]$CPRS_sd_3_Q20=3
table(female$CPRS_sd_3_Q20)
#1      2      3 
#47920 143758  47920 

#生成新变量PRS_CRC_FH，组合PRS及家族史讲人群分为6层
female$CPRS_FH=0
female[female$CPRS_sd_3_Q20==1&female$FH_total_2==0,]$CPRS_FH=1  #对照（1，2，3，4）（1，2，5，6）
female[female$CPRS_sd_3_Q20==1&female$FH_total_2==1,]$CPRS_FH=2  #家族史效应
female[female$CPRS_sd_3_Q20==2&female$FH_total_2==0,]$CPRS_FH=3  #PRS中位效应
female[female$CPRS_sd_3_Q20==2&female$FH_total_2==1,]$CPRS_FH=4  #PRS中位联合家族史效应
female[female$CPRS_sd_3_Q20==3&female$FH_total_2==0,]$CPRS_FH=5  #PRS高位效应
female[female$CPRS_sd_3_Q20==3&female$FH_total_2==1,]$CPRS_FH=6  #PRS高位联合家族史效应
table(female$CPRS_FH)
#1      2      3      4      5      6 
#35545 12375 99936 43822 32053 15867  


#case/人年
fit2 <- pyears(Surv(survival_time,cancer_total_20)~1,data=female[female$CPRS_FH==1,],scale = 1)
fit2$pyears#
fit2$event#
#round(1963/379576.8*100000,2)=517.15

fit2 <- pyears(Surv(survival_time,cancer_total_20)~1,data=female[female$CPRS_FH==2,],scale = 1)
fit2$pyears#
fit2$event#  
#round(862/131911.5*100000,2)=653.47

fit2 <- pyears(Surv(survival_time,cancer_total_20)~1,data=female[female$CPRS_FH==3,],scale = 1)
fit2$pyears#
fit2$event#   
#round(6917/1067962*100000,2)=647.68

fit2 <- pyears(Surv(survival_time,cancer_total_20)~1,data=female[female$CPRS_FH==4,],scale = 1)
fit2$pyears#
fit2$event#  
#round(3735/464855.4*100000,2)=803.48

fit2 <- pyears(Surv(survival_time,cancer_total_20)~1,data=female[female$CPRS_FH==5,],scale = 1)
fit2$pyears#
fit2$event#   
#round(2857/339022.2*100000,2)=842.72

fit2 <- pyears(Surv(survival_time,cancer_total_20)~1,data=female[female$CPRS_FH==6,],scale = 1)
fit2$pyears#
fit2$event#  
#round(1689/166655.9*100000,2)=1013.47


model1 <- coxph(Surv(survival_time,cancer_total_20)~as.factor(CPRS_FH)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)
summary(model1)
res <- summary(model1)
print(res$coefficients,2)
print(res$conf.int,digit=2)

#                                     coef  exp(coef)   se(coef)       z Pr(>|z|)    
#as.factor(CPRS_FH)2                  1.6e-01      1.17  0.04109   3.880  1.0e-04
#as.factor(CPRS_FH)3                  2.1e-01      1.23  0.02648   7.802  6.1e-15
#as.factor(CPRS_FH)4                  3.6e-01      1.43  0.02888  12.474  1.0e-35
#as.factor(CPRS_FH)5                  4.7e-01      1.60  0.03033  15.585  9.2e-55
#as.factor(CPRS_FH)6                  6.0e-01      1.82  0.03413  17.578  3.6e-69

#                                      exp(coef) exp(-coef) lower .95 upper .95
#as.factor(CPRS_FH)2                      1.17       0.85      1.08      1.27
#as.factor(CPRS_FH)3                      1.23       0.81      1.17      1.29
#as.factor(CPRS_FH)4                      1.43       0.70      1.35      1.52
#as.factor(CPRS_FH)5                      1.60       0.62      1.51      1.70
#as.factor(CPRS_FH)6                      1.82       0.55      1.70      1.95


model2 <- coxph(Surv(survival_time,cancer_total_20)~as.numeric(CPRS_FH)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)
res=summary(model2)
print(res$coefficients,2)
print(res$conf.int,digit=2)
#                                      coef  exp(coef)   se(coef)       z Pr(>|z|)    
#as.numeric(CPRS_FH)                  1.2e-01      1.13  0.00564  21.5807 2.7e-103
#                                         exp(coef) exp(-coef) lower .95 upper .95
#as.numeric(CPRS_FH)                      1.13       0.89      1.12      1.14



#计算家族史（有/无）和遗传风险（低/中/高）之间的相加交互作用
#设HR11表示2 个危险因素同时存在的HR值，HR10,HR01分别表示存在1个因素时HR值。相加交互作用的评价指标：RERI=HR1-HR01-HR10+1,AP=RERI/HR11

##  相加交互作用检验
#coef的2，3，4，5，6分别对应分组时的2-6组，但是却是生成系数的第1-5行

RERI1=epi.interaction(model=model1,coef=c(1,2,3),param = c("dummy"),conf.level=0.95) 
RERI1  

#$reri
#est       lower     upper
#1 0.03131968 -0.07331608 0.1359554

#$apab
#est       lower     upper
#1 0.0218466 -0.05125311 0.0949463

#$s
#est     lower    upper
#1 1.077852 0.8296668 1.400278

#$multiplicative
#est    lower    upper
#1 1.433618 1.354732 1.517098

RERI2=epi.interaction(model=model1,coef=c(1,4,5),param = c("dummy"),conf.level=0.95)
RERI2
#$reri
#est       lower     upper
#1 0.04476267 -0.09296329 0.1824886

#$apab
#est      lower      upper
#1 0.02456916 -0.0506508 0.09978912

#$s
#est    lower    upper
#1 1.057599 0.887469 1.260344

#$multiplicative
##1 1.821905 1.704027 1.947937


###森林图
   HR <- c(1.0,1.17,NA,1.23,1.43,NA,1.60,1.82)
lower <- c(1.0,1.08,NA,1.17,1.35,NA,1.51,1.70)
upper <- c(1.0,1.27,NA,1.29,1.52,NA,1.70,1.95)

foredata <- data.frame(HR,lower,upper)
foredata
class(foredata)
write.csv(foredata,file = "CPRS_FH_foredata_overall_female.csv",row.names = F)

pdata <- read.csv('CPRS_FH_foredata_overall_female.csv',header=T)
pdata
colnames(pdata)

labeltext <- as.matrix(pdata[,c(1:3)])#将数据框的XX列转换成矩阵
attach(pdata)#绑定数据框
names(pdata)
##森林图
forestplot <- forestplot(
  labeltext,
  new_page = TRUE,
  graph.pos = 4,
  mean = HR,
  lower = lower,
  upper = upper,
  zero = 1,
  xlog = FALSE,
  fn.ci_norm = fpDrawNormalCI,
  col = fpColors(box = "#00468BFF", line = "#1B1919FF", zero = "gray50"),
  boxsize = 0.4,
  lty.ci = 1,
  lwd.ci = 1,
  ci.vertices = FALSE,
  txt_gp = fpTxtGp(label=gpar(fontfamily="sans"),
                   ticks = gpar(cex = 1),
                   xlab = gpar(cex = 1),
                   cex = 1),
  lineheight = "auto",
  line.margin = .1,
  clip = c(0.95, 1.4),
  xticks = c(0.5,1.0,1.5,2.0),
)
forestplot
topptx(forestplot,"forestplot_overall_female_FH_CPRS.pptx",width = 5,height = 4)
