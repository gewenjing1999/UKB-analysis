###sfigure 5 & sfigure 11###

#敏感性分析--剔除随访第一年患肿瘤的个体

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


#####Sfigure 5#####

#肠癌

#男性
###生成新变量来表示随访第一年内的病例，0：否，1：是
male$CRC<-ifelse(male$CRC_difftime_new<=1&male$CRC_total.y==1,1,0)
###剔除随访第一年内发生的病例
dat1<-subset(male,male$CRC==0) #余下202571人，即在随访第一年患肠癌的男性有202801-202571=230人

#分层分析
control <- dat1[which(dat1$FH_Bowel_2==0),]
summary(control$COL_PRS)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.718   2.910   3.167   3.180   3.437   5.035 

sd(control$COL_PRS)  ##0.3905702
dat1$COL_PRS_sd<-dat1$COL_PRS/0.3905702

Q20=quantile(dat1$COL_PRS_sd,seq(0.05,1,0.05))[[4]]
Q80=quantile(dat1$COL_PRS_sd,seq(0.05,1,0.05))[[16]]

####对原始数据的PRS三分位，<=Q20，1；>=Q80，3
dat1$PRS_1_sd_3_Q20=0
dat1[dat1$COL_PRS_sd<=Q20,]$PRS_1_sd_3_Q20=1
dat1[dat1$COL_PRS_sd>Q20&dat1$COL_PRS_sd<Q80,]$PRS_1_sd_3_Q20=2
dat1[dat1$COL_PRS_sd>=Q80,]$PRS_1_sd_3_Q20=3

table(dat1$PRS_1_sd_3_Q20)
#1      2      3 
#40515 121541  40515

#生成新变量COL_PRS_FH，组合PRS及家族史将人群分为6层
dat1$COL_PRS_FH=0
dat1[dat1$PRS_1_sd_3_Q20==1&dat1$FH_Bowel_2==0,]$COL_PRS_FH=1  #对照（1，2，3，4）（1，2，5，6）
dat1[dat1$PRS_1_sd_3_Q20==1&dat1$FH_Bowel_2==1,]$COL_PRS_FH=2  #家族史效应
dat1[dat1$PRS_1_sd_3_Q20==2&dat1$FH_Bowel_2==0,]$COL_PRS_FH=3  #PRS中位效应
dat1[dat1$PRS_1_sd_3_Q20==2&dat1$FH_Bowel_2==1,]$COL_PRS_FH=4  #PRS中位联合家族史效应
dat1[dat1$PRS_1_sd_3_Q20==3&dat1$FH_Bowel_2==0,]$COL_PRS_FH=5  #PRS高位效应
dat1[dat1$PRS_1_sd_3_Q20==3&dat1$FH_Bowel_2==1,]$COL_PRS_FH=6  #PRS高位联合家族史效应
table(dat1$COL_PRS_FH)
#1      2      3      4      5      6 
#36720   3795 108056  13485  35087   5428 


#case/人年
fit2 <- pyears(Surv(CRC_difftime_new,CRC_total.y)~1,data=dat1[dat1$COL_PRS_FH==1,],scale = 1)
fit2$pyears#
fit2$event#
#round(227/396789.1*100000,2)=57.21

fit2 <- pyears(Surv(CRC_difftime_new,CRC_total.y)~1,data=dat1[dat1$COL_PRS_FH==2,],scale = 1)
fit2$pyears#
fit2$event#  
#round(32/40931.14*100000,2)=78.18

fit2 <- pyears(Surv(CRC_difftime_new,CRC_total.y)~1,data=dat1[dat1$COL_PRS_FH==3,],scale = 1)
fit2$pyears#
fit2$event#  
#round(1267/1165953*100000,2)=108.67

fit2 <- pyears(Surv(CRC_difftime_new,CRC_total.y)~1,data=dat1[dat1$COL_PRS_FH==4,],scale = 1)
fit2$pyears#
fit2$event#  
#round(241/145200*100000,2)=165.98

fit2 <- pyears(Surv(CRC_difftime_new,CRC_total.y)~1,data=dat1[dat1$COL_PRS_FH==5,],scale = 1)
fit2$pyears#377284.1
fit2$event#713   
#round(705/376876.8*100000,2)=187.06

fit2 <- pyears(Surv(CRC_difftime_new,CRC_total.y)~1,data=dat1[dat1$COL_PRS_FH==6,],scale = 1)
fit2$pyears#
fit2$event#   
#round(136/58214.14*100000,2)=233.62


model1 <- coxph(Surv(CRC_difftime_new,CRC_total.y)~as.factor(COL_PRS_FH)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=dat1)
res <- summary(model1)
print(res$coefficients,digits = 2)
print(res$conf.int,digits = 2)
#                                        coef exp(coef) se(coef)      z Pr(>|z|)
#as.factor(COL_PRS_FH)2               0.17741      1.19  0.18890  0.939  3.5e-01
#as.factor(COL_PRS_FH)3               0.62221      1.86  0.07216  8.622  6.6e-18
#as.factor(COL_PRS_FH)4               0.92522      2.52  0.09270  9.980  1.9e-23
#as.factor(COL_PRS_FH)5               1.15475      3.17  0.07654 15.087  2.0e-51
#as.factor(COL_PRS_FH)6               1.28947      3.63  0.10867 11.866  1.8e-32

#                                      exp(coef) exp(-coef) lower .95 upper .95
#as.factor(COL_PRS_FH)2                   1.19       0.84      0.82      1.73
#as.factor(COL_PRS_FH)3                   1.86       0.54      1.62      2.15
#as.factor(COL_PRS_FH)4                   2.52       0.40      2.10      3.03
#as.factor(COL_PRS_FH)5                   3.17       0.32      2.73      3.69
#as.factor(COL_PRS_FH)6                   3.63       0.28      2.93      4.49


#计算家族史（无/有）和遗传风险（低/中/高）之间的相加交互作用
#设HR11表示2 个危险因素同时存在的HR值，HR10,HR01分别表示存在1个因素时HR值。相加交互作用的评价指标：RERI=HR1-HR01-HR10+1,AP=RERI/HR11
##  相加交互作用检验
#coef的2，3，4，5，6分别对应分组时的2-6组，但是却是生成系数的第1-5行


RERI1=epi.interaction(model=model1,coef=c(1,2,3),param = c("dummy"),conf.level=0.95) 
RERI1  
#$reri
#est       lower   upper
#1 0.465269 -0.07144224 1.00198

#$apab
#est       lower     upper
#1 0.1844526 -0.02006019 0.3889654

#$s
#est     lower    upper
#1 1.440111 0.8894637 2.331653

#$multiplicative
#est    lower    upper
#1 2.522431 2.103338 3.025029

RERI2=epi.interaction(model=model1,coef=c(1,4,5),param = c("dummy"),conf.level=0.95)
RERI2
#$reri
#est      lower  upper
#1 0.263507 -0.5160857 1.0431

#$apab
#est      lower     upper
#1 0.07257424 -0.1341515 0.2792999

#$s
#est     lower    upper
#1 1.111309 0.8128484 1.519357

#$multiplicative
#est    lower    upper
#1 3.630861 2.934318 4.492747


###森林图
setwd("Sfigure 5 & 11/中间文件/PRS")

   HR <- c(1.0,1.19,NA,1.86,2.52,NA,3.17,3.63)
lower <- c(1.0,0.82,NA,1.62,2.10,NA,2.73,2.93)
upper <- c(1.0,1.73,NA,2.15,3.03,NA,3.69,4.49)

foredata <- data.frame(HR,lower,upper)
foredata
class(foredata)
write.csv(foredata,file = "Sensitivity_analysis_PRS_FH_foredata_colorectal_male.csv",row.names = F)

pdata <- read.csv('Sensitivity_analysis_PRS_FH_foredata_colorectal_male.csv',header=T)
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
  xticks = c(0.5,1.0,2.0,3.0,4.0,5.0),
)
forestplot

topptx(forestplot,"Sensitivity_analysis_forestplot_colorectal_male_FH_PRS.pptx",width = 5,height = 4)



#女性
###生成新变量来表示随访第一年内的病例，0：否，1：是
female$CRC<-ifelse(female$CRC_difftime_new<=1&female$CRC_total.y==1,1,0)
###剔除随访第一年内发生的病例
dat1<-subset(female,female$CRC==0) #余下239466人，即在随访第一年患肠癌的女性有239598-239466=132人

#分层分析
control <- dat1[which(dat1$FH_Bowel_2==0),]
summary(control$COL_PRS)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.595   2.911   3.168   3.181   3.437   5.200 

sd(control$COL_PRS)  ##0.3912937
dat1$COL_PRS_sd<-dat1$COL_PRS/0.3912937

Q20=quantile(dat1$COL_PRS_sd,seq(0.05,1,0.05))[[4]]
Q80=quantile(dat1$COL_PRS_sd,seq(0.05,1,0.05))[[16]]

####对原始数据的PRS三分位，<=Q20，1；>=Q80，3
dat1$PRS_1_sd_3_Q20=0
dat1[dat1$COL_PRS_sd<=Q20,]$PRS_1_sd_3_Q20=1
dat1[dat1$COL_PRS_sd>Q20&dat1$COL_PRS_sd<Q80,]$PRS_1_sd_3_Q20=2
dat1[dat1$COL_PRS_sd>=Q80,]$PRS_1_sd_3_Q20=3

table(dat1$PRS_1_sd_3_Q20)
#1      2      3 
#47894 143678  47894  

dat1$COL_PRS_FH=0
dat1[dat1$PRS_1_sd_3_Q20==1&dat1$FH_Bowel_2==0,]$COL_PRS_FH=1  #对照（1，2，3，4）（1，2，5，6）
dat1[dat1$PRS_1_sd_3_Q20==1&dat1$FH_Bowel_2==1,]$COL_PRS_FH=2  #家族史效应
dat1[dat1$PRS_1_sd_3_Q20==2&dat1$FH_Bowel_2==0,]$COL_PRS_FH=3  #PRS中位效应
dat1[dat1$PRS_1_sd_3_Q20==2&dat1$FH_Bowel_2==1,]$COL_PRS_FH=4  #PRS中位联合家族史效应
dat1[dat1$PRS_1_sd_3_Q20==3&dat1$FH_Bowel_2==0,]$COL_PRS_FH=5  #PRS高位效应
dat1[dat1$PRS_1_sd_3_Q20==3&dat1$FH_Bowel_2==1,]$COL_PRS_FH=6  #PRS高位联合家族史效应
table(dat1$COL_PRS_FH)
#1      2      3      4      5      6 
#43702   4192 128346  15332  41641   6253  


#case/人年
fit2 <- pyears(Surv(CRC_difftime_new,CRC_total.y)~1,data=dat1[dat1$COL_PRS_FH==1,],scale = 1)
fit2$pyears#
fit2$event#
#round(179/478226*100000,2)=37.43

fit2 <- pyears(Surv(CRC_difftime_new,CRC_total.y)~1,data=dat1[dat1$COL_PRS_FH==2,],scale = 1)
fit2$pyears#
fit2$event#  
#round(26/45805.97*100000,2)=56.76

fit2 <- pyears(Surv(CRC_difftime_new,CRC_total.y)~1,data=dat1[dat1$COL_PRS_FH==3,],scale = 1)
fit2$pyears#
fit2$event#   
#round(973/1406263*100000,2)=69.19

fit2 <- pyears(Surv(CRC_difftime_new,CRC_total.y)~1,data=dat1[dat1$COL_PRS_FH==4,],scale = 1)
fit2$pyears#
fit2$event# 
#round(172/167631.1*100000,2)=102.61

fit2 <- pyears(Surv(CRC_difftime_new,CRC_total.y)~1,data=dat1[dat1$COL_PRS_FH==5,],scale = 1)
fit2$pyears#
fit2$event# 
#round(590/455143.6*100000,2)=129.63

fit2 <- pyears(Surv(CRC_difftime_new,CRC_total.y)~1,data=dat1[dat1$COL_PRS_FH==6,],scale = 1)
fit2$pyears#
fit2$event#  
#round(95/68551.81*100000,2)=138.58


model1 <- coxph(Surv(CRC_difftime_new,CRC_total.y)~as.factor(COL_PRS_FH)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=dat1)
res <- summary(model1)
print(res$coefficients,digits = 2)
print(res$conf.int,digits = 2)
#                                             coef exp(coef) se(coef)      z Pr(>|z|)
#as.factor(COL_PRS_FH)2               3.1e-01      1.36  0.20997  1.457  1.5e-01
#as.factor(COL_PRS_FH)3               6.1e-01      1.84  0.08146  7.478  7.6e-14
#as.factor(COL_PRS_FH)4               9.1e-01      2.48  0.10707  8.480  2.3e-17
#as.factor(COL_PRS_FH)5               1.2e+00      3.46  0.08569 14.472  1.8e-47
#as.factor(COL_PRS_FH)6               1.2e+00      3.36  0.12732  9.528  1.6e-21
#                                      exp(coef) exp(-coef) lower .95 upper .95
#as.factor(COL_PRS_FH)2                   1.36       0.74      0.90      2.05
#as.factor(COL_PRS_FH)3                   1.84       0.54      1.57      2.16
#as.factor(COL_PRS_FH)4                   2.48       0.40      2.01      3.06
#as.factor(COL_PRS_FH)5                   3.46       0.29      2.92      4.09
#as.factor(COL_PRS_FH)6                   3.36       0.30      2.62      4.32


##  相加交互作用检验
#coef的2，3，4，5，6分别对应分组时的2-6组，但是却是生成系数的第1-5行

RERI1=epi.interaction(model=model1,coef=c(1,2,3),param = c("dummy"),conf.level=0.95) 
RERI1  
#$reri
#est      lower     upper
#1 0.282442 -0.3764871 0.9413711

#$apab
#est      lower     upper
#1 0.1139282 -0.1459374 0.3737939

#$s
#est     lower    upper
#1 1.236021 0.7257038 2.105196

#$multiplicative
#est    lower   upper
#1 2.479122 2.009839 3.05798


RERI2=epi.interaction(model=model1,coef=c(1,4,5),param = c("dummy"),conf.level=0.95)
RERI2
#$reri
#est     lower     upper
#1 -0.4497865 -1.373486 0.4739127

#$apab
#est      lower    upper
#1 -0.133703 -0.4249939 0.157588

#$s
#est     lower    upper
#1 0.8401532 0.5875015 1.201456

#$multiplicative
#est   lower    upper
#1 3.364073 2.62112 4.317616


###森林图
setwd("Sfigure 5 & 11/中间文件/PRS")

   HR <- c(1.0,1.36,NA,1.84,2.48,NA,3.46,3.36)
lower <- c(1.0,0.90,NA,1.57,2.01,NA,2.92,2.62)
upper <- c(1.0,2.05,NA,2.16,3.06,NA,4.09,4.32)

foredata <- data.frame(HR,lower,upper)
foredata
class(foredata)
write.csv(foredata,file = "Sensitivity_analysis_PRS_FH_foredata_colorectal_female.csv",row.names = F)

pdata <- read.csv('Sensitivity_analysis_PRS_FH_foredata_colorectal_female.csv',header=T)
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
  xticks = c(0.5,1.0,2.0,3.0,4.0,4.5),
)
forestplot
topptx(forestplot,"Sensitivity_analysis_forestplot_colorectal_female_FH_PRS.pptx",width = 5,height = 4)




#肺癌

#男性
###生成新变量来表示随访第一年内的病例，0：否，1：是
male$LC<-ifelse(male$LC_difftime_new<=1&male$LC_total.y==1,1,0)
###剔除随访第一年内发生的病例
dat1<-subset(male,male$LC==0) #余下202695人，即在随访第一年患肺癌的男性有202801-202695=106人

#分层分析
control <- dat1[which(dat1$FH_Lung_2==0),]
summary(control$LC_PRS)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-0.378934  0.009061  0.065268  0.083836  0.139757  1.019890 

sd(control$LC_PRS)  ##0.1093155
dat1$LC_PRS_sd<-dat1$LC_PRS/0.1093155

Q20=quantile(dat1$LC_PRS_sd,seq(0.05,1,0.05))[[4]]
Q80=quantile(dat1$LC_PRS_sd,seq(0.05,1,0.05))[[16]]

####对原始数据的PRS三分位，<=Q20，1；>=Q80，3
dat1$PRS_7_sd_3_Q20=0
dat1[dat1$LC_PRS_sd<=Q20,]$PRS_7_sd_3_Q20=1
dat1[dat1$LC_PRS_sd>Q20&dat1$LC_PRS_sd<Q80,]$PRS_7_sd_3_Q20=2
dat1[dat1$LC_PRS_sd>=Q80,]$PRS_7_sd_3_Q20=3

table(dat1$PRS_7_sd_3_Q20)
#1      2      3 
#40539 121617  40539  

dat1$LC_PRS_FH=0
dat1[dat1$PRS_7_sd_3_Q20==1&dat1$FH_Lung_2==0,]$LC_PRS_FH=1  #对照（1，2，3，4）（1，2，5，6）
dat1[dat1$PRS_7_sd_3_Q20==1&dat1$FH_Lung_2==1,]$LC_PRS_FH=2  #家族史效应
dat1[dat1$PRS_7_sd_3_Q20==2&dat1$FH_Lung_2==0,]$LC_PRS_FH=3  #PRS中位效应
dat1[dat1$PRS_7_sd_3_Q20==2&dat1$FH_Lung_2==1,]$LC_PRS_FH=4  #PRS中位联合家族史效应
dat1[dat1$PRS_7_sd_3_Q20==3&dat1$FH_Lung_2==0,]$LC_PRS_FH=5  #PRS高位效应
dat1[dat1$PRS_7_sd_3_Q20==3&dat1$FH_Lung_2==1,]$LC_PRS_FH=6  #PRS高位联合家族史效应
table(dat1$LC_PRS_FH)
#1      2      3      4      5      6 
#36063   4476 107026  14591  34972   5567


#case/人年
fit2 <- pyears(Surv(LC_difftime_new,LC_total.y)~1,data=dat1[dat1$LC_PRS_FH==1,],scale = 1)
fit2$pyears#
fit2$event#
#round(166/390528.7*100000,2)=42.51

fit2 <- pyears(Surv(LC_difftime_new,LC_total.y)~1,data=dat1[dat1$LC_PRS_FH==2,],scale = 1)
fit2$pyears#
fit2$event#  
#round(45/48298.28*100000,2)=93.17

fit2 <- pyears(Surv(LC_difftime_new,LC_total.y)~1,data=dat1[dat1$LC_PRS_FH==3,],scale = 1)
fit2$pyears#
fit2$event#   
#round(819/1157948*100000,2)=70.73

fit2 <- pyears(Surv(LC_difftime_new,LC_total.y)~1,data=dat1[dat1$LC_PRS_FH==4,],scale = 1)
fit2$pyears#
fit2$event#  
#round(221/157495.7*100000,2)=140.32

fit2 <- pyears(Surv(LC_difftime_new,LC_total.y)~1,data=dat1[dat1$LC_PRS_FH==5,],scale = 1)
fit2$pyears#
fit2$event#   
#round(356/378456.2*100000,2)=94.07

fit2 <- pyears(Surv(LC_difftime_new,LC_total.y)~1,data=dat1[dat1$LC_PRS_FH==6,],scale = 1)
fit2$pyears#
fit2$event#  
#round(94/59977.42*100000,2)=156.73


model1 <- coxph(Surv(LC_difftime_new,LC_total.y)~as.factor(LC_PRS_FH)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=dat1)
res <- summary(model1)
print(res$coefficients,2)
print(res$conf.int,digits = 2)
#                                             coef exp(coef) se(coef)        z Pr(>|z|)
#as.factor(LC_PRS_FH)2                5.9e-01      1.80  0.16821   3.4988  4.7e-04
#as.factor(LC_PRS_FH)3                5.3e-01      1.69  0.08513   6.1721  6.7e-10
#as.factor(LC_PRS_FH)4                1.0e+00      2.80  0.10290  10.0103  1.4e-23
#as.factor(LC_PRS_FH)5                8.2e-01      2.27  0.09408   8.7235  2.7e-18
#as.factor(LC_PRS_FH)6                1.1e+00      3.16  0.12926   8.8892  6.2e-19

#                                      exp(coef) exp(-coef) lower .95 upper .95
#as.factor(LC_PRS_FH)2                    1.80       0.56      1.30      2.50
#as.factor(LC_PRS_FH)3                    1.69       0.59      1.43      2.00
#as.factor(LC_PRS_FH)4                    2.80       0.36      2.29      3.43
#as.factor(LC_PRS_FH)5                    2.27       0.44      1.89      2.73
#as.factor(LC_PRS_FH)6                    3.16       0.32      2.45      4.06

##  相加交互作用检验
#coef的2，3，4，5，6分别对应分组时的2-6组，但是却是生成系数的第1-5行


RERI1=epi.interaction(model=model1,coef=c(1,2,3),param = c("dummy"),conf.level=0.95) 
RERI1  
#$reri
#est      lower     upper
#1 0.3085313 -0.3533075 0.9703701

#$apab
#est      lower    upper
#1 0.1101469 -0.1216732 0.341967

#$s
#est     lower    upper
#1 1.206713 0.7833457 1.858893

#$multiplicative
#est    lower    upper
#1 2.80109 2.289509 3.426982

RERI2=epi.interaction(model=model1,coef=c(1,4,5),param = c("dummy"),conf.level=0.95)
RERI2
#$reri
#est      lower     upper
#1 0.08167633 -0.7898036 0.9531563

#$apab
#est      lower     upper
#1 0.02588636 -0.2471747 0.2989474

#$s
#est     lower    upper
#1 1.03939 0.6874846 1.571428
#
#$multiplicative
#est    lower    upper
#1 3.155188 2.449041 4.064942


###森林图
setwd("Sfigure 5 & 11/中间文件/PRS")

   HR <- c(1.0,1.80,NA,1.69,2.80,NA,2.27,3.16)
lower <- c(1.0,1.30,NA,1.43,2.29,NA,1.89,2.45)
upper <- c(1.0,2.50,NA,2.00,3.43,NA,2.73,4.06)

foredata <- data.frame(HR,lower,upper)
foredata
class(foredata)
write.csv(foredata,file = "Sensitivity_analysis_PRS_FH_foredata_lung_male.csv",row.names = F)

pdata <- read.csv('Sensitivity_analysis_PRS_FH_foredata_lung_male.csv',header=T)
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
  xticks = c(0.5,1.0,2.0,3.0,4.0,5.0),
)
forestplot
topptx(forestplot,"Sensitivity_analysis_forestplot_lung_male_FH_PRS.pptx",width = 5,height = 4)



#女性
###生成新变量来表示随访第一年内的病例，0：否，1：是
female$LC<-ifelse(female$LC_difftime_new<=1&female$LC_total.y==1,1,0)
###剔除随访第一年内发生的病例
dat1<-subset(female,female$LC==0) #余下239511人，即在随访第一年患肺癌的女性有239598-239511=87人

###lung-female###
control <- dat1[which(dat1$FH_Lung_2==0),]
summary(control$LC_PRS)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-0.401799  0.009541  0.065647  0.084038  0.139876  1.199014 

sd(control$LC_PRS)  ##0.109025
dat1$LC_PRS_sd<-dat1$LC_PRS/0.109025

Q20=quantile(dat1$LC_PRS_sd,seq(0.05,1,0.05))[[4]]
Q80=quantile(dat1$LC_PRS_sd,seq(0.05,1,0.05))[[16]]

####对原始数据的PRS三分位，<=Q20，1；>=Q80，3
dat1$PRS_7_sd_3_Q20=0
dat1[dat1$LC_PRS_sd<=Q20,]$PRS_7_sd_3_Q20=1
dat1[dat1$LC_PRS_sd>Q20&dat1$LC_PRS_sd<Q80,]$PRS_7_sd_3_Q20=2
dat1[dat1$LC_PRS_sd>=Q80,]$PRS_7_sd_3_Q20=3

table(dat1$PRS_7_sd_3_Q20)
#1      2      3 
#47903 143705  47903  

dat1$LC_PRS_FH=0
dat1[dat1$PRS_7_sd_3_Q20==1&dat1$FH_Lung_2==0,]$LC_PRS_FH=1  #对照（1，2，3，4）（1，2，5，6）
dat1[dat1$PRS_7_sd_3_Q20==1&dat1$FH_Lung_2==1,]$LC_PRS_FH=2  #家族史效应
dat1[dat1$PRS_7_sd_3_Q20==2&dat1$FH_Lung_2==0,]$LC_PRS_FH=3  #PRS中位效应
dat1[dat1$PRS_7_sd_3_Q20==2&dat1$FH_Lung_2==1,]$LC_PRS_FH=4  #PRS中位联合家族史效应
dat1[dat1$PRS_7_sd_3_Q20==3&dat1$FH_Lung_2==0,]$LC_PRS_FH=5  #PRS高位效应
dat1[dat1$PRS_7_sd_3_Q20==3&dat1$FH_Lung_2==1,]$LC_PRS_FH=6  #PRS高位联合家族史效应
table(dat1$LC_PRS_FH)
#1      2      3      4      5      6 
#42412   5491 125504  18201  41272   6631 


#case/人年
fit2 <- pyears(Surv(LC_difftime_new,LC_total.y)~1,data=dat1[dat1$LC_PRS_FH==1,],scale = 1)
fit2$pyears#
fit2$event#
round(181/465367.8*100000,2) #38.89

fit2 <- pyears(Surv(LC_difftime_new,LC_total.y)~1,data=dat1[dat1$LC_PRS_FH==2,],scale = 1)
fit2$pyears#
fit2$event#  
round(41/60500.59*100000,2)#67.77

fit2 <- pyears(Surv(LC_difftime_new,LC_total.y)~1,data=dat1[dat1$LC_PRS_FH==3,],scale = 1)
fit2$pyears#
fit2$event#   
round(746/1376662*100000,2)#54.19

fit2 <- pyears(Surv(LC_difftime_new,LC_total.y)~1,data=dat1[dat1$LC_PRS_FH==4,],scale = 1)
fit2$pyears#
fit2$event#  
round(218/199518.4*100000,2)#109.26

fit2 <- pyears(Surv(LC_difftime_new,LC_total.y)~1,data=dat1[dat1$LC_PRS_FH==5,],scale = 1)
fit2$pyears#
fit2$event#   
round(312/452406.5*100000,2)#68.96

fit2 <- pyears(Surv(LC_difftime_new,LC_total.y)~1,data=dat1[dat1$LC_PRS_FH==6,],scale = 1)
fit2$pyears#
fit2$event#  
round(97/72873.27*100000,2)#133.11


model1 <- coxph(Surv(LC_difftime_new,LC_total.y)~as.factor(LC_PRS_FH)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=dat1)
res <- summary(model1)
print(res$coefficients,digits = 2)
print(res$conf.int,digits = 2)
#                                             coef exp(coef) se(coef)        z Pr(>|z|)
#as.factor(LC_PRS_FH)2                3.2e-01      1.38  0.17311   1.8706  6.1e-02
#as.factor(LC_PRS_FH)3                3.4e-01      1.41  0.08288   4.1532  3.3e-05
#as.factor(LC_PRS_FH)4                8.2e-01      2.26  0.10077   8.1129  4.9e-16
#as.factor(LC_PRS_FH)5                6.0e-01      1.82  0.09350   6.4052  1.5e-10
#as.factor(LC_PRS_FH)6                1.1e+00      2.89  0.12598   8.4145  3.9e-17

#                                      exp(coef) exp(-coef) lower .95 upper .95
#as.factor(LC_PRS_FH)2                    1.38       0.72      0.98      1.94
#as.factor(LC_PRS_FH)3                    1.41       0.71      1.20      1.66
#as.factor(LC_PRS_FH)4                    2.26       0.44      1.86      2.76
#as.factor(LC_PRS_FH)5                    1.82       0.55      1.52      2.19
#as.factor(LC_PRS_FH)6                    2.89       0.35      2.25      3.69

##  相加交互作用检验
#coef的2，3，4，5，6分别对应分组时的2-6组，但是却是生成系数的第1-5行

RERI1=epi.interaction(model=model1,coef=c(1,2,3),param = c("dummy"),conf.level=0.95) 
RERI1  

#$reri
#est       lower    upper
#1 0.4716645 -0.06292777 1.006257

#$apab
#est       lower     upper
#1 0.2082469 -0.02000953 0.4365032

#$s
#est     lower    upper
#1 1.594586 0.8334118 3.050958

#$multiplicative
#est    lower    upper
#1 2.26493 1.858997 2.759502

RERI2=epi.interaction(model=model1,coef=c(1,4,5),param = c("dummy"),conf.level=0.95)
RERI2

#$reri
#est       lower    upper
#1 0.6839875 -0.05911348 1.427089

#$apab
#est       lower     upper
#1 0.2369659 0.008911911 0.4650199

#$s
#est     lower    upper
#1 1.568828 0.9311815 2.643116

#$multiplicative
#est   lower    upper
#1 2.886438 2.25492 3.694821

###森林图
setwd("Sfigure 5 & 11/中间文件/PRS")

   HR <- c(1.0,1.38,NA,1.41,2.26,NA,1.82,2.89)
lower <- c(1.0,0.98,NA,1.20,1.86,NA,1.52,2.25)
upper <- c(1.0,1.94,NA,1.66,2.76,NA,2.19,3.69)

foredata <- data.frame(HR,lower,upper)
foredata
class(foredata)
write.csv(foredata,file = "Sensitivity_analysis_PRS_FH_foredata_lung_female.csv",row.names = F)

pdata <- read.csv('Sensitivity_analysis_PRS_FH_foredata_lung_female.csv',header=T)
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
  xticks = c(0.5,1.0,2.0,3.0,4.0),
)
forestplot
topptx(forestplot,"Sensitivity_analysis_forestplot_lung_female_FH_PRS.pptx",width = 5,height = 4)




#前列腺癌

#男性
###生成新变量来表示随访第一年内的病例，0：否，1：是
male$PRC<-ifelse(male$PRC_difftime_new<=1&male$PRC_total.y==1,1,0)
###剔除随访第一年内发生的病例
dat1<-subset(male,male$PRC==0) #余下202186人，即在随访第一年患前列腺癌的男性有202801-202186=615人

#分层分析
control <- dat1[which(dat1$FH_Prostate_2==0),]
summary(control$PRO_PRS)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#5.459   7.735   8.198   8.213   8.673  11.766  

sd(control$PRO_PRS)  ##0.6948325
dat1$PRO_PRS_sd<-dat1$PRO_PRS/0.6948325

Q20=quantile(dat1$PRO_PRS_sd,seq(0.05,1,0.05))[[4]]
Q80=quantile(dat1$PRO_PRS_sd,seq(0.05,1,0.05))[[16]]

####对原始数据的PRS三分位，<=Q20，1；>=Q80，3
dat1$PRS_8_sd_3_Q20=0
dat1[dat1$PRO_PRS_sd<=Q20,]$PRS_8_sd_3_Q20=1
dat1[dat1$PRO_PRS_sd>Q20&dat1$PRO_PRS_sd<Q80,]$PRS_8_sd_3_Q20=2
dat1[dat1$PRO_PRS_sd>=Q80,]$PRS_8_sd_3_Q20=3

table(dat1$PRS_8_sd_3_Q20)
#1      2      3 
#40438 121310  40438  

dat1$PRO_PRS_FH=0
dat1[dat1$PRS_8_sd_3_Q20==1&dat1$FH_Prostate_2==0,]$PRO_PRS_FH=1  #对照（1，2，3，4）（1，2，5，6）
dat1[dat1$PRS_8_sd_3_Q20==1&dat1$FH_Prostate_2==1,]$PRO_PRS_FH=2  #家族史效应
dat1[dat1$PRS_8_sd_3_Q20==2&dat1$FH_Prostate_2==0,]$PRO_PRS_FH=3  #PRS中位效应
dat1[dat1$PRS_8_sd_3_Q20==2&dat1$FH_Prostate_2==1,]$PRO_PRS_FH=4  #PRS中位联合家族史效应
dat1[dat1$PRS_8_sd_3_Q20==3&dat1$FH_Prostate_2==0,]$PRO_PRS_FH=5  #PRS高位效应
dat1[dat1$PRS_8_sd_3_Q20==3&dat1$FH_Prostate_2==1,]$PRO_PRS_FH=6  #PRS高位联合家族史效应
table(dat1$PRO_PRS_FH)
#1      2      3      4      5      6 
#38037   2401 111937   9373  36519   3919 


#case/人年
fit2 <- pyears(Surv(PRC_difftime_new,PRC_total.y)~1,data=dat1[dat1$PRO_PRS_FH==1,],scale = 1)
fit2$pyears#
fit2$event#
round(609/408766.3*100000,2) #148.98

fit2 <- pyears(Surv(PRC_difftime_new,PRC_total.y)~1,data=dat1[dat1$PRO_PRS_FH==2,],scale = 1)
fit2$pyears#
fit2$event#  
round(50/25914.26*100000,2) #192.94

fit2 <- pyears(Surv(PRC_difftime_new,PRC_total.y)~1,data=dat1[dat1$PRO_PRS_FH==3,],scale = 1)
fit2$pyears#
fit2$event#   
round(4190/1194487*100000,2) #350.78

fit2 <- pyears(Surv(PRC_difftime_new,PRC_total.y)~1,data=dat1[dat1$PRO_PRS_FH==4,],scale = 1)
fit2$pyears#
fit2$event#  
round(600/98698.76*100000,2) #607.91

fit2 <- pyears(Surv(PRC_difftime_new,PRC_total.y)~1,data=dat1[dat1$PRO_PRS_FH==5,],scale = 1)
fit2$pyears#
fit2$event#   
round(3013/382350*100000,2) #788.02

fit2 <- pyears(Surv(PRC_difftime_new,PRC_total.y)~1,data=dat1[dat1$PRO_PRS_FH==6,],scale = 1)
fit2$pyears#
fit2$event#  
round(508/40101*100000,2) #1266.81

table(dat1$PRC_total.y)
#0      1 
#193216   8970 


model1 <- coxph(Surv(PRC_difftime_new,PRC_total.y)~as.factor(PRO_PRS_FH)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=dat1)
res <- summary(model1)
print(res$coefficients,digits=2)
print(res$conf.int,digits=2)
#                                             coef exp(coef) se(coef)      z Pr(>|z|)
#as.factor(PRO_PRS_FH)2               2.1e-01      1.23  0.14712  1.411  1.6e-01
#as.factor(PRO_PRS_FH)3               8.7e-01      2.39  0.04338 20.063  1.5e-89
#as.factor(PRO_PRS_FH)4               1.4e+00      4.01  0.05755 24.115 1.8e-128
#as.factor(PRO_PRS_FH)5               1.7e+00      5.57  0.04454 38.539  0.0e+00
#as.factor(PRO_PRS_FH)6               2.2e+00      8.80  0.06018 36.145 4.5e-286
#                                      exp(coef) exp(-coef) lower .95 upper .95
#as.factor(PRO_PRS_FH)2                   1.23       0.81      0.92      1.64
#as.factor(PRO_PRS_FH)3                   2.39       0.42      2.19      2.60
#as.factor(PRO_PRS_FH)4                   4.01       0.25      3.58      4.48
#as.factor(PRO_PRS_FH)5                   5.57       0.18      5.10      6.07
#as.factor(PRO_PRS_FH)6                   8.80       0.11      7.82      9.91

##  相加交互作用检验
#coef的2，3，4，5，6分别对应分组时的2-6组，但是却是生成系数的第1-5行

RERI1=epi.interaction(model=model1,coef=c(1,2,3),param = c("dummy"),conf.level=0.95) 
RERI1  
#$reri
#est     lower    upper
#1 1.387332 0.9125686 1.862096

#$apab
#est     lower     upper
#1 0.3463149 0.2428327 0.4497972

#$s
#est    lower    upper
#1 1.857091 1.445048 2.386624

#$multiplicative
#est   lower    upper
#1 4.005984 3.57868 4.484309

RERI2=epi.interaction(model=model1,coef=c(1,4,5),param = c("dummy"),conf.level=0.95)
RERI2
#$reri
#est    lower    upper
#1 3.007519 2.130878 3.884159

#$apab
#est     lower     upper
#1 0.3416244 0.2683063 0.4149426

#$s
#est    lower   upper
#1 1.627081 1.428497 1.85327

#$multiplicative
#est    lower    upper
#1 8.803582 7.824112 9.905668


###森林图
setwd("Sfigure 5 & 11/中间文件/PRS")

   HR <- c(1.0,1.23,NA,2.39,4.01,NA,5.57,8.80)
lower <- c(1.0,0.92,NA,2.19,3.58,NA,5.10,7.82)
upper <- c(1.0,1.64,NA,2.60,4.48,NA,6.07,9.91)

foredata <- data.frame(HR,lower,upper)
foredata
class(foredata)
write.csv(foredata,file = "Sensitiviry_analysis_PRS_FH_foredata_prostate_male.csv",row.names = F)

pdata <- read.csv('Sensitiviry_analysisPRS_FH_foredata_prostate_male.csv',header=T)
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
  xticks = c(0.5,1.0,3.0,5.0,7.0,9.0,11.0),
)
forestplot
topptx(forestplot,"Sensitiviry_analysis_forestplot_prostate_male_FH_PRS.pptx",width = 5,height = 4)



#女性
###生成新变量来表示随访第一年内的病例，0：否，1：是
female$BC<-ifelse(female$BC_difftime_new<=1&female$BC_total==1,1,0)
###剔除随访第一年内发生的病例
dat1<-subset(female,female$BC==0) #余下238819人，即在随访第一年患乳腺癌的女性有239598-238819=779人

control <- dat1[which(dat1$FH_Breast_2==0),]
summary(control$BC_PRS)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-0.74588 -0.18928 -0.10448 -0.10545 -0.02075  0.50938

sd(control$BC_PRS)  ##0.1271475
dat1$BC_PRS_sd<-dat1$BC_PRS/0.1271475

Q20=quantile(dat1$BC_PRS_sd,seq(0.05,1,0.05))[[4]]
Q80=quantile(dat1$BC_PRS_sd,seq(0.05,1,0.05))[[16]]
####对原始数据的PRS三分位，<=Q20，1；>=Q80，3

dat1$PRS_28_sd_3_Q20=0
dat1[dat1$BC_PRS_sd<=Q20,]$PRS_28_sd_3_Q20=1
dat1[dat1$BC_PRS_sd>Q20&dat1$BC_PRS_sd<Q80,]$PRS_28_sd_3_Q20=2
dat1[dat1$BC_PRS_sd>=Q80,]$PRS_28_sd_3_Q20=3

table(dat1$PRS_28_sd_3_Q20)
#1      2      3 
#47764 143291  47764  

dat1$BC_PRS_FH=0
dat1[dat1$PRS_28_sd_3_Q20==1&dat1$FH_Breast_2==0,]$BC_PRS_FH=1  #对照（1，2，3，4）（1，2，5，6）
dat1[dat1$PRS_28_sd_3_Q20==1&dat1$FH_Breast_2==1,]$BC_PRS_FH=2  #家族史效应
dat1[dat1$PRS_28_sd_3_Q20==2&dat1$FH_Breast_2==0,]$BC_PRS_FH=3  #PRS中位效应
dat1[dat1$PRS_28_sd_3_Q20==2&dat1$FH_Breast_2==1,]$BC_PRS_FH=4  #PRS中位联合家族史效应
dat1[dat1$PRS_28_sd_3_Q20==3&dat1$FH_Breast_2==0,]$BC_PRS_FH=5  #PRS高位效应
dat1[dat1$PRS_28_sd_3_Q20==3&dat1$FH_Breast_2==1,]$BC_PRS_FH=6  #PRS高位联合家族史效应
table(dat1$BC_PRS_FH)
#1      2      3      4      5      6 
#43437   4327 127760  15531  41538   6226 


#case/人年
fit2 <- pyears(Surv(BC_difftime_new,BC_total)~1,data=dat1[dat1$BC_PRS_FH==1,],scale = 1)
fit2$pyears#
fit2$event#
round(784/473958.8*100000,2) #165.42

fit2 <- pyears(Surv(BC_difftime_new,BC_total)~1,data=dat1[dat1$BC_PRS_FH==2,],scale = 1)
fit2$pyears#
fit2$event#  
round(133/46872.67*100000,2) #283.75

fit2 <- pyears(Surv(BC_difftime_new,BC_total)~1,data=dat1[dat1$BC_PRS_FH==3,],scale = 1)
fit2$pyears#
fit2$event#   
round(3647/1385944*100000,2) #263.14

fit2 <- pyears(Surv(BC_difftime_new,BC_total)~1,data=dat1[dat1$BC_PRS_FH==4,],scale = 1)
fit2$pyears#
fit2$event#  
round(659/167216.7*100000,2) #394.10

fit2 <- pyears(Surv(BC_difftime_new,BC_total)~1,data=dat1[dat1$BC_PRS_FH==5,],scale = 1)
fit2$pyears#
fit2$event#   
round(1798/446122.3*100000,2) #403.03

fit2 <- pyears(Surv(BC_difftime_new,BC_total)~1,data=dat1[dat1$BC_PRS_FH==6,],scale = 1)
fit2$pyears#
fit2$event#  
round(396/66368.3*100000,2) #596.67
 
table(dat1$BC_total)
#0      1 
#231402   7417

model1 <- coxph(Surv(BC_difftime_new,BC_total)~as.factor(BC_PRS_FH)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=dat1)
res <- summary(model1)
print(res$coefficients,2)
print(res$conf.int,2)
#                                       coef exp(coef) se(coef)      z Pr(>|z|)
#as.factor(BC_PRS_FH)2                5.2e-01      1.68  0.09379  5.561  2.7e-08
#as.factor(BC_PRS_FH)3                4.8e-01      1.61  0.03939 12.097  1.1e-33
#as.factor(BC_PRS_FH)4                8.6e-01      2.36  0.05287 16.207  4.5e-59
#as.factor(BC_PRS_FH)5                9.2e-01      2.51  0.04293 21.467 3.2e-102
#as.factor(BC_PRS_FH)6                1.3e+00      3.61  0.06168 20.816  3.1e-96
#                                      exp(coef) exp(-coef) lower .95 upper .95
#as.factor(BC_PRS_FH)2                    1.68       0.59      1.40      2.02
#as.factor(BC_PRS_FH)3                    1.61       0.62      1.49      1.74
#as.factor(BC_PRS_FH)4                    2.36       0.42      2.12      2.61
#as.factor(BC_PRS_FH)5                    2.51       0.40      2.31      2.73
#as.factor(BC_PRS_FH)6                    3.61       0.28      3.20      4.07


##  相加交互作用检验
#coef的2，3，4，5，6分别对应分组时的2-6组，但是却是生成系数的第1-5行

RERI1=epi.interaction(model=model1,coef=c(1,2,3),param = c("dummy"),conf.level=0.95) 
RERI1  
#$reri
#est      lower     upper
#1 0.0606248 -0.2877713 0.4090209

#$apab
#est      lower     upper
#1 0.0257366 -0.1215134 0.1729866

#$s
#est     lower    upper
#1 1.046816 0.8014933 1.367227

#$multiplicative
#est    lower    upper
#1 2.355587 2.123733 2.612753



RERI2=epi.interaction(model=model1,coef=c(1,4,5),param = c("dummy"),conf.level=0.95)
RERI2

#reri
#est       lower     upper
#1 0.4129853 -0.05999068 0.8859613

#$apab
#est        lower     upper
#1 0.1143735 -0.009369441 0.2381165

#$s
#est     lower  upper
#1 1.187903 0.9733165 1.4498

#$multiplicative
#est    lower   upper
#1 3.610846 3.199687 4.07484

###森林图
setwd("Sfigure 5 & 11/中间文件/PRS")

   HR <- c(1.0,1.68,NA,1.61,2.36,NA,2.51,3.61)
lower <- c(1.0,1.40,NA,1.49,2.12,NA,2.31,3.20)
upper <- c(1.0,2.02,NA,1.74,2.61,NA,2.73,4.07)

foredata <- data.frame(HR,lower,upper)
foredata
class(foredata)
write.csv(foredata,file = "Sensitiviry_analysis_PRS_FH_foredata_breast_female.csv",row.names = F)

pdata <- read.csv('Sensitiviry_analysis_PRS_FH_foredata_breast_female.csv',header=T)
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
  xticks = c(0.5,1.0,2.0,3.0,4.0,5.0),
)
forestplot
topptx(forestplot,"Sensitiviry_analysis_forestplot_breast_female_FH_PRS.pptx",width = 5,height = 4)





#####Sfigure 11#####

#多肿瘤家族史及CPRS
#男性
###生成新变量来表示随访第一年内的病例，0：否，1：是
male$Overall<-ifelse(male$survival_time<=1&male$cancer_total_20==1,1,0)
###剔除随访第一年内发生的病例
dat1<-subset(male,male$Overall==0) #余下201440人，即在随访第一年患全肿瘤的男性有202801-201440=1361人

control <- dat1[which(dat1$FH_total_2==0),]
summary(control$CPRS)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#2.878   3.603   3.735   3.738   3.870   4.801 
sd(control$CPRS)  ##0.198777
dat1$CPRS_sd<-dat1$CPRS/0.198777


Q20=quantile(dat1$CPRS_sd,seq(0.05,1,0.05))[[4]]
Q80=quantile(dat1$CPRS_sd,seq(0.05,1,0.05))[[16]]

####对原始数据的PRS三分位，<=Q20，1；>=Q80，3
dat1$CPRS_sd_3_Q20=0
dat1[dat1$CPRS_sd<=Q20,]$CPRS_sd_3_Q20=1
dat1[dat1$CPRS_sd>Q20&dat1$CPRS_sd<Q80,]$CPRS_sd_3_Q20=2
dat1[dat1$CPRS_sd>=Q80,]$CPRS_sd_3_Q20=3

table(dat1$CPRS_sd_3_Q20)
#1      2      3 
#40288 120864  40288 

#生成新变量COL_PRS_FH，组合PRS及家族史讲人群分为6层
dat1$CPRS_FH=0
dat1[dat1$CPRS_sd_3_Q20==1&dat1$FH_total_2==0,]$CPRS_FH=1  #对照（1，2，3，4）（1，2，5，6）
dat1[dat1$CPRS_sd_3_Q20==1&dat1$FH_total_2==1,]$CPRS_FH=2  #家族史效应
dat1[dat1$CPRS_sd_3_Q20==2&dat1$FH_total_2==0,]$CPRS_FH=3  #PRS中位效应
dat1[dat1$CPRS_sd_3_Q20==2&dat1$FH_total_2==1,]$CPRS_FH=4  #PRS中位联合家族史效应
dat1[dat1$CPRS_sd_3_Q20==3&dat1$FH_total_2==0,]$CPRS_FH=5  #PRS高位效应
dat1[dat1$CPRS_sd_3_Q20==3&dat1$FH_total_2==1,]$CPRS_FH=6  #PRS高位联合家族史效应
table(dat1$CPRS_FH)
#1      2      3      4      5      6 
#30237 10051 87432 33432 27964 12324


#case/人年
fit2 <- pyears(Surv(survival_time,cancer_total_20)~1,data=dat1[dat1$CPRS_FH==1,],scale = 1)
fit2$pyears#
fit2$event#
round(1720/322001.6*100000,2) #534.16

fit2 <- pyears(Surv(survival_time,cancer_total_20)~1,data=dat1[dat1$CPRS_FH==2,],scale = 1)
fit2$pyears#
fit2$event#  
round(742/106149*100000,2) #699.02

fit2 <- pyears(Surv(survival_time,cancer_total_20)~1,data=dat1[dat1$CPRS_FH==3,],scale = 1)
fit2$pyears#
fit2$event#   
round(7132/922343.1*100000,2) #773.25

fit2 <- pyears(Surv(survival_time,cancer_total_20)~1,data=dat1[dat1$CPRS_FH==4,],scale = 1)
fit2$pyears#
fit2$event#  
round(3297/350563*100000,2) #940.49

fit2 <- pyears(Surv(survival_time,cancer_total_20)~1,data=dat1[dat1$CPRS_FH==5,],scale = 1)
fit2$pyears#
fit2$event#   
round(3512/289608*100000,2) #1234.71

fit2 <- pyears(Surv(survival_time,cancer_total_20)~1,data=dat1[dat1$CPRS_FH==6,],scale = 1)
fit2$pyears#
fit2$event#  
round(1940/125578*100000,2) #1212.67


table(dat1$cancer_total_20)
#0      1 
#183133  18307

model1 <- coxph(Surv(survival_time,cancer_total_20)~as.factor(CPRS_FH)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=dat1)
res <- summary(model1)
print(res$coefficients,digit=2)
print(res$conf.int,digit=2)
#                                         coef exp(coef) se(coef)      z Pr(>|z|)
#as.factor(CPRS_FH)2                  0.14178      1.15  0.04397   3.22  1.3e-03
#as.factor(CPRS_FH)3                  0.36519      1.44  0.02694  13.56  7.4e-42
#as.factor(CPRS_FH)4                  0.46524      1.59  0.02988  15.57  1.2e-54
#as.factor(CPRS_FH)5                  0.83827      2.31  0.02957  28.35 7.9e-177
#as.factor(CPRS_FH)6                  0.97267      2.64  0.03343  29.09 4.6e-186

#                                     exp(coef) exp(-coef) lower .95 upper .95
#as.factor(CPRS_FH)2                      1.15       0.87      1.06      1.26
#as.factor(CPRS_FH)3                      1.44       0.69      1.37      1.52
#as.factor(CPRS_FH)4                      1.59       0.63      1.50      1.69
#as.factor(CPRS_FH)5                      2.31       0.43      2.18      2.45
#as.factor(CPRS_FH)6                      2.64       0.38      2.48      2.82


#计算家族史（无/有）和遗传风险（低/中/高）之间的相加交互作用
#设HR11表示2 个危险因素同时存在的HR值，HR10,HR01分别表示存在1个因素时HR值。相加交互作用的评价指标：RERI=HR1-HR01-HR10+1,AP=RERI/HR11

##  相加交互作用检验
#coef的2，3，4，5，6分别对应分组时的2-6组，但是却是生成系数的第1-5行

RERI1=epi.interaction(model=model1,coef=c(1,2,3),param = c("dummy"),conf.level=0.95) 
RERI1  

#$reri
#est      lower     upper
#1 -0.0007089394 -0.1155403 0.1141225

#$apab
#est       lower      upper
#1 -0.0004452026 -0.07255603 0.07166562

#$s
#est     lower   upper
#1 0.9988047 0.8230028 1.21216

#$multiplicative
#est    lower    upper
#1 1.592397 1.501804 1.688455


RERI2=epi.interaction(model=model1,coef=c(1,4,5),param = c("dummy"),conf.level=0.95)
RERI2
#$reri
#est      lower     upper
#1 0.1803242 0.01195377 0.3486946

#$apab
#est       lower     upper
#1 0.0681755 0.005798831 0.1305522

#$s
#est    lower    upper
#1 1.123115 1.004616 1.255593

#$multiplicative
#est    lower    upper
#1 2.644999 2.477227 2.824134

###森林图
setwd("Sfigure 5 & 11/中间文件/CPRS")

   HR <- c(1.0,1.15,NA,1.44,1.59,NA,2.31,2.64)
lower <- c(1.0,1.06,NA,1.37,1.50,NA,2.18,2.48)
upper <- c(1.0,1.26,NA,1.52,1.69,NA,2.45,2.82)

foredata <- data.frame(HR,lower,upper)
foredata
class(foredata)
write.csv(foredata,file = "Sensitivity_analysis_CPRS_FH_foredata_overall_male.csv",row.names = F)

pdata <- read.csv('Sensitivity_analysis_CPRS_FH_foredata_overall_male.csv',header=T)
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
topptx(forestplot,"Sensitivity_analysis_forestplot_overall_male_FH_CPRS.pptx",width = 5,height = 4)



#女性
###生成新变量来表示随访第一年内的病例，0：否，1：是
female$Overall<-ifelse(female$survival_time<=1&female$cancer_total_20==1,1,0)
###剔除随访第一年内发生的病例
dat1<-subset(female,female$Overall==0) #余下238134人，即在随访第一年患全肿瘤的女性有239598-238134=1464人

#分层分析
control <- dat1[which(dat1$FH_total_2==0),]
summary(control$CPRS)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.070   1.385   1.438   1.437   1.489   1.825 

sd(control$CPRS)  ##0.07823804
dat1$CPRS_sd<-dat1$CPRS/0.07823804

Q20=quantile(dat1$CPRS_sd,seq(0.05,1,0.05))[[4]]
Q80=quantile(dat1$CPRS_sd,seq(0.05,1,0.05))[[16]]

####对原始数据的PRS三分位，<=Q20，1；>=Q80，3
dat1$CPRS_sd_3_Q20=0
dat1[dat1$CPRS_sd<=Q20,]$CPRS_sd_3_Q20=1
dat1[dat1$CPRS_sd>Q20&dat1$CPRS_sd<Q80,]$CPRS_sd_3_Q20=2
dat1[dat1$CPRS_sd>=Q80,]$CPRS_sd_3_Q20=3

table(dat1$CPRS_sd_3_Q20)
#1      2      3 
#47627 142880  47627 

#生成新变量COL_PRS_FH，组合PRS及家族史讲人群分为6层
dat1$CPRS_FH=0
dat1[dat1$CPRS_sd_3_Q20==1&dat1$FH_total_2==0,]$CPRS_FH=1  #对照（1，2，3，4）（1，2，5，6）
dat1[dat1$CPRS_sd_3_Q20==1&dat1$FH_total_2==1,]$CPRS_FH=2  #家族史效应
dat1[dat1$CPRS_sd_3_Q20==2&dat1$FH_total_2==0,]$CPRS_FH=3  #PRS中位效应
dat1[dat1$CPRS_sd_3_Q20==2&dat1$FH_total_2==1,]$CPRS_FH=4  #PRS中位联合家族史效应
dat1[dat1$CPRS_sd_3_Q20==3&dat1$FH_total_2==0,]$CPRS_FH=5  #PRS高位效应
dat1[dat1$CPRS_sd_3_Q20==3&dat1$FH_total_2==1,]$CPRS_FH=6  #PRS高位联合家族史效应
table(dat1$CPRS_FH)
#1      2      3      4      5      6 
#35346 12281 99382 43498 31887 15740   


#case/人年
fit2 <- pyears(Surv(survival_time,cancer_total_20)~1,data=dat1[dat1$CPRS_FH==1,],scale = 1)
fit2$pyears#
fit2$event#
round(1825/378821*100000,2) #481.76

fit2 <- pyears(Surv(survival_time,cancer_total_20)~1,data=dat1[dat1$CPRS_FH==2,],scale = 1)
fit2$pyears#
fit2$event#  
round(792/131599*100000,2) #601.83

fit2 <- pyears(Surv(survival_time,cancer_total_20)~1,data=dat1[dat1$CPRS_FH==3,],scale = 1)
fit2$pyears#
fit2$event#   
round(6372/1067585*100000,2) #596.86

fit2 <- pyears(Surv(survival_time,cancer_total_20)~1,data=dat1[dat1$CPRS_FH==4,],scale = 1)
fit2$pyears#
fit2$event#  
round(3407/464712*100000,2) #733.14

fit2 <- pyears(Surv(survival_time,cancer_total_20)~1,data=dat1[dat1$CPRS_FH==5,],scale = 1)
fit2$pyears#
fit2$event#   
round(2621/339677*100000,2) #771.62

fit2 <- pyears(Surv(survival_time,cancer_total_20)~1,data=dat1[dat1$CPRS_FH==6,],scale = 1)
fit2$pyears#
fit2$event#  
round(1542/166822*100000,2) #924.34


model1 <- coxph(Surv(survival_time,cancer_total_20)~as.factor(CPRS_FH)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=dat1)
res <- summary(model1)
print(res$coefficients,digit=2)
print(res$conf.int,digit=2)

#                                         coef exp(coef) se(coef)      z Pr(>|z|)
#as.factor(CPRS_FH)2                  1.5e-01      1.16  0.04279   3.52  4.3e-04
#as.factor(CPRS_FH)3                  2.0e-01      1.22  0.02750   7.19  6.4e-13
#as.factor(CPRS_FH)4                  3.4e-01      1.41  0.03006  11.40  4.1e-30
#as.factor(CPRS_FH)5                  4.6e-01      1.58  0.03155  14.52  8.5e-48
#as.factor(CPRS_FH)6                  5.8e-01      1.79  0.03557  16.38  2.5e-60

#                                   exp(coef) exp(-coef) lower .95 upper .95
#as.factor(CPRS_FH)2                      1.16       0.86      1.07      1.26
#as.factor(CPRS_FH)3                      1.22       0.82      1.15      1.29
#as.factor(CPRS_FH)4                      1.41       0.71      1.33      1.49
#as.factor(CPRS_FH)5                      1.58       0.63      1.49      1.68
#as.factor(CPRS_FH)6                      1.79       0.56      1.67      1.92


#计算家族史（有/无）和遗传风险（低/中/高）之间的相加交互作用
#设HR11表示2 个危险因素同时存在的HR值，HR10,HR01分别表示存在1个因素时HR值。相加交互作用的评价指标：RERI=HR1-HR01-HR10+1,AP=RERI/HR11

##  相加交互作用检验
#coef的2，3，4，5，6分别对应分组时的2-6组，但是却是生成系数的第1-5行

RERI1=epi.interaction(model=model1,coef=c(1,2,3),param = c("dummy"),conf.level=0.95) 
RERI1  

#$reri
#est       lower     upper
#1 0.02734058 -0.08089522 0.1355764

#$apab
#est      lower      upper
#1 0.01940749 -0.0575273 0.09634227

#$s
#est     lower    upper
#1 1.07168 0.8056295 1.425592

#$multiplicative
#est    lower    upper
#1 1.408764 1.328167 1.494253

RERI2=epi.interaction(model=model1,coef=c(1,4,5),param = c("dummy"),conf.level=0.95)
RERI2
#$reri
#est       lower     upper
#1 0.04707322 -0.09495847 0.1891049

#$apab
#est       lower     upper
#1 0.02628105 -0.05259741 0.1051595

#$s
#est     lower   upper
#1 1.063264 0.8803798 1.28414

#$multiplicative
#est    lower    upper
#1 1.791147 1.670516 1.920489



###森林图
setwd("Sfigure 5 & 11/中间文件/CPRS")

   HR <- c(1.0,1.16,NA,1.22,1.41,NA,1.58,1.79)
lower <- c(1.0,1.07,NA,1.15,1.33,NA,1.49,1.67)
upper <- c(1.0,1.26,NA,1.29,1.49,NA,1.68,1.92)

foredata <- data.frame(HR,lower,upper)
foredata
class(foredata)
write.csv(foredata,file = "Sensitivity_analysis_CPRS_FH_foredata_overall_female.csv",row.names = F)

pdata <- read.csv('Sensitivity_analysis_CPRS_FH_foredata_overall_female.csv',header=T)
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
topptx(forestplot,"Sensitivity_analysis_forestplot_overall_female_FH_CPRS.pptx",width = 5,height = 4)

