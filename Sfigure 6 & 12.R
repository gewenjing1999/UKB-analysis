####Sfigure6 & Sfigure 12####

#敏感性分析-排除随访时间小于两年的个体excluding participants who with <2 y of follow-up


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


####Sfigure6

#肠癌
#男性
###剔除随访时间小于两年
dat1<-subset(male,male$CRC_difftime_new>2) #余下201223人，即随访时间小于两年的男性有202801-201223=1578人

control <- dat1[which(dat1$FH_Bowel_2==0),]
summary(control$COL_PRS)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.718   2.910   3.167   3.180   3.437   5.035

sd(control$COL_PRS)  ##0.390637
dat1$COL_PRS_sd<-dat1$COL_PRS/0.390637

Q20=quantile(dat1$COL_PRS_sd,seq(0.05,1,0.05))[[4]]
Q80=quantile(dat1$COL_PRS_sd,seq(0.05,1,0.05))[[16]]

####对原始数据的PRS三分位，<=Q20，1；>=Q80，3
dat1$PRS_1_sd_3_Q20=0
dat1[dat1$COL_PRS_sd<=Q20,]$PRS_1_sd_3_Q20=1
dat1[dat1$COL_PRS_sd>Q20&dat1$COL_PRS_sd<Q80,]$PRS_1_sd_3_Q20=2
dat1[dat1$COL_PRS_sd>=Q80,]$PRS_1_sd_3_Q20=3

table(dat1$PRS_1_sd_3_Q20)
#1      2      3 
#40245 120733  40245  

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
#36481   3764 107334  13399  34855   5390  


#case/人年
fit2 <- pyears(Surv(CRC_difftime_new,CRC_total.y)~1,data=dat1[dat1$COL_PRS_FH==1,],scale = 1)
fit2$pyears#
fit2$event#
round(199/396100*100000,2) #50.24

fit2 <- pyears(Surv(CRC_difftime_new,CRC_total.y)~1,data=dat1[dat1$COL_PRS_FH==2,],scale = 1)
fit2$pyears#
fit2$event# 
round(30/40844*100000,2) #73.45

fit2 <- pyears(Surv(CRC_difftime_new,CRC_total.y)~1,data=dat1[dat1$COL_PRS_FH==3,],scale = 1)
fit2$pyears#
fit2$event# 
round(1159/1165279*100000,2) #99.46

fit2 <- pyears(Surv(CRC_difftime_new,CRC_total.y)~1,data=dat1[dat1$COL_PRS_FH==4,],scale = 1)
fit2$pyears#
fit2$event#  
round(219/145134*100000,2) #150.89

fit2 <- pyears(Surv(CRC_difftime_new,CRC_total.y)~1,data=dat1[dat1$COL_PRS_FH==5,],scale = 1)
fit2$pyears#
fit2$event# 
round(652/376816*100000,2) #173.03

fit2 <- pyears(Surv(CRC_difftime_new,CRC_total.y)~1,data=dat1[dat1$COL_PRS_FH==6,],scale = 1)
fit2$pyears#
fit2$event#   
round(119/58172*100000,2) #204.56

table(dat1$CRC_total.y)
#0      1 
#198845   2378 


model1 <- coxph(Surv(CRC_difftime_new,CRC_total.y)~as.factor(COL_PRS_FH)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=dat1)
res <- summary(model1)
print(res$coefficients,digits = 2)
print(res$conf.int,digits = 2)
#                                             coef exp(coef) se(coef)      z Pr(>|z|)
#as.factor(COL_PRS_FH)2               0.24340      1.28  0.19593  1.242  2.1e-01
#as.factor(COL_PRS_FH)3               0.66238      1.94  0.07682  8.622  6.6e-18
#as.factor(COL_PRS_FH)4               0.95693      2.60  0.09815  9.750  1.8e-22
#as.factor(COL_PRS_FH)5               1.20559      3.34  0.08121 14.845  7.5e-50
#as.factor(COL_PRS_FH)6               1.28462      3.61  0.11612 11.063  1.9e-28

#                                      exp(coef) exp(-coef) lower .95 upper .95
#as.factor(COL_PRS_FH)2                   1.28       0.78      0.87      1.87
#as.factor(COL_PRS_FH)3                   1.94       0.52      1.67      2.25
#as.factor(COL_PRS_FH)4                   2.60       0.38      2.15      3.16
#as.factor(COL_PRS_FH)5                   3.34       0.30      2.85      3.91
#as.factor(COL_PRS_FH)6                   3.61       0.28      2.88      4.54


#计算家族史（无/有）和遗传风险（低/中/高）之间的相加交互作用
#设HR11表示2 个危险因素同时存在的HR值，HR10,HR01分别表示存在1个因素时HR值。相加交互作用的评价指标：RERI=HR1-HR01-HR10+1,AP=RERI/HR11
##  相加交互作用检验
#coef的2，3，4，5，6分别对应分组时的2-6组，但是却是生成系数的第1-5行


RERI1=epi.interaction(model=model1,coef=c(1,2,3),param = c("dummy"),conf.level=0.95) 
RERI1  
#$reri
#est      lower     upper
#1 0.3887236 -0.2003849 0.9778321

#$apab
#est       lower     upper
#1 0.149297 -0.07007024 0.3686642

#$s
#est     lower    upper
#1 1.319945 0.8292685 2.100954

#$multiplicative
#est    lower   upper
#1 2.603694 2.148056 3.15598



RERI2=epi.interaction(model=model1,coef=c(1,4,5),param = c("dummy"),conf.level=0.95)
RERI2
#$reri
#est      lower     upper
#1 -0.001002753 -0.8465807 0.8445752

#$apab
#est      lower     upper
#1 -0.0002775171 -0.2343277 0.2337726

#$s
#est     lower    upper
#1 0.9996164 0.7233315 1.381432

#$multiplicative
#est    lower    upper
#1 3.6133 2.877826 4.536737


###森林图
setwd("Sfigure 6 & 12/中间文件/PRS")

   HR <- c(1.0,1.28,NA,1.94,2.60,NA,3.34,3.61)
lower <- c(1.0,0.87,NA,1.67,2.15,NA,2.85,2.88)
upper <- c(1.0,1.87,NA,2.25,3.16,NA,3.91,4.54)

foredata <- data.frame(HR,lower,upper)
foredata
class(foredata)
write.csv(foredata,file = "Sensitivity_analysis_excluding_participants_with_less_than_two_years_follow_up_PRS_FH_foredata_colorectal_male.csv",row.names = F)

pdata <- read.csv('Sensitivity_analysis_excluding_participants_with_less_than_two_years_follow_up_PRS_FH_foredata_colorectal_male.csv',header=T)
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
  xticks = c(0.5,1.0,2.0,3.0,4.0,5.0)
)
forestplot

topptx(forestplot,"Sensitivity_analysis_excluding_participants_with_less_than_two_years_follow_up_PRS_FH_foredata_colorectal_male.pptx",width = 5,height = 4)



###colorectal-female###
#女性
###剔除随访时间小于两年
dat1<-subset(female,female$CRC_difftime_new>2) #余下238673，即随访时间小于两年的女性有239598-238673=925人

control <- dat1[which(dat1$FH_Bowel_2==0),]
summary(control$COL_PRS)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.595   2.911   3.168   3.181   3.437   5.200

sd(control$COL_PRS)  ##0.3911737
dat1$COL_PRS_sd<-dat1$COL_PRS/0.3911737

Q20=quantile(dat1$COL_PRS_sd,seq(0.05,1,0.05))[[4]]
Q80=quantile(dat1$COL_PRS_sd,seq(0.05,1,0.05))[[16]]

####对原始数据的PRS三分位，<=Q20，1；>=Q80，3
dat1$PRS_1_sd_3_Q20=0
dat1[dat1$COL_PRS_sd<=Q20,]$PRS_1_sd_3_Q20=1
dat1[dat1$COL_PRS_sd>Q20&dat1$COL_PRS_sd<Q80,]$PRS_1_sd_3_Q20=2
dat1[dat1$COL_PRS_sd>=Q80,]$PRS_1_sd_3_Q20=3

table(dat1$PRS_1_sd_3_Q20)
#1      2      3 
#47735 143203  47735  

dat1$COL_PRS_FH=0
dat1[dat1$PRS_1_sd_3_Q20==1&dat1$FH_Bowel_2==0,]$COL_PRS_FH=1  #对照（1，2，3，4）（1，2，5，6）
dat1[dat1$PRS_1_sd_3_Q20==1&dat1$FH_Bowel_2==1,]$COL_PRS_FH=2  #家族史效应
dat1[dat1$PRS_1_sd_3_Q20==2&dat1$FH_Bowel_2==0,]$COL_PRS_FH=3  #PRS中位效应
dat1[dat1$PRS_1_sd_3_Q20==2&dat1$FH_Bowel_2==1,]$COL_PRS_FH=4  #PRS中位联合家族史效应
dat1[dat1$PRS_1_sd_3_Q20==3&dat1$FH_Bowel_2==0,]$COL_PRS_FH=5  #PRS高位效应
dat1[dat1$PRS_1_sd_3_Q20==3&dat1$FH_Bowel_2==1,]$COL_PRS_FH=6  #PRS高位联合家族史效应
table(dat1$COL_PRS_FH)
#1      2      3      4      5      6 
#43553   4182 127929  15274  41497   6238 


#case/人年
fit2 <- pyears(Surv(CRC_difftime_new,CRC_total.y)~1,data=dat1[dat1$COL_PRS_FH==1,],scale = 1)
fit2$pyears#
fit2$event#
round(170/477948*100000,2) #35.57

fit2 <- pyears(Surv(CRC_difftime_new,CRC_total.y)~1,data=dat1[dat1$COL_PRS_FH==2,],scale = 1)
fit2$pyears#
fit2$event#  
round(23/45773*100000,2) #50.25

fit2 <- pyears(Surv(CRC_difftime_new,CRC_total.y)~1,data=dat1[dat1$COL_PRS_FH==3,],scale = 1)
fit2$pyears#
fit2$event#  
round(890/1405436*100000,2) #63.33

fit2 <- pyears(Surv(CRC_difftime_new,CRC_total.y)~1,data=dat1[dat1$COL_PRS_FH==4,],scale = 1)
fit2$pyears#
fit2$event#  
round(164/167509*100000,2) #97.91

fit2 <- pyears(Surv(CRC_difftime_new,CRC_total.y)~1,data=dat1[dat1$COL_PRS_FH==5,],scale = 1)
fit2$pyears#
fit2$event#
round(527/455364*100000,2) #115.73

fit2 <- pyears(Surv(CRC_difftime_new,CRC_total.y)~1,data=dat1[dat1$COL_PRS_FH==6,],scale = 1)
fit2$pyears#
fit2$event#  
round(88/68601*100000,2) #128.28


model1 <- coxph(Surv(CRC_difftime_new,CRC_total.y)~as.factor(COL_PRS_FH)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=dat1)
res <- summary(model1)
print(res$coefficients,digits = 2)
print(res$conf.int,digits = 2)
#                                             coef exp(coef) se(coef)      z Pr(>|z|)
#as.factor(COL_PRS_FH)2               0.23849      1.27  0.22226  1.073  2.8e-01
#as.factor(COL_PRS_FH)3               0.57328      1.77  0.08384  6.838  8.1e-12
#as.factor(COL_PRS_FH)4               0.91572      2.50  0.10977  8.342  7.3e-17
#as.factor(COL_PRS_FH)5               1.18099      3.26  0.08859 13.331  1.5e-40
#as.factor(COL_PRS_FH)6               1.19208      3.29  0.13174  9.049  1.4e-19
#                                     exp(coef) exp(-coef) lower .95 upper .95
#as.factor(COL_PRS_FH)2                   1.27       0.79      0.82      1.96
#as.factor(COL_PRS_FH)3                   1.77       0.56      1.51      2.09
#as.factor(COL_PRS_FH)4                   2.50       0.40      2.01      3.10
#as.factor(COL_PRS_FH)5                   3.26       0.31      2.74      3.88
#as.factor(COL_PRS_FH)6                   3.29       0.30      2.54      4.26


##  相加交互作用检验
#coef的2，3，4，5，6分别对应分组时的2-6组，但是却是生成系数的第1-5行

RERI1=epi.interaction(model=model1,coef=c(1,2,3),param = c("dummy"),conf.level=0.95) 
RERI1  
#$reri
#est      lower   upper
#1 0.455164 -0.2049821 1.11531

#$apab
#est       lower     upper
#1 0.1821694 -0.07194945 0.4362883

#$s
#est     lower    upper
#1 1.436227 0.7881635 2.617157

#$multiplicative
#est    lower    upper
#1 2.498575 2.014899 3.098357


RERI2=epi.interaction(model=model1,coef=c(1,4,5),param = c("dummy"),conf.level=0.95)
RERI2
#$reri
#est     lower     upper
#1 -0.2329857 -1.157528 0.6915566

#$apab
#est      lower     upper
#1 -0.07073172 -0.3607292 0.2192658

#$s
#est     lower    upper
#1 0.9077986 0.6185408 1.332326

#$multiplicative
#est    lower    upper
#1 3.293936 2.544342 4.264369

###森林图
setwd("Sfigure 6 & 12/中间文件/PRS")

   HR <- c(1.0,1.27,NA,1.77,2.50,NA,3.26,3.29)
lower <- c(1.0,0.82,NA,1.51,2.01,NA,2.74,2.54)
upper <- c(1.0,1.96,NA,2.09,3.10,NA,3.88,4.26)

foredata <- data.frame(HR,lower,upper)
foredata
class(foredata)
write.csv(foredata,file = "Sensitivity_analysis_excluding_participants_with_less_than_two_years_follow_up_PRS_FH_foredata_colorectal_female.csv",row.names = F)

pdata <- read.csv('Sensitivity_analysis_excluding_participants_with_less_than_two_years_follow_up_PRS_FH_foredata_colorectal_female.csv',header=T)
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
topptx(forestplot,"Sensitivity_analysis_excluding_participants_with_less_than_two_years_follow_up_PRS_FH_foredata_colorectal_female.pptx",width = 5,height = 4)



#肺癌

#男性
###剔除随访时间小于两年
dat1<-subset(male,male$LC_difftime_new>2) #余下201515人，即随访时间小于两年的男性有202801-201515=1286人

#分层分析
control <- dat1[which(dat1$FH_Lung_2==0),]
summary(control$LC_PRS)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-0.378934  0.009052  0.065273  0.083825  0.139754  1.019890 

sd(control$LC_PRS)  ##0.1093029
dat1$LC_PRS_sd<-dat1$LC_PRS/0.1093029

Q20=quantile(dat1$LC_PRS_sd,seq(0.05,1,0.05))[[4]]
Q80=quantile(dat1$LC_PRS_sd,seq(0.05,1,0.05))[[16]]

####对原始数据的PRS三分位，<=Q20，1；>=Q80，3
dat1$PRS_7_sd_3_Q20=0
dat1[dat1$LC_PRS_sd<=Q20,]$PRS_7_sd_3_Q20=1
dat1[dat1$LC_PRS_sd>Q20&dat1$LC_PRS_sd<Q80,]$PRS_7_sd_3_Q20=2
dat1[dat1$LC_PRS_sd>=Q80,]$PRS_7_sd_3_Q20=3

table(dat1$PRS_7_sd_3_Q20)
#1      2      3 
#40303 120909  40303  

dat1$LC_PRS_FH=0
dat1[dat1$PRS_7_sd_3_Q20==1&dat1$FH_Lung_2==0,]$LC_PRS_FH=1  #对照（1，2，3，4）（1，2，5，6）
dat1[dat1$PRS_7_sd_3_Q20==1&dat1$FH_Lung_2==1,]$LC_PRS_FH=2  #家族史效应
dat1[dat1$PRS_7_sd_3_Q20==2&dat1$FH_Lung_2==0,]$LC_PRS_FH=3  #PRS中位效应
dat1[dat1$PRS_7_sd_3_Q20==2&dat1$FH_Lung_2==1,]$LC_PRS_FH=4  #PRS中位联合家族史效应
dat1[dat1$PRS_7_sd_3_Q20==3&dat1$FH_Lung_2==0,]$LC_PRS_FH=5  #PRS高位效应
dat1[dat1$PRS_7_sd_3_Q20==3&dat1$FH_Lung_2==1,]$LC_PRS_FH=6  #PRS高位联合家族史效应
table(dat1$LC_PRS_FH)
#1      2      3      4      5      6 
#35865   4438 106416  14493  34768   5535


#case/人年
fit2 <- pyears(Surv(LC_difftime_new,LC_total.y)~1,data=dat1[dat1$LC_PRS_FH==1,],scale = 1)
fit2$pyears#
fit2$event#
round(157/390383*100000,2) #40.22

fit2 <- pyears(Surv(LC_difftime_new,LC_total.y)~1,data=dat1[dat1$LC_PRS_FH==2,],scale = 1)
fit2$pyears#
fit2$event#  
round(39/48268*100000,2) #80.80

fit2 <- pyears(Surv(LC_difftime_new,LC_total.y)~1,data=dat1[dat1$LC_PRS_FH==3,],scale = 1)
fit2$pyears#
fit2$event#   
round(753/1157233*100000,2) #65.07

fit2 <- pyears(Surv(LC_difftime_new,LC_total.y)~1,data=dat1[dat1$LC_PRS_FH==4,],scale = 1)
fit2$pyears#
fit2$event#  
round(213/157372*100000,2) #135.35

fit2 <- pyears(Surv(LC_difftime_new,LC_total.y)~1,data=dat1[dat1$LC_PRS_FH==5,],scale = 1)
fit2$pyears#
fit2$event#   
round(338/378151*100000,2) #89.38

fit2 <- pyears(Surv(LC_difftime_new,LC_total.y)~1,data=dat1[dat1$LC_PRS_FH==6,],scale = 1)
fit2$pyears#
fit2$event#  
round(89/59932*100000,2) #148.50


model1 <- coxph(Surv(LC_difftime_new,LC_total.y)~as.factor(LC_PRS_FH)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=dat1)
res <- summary(model1)
print(res$coefficients,digits = 2)
print(res$conf.int,digits = 2)
#                                             coef exp(coef) se(coef)       z Pr(>|z|)
#as.factor(LC_PRS_FH)2                0.50565      1.66  0.17906   2.824  4.7e-03
#as.factor(LC_PRS_FH)3                0.49715      1.64  0.08775   5.666  1.5e-08
#as.factor(LC_PRS_FH)4                1.05363      2.87  0.10539   9.998  1.6e-23
#as.factor(LC_PRS_FH)5                0.82449      2.28  0.09669   8.528  1.5e-17
#as.factor(LC_PRS_FH)6                1.15532      3.18  0.13288   8.695  3.5e-18

#                                      exp(coef) exp(-coef) lower .95 upper .95
#as.factor(LC_PRS_FH)2                    1.66       0.60      1.17      2.36
#as.factor(LC_PRS_FH)3                    1.64       0.61      1.38      1.95
#as.factor(LC_PRS_FH)4                    2.87       0.35      2.33      3.53
#as.factor(LC_PRS_FH)5                    2.28       0.44      1.89      2.76
#as.factor(LC_PRS_FH)6                    3.18       0.31      2.45      4.12

##  相加交互作用检验
#coef的2，3，4，5，6分别对应分组时的2-6组，但是却是生成系数的第1-5行


RERI1=epi.interaction(model=model1,coef=c(1,2,3),param = c("dummy"),conf.level=0.95) 
RERI1  
#$reri
#est       lower    upper
#1 0.5659351 -0.09571505 1.227585

#$apab
#est       lower     upper
#1 0.1973254 -0.02449144 0.4191421

#$s
#est     lower    upper
#1 1.434634 0.8835998 2.329307

#$multiplicative
#est    lower    upper
#1 2.86803 2.332801 3.526061



RERI2=epi.interaction(model=model1,coef=c(1,4,5),param = c("dummy"),conf.level=0.95)
RERI2
#$reri
#est      lower    upper
#1 0.2362579 -0.6468112 1.119327

#$apab
#est      lower     upper
#1 0.07441099 -0.1942386 0.3430606

#$s
#est     lower    upper
#1 1.121859 0.7272953 1.730477

#$multiplicative
#est    lower    upper
#1 3.17504 2.447049 4.119606




###森林图
setwd("Sfigure 6 & 12/中间文件/PRS")


   HR <- c(1.0,1.66,NA,1.64,2.87,NA,2.28,3.18)
lower <- c(1.0,1.17,NA,1.38,2.33,NA,1.89,2.45)
upper <- c(1.0,2.36,NA,1.95,3.53,NA,2.76,4.12)

foredata <- data.frame(HR,lower,upper)
foredata
class(foredata)
write.csv(foredata,file = "Sensitivity_analysis_excluding_participants_with_less_than_two_years_follow_up_PRS_FH_foredata_lung_male.csv",row.names = F)

pdata <- read.csv('Sensitivity_analysis_excluding_participants_with_less_than_two_years_follow_up_PRS_FH_foredata_lung_male.csv',header=T)
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
topptx(forestplot,"Sensitivity_analysis_excluding_participants_with_less_than_two_years_follow_up_PRS_FH_foredata_lung_male.pptx",width = 5,height = 4)



#女性
###剔除随访时间小于两年
dat1<-subset(female,female$LC_difftime_new>2) #余下238850人，即随访时间小于两年的女性有239598-238850=748人

###分层分析
control <- dat1[which(dat1$FH_Lung_2==0),]
summary(control$LC_PRS)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.401799  0.009532  0.065641  0.084031  0.139869  1.199014

sd(control$LC_PRS)  ##0.109038
dat1$LC_PRS_sd<-dat1$LC_PRS/0.109038

Q20=quantile(dat1$LC_PRS_sd,seq(0.05,1,0.05))[[4]]
Q80=quantile(dat1$LC_PRS_sd,seq(0.05,1,0.05))[[16]]

####对原始数据的PRS三分位，<=Q20，1；>=Q80，3
dat1$PRS_7_sd_3_Q20=0
dat1[dat1$LC_PRS_sd<=Q20,]$PRS_7_sd_3_Q20=1
dat1[dat1$LC_PRS_sd>Q20&dat1$LC_PRS_sd<Q80,]$PRS_7_sd_3_Q20=2
dat1[dat1$LC_PRS_sd>=Q80,]$PRS_7_sd_3_Q20=3

table(dat1$PRS_7_sd_3_Q20)
#1      2      3 
#47770 143310  47770  

dat1$LC_PRS_FH=0
dat1[dat1$PRS_7_sd_3_Q20==1&dat1$FH_Lung_2==0,]$LC_PRS_FH=1  #对照（1，2，3，4）（1，2，5，6）
dat1[dat1$PRS_7_sd_3_Q20==1&dat1$FH_Lung_2==1,]$LC_PRS_FH=2  #家族史效应
dat1[dat1$PRS_7_sd_3_Q20==2&dat1$FH_Lung_2==0,]$LC_PRS_FH=3  #PRS中位效应
dat1[dat1$PRS_7_sd_3_Q20==2&dat1$FH_Lung_2==1,]$LC_PRS_FH=4  #PRS中位联合家族史效应
dat1[dat1$PRS_7_sd_3_Q20==3&dat1$FH_Lung_2==0,]$LC_PRS_FH=5  #PRS高位效应
dat1[dat1$PRS_7_sd_3_Q20==3&dat1$FH_Lung_2==1,]$LC_PRS_FH=6  #PRS高位联合家族史效应
table(dat1$LC_PRS_FH)
#1      2      3      4      5      6 
#42295   5475 125171  18139  41153   6617


#case/人年
fit2 <- pyears(Surv(LC_difftime_new,LC_total.y)~1,data=dat1[dat1$LC_PRS_FH==1,],scale = 1)
fit2$pyears#
fit2$event#
round(171/465204*100000,2) #36.76

fit2 <- pyears(Surv(LC_difftime_new,LC_total.y)~1,data=dat1[dat1$LC_PRS_FH==2,],scale = 1)
fit2$pyears#
fit2$event#  
round(37/60483*100000,2) #61.17

fit2 <- pyears(Surv(LC_difftime_new,LC_total.y)~1,data=dat1[dat1$LC_PRS_FH==3,],scale = 1)
fit2$pyears#
fit2$event#   
round(688/1376368*100000,2) #49.99

fit2 <- pyears(Surv(LC_difftime_new,LC_total.y)~1,data=dat1[dat1$LC_PRS_FH==4,],scale = 1)
fit2$pyears#
fit2$event#  
round(205/199447*100000,2) #102.78

fit2 <- pyears(Surv(LC_difftime_new,LC_total.y)~1,data=dat1[dat1$LC_PRS_FH==5,],scale = 1)
fit2$pyears#
fit2$event#   
round(299/452164*100000,2) #66.13

fit2 <- pyears(Surv(LC_difftime_new,LC_total.y)~1,data=dat1[dat1$LC_PRS_FH==6,],scale = 1)
fit2$pyears#
fit2$event#  
round(92/72848*100000,2) #126.29


model1 <- coxph(Surv(LC_difftime_new,LC_total.y)~as.factor(LC_PRS_FH)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=dat1)
res <- summary(model1)
print(res$coefficients,digits = 2)
print(res$conf.int,digits = 2)
#                                             coef exp(coef) se(coef)       z Pr(>|z|)
#as.factor(LC_PRS_FH)2                0.27601      1.32  0.18147   1.521  1.3e-01
#as.factor(LC_PRS_FH)3                0.32014      1.38  0.08547   3.746  1.8e-04
#as.factor(LC_PRS_FH)4                0.81209      2.25  0.10379   7.824  5.1e-15
#as.factor(LC_PRS_FH)5                0.61562      1.85  0.09594   6.417  1.4e-10
#as.factor(LC_PRS_FH)6                1.06388      2.90  0.12945   8.219  2.1e-16

#                                      exp(coef) exp(-coef) lower .95 upper .95
#as.factor(LC_PRS_FH)2                    1.32       0.76      0.92      1.88
#as.factor(LC_PRS_FH)3                    1.38       0.73      1.16      1.63
#as.factor(LC_PRS_FH)4                    2.25       0.44      1.84      2.76
#as.factor(LC_PRS_FH)5                    1.85       0.54      1.53      2.23
#as.factor(LC_PRS_FH)6                    2.90       0.35      2.25      3.73

##  相加交互作用检验
#coef的2，3，4，5，6分别对应分组时的2-6组，但是却是生成系数的第1-5行

RERI1=epi.interaction(model=model1,coef=c(1,2,3),param = c("dummy"),conf.level=0.95) 
RERI1  

#$reri
#est     lower  upper
#1 0.557445 0.0183903 1.0965

#$apab
#est      lower     upper
#1 0.2474654 0.01795227 0.4769786

#$s
#est     lower    upper
#1 1.80188 0.8566727 3.789978

#$multiplicative
#est   lower    upper
#1 2.252618 1.83798 2.760795


RERI2=epi.interaction(model=model1,coef=c(1,4,5),param = c("dummy"),conf.level=0.95)
RERI2

#$reri
#est       lower    upper
#1 0.7289213 -0.03060201 1.488445

#$apab
#est      lower     upper
#1 0.2515623 0.02177738 0.4813473

#$s
#est     lower    upper
#1 1.623726 0.9448159 2.790476

#$multiplicative
#est    lower    upper
#1 2.897577 2.248271 3.734405



###森林图
setwd("Sfigure 6 & 12/中间文件/PRS")

   HR <- c(1.0,1.32,NA,1.38,2.25,NA,1.85,2.90)
lower <- c(1.0,0.92,NA,1.16,1.84,NA,1.53,2.25)
upper <- c(1.0,1.88,NA,1.63,2.76,NA,2.23,3.73)

foredata <- data.frame(HR,lower,upper)
foredata
class(foredata)
write.csv(foredata,file = "Sensitivity_analysis_excluding_participants_with_less_than_two_years_follow_up_PRS_FH_foredata_lung_female.csv",row.names = F)

pdata <- read.csv('Sensitivity_analysis_excluding_participants_with_less_than_two_years_follow_up_PRS_FH_foredata_lung_female.csv',header=T)
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
topptx(forestplot,"Sensitivity_analysis_excluding_participants_with_less_than_two_years_follow_up_PRS_FH_foredata_lung_female.pptx",width = 5,height = 4)




#前列腺癌

#男性
###剔除随访时间小于两年
dat1<-subset(male,male$PRC_difftime_new>2) #余下200347人，即随访时间小于两年的男性有202801-200347=2454人

#分层分析
control <- dat1[which(dat1$FH_Prostate_2==0),]
summary(control$PRO_PRS)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#5.459   7.734   8.197   8.211   8.672  11.766 

sd(control$PRO_PRS)  ##0.6943598
dat1$PRO_PRS_sd<-dat1$PRO_PRS/0.6943598

Q20=quantile(dat1$PRO_PRS_sd,seq(0.05,1,0.05))[[4]]
Q80=quantile(dat1$PRO_PRS_sd,seq(0.05,1,0.05))[[16]]

####对原始数据的PRS三分位，<=Q20，1；>=Q80，3
dat1$PRS_8_sd_3_Q20=0
dat1[dat1$PRO_PRS_sd<=Q20,]$PRS_8_sd_3_Q20=1
dat1[dat1$PRO_PRS_sd>Q20&dat1$PRO_PRS_sd<Q80,]$PRS_8_sd_3_Q20=2
dat1[dat1$PRO_PRS_sd>=Q80,]$PRS_8_sd_3_Q20=3

table(dat1$PRS_8_sd_3_Q20)
#1      2      3 
#40070 120207  40070  

dat1$PRO_PRS_FH=0
dat1[dat1$PRS_8_sd_3_Q20==1&dat1$FH_Prostate_2==0,]$PRO_PRS_FH=1  #对照（1，2，3，4）（1，2，5，6）
dat1[dat1$PRS_8_sd_3_Q20==1&dat1$FH_Prostate_2==1,]$PRO_PRS_FH=2  #家族史效应
dat1[dat1$PRS_8_sd_3_Q20==2&dat1$FH_Prostate_2==0,]$PRO_PRS_FH=3  #PRS中位效应
dat1[dat1$PRS_8_sd_3_Q20==2&dat1$FH_Prostate_2==1,]$PRO_PRS_FH=4  #PRS中位联合家族史效应
dat1[dat1$PRS_8_sd_3_Q20==3&dat1$FH_Prostate_2==0,]$PRO_PRS_FH=5  #PRS高位效应
dat1[dat1$PRS_8_sd_3_Q20==3&dat1$FH_Prostate_2==1,]$PRO_PRS_FH=6  #PRS高位联合家族史效应
table(dat1$PRO_PRS_FH)
#1      2      3      4      5      6 
#37716   2354 111052   9155  36060   4010


#case/人年
fit2 <- pyears(Surv(PRC_difftime_new,PRC_total.y)~1,data=dat1[dat1$PRO_PRS_FH==1,],scale = 1)
fit2$pyears#
fit2$event#
round(567/407578*100000,2) #139.11

fit2 <- pyears(Surv(PRC_difftime_new,PRC_total.y)~1,data=dat1[dat1$PRO_PRS_FH==2,],scale = 1)
fit2$pyears#
fit2$event#  
round(46/25864*100000,2) #177.85

fit2 <- pyears(Surv(PRC_difftime_new,PRC_total.y)~1,data=dat1[dat1$PRO_PRS_FH==3,],scale = 1)
fit2$pyears#
fit2$event#   
round(3879/1193021*100000,2) #325.14

fit2 <- pyears(Surv(PRC_difftime_new,PRC_total.y)~1,data=dat1[dat1$PRO_PRS_FH==4,],scale = 1)
fit2$pyears#
fit2$event#  
round(540/98433*100000,2) #548.60

fit2 <- pyears(Surv(PRC_difftime_new,PRC_total.y)~1,data=dat1[dat1$PRO_PRS_FH==5,],scale = 1)
fit2$pyears#
fit2$event#   
round(2794/382899*100000,2) #729.70

fit2 <- pyears(Surv(PRC_difftime_new,PRC_total.y)~1,data=dat1[dat1$PRO_PRS_FH==6,],scale = 1)
fit2$pyears#
fit2$event#  
round(461/40173*100000,2) #1147.53


model1 <- coxph(Surv(PRC_difftime_new,PRC_total.y)~as.factor(PRO_PRS_FH)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=dat1)
res <- summary(model1)
print(res$coefficients,digits=2)
print(res$conf.int,digits=2)
#                                             coef exp(coef) se(coef)      z Pr(>|z|)
#as.factor(PRO_PRS_FH)2               0.19543      1.22  0.15331  1.275  2.0e-01
#as.factor(PRO_PRS_FH)3               0.86458      2.37  0.04498 19.223  2.4e-82
#as.factor(PRO_PRS_FH)4               1.35548      3.88  0.06016 22.532 2.0e-112
#as.factor(PRO_PRS_FH)5               1.71176      5.54  0.04617 37.072 8.1e-301
#as.factor(PRO_PRS_FH)6               2.15066      8.59  0.06280 34.244 5.4e-257
#                                      exp(coef) exp(-coef) lower .95 upper .95
#as.factor(PRO_PRS_FH)2                   1.22       0.82      0.90      1.64
#as.factor(PRO_PRS_FH)3                   2.37       0.42      2.17      2.59
#as.factor(PRO_PRS_FH)4                   3.88       0.26      3.45      4.36
#as.factor(PRO_PRS_FH)5                   5.54       0.18      5.06      6.06
#as.factor(PRO_PRS_FH)6                   8.59       0.12      7.60      9.72


##  相加交互作用检验
#coef的2，3，4，5，6分别对应分组时的2-6组，但是却是生成系数的第1-5行

RERI1=epi.interaction(model=model1,coef=c(1,2,3),param = c("dummy"),conf.level=0.95) 
RERI1  
#$reri
#est    lower    upper
#1 1.288769 0.802234 1.775303

#$apab
#est     lower     upper
#1 0.3322753 0.2218146 0.4427361

#$s
#est   lower    upper
#1 1.810624 1.39079 2.357192

#$multiplicative
#est   lower    upper
#1 3.878616 3.44723 4.363986


RERI2=epi.interaction(model=model1,coef=c(1,4,5),param = c("dummy"),conf.level=0.95)
RERI2
#$reri
#est    lower    upper
#1 2.835972 1.938997 3.732946

#$apab
#est    lower     upper
#1 0.3301286 0.252145 0.4081121

#$s
#est    lower    upper
#1 1.596477 1.393026 1.829642

#$multiplicative
#est    lower    upper
#1 8.590507 7.595552 9.715793


###森林图
setwd("Sfigure 6 & 12/中间文件/PRS")

   HR <- c(1.0,1.22,NA,2.37,3.88,NA,5.54,8.59)
lower <- c(1.0,0.90,NA,2.17,3.45,NA,5.06,7.60)
upper <- c(1.0,1.64,NA,2.59,4.36,NA,6.06,9.72)

foredata <- data.frame(HR,lower,upper)
foredata
class(foredata)
write.csv(foredata,file = "Sensitivity_analysis_excluding_participants_with_less_than_two_years_follow_up_PRS_FH_foredata_prostate_male.csv",row.names = F)

pdata <- read.csv('Sensitivity_analysis_excluding_participants_with_less_than_two_years_follow_up_PRS_FH_foredata_prostate_male.csv',header=T)
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
topptx(forestplot,"Sensitivity_analysis_excluding_participants_with_less_than_two_years_follow_up_PRS_FH_foredata_prostate_male.pptx",width = 5,height = 4)



#女性
###剔除随访时间小于两年
dat1<-subset(female,female$BC_difftime_new>2) #余下237559人，即随访时间小于两年女性有239598-237559=2039人

control <- dat1[which(dat1$FH_Breast_2==0),]
summary(control$BC_PRS)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-0.74588 -0.18938 -0.10454 -0.10553 -0.02081  0.50938

sd(control$BC_PRS)  ##0.1271452
dat1$BC_PRS_sd<-dat1$BC_PRS/0.1271452

Q20=quantile(dat1$BC_PRS_sd,seq(0.05,1,0.05))[[4]]
Q80=quantile(dat1$BC_PRS_sd,seq(0.05,1,0.05))[[16]]
####对原始数据的PRS三分位，<=Q20，1；>=Q80，3

dat1$PRS_28_sd_3_Q20=0
dat1[dat1$BC_PRS_sd<=Q20,]$PRS_28_sd_3_Q20=1
dat1[dat1$BC_PRS_sd>Q20&dat1$BC_PRS_sd<Q80,]$PRS_28_sd_3_Q20=2
dat1[dat1$BC_PRS_sd>=Q80,]$PRS_28_sd_3_Q20=3

table(dat1$PRS_28_sd_3_Q20)
#1      2      3 
#47512 142535  47512  

dat1$BC_PRS_FH=0
dat1[dat1$PRS_28_sd_3_Q20==1&dat1$FH_Breast_2==0,]$BC_PRS_FH=1  #对照（1，2，3，4）（1，2，5，6）
dat1[dat1$PRS_28_sd_3_Q20==1&dat1$FH_Breast_2==1,]$BC_PRS_FH=2  #家族史效应
dat1[dat1$PRS_28_sd_3_Q20==2&dat1$FH_Breast_2==0,]$BC_PRS_FH=3  #PRS中位效应
dat1[dat1$PRS_28_sd_3_Q20==2&dat1$FH_Breast_2==1,]$BC_PRS_FH=4  #PRS中位联合家族史效应
dat1[dat1$PRS_28_sd_3_Q20==3&dat1$FH_Breast_2==0,]$BC_PRS_FH=5  #PRS高位效应
dat1[dat1$PRS_28_sd_3_Q20==3&dat1$FH_Breast_2==1,]$BC_PRS_FH=6  #PRS高位联合家族史效应
table(dat1$BC_PRS_FH)
#1      2      3      4      5      6 
#43215   4297 127093  15442  41332   6180 


#case/人年
fit2 <- pyears(Surv(BC_difftime_new,BC_total)~1,data=dat1[dat1$BC_PRS_FH==1,],scale = 1)
fit2$pyears#
fit2$event#
round(718/473230*100000,2) #151.72

fit2 <- pyears(Surv(BC_difftime_new,BC_total)~1,data=dat1[dat1$BC_PRS_FH==2,],scale = 1)
fit2$pyears#
fit2$event#  
round(121/46725*100000,2) #258.96

fit2 <- pyears(Surv(BC_difftime_new,BC_total)~1,data=dat1[dat1$BC_PRS_FH==3,],scale = 1)
fit2$pyears#
fit2$event#   
round(3321/1385043*100000,2) #239.78

fit2 <- pyears(Surv(BC_difftime_new,BC_total)~1,data=dat1[dat1$BC_PRS_FH==4,],scale = 1)
fit2$pyears#
fit2$event#  
round(611/167122*100000,2) #365.60

fit2 <- pyears(Surv(BC_difftime_new,BC_total)~1,data=dat1[dat1$BC_PRS_FH==5,],scale = 1)
fit2$pyears#
fit2$event#   
round(1655/446283*100000,2) #370.84

fit2 <- pyears(Surv(BC_difftime_new,BC_total)~1,data=dat1[dat1$BC_PRS_FH==6,],scale = 1)
fit2$pyears#
fit2$event#  
round(354/66382*100000,2) #533.28


model1 <- coxph(Surv(BC_difftime_new,BC_total)~as.factor(BC_PRS_FH)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=dat1)
res <- summary(model1)
print(res$coefficients,digits = 2)
print(res$conf.int,digits = 2)
#                                             coef exp(coef) se(coef)     z Pr(>|z|)
#as.factor(BC_PRS_FH)2                0.51788      1.68  0.09829  5.27  1.4e-07
#as.factor(BC_PRS_FH)3                0.46922      1.60  0.04118 11.39  4.4e-30
#as.factor(BC_PRS_FH)4                0.86914      2.38  0.05506 15.79  3.9e-56
#as.factor(BC_PRS_FH)5                0.92332      2.52  0.04483 20.59  3.0e-94
#as.factor(BC_PRS_FH)6                1.25847      3.52  0.06497 19.37  1.4e-83
#                                      exp(coef) exp(-coef) lower .95 upper .95
#as.factor(BC_PRS_FH)2                    1.68       0.60      1.38      2.04
#as.factor(BC_PRS_FH)3                    1.60       0.63      1.47      1.73
#as.factor(BC_PRS_FH)4                    2.38       0.42      2.14      2.66
#as.factor(BC_PRS_FH)5                    2.52       0.40      2.31      2.75
#as.factor(BC_PRS_FH)6                    3.52       0.28      3.10      4.00


##  相加交互作用检验
#coef的2，3，4，5，6分别对应分组时的2-6组，但是却是生成系数的第1-5行

RERI1=epi.interaction(model=model1,coef=c(1,2,3),param = c("dummy"),conf.level=0.95) 
RERI1  

#$reri
#est      lower     upper
#1 0.1076362 -0.2562721 0.4715444

#$apab
#est      lower   upper
#1 0.04513332 -0.1062634 0.19653

#$s
#est     lower    upper
#1 1.084274 0.8184508 1.436434

#$multiplicative
#est    lower    upper
#1 2.38485 2.140896 2.656603



RERI2=epi.interaction(model=model1,coef=c(1,4,5),param = c("dummy"),conf.level=0.95)
RERI2

#$reri
#est     lower     upper
#1 0.3239666 -0.167099 0.8150323

#$apab
#est       lower     upper
#1 0.0920347 -0.04134844 0.2254178

#$s
#est     lower    upper
#1 1.14752 0.9299375 1.416012

#$multiplicative
#est    lower    upper
#1 3.520049 3.099149 3.998111



###森林图
setwd("Sfigure 6 & 12/中间文件/PRS")

   HR <- c(1.0,1.68,NA,1.60,2.38,NA,2.52,3.52)
lower <- c(1.0,1.38,NA,1.47,2.14,NA,2.31,3.10)
upper <- c(1.0,2.04,NA,1.73,2.66,NA,2.75,4.00)

foredata <- data.frame(HR,lower,upper)
foredata
class(foredata)
write.csv(foredata,file = "Sensitivity_analysis_excluding_participants_with_less_than_two_years_follow_up_PRS_FH_foredata_breast_female.csv",row.names = F)

pdata <- read.csv('Sensitivity_analysis_excluding_participants_with_less_than_two_years_follow_up_PRS_FH_foredata_breast_female.csv',header=T)
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
topptx(forestplot,"Sensitivity_analysis_excluding_participants_with_less_than_two_years_follow_up_PRS_FH_foredata_breast_female.pptx",width = 5,height = 4)






####Sfigure12


#多肿瘤家族史及CPRS
#男性
###剔除随访时间小于两年
dat1<-subset(male,male$survival_time>2) #余下199185人，即随访时间小于两年男性有202801-199185=3616人

control <- dat1[which(dat1$FH_total_2==0),]
summary(control$CPRS)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#2.878   3.602   3.735   3.737   3.869   4.801 
sd(control$CPRS)  ##0.1986839
dat1$CPRS_sd<-dat1$CPRS/0.1986839


Q20=quantile(dat1$CPRS_sd,seq(0.05,1,0.05))[[4]]
Q80=quantile(dat1$CPRS_sd,seq(0.05,1,0.05))[[16]]

####对原始数据的PRS三分位，<=Q20，1；>=Q80，3
dat1$CPRS_sd_3_Q20=0
dat1[dat1$CPRS_sd<=Q20,]$CPRS_sd_3_Q20=1
dat1[dat1$CPRS_sd>Q20&dat1$CPRS_sd<Q80,]$CPRS_sd_3_Q20=2
dat1[dat1$CPRS_sd>=Q80,]$CPRS_sd_3_Q20=3

table(dat1$CPRS_sd_3_Q20)
#1      2      3 
#39837 119511  39837 

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
#29929  9908 86478 33033 27663 12174


#case/人年
fit2 <- pyears(Surv(survival_time,cancer_total_20)~1,data=dat1[dat1$CPRS_FH==1,],scale = 1)
fit2$pyears#
fit2$event#
round(1591/321025*100000,2) #495.60

fit2 <- pyears(Surv(survival_time,cancer_total_20)~1,data=dat1[dat1$CPRS_FH==2,],scale = 1)
fit2$pyears#
fit2$event#  
round(669/105677*100000,2) #633.06

fit2 <- pyears(Surv(survival_time,cancer_total_20)~1,data=dat1[dat1$CPRS_FH==3,],scale = 1)
fit2$pyears#
fit2$event#   
round(6555/920853*100000,2) #711.84

fit2 <- pyears(Surv(survival_time,cancer_total_20)~1,data=dat1[dat1$CPRS_FH==4,],scale = 1)
fit2$pyears#
fit2$event#  
round(3034/349932*100000,2) #867.03

fit2 <- pyears(Surv(survival_time,cancer_total_20)~1,data=dat1[dat1$CPRS_FH==5,],scale = 1)
fit2$pyears#
fit2$event#   
round(3259/289975*100000,2) #1123.89

fit2 <- pyears(Surv(survival_time,cancer_total_20)~1,data=dat1[dat1$CPRS_FH==6,],scale = 1)
fit2$pyears#
fit2$event#  
round(1759/126397*100000,2) #1391.64


model1 <- coxph(Surv(survival_time,cancer_total_20)~as.factor(CPRS_FH)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=dat1)
res <- summary(model1)
print(res$coefficients,digit=2)
print(res$conf.int,digit=2)
#                                         coef exp(coef) se(coef)      z Pr(>|z|)
#as.factor(CPRS_FH)2                  0.11948      1.13  0.04613   2.59  9.6e-03
#as.factor(CPRS_FH)3                  0.35886      1.43  0.02803  12.80  1.6e-37
#as.factor(CPRS_FH)4                  0.46150      1.59  0.03111  14.84  8.5e-50
#as.factor(CPRS_FH)5                  0.83987      2.32  0.03073  27.33 1.9e-164
#as.factor(CPRS_FH)6                  0.97211      2.64  0.03478  27.95 6.5e-172

#                                     exp(coef) exp(-coef) lower .95 upper .95
#as.factor(CPRS_FH)2                      1.13       0.89      1.03      1.23
#as.factor(CPRS_FH)3                      1.43       0.70      1.36      1.51
#as.factor(CPRS_FH)4                      1.59       0.63      1.49      1.69
#as.factor(CPRS_FH)5                      2.32       0.43      2.18      2.46
#as.factor(CPRS_FH)6                      2.64       0.38      2.47      2.83

#计算家族史（无/有）和遗传风险（低/中/高）之间的相加交互作用
#设HR11表示2 个危险因素同时存在的HR值，HR10,HR01分别表示存在1个因素时HR值。相加交互作用的评价指标：RERI=HR1-HR01-HR10+1,AP=RERI/HR11

##  相加交互作用检验
#coef的2，3，4，5，6分别对应分组时的2-6组，但是却是生成系数的第1-5行

RERI1=epi.interaction(model=model1,coef=c(1,2,3),param = c("dummy"),conf.level=0.95) 
RERI1  

#$reri
#est       lower     upper
#1 0.02784417 -0.09032784 0.1460162

#$apab
#est       lower      upper
#1 0.01755124 -0.05699508 0.09209757

#$s
#est     lower    upper
#1 1.049846 0.8492085 1.297886

#$multiplicative
#est   lower    upper
#1 1.58645 1.49262 1.686179


RERI2=epi.interaction(model=model1,coef=c(1,4,5),param = c("dummy"),conf.level=0.95)
RERI2
#$reri
#est      lower     upper
#1 0.2005499 0.02610607 0.3749937

#$apab
#est      lower    upper
#1 0.07586477 0.01136558 0.140364

#$s
#est    lower    upper
#1 1.138984 1.013379 1.280158

#$multiplicative
#est    lower    upper
#1 2.643518 2.469323 2.830002



###森林图
setwd("Sfigure 6 & 12/中间文件/CPRS")

   HR <- c(1.0,1.13,NA,1.43,1.59,NA,2.32,2.64)
lower <- c(1.0,1.03,NA,1.36,1.49,NA,2.18,2.47)
upper <- c(1.0,1.23,NA,1.51,1.69,NA,2.46,2.83)

foredata <- data.frame(HR,lower,upper)
foredata
class(foredata)
write.csv(foredata,file = "Sensitivity_analysis_excluding_participants_with_less_than_two_years_follow_up_CPRS_FH_foredata_overall_male.csv",row.names = F)

pdata <- read.csv('Sensitivity_analysis_excluding_participants_with_less_than_two_years_follow_up_CPRS_FH_foredata_overall_male.csv',header=T)
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
topptx(forestplot,"Sensitivity_analysis_excluding_participants_with_less_than_two_years_follow_up_CPRS_FH_foredata_overall_male.pptx",width = 5,height = 4)



#女性
###剔除随访时间小于两年
dat1<-subset(female,female$survival_time>2) #余下236301人，即随访时间小于两年女性有239598-236301=3297人

#分层分析
control <- dat1[which(dat1$FH_total_2==0),]
summary(control$CPRS)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.070   1.385   1.438   1.437   1.489   1.825

sd(control$CPRS)  ##0.0782152
dat1$CPRS_sd<-dat1$CPRS/0.0782152

Q20=quantile(dat1$CPRS_sd,seq(0.05,1,0.05))[[4]]
Q80=quantile(dat1$CPRS_sd,seq(0.05,1,0.05))[[16]]

####对原始数据的PRS三分位，<=Q20，1；>=Q80，3
dat1$CPRS_sd_3_Q20=0
dat1[dat1$CPRS_sd<=Q20,]$CPRS_sd_3_Q20=1
dat1[dat1$CPRS_sd>Q20&dat1$CPRS_sd<Q80,]$CPRS_sd_3_Q20=2
dat1[dat1$CPRS_sd>=Q80,]$CPRS_sd_3_Q20=3

table(dat1$CPRS_sd_3_Q20)
#1      2      3 
#47261 141779  47261  

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
#35084 12177 98638 43141 31646 15615   


#case/人年
fit2 <- pyears(Surv(survival_time,cancer_total_20)~1,data=dat1[dat1$CPRS_FH==1,],scale = 1)
fit2$pyears#
fit2$event#
round(1656/378075*100000,2) #438.01

fit2 <- pyears(Surv(survival_time,cancer_total_20)~1,data=dat1[dat1$CPRS_FH==2,],scale = 1)
fit2$pyears#
fit2$event#  
round(724/131305*100000,2) #551.39

fit2 <- pyears(Surv(survival_time,cancer_total_20)~1,data=dat1[dat1$CPRS_FH==3,],scale = 1)
fit2$pyears#
fit2$event#   
round(5827/1066267*100000,2) #546.49

fit2 <- pyears(Surv(survival_time,cancer_total_20)~1,data=dat1[dat1$CPRS_FH==4,],scale = 1)
fit2$pyears#
fit2$event#  
round(3133/464043*100000,2) #675.15

fit2 <- pyears(Surv(survival_time,cancer_total_20)~1,data=dat1[dat1$CPRS_FH==5,],scale = 1)
fit2$pyears#
fit2$event#   
round(2385/339974*100000,2) #701.52

fit2 <- pyears(Surv(survival_time,cancer_total_20)~1,data=dat1[dat1$CPRS_FH==6,],scale = 1)
fit2$pyears#
fit2$event#  
round(1409/166955*100000,2) #843.94


model1 <- coxph(Surv(survival_time,cancer_total_20)~as.factor(CPRS_FH)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=dat1)
res <- summary(model1)
print(res$coefficients,digit=2)
print(res$conf.int,digit=2)

#                                          coef exp(coef) se(coef)      z Pr(>|z|)
#as.factor(CPRS_FH)2                  1.6e-01      1.18  0.04482   3.64  2.7e-04
#as.factor(CPRS_FH)3                  2.1e-01      1.24  0.02887   7.33  2.2e-13
#as.factor(CPRS_FH)4                  3.6e-01      1.44  0.03152  11.56  6.6e-31
#as.factor(CPRS_FH)5                  4.7e-01      1.59  0.03314  14.09  4.5e-45
#as.factor(CPRS_FH)6                  6.0e-01      1.82  0.03731  16.01  1.2e-57

#                                   exp(coef) exp(-coef) lower .95 upper .95
#as.factor(CPRS_FH)2                      1.18       0.85      1.08      1.29
#as.factor(CPRS_FH)3                      1.24       0.81      1.17      1.31
#as.factor(CPRS_FH)4                      1.44       0.69      1.35      1.53
#as.factor(CPRS_FH)5                      1.59       0.63      1.49      1.70
#as.factor(CPRS_FH)6                      1.82       0.55      1.69      1.95


#计算家族史（有/无）和遗传风险（低/中/高）之间的相加交互作用
#设HR11表示2 个危险因素同时存在的HR值，HR10,HR01分别表示存在1个因素时HR值。相加交互作用的评价指标：RERI=HR1-HR01-HR10+1,AP=RERI/HR11

##  相加交互作用检验
#coef的2，3，4，5，6分别对应分组时的2-6组，但是却是生成系数的第1-5行

RERI1=epi.interaction(model=model1,coef=c(1,2,3),param = c("dummy"),conf.level=0.95) 
RERI1  

#$reri
#est       lower     upper
#1 0.02636144 -0.08825972 0.1409826

#$apab
#est      lower      upper
#1 0.01831163 -0.0614088 0.09803206

#$s
#est     lower    upper
#1 1.063792 0.8050243 1.405739

#$multiplicative
#est   lower    upper
#1 1.439601 1.35336 1.531338


RERI2=epi.interaction(model=model1,coef=c(1,4,5),param = c("dummy"),conf.level=0.95)
RERI2
#$reri
#est      lower     upper
#1 0.044586 -0.1058191 0.1949911

#$apab
#est       lower     upper
#1 0.02453976 -0.05782991 0.1069094

#$s
#est     lower    upper
#1 1.057731 0.8723313 1.282535

#$multiplicative
#est    lower    upper
#1 1.816888 1.688772 1.954723

###森林图
setwd("Sfigure 6 & 12/中间文件/CPRS")


   HR <- c(1.0,1.18,NA,1.24,1.44,NA,1.59,1.82)
lower <- c(1.0,1.08,NA,1.17,1.35,NA,1.49,1.69)
upper <- c(1.0,1.29,NA,1.31,1.53,NA,1.70,1.95)

foredata <- data.frame(HR,lower,upper)
foredata
class(foredata)
write.csv(foredata,file = "Sensitivity_analysis_excluding_participants_with_less_than_two_years_follow_up_CPRS_FH_foredata_overall_female.csv",row.names = F)

pdata <- read.csv('Sensitivity_analysis_excluding_participants_with_less_than_two_years_follow_up_CPRS_FH_foredata_overall_female.csv',header=T)
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
topptx(forestplot,"Sensitivity_analysis_excluding_participants_with_less_than_two_years_follow_up_CPRS_FH_foredata_overall_female.pptx",width = 5,height = 4)

