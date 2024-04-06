####Sfigure 9 & 15####

#敏感性分析-局限于Unrelated white British population
#first-,second-, and third-degree relatives and population other than White British were excluded (359,941 individuals retained)


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
library(eoffice) 

#定义亲缘关系和种族
relate=read.table('genetic_relatedness.tab',h=T)
names(relate)[1]='f.eid'
names(relate)[2]='ethnic_background'
names(relate)[3]='genetic_sex'
names(relate)[4]='recommend_exclusion'
names(relate)[5]='genetic_relate'
names(relate)[6]='kinship'

mer=merge(male,relate,by.x=c("f.eid"),by.y=c("f.eid"))
male1=subset(mer,mer$kinship==0&mer$ethnic_background=='1001')  #剩余123649人，排除202801-123649=79152人
 
mer=merge(female,relate,by.x=c("f.eid"),by.y=c("f.eid"))
female1=subset(mer,mer$kinship==0&mer$ethnic_background=='1001') #剩余141133， 排除239598-141133=98465人



####Sfigure 9

########肠癌 CRC  COL_PRS

#male  COL_PRS
control <- male1[which(male1$FH_Bowel_2==0),]
summary(control$COL_PRS)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.718   2.922   3.178   3.190   3.446   5.035

sd(control$COL_PRS)  ##0.389234
male1$COL_PRS_sd<-male1$COL_PRS/0.389234

Q20=quantile(male1$COL_PRS_sd,seq(0.05,1,0.05))[[4]]
Q80=quantile(male1$COL_PRS_sd,seq(0.05,1,0.05))[[16]]

####对原始数据的PRS三分位，<=Q20，1；>=Q80，3
male1$PRS_1_sd_3_Q20=0
male1[male1$COL_PRS_sd<=Q20,]$PRS_1_sd_3_Q20=1
male1[male1$COL_PRS_sd>Q20&male1$COL_PRS_sd<Q80,]$PRS_1_sd_3_Q20=2
male1[male1$COL_PRS_sd>=Q80,]$PRS_1_sd_3_Q20=3

table(male1$PRS_1_sd_3_Q20)
#1      2      3 
#24730 74189 24730

#生成新变量COL_PRS_FH，组合PRS及家族史将人群分为6层
male1$COL_PRS_FH=0
male1[male1$PRS_1_sd_3_Q20==1&male1$FH_Bowel_2==0,]$COL_PRS_FH=1  #对照（1，2，3，4）（1，2，5，6）
male1[male1$PRS_1_sd_3_Q20==1&male1$FH_Bowel_2==1,]$COL_PRS_FH=2  #家族史效应
male1[male1$PRS_1_sd_3_Q20==2&male1$FH_Bowel_2==0,]$COL_PRS_FH=3  #PRS中位效应
male1[male1$PRS_1_sd_3_Q20==2&male1$FH_Bowel_2==1,]$COL_PRS_FH=4  #PRS中位联合家族史效应
male1[male1$PRS_1_sd_3_Q20==3&male1$FH_Bowel_2==0,]$COL_PRS_FH=5  #PRS高位效应
male1[male1$PRS_1_sd_3_Q20==3&male1$FH_Bowel_2==1,]$COL_PRS_FH=6  #PRS高位联合家族史效应
table(male1$COL_PRS_FH)
#1      2      3      4      5      6 
#22370  2360 65833  8356 21406  3324  


#case/人年
fit2 <- pyears(Surv(CRC_difftime_new,CRC_total.y)~1,data=male1[male1$COL_PRS_FH==1,],scale = 1)
fit2$pyears#
fit2$event#
round(148/242128*100000,2) #61.12

fit2 <- pyears(Surv(CRC_difftime_new,CRC_total.y)~1,data=male1[male1$COL_PRS_FH==2,],scale = 1)
fit2$pyears#
fit2$event#  
round(20/25389*100000,2) #78.78

fit2 <- pyears(Surv(CRC_difftime_new,CRC_total.y)~1,data=male1[male1$COL_PRS_FH==3,],scale = 1)
fit2$pyears#
fit2$event#  
round(868/710166*100000,2) #122.23

fit2 <- pyears(Surv(CRC_difftime_new,CRC_total.y)~1,data=male1[male1$COL_PRS_FH==4,],scale = 1)
fit2$pyears#
fit2$event#   
round(157/89740*100000,2) #174.95

fit2 <- pyears(Surv(CRC_difftime_new,CRC_total.y)~1,data=male1[male1$COL_PRS_FH==5,],scale = 1)
fit2$pyears#
fit2$event#   
round(477/229465*100000,2) #207.87

fit2 <- pyears(Surv(CRC_difftime_new,CRC_total.y)~1,data=male1[male1$COL_PRS_FH==6,],scale = 1)
fit2$pyears#
fit2$event#   
round(89/35562*100000,2) #250.27


model1 <- coxph(Surv(CRC_difftime_new,CRC_total.y)~as.factor(COL_PRS_FH)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male1)
res <- summary(model1)
print(res$coefficients,digits = 2)
print(res$conf.int,digits = 2)
#                                             coef exp(coef) se(coef)      z Pr(>|z|)
#as.factor(COL_PRS_FH)2               0.14765      1.16   0.2383  0.620  5.4e-01
#as.factor(COL_PRS_FH)3               0.69610      2.01   0.0889  7.826  5.0e-15
#as.factor(COL_PRS_FH)4               0.95058      2.59   0.1146  8.292  1.1e-16
#as.factor(COL_PRS_FH)5               1.23568      3.44   0.0941 13.127  2.3e-39
#as.factor(COL_PRS_FH)6               1.33868      3.81   0.1342  9.974  2.0e-23

#                                      exp(coef) exp(-coef) lower .95 upper .95
#as.factor(COL_PRS_FH)2                   1.16       0.86      0.73      1.85
#as.factor(COL_PRS_FH)3                   2.01       0.50      1.69      2.39
#as.factor(COL_PRS_FH)4                   2.59       0.39      2.07      3.24
#as.factor(COL_PRS_FH)5                   3.44       0.29      2.86      4.14
#as.factor(COL_PRS_FH)6                   3.81       0.26      2.93      4.96


#计算家族史（无/有）和遗传风险（低/中/高）之间的相加交互作用
#设HR11表示2 个危险因素同时存在的HR值，HR10,HR01分别表示存在1个因素时HR值。相加交互作用的评价指标：RERI=HR1-HR01-HR10+1,AP=RERI/HR11
##  相加交互作用检验
#coef的2，3，4，5，6分别对应分组时的2-6组，但是却是生成系数的第1-5行


RERI1=epi.interaction(model=model1,coef=c(1,2,3),param = c("dummy"),conf.level=0.95) 
RERI1  
#$reri
#est      lower    upper
#1 0.422198 -0.2474631 1.091859

#$apab
#est      lower     upper
#1 0.1631861 -0.0862862 0.4126585

#$s
#est     lower    upper
#1 1.362396 0.7885782 2.353758

#$multiplicative
#est    lower    upper
#1 2.587218 2.066594 3.238998


RERI2=epi.interaction(model=model1,coef=c(1,4,5),param = c("dummy"),conf.level=0.95)
RERI2
#$reri
#est      lower    upper
#1 0.2142026 -0.7847094 1.213115

#$apab
#est      lower     upper
#1 0.05616185 -0.1979246 0.3102483

#$s
#est     lower    upper
#1 1.082391 0.7484692 1.565289

#$multiplicative
#est    lower    upper
#1 3.814023 2.931795 4.961727



###森林图
setwd("Sfigure 9 & 15/中间文件/PRS")

   HR <- c(1.0,1.16,NA,2.01,2.59,NA,3.44,3.81)
lower <- c(1.0,0.73,NA,1.69,2.07,NA,2.86,2.93)
upper <- c(1.0,1.85,NA,2.39,3.24,NA,4.14,4.96)

foredata <- data.frame(HR,lower,upper)
foredata
class(foredata)
write.csv(foredata,file = "Sensitivity_analysis_british_population_PRS_FH_foredata_colorectal_male.csv",row.names = F)

pdata <- read.csv('Sensitivity_analysis_british_population_PRS_FH_foredata_colorectal_male.csv',header=T)
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

topptx(forestplot,"Sensitivity_analysis_british_population_forestplot_colorectal_male_FH_PRS.pptx",width = 5,height = 4)



###colorectal-female###
control <- female1[which(female1$FH_Bowel_2==0),]
summary(control$COL_PRS)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.595   2.922   3.178   3.191   3.448   5.200

sd(control$COL_PRS)  ##0.3902893
female1$COL_PRS_sd<-female1$COL_PRS/0.3902893

Q20=quantile(female1$COL_PRS_sd,seq(0.05,1,0.05))[[4]]
Q80=quantile(female1$COL_PRS_sd,seq(0.05,1,0.05))[[16]]

####对原始数据的PRS三分位，<=Q20，1；>=Q80，3
female1$PRS_1_sd_3_Q20=0
female1[female1$COL_PRS_sd<=Q20,]$PRS_1_sd_3_Q20=1
female1[female1$COL_PRS_sd>Q20&female1$COL_PRS_sd<Q80,]$PRS_1_sd_3_Q20=2
female1[female1$COL_PRS_sd>=Q80,]$PRS_1_sd_3_Q20=3

table(female1$PRS_1_sd_3_Q20)
#1      2      3 
#28227 84679 28227   

female1$COL_PRS_FH=0
female1[female1$PRS_1_sd_3_Q20==1&female1$FH_Bowel_2==0,]$COL_PRS_FH=1  #对照（1，2，3，4）（1，2，5，6）
female1[female1$PRS_1_sd_3_Q20==1&female1$FH_Bowel_2==1,]$COL_PRS_FH=2  #家族史效应
female1[female1$PRS_1_sd_3_Q20==2&female1$FH_Bowel_2==0,]$COL_PRS_FH=3  #PRS中位效应
female1[female1$PRS_1_sd_3_Q20==2&female1$FH_Bowel_2==1,]$COL_PRS_FH=4  #PRS中位联合家族史效应
female1[female1$PRS_1_sd_3_Q20==3&female1$FH_Bowel_2==0,]$COL_PRS_FH=5  #PRS高位效应
female1[female1$PRS_1_sd_3_Q20==3&female1$FH_Bowel_2==1,]$COL_PRS_FH=6  #PRS高位联合家族史效应
table(female1$COL_PRS_FH)
#1      2      3      4      5      6 
#25656  2571 75692  8987 24561  3666  


#case/人年
fit2 <- pyears(Surv(CRC_difftime_new,CRC_total.y)~1,data=female1[female1$COL_PRS_FH==1,],scale = 1)
fit2$pyears#
fit2$event#
round(114/281563*100000,2) #40.49

fit2 <- pyears(Surv(CRC_difftime_new,CRC_total.y)~1,data=female1[female1$COL_PRS_FH==2,],scale = 1)
fit2$pyears#
fit2$event#  
round(19/28150*100000,2) #67.50

fit2 <- pyears(Surv(CRC_difftime_new,CRC_total.y)~1,data=female1[female1$COL_PRS_FH==3,],scale = 1)
fit2$pyears#
fit2$event#   
round(637/829920*100000,2) #76.75

fit2 <- pyears(Surv(CRC_difftime_new,CRC_total.y)~1,data=female1[female1$COL_PRS_FH==4,],scale = 1)
fit2$pyears#
fit2$event#  
round(104/98465*100000,2) #105.62

fit2 <- pyears(Surv(CRC_difftime_new,CRC_total.y)~1,data=female1[female1$COL_PRS_FH==5,],scale = 1)
fit2$pyears#
fit2$event#  
round(362/268535*100000,2) #134.81

fit2 <- pyears(Surv(CRC_difftime_new,CRC_total.y)~1,data=female1[female1$COL_PRS_FH==6,],scale = 1)
fit2$pyears#
fit2$event#  
round(54/40133*100000,2) #134.55


model1 <- coxph(Surv(CRC_difftime_new,CRC_total.y)~as.factor(COL_PRS_FH)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female1)
res <- summary(model1)
print(res$coefficients,digits = 2)
print(res$conf.int,digits = 2)
#                                             coef exp(coef) se(coef)     z Pr(>|z|)
#as.factor(COL_PRS_FH)2               0.4260      1.53   0.2479  1.72  8.6e-02
#as.factor(COL_PRS_FH)3               0.6377      1.89   0.1017  6.27  3.6e-10
#as.factor(COL_PRS_FH)4               0.8762      2.40   0.1357  6.46  1.1e-10
#as.factor(COL_PRS_FH)5               1.2072      3.34   0.1074 11.24  2.7e-29
#as.factor(COL_PRS_FH)6               1.1356      3.11   0.1653  6.87  6.4e-12
#                                      exp(coef) exp(-coef) lower .95 upper .95
#as.factor(COL_PRS_FH)2                   1.53       0.65      0.94      2.49
#as.factor(COL_PRS_FH)3                   1.89       0.53      1.55      2.31
#as.factor(COL_PRS_FH)4                   2.40       0.42      1.84      3.13
#as.factor(COL_PRS_FH)5                   3.34       0.30      2.71      4.13
#as.factor(COL_PRS_FH)6                   3.11       0.32      2.25      4.30


##  相加交互作用检验
#coef的2，3，4，5，6分别对应分组时的2-6组，但是却是生成系数的第1-5行

RERI1=epi.interaction(model=model1,coef=c(1,2,3),param = c("dummy"),conf.level=0.95) 
RERI1  
#$reri
#est      lower     upper
#1 -0.02156388 -0.8841259 0.8409982

#$apab
#est      lower     upper
#1 -0.00897859 -0.3686985 0.3507413

#$s
#est     lower    upper
#1 0.984849 0.5363409 1.808416

#$multiplicative
#est    lower    upper
#1 2.4017 1.840867 3.133395


RERI2=epi.interaction(model=model1,coef=c(1,4,5),param = c("dummy"),conf.level=0.95)
RERI2
#reri
#est     lower     upper
#1 -0.7623274 -1.940043 0.4153879

#$apab
#est      lower    upper
#1 -0.2448818 -0.6631287 0.173365

#$s
#est     lower    upper
#1 0.7348767 0.4559511 1.184434

#$multiplicative
#est   lower    upper
#1 3.113042 2.25161 4.304044



###森林图
setwd("Sfigure 9 & 15/中间文件/PRS")

   HR <- c(1.0,1.53,NA,1.89,2.40,NA,3.34,3.11)
lower <- c(1.0,0.94,NA,1.55,1.84,NA,2.71,2.25)
upper <- c(1.0,2.49,NA,2.31,3.13,NA,4.13,4.30)

foredata <- data.frame(HR,lower,upper)
foredata
class(foredata)
write.csv(foredata,file = "Sensitivity_analysis_british_population_PRS_FH_foredata_colorectal_female.csv",row.names = F)

pdata <- read.csv('Sensitivity_analysis_british_population_PRS_FH_foredata_colorectal_female.csv',header=T)
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
topptx(forestplot,"Sensitivity_analysis_british_population_forestplot_colorectal_female_FH_PRS.pptx",width = 5,height = 4)




########肺癌 LC LC_PRS########

#male
control <- male1[which(male1$FH_Lung_2==0),]
summary(control$LC_PRS)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-0.37893  0.00933  0.06581  0.08419  0.14042  1.01989

sd(control$LC_PRS)  ##0.1094446
male1$LC_PRS_sd<-male1$LC_PRS/0.1094446

Q20=quantile(male1$LC_PRS_sd,seq(0.05,1,0.05))[[4]]
Q80=quantile(male1$LC_PRS_sd,seq(0.05,1,0.05))[[16]]

####对原始数据的PRS三分位，<=Q20，1；>=Q80，3
male1$PRS_7_sd_3_Q20=0
male1[male1$LC_PRS_sd<=Q20,]$PRS_7_sd_3_Q20=1
male1[male1$LC_PRS_sd>Q20&male1$LC_PRS_sd<Q80,]$PRS_7_sd_3_Q20=2
male1[male1$LC_PRS_sd>=Q80,]$PRS_7_sd_3_Q20=3

table(male1$PRS_7_sd_3_Q20)
#1      2      3 
#24730 74189 24730  

male1$LC_PRS_FH=0
male1[male1$PRS_7_sd_3_Q20==1&male1$FH_Lung_2==0,]$LC_PRS_FH=1  #对照（1，2，3，4）（1，2，5，6）
male1[male1$PRS_7_sd_3_Q20==1&male1$FH_Lung_2==1,]$LC_PRS_FH=2  #家族史效应
male1[male1$PRS_7_sd_3_Q20==2&male1$FH_Lung_2==0,]$LC_PRS_FH=3  #PRS中位效应
male1[male1$PRS_7_sd_3_Q20==2&male1$FH_Lung_2==1,]$LC_PRS_FH=4  #PRS中位联合家族史效应
male1[male1$PRS_7_sd_3_Q20==3&male1$FH_Lung_2==0,]$LC_PRS_FH=5  #PRS高位效应
male1[male1$PRS_7_sd_3_Q20==3&male1$FH_Lung_2==1,]$LC_PRS_FH=6  #PRS高位联合家族史效应
table(male1$LC_PRS_FH)
#1      2      3      4      5      6 
#21946  2784 65224  8965 21288  3442


#case/人年
fit2 <- pyears(Surv(LC_difftime_new,LC_total.y)~1,data=male1[male1$LC_PRS_FH==1,],scale = 1)
fit2$pyears#
fit2$event#
round(102/237694*100000,2) #42.91

fit2 <- pyears(Surv(LC_difftime_new,LC_total.y)~1,data=male1[male1$LC_PRS_FH==2,],scale = 1)
fit2$pyears#
fit2$event#  
round(23/30039*100000,2) #76.57

fit2 <- pyears(Surv(LC_difftime_new,LC_total.y)~1,data=male1[male1$LC_PRS_FH==3,],scale = 1)
fit2$pyears#
fit2$event#   
round(524/705998*100000,2) #74.22

fit2 <- pyears(Surv(LC_difftime_new,LC_total.y)~1,data=male1[male1$LC_PRS_FH==4,],scale = 1)
fit2$pyears#
fit2$event#  
round(138/96592*100000,2) #142.87

fit2 <- pyears(Surv(LC_difftime_new,LC_total.y)~1,data=male1[male1$LC_PRS_FH==5,],scale = 1)
fit2$pyears#
fit2$event#   
round(230/230455*100000,2) #99.80

fit2 <- pyears(Surv(LC_difftime_new,LC_total.y)~1,data=male1[male1$LC_PRS_FH==6,],scale = 1)
fit2$pyears#
fit2$event#  
round(62/37034*100000,2) #167.42


model1 <- coxph(Surv(LC_difftime_new,LC_total.y)~as.factor(LC_PRS_FH)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male1)
res <- summary(model1)
print(res$coefficients,digits = 2)
print(res$conf.int,digits = 2)
#                                             coef exp(coef) se(coef)       z Pr(>|z|)
#as.factor(LC_PRS_FH)2                0.38787      1.47   0.2310   1.679  9.3e-02
#as.factor(LC_PRS_FH)3                0.55769      1.75   0.1082   5.152  2.6e-07
#as.factor(LC_PRS_FH)4                1.05728      2.88   0.1307   8.089  6.0e-16
#as.factor(LC_PRS_FH)5                0.86506      2.38   0.1190   7.269  3.6e-13
#as.factor(LC_PRS_FH)6                1.21092      3.36   0.1612   7.513  5.8e-14

#                                      exp(coef) exp(-coef) lower .95 upper .95
#as.factor(LC_PRS_FH)2                    1.47       0.68      0.94      2.32
#as.factor(LC_PRS_FH)3                    1.75       0.57      1.41      2.16
#as.factor(LC_PRS_FH)4                    2.88       0.35      2.23      3.72
#as.factor(LC_PRS_FH)5                    2.38       0.42      1.88      3.00
#as.factor(LC_PRS_FH)6                    3.36       0.30      2.45      4.60

##  相加交互作用检验
#coef的2，3，4，5，6分别对应分组时的2-6组，但是却是生成系数的第1-5行


RERI1=epi.interaction(model=model1,coef=c(1,2,3),param = c("dummy"),conf.level=0.95) 
RERI1  
#$reri
#est      lower    upper
#1 0.6580518 -0.1296745 1.445778

#$apab
#est       lower     upper
#1 0.2286075 -0.03142378 0.4886388

#$s
#est     lower    upper
#1 1.539179 0.8376479 2.828243

#$multiplicative
#est    lower    upper
#1 2.878522 2.227984 3.719006

RERI2=epi.interaction(model=model1,coef=c(1,4,5),param = c("dummy"),conf.level=0.95)
RERI2
#$reri
#est    lower   upper
#1 0.5076051 -0.57172 1.58693

#$apab
#est      lower     upper
#1 0.1512266 -0.1464021 0.4488553

#$s
#est   lower   upper
#1 1.274532 0.75386 2.15482

#$multiplicative
#est    lower    upper
#1 3.356585 2.447421 4.603485


###森林图
setwd("Sfigure 9 & 15/中间文件/PRS")

   HR <- c(1.0,1.47,NA,1.75,2.88,NA,2.38,3.36)
lower <- c(1.0,0.94,NA,1.41,2.23,NA,1.88,2.45)
upper <- c(1.0,2.32,NA,2.16,3.72,NA,3.00,4.60)

foredata <- data.frame(HR,lower,upper)
foredata
class(foredata)
write.csv(foredata,file = "Sensitivity_analysis_british_population_PRS_FH_foredata_lung_male.csv",row.names = F)

pdata <- read.csv('Sensitivity_analysis_british_population_PRS_FH_foredata_lung_male.csv',header=T)
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
topptx(forestplot,"Sensitivity_analysis_british_population_forestplot_lung_male_FH_PRS.pptx",width = 5,height = 4)



###lung-female###
control <- female1[which(female1$FH_Lung_2==0),]
summary(control$LC_PRS)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-0.31775  0.01027  0.06676  0.08487  0.14145  1.19901 

sd(control$LC_PRS)  ##0.1088961
female1$LC_PRS_sd<-female1$LC_PRS/0.1088961

Q20=quantile(female1$LC_PRS_sd,seq(0.05,1,0.05))[[4]]
Q80=quantile(female1$LC_PRS_sd,seq(0.05,1,0.05))[[16]]

####对原始数据的PRS三分位，<=Q20，1；>=Q80，3
female1$PRS_7_sd_3_Q20=0
female1[female1$LC_PRS_sd<=Q20,]$PRS_7_sd_3_Q20=1
female1[female1$LC_PRS_sd>Q20&female1$LC_PRS_sd<Q80,]$PRS_7_sd_3_Q20=2
female1[female1$LC_PRS_sd>=Q80,]$PRS_7_sd_3_Q20=3

table(female1$PRS_7_sd_3_Q20)
#1      2      3 
#28227 84679 28227  

female1$LC_PRS_FH=0
female1[female1$PRS_7_sd_3_Q20==1&female1$FH_Lung_2==0,]$LC_PRS_FH=1  #对照（1，2，3，4）（1，2，5，6）
female1[female1$PRS_7_sd_3_Q20==1&female1$FH_Lung_2==1,]$LC_PRS_FH=2  #家族史效应
female1[female1$PRS_7_sd_3_Q20==2&female1$FH_Lung_2==0,]$LC_PRS_FH=3  #PRS中位效应
female1[female1$PRS_7_sd_3_Q20==2&female1$FH_Lung_2==1,]$LC_PRS_FH=4  #PRS中位联合家族史效应
female1[female1$PRS_7_sd_3_Q20==3&female1$FH_Lung_2==0,]$LC_PRS_FH=5  #PRS高位效应
female1[female1$PRS_7_sd_3_Q20==3&female1$FH_Lung_2==1,]$LC_PRS_FH=6  #PRS高位联合家族史效应
table(female1$LC_PRS_FH)
#1      2      3      4      5      6 
#25001  3226 73980 10699 24293  3934 


#case/人年
fit2 <- pyears(Surv(LC_difftime_new,LC_total.y)~1,data=female1[female1$LC_PRS_FH==1,],scale = 1)
fit2$pyears#
fit2$event#
round(107/274808*100000,2) #38.94

fit2 <- pyears(Surv(LC_difftime_new,LC_total.y)~1,data=female1[female1$LC_PRS_FH==2,],scale = 1)
fit2$pyears#
fit2$event#  
round(20/35599*100000,2) #56.18

fit2 <- pyears(Surv(LC_difftime_new,LC_total.y)~1,data=female1[female1$LC_PRS_FH==3,],scale = 1)
fit2$pyears#
fit2$event#   
round(468/812629*100000,2) #57.59

fit2 <- pyears(Surv(LC_difftime_new,LC_total.y)~1,data=female1[female1$LC_PRS_FH==4,],scale = 1)
fit2$pyears#
fit2$event#  
round(118/117418*100000,2) #100.50

fit2 <- pyears(Surv(LC_difftime_new,LC_total.y)~1,data=female1[female1$LC_PRS_FH==5,],scale = 1)
fit2$pyears#
fit2$event#   
round(194/266562*100000,2) #72.78

fit2 <- pyears(Surv(LC_difftime_new,LC_total.y)~1,data=female1[female1$LC_PRS_FH==6,],scale = 1)
fit2$pyears#
fit2$event#  
round(59/43154*100000,2) #136.72


model1 <- coxph(Surv(LC_difftime_new,LC_total.y)~as.factor(LC_PRS_FH)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female1)
res <- summary(model1)
print(res$coefficients,digits = 2)
print(res$conf.int,digits = 2)
#                                             coef exp(coef) se(coef)       z Pr(>|z|)
#as.factor(LC_PRS_FH)2                0.15937      1.17   0.2438   0.654  5.1e-01
#as.factor(LC_PRS_FH)3                0.40117      1.49   0.1072   3.743  1.8e-04
#as.factor(LC_PRS_FH)4                0.76089      2.14   0.1337   5.691  1.3e-08
#as.factor(LC_PRS_FH)5                0.65510      1.93   0.1205   5.438  5.4e-08
#as.factor(LC_PRS_FH)6                1.11287      3.04   0.1623   6.855  7.1e-12

#                                      exp(coef) exp(-coef) lower .95 upper .95
#as.factor(LC_PRS_FH)2                    1.17       0.85      0.73      1.89
#as.factor(LC_PRS_FH)3                    1.49       0.67      1.21      1.84
#as.factor(LC_PRS_FH)4                    2.14       0.47      1.65      2.78
#as.factor(LC_PRS_FH)5                    1.93       0.52      1.52      2.44
#as.factor(LC_PRS_FH)6                    3.04       0.33      2.21      4.18

##  相加交互作用检验
#coef的2，3，4，5，6分别对应分组时的2-6组，但是却是生成系数的第1-5行

RERI1=epi.interaction(model=model1,coef=c(1,2,3),param = c("dummy"),conf.level=0.95) 
RERI1  

#$reri
#est      lower   upper
#1 0.4738466 -0.1908269 1.13852

#$apab
#est      lower    upper
#1 0.2214044 -0.0773873 0.520196

#$s
#est     lower    upper
#1 1.711119 0.6580471 4.449419

#$multiplicative
#est    lower    upper
#1 2.140186 1.646804 2.781386



RERI2=epi.interaction(model=model1,coef=c(1,4,5),param = c("dummy"),conf.level=0.95)
RERI2

#$reri
#est       lower    upper
#1 0.9449536 -0.02541483 1.915322

#$apab
#est      lower     upper
#1 0.310527 0.04345424 0.5775997

#$s
#est     lower    upper
#1 1.860527 0.9250279 3.742113

#$multiplicative
#est   lower    upper
#1 3.043064 2.21376 4.183037



###森林图
setwd("Sfigure 9 & 15/中间文件/PRS")

   HR <- c(1.0,1.17,NA,1.49,2.14,NA,1.93,3.04)
lower <- c(1.0,0.73,NA,1.21,1.65,NA,1.52,2.21)
upper <- c(1.0,1.89,NA,1.84,2.78,NA,2.44,4.18)

foredata <- data.frame(HR,lower,upper)
foredata
class(foredata)
write.csv(foredata,file = "Sensitivity_analysis_british_population_PRS_FH_foredata_lung_female.csv",row.names = F)

pdata <- read.csv('Sensitivity_analysis_british_population_PRS_FH_foredata_lung_female.csv',header=T)
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
topptx(forestplot,"Sensitivity_analysis_british_population_forestplot_lung_female_FH_PRS.pptx",width = 5,height = 4)



########前列腺癌 PRC PRO_PRS########

#male#
control <- male1[which(male1$FH_Prostate_2==0),]
summary(control$PRO_PRS)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#5.459   7.729   8.188   8.204   8.661  11.561

sd(control$PRO_PRS)  ##0.6924959
male1$PRO_PRS_sd<-male1$PRO_PRS/0.6924959

Q20=quantile(male1$PRO_PRS_sd,seq(0.05,1,0.05))[[4]]
Q80=quantile(male1$PRO_PRS_sd,seq(0.05,1,0.05))[[16]]

####对原始数据的PRS三分位，<=Q20，1；>=Q80，3
male1$PRS_8_sd_3_Q20=0
male1[male1$PRO_PRS_sd<=Q20,]$PRS_8_sd_3_Q20=1
male1[male1$PRO_PRS_sd>Q20&male1$PRO_PRS_sd<Q80,]$PRS_8_sd_3_Q20=2
male1[male1$PRO_PRS_sd>=Q80,]$PRS_8_sd_3_Q20=3

table(male1$PRS_8_sd_3_Q20)
#1      2      3 
#24730 74189 24730   

male1$PRO_PRS_FH=0
male1[male1$PRS_8_sd_3_Q20==1&male1$FH_Prostate_2==0,]$PRO_PRS_FH=1  #对照（1，2，3，4）（1，2，5，6）
male1[male1$PRS_8_sd_3_Q20==1&male1$FH_Prostate_2==1,]$PRO_PRS_FH=2  #家族史效应
male1[male1$PRS_8_sd_3_Q20==2&male1$FH_Prostate_2==0,]$PRO_PRS_FH=3  #PRS中位效应
male1[male1$PRS_8_sd_3_Q20==2&male1$FH_Prostate_2==1,]$PRO_PRS_FH=4  #PRS中位联合家族史效应
male1[male1$PRS_8_sd_3_Q20==3&male1$FH_Prostate_2==0,]$PRO_PRS_FH=5  #PRS高位效应
male1[male1$PRS_8_sd_3_Q20==3&male1$FH_Prostate_2==1,]$PRO_PRS_FH=6  #PRS高位联合家族史效应
table(male1$PRO_PRS_FH)
#1      2      3      4      5      6 
#23265  1465 68403  5786 22223  2507 


#case/人年
fit2 <- pyears(Surv(PRC_difftime_new,PRC_total.y)~1,data=male1[male1$PRO_PRS_FH==1,],scale = 1)
fit2$pyears#
fit2$event#
round(377/249786*100000,2) #150.93

fit2 <- pyears(Surv(PRC_difftime_new,PRC_total.y)~1,data=male1[male1$PRO_PRS_FH==2,],scale = 1)
fit2$pyears#
fit2$event#  
round(34/15830*100000,2) #214.78

fit2 <- pyears(Surv(PRC_difftime_new,PRC_total.y)~1,data=male1[male1$PRO_PRS_FH==3,],scale = 1)
fit2$pyears#
fit2$event#   
round(2792/728222*100000,2) #383.40

fit2 <- pyears(Surv(PRC_difftime_new,PRC_total.y)~1,data=male1[male1$PRO_PRS_FH==4,],scale = 1)
fit2$pyears#
fit2$event#  
round(393/60759*100000,2) #646.82

fit2 <- pyears(Surv(PRC_difftime_new,PRC_total.y)~1,data=male1[male1$PRO_PRS_FH==5,],scale = 1)
fit2$pyears#
fit2$event#   
round(1977/231641*100000,2) #853.48

fit2 <- pyears(Surv(PRC_difftime_new,PRC_total.y)~1,data=male1[male1$PRO_PRS_FH==6,],scale = 1)
fit2$pyears#
fit2$event#  
round(359/25440*100000,2) #1411.16


table(male1$PRC_total.y)
#0      1 
#117717   5932


model1 <- coxph(Surv(PRC_difftime_new,PRC_total.y)~as.factor(PRO_PRS_FH)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male1)
res <- summary(model1)
print(res$coefficients,digits = 2)
print(res$conf.int,digits = 2)
#                                             coef exp(coef) se(coef)       z Pr(>|z|)
#as.factor(PRO_PRS_FH)2               3.1e-01      1.36   0.1791  1.7076  8.8e-02
#as.factor(PRO_PRS_FH)3               9.4e-01      2.56   0.0549 17.1386  7.6e-66
#as.factor(PRO_PRS_FH)4               1.4e+00      4.24   0.0721 20.0466  2.2e-89
#as.factor(PRO_PRS_FH)5               1.8e+00      5.94   0.0562 31.7018 1.5e-220
#as.factor(PRO_PRS_FH)6               2.3e+00      9.86   0.0738 31.0088 4.1e-211
#                                      exp(coef) exp(-coef) lower .95 upper .95
#as.factor(PRO_PRS_FH)2                   1.36       0.74      0.96      1.93
#as.factor(PRO_PRS_FH)3                   2.56       0.39      2.30      2.85
#as.factor(PRO_PRS_FH)4                   4.24       0.24      3.68      4.89
#as.factor(PRO_PRS_FH)5                   5.94       0.17      5.32      6.64
#as.factor(PRO_PRS_FH)6                   9.86       0.10      8.53     11.39

##  相加交互作用检验
#coef的2，3，4，5，6分别对应分组时的2-6组，但是却是生成系数的第1-5行

RERI1=epi.interaction(model=model1,coef=c(1,2,3),param = c("dummy"),conf.level=0.95) 
RERI1  
#$reri
#est     lower    upper
#1 1.325267 0.6971083 1.953426

#$apab
#est     lower     upper
#1 0.3122582 0.1808938 0.4436227

#$s
#est    lower    upper
#1 1.690649 1.271334 2.248264

#$multiplicative
#est   lower    upper
#1 4.244138 3.68477 4.888422



RERI2=epi.interaction(model=model1,coef=c(1,4,5),param = c("dummy"),conf.level=0.95)
RERI2
#$reri
#est    lower   upper
#1 3.554617 2.377954 4.73128

#$apab
#est    lower     upper
#1 0.3606743 0.275021 0.4463276

#$s
#est    lower    upper
#1 1.670574 1.431553 1.949503

#$multiplicative
#est    lower    upper
#1 9.855475 8.528458 11.38897



###森林图
setwd("Sfigure 9 & 15/中间文件/PRS")

   HR <- c(1.0,1.36,NA,2.56,4.24,NA,5.94,9.86)
lower <- c(1.0,0.96,NA,2.30,3.68,NA,5.32,8.53)
upper <- c(1.0,1.93,NA,2.85,4.89,NA,6.64,11.39)

foredata <- data.frame(HR,lower,upper)
foredata
class(foredata)
write.csv(foredata,file = "Sensitivity_analysis_british_population_PRS_FH_foredata_prostate_male.csv",row.names = F)

pdata <- read.csv('Sensitivity_analysis_british_population_PRS_FH_foredata_prostate_male.csv',header=T)
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
  xticks = c(0.5,1.0,3.0,5.0,7.0,9.0,11.0,13.0),
)
forestplot
topptx(forestplot,"Sensitivity_analysis_british_population_forestplot_prostate_male_FH_PRS.pptx",width = 5,height = 4)



####breast-female  BC_PRS####
control <- female1[which(female1$FH_Breast_2==0),]
summary(control$BC_PRS)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-0.74588 -0.19206 -0.10766 -0.10835 -0.02424  0.46687 

sd(control$BC_PRS)  ##0.126784
female1$BC_PRS_sd<-female1$BC_PRS/0.126784

Q20=quantile(female1$BC_PRS_sd,seq(0.05,1,0.05))[[4]]
Q80=quantile(female1$BC_PRS_sd,seq(0.05,1,0.05))[[16]]
####对原始数据的PRS三分位，<=Q20，1；>=Q80，3

female1$PRS_28_sd_3_Q20=0
female1[female1$BC_PRS_sd<=Q20,]$PRS_28_sd_3_Q20=1
female1[female1$BC_PRS_sd>Q20&female1$BC_PRS_sd<Q80,]$PRS_28_sd_3_Q20=2
female1[female1$BC_PRS_sd>=Q80,]$PRS_28_sd_3_Q20=3

table(female1$PRS_28_sd_3_Q20)
#1      2      3 
#28227 84679 28227 

female1$BC_PRS_FH=0
female1[female1$PRS_28_sd_3_Q20==1&female1$FH_Breast_2==0,]$BC_PRS_FH=1  #对照（1，2，3，4）（1，2，5，6）
female1[female1$PRS_28_sd_3_Q20==1&female1$FH_Breast_2==1,]$BC_PRS_FH=2  #家族史效应
female1[female1$PRS_28_sd_3_Q20==2&female1$FH_Breast_2==0,]$BC_PRS_FH=3  #PRS中位效应
female1[female1$PRS_28_sd_3_Q20==2&female1$FH_Breast_2==1,]$BC_PRS_FH=4  #PRS中位联合家族史效应
female1[female1$PRS_28_sd_3_Q20==3&female1$FH_Breast_2==0,]$BC_PRS_FH=5  #PRS高位效应
female1[female1$PRS_28_sd_3_Q20==3&female1$FH_Breast_2==1,]$BC_PRS_FH=6  #PRS高位联合家族史效应
table(female1$BC_PRS_FH)
#1      2      3      4      5      6 
#25701  2526 75430  9249 24429  3798


#case/人年
fit2 <- pyears(Surv(BC_difftime_new,BC_total)~1,data=female1[female1$BC_PRS_FH==1,],scale = 1)
fit2$pyears#
fit2$event#
round(503/280183*100000,2) #179.53

fit2 <- pyears(Surv(BC_difftime_new,BC_total)~1,data=female1[female1$BC_PRS_FH==2,],scale = 1)
fit2$pyears#
fit2$event#  
round(85/27278*100000,2) #311.61

fit2 <- pyears(Surv(BC_difftime_new,BC_total)~1,data=female1[female1$BC_PRS_FH==3,],scale = 1)
fit2$pyears#
fit2$event#   
round(2368/817035*100000,2) #289.83

fit2 <- pyears(Surv(BC_difftime_new,BC_total)~1,data=female1[female1$BC_PRS_FH==4,],scale = 1)
fit2$pyears#
fit2$event#  
round(441/99323*100000,2) #444.01

fit2 <- pyears(Surv(BC_difftime_new,BC_total)~1,data=female1[female1$BC_PRS_FH==5,],scale = 1)
fit2$pyears#
fit2$event#   
round(1251/261438*100000,2) #478.51

fit2 <- pyears(Surv(BC_difftime_new,BC_total)~1,data=female1[female1$BC_PRS_FH==6,],scale = 1)
fit2$pyears#
fit2$event#  
round(265/40296*100000,2) #657.63


table(female1$BC_total)
#0      1 
#136220   4913

model1 <- coxph(Surv(BC_difftime_new,BC_total)~as.factor(BC_PRS_FH)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female1)
res <- summary(model1)
print(res$coefficients,digits = 2)
print(res$conf.int,digits = 2)
#                                            coef exp(coef) se(coef)     z Pr(>|z|)
#as.factor(BC_PRS_FH)2                0.53817      1.71   0.1173  4.59  4.5e-06
#as.factor(BC_PRS_FH)3                0.48126      1.62   0.0491  9.80  1.1e-22
#as.factor(BC_PRS_FH)4                0.89282      2.44   0.0653 13.68  1.3e-42
#as.factor(BC_PRS_FH)5                0.98289      2.67   0.0528 18.61  2.6e-77
#as.factor(BC_PRS_FH)6                1.28973      3.63   0.0759 16.99  1.0e-64
#                                      exp(coef) exp(-coef) lower .95 upper .95
#as.factor(BC_PRS_FH)2                    1.71       0.58      1.36      2.16
#as.factor(BC_PRS_FH)3                    1.62       0.62      1.47      1.78
#as.factor(BC_PRS_FH)4                    2.44       0.41      2.15      2.78
#as.factor(BC_PRS_FH)5                    2.67       0.37      2.41      2.96
#as.factor(BC_PRS_FH)6                    3.63       0.28      3.13      4.21


##  相加交互作用检验
#coef的2，3，4，5，6分别对应分组时的2-6组，但是却是生成系数的第1-5行

RERI1=epi.interaction(model=model1,coef=c(1,2,3),param = c("dummy"),conf.level=0.95) 
RERI1  

#$reri
#est      lower    upper
#1 0.1110357 -0.3303857 0.552457

#$apab
#est      lower     upper
#1 0.04546904 -0.1338884 0.2248265

#$s
#est     lower   upper
#1 1.083425 0.7808629 1.50322

#$multiplicative
#est    lower   upper
#1 2.442006 2.148838 2.77517




RERI2=epi.interaction(model=model1,coef=c(1,4,5),param = c("dummy"),conf.level=0.95)
RERI2

#$reri
#est      lower     upper
#1 0.2467997 -0.3449448 0.8385441

#$apab
#est      lower     upper
#1 0.06795475 -0.0898226 0.2257321

#$s
#est     lower    upper
#1 1.103479 0.8702278 1.399249

#$multiplicative
#est    lower    upper
#1 3.631824 3.129647 4.214578



###森林图
setwd("Sfigure 9 & 15/中间文件/PRS")

   HR <- c(1.0,1.71,NA,1.62,2.44,NA,2.67,3.63)
lower <- c(1.0,1.36,NA,1.47,2.15,NA,2.41,3.13)
upper <- c(1.0,2.16,NA,1.78,2.78,NA,2.96,4.21)

foredata <- data.frame(HR,lower,upper)
foredata
class(foredata)
write.csv(foredata,file = "Sensitivity_analysis_british_population_PRS_FH_foredata_breast_female.csv",row.names = F)

pdata <- read.csv('Sensitivity_analysis_british_population_PRS_FH_foredata_breast_female.csv',header=T)
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
topptx(forestplot,"Sensitivity_analysis_british_population_forestplot_breast_female_FH_PRS.pptx",width = 5,height = 4)




####Sfigure 15

############多肿瘤家族史与CPRS联合效应############

#male#

control <- male1[which(male1$FH_total_2==0),]
summary(control$CPRS)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#2.964   3.609   3.740   3.743   3.874   4.801

sd(control$CPRS)  ##0.1973115
male1$CPRS_sd<-male1$CPRS/0.1973115


Q20=quantile(male1$CPRS_sd,seq(0.05,1,0.05))[[4]]
Q80=quantile(male1$CPRS_sd,seq(0.05,1,0.05))[[16]]

####对原始数据的PRS三分位，<=Q20，1；>=Q80，3
male1$CPRS_sd_3_Q20=0
male1[male1$CPRS_sd<=Q20,]$CPRS_sd_3_Q20=1
male1[male1$CPRS_sd>Q20&male1$CPRS_sd<Q80,]$CPRS_sd_3_Q20=2
male1[male1$CPRS_sd>=Q80,]$CPRS_sd_3_Q20=3

table(male1$CPRS_sd_3_Q20)
#1      2      3 
#24730 74189 24730

#生成新变量COL_PRS_FH，组合PRS及家族史讲人群分为6层
male1$CPRS_FH=0
male1[male1$CPRS_sd_3_Q20==1&male1$FH_total_2==0,]$CPRS_FH=1  #对照（1，2，3，4）（1，2，5，6）
male1[male1$CPRS_sd_3_Q20==1&male1$FH_total_2==1,]$CPRS_FH=2  #家族史效应
male1[male1$CPRS_sd_3_Q20==2&male1$FH_total_2==0,]$CPRS_FH=3  #PRS中位效应
male1[male1$CPRS_sd_3_Q20==2&male1$FH_total_2==1,]$CPRS_FH=4  #PRS中位联合家族史效应
male1[male1$CPRS_sd_3_Q20==3&male1$FH_total_2==0,]$CPRS_FH=5  #PRS高位效应
male1[male1$CPRS_sd_3_Q20==3&male1$FH_total_2==1,]$CPRS_FH=6  #PRS高位联合家族史效应
table(male1$CPRS_FH)
#1      2      3      4      5      6 
#18369  6361 53395 20794 17116  7614  


#case/人年
fit2 <- pyears(Surv(survival_time,cancer_total_20)~1,data=male1[male1$CPRS_FH==1,],scale = 1)
fit2$pyears#
fit2$event#
round(1171/194982*100000,2) #600.57

fit2 <- pyears(Surv(survival_time,cancer_total_20)~1,data=male1[male1$CPRS_FH==2,],scale = 1)
fit2$pyears#
fit2$event#  
round(497/66773*100000,2) #744.31

fit2 <- pyears(Surv(survival_time,cancer_total_20)~1,data=male1[male1$CPRS_FH==3,],scale = 1)
fit2$pyears#
fit2$event#   
round(4764/560022*100000,2) #850.68

fit2 <- pyears(Surv(survival_time,cancer_total_20)~1,data=male1[male1$CPRS_FH==4,],scale = 1)
fit2$pyears#
fit2$event#  
round(2178/216518*100000,2) #1005.92

fit2 <- pyears(Surv(survival_time,cancer_total_20)~1,data=male1[male1$CPRS_FH==5,],scale = 1)
fit2$pyears#
fit2$event#   
round(2289/176115*100000,2) #1299.72

fit2 <- pyears(Surv(survival_time,cancer_total_20)~1,data=male1[male1$CPRS_FH==6,],scale = 1)
fit2$pyears#
fit2$event#  
round(1263/76940*100000,2) #1641.54


table(male1$cancer_total_20)
#0      1 
#111487  12162 

model1 <- coxph(Surv(survival_time,cancer_total_20)~as.factor(CPRS_FH)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male1)
res <- summary(model1)
print(res$coefficients,digit=2)
print(res$conf.int,digit=2)
#                                         coef exp(coef) se(coef)       z Pr(>|z|)
#as.factor(CPRS_FH)2                  1.3e-01      1.14   0.0535  2.3836  1.7e-02
#as.factor(CPRS_FH)3                  3.7e-01      1.44   0.0326 11.2006  4.0e-29
#as.factor(CPRS_FH)4                  4.5e-01      1.57   0.0363 12.4321  1.8e-35
#as.factor(CPRS_FH)5                  8.2e-01      2.28   0.0360 22.8845 6.6e-116
#as.factor(CPRS_FH)6                  9.7e-01      2.65   0.0406 23.9624 6.9e-127

#                                     exp(coef) exp(-coef) lower .95 upper .95
#as.factor(CPRS_FH)2                      1.14       0.88      1.02      1.26
#as.factor(CPRS_FH)3                      1.44       0.69      1.35      1.54
#as.factor(CPRS_FH)4                      1.57       0.64      1.46      1.69
#as.factor(CPRS_FH)5                      2.28       0.44      2.12      2.44
#as.factor(CPRS_FH)6                      2.65       0.38      2.44      2.87


#计算家族史（无/有）和遗传风险（低/中/高）之间的相加交互作用
#设HR11表示2 个危险因素同时存在的HR值，HR10,HR01分别表示存在1个因素时HR值。相加交互作用的评价指标：RERI=HR1-HR01-HR10+1,AP=RERI/HR11

##  相加交互作用检验
#coef的2，3，4，5，6分别对应分组时的2-6组，但是却是生成系数的第1-5行

RERI1=epi.interaction(model=model1,coef=c(1,2,3),param = c("dummy"),conf.level=0.95) 
RERI1  

#$reri
#est      lower     upper
#1 -0.007742566 -0.1467766 0.1312915

#$apab
#est       lower      upper
#1 -0.004933024 -0.09349751 0.08363146

#$s
#est     lower    upper
#1 0.9865878 0.7755228 1.255096

#$multiplicative
#est    lower    upper
#1 1.569537 1.461866 1.685139


RERI2=epi.interaction(model=model1,coef=c(1,4,5),param = c("dummy"),conf.level=0.95)
RERI2
#$reri
#est      lower     upper
#1 0.2326437 0.02790747 0.4373798

#$apab
#est      lower     upper
#1 0.08792737 0.01263748 0.1632173

#$s
#est    lower    upper
#1 1.16462 1.013351 1.338469

#$multiplicative
#est    lower    upper
#1 2.645862 2.443453 2.865038


###森林图
setwd("Sfigure 9 & 15/中间文件/CPRS")

   HR <- c(1.0,1.14,NA,1.44,1.57,NA,2.28,2.65)
lower <- c(1.0,1.02,NA,1.35,1.46,NA,2.12,2.44)
upper <- c(1.0,1.26,NA,1.54,1.69,NA,2.44,2.87)

foredata <- data.frame(HR,lower,upper)
foredata
class(foredata)
write.csv(foredata,file = "Sensitivity_analysis_british_population_CPRS_FH_foredata_overall_male.csv",row.names = F)

pdata <- read.csv('Sensitivity_analysis_british_population_CPRS_FH_foredata_overall_male.csv',header=T)
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
topptx(forestplot,"Sensitivity_analysis_british_population_forestplot_overall_male_FH_CPRS.pptx",width = 5,height = 4)


###female###
control <- female1[which(female1$FH_total_2==0),]
summary(control$CPRS)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.126   1.394   1.444   1.444   1.494   1.825 

sd(control$CPRS)  ##0.07382597
female1$CPRS_sd<-female1$CPRS/0.07382597

Q20=quantile(female1$CPRS_sd,seq(0.05,1,0.05))[[4]]
Q80=quantile(female1$CPRS_sd,seq(0.05,1,0.05))[[16]]

####对原始数据的PRS三分位，<=Q20，1；>=Q80，3
female1$CPRS_sd_3_Q20=0
female1[female1$CPRS_sd<=Q20,]$CPRS_sd_3_Q20=1
female1[female1$CPRS_sd>Q20&female1$CPRS_sd<Q80,]$CPRS_sd_3_Q20=2
female1[female1$CPRS_sd>=Q80,]$CPRS_sd_3_Q20=3

table(female1$CPRS_sd_3_Q20)
#1      2      3 
#28227 84679 28227 

#生成新变量COL_PRS_FH，组合PRS及家族史讲人群分为6层
female1$CPRS_FH=0
female1[female1$CPRS_sd_3_Q20==1&female1$FH_total_2==0,]$CPRS_FH=1  #对照（1，2，3，4）（1，2，5，6）
female1[female1$CPRS_sd_3_Q20==1&female1$FH_total_2==1,]$CPRS_FH=2  #家族史效应
female1[female1$CPRS_sd_3_Q20==2&female1$FH_total_2==0,]$CPRS_FH=3  #PRS中位效应
female1[female1$CPRS_sd_3_Q20==2&female1$FH_total_2==1,]$CPRS_FH=4  #PRS中位联合家族史效应
female1[female1$CPRS_sd_3_Q20==3&female1$FH_total_2==0,]$CPRS_FH=5  #PRS高位效应
female1[female1$CPRS_sd_3_Q20==3&female1$FH_total_2==1,]$CPRS_FH=6  #PRS高位联合家族史效应
table(female1$CPRS_FH)
#1      2      3      4      5      6 
#20248  7979 59109 25570 19009  9218 


#case/人年
fit2 <- pyears(Surv(survival_time,cancer_total_20)~1,data=female1[female1$CPRS_FH==1,],scale = 1)
fit2$pyears#
fit2$event#
round(1164/217669*100000,2) #534.76

fit2 <- pyears(Surv(survival_time,cancer_total_20)~1,data=female1[female1$CPRS_FH==2,],scale = 1)
fit2$pyears#
fit2$event#  
round(554/85277*100000,2) #649.65

fit2 <- pyears(Surv(survival_time,cancer_total_20)~1,data=female1[female1$CPRS_FH==3,],scale = 1)
fit2$pyears#
fit2$event#   
round(4183/631428*100000,2) #662.47

fit2 <- pyears(Surv(survival_time,cancer_total_20)~1,data=female1[female1$CPRS_FH==4,],scale = 1)
fit2$pyears#
fit2$event#  
round(2198/271530*100000,2) #809.49

fit2 <- pyears(Surv(survival_time,cancer_total_20)~1,data=female1[female1$CPRS_FH==5,],scale = 1)
fit2$pyears#
fit2$event#   
round(1702/200916*100000,2) #847.12

fit2 <- pyears(Surv(survival_time,cancer_total_20)~1,data=female1[female1$CPRS_FH==6,],scale = 1)
fit2$pyears#
fit2$event#  
round(969/96945*100000,2) #999.54


model1 <- coxph(Surv(survival_time,cancer_total_20)~as.factor(CPRS_FH)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female1)
res <- summary(model1)
print(res$coefficients,digit=2)
print(res$conf.int,digit=2)

#                                          coef exp(coef) se(coef)       z Pr(>|z|)
#as.factor(CPRS_FH)2                  0.14102      1.15   0.0516   2.731  6.3e-03
#as.factor(CPRS_FH)3                  0.21845      1.24   0.0331   6.591  4.4e-11
#as.factor(CPRS_FH)4                  0.36641      1.44   0.0363  10.096  5.8e-24
#as.factor(CPRS_FH)5                  0.47215      1.60   0.0381  12.406  2.4e-35
#as.factor(CPRS_FH)6                  0.58967      1.80   0.0435  13.543  8.7e-42

#                                   exp(coef) exp(-coef) lower .95 upper .95
#as.factor(CPRS_FH)2                      1.15       0.87      1.04      1.27
#as.factor(CPRS_FH)3                      1.24       0.80      1.17      1.33
#as.factor(CPRS_FH)4                      1.44       0.69      1.34      1.55
#as.factor(CPRS_FH)5                      1.60       0.62      1.49      1.73
#as.factor(CPRS_FH)6                      1.80       0.55      1.66      1.96


#计算家族史（有/无）和遗传风险（低/中/高）之间的相加交互作用
#设HR11表示2 个危险因素同时存在的HR值，HR10,HR01分别表示存在1个因素时HR值。相加交互作用的评价指标：RERI=HR1-HR01-HR10+1,AP=RERI/HR11

##  相加交互作用检验
#coef的2，3，4，5，6分别对应分组时的2-6组，但是却是生成系数的第1-5行

RERI1=epi.interaction(model=model1,coef=c(1,2,3),param = c("dummy"),conf.level=0.95) 
RERI1  

#$reri
#est       lower     upper
#1 0.04695958 -0.08438924 0.1783084

#$apab
#est       lower     upper
#1 0.03255311 -0.05868299 0.1237892

#$s
#est     lower    upper
#1 1.118707 0.8002649 1.563863

#$multiplicative
#est    lower    upper
#1 1.442553 1.343501 1.548907



RERI2=epi.interaction(model=model1,coef=c(1,4,5),param = c("dummy"),conf.level=0.95)
RERI2
#$reri
#est      lower     upper
#1 0.04851606 -0.1272282 0.2242604

#$apab
#est       lower     upper
#1 0.02690262 -0.06995933 0.1237646

#$s
#est     lower    upper
#1 1.06427 0.8457519 1.339247

#$multiplicative
#est    lower    upper
#1 1.803395 1.655878 1.964054



###森林图
setwd("Sfigure 9 & 15/中间文件/CPRS")


   HR <- c(1.0,1.15,NA,1.24,1.44,NA,1.60,1.80)
lower <- c(1.0,1.04,NA,1.17,1.34,NA,1.49,1.66)
upper <- c(1.0,1.27,NA,1.33,1.55,NA,1.73,1.96)

foredata <- data.frame(HR,lower,upper)
foredata
class(foredata)
write.csv(foredata,file = "Sensitivity_analysis_british_population_CPRS_FH_foredata_overall_female.csv",row.names = F)

pdata <- read.csv('Sensitivity_analysis_british_population_CPRS_FH_foredata_overall_female.csv',header=T)
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
topptx(forestplot,"Sensitivity_analysis_british_population_forestplot_overall_female_FH_CPRS.pptx",width = 5,height = 4)

