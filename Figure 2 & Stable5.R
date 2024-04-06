#Figure 2 & Stable 5

#肿瘤家族史与PRS(CPRS)联合效应及交互作用--更新7种肿瘤PRS（单个肿瘤家族史、多肿瘤家族史）


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
library("ggplot2")



########肠癌 CRC  COL_PRS ########

#male  COL_PRS
control <- male[which(male$FH_Bowel_2==0),]
summary(control$COL_PRS)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.718   2.910   3.168   3.180   3.437   5.035  

sd(control$COL_PRS)  ##0.390638
male$COL_PRS_sd<-male$COL_PRS/0.390638

Q20=quantile(male$COL_PRS_sd,seq(0.05,1,0.05))[[4]]
Q80=quantile(male$COL_PRS_sd,seq(0.05,1,0.05))[[16]]

####对原始数据的PRS三分位，<=Q20，1；>=Q80，3
male$PRS_1_sd_3_Q20=0
male[male$COL_PRS_sd<=Q20,]$PRS_1_sd_3_Q20=1
male[male$COL_PRS_sd>Q20&male$COL_PRS_sd<Q80,]$PRS_1_sd_3_Q20=2
male[male$COL_PRS_sd>=Q80,]$PRS_1_sd_3_Q20=3

table(male$PRS_1_sd_3_Q20)
#1      2      3 
#40561 121679  40561 

#生成新变量COL_PRS_FH，组合PRS及家族史将人群分为6层
male$COL_PRS_FH=0
male[male$PRS_1_sd_3_Q20==1&male$FH_Bowel_2==0,]$COL_PRS_FH=1  #对照（1，2，3，4）（1，2，5，6）
male[male$PRS_1_sd_3_Q20==1&male$FH_Bowel_2==1,]$COL_PRS_FH=2  #家族史效应
male[male$PRS_1_sd_3_Q20==2&male$FH_Bowel_2==0,]$COL_PRS_FH=3  #PRS中位效应
male[male$PRS_1_sd_3_Q20==2&male$FH_Bowel_2==1,]$COL_PRS_FH=4  #PRS中位联合家族史效应
male[male$PRS_1_sd_3_Q20==3&male$FH_Bowel_2==0,]$COL_PRS_FH=5  #PRS高位效应
male[male$PRS_1_sd_3_Q20==3&male$FH_Bowel_2==1,]$COL_PRS_FH=6  #PRS高位联合家族史效应
table(male$COL_PRS_FH)
#1      2      3      4      5      6 
#36754   3807 108173  13506  35121   5440  
36754/202801*100
#18.12%


#case/人年
fit2 <- pyears(Surv(CRC_difftime_new,CRC_total.y)~1,data=male[male$COL_PRS_FH==1,],scale = 1)
fit2$pyears#
fit2$event#
#246/396946.8*100000=61.97

fit2 <- pyears(Surv(CRC_difftime_new,CRC_total.y)~1,data=male[male$COL_PRS_FH==2,],scale = 1)
fit2$pyears#
fit2$event# 
#39/40990.3*100000=95.14

fit2 <- pyears(Surv(CRC_difftime_new,CRC_total.y)~1,data=male[male$COL_PRS_FH==3,],scale = 1)
fit2$pyears#
fit2$event# 
#1369/1166179*100000=117.39

fit2 <- pyears(Surv(CRC_difftime_new,CRC_total.y)~1,data=male[male$COL_PRS_FH==4,],scale = 1)
fit2$pyears#
fit2$event#  
#265/145178.7*100000=182.53

fit2 <- pyears(Surv(CRC_difftime_new,CRC_total.y)~1,data=male[male$COL_PRS_FH==5,],scale = 1)
fit2$pyears#
fit2$event# 
#round(769/376592.3*100000,2)=204.20

fit2 <- pyears(Surv(CRC_difftime_new,CRC_total.y)~1,data=male[male$COL_PRS_FH==6,],scale = 1)
fit2$pyears#
fit2$event#   
#150/58199.61*100000=257.73



model1 <- coxph(Surv(CRC_difftime_new,CRC_total.y)~as.factor(COL_PRS_FH)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)
res=summary(model1,digit=10)
round(res$conf.int,2)
#                                        exp(coef) exp(-coef) lower .95 upper .95
#as.factor(COL_PRS_FH)2                   1.34       0.75      0.96      1.88
#as.factor(COL_PRS_FH)3                   1.86       0.54      1.62      2.13
#as.factor(COL_PRS_FH)4                   2.56       0.39      2.15      3.04
#as.factor(COL_PRS_FH)5                   3.20       0.31      2.77      3.69
#as.factor(COL_PRS_FH)6                   3.69       0.27      3.01      4.52




dat1=male

dat_df <- with(dat1,
               data.frame(COL_PRS_FH = c(1, 2, 3, 4, 5, 6), 
                          Age=rep(mean(Age,na.rm = TRUE),6),
                          Height_imp=rep(mean(Height_imp,na.rm = TRUE),6),
                          Townsend_deprivation_index_imp_cat=c(2,2,2,2,2,2),
                          smoking_health=rep(mean(smoking_health,na.rm = TRUE),6),
                          Alcohol_health=rep(mean(Alcohol_health,na.rm = TRUE),6),
                          BMI_health=rep(mean(BMI_health,na.rm = TRUE),6),
                          physical_health=rep(mean(physical_health,na.rm = TRUE),6),
                          PCA1=rep(mean(PCA1,na.rm = TRUE),6),
                          PCA2=rep(mean(PCA2,na.rm = TRUE),6),
                          PCA3=rep(mean(PCA3,na.rm = TRUE),6),
                          PCA4=rep(mean(PCA4,na.rm = TRUE),6),
                          PCA5=rep(mean(PCA5,na.rm = TRUE),6),
                          PCA6=rep(mean(PCA6,na.rm = TRUE),6),
                          PCA7=rep(mean(PCA7,na.rm = TRUE),6),
                          PCA8=rep(mean(PCA8,na.rm = TRUE),6),
                          PCA9=rep(mean(PCA9,na.rm = TRUE),6),
                          PCA10=rep(mean(PCA10,na.rm = TRUE),6)			  
               )
)



dat_df$COL_PRS_FH=as.factor(dat_df$COL_PRS_FH)
dat_df$Townsend_deprivation_index_imp_cat=as.factor(dat_df$Townsend_deprivation_index_imp_cat)


fit2 <- survfit(model1, newdata = dat_df)#survfit()：使用公式或已构建的Cox模型拟合生存曲线

p2=ggsurvplot(fit2,
              linetype=c(1,1,1,1,1,1),
              xlab=c("Years of Follow-up"),
              xlim=c(0,8),
              ylim=c(0,0.02),
              ylab=c("Standardized Colorectal cancer Event Rate"),
              legend=c(0.30,0.8),
              legend.title="FH combined PRS",
              legend.labs=c("Group1 (Reference)","Group2 (HR):1.34 (0.96-1.88)","Group3 (HR):1.86 (1.62-2.13)",
                            "Group4 (HR):2.56 (2.15-3.04)","Group5 (HR):3.20 (2.77-3.69)","Group6 (HR):3.69 (3.01-4.52)"),			
              conf.int=TRUE,
              conf.int.fill="gray",
              conf.int.style="ribbon",
              conf.int.alpha=0.2,
              risk.table.col="strata",
              palette=c("#00468B99","#42B540E5","#0099B4FF","#FDAF91E5","#925E9FE5","#AD002AFF"),
              ggtheme=theme_bw(),
              fun="event",
              font.xtickslab=c(12,"bold"),
              font.ytickslab=c(12,"bold"),
              font.x=c(14,"bold"),
              font.y=c(14,"bold"),
              font.legend=c(10,"bold"),
              censor=FALSE,
              dat=dat_df)
p2

#7.08×5.84

ggsave("Male_75_years_old_Colorectal_Cancer_risk_predtion_of_PRS_combined_FH.pdf", plot = p2, width = 8, height = 6, dpi = 600)

#计算家族史（无/有）和遗传风险（低/中/高）之间的相加交互作用
#设HR11表示2 个危险因素同时存在的HR值，HR10,HR01分别表示存在1个因素时HR值。相加交互作用的评价指标：RERI=HR1-HR01-HR10+1,AP=RERI/HR11
##  相加交互作用检验
#coef的2，3，4，5，6分别对应分组时的2-6组，但是却是生成系数的第1-5行


RERI1=epi.interaction(model=model1,coef=c(1,2,3),param = c("dummy"),conf.level=0.95) 
RERI1  
#$reri
#est    lower     upper
#1 0.3605913 -0.17582 0.8970026

#$apab
#est       lower     upper
#1 0.1409364 -0.06293004 0.3448029

#$s
#est     lower    upper
#1 1.301008 0.8457625 2.001296

#$multiplicative
#est   lower    upper
#1 2.558538 2.15007 3.044607

RERI2=epi.interaction(model=model1,coef=c(1,4,5),param = c("dummy"),conf.level=0.95)
RERI2
#$reri
#est      lower     upper
#1 0.1534309 -0.6134307 0.9202925

#$apab
#est      lower     upper
#1 0.04156593 -0.1619014 0.2450332

#$s
#est     lower    upper
#1 1.060457 0.7905514 1.422513

#$multiplicative
#est    lower    upper
#1 3.691265 3.011603 4.524314







###colorectal-female###
control <- female[which(female$FH_Bowel_2==0),]
summary(control$COL_PRS)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.595   2.911   3.168   3.181   3.437   5.200 

sd(control$COL_PRS)  ##0.3913549
female$COL_PRS_sd<-female$COL_PRS/0.3913549

Q20=quantile(female$COL_PRS_sd,seq(0.05,1,0.05))[[4]]
Q80=quantile(female$COL_PRS_sd,seq(0.05,1,0.05))[[16]]

####对原始数据的PRS三分位，<=Q20，1；>=Q80，3
female$PRS_1_sd_3_Q20=0
female[female$COL_PRS_sd<=Q20,]$PRS_1_sd_3_Q20=1
female[female$COL_PRS_sd>Q20&female$COL_PRS_sd<Q80,]$PRS_1_sd_3_Q20=2
female[female$COL_PRS_sd>=Q80,]$PRS_1_sd_3_Q20=3

table(female$PRS_1_sd_3_Q20)
#1      2      3 
#47920 143758  47920  

female$COL_PRS_FH=0
female[female$PRS_1_sd_3_Q20==1&female$FH_Bowel_2==0,]$COL_PRS_FH=1  #对照（1，2，3，4）（1，2，5，6）
female[female$PRS_1_sd_3_Q20==1&female$FH_Bowel_2==1,]$COL_PRS_FH=2  #家族史效应
female[female$PRS_1_sd_3_Q20==2&female$FH_Bowel_2==0,]$COL_PRS_FH=3  #PRS中位效应
female[female$PRS_1_sd_3_Q20==2&female$FH_Bowel_2==1,]$COL_PRS_FH=4  #PRS中位联合家族史效应
female[female$PRS_1_sd_3_Q20==3&female$FH_Bowel_2==0,]$COL_PRS_FH=5  #PRS高位效应
female[female$PRS_1_sd_3_Q20==3&female$FH_Bowel_2==1,]$COL_PRS_FH=6  #PRS高位联合家族史效应
table(female$COL_PRS_FH)
#1      2      3      4      5      6 
#43725   4195 128413  15345  41662   6258  


#case/人年
fit2 <- pyears(Surv(CRC_difftime_new,CRC_total.y)~1,data=female[female$COL_PRS_FH==1,],scale = 1)
fit2$pyears#
fit2$event#
#round(194/478324.8*100000,2)=40.56

fit2 <- pyears(Surv(CRC_difftime_new,CRC_total.y)~1,data=female[female$COL_PRS_FH==2,],scale = 1)
fit2$pyears#
fit2$event#  
#round(27/45826.83*100000,2)=58.92

fit2 <- pyears(Surv(CRC_difftime_new,CRC_total.y)~1,data=female[female$COL_PRS_FH==3,],scale = 1)
fit2$pyears#
fit2$event#  
#round(1034/1406357*100000,2)=73.52

fit2 <- pyears(Surv(CRC_difftime_new,CRC_total.y)~1,data=female[female$COL_PRS_FH==4,],scale = 1)
fit2$pyears#
fit2$event#  
#round(185/167639.3*100000,2)=110.36

fit2 <- pyears(Surv(CRC_difftime_new,CRC_total.y)~1,data=female[female$COL_PRS_FH==5,],scale = 1)
fit2$pyears#
fit2$event#
#round(625/455008.6*100000,2)=137.36

fit2 <- pyears(Surv(CRC_difftime_new,CRC_total.y)~1,data=female[female$COL_PRS_FH==6,],scale = 1)
fit2$pyears#
fit2$event#  
#round(102/68533.01*100000,2)=148.83



model1 <- coxph(Surv(CRC_difftime_new,CRC_total.y)~as.factor(COL_PRS_FH)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)
res=summary(model1,digit=10)
round(res$conf.int,2)
#                                       exp(coef) exp(-coef) lower .95 upper .95
#as.factor(COL_PRS_FH)2                   1.30       0.77      0.87      1.94
#as.factor(COL_PRS_FH)3                   1.80       0.55      1.55      2.10
#as.factor(COL_PRS_FH)4                   2.46       0.41      2.01      3.01
#as.factor(COL_PRS_FH)5                   3.38       0.30      2.87      3.97
#as.factor(COL_PRS_FH)6                   3.33       0.30      2.62      4.23



dat1=female

dat_df <- with(dat1,
               data.frame(COL_PRS_FH = c(1, 2, 3, 4, 5, 6), 
                          Age=rep(mean(Age,na.rm = TRUE),6),
                          Height_imp=rep(mean(Height_imp,na.rm = TRUE),6),
                          Townsend_deprivation_index_imp_cat=c(2,2,2,2,2,2),
                          smoking_health=rep(mean(smoking_health,na.rm = TRUE),6),
                          Alcohol_health=rep(mean(Alcohol_health,na.rm = TRUE),6),
                          BMI_health=rep(mean(BMI_health,na.rm = TRUE),6),
                          physical_health=rep(mean(physical_health,na.rm = TRUE),6),
                          PCA1=rep(mean(PCA1,na.rm = TRUE),6),
                          PCA2=rep(mean(PCA2,na.rm = TRUE),6),
                          PCA3=rep(mean(PCA3,na.rm = TRUE),6),
                          PCA4=rep(mean(PCA4,na.rm = TRUE),6),
                          PCA5=rep(mean(PCA5,na.rm = TRUE),6),
                          PCA6=rep(mean(PCA6,na.rm = TRUE),6),
                          PCA7=rep(mean(PCA7,na.rm = TRUE),6),
                          PCA8=rep(mean(PCA8,na.rm = TRUE),6),
                          PCA9=rep(mean(PCA9,na.rm = TRUE),6),
                          PCA10=rep(mean(PCA10,na.rm = TRUE),6)			  
               )
)



dat_df$COL_PRS_FH=as.factor(dat_df$COL_PRS_FH)
dat_df$Townsend_deprivation_index_imp_cat=as.factor(dat_df$Townsend_deprivation_index_imp_cat)


fit2 <- survfit(model1, newdata = dat_df)#survfit()：使用公式或已构建的Cox模型拟合生存曲线

p2=ggsurvplot(fit2,
              linetype=c(1,1,1,1,1,1),
              xlab=c("Years of Follow-up"),
              xlim=c(0,8),
              ylim=c(0,0.01),
              ylab=c("Standardized Colorectal cancer Event Rate"),
              legend=c(0.30,0.8),
              legend.title="FH combined PRS",
              legend.labs=c("Group1 (Reference)","Group2 (HR):1.30 (0.87-1.94)","Group3 (HR):1.80 (1.55-2.10)",
                            "Group4 (HR):2.46 (2.01-3.01)","Group5 (HR):3.38 (2.87-3.97)","Group6 (HR):3.33 (2.62-4.23)"),			
              conf.int=TRUE,
              conf.int.fill="gray",
              conf.int.style="ribbon",
              conf.int.alpha=0.2,
              risk.table.col="strata",
              palette=c("#00468B99","#42B540E5","#0099B4FF","#FDAF91E5","#925E9FE5","#AD002AFF"),
              ggtheme=theme_bw(),
              fun="event",
              font.xtickslab=c(12,"bold"),
              font.ytickslab=c(12,"bold"),
              font.x=c(14,"bold"),
              font.y=c(14,"bold"),
              font.legend=c(10,"bold"),
              censor=FALSE,
              dat=dat_df)
p2

#7.08×5.84

ggsave("Female_75_years_old_Colorectal_Cancer_risk_predtion_of_PRS_combined_FH.pdf", plot = p2, width = 8, height = 6, dpi = 600)


##  相加交互作用检验
#coef的2，3，4，5，6分别对应分组时的2-6组，但是却是生成系数的第1-5行

RERI1=epi.interaction(model=model1,coef=c(1,2,3),param = c("dummy"),conf.level=0.95) 
RERI1  
#$reri
#est      lower     upper
#1 0.3550867 -0.2664712 0.9766446

#$apab
#est      lower     upper
#1 0.1444498 -0.1010502 0.3899497

#$s
#est     lower  upper
#1 1.321894 0.7702571 2.2686

#$multiplicative
#est    lower    upper
#1 2.458202 2.008626 3.008404

RERI2=epi.interaction(model=model1,coef=c(1,4,5),param = c("dummy"),conf.level=0.95)
RERI2
#$reri
#est     lower     upper
#1 -0.3473025 -1.221454 0.5268492

#$apab
#est      lower     upper
#1 -0.104323 -0.3795808 0.1709347

#$s
#est     lower    upper
#1 0.8702357 0.6133717 1.234667

#$multiplicative
#est    lower    upper
#1 3.329107 2.617586 4.234038





########肺癌 LC LC_PRS########

#male
control <- male[which(male$FH_Lung_2==0),]
summary(control$LC_PRS)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-0.378934  0.009072  0.065275  0.083848  0.139785  1.019890 

sd(control$LC_PRS)  ##0.1093194
male$LC_PRS_sd<-male$LC_PRS/0.1093194

Q20=quantile(male$LC_PRS_sd,seq(0.05,1,0.05))[[4]]
Q80=quantile(male$LC_PRS_sd,seq(0.05,1,0.05))[[16]]

####对原始数据的PRS三分位，<=Q20，1；>=Q80，3
male$PRS_7_sd_3_Q20=0
male[male$LC_PRS_sd<=Q20,]$PRS_7_sd_3_Q20=1
male[male$LC_PRS_sd>Q20&male$LC_PRS_sd<Q80,]$PRS_7_sd_3_Q20=2
male[male$LC_PRS_sd>=Q80,]$PRS_7_sd_3_Q20=3
table(male$PRS_7_sd_3_Q20)
#1      2      3 
#40561 121679  40561  

male$LC_PRS_FH=0
male[male$PRS_7_sd_3_Q20==1&male$FH_Lung_2==0,]$LC_PRS_FH=1  #对照（1，2，3，4）（1，2，5，6）
male[male$PRS_7_sd_3_Q20==1&male$FH_Lung_2==1,]$LC_PRS_FH=2  #家族史效应
male[male$PRS_7_sd_3_Q20==2&male$FH_Lung_2==0,]$LC_PRS_FH=3  #PRS中位效应
male[male$PRS_7_sd_3_Q20==2&male$FH_Lung_2==1,]$LC_PRS_FH=4  #PRS中位联合家族史效应
male[male$PRS_7_sd_3_Q20==3&male$FH_Lung_2==0,]$LC_PRS_FH=5  #PRS高位效应
male[male$PRS_7_sd_3_Q20==3&male$FH_Lung_2==1,]$LC_PRS_FH=6  #PRS高位联合家族史效应
table(male$LC_PRS_FH)
#1      2      3      4      5      6 
#36083   4478 107075  14604  34988   5573 


#case/人年
fit2 <- pyears(Surv(LC_difftime_new,LC_total.y)~1,data=male[male$LC_PRS_FH==1,],scale = 1)
fit2$pyears#
fit2$event#
#round(180/390603.8*100000,2)=46.08

fit2 <- pyears(Surv(LC_difftime_new,LC_total.y)~1,data=male[male$LC_PRS_FH==2,],scale = 1)
fit2$pyears#
fit2$event#  
#round(46/48311.66*100000,2)=95.22

fit2 <- pyears(Surv(LC_difftime_new,LC_total.y)~1,data=male[male$LC_PRS_FH==3,],scale = 1)
fit2$pyears#
fit2$event#   
#round(863/1158021*100000,2)=74.52

fit2 <- pyears(Surv(LC_difftime_new,LC_total.y)~1,data=male[male$LC_PRS_FH==4,],scale = 1)
fit2$pyears#
fit2$event#  
#round(234/157498.6*100000,2)=148.57

fit2 <- pyears(Surv(LC_difftime_new,LC_total.y)~1,data=male[male$LC_PRS_FH==5,],scale = 1)
fit2$pyears#
fit2$event#   
#round(383/378354.8*100000,2)=101.23

fit2 <- pyears(Surv(LC_difftime_new,LC_total.y)~1,data=male[male$LC_PRS_FH==6,],scale = 1)
fit2$pyears#
fit2$event#  
#round(101/59971.26*100000,2)=168.41



model1 <- coxph(Surv(LC_difftime_new,LC_total.y)~as.factor(LC_PRS_FH)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)
res=summary(model1,digit=10)
round(res$conf.int,2)
#                                       exp(coef) exp(-coef) lower .95 upper .95
#as.factor(LC_PRS_FH)2                    1.69       0.59      1.22      2.34
#as.factor(LC_PRS_FH)3                    1.64       0.61      1.40      1.93
#as.factor(LC_PRS_FH)4                    2.73       0.37      2.24      3.31
#as.factor(LC_PRS_FH)5                    2.26       0.44      1.89      2.69
#as.factor(LC_PRS_FH)6                    3.11       0.32      2.44      3.98




dat1=male

dat_df <- with(dat1,
               data.frame(LC_PRS_FH = c(1, 2, 3, 4, 5, 6), 
                          Age=rep(mean(Age,na.rm = TRUE),6),
                          Height_imp=rep(mean(Height_imp,na.rm = TRUE),6),
                          Townsend_deprivation_index_imp_cat=c(2,2,2,2,2,2),
                          smoking_health=rep(mean(smoking_health,na.rm = TRUE),6),
                          Alcohol_health=rep(mean(Alcohol_health,na.rm = TRUE),6),
                          BMI_health=rep(mean(BMI_health,na.rm = TRUE),6),
                          physical_health=rep(mean(physical_health,na.rm = TRUE),6),
                          PCA1=rep(mean(PCA1,na.rm = TRUE),6),
                          PCA2=rep(mean(PCA2,na.rm = TRUE),6),
                          PCA3=rep(mean(PCA3,na.rm = TRUE),6),
                          PCA4=rep(mean(PCA4,na.rm = TRUE),6),
                          PCA5=rep(mean(PCA5,na.rm = TRUE),6),
                          PCA6=rep(mean(PCA6,na.rm = TRUE),6),
                          PCA7=rep(mean(PCA7,na.rm = TRUE),6),
                          PCA8=rep(mean(PCA8,na.rm = TRUE),6),
                          PCA9=rep(mean(PCA9,na.rm = TRUE),6),
                          PCA10=rep(mean(PCA10,na.rm = TRUE),6)			  
               )
)



dat_df$LC_PRS_FH=as.factor(dat_df$LC_PRS_FH)
dat_df$Townsend_deprivation_index_imp_cat=as.factor(dat_df$Townsend_deprivation_index_imp_cat)


fit2 <- survfit(model1, newdata = dat_df)#survfit()：使用公式或已构建的Cox模型拟合生存曲线

p2=ggsurvplot(fit2,
              linetype=c(1,1,1,1,1,1),
              xlab=c("Years of Follow-up"),
              xlim=c(0,7),
              ylim=c(0,0.005),
              ylab=c("Standardized Lung cancer Event Rate"),
              legend=c(0.30,0.8),
              legend.title="FH combined PRS",
              legend.labs=c("Group1 (Reference)","Group2 (HR):1.69 (1.22-2.34)","Group3 (HR):1.64 (1.40-1.93)",
                            "Group4 (HR):2.73 (2.24-3.31)","Group5 (HR):2.26 (1.89-2.69)","Group6 (HR):3.11 (2.44-3.98)"),			
              conf.int=TRUE,
              conf.int.fill="gray",
              conf.int.style="ribbon",
              conf.int.alpha=0.2,
              risk.table.col="strata",
              palette=c("#00468B99","#42B540E5","#0099B4FF","#FDAF91E5","#925E9FE5","#AD002AFF"),
              ggtheme=theme_bw(),
              fun="event",
              font.xtickslab=c(12,"bold"),
              font.ytickslab=c(12,"bold"),
              font.x=c(14,"bold"),
              font.y=c(14,"bold"),
              font.legend=c(10,"bold"),
              censor=FALSE,
              dat=dat_df)
p2

#7.08×5.84

ggsave("Male_75_years_old_Lung_Cancer_risk_predtion_of_PRS_combined_FH.pdf", plot = p2, width = 8, height = 6, dpi = 600)


##  相加交互作用检验
#coef的2，3，4，5，6分别对应分组时的2-6组，但是却是生成系数的第1-5行


RERI1=epi.interaction(model=model1,coef=c(1,2,3),param = c("dummy"),conf.level=0.95) 
RERI1  
#reri
#est      lower    upper
#1 0.3912444 -0.2259501 1.008439

#$apab
#est       lower     upper
#1 0.1435207 -0.07712207 0.3641634

#$s
#est     lower    upper
#1 1.29311 0.8270027 2.021919

#$multiplicative
#est    lower    upper
#1 2.726049 2.243839 3.311889

RERI2=epi.interaction(model=model1,coef=c(1,4,5),param = c("dummy"),conf.level=0.95)
RERI2
#$reri
#est      lower     upper
#1 0.168027 -0.6529164 0.9889704

#$apab
#est      lower     upper
#1 0.05395067 -0.2032625 0.3111639

#$s
#est     lower    upper
#1 1.086326 0.7230864 1.632037

#$multiplicative
#est    lower    upper
#1 3.114457 2.440119 3.975149





###lung-female###
control <- female[which(female$FH_Lung_2==0),]
summary(control$LC_PRS)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-0.401799  0.009546  0.065659  0.084044  0.139884  1.199014 

sd(control$LC_PRS)  ##0.109029
female$LC_PRS_sd<-female$LC_PRS/0.109029

Q20=quantile(female$LC_PRS_sd,seq(0.05,1,0.05))[[4]]
Q80=quantile(female$LC_PRS_sd,seq(0.05,1,0.05))[[16]]

####对原始数据的PRS三分位，<=Q20，1；>=Q80，3
female$PRS_7_sd_3_Q20=0
female[female$LC_PRS_sd<=Q20,]$PRS_7_sd_3_Q20=1
female[female$LC_PRS_sd>Q20&female$LC_PRS_sd<Q80,]$PRS_7_sd_3_Q20=2
female[female$LC_PRS_sd>=Q80,]$PRS_7_sd_3_Q20=3
table(female$PRS_7_sd_3_Q20)
#1      2      3 
#47920 143758  47920  

female$LC_PRS_FH=0
female[female$PRS_7_sd_3_Q20==1&female$FH_Lung_2==0,]$LC_PRS_FH=1  #对照（1，2，3，4）（1，2，5，6）
female[female$PRS_7_sd_3_Q20==1&female$FH_Lung_2==1,]$LC_PRS_FH=2  #家族史效应
female[female$PRS_7_sd_3_Q20==2&female$FH_Lung_2==0,]$LC_PRS_FH=3  #PRS中位效应
female[female$PRS_7_sd_3_Q20==2&female$FH_Lung_2==1,]$LC_PRS_FH=4  #PRS中位联合家族史效应
female[female$PRS_7_sd_3_Q20==3&female$FH_Lung_2==0,]$LC_PRS_FH=5  #PRS高位效应
female[female$PRS_7_sd_3_Q20==3&female$FH_Lung_2==1,]$LC_PRS_FH=6  #PRS高位联合家族史效应
table(female$LC_PRS_FH)
#1      2      3      4      5      6 
#42427   5493 125539  18219  41283   6637 


#case/人年
fit2 <- pyears(Surv(LC_difftime_new,LC_total.y)~1,data=female[female$LC_PRS_FH==1,],scale = 1)
fit2$pyears#
fit2$event#
#round(188/465458.5*100000,2)=40.39

fit2 <- pyears(Surv(LC_difftime_new,LC_total.y)~1,data=female[female$LC_PRS_FH==2,],scale = 1)
fit2$pyears#
fit2$event#  
#round(43/60501.15*100000,2)=71.07

fit2 <- pyears(Surv(LC_difftime_new,LC_total.y)~1,data=female[female$LC_PRS_FH==3,],scale = 1)
fit2$pyears#
fit2$event#   
#round(785/1376645*100000,2)=57.02

fit2 <- pyears(Surv(LC_difftime_new,LC_total.y)~1,data=female[female$LC_PRS_FH==4,],scale = 1)
fit2$pyears#
fit2$event#  
#round(235/199539.1*100000,2)=117.77

fit2 <- pyears(Surv(LC_difftime_new,LC_total.y)~1,data=female[female$LC_PRS_FH==5,],scale = 1)
fit2$pyears#
fit2$event#   
#round(327/452367.6*100000,2)=72.29

fit2 <- pyears(Surv(LC_difftime_new,LC_total.y)~1,data=female[female$LC_PRS_FH==6,],scale = 1)
fit2$pyears#
fit2$event#  
#round(104/72867.82*100000,2)=142.72



model1 <- coxph(Surv(LC_difftime_new,LC_total.y)~as.factor(LC_PRS_FH)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)
res=summary(model1,digit=10)
round(res$conf.int,2)
#                                        exp(coef) exp(-coef) lower .95 upper .95
#as.factor(LC_PRS_FH)2                    1.40       0.72      1.00      1.94
#as.factor(LC_PRS_FH)3                    1.43       0.70      1.22      1.67
#as.factor(LC_PRS_FH)4                    2.35       0.43      1.94      2.85
#as.factor(LC_PRS_FH)5                    1.83       0.55      1.53      2.20
#as.factor(LC_PRS_FH)6                    2.98       0.34      2.34      3.79



dat1=female

dat_df <- with(dat1,
               data.frame(LC_PRS_FH = c(1, 2, 3, 4, 5, 6), 
                          Age=rep(mean(Age,na.rm = TRUE),6),
                          Height_imp=rep(mean(Height_imp,na.rm = TRUE),6),
                          Townsend_deprivation_index_imp_cat=c(2,2,2,2,2,2),
                          smoking_health=rep(mean(smoking_health,na.rm = TRUE),6),
                          Alcohol_health=rep(mean(Alcohol_health,na.rm = TRUE),6),
                          BMI_health=rep(mean(BMI_health,na.rm = TRUE),6),
                          physical_health=rep(mean(physical_health,na.rm = TRUE),6),
                          PCA1=rep(mean(PCA1,na.rm = TRUE),6),
                          PCA2=rep(mean(PCA2,na.rm = TRUE),6),
                          PCA3=rep(mean(PCA3,na.rm = TRUE),6),
                          PCA4=rep(mean(PCA4,na.rm = TRUE),6),
                          PCA5=rep(mean(PCA5,na.rm = TRUE),6),
                          PCA6=rep(mean(PCA6,na.rm = TRUE),6),
                          PCA7=rep(mean(PCA7,na.rm = TRUE),6),
                          PCA8=rep(mean(PCA8,na.rm = TRUE),6),
                          PCA9=rep(mean(PCA9,na.rm = TRUE),6),
                          PCA10=rep(mean(PCA10,na.rm = TRUE),6)			  
               )
)



dat_df$LC_PRS_FH=as.factor(dat_df$LC_PRS_FH)
dat_df$Townsend_deprivation_index_imp_cat=as.factor(dat_df$Townsend_deprivation_index_imp_cat)


fit2 <- survfit(model1, newdata = dat_df)#survfit()：使用公式或已构建的Cox模型拟合生存曲线

p2=ggsurvplot(fit2,
              linetype=c(1,1,1,1,1,1),
              xlab=c("Years of Follow-up"),
              xlim=c(0,7),
              ylim=c(0,0.005),
              ylab=c("Standardized Lung cancer Event Rate"),
              legend=c(0.30,0.8),
              legend.title="FH combined PRS",
              legend.labs=c("Group1 (Reference)","Group2 (HR):1.40 (1.00-1.94)","Group3 (HR):1.43 (1.22-1.67)",
                            "Group4 (HR):2.35 (1.94-2.85)","Group5 (HR):1.83 (1.53-2.20)","Group6 (HR):2.98 (2.34-3.79)"),			
              conf.int=TRUE,
              conf.int.fill="gray",
              conf.int.style="ribbon",
              conf.int.alpha=0.2,
              risk.table.col="strata",
              palette=c("#00468B99","#42B540E5","#0099B4FF","#FDAF91E5","#925E9FE5","#AD002AFF"),
              ggtheme=theme_bw(),
              fun="event",
              font.xtickslab=c(12,"bold"),
              font.ytickslab=c(12,"bold"),
              font.x=c(14,"bold"),
              font.y=c(14,"bold"),
              font.legend=c(10,"bold"),
              censor=FALSE,
              dat=dat_df)
p2

#7.08×5.84

ggsave("Female_75_years_old_Lung_Cancer_risk_predtion_of_PRS_combined_FH.pdf", plot = p2, width = 8, height = 6, dpi = 600)

##  相加交互作用检验
#coef的2，3，4，5，6分别对应分组时的2-6组，但是却是生成系数的第1-5行

RERI1=epi.interaction(model=model1,coef=c(1,2,3),param = c("dummy"),conf.level=0.95) 
RERI1  

#$reri
#est         lower    upper
#1 0.5276855 -0.0005259251 1.055897

#$apab
#est      lower     upper
#1 0.2244439 0.00807556 0.4408122

#$s
#est     lower   upper
#1 1.640866 0.8864619 3.03729

#$multiplicative
#est    lower    upper
#1 2.35108 1.939994 2.849276


RERI2=epi.interaction(model=model1,coef=c(1,4,5),param = c("dummy"),conf.level=0.95)
RERI2

#$reri
#est      lower    upper
#1 0.7500371 0.01315726 1.486917

#$apab
#est      lower     upper
#1 0.2517154 0.03482424 0.4686065

#$s
#est     lower    upper
#1 1.609952 0.9758817 2.656003

#$multiplicative
#est    lower    upper
#1 2.979703 2.344394 3.787175






########前列腺癌 PRC PRO_PRS########

#male#
control <- male[which(male$FH_Prostate_2==0),]
summary(control$PRO_PRS)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#5.459   7.736   8.199   8.214   8.675  11.766 

sd(control$PRO_PRS)  ##0.6953242
male$PRO_PRS_sd<-male$PRO_PRS/0.6953242

Q20=quantile(male$PRO_PRS_sd,seq(0.05,1,0.05))[[4]]
Q80=quantile(male$PRO_PRS_sd,seq(0.05,1,0.05))[[16]]

####对原始数据的PRS三分位，<=Q20，1；>=Q80，3
male$PRS_8_sd_3_Q20=0
male[male$PRO_PRS_sd<=Q20,]$PRS_8_sd_3_Q20=1
male[male$PRO_PRS_sd>Q20&male$PRO_PRS_sd<Q80,]$PRS_8_sd_3_Q20=2
male[male$PRO_PRS_sd>=Q80,]$PRS_8_sd_3_Q20=3
table(male$PRS_8_sd_3_Q20)
#1      2      3 
#40561 121679  40561   

male$PRO_PRS_FH=0
male[male$PRS_8_sd_3_Q20==1&male$FH_Prostate_2==0,]$PRO_PRS_FH=1  #对照（1，2，3，4）（1，2，5，6）
male[male$PRS_8_sd_3_Q20==1&male$FH_Prostate_2==1,]$PRO_PRS_FH=2  #家族史效应
male[male$PRS_8_sd_3_Q20==2&male$FH_Prostate_2==0,]$PRO_PRS_FH=3  #PRS中位效应
male[male$PRS_8_sd_3_Q20==2&male$FH_Prostate_2==1,]$PRO_PRS_FH=4  #PRS中位联合家族史效应
male[male$PRS_8_sd_3_Q20==3&male$FH_Prostate_2==0,]$PRO_PRS_FH=5  #PRS高位效应
male[male$PRS_8_sd_3_Q20==3&male$FH_Prostate_2==1,]$PRO_PRS_FH=6  #PRS高位联合家族史效应
table(male$PRO_PRS_FH)
#1      2      3      4      5      6 
#38147   2414 112266   9413  36615   3946  


#case/人年
fit2 <- pyears(Surv(PRC_difftime_new,PRC_total.y)~1,data=male[male$PRO_PRS_FH==1,],scale = 1)
fit2$pyears#
fit2$event#
#round(650/409560*100000,2)=158.71

fit2 <- pyears(Surv(PRC_difftime_new,PRC_total.y)~1,data=male[male$PRO_PRS_FH==2,],scale = 1)
fit2$pyears#
fit2$event#  
#round(56/25995.15*100000,2)=215.42

fit2 <- pyears(Surv(PRC_difftime_new,PRC_total.y)~1,data=male[male$PRO_PRS_FH==3,],scale = 1)
fit2$pyears#
fit2$event#   
#round(4472/1195184*100000,2)=374.17

fit2 <- pyears(Surv(PRC_difftime_new,PRC_total.y)~1,data=male[male$PRO_PRS_FH==4,],scale = 1)
fit2$pyears#
fit2$event#  
#round(636/98747.46*100000,2)=644.07

fit2 <- pyears(Surv(PRC_difftime_new,PRC_total.y)~1,data=male[male$PRO_PRS_FH==5,],scale = 1)
fit2$pyears#
fit2$event#   
#round(3225/381134.6*100000,2)=846.16

fit2 <- pyears(Surv(PRC_difftime_new,PRC_total.y)~1,data=male[male$PRO_PRS_FH==6,],scale = 1)
fit2$pyears#
fit2$event#  
#round(546/40009.64*100000,2)=1364.67





model1 <- coxph(Surv(PRC_difftime_new,PRC_total.y)~as.factor(PRO_PRS_FH)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)
res=summary(model1,digit=10)
round(res$conf.int,2)
#                                        exp(coef) exp(-coef) lower .95 upper .95
#as.factor(PRO_PRS_FH)2                   1.29       0.77      0.98      1.70
#as.factor(PRO_PRS_FH)3                   2.39       0.42      2.20      2.60
#as.factor(PRO_PRS_FH)4                   3.98       0.25      3.56      4.44
#as.factor(PRO_PRS_FH)5                   5.60       0.18      5.14      6.09
#as.factor(PRO_PRS_FH)6                   8.87       0.11      7.92      9.94




dat1=male

dat_df <- with(dat1,
               data.frame(PRO_PRS_FH = c(1, 2, 3, 4, 5, 6), 
                          Age=rep(mean(Age,na.rm = TRUE),6),
                          Height_imp=rep(mean(Height_imp,na.rm = TRUE),6),
                          Townsend_deprivation_index_imp_cat=c(2,2,2,2,2,2),
                          smoking_health=rep(mean(smoking_health,na.rm = TRUE),6),
                          Alcohol_health=rep(mean(Alcohol_health,na.rm = TRUE),6),
                          BMI_health=rep(mean(BMI_health,na.rm = TRUE),6),
                          physical_health=rep(mean(physical_health,na.rm = TRUE),6),
                          PCA1=rep(mean(PCA1,na.rm = TRUE),6),
                          PCA2=rep(mean(PCA2,na.rm = TRUE),6),
                          PCA3=rep(mean(PCA3,na.rm = TRUE),6),
                          PCA4=rep(mean(PCA4,na.rm = TRUE),6),
                          PCA5=rep(mean(PCA5,na.rm = TRUE),6),
                          PCA6=rep(mean(PCA6,na.rm = TRUE),6),
                          PCA7=rep(mean(PCA7,na.rm = TRUE),6),
                          PCA8=rep(mean(PCA8,na.rm = TRUE),6),
                          PCA9=rep(mean(PCA9,na.rm = TRUE),6),
                          PCA10=rep(mean(PCA10,na.rm = TRUE),6)			  
               )
)



dat_df$PRO_PRS_FH=as.factor(dat_df$PRO_PRS_FH)
dat_df$Townsend_deprivation_index_imp_cat=as.factor(dat_df$Townsend_deprivation_index_imp_cat)


fit2 <- survfit(model1, newdata = dat_df)#survfit()：使用公式或已构建的Cox模型拟合生存曲线

p2=ggsurvplot(fit2,
              linetype=c(1,1,1,1,1,1),
              xlab=c("Years of Follow-up"),
              xlim=c(0,7),
              ylim=c(0,0.08),
              ylab=c("Standardized Prostate cancer Event Rate"),
              legend=c(0.30,0.8),
              legend.title="FH combined PRS",
              legend.labs=c("Group1 (Reference)","Group2 (HR):1.29 (0.98-1.70)","Group3 (HR):2.39 (2.20-2.60)",
                            "Group4 (HR):3.98 (3.56-4.44)","Group5 (HR):5.60 (5.14-6.09)","Group6 (HR):8.87 (7.92-9.94)"),			
              conf.int=TRUE,
              conf.int.fill="gray",
              conf.int.style="ribbon",
              conf.int.alpha=0.2,
              risk.table.col="strata",
              palette=c("#00468B99","#42B540E5","#0099B4FF","#FDAF91E5","#925E9FE5","#AD002AFF"),
              ggtheme=theme_bw(),
              fun="event",
              font.xtickslab=c(12,"bold"),
              font.ytickslab=c(12,"bold"),
              font.x=c(14,"bold"),
              font.y=c(14,"bold"),
              font.legend=c(10,"bold"),
              censor=FALSE,
              dat=dat_df)
p2

#7.08×5.84

ggsave("Male_75_years_old_Prostate_Cancer_risk_predtion_of_PRS_combined_FH.pdf", plot = p2, width = 8, height = 6, dpi = 600)


##  相加交互作用检验
#coef的2，3，4，5，6分别对应分组时的2-6组，但是却是生成系数的第1-5行

RERI1=epi.interaction(model=model1,coef=c(1,2,3),param = c("dummy"),conf.level=0.95) 
RERI1  
#$reri
#est     lower    upper
#1 1.295007 0.8310488 1.758965

#$apab
#est     lower     upper
#1 0.3256993 0.2224539 0.4289447

#$s
#est    lower   upper
#1 1.770345 1.393587 2.24896

#$multiplicative
#est    lower    upper
#1 3.976082 3.564156 4.435616

RERI2=epi.interaction(model=model1,coef=c(1,4,5),param = c("dummy"),conf.level=0.95)
RERI2
#$reri
#est    lower    upper
#1 2.982429 2.128257 3.836601

#$apab
#est     lower     upper
#1 0.3361904 0.2647315 0.4076492

#$s
#est    lower    upper
#1 1.610051 1.420077 1.825439

#$multiplicative
#est    lower    upper
#1 8.871251 7.915796 9.942032





####breast-female  BC_PRS####
control <- female[which(female$FH_Breast_2==0),]
summary(control$BC_PRS)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-0.7459 -0.1892 -0.1044 -0.1053 -0.0206  0.5094 

sd(control$BC_PRS)  ##0.1271848
female$BC_PRS_sd<-female$BC_PRS/0.1271848

Q20=quantile(female$BC_PRS_sd,seq(0.05,1,0.05))[[4]]
Q80=quantile(female$BC_PRS_sd,seq(0.05,1,0.05))[[16]]
####对原始数据的PRS三分位，<=Q20，1；>=Q80，3

female$PRS_28_sd_3_Q20=0
female[female$BC_PRS_sd<=Q20,]$PRS_28_sd_3_Q20=1
female[female$BC_PRS_sd>Q20&female$BC_PRS_sd<Q80,]$PRS_28_sd_3_Q20=2
female[female$BC_PRS_sd>=Q80,]$PRS_28_sd_3_Q20=3
table(female$PRS_28_sd_3_Q20)
#1      2      3 
#47920 143758  47920  

female$BC_PRS_FH=0
female[female$PRS_28_sd_3_Q20==1&female$FH_Breast_2==0,]$BC_PRS_FH=1  #对照（1，2，3，4）（1，2，5，6）
female[female$PRS_28_sd_3_Q20==1&female$FH_Breast_2==1,]$BC_PRS_FH=2  #家族史效应
female[female$PRS_28_sd_3_Q20==2&female$FH_Breast_2==0,]$BC_PRS_FH=3  #PRS中位效应
female[female$PRS_28_sd_3_Q20==2&female$FH_Breast_2==1,]$BC_PRS_FH=4  #PRS中位联合家族史效应
female[female$PRS_28_sd_3_Q20==3&female$FH_Breast_2==0,]$BC_PRS_FH=5  #PRS高位效应
female[female$PRS_28_sd_3_Q20==3&female$FH_Breast_2==1,]$BC_PRS_FH=6  #PRS高位联合家族史效应
table(female$BC_PRS_FH)
#1      2      3      4      5      6 
#43578   4342 128141  15617  41654   6266 



#case/人年
fit2 <- pyears(Surv(BC_difftime_new,BC_total)~1,data=female[female$BC_PRS_FH==1,],scale = 1)
fit2$pyears#
fit2$event#
#round(866/474664.3*100000,2)=182.44

fit2 <- pyears(Surv(BC_difftime_new,BC_total)~1,data=female[female$BC_PRS_FH==2,],scale = 1)
fit2$pyears#
fit2$event#  
#round(142/46947.72*100000,2)=302.46

fit2 <- pyears(Surv(BC_difftime_new,BC_total)~1,data=female[female$BC_PRS_FH==3,],scale = 1)
fit2$pyears#
fit2$event#   
#round(4014/1386300*100000,2)=289.55

fit2 <- pyears(Surv(BC_difftime_new,BC_total)~1,data=female[female$BC_PRS_FH==4,],scale = 1)
fit2$pyears#
fit2$event#  
#round(744/167271.9*100000,2)=444.78

fit2 <- pyears(Surv(BC_difftime_new,BC_total)~1,data=female[female$BC_PRS_FH==5,],scale = 1)
fit2$pyears#
fit2$event#   
#round(1987/445389.4*100000,2)=446.13

fit2 <- pyears(Surv(BC_difftime_new,BC_total)~1,data=female[female$BC_PRS_FH==6,],scale = 1)
fit2$pyears#
fit2$event#  
#round(443/66316.58*100000,2)=668.01



model1 <- coxph(Surv(BC_difftime_new,BC_total)~as.factor(BC_PRS_FH)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)
res=summary(model1,digit=10)
round(res$conf.int,2)
#                                       exp(coef) exp(-coef) lower .95 upper .95
#as.factor(BC_PRS_FH)2                    1.62       0.62      1.36      1.94
#as.factor(BC_PRS_FH)3                    1.61       0.62      1.49      1.73
#as.factor(BC_PRS_FH)4                    2.40       0.42      2.18      2.65
#as.factor(BC_PRS_FH)5                    2.52       0.40      2.33      2.73
#as.factor(BC_PRS_FH)6                    3.65       0.27      3.26      4.10


dat1=female

dat_df <- with(dat1,
               data.frame(BC_PRS_FH = c(1, 2, 3, 4, 5, 6), 
                          Age=rep(mean(Age,na.rm = TRUE),6),
                          Height_imp=rep(mean(Height_imp,na.rm = TRUE),6),
                          Townsend_deprivation_index_imp_cat=c(2,2,2,2,2,2),
                          smoking_health=rep(mean(smoking_health,na.rm = TRUE),6),
                          Alcohol_health=rep(mean(Alcohol_health,na.rm = TRUE),6),
                          BMI_health=rep(mean(BMI_health,na.rm = TRUE),6),
                          physical_health=rep(mean(physical_health,na.rm = TRUE),6),
                          PCA1=rep(mean(PCA1,na.rm = TRUE),6),
                          PCA2=rep(mean(PCA2,na.rm = TRUE),6),
                          PCA3=rep(mean(PCA3,na.rm = TRUE),6),
                          PCA4=rep(mean(PCA4,na.rm = TRUE),6),
                          PCA5=rep(mean(PCA5,na.rm = TRUE),6),
                          PCA6=rep(mean(PCA6,na.rm = TRUE),6),
                          PCA7=rep(mean(PCA7,na.rm = TRUE),6),
                          PCA8=rep(mean(PCA8,na.rm = TRUE),6),
                          PCA9=rep(mean(PCA9,na.rm = TRUE),6),
                          PCA10=rep(mean(PCA10,na.rm = TRUE),6)			  
               )
)



dat_df$BC_PRS_FH=as.factor(dat_df$BC_PRS_FH)
dat_df$Townsend_deprivation_index_imp_cat=as.factor(dat_df$Townsend_deprivation_index_imp_cat)


fit2 <- survfit(model1, newdata = dat_df)#survfit()：使用公式或已构建的Cox模型拟合生存曲线

p2=ggsurvplot(fit2,
              linetype=c(1,1,1,1,1,1),
              xlab=c("Years of Follow-up"),
              xlim=c(0,7),
              ylim=c(0,0.05),
              ylab=c("Standardized Breast cancer Event Rate"),
              legend=c(0.30,0.8),
              legend.title="FH combined PRS",
              legend.labs=c("Group1 (Reference)","Group2 (HR):1.62 (1.36-1.94)","Group3 (HR):1.61 (1.49-1.73)",
                            "Group4 (HR):2.40 (2.18-2.65)","Group5 (HR):2.52 (2.33-2.73)","Group6 (HR):3.65 (3.26-4.10)"),			
              conf.int=TRUE,
              conf.int.fill="gray",
              conf.int.style="ribbon",
              conf.int.alpha=0.2,
              risk.table.col="strata",
              palette=c("#00468B99","#42B540E5","#0099B4FF","#FDAF91E5","#925E9FE5","#AD002AFF"),
              ggtheme=theme_bw(),
              fun="event",
              font.xtickslab=c(12,"bold"),
              font.ytickslab=c(12,"bold"),
              font.x=c(14,"bold"),
              font.y=c(14,"bold"),
              font.legend=c(10,"bold"),
              censor=FALSE,
              dat=dat_df)
p2

#7.08×5.84

ggsave("Female_75_years_old_Breast_Cancer_risk_predtion_of_PRS_combined_FH.pdf", plot = p2, width = 8, height = 6, dpi = 600)



##  相加交互作用检验
#coef的2，3，4，5，6分别对应分组时的2-6组，但是却是生成系数的第1-5行

RERI1=epi.interaction(model=model1,coef=c(1,2,3),param = c("dummy"),conf.level=0.95) 
RERI1  

#$reri
#est      lower     upper
#1 0.173497 -0.1530246 0.5000187

#$apab
#est       lower     upper
#1 0.07217267 -0.06189515 0.2062405

#$s
#est     lower    upper
#1 1.141006 0.8798602 1.479662

#$multiplicative
#est    lower    upper
#1 2.403916 2.179493 2.651447


RERI2=epi.interaction(model=model1,coef=c(1,4,5),param = c("dummy"),conf.level=0.95)
RERI2

#$reri
#est      lower     upper
#1 0.5075735 0.05976357 0.9553835

#$apab
#est      lower     upper
#1 0.138891 0.02475991 0.2530221

#$s
#est   lower    upper
#1 1.236422 1.02287 1.494557

#$multiplicative
#est    lower    upper
#1 3.654474 3.258972 4.097974
