###sfigure 4 & sfigure10###

#75岁时不同风险人群的肿瘤发病

#加载数据
load("imputed_final_20240120_UKB_database_for_analysis_202801_without_prevalent_cancer_male.rdata")
load("imputed_final_20240120_UKB_database_for_analysis_239598_without_prevalent_cancer_female.rdata")

#加载R包
library("ggplot2")
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
library("Cairo") 


#sfigure 4

#计算出现结局的年龄_cancer_age   

#CRC_total.y  是否出现肿瘤结局
#CRC_endpoint.y  随访截止时间
#birth_date   出生日期

########肠癌 CRC COL_PRS   COL_PRS

#male#

#更新男性出现肠癌结局的年龄
male$CRC_cancer_age=as.numeric(as.Date(male$CRC_endpoint.y)-as.Date(male$birth_date))/365.25
summary(male$CRC_cancer_age)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#40.38   61.04   68.67   67.57   74.12   85.96 

control <- male[which(male$FH_Bowel_2==0),]
summary(control$COL_PRS)
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

male$COL_PRS_FH=0
male[male$PRS_1_sd_3_Q20==1&male$FH_Bowel_2==0,]$COL_PRS_FH=1  #对照（1，2，3，4）（1，2，5，6）
male[male$PRS_1_sd_3_Q20==1&male$FH_Bowel_2==1,]$COL_PRS_FH=2  #家族史效应
male[male$PRS_1_sd_3_Q20==2&male$FH_Bowel_2==0,]$COL_PRS_FH=3  #PRS中位效应
male[male$PRS_1_sd_3_Q20==2&male$FH_Bowel_2==1,]$COL_PRS_FH=4  #PRS中位联合家族史效应
male[male$PRS_1_sd_3_Q20==3&male$FH_Bowel_2==0,]$COL_PRS_FH=5  #PRS高位效应
male[male$PRS_1_sd_3_Q20==3&male$FH_Bowel_2==1,]$COL_PRS_FH=6  #PRS高位联合家族史效应
table(male$COL_PRS_FH)
male$COL_PRS_FH=as.factor(male$COL_PRS_FH)
#1      2      3      4      5      6 
#36754   3807 108173  13506  35121   5440 


##########male:年龄 vs 肿瘤累积风险图##########

dat1=male
model2 <-coxph(Surv(CRC_cancer_age,CRC_total.y)~COL_PRS_FH+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=dat1)
summary((model2),data=dat1)


dat_df <- with(dat1,
               data.frame(CRC_cancer_age=c(75,75,75,75,75,75),
                          CRC_total.y=c(1,1,1,1,1,1),
                          COL_PRS_FH = c(1, 2, 3,4,5,6), 
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

fit2 <- survfit(model2, newdata = dat_df)

fit2_1 <- predict(model2, newdata = dat_df, type="expected", se=T)

fit2_1 #by 75 age#
#$fit
#[1] 0.02226429 0.02967297 0.04136553 0.05671569 0.07120192 0.08211633

#$se.fit
#[1] 0.001625680 0.004873778 0.001848845 0.004056476 0.003619276 0.007355179

x=c(0.02226429, 0.02967297, 0.04136553, 0.05671569, 0.07120192, 0.08211633)
y=c(0.001625680, 0.004873778, 0.001848845, 0.004056476, 0.003619276, 0.007355179)

X=round(x*100,2)
X
#2.23 2.97 4.14 5.67 7.12 8.21
Y=round(y*100,2)
Y

LCI=X-1.96*Y
HCI=X+1.96*Y

LCI
LCI=round(LCI,2)
LCI
#1.92 2.01 3.79 4.87 6.41 6.76

HCI
HCI=round(HCI,2)
HCI
#2.54 3.93 4.49 6.47 7.83 9.66


ggsurvplot(fit2,
              linetype=c(1,1,1,1,1,1),
              xlab=c("Age (Years)"),
              xlim=c(50,75),
              ylim=c(0,0.15),
              break.x.by=10,#X轴的刻度
              ylab=c("Cumulative Cancer Event Rate"),
              legend=c(0.35,0.8),  #设置标签位置
              legend.title="FH combined PRS",
              legend.labs=c("Group1: 2.23 (95%CI:1.92-2.54)","Group2: 2.97 (95%CI:2.01-3.93)","Group3: 4.14 (95%CI:3.79-4.49)",
                            "Group4: 5.67 (95%CI:4.87-6.47)","Group5: 7.12 (95%CI:6.41-7.83)","Group6: 8.21 (95%CI:6.76-9.66)"),			
              conf.int=TRUE,#设置曲线置信区间
              conf.int.fill="gray",
              conf.int.style="ribbon",#控制置信区间的显示方式，绘制一条带状区域
              conf.int.alpha=0.2,#置信区间颜色深浅，数值越大颜色越深
              risk.table.col="strata",
              palette=c("#00468B99","#42B540E5","#0099B4FF","#FDAF91E5","#925E9FE5","#AD002AFF"),
              ggtheme=theme_bw(),
              fun="cumhaz",
              font.xtickslab=c(14,"bold"),
              font.ytickslab=c(14,"bold"),
              font.x=c(16,"bold"),
              font.y=c(16,"bold"),
              font.legend=c(14,"bold"),  #设置标签字体大小
              censor=FALSE,
              dat=dat_df)



ggsave("Male_75_years_old_lifetime_risk_Colorectal_Cancer_risk_predtion_of_PRS_combined_FH.pdf",height=8,width=8)



#female#

#更新女性出现肠癌结局的年龄

female$CRC_cancer_age=as.numeric(as.Date(female$CRC_endpoint.y)-as.Date(female$birth_date))/365.25
summary(female$CRC_cancer_age)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#41.17   61.04   68.38   67.47   73.88   83.71 

control <- female[which(female$FH_Bowel_2==0),]
summary(control$COL_PRS)
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
female$COL_PRS_FH=as.factor(female$COL_PRS_FH)
#1      2      3      4      5      6 
#43725   4195 128413  15345  41662   6258 


##########male:年龄 vs 肿瘤累积风险图##########
dat1=female
model2 <-coxph(Surv(CRC_cancer_age,CRC_total.y)~COL_PRS_FH+Age+Townsend_deprivation_index_imp_cat+Height_imp+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=dat1)
summary((model2),data=dat1)

#all#
dat_df <- with(dat1,
               data.frame(CRC_cancer_age=c(75,75,75,75,75,75),
                          CRC_total.y=c(1,1,1,1,1,1),
                          COL_PRS_FH = c(1, 2, 3,4,5,6), 
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

fit2 <- survfit(model2, newdata = dat_df)

fit2_1 <- predict(model2, newdata = dat_df, type="expected", se=T)

fit2_1 #by 75 age#
#$fit
#[1] 0.01599919 0.02077227 0.02884814 0.03924148 0.05399645 0.05297023

#$se.fit
#[1] 0.001321546 0.004091352 0.001477996 0.003319724 0.003086993 0.005697846

x=c(0.01599919, 0.02077227, 0.02884814, 0.03924148, 0.05399645, 0.05297023)
y=c(0.001321546, 0.004091352, 0.001477996, 0.003319724, 0.003086993, 0.005697846)

X=round(x*100,2)
X
#1.60 2.08 2.88 3.92 5.40 5.30

Y=round(y*100,2)
Y

LCI=X-1.96*Y
HCI=X+1.96*Y

LCI=round(LCI,2)
LCI
#1.35 1.28 2.59 3.27 4.79 4.18

HCI=round(HCI,2)
HCI
#1.85 2.88 3.17 4.57 6.01 6.42



ggsurvplot(fit2,
              linetype=c(1,1,1,1,1,1),
              xlab=c("Age (Years)"),
              xlim=c(50,75),
              ylim=c(0,0.10),
              break.x.by=10,
              ylab=c("Cumulative Cancer Event Rate"),
              legend=c(0.35,0.8),
              legend.title="FH combined PRS",
              legend.labs=c("Group1: 1.60 (95%CI:1.35-1.85)","Group2: 2.08 (95%CI:1.28-2.88)","Group3: 2.88 (95%CI:2.59-3.17)",
                            "Group4: 3.92 (95%CI:3.27-4.57)","Group5: 5.40 (95%CI:4.79-6.01)","Group6: 5.30 (95%CI:4.18-6.42)"),			
              conf.int=TRUE,
              conf.int.fill="gray",
              conf.int.style="ribbon",
              conf.int.alpha=0.2,
              risk.table.col="strata",
              palette=c("#00468B99","#42B540E5","#0099B4FF","#FDAF91E5","#925E9FE5","#AD002AFF"),
              ggtheme=theme_bw(),
              fun="cumhaz",
              font.xtickslab=c(14,"bold"),
              font.ytickslab=c(14,"bold"),
              font.x=c(16,"bold"),
              font.y=c(16,"bold"),
              font.legend=c(14,"bold"),  #设置标签字体大小
              censor=FALSE,
              dat=dat_df)

ggsave("Female_75_years_old_lifetime_risk_Colorectal_Cancer_risk_predtion_of_PRS_combined_FH.pdf",height=8,width=8)






########肺癌 LC LC_PRS########

#male

#更新男性出现肺癌结局的年龄
male$LC_cancer_age=as.numeric(as.Date(male$LC_endpoint.y)-as.Date(male$birth_date))/365.25
summary(male$LC_cancer_age)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#40.38   61.04   68.71   67.61   74.21   85.96

control <- male[which(male$FH_Lung_2==0),]
summary(control$LC_PRS)
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
male$LC_PRS_FH=as.factor(male$LC_PRS_FH)
#1      2      3      4      5      6 
#36083   4478 107075  14604  34988   5573 


##########male:年龄 vs 肿瘤累积风险图##########
dat1=male
model2 <-coxph(Surv(LC_cancer_age,LC_total.y)~LC_PRS_FH+Age+Townsend_deprivation_index_imp_cat+Height_imp+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=dat1)
summary((model2),data=dat1)

#all#
dat_df <- with(dat1,
               data.frame(LC_cancer_age=c(75,75,75,75,75,75),
                          LC_total.y=c(1,1,1,1,1,1),
                          LC_PRS_FH = c(1, 2, 3,4,5,6), 
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

fit2 <- survfit(model2, newdata = dat_df)

fit2_1 <- predict(model2, newdata = dat_df, type="expected", se=T)

fit2_1 #by 75 age#
#$fit
#[1] 0.009105593 0.015241196 0.014964134 0.024668183 0.020513240 0.028138763

#$se.fit
#[1] 0.0008610621 0.0024267358 0.0010056731 0.0021843971 0.0015828377 0.0032672129

x=c(0.009105593, 0.015241196, 0.014964134, 0.024668183, 0.020513240, 0.028138763)
y=c(0.0008610621, 0.0024267358, 0.0010056731, 0.0021843971, 0.0015828377, 0.0032672129)

X=round(x*100,2)
X
#0.91 1.52 1.50 2.47 2.05 2.81
Y=round(y*100,2)
Y

LCI=X-1.96*Y
HCI=X+1.96*Y

LCI
LCI=round(LCI,2)
LCI
#0.73 1.05 1.30 2.04 1.74 2.16

HCI
HCI=round(HCI,2)
HCI
#1.09 1.99 1.70 2.90 2.36 3.46



ggsurvplot(fit2,
              linetype=c(1,1,1,1,1,1),
              xlab=c("Age (Years)"),
              xlim=c(50,75),
              ylim=c(0,0.05),
              break.x.by=10,
              ylab=c("Cumulative Cancer Event Rate"),
              legend=c(0.35,0.8),
              legend.title="FH combined PRS",
              legend.labs=c("Group1: 0.91 (95%CI:0.73-1.09)","Group2: 1.52 (95%CI:1.05-1.99)","Group3: 1.50 (95%CI:1.30-1.70)",
                            "Group4: 2.47 (95%CI:2.04-2.90)","Group5: 2.05 (95%CI:1.74-2.36)","Group6: 2.81 (95%CI:2.16-3.46)"),			
              conf.int=TRUE,
              conf.int.fill="gray",
              conf.int.style="ribbon",
              conf.int.alpha=0.2,
              risk.table.col="strata",
              palette=c("#00468B99","#42B540E5","#0099B4FF","#FDAF91E5","#925E9FE5","#AD002AFF"),
              ggtheme=theme_bw(),
              fun="cumhaz",
              font.xtickslab=c(14,"bold"),
              font.ytickslab=c(14,"bold"),
              font.x=c(16,"bold"),
              font.y=c(16,"bold"),
              font.legend=c(14,"bold"),  #设置标签字体大小
              censor=FALSE,
              dat=dat_df)

ggsave("Male_75_years_old_lifetime_risk_Lung_Cancer_risk_predtion_of_PRS_combined_FH.pdf",height=8,width=8)



###female###

#更新女性出现肺癌结局的年龄

female$LC_cancer_age=as.numeric(as.Date(female$LC_endpoint.y)-as.Date(female$birth_date))/365.25
summary(female$LC_cancer_age)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#41.17   61.12   68.38   67.49   73.88   83.71

control <- female[which(female$FH_Lung_2==0),]
summary(control$LC_PRS)
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
female$LC_PRS_FH=as.factor(female$LC_PRS_FH)
#1      2      3      4      5      6 
#42427   5493 125539  18219  41283   6637 



##########female:年龄 vs 肿瘤累积风险图##########
dat1=female
model2 <-coxph(Surv(LC_cancer_age,LC_total.y)~LC_PRS_FH+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=dat1)
summary((model2),data=dat1)

#all#
dat_df <- with(dat1,
               data.frame(LC_cancer_age=c(75,75,75,75,75,75),
                          LC_total.y=c(1,1,1,1,1,1),
                          LC_PRS_FH = c(1, 2, 3,4,5,6), 
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

fit2 <- survfit(model2, newdata = dat_df)

fit2_1 <- predict(model2, newdata = dat_df, type="expected", se=T)

fit2_1 #by 75 age#
#$fit
#[1] 0.01000375 0.01392801 0.01431521 0.02339344 0.01840344 0.02974707

#$se.fit
#[1] 0.0009224575 0.0022796032 0.0009552191 0.0020561787 0.0014525365 0.0033969697

x=c(0.01000375, 0.01392801, 0.01431521, 0.02339344, 0.01840344, 0.02974707)
y=c(0.0009224575, 0.0022796032, 0.0009552191, 0.0020561787, 0.0014525365, 0.0033969697)

X=round(x*100,2)
X
#1.00 1.39 1.43 2.34 1.84 2.97
Y=round(y*100,2)
Y

LCI=X-1.96*Y
HCI=X+1.96*Y

LCI=round(LCI,2)
LCI
#0.82 0.94 1.23 1.93 1.55 2.30

HCI=round(HCI,2)
HCI
#1.18 1.84 1.63 2.75 2.13 3.64

ggsurvplot(fit2,
              linetype=c(1,1,1,1,1,1),
              xlab=c("Age (Years)"),
              xlim=c(50,75),
              ylim=c(0,0.05),
              break.x.by=10,
              ylab=c("Cumulative Cancer Event Rate"),
              legend=c(0.35,0.8),
              legend.title="FH combined PRS",
              legend.labs=c("Group1: 1.00 (95%CI:0.82-1.18)","Group2: 1.39 (95%CI:0.94-1.84)","Group3: 1.43 (95%CI:1.23-1.63)",
                            "Group4: 2.34 (95%CI:1.93-2.75)","Group5: 1.84 (95%CI:1.55-2.13)","Group6: 2.97 (95%CI:2.30-3.64)"),			
              conf.int=TRUE,
              conf.int.fill="gray",
              conf.int.style="ribbon",
              conf.int.alpha=0.2,
              risk.table.col="strata",
              palette=c("#00468B99","#42B540E5","#0099B4FF","#FDAF91E5","#925E9FE5","#AD002AFF"),
              ggtheme=theme_bw(),
              fun="cumhaz",
              font.xtickslab=c(14,"bold"),
              font.ytickslab=c(14,"bold"),
              font.x=c(16,"bold"),
              font.y=c(16,"bold"),
              font.legend=c(14,"bold"),  #设置标签字体大小
              censor=FALSE,
              dat=dat_df)


ggsave("Female_75_years_old_lifetime_risk_Lung_Cancer_risk_predtion_of_PRS_combined_FH.pdf",height=8,width=8)






########前列腺癌 PRC PRO_PRS ########

#male

#更新男性出现前列腺癌结局的年龄
male$PRC_cancer_age=as.numeric(as.Date(male$PRC_endpoint.y)-as.Date(male$birth_date))/365.25
summary(male$PRC_cancer_age)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#40.38   60.96   68.38   67.40   73.88   85.96 

control <- male[which(male$FH_Prostate_2==0),]
summary(control$PRO_PRS)
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
male$PRO_PRS_FH=as.factor(male$PRO_PRS_FH)
#1      2      3      4      5      6 
#38147   2414 112266   9413  36615   3946



##########male:年龄 vs 肿瘤累积风险图##########
dat1=male
model2 <-coxph(Surv(PRC_cancer_age,PRC_total.y)~PRO_PRS_FH+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=dat1)
summary((model2),data=dat1)

#all#
dat_df <- with(dat1,
               data.frame(PRC_cancer_age=c(75,75,75,75,75,75),
                          PRC_total.y=c(1,1,1,1,1,1),
                          PRO_PRS_FH = c(1, 2, 3,4,5,6), 
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

fit2 <- survfit(model2, newdata = dat_df)

fit2_1 <- predict(model2, newdata = dat_df, type="expected", se=T)

fit2_1 #by 75 age#
#$fit
#[1] 0.06154390 0.07846314 0.14684051 0.24132957 0.34419727 0.54465537

#$se.fit
#[1] 0.002719718 0.010606071 0.003704176 0.010754213 0.009281848 0.025853132

x=c(0.06154390, 0.07846314, 0.14684051, 0.24132957, 0.34419727, 0.54465537)
y=c(0.002719718, 0.010606071, 0.003704176, 0.010754213, 0.009281848, 0.025853132)

X=round(x*100,2)
X
#6.15  7.85 14.68 24.13 34.42 54.47
Y=round(y*100,2)
Y

LCI=X-1.96*Y
HCI=X+1.96*Y

LCI
LCI=round(LCI,2)
LCI
#5.62  5.77 13.95 22.01 32.60 49.39

HCI
HCI=round(HCI,2)
HCI
#6.68  9.93 15.41 26.25 36.24 59.55



ggsurvplot(fit2,
              linetype=c(1,1,1,1,1,1),
              xlab=c("Age (Years)"),
              xlim=c(50,75),
              ylim=c(0,0.65),
              break.x.by=10,
              ylab=c("Cumulative Cancer Event Rate"),
              legend=c(0.35,0.8),
              legend.title="FH combined PRS",
              legend.labs=c("Group1: 6.15 (95%CI:5.62-6.68)","Group2: 7.85 (95%CI:5.77-9.93)","Group3: 14.68 (95%CI:13.95-15.41)",
                            "Group4: 24.13 (95%CI:22.01-26.25)","Group5: 34.42 (95%CI:32.60-36.24)","Group6: 54.47 (95%CI:49.39-59.55)"),			
              conf.int=TRUE,
              conf.int.fill="gray",
              conf.int.style="ribbon",
              conf.int.alpha=0.2,
              risk.table.col="strata",
              palette=c("#00468B99","#42B540E5","#0099B4FF","#FDAF91E5","#925E9FE5","#AD002AFF"),
              ggtheme=theme_bw(),
              fun="cumhaz",
              font.xtickslab=c(14,"bold"),
              font.ytickslab=c(14,"bold"),
              font.x=c(16,"bold"),
              font.y=c(16,"bold"),
              font.legend=c(14,"bold"),  #设置标签字体大小
              censor=FALSE,
              dat=dat_df)


ggsave("Male_75_years_old_lifetime_risk_Prostate_Cancer_risk_predtion_of_PRS_combined_FH.pdf",height = 8,width = 8)





#########female  BRC  BC_PRS#########

#更新女性出现乳腺癌结局的年龄

female$BRC_cancer_age=as.numeric(as.Date(female$BC_endpoint)-as.Date(female$birth_date))/365.25
summary(female$BRC_cancer_age)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#40.78   60.96   68.21   67.32   73.72   83.71 

control <- female[which(female$FH_Breast_2==0),]
summary(control$BC_PRS)
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
female$BC_PRS_FH=as.factor(female$BC_PRS_FH)
#1      2      3      4      5      6 
#43578   4342 128141  15617  41654   6266



##########female:年龄 vs 肿瘤累积风险图##########

#by 75 age#
#$fit
#[1] 0.06629839 0.10753745 0.10655944 0.15875790 0.16748440 0.24307692

#$se.fit
#[1] 0.002597417 0.009269198 0.002668938 0.006604128 0.004978175 0.012501546

x=c(0.06629839, 0.10753745, 0.10655944, 0.15875790, 0.16748440, 0.24307692)
y=c(0.002597417, 0.009269198, 0.002668938, 0.006604128, 0.004978175, 0.012501546)

X=round(x*100,2)
X
#6.63 10.75 10.66 15.88 16.75 24.31
Y=round(y*100,2)
Y

LCI=X-1.96*Y
HCI=X+1.96*Y

LCI=round(LCI,2)
LCI
# 6.12  8.93 10.13 14.59 15.77 21.86

HCI=round(HCI,2)
HCI
#7.14 12.57 11.19 17.17 17.73 26.76



ggsurvplot(fit2,
              linetype=c(1,1,1,1,1,1),
              xlab=c("Age (Years)"),
              xlim=c(50,75),
              ylim=c(0,0.3),
              break.x.by=10,
              ylab=c("Cumulative Cancer Event Rate"),
              legend=c(0.35,0.8),
              legend.title="FH combined PRS",
              legend.labs=c("Group1: 6.63 (95%CI:6.12-7.14)","Group2: 10.75 (95%CI:8.93-12.57)","Group3: 10.66 (95%CI:10.13-11.19)",
                            "Group4: 15.88 (95%CI:14.59-17.17)","Group5: 16.75 (95%CI:15.77-17.73)","Group6: 24.31 (95%CI:21.86-26.76)"),			
              conf.int=TRUE,
              conf.int.fill="gray",
              conf.int.style="ribbon",
              conf.int.alpha=0.2,
              risk.table.col="strata", 
              palette=c("#00468B99","#42B540E5","#0099B4FF","#FDAF91E5","#925E9FE5","#AD002AFF"),
              ggtheme=theme_bw(),
              fun="cumhaz",
              font.xtickslab=c(14,"bold"),
              font.ytickslab=c(14,"bold"),
              font.x=c(16,"bold"),
              font.y=c(16,"bold"),
              font.legend=c(14,"bold"),  #设置标签字体大小
              censor=FALSE,
              dat=dat_df)


ggsave("Female_75_years_old_lifetime_risk_Breast_Cancer_risk_predtion_of_PRS_combined_FH.pdf",height = 8,width = 8)



#sfigure 10

####家族史联合CPRS预测75岁时全肿瘤发生风险########

#male
male$Overall_cancer_age=as.numeric(as.Date(male$cancer_endpoint_20)-as.Date(male$birth_date))/365.25
summary(male$Overall_cancer_age)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#40.38   60.79   68.12   67.24   73.71   85.96

#male    
table(male$FH_total_2)
#0      1 
#146524  56277 

male$cprs=male$CPRS
control <- male[which(male$FH_total_2==0),]
summary(control$cprs)
sd(control$cprs)  ##0.1988989
male$cprs_sd<-male$cprs/0.1988989

Q20=quantile(male$cprs_sd,seq(0.05,1,0.05))[[4]]
Q80=quantile(male$cprs_sd,seq(0.05,1,0.05))[[16]]

####对原始数据的PRS三分位，<=Q20，1；>=Q80，3
male$cprs_sd_3=0
male[male$cprs_sd<=Q20,]$cprs_sd_3=1
male[male$cprs_sd>Q20&male$cprs_sd<Q80,]$cprs_sd_3=2
male[male$cprs_sd>=Q80,]$cprs_sd_3=3

table(male$cprs_sd_3)
#1      2      3 
#40561 121679  40561 

male$cprs_FH=0
male[male$cprs_sd_3==1&male$FH_total_2==0,]$cprs_FH=1  #对照（1，2，3，4）（1，2，5，6）
male[male$cprs_sd_3==1&male$FH_total_2==1,]$cprs_FH=2  #家族史效应
male[male$cprs_sd_3==2&male$FH_total_2==0,]$cprs_FH=3  #CPRS中位效应
male[male$cprs_sd_3==2&male$FH_total_2==1,]$cprs_FH=4  #CPRS中位联合家族史效应
male[male$cprs_sd_3==3&male$FH_total_2==0,]$cprs_FH=5  #CPRS高位效应
male[male$cprs_sd_3==3&male$FH_total_2==1,]$cprs_FH=6  #CPRS高位联合家族史效应
table(male$cprs_FH)
male$cprs_FH=as.factor(male$cprs_FH)
#1      2      3      4      5      6 
#30414 10147 87988 33691 28122 12439


##########male:年龄 vs 肿瘤累积风险图##########
dat1=male   #coxph（出现结局的年龄,是否出现结局）
model2 <-coxph(Surv(Overall_cancer_age,cancer_total_20)~cprs_FH+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=dat1)
summary((model2),data=dat1)

#all#
dat_df <- with(dat1,
               data.frame(Overall_cancer_age=c(75,75,75,75,75,75),
                          cancer_total_20=c(1,1,1,1,1,1),
                          cprs_FH = c(1, 2, 3,4,5,6), 
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


dat_df$cprs_FH=as.factor(dat_df$cprs_FH)
dat_df$Townsend_deprivation_index_imp_cat=as.factor(dat_df$Townsend_deprivation_index_imp_cat)

fit2 <- survfit(model2, newdata = dat_df)

fit2_1 <- predict(model2, newdata = dat_df, type="expected", se=T)


fit2_1  #75岁全肿瘤发生率
#$fit
#[1] 0.2329762 0.2712979 0.3368872 0.3695155 0.5400025 0.6158675

#$se.fit
#[1] 0.006326544 0.010281955 0.006058584 0.008101597 0.011604129 0.016195489


x=c(0.2329762, 0.2712979, 0.3368872, 0.3695155, 0.5400025, 0.6158675)
y=c(0.006326544, 0.010281955, 0.006058584, 0.008101597, 0.011604129, 0.016195489)

X=round(x*100,2)
X
#23.30 27.13 33.69 36.95 54.00 61.59
Y=round(y*100,2)
Y

LCI=X-1.96*Y
HCI=X+1.96*Y

LCI=round(LCI,2)
LCI
# 22.07 25.11 32.49 35.36 51.73 58.41

HCI=round(HCI,2)
HCI
#24.53 29.15 34.89 38.54 56.27 64.77


ggsurvplot(fit2,
              linetype=c(1,1,1,1,1,1),
              xlab=c("Age (Years)"),
              xlim=c(50,75),
              ylim=c(0,0.7),
              break.x.by=10,
              ylab=c("Cumulative Cancer Event Rate"),
              legend=c(0.35,0.8),
              legend.title="FH of multiple cancer combined CPRS",
              legend.labs=c("Group1: 23.30 (95%CI:22.07-24.53)","Group2: 27.13 (95%CI:25.11-29.15)","Group3: 33.69 (95%CI:32.49-34.89)",
                            "Group4: 36.95 (95%CI:35.36-38.54)","Group5: 54.00 (95%CI:51.73-56.27)","Group6: 61.59 (95%CI:58.41-64.77)"),			
              conf.int=TRUE,
              conf.int.fill="gray",
              conf.int.style="ribbon",
              conf.int.alpha=0.2,
              risk.table.col="strata", 
              palette=c("#00468B99","#42B540E5","#0099B4FF","#FDAF91E5","#925E9FE5","#AD002AFF"),
              ggtheme=theme_bw(),
              fun="cumhaz",
              font.xtickslab=c(14,"bold"),
              font.ytickslab=c(14,"bold"),
              font.x=c(16,"bold"),
              font.y=c(16,"bold"),
              font.legend=c(14,"bold"),  #设置标签字体大小
              censor=FALSE,
              dat=dat_df)


ggsave("Male_75_years_old_lifetime_overall_Cancer_risk_predtion_of_PRS_combined_FH_of_multiple_cancer.pdf",height = 8,width = 8)





###female###
female$Overall_cancer_age=as.numeric(as.Date(female$cancer_endpoint_20)-as.Date(female$birth_date))/365.25
summary(female$Overall_cancer_age)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#40.78   60.79   67.96   67.17   73.63   83.71

table(female$FH_total_2)
#0      1 
#167534  72064

female$cprs=female$CPRS
control <- female[which(female$FH_total_2==0),]
summary(control$cprs)
sd(control$cprs)  ##0.0782589
female$cprs_sd<-female$cprs/0.0782589

Q20=quantile(female$cprs_sd,seq(0.05,1,0.05))[[4]]
Q80=quantile(female$cprs_sd,seq(0.05,1,0.05))[[16]]

####对原始数据的PRS三分位，<=Q20，1；>=Q80，3
female$cprs_sd_3=0
female[female$cprs_sd<=Q20,]$cprs_sd_3=1
female[female$cprs_sd>Q20&female$cprs_sd<Q80,]$cprs_sd_3=2
female[female$cprs_sd>=Q80,]$cprs_sd_3=3

table(female$cprs_sd_3)
#1      2      3 
#47920 143758  47920

female$cprs_FH=0
female[female$cprs_sd_3==1&female$FH_total_2==0,]$cprs_FH=1  #对照（1，2，3，4）（1，2，5，6）
female[female$cprs_sd_3==1&female$FH_total_2==1,]$cprs_FH=2  #家族史效应
female[female$cprs_sd_3==2&female$FH_total_2==0,]$cprs_FH=3  #CPRS中位效应
female[female$cprs_sd_3==2&female$FH_total_2==1,]$cprs_FH=4  #CPRS中位联合家族史效应
female[female$cprs_sd_3==3&female$FH_total_2==0,]$cprs_FH=5  #CPRS高位效应
female[female$cprs_sd_3==3&female$FH_total_2==1,]$cprs_FH=6  #CPRS高位联合家族史效应
table(female$cprs_FH)
female$cprs_FH=as.factor(female$cprs_FH)
#1      2      3      4      5      6 
#35545 12375 99936 43822 32053 15867


##########male:年龄 vs 肿瘤累积风险图##########
dat1=female    #coxph（出现结局的年龄,是否出现结局）
model2 <-coxph(Surv(Overall_cancer_age,cancer_total_20)~cprs_FH+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=dat1)
summary((model2),data=dat1)

#all#
dat_df <- with(dat1,
               data.frame(Overall_cancer_age=c(75,75,75,75,75,75),
                          cancer_total_20=c(1,1,1,1,1,1),
                          cprs_FH = c(1, 2, 3,4,5,6), 
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


dat_df$cprs_FH=as.factor(dat_df$cprs_FH)
dat_df$Townsend_deprivation_index_imp_cat=as.factor(dat_df$Townsend_deprivation_index_imp_cat)

fit2 <- survfit(model2, newdata = dat_df)

fit2_1 <- predict(model2, newdata = dat_df, type="expected", se=T)


fit2_1  #75岁全肿瘤发生率
#$fit
#[1] 0.2092829 0.2441920 0.2573067 0.2985049 0.3361797 0.3790893

#$se.fit
#[1] 0.005647017 0.008992674 0.004686606 0.006438469 0.007825995 0.010659797


x=c(0.2092829, 0.2441920, 0.2573067, 0.2985049, 0.3361797, 0.3790893)
y=c(0.005647017, 0.008992674, 0.004686606, 0.006438469, 0.007825995, 0.010659797)

X=round(x*100,2)
X
#20.93 24.42 25.73 29.85 33.62 37.91
Y=round(y*100,2)
Y

LCI=X-1.96*Y
HCI=X+1.96*Y

LCI=round(LCI,2)
LCI
#19.83 22.66 24.81 28.60 32.09 35.81

HCI=round(HCI,2)
HCI
#22.03 26.18 26.65 31.10 35.15 40.01


ggsurvplot(fit2,
              linetype=c(1,1,1,1,1,1),
              xlab=c("Age (Years)"),
              xlim=c(50,75),
              ylim=c(0,0.5),
              break.x.by=10,
              ylab=c("Cumulative Cancer Event Rate"),
              legend=c(0.35,0.8),
              legend.title="FH of multiple cancer combined CPRS",
              legend.labs=c("Group1: 20.93 (95%CI:19.83-22.03)","Group2: 24.42 (95%CI:22.66-26.18)","Group3: 25.73 (95%CI:24.81-26.65)",
                            "Group4: 29.85 (95%CI:28.60-31.10)","Group5: 33.62 (95%CI:32.09-35.15)","Group6: 37.91 (95%CI:35.81-40.01)"),			
              conf.int=TRUE,
              conf.int.fill="gray",
              conf.int.style="ribbon",
              conf.int.alpha=0.2,
              risk.table.col="strata", 
              palette=c("#00468B99","#42B540E5","#0099B4FF","#FDAF91E5","#925E9FE5","#AD002AFF"),
              ggtheme=theme_bw(),
              fun="cumhaz",
              font.xtickslab=c(14,"bold"),
              font.ytickslab=c(14,"bold"),
              font.x=c(16,"bold"),
              font.y=c(16,"bold"),
              font.legend=c(14,"bold"),  #设置标签字体大小
              censor=FALSE,
              dat=dat_df)


ggsave("Female_75_years_old_lifetime_overall_Cancer_risk_predtion_of_PRS_combined_FH_of_multiple_cancer.pdf",height = 8,width = 8)


