#肿瘤家族史与其他肿瘤发生风险   森林图

#4种肿瘤家族史与其他20肿瘤发生风险
setwd("D://2022.11.10   PRS & 家族史/20240120最新PRS/data")


#加载数据
load("imputed_final_20240120_UKB_database_for_analysis_202801_without_prevalent_cancer_male.rdata")
load("imputed_final_20240120_UKB_database_for_analysis_239598_without_prevalent_cancer_female.rdata")


#加载R包
library("survival")
library("survminer")
library("tidyverse")
library("psych")
library("reshape")
library("openxlsx")
library("data.table")
library("writexl")
library("eoffice")



#定义总结函数结果的代码
HR_CI <- function(model1) {
  a <- as.data.frame(summary(model1)$coefficients)
  P.value <- format.pval(a$`Pr(>|z|)`[1])  
  HR <- round(a$`exp(coef)`, 2)[1]
  b <- as.data.frame(summary(model1)$conf.int)[1, ]
  lower <- round((b$`lower .95`), 2)
  upper <- round((b$`upper .95`), 2)
  CI95 <- paste("(", round((b$`lower .95`), 2), "-", round((b$`upper .95`), 2), ")", sep = "")
  dat <- data.frame(P.value = (P.value), HR = c(HR), lower = lower, upper = upper, CI95 = c(CI95))
  return(dat)
}




#男性  16种
#人年 肠癌 肺癌 前列腺癌 胃癌 胰腺癌  睾丸癌  食管癌  霍奇金  多发性骨髓瘤 口腔癌  皮肤黑色素瘤  膀胱癌 甲状腺 脑神经 肾脏 淋巴细胞白血病
fit2 <- pyears(Surv(CRC_difftime_new,CRC_total.y)~1,data=male,scale = 1)
fit2$pyears#2184087
fit2$event#2838
fit2 <- pyears(Surv(LC_difftime_new,LC_total.y)~1,data=male,scale = 1)
fit2$pyears#2192761
fit2$event#1807
fit2 <- pyears(Surv(PRC_difftime_new,PRC_total.y)~1,data=male,scale = 1)
fit2$pyears#2150631
fit2$event#9585
fit2 <- pyears(Surv(GC_difftime_new,GC_total.y)~1,data=male,scale = 1)
fit2$pyears#2194663
fit2$event#421
fit2 <- pyears(Surv(PAC_difftime_new,PAC_total.y)~1,data=male,scale = 1)
fit2$pyears#2195011
fit2$event#594
fit2 <- pyears(Surv(TEC_difftime_new,TEC_total)~1,data=male,scale = 1)
fit2$pyears#2195002
fit2$event#80
fit2 <- pyears(Surv(OEC_difftime_new,ESC_total)~1,data=male,scale = 1)
fit2$pyears#2194025
fit2$event#605
fit2 <- pyears(Surv(Hodgkin_difftime_new,Hodgkin_total.y)~1,data=male,scale = 1)
fit2$pyears#2195243
fit2$event#65
fit2 <- pyears(Surv(MM_difftime_new,MM_total)~1,data=male,scale = 1)
fit2$pyears#2194001
fit2$event#448
fit2 <- pyears(Surv(LOCPC_difftime_new,LOCPC_total)~1,data=male,scale = 1)
fit2$pyears#2193258
fit2$event#582
fit2 <- pyears(Surv(SKM_difftime_new,SKM_total)~1,data=male,scale = 1)
fit2$pyears#2190486
fit2$event#1112
fit2 <- pyears(Surv(BLC_difftime_new,BLC_total.y)~1,data=male,scale = 1)
fit2$pyears#2192408
fit2$event#789
fit2 <- pyears(Surv(THC_difftime_new,THC_total)~1,data=male,scale = 1)
fit2$pyears#2195069
fit2$event#98
fit2 <- pyears(Surv(BCNC_difftime_new,BCNC_total)~1,data=male,scale = 1)
fit2$pyears#2195017
fit2$event#426
fit2 <- pyears(Surv(KIC_difftime_new,KIC_total.y)~1,data=male,scale = 1)
fit2$pyears#2192565
fit2$event#827
fit2 <- pyears(Surv(LYLE_difftime_new,LYLE_total)~1,data=male,scale = 1)
fit2$pyears#2193794
fit2$event#389


#女性 18种
#人年 肠癌 肺癌 乳腺癌 胃癌 胰腺癌  食管癌 子宫颈  霍奇金  多发性骨髓瘤 口腔癌  皮肤黑色素瘤 子宫体 膀胱癌 甲状腺 脑神经 肾脏 淋巴细胞白血病 卵巢
fit2 <- pyears(Surv(CRC_difftime_new,CRC_total.y)~1,data=female,scale = 1)
fit2$pyears#2621689
fit2$event#2167
fit2 <- pyears(Surv(LC_difftime_new,LC_total.y)~1,data=female,scale = 1)
fit2$pyears#2627379
fit2$event#1682
fit2 <- pyears(Surv(BC_difftime_new,BC_total)~1,data=female,scale = 1)
fit2$pyears#2586890
fit2$event#8196
fit2 <- pyears(Surv(GC_difftime_new,GC_total.y)~1,data=female,scale = 1)
fit2$pyears#2630309
fit2$event#195
fit2 <- pyears(Surv(PAC_difftime_new,PAC_total.y)~1,data=female,scale = 1)
fit2$pyears#2630123
fit2$event#507
fit2 <- pyears(Surv(OEC_difftime_new,ESC_total)~1,data=female,scale = 1)
fit2$pyears#2630123
fit2$event#222
fit2 <- pyears(Surv(CEUC_difftime_new,CEUC_total)~1,data=female,scale = 1)
fit2$pyears#2630197
fit2$event# 107
fit2 <- pyears(Surv(Hodgkin_difftime_new,Hodgkin_total.y)~1,data=female,scale = 1)
fit2$pyears#2630403
fit2$event#60
fit2 <- pyears(Surv(MM_difftime_new,MM_total)~1,data=female,scale = 1)
fit2$pyears#2629563
fit2$event#325
fit2 <- pyears(Surv(LOCPC_difftime_new,LOCPC_total)~1,data=female,scale = 1)
fit2$pyears#2629402
fit2$event#305
fit2 <- pyears(Surv(SKM_difftime_new,SKM_total)~1,data=female,scale = 1)
fit2$pyears#2624596
fit2$event#1189
fit2 <- pyears(Surv(COUC_difftime_new,COUC_total)~1,data=female,scale = 1)
fit2$pyears# 2624194
fit2$event# 1330
fit2 <- pyears(Surv(BLC_difftime_new,BLC_total.y)~1,data=female,scale = 1)
fit2$pyears#2629797
fit2$event#232
fit2 <- pyears(Surv(THC_difftime_new,THC_total)~1,data=female,scale = 1)
fit2$pyears#2629256
fit2$event#277
fit2 <- pyears(Surv(BCNC_difftime_new,BCNC_total)~1,data=female,scale = 1)
fit2$pyears#2630291
fit2$event#303
fit2 <- pyears(Surv(KIC_difftime_new,KIC_total.y)~1,data=female,scale = 1)
fit2$pyears#2628903
fit2$event#453
fit2 <- pyears(Surv(LYLE_difftime_new,LYLE_total)~1,data=female,scale = 1)
fit2$pyears#2629521
fit2$event# 252
fit2 <- pyears(Surv(OVC_difftime_new,OVC_total.y)~1,data=female,scale = 1)
fit2$pyears#2627518
fit2$event#870





#男性3种肿瘤家族史与其他肿瘤发生风险

#colorectal#
model1 <- coxph(Surv(LC_difftime_new,LC_total.y)~as.factor(FH_Bowel_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)

model2 <- coxph(Surv(PRC_difftime_new,PRC_total.y)~as.factor(FH_Bowel_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)

model3 <- coxph(Surv(GC_difftime_new,GC_total.y)~as.factor(FH_Bowel_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)

model4 <- coxph(Surv(PAC_difftime_new,PAC_total.y)~as.factor(FH_Bowel_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)

model5 <- coxph(Surv(TEC_difftime_new,TEC_total)~as.factor(FH_Bowel_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)

model6 <- coxph(Surv(OEC_difftime_new,ESC_total)~as.factor(FH_Bowel_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)

model7 <- coxph(Surv(Hodgkin_difftime_new,Hodgkin_total.y)~as.factor(FH_Bowel_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)

model8 <- coxph(Surv(MM_difftime_new,MM_total)~as.factor(FH_Bowel_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)

model9 <- coxph(Surv(LOCPC_difftime_new,LOCPC_total)~as.factor(FH_Bowel_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)

model10 <- coxph(Surv(SKM_difftime_new,SKM_total)~as.factor(FH_Bowel_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)

model11 <- coxph(Surv(BLC_difftime_new,BLC_total.y)~as.factor(FH_Bowel_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)

model12 <- coxph(Surv(THC_difftime_new,THC_total)~as.factor(FH_Bowel_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)

model13 <- coxph(Surv(BCNC_difftime_new,BCNC_total)~as.factor(FH_Bowel_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)

model14 <- coxph(Surv(KIC_difftime_new,KIC_total.y)~as.factor(FH_Bowel_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)

model15 <- coxph(Surv(LYLE_difftime_new,LYLE_total)~as.factor(FH_Bowel_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)

model16 <- coxph(Surv(survival_time,cancer_total_20)~as.factor(FH_Bowel_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)


summary_HR_CI <- data.frame(rbind(
  HR_CI(model1),
  HR_CI(model2),
  HR_CI(model3),
  HR_CI(model4),
  HR_CI(model5),
  HR_CI(model6),
  HR_CI(model7),
  HR_CI(model8),
  HR_CI(model9),
  HR_CI(model10),
  HR_CI(model11),
  HR_CI(model12),
  HR_CI(model13),
  HR_CI(model14),
  HR_CI(model15),
  HR_CI(model16)
))

summary_HR_CI


write.csv(summary_HR_CI,file = "male_colorectal_FH_with_other_cancer_risk_included_overall_cancer_risk_summary.csv",row.names = F)

#FDR校正
# 创建一个包含原始P值的向量
p_values <- summary_HR_CI$P.value
#[1] "0.82289" "0.0494"  "0.40453" "0.47301" "0.06907" "0.02325" "0.19572" "0.67229" "0.04332" "0.67983" "0.955"   "0.79915"
#[13] "0.29453" "0.37976" "0.48484" "0.00083"

# 使用p.adjust()函数进行FDR校正
fdr_adjusted <- p.adjust(p_values, method = "fdr")

# 输出校正后的p值
print(fdr_adjusted)

#[1] 0.8777493 0.1976000 0.7052218 0.7052218 0.2210240 0.1860000 0.5219200 0.8367138 0.1976000 0.8367138 0.9550000 0.8777493
#[13] 0.6732114 0.7052218 0.7052218 0.0132800




#lung#
model1 <- coxph(Surv(CRC_difftime_new,CRC_total.y)~as.factor(FH_Lung_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)

model2 <- coxph(Surv(PRC_difftime_new,PRC_total.y)~as.factor(FH_Lung_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)

model3 <- coxph(Surv(GC_difftime_new,GC_total.y)~as.factor(FH_Lung_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)

model4 <- coxph(Surv(PAC_difftime_new,PAC_total.y)~as.factor(FH_Lung_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)

model5 <- coxph(Surv(TEC_difftime_new,TEC_total)~as.factor(FH_Lung_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)

model6 <- coxph(Surv(OEC_difftime_new,ESC_total)~as.factor(FH_Lung_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)

model7 <- coxph(Surv(Hodgkin_difftime_new,Hodgkin_total.y)~as.factor(FH_Lung_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)

model8 <- coxph(Surv(MM_difftime_new,MM_total)~as.factor(FH_Lung_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)

model9 <- coxph(Surv(LOCPC_difftime_new,LOCPC_total)~as.factor(FH_Lung_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)

model10 <- coxph(Surv(SKM_difftime_new,SKM_total)~as.factor(FH_Lung_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)

model11 <- coxph(Surv(BLC_difftime_new,BLC_total.y)~as.factor(FH_Lung_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)

model12 <- coxph(Surv(THC_difftime_new,THC_total)~as.factor(FH_Lung_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)

model13 <- coxph(Surv(BCNC_difftime_new,BCNC_total)~as.factor(FH_Lung_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)

model14 <- coxph(Surv(KIC_difftime_new,KIC_total.y)~as.factor(FH_Lung_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)

model15 <- coxph(Surv(LYLE_difftime_new,LYLE_total)~as.factor(FH_Lung_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)

model16 <- coxph(Surv(survival_time,cancer_total_20)~as.factor(FH_Lung_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)


summary_HR_CI <- data.frame(rbind(
  HR_CI(model1),
  HR_CI(model2),
  HR_CI(model3),
  HR_CI(model4),
  HR_CI(model5),
  HR_CI(model6),
  HR_CI(model7),
  HR_CI(model8),
  HR_CI(model9),
  HR_CI(model10),
  HR_CI(model11),
  HR_CI(model12),
  HR_CI(model13),
  HR_CI(model14),
  HR_CI(model15),
  HR_CI(model16)
))

summary_HR_CI


write.csv(summary_HR_CI,file = "male_colorectal_FH_with_other_cancer_risk_included_overall_cancer_risk_summary.csv",row.names = F)


#FDR校正
# 创建一个包含原始P值的向量
p_values <- summary_HR_CI$P.value
#[1] "0.11091" "0.02514" "0.79413" "0.1439"  "0.38258" "0.00043" "0.59352" "0.36122" "0.15054" "0.53363" "0.32679" "0.73527"
#[13] "0.69098" "0.06045" "0.76514" "0.00656"

# 使用p.adjust()函数进行FDR校正
fdr_adjusted <- p.adjust(p_values, method = "fdr")

# 输出校正后的p值
print(fdr_adjusted)

# [1] 0.3440914 0.1340800 0.7941300 0.3440914 0.6121280 0.0068800 0.7913600 0.6121280 0.3440914 0.7761891 0.6121280 0.7941300
#[13] 0.7941300 0.2418000 0.7941300 0.0524800








#prostate#
model1 <- coxph(Surv(CRC_difftime_new,CRC_total.y)~as.factor(FH_Prostate_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)

model2 <- coxph(Surv(LC_difftime_new,LC_total.y)~as.factor(FH_Prostate_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)

model3 <- coxph(Surv(GC_difftime_new,GC_total.y)~as.factor(FH_Prostate_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)

model4 <- coxph(Surv(PAC_difftime_new,PAC_total.y)~as.factor(FH_Prostate_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)

model5 <- coxph(Surv(TEC_difftime_new,TEC_total)~as.factor(FH_Prostate_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)

model6 <- coxph(Surv(OEC_difftime_new,ESC_total)~as.factor(FH_Prostate_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)

model7 <- coxph(Surv(Hodgkin_difftime_new,Hodgkin_total.y)~as.factor(FH_Prostate_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)

model8 <- coxph(Surv(MM_difftime_new,MM_total)~as.factor(FH_Prostate_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)

model9 <- coxph(Surv(LOCPC_difftime_new,LOCPC_total)~as.factor(FH_Prostate_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)

model10 <- coxph(Surv(SKM_difftime_new,SKM_total)~as.factor(FH_Prostate_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)

model11 <- coxph(Surv(BLC_difftime_new,BLC_total.y)~as.factor(FH_Prostate_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)

model12 <- coxph(Surv(THC_difftime_new,THC_total)~as.factor(FH_Prostate_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)

model13 <- coxph(Surv(BCNC_difftime_new,BCNC_total)~as.factor(FH_Prostate_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)

model14 <- coxph(Surv(KIC_difftime_new,KIC_total.y)~as.factor(FH_Prostate_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)

model15 <- coxph(Surv(LYLE_difftime_new,LYLE_total)~as.factor(FH_Prostate_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)

model16 <- coxph(Surv(survival_time,cancer_total_20)~as.factor(FH_Prostate_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)


summary_HR_CI <- data.frame(rbind(
  HR_CI(model1),
  HR_CI(model2),
  HR_CI(model3),
  HR_CI(model4),
  HR_CI(model5),
  HR_CI(model6),
  HR_CI(model7),
  HR_CI(model8),
  HR_CI(model9),
  HR_CI(model10),
  HR_CI(model11),
  HR_CI(model12),
  HR_CI(model13),
  HR_CI(model14),
  HR_CI(model15),
  HR_CI(model16)
))

summary_HR_CI

write.csv(summary_HR_CI,file = "male_prostate_FH_with_other_cancer_risk_included_overall_cancer_risk_summary.csv",row.names = F)

#FDR校正
# 创建一个包含原始P值的向量
p_values <- summary_HR_CI$P.value
#[1] "0.80851"    "0.33206"    "0.12632"    "0.50782"    "0.37531"    "0.93866"    "0.18178"    "0.84615"    "0.69222"   
#[10] "0.98648"    "0.15023"    "0.81972"    "0.01465"    "0.4699"     "0.19687"    "< 2.22e-16"

# 使用p.adjust()函数进行FDR校正
fdr_adjusted <- p.adjust(p_values, method = "fdr")

# 输出校正后的p值
print(fdr_adjusted)

# [1] 0.9763269 0.8042357 0.5906100 0.8463667 0.8042357 0.9864800 0.5906100 0.9763269 0.9763269 0.9864800 0.5906100 0.9763269
#[13] 0.2197500 0.8463667 0.5906100        NA



#女性3种肿瘤家族史与其他肿瘤发生风险#

#colorectal#
model1 <- coxph(Surv(LC_difftime_new,LC_total.y)~as.factor(FH_Bowel_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)

model2 <- coxph(Surv(BC_difftime_new,BC_total)~as.factor(FH_Bowel_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)

model3 <- coxph(Surv(GC_difftime_new,GC_total.y)~as.factor(FH_Bowel_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)

model4 <- coxph(Surv(PAC_difftime_new,PAC_total.y)~as.factor(FH_Bowel_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)

model5 <- coxph(Surv(OEC_difftime_new,ESC_total)~as.factor(FH_Bowel_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)

model6 <- coxph(Surv(CEUC_difftime_new,CEUC_total)~as.factor(FH_Bowel_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)

model7 <- coxph(Surv(Hodgkin_difftime_new,Hodgkin_total.y)~as.factor(FH_Bowel_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)

model8 <- coxph(Surv(MM_difftime_new,MM_total)~as.factor(FH_Bowel_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)

model9 <- coxph(Surv(LOCPC_difftime_new,LOCPC_total)~as.factor(FH_Bowel_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)

model10 <- coxph(Surv(SKM_difftime_new,SKM_total)~as.factor(FH_Bowel_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)

model11 <- coxph(Surv(COUC_difftime_new,COUC_total)~as.factor(FH_Bowel_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)

model12 <- coxph(Surv(BLC_difftime_new,BLC_total.y)~as.factor(FH_Bowel_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)

model13 <- coxph(Surv(THC_difftime_new,THC_total)~as.factor(FH_Bowel_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)

model14 <- coxph(Surv(BCNC_difftime_new,BCNC_total)~as.factor(FH_Bowel_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)

model15 <- coxph(Surv(KIC_difftime_new,KIC_total.y)~as.factor(FH_Bowel_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)

model16 <- coxph(Surv(LYLE_difftime_new,LYLE_total)~as.factor(FH_Bowel_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)

model17 <- coxph(Surv(OVC_difftime_new,OVC_total.y)~as.factor(FH_Bowel_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)

model18 <- coxph(Surv(survival_time,cancer_total_20)~as.factor(FH_Bowel_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)


summary_HR_CI <- data.frame(rbind(
  HR_CI(model1),
  HR_CI(model2),
  HR_CI(model3),
  HR_CI(model4),
  HR_CI(model5),
  HR_CI(model6),
  HR_CI(model7),
  HR_CI(model8),
  HR_CI(model9),
  HR_CI(model10),
  HR_CI(model11),
  HR_CI(model12),
  HR_CI(model13),
  HR_CI(model14),
  HR_CI(model15),
  HR_CI(model16),
  HR_CI(model17),
  HR_CI(model18)
))

summary_HR_CI


write.csv(summary_HR_CI,file = "female_colorectal_FH_with_other_cancer_risk_included_overall_cancer_risk_summary.csv",row.names = F)


#FDR校正
# 创建一个包含原始P值的向量
p_values <- summary_HR_CI$P.value

#[1] "0.30008" "0.00227" "0.66242" "0.46463" "0.85495" "0.80629" "0.76992" "0.47864" "0.85482" "0.41948" "0.29428" "0.62103"
#[13] "0.22994" "0.47412" "0.02992" "0.83479" "0.07714" "0"

# 使用p.adjust()函数进行FDR校正
fdr_adjusted <- p.adjust(p_values, method = "fdr")

# 输出校正后的p值  即q值
print(fdr_adjusted)

# [1] 0.7716343 0.0204300 0.8549500 0.7832291 0.8549500 0.8549500 0.8549500 0.7832291 0.8549500 0.7832291 0.7716343 0.8549500
#[13] 0.7716343 0.7832291 0.1795200 0.8549500 0.3471300 0.0000000




#lung#
model1 <- coxph(Surv(CRC_difftime_new,CRC_total.y)~as.factor(FH_Lung_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)

model2 <- coxph(Surv(BC_difftime_new,BC_total)~as.factor(FH_Lung_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)

model3 <- coxph(Surv(GC_difftime_new,GC_total.y)~as.factor(FH_Lung_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)

model4 <- coxph(Surv(PAC_difftime_new,PAC_total.y)~as.factor(FH_Lung_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)

model5 <- coxph(Surv(OEC_difftime_new,ESC_total)~as.factor(FH_Lung_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)

model6 <- coxph(Surv(CEUC_difftime_new,CEUC_total)~as.factor(FH_Lung_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)

model7 <- coxph(Surv(Hodgkin_difftime_new,Hodgkin_total.y)~as.factor(FH_Lung_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)

model8 <- coxph(Surv(MM_difftime_new,MM_total)~as.factor(FH_Lung_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)

model9 <- coxph(Surv(LOCPC_difftime_new,LOCPC_total)~as.factor(FH_Lung_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)

model10 <- coxph(Surv(SKM_difftime_new,SKM_total)~as.factor(FH_Lung_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)

model11 <- coxph(Surv(COUC_difftime_new,COUC_total)~as.factor(FH_Lung_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)

model12 <- coxph(Surv(BLC_difftime_new,BLC_total.y)~as.factor(FH_Lung_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)

model13 <- coxph(Surv(THC_difftime_new,THC_total)~as.factor(FH_Lung_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)

model14 <- coxph(Surv(BCNC_difftime_new,BCNC_total)~as.factor(FH_Lung_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)

model15 <- coxph(Surv(KIC_difftime_new,KIC_total.y)~as.factor(FH_Lung_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)

model16 <- coxph(Surv(LYLE_difftime_new,LYLE_total)~as.factor(FH_Lung_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)

model17 <- coxph(Surv(OVC_difftime_new,OVC_total.y)~as.factor(FH_Lung_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)

model18 <- coxph(Surv(survival_time,cancer_total_20)~as.factor(FH_Lung_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)


summary_HR_CI <- data.frame(rbind(
  HR_CI(model1),
  HR_CI(model2),
  HR_CI(model3),
  HR_CI(model4),
  HR_CI(model5),
  HR_CI(model6),
  HR_CI(model7),
  HR_CI(model8),
  HR_CI(model9),
  HR_CI(model10),
  HR_CI(model11),
  HR_CI(model12),
  HR_CI(model13),
  HR_CI(model14),
  HR_CI(model15),
  HR_CI(model16),
  HR_CI(model17),
  HR_CI(model18)
))

summary_HR_CI

write.csv(summary_HR_CI,file = "female_lung_FH_with_other_cancer_risk_included_overall_cancer_risk_summary.csv",row.names = F)

#FDR校正
# 创建一个包含原始P值的向量
p_values <- summary_HR_CI$P.value

#[1] "0.87151" "0.06084" "0.76955" "0.54315" "0.31297" "0.65954" "0.78788" "0.81084" "0.96461" "0.38346" "0.52272" "0.41321"
#[13] "0.75716" "0.83421" "0.86302" "0.67261" "0.5476"  "0.00044"

# 使用p.adjust()函数进行FDR校正
fdr_adjusted <- p.adjust(p_values, method = "fdr")

# 输出校正后的p值  即q值
print(fdr_adjusted)

# [1] 0.9227753 0.5475600 0.9227753 0.9227753 0.9227753 0.9227753 0.9227753 0.9227753 0.9646100 0.9227753 0.9227753 0.9227753
#[13] 0.9227753 0.9227753 0.9227753 0.9227753 0.9227753 0.0079200







#breast#
model1 <- coxph(Surv(CRC_difftime_new,CRC_total.y)~as.factor(FH_Breast_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)

model2 <- coxph(Surv(LC_difftime_new,LC_total.y)~as.factor(FH_Breast_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)

model3 <- coxph(Surv(GC_difftime_new,GC_total.y)~as.factor(FH_Breast_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)

model4 <- coxph(Surv(PAC_difftime_new,PAC_total.y)~as.factor(FH_Breast_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)

model5 <- coxph(Surv(OEC_difftime_new,ESC_total)~as.factor(FH_Breast_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)

model6 <- coxph(Surv(CEUC_difftime_new,CEUC_total)~as.factor(FH_Breast_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)

model7 <- coxph(Surv(Hodgkin_difftime_new,Hodgkin_total.y)~as.factor(FH_Breast_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)

model8 <- coxph(Surv(MM_difftime_new,MM_total)~as.factor(FH_Breast_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)

model9 <- coxph(Surv(LOCPC_difftime_new,LOCPC_total)~as.factor(FH_Breast_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)

model10 <- coxph(Surv(SKM_difftime_new,SKM_total)~as.factor(FH_Breast_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)

model11 <- coxph(Surv(COUC_difftime_new,COUC_total)~as.factor(FH_Breast_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)

model12 <- coxph(Surv(BLC_difftime_new,BLC_total.y)~as.factor(FH_Breast_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)

model13 <- coxph(Surv(THC_difftime_new,THC_total)~as.factor(FH_Breast_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)

model14 <- coxph(Surv(BCNC_difftime_new,BCNC_total)~as.factor(FH_Breast_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)

model15 <- coxph(Surv(KIC_difftime_new,KIC_total.y)~as.factor(FH_Breast_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)

model16 <- coxph(Surv(LYLE_difftime_new,LYLE_total)~as.factor(FH_Breast_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)

model17 <- coxph(Surv(OVC_difftime_new,OVC_total.y)~as.factor(FH_Breast_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)

model18 <- coxph(Surv(survival_time,cancer_total_20)~as.factor(FH_Breast_2)+Age+Height_imp+Townsend_deprivation_index_imp_cat+smoking_health+Alcohol_health+BMI_health+physical_health+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)


summary_HR_CI <- data.frame(rbind(
  HR_CI(model1),
  HR_CI(model2),
  HR_CI(model3),
  HR_CI(model4),
  HR_CI(model5),
  HR_CI(model6),
  HR_CI(model7),
  HR_CI(model8),
  HR_CI(model9),
  HR_CI(model10),
  HR_CI(model11),
  HR_CI(model12),
  HR_CI(model13),
  HR_CI(model14),
  HR_CI(model15),
  HR_CI(model16),
  HR_CI(model17),
  HR_CI(model18)
))

summary_HR_CI

write.csv(summary_HR_CI,file = "female_breast_FH_with_other_cancer_risk_included_overall_cancer_risk_summary.csv",row.names = F)

#FDR校正
# 创建一个包含原始P值的向量
p_values <- summary_HR_CI$P.value

#[1] "0.67081"    "0.50278"    "0.67064"    "0.00606"    "0.9158"     "0.44135"    "0.13463"    "0.74594"    "0.68939"   
#[10] "0.57655"    "0.18748"    "0.90307"    "0.72563"    "0.9517"     "0.68024"    "0.83137"    "0.01172"    "< 2.22e-16"

# 使用p.adjust()函数进行FDR校正
fdr_adjusted <- p.adjust(p_values, method = "fdr")

# 输出校正后的p值  即q值
print(fdr_adjusted)

# [1] 0.9517000 0.9517000 0.9517000 0.0996200 0.9517000 0.9517000 0.7629033 0.9517000 0.9517000 0.9517000 0.7967900 0.9517000
#[13] 0.9517000 0.9517000 0.9517000 0.9517000 0.0996200        NA





##森林图##

library("forestplot")

setwd("Figure 2/中间文件")

#male 肠癌

###森林图
   HR <- c(NA, 1.02, 1.06, 1.13, 1.09, 0.34, 0.73, 1.54, 1.06, 0.75, 1.04, 0.99, 1.08, 0.85, 1.09, 0.89, 1.07)
lower <- c(NA, 0.88, 1.00, 0.85, 0.86, 0.11, 0.56, 0.80, 0.81, 0.56, 0.87, 0.81, 0.59, 0.62, 0.89, 0.65, 1.03)
upper <- c(NA, 1.17, 1.13, 1.49, 1.38, 1.09, 0.96, 2.95, 1.40, 0.99, 1.24, 1.23, 1.98, 1.16, 1.34, 1.23, 1.12)

foredata <- data.frame(HR,lower,upper)
foredata
class(foredata)
write.csv(foredata,file = "male_colorectal_FH_with_other_cancer_risk_included_overall_cancer_risk_summary.csv",row.names = F)

pdata <- read.csv('male_colorectal_FH_with_other_cancer_risk_included_overall_cancer_risk_summary.csv',header=T)
pdata
colnames(pdata)

labeltext <- as.matrix(pdata[,c(1:3)])#将数据框的XX列转换成矩阵
attach(pdata)#绑定数据框
names(pdata)

##森林图
pdf("Forestplot_of_male_colorectal_FH_with_other_cancer_risk_included_overall_cancer_risk.pdf", width = 8, height = 6)

#p2
forestplot(
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
  boxsize = 0.5,
  lty.ci = 1,
  lwd.ci = 2, #修改置信区间线粗细
  ci.vertices = FALSE,
  txt_gp = fpTxtGp(label=gpar(fontfamily="sans"),
                   ticks = gpar(cex = 1),
                   xlab = gpar(cex = 1.5), #修改X坐标的字体大小
                   cex = 1),
  lineheight = "auto",
  line.margin = .1,
  clip = c(0.95, 1.4),
  xticks = c(0.10,0.5,1.0,1.5,2.0,2.5,3.0),
)
dev.off()


topptx(p2,"Forestplot_of_male_colorectal_FH_with_other_cancer_risk_included_overall_cancer_risk.pptx",width = 4,height = 6)





#male 肺癌

###森林图
   HR <- c(1.09, NA, 0.93, 1.04, 1.18, 0.69, 1.44, 1.20, 1.13, 1.18, 1.06, 1.10, 0.90, 1.06, 0.81, 1.05, 1.06)
lower <- c(0.98, NA, 0.88, 0.79, 0.95, 0.30, 1.18, 0.61, 0.87, 0.94, 0.89, 0.91, 0.48, 0.80, 0.65, 0.78, 1.02)
upper <- c(1.21, NA, 0.99, 1.36, 1.47, 1.59, 1.77, 2.37, 1.47, 1.49, 1.25, 1.34, 1.68, 1.40, 1.01, 1.40, 1.10)

foredata <- data.frame(HR,lower,upper)
foredata
class(foredata)
write.csv(foredata,file = "male_lung_FH_with_other_cancer_risk_included_overall_cancer_risk_summary.csv",row.names = F)

pdata <- read.csv('male_lung_FH_with_other_cancer_risk_included_overall_cancer_risk_summary.csv',header=T)
pdata
colnames(pdata)

labeltext <- as.matrix(pdata[,c(1:3)])#将数据框的XX列转换成矩阵
attach(pdata)#绑定数据框
names(pdata)

##森林图
pdf("Forestplot_of_male_lung_FH_with_other_cancer_risk_included_overall_cancer_risk.pdf", width = 8, height = 6)

#p2
forestplot(
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
  boxsize = 0.5,
  lty.ci = 1,
  lwd.ci = 2, #修改置信区间线粗细
  ci.vertices = FALSE,
  txt_gp = fpTxtGp(label=gpar(fontfamily="sans"),
                   ticks = gpar(cex = 1),
                   xlab = gpar(cex = 1.5), #修改X坐标的字体大小
                   cex = 1),
  lineheight = "auto",
  line.margin = .1,
  clip = c(0.95, 1.4),
  xticks = c(0.10,0.5,1.0,1.5,2.0,2.5),
)
dev.off()


topptx(p2,"Forestplot_of_male_lung_FH_with_other_cancer_risk_included_overall_cancer_risk.pptx",width = 4,height = 6)





#male 前列腺癌

###森林图
   HR <- c(1.02, 1.09, NA, 0.73, 1.10, 0.63, 1.01, 0.38, 0.97, 0.94, 1.00, 0.81, 0.91, 0.57, 1.09, 1.25, 1.06)
lower <- c(0.89, 0.92, NA, 0.48, 0.83, 0.23, 0.75, 0.09, 0.68, 0.68, 0.80, 0.61, 0.42, 0.36, 0.86, 0.89, 1.02)
upper <- c(1.16, 1.29, NA, 1.09, 1.47, 1.74, 1.36, 1.57, 1.36, 1.29, 1.24, 1.08, 1.97, 0.90, 1.40, 1.74, 1.10)

foredata <- data.frame(HR,lower,upper)
foredata
class(foredata)
write.csv(foredata,file = "male_prostate_FH_with_other_cancer_risk_included_overall_cancer_risk_summary.csv",row.names = F)

pdata <- read.csv('male_prostate_FH_with_other_cancer_risk_included_overall_cancer_risk_summary.csv',header=T)
pdata
colnames(pdata)

labeltext <- as.matrix(pdata[,c(1:3)])#将数据框的XX列转换成矩阵
attach(pdata)#绑定数据框
names(pdata)

##森林图
pdf("Forestplot_of_male_prostate_FH_with_other_cancer_risk_included_overall_cancer_risk.pdf", width = 8, height = 6)

#p2
forestplot(
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
  boxsize = 0.5,
  lty.ci = 1,
  lwd.ci = 2, #修改置信区间线粗细
  ci.vertices = FALSE,
  txt_gp = fpTxtGp(label=gpar(fontfamily="sans"),
                   ticks = gpar(cex = 1),
                   xlab = gpar(cex = 1.5), #修改X坐标的字体大小
                   cex = 1),
  lineheight = "auto",
  line.margin = .1,
  clip = c(0.95, 1.4),
  xticks = c(0.05,0.5,1.0,1.5,2.0),
)
dev.off()

topptx(p2,"Forestplot_of_male_prostate_FH_with_other_cancer_risk_included_overall_cancer_risk.pptx",width = 4,height = 6)




#female 肠癌

###森林图
   HR <- c(NA, 1.08, 1.11, 1.10, 1.10, 1.04, 1.08, 0.88, 1.13, 1.03, 1.07, 1.09, 0.90, 1.24, 0.87, 1.33, 0.96, 1.19, 1.12)
lower <- c(NA, 0.94, 1.04, 0.72, 0.85, 0.70, 0.59, 0.38, 0.81, 0.73, 0.90, 0.93, 0.60, 0.87, 0.60, 1.03, 0.65, 0.98, 1.07)
upper <- c(NA, 1.24, 1.19, 1.67, 1.43, 1.54, 1.97, 2.05, 1.56, 1.46, 1.28, 1.29, 1.36, 1.78, 1.27, 1.72, 1.41, 1.45, 1.17)

foredata <- data.frame(HR,lower,upper)
foredata
class(foredata)
write.csv(foredata,file = "female_colorectal_FH_with_other_cancer_risk_included_overall_cancer_risk_summary.csv",row.names = F)

pdata <- read.csv('female_colorectal_FH_with_other_cancer_risk_included_overall_cancer_risk_summary.csv',header=T)
pdata
colnames(pdata)

labeltext <- as.matrix(pdata[,c(1:3)])#将数据框的XX列转换成矩阵
attach(pdata)#绑定数据框
names(pdata)

##森林图
pdf("Forestplot_of_female_colorectal_FH_with_other_cancer_risk_included_overall_cancer_risk.pdf", width = 8, height = 6)

#p2
forestplot(
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
  boxsize = 0.5,
  lty.ci = 1,
  lwd.ci = 2, #修改置信区间线粗细
  ci.vertices = FALSE,
  txt_gp = fpTxtGp(label=gpar(fontfamily="sans"),
                   ticks = gpar(cex = 1),
                   xlab = gpar(cex = 1.5), #修改X坐标的字体大小
                   cex = 1),
  lineheight = "auto",
  line.margin = .1,
  clip = c(0.95, 1.4),
  xticks = c(0.1,0.5,1.0,1.5,2.0,2.5),
)
dev.off()


topptx(p2,"Forestplot_of_female_colorectal_FH_with_other_cancer_risk_included_overall_cancer_risk.pptx",width = 4,height = 6)




#female 肺癌

###森林图
   HR <- c(0.99, NA, 1.06, 1.06, 0.92, 1.20, 1.13, 0.90, 0.96, 0.99, 0.93, 0.95, 1.16, 0.94, 0.96, 0.98, 0.92, 1.06, 1.12)
lower <- c(0.88, NA, 1.00, 0.71, 0.71, 0.84, 0.65, 0.41, 0.69, 0.72, 0.78, 0.81, 0.82, 0.65, 0.69, 0.75, 0.64, 0.88, 1.07)
upper <- c(1.12, NA, 1.13, 1.58, 1.19, 1.70, 1.96, 1.98, 1.33, 1.38, 1.10, 1.11, 1.63, 1.37, 1.35, 1.28, 1.33, 1.28, 1.17)

foredata <- data.frame(HR,lower,upper)
foredata
class(foredata)
write.csv(foredata,file = "female_lung_FH_with_other_cancer_risk_included_overall_cancer_risk_summary.csv",row.names = F)

pdata <- read.csv('female_lung_FH_with_other_cancer_risk_included_overall_cancer_risk_summary.csv',header=T)
pdata
colnames(pdata)

labeltext <- as.matrix(pdata[,c(1:3)])#将数据框的XX列转换成矩阵
attach(pdata)#绑定数据框
names(pdata)

##森林图
pdf("Forestplot_of_female_lung_FH_with_other_cancer_risk_included_overall_cancer_risk.pdf", width = 8, height = 6)

#p2
forestplot(
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
  boxsize = 0.5,
  lty.ci = 1,
  lwd.ci = 2, #修改置信区间线粗细
  ci.vertices = FALSE,
  txt_gp = fpTxtGp(label=gpar(fontfamily="sans"),
                   ticks = gpar(cex = 1),
                   xlab = gpar(cex = 1.5), #修改X坐标的字体大小
                   cex = 1),
  lineheight = "auto",
  line.margin = .1,
  clip = c(0.95, 1.4),
  xticks = c(0.1,0.5,1.0,1.5,2.0),
)
dev.off()

topptx(p2,"Forestplot_of_female_lung_FH_with_other_cancer_risk_included_overall_cancer_risk.pptx",width = 4,height = 6)




#female 乳腺癌

###森林图
   HR <- c(1.03, 1.05, NA, 0.90, 1.40, 0.98, 0.76, 0.41, 0.94, 0.93, 0.95, 1.12, 0.97, 0.93, 0.99, 1.06, 1.04, 1.28, 1.26)
lower <- c(0.90, 0.91, NA, 0.57, 1.10, 0.64, 0.39, 0.13, 0.66, 0.64, 0.79, 0.95, 0.65, 0.63, 0.69, 0.80, 0.71, 1.06, 1.21)
upper <- c(1.17, 1.22, NA, 1.44, 1.79, 1.48, 1.51, 1.32, 1.35, 1.34, 1.14, 1.32, 1.47, 1.38, 1.41, 1.41, 1.53, 1.55, 1.31)

foredata <- data.frame(HR,lower,upper)
foredata
class(foredata)
write.csv(foredata,file = "female_breast_FH_with_other_cancer_risk_included_overall_cancer_risk_summary.csv",row.names = F)

pdata <- read.csv('female_breast_FH_with_other_cancer_risk_included_overall_cancer_risk_summary.csv',header=T)
pdata
colnames(pdata)

labeltext <- as.matrix(pdata[,c(1:3)])#将数据框的XX列转换成矩阵
attach(pdata)#绑定数据框
names(pdata)

##森林图
pdf("Forestplot_of_female_breast_FH_with_other_cancer_risk_included_overall_cancer_risk.pdf", width = 8, height = 6)

#p2
forestplot(
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
  boxsize = 0.5,
  lty.ci = 1,
  lwd.ci = 2, #修改置信区间线粗细
  ci.vertices = FALSE,
  txt_gp = fpTxtGp(label=gpar(fontfamily="sans"),
                   ticks = gpar(cex = 1),
                   xlab = gpar(cex = 1.5), #修改X坐标的字体大小
                   cex = 1),
  lineheight = "auto",
  line.margin = .1,
  clip = c(0.95, 1.4),
  xticks = c(0.1,0.5,1.0,1.5,2.0),
)
dev.off()


topptx(p2,"Forestplot_of_female_breast_FH_with_other_cancer_risk_included_overall_cancer_risk.pptx",width = 4,height = 6)



