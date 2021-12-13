#### Load packages ####
library(ggpubr) 
library(lubridate)
library(gtsummary)
library(MASS)
library(sandwich)
library(lmtest)
library(broom)
library(mice)
library(tidyverse)
library(lavaan)
library(snow)

#### Load Data and format variables ####
answers <- read_csv("c:/users/atala/sync/research/projects/to publish/answers study/phase I survey/ANSWERSI_Data_Scored.csv") %>% 
  filter(Age<26) %>%
  mutate(
   "ID" = factor(ID),
   "Sex" = factor(Sex, levels = c("Male","Female")),
   "Race" = factor(Race,levels = c("White","Black","Asian","Native","Multiracial")),
   "Ethnicity" = factor(Ethnicity, levels = c("Non-Hispanic","Hispanic")),
   "WD_SQ" = as.numeric(WD_SQ),
   "WE_SQ" = as.numeric(WE_SQ),
   "WD_TST" = round(WD_TST/60,2),
   "WE_TST" = round(WE_TST/60,2),
   "TST" = (5*WD_TST + 2*WE_TST)/7,
   "WD_TWT" = WD_TWT/60, 
   "WE_TWT" = WE_TWT/60,
   "TWT" = (5*WD_TWT+2*WE_TWT)/7,
   "SE" = (5*WD_SE+2*WE_SE)/7,
   "MSFsc2" = MSFsc2*24,
   "MSFsc1" = MSFsc1*24,
   "SocialJetlag" = SocialJetlag/60,
   "absSocialJetlag" = abs(SocialJetlag),
   "isi_cat" = factor(case_when(isi_total<8 ~ "No Insomnia",
                                isi_total>7 & isi_total < 15 ~ "Mild Insomnia",
                                isi_total>14 & isi_total < 22 ~ "Moderate Insomnia",
                                isi_total>21 ~ "Severe Insomnia"),
                      levels = c("No Insomnia","Mild Insomnia",'Moderate Insomnia',"Severe Insomnia")),
   "isi_cat2" = factor(case_when(isi_total<15 ~ "No Insomnia",
                                 isi_total>=15 ~ "Insomnia"),
                       levels = c("No Insomnia","Insomnia")),
   "isi_cat2n" = as.numeric(isi_cat2)-1,
   "ddnsi_cat" = factor(case_when(ddnsi_total<=10 ~ "No Nightmare Disorder",
                                  ddnsi_total>10 ~ "Nightmare Disorder"),
                        levels = c("No Nightmare Disorder","Nightmare Disorder")),
   "ddnsi_catn" = as.numeric(ddnsi_cat)-1,
   "LifetimeSI" = case_when(SI_lf_severity<2&SI_3mo_severity<2~0,
                             SI_lf_severity>1|SI_3mo_severity>1~1),
   "RecentSI" = case_when(LifetimeSI==0~NA_real_,
                           LifetimeSI==1&SI_3mo_severity<2~0,
                           LifetimeSI==1&SI_3mo_severity>1~1),
   "LifetimeSA" = case_when(LifetimeSI==0 & SA_lf=="No" & SA_3mo=="No" ~ NA_real_,
                            LifetimeSI==1 & SA_lf=="No" & SA_3mo=="No" ~ 0,
                            LifetimeSI==1 & SA_lf=="Yes" | SA_3mo=="Yes" ~ 1),
   "RecentSA" = case_when(LifetimeSA==0 ~ NA_real_,
                          LifetimeSA==1 & SA_3mo=="No" ~ 0,
                          LifetimeSA==1 & SA_3mo=="Yes" ~ 1),
   "SIcat" = factor(case_when(LifetimeSI==0 ~ "None",
                              LifetimeSI==1 & RecentSI==0 ~ "Lifetime", 
                              LifetimeSI==1 & RecentSI==1 ~ "Recent"),
                    levels = c("None",'Lifetime',"Recent")),
   "SAcat" = factor(case_when(LifetimeSA==0 ~ "None",
                              LifetimeSA==1 & RecentSA==0 ~ "Lifetime", 
                              LifetimeSA==1 & RecentSA==1 ~ "Recent"),
                    levels = c("None",'Lifetime',"Recent")),
   "CESD_sad" = if_else(cesd_2==3|abs(cesd_4-3)==3|cesd_6==3,1,0),
   "CESD_anh" = if_else(abs(cesd_8-3)==3|cesd_10==3,1,0),
   "CESD_app" = if_else(cesd_1==3|cesd_18==3,1,0),
   "CESD_slp" = if_else(cesd_5==3|cesd_11==3|cesd_19==3,1,0),
   "CESD_con" = if_else(cesd_3==3|cesd_20==3,1,0),
   "CESD_glt" = if_else(cesd_9==3|cesd_17==3,1,0),
   "CESD_trd" = if_else(cesd_7==3|abs(cesd_16-3)==3,1,0),
   "CESD_mvt" = if_else(cesd_13==3|abs(cesd_12-3)==3,1,0),
   "CESD_scd" = if_else(cesd_14==3|cesd_15==3,1,0),
   "CESDcat" = case_when(cesd_total<16 ~ "No Clinical Depression",
                         CESD_sad==1&CESD_anh==1&(CESD_app+CESD_slp+CESD_con+CESD_glt+CESD_trd+CESD_mvt+CESD_scd)>=4 ~ "Clinical Depression",
                         CESD_sad==1&CESD_anh==1&(CESD_app+CESD_slp+CESD_con+CESD_glt+CESD_trd+CESD_mvt+CESD_scd)==3 ~ "Clinical Depression",
                         CESD_sad==1&CESD_anh==1&(CESD_app+CESD_slp+CESD_con+CESD_glt+CESD_trd+CESD_mvt+CESD_scd)==2 ~ "Clinical Depression",
                         CESD_sad==1&CESD_anh==1&(CESD_app+CESD_slp+CESD_con+CESD_glt+CESD_trd+CESD_mvt+CESD_scd)<=1 ~ "Subclinical Depression",
                         CESD_sad==0|CESD_anh==0&cesd_total>=16 ~ "Subclinical Depression") %>%
     factor(levels = c("No Clinical Depression","Subclinical Depression","Clinical Depression")),
   "GADcat" = if_else(gad_total>9,1,0) %>% factor(0:1),
   "subs_caff_1" = factor(subs_caff_1,c("Never","Once a month or less","Once a week or less",
                                                     "A few times a week","Daily","Multiple times daily")),
   "subs_alc_1" = factor(subs_alc_1,
                         c("Never","Once a month or less","Once a week or less","A few times a week","Daily","Multiple times daily"),
                         c("Never","Once a week or less","Once a week or less",
                           "More than once a week","More than once a week","More than once a week")),
   "subs_alc_3" = factor(subs_alc_3,c("Never","Rarely","Often")),
   "subs_can_1" = factor(subs_can_1,
                         c("Never","Once a month or less","Once a week or less","A few times a week","Daily","Multiple times daily"),
                         c("Never","Once a week or less","Once a week or less",
                           "More than once a week","More than once a week","More than once a week")),
   "subs_can_3" = factor(subs_can_3,c("Never","Rarely","Often")),
   "BedPartner" = factor(BedPartner,c("No partner","Partner in other room",
                                       "Partner in other bed","Partner in same bed"))
 )
## Impute missing MSFsc2
answers.id <- answers$ID
tempmids <- mice(dplyr::select(answers[,2:193],-MSFsc1,-psqi_4,-starts_with("cssrs_")),maxit=0)
pmatrix <- tempmids$predictorMatrix
answers.mids <- mice(dplyr::select(answers[,2:193],-MSFsc1,-psqi_4,-starts_with("cssrs_")),m=10,seed=123,predictorMatrix = pmatrix)
answers.imp <- complete(answers.mids)
answers <- answers %>% mutate(
  "MSFsc2" = answers.imp$MSFsc2,
  "MSFsccat" = factor(case_when(MSFsc2<=quantile(MSFsc2,.1,na.rm=T)~"Early",
                                MSFsc2>=quantile(MSFsc2,.9,na.rm=T)~"Late",
                                MSFsc2>quantile(MSFsc2,.1,na.rm=T)&MSFsc2<quantile(MSFsc2,.9,na.rm=T)~"Typical"),
                      levels = c("Typical","Early","Late"))
)

#### Begin Analyses ####
### Outcome 1: Suicidal Ideation
## Table 1: Sample characteristics by suicidal ideation 
tbl1 <- answers %>% 
  dplyr::select(Age,Sex,Race,Ethnicity,psqi_total,CESDcat,GADcat,
                subs_alc_1,subs_can_1,SIcat) %>%
  tbl_summary(statistic= list(all_continuous() ~ "{mean} ({sd})",
                              all_categorical() ~ "{n} ({p}%)"),
              digits = list(all_continuous() ~ c(1,2)),
              missing = "no",
              by=c("SIcat")) %>% 
  add_p(test = list(all_continuous() ~ "aov",all_categorical() ~ "chisq.test"),
        pvalue_fun = function(x) style_pvalue(x, digits=3))
## Part 1: Lifetime SI differences in sleep parameters
# Create subset dataset with desired variables
lifetimesi.data <- answers %>% 
  transmute(ID,LifetimeSI,
            "wTIB" = (5*WD_TIB+2*WE_TIB)/7,
            "wTST" = 60*(5*WD_TST+2*WE_TST)/7,
            "wTWT" = 60*(5*WD_TWT+2*WE_TWT)/7,
            "wSE" = (5*WD_SE+2*WE_SE)/7,
            "wSQ" = (5*WD_SQ+2*WE_SQ)/7,
            MSFsc2,MSFsccat,SocialJetlag,absSocialJetlag,brisc_total,
            isi_cat2n,ddnsi_catn,
            Age,Sex,Race,Ethnicity,cesd_total,gad_total,CESDcat,GADcat,isi_total,ddnsi_total,
            "Alcohol" = subs_alc_1, "Cannabis" = subs_can_1)
# Run iterative linear modeling
outcomes <- c("wTIB","wTST","wTWT","wSE","wSQ","SocialJetlag","absSocialJetlag","MSFsc2","brisc_total")
covariates <- c("LifetimeSI","LifetimeSI+Sex+Race+Ethnicity+CESDcat+GADcat+Alcohol+Cannabis")
modelnames <- c("Unadjusted","Adjusted")
simodels.lm <- tibble()
for(idx in 1:length(outcomes)){
  for(jdx in 1:length(covariates)){
    simodels.lm <- bind_rows(
      simodels.lm,
      tidy(lm(as.formula(paste(outcomes[idx],"~",covariates[jdx])),lifetimesi.data),conf.int=T) %>%
        filter(term == "LifetimeSI") %>%
        transmute("Outcome" = outcomes[idx],
                  "Model" = modelnames[jdx], 
                  estimate = estimate, 
                  "Lower" = (conf.low), 
                  "Upper" = (conf.high), 
                  "P" = round(p.value,4))
    )
  }
}
# Run iterative robust Poisson modeling
outcomes <- c("isi_cat2n","ddnsi_catn")
simodels.glm <- tibble()
for(idx in 1:length(outcomes)){
  for(jdx in 1:length(covariates)){
    simodel.temp <- glm(as.formula(paste(outcomes[idx],"~",covariates[jdx])),lifetimesi.data,family="poisson")
    simodel.sandwich <- sandwich(simodel.temp)
    simodel.coef <- coeftest(simodel.temp,vcov. = simodel.sandwich) %>% tidy(conf.int=T) %>% 
      filter(term=="LifetimeSI") %>% 
      transmute("Outcome" = outcomes[idx],
                "Model" = modelnames[jdx], 
                estimate = exp(estimate), 
                "Lower" = exp(conf.low), 
                "Upper" = exp(conf.high), 
                "P" = round(p.value,4))
    simodels.glm <- bind_rows(simodels.glm,simodel.coef)
  }  
}
# Combine linear and Poisson models
simodels <- bind_rows(simodels.lm,simodels.glm) %>%
  mutate("95% CI" = paste("[",round(Lower,2),", ",round(Upper,2),"]",sep=""))
# Table 2
# write_csv(simodels, "c:/users/atala/sync/research/projects/to publish/answers study/phase I survey/suicidality/table2.csv")

## Part 2: Variable selection for Recent SI
recentsi.data <- answers %>% filter(is.na(RecentSI)==F) %>%
  dplyr::select("RecentSI","WD_TIB","WD_TST","WD_TWT","WD_SE","WD_SQ",
                "WE_TIB","WE_TST","WE_TWT","WE_SE","WE_SQ",
                "MSFsc2","SocialJetlag","absSocialJetlag",
                "brisc_total","isi_total","ddnsi_total")
recentsi.step <- stepAIC(glm(RecentSI ~ ., recentsi.data,family="binomial"),
                         trace=F,direction = "both")
recentsi.m <- glm(RecentSI ~ isi_total,family="poisson",recentsi.data)
recentsi.m.final <- coeftest(recentsi.m,vcov. = sandwich(recentsi.m)) %>% tidy(conf.int=T) %>%
  filter(term!="(Intercept)") %>% 
  transmute("Predictor" = term, 
            "PR" = round(exp(estimate),2),
            "Lower" = round(exp(conf.low),2),
            "Upper" = round(exp(conf.high),2),
            "95% CI" = paste("[",Lower,", ",Upper,"]",sep=""),
            "pvalue" = round(p.value,4))
# Robust Poisson Modeling of Recent SI
covariates <-  c("",
                 "+Sex+Race+Ethnicity+CESDcat+GADcat+subs_alc_1+subs_can_1")
modelnames <- c("Unadjusted","Adjusted")
si.modeltable <- tibble()
for(idx in 1:length(modelnames)){
  modelform <- as.formula(paste("RecentSI~isi_total",covariates[idx]))
  si.m <- coeftest(glm(modelform,filter(answers,is.na(RecentSI)==F),family="poisson"),
                    sandwich(glm(modelform,filter(answers,is.na(RecentSI)==F),family="poisson"))) %>%
    tidy(conf.int=T) %>% filter(term!="(Intercept)") %>%
    transmute("Model" = modelnames[idx],
              "Predictors" = term,
              "PR" = round(exp(estimate),2),
              "Lower" = round(exp(conf.low),2),
              "Upper" = round(exp(conf.high),2),
              "95% CI" = paste("[",Lower,", ",Upper,"]",sep=""),
              "pvalue" = round(p.value,4))
  si.modeltable <- bind_rows(si.modeltable,si.m)
}

### Outcome 2: Suicide Attempt
tbl3 <- answers %>% 
  dplyr::select(Age,Sex,Race,Ethnicity,psqi_total,CESDcat,GADcat,
                subs_alc_1,subs_can_1,SAcat) %>%
  tbl_summary(statistic= list(all_continuous() ~ "{mean} ({sd})",
                              all_categorical() ~ "{n} ({p}%)"),
              digits = list(all_continuous() ~ c(1,2)),
              missing = "no",
              by=c("SAcat")) %>% 
  add_p(test = list(all_continuous() ~ "aov",all_categorical() ~ "chisq.test"),
        pvalue_fun = function(x) style_pvalue(x, digits=3))

## Part 1: Lifetime SA didfferences on sleep parameters
# Create data subset with desired variables
lifetimesa.data <- answers %>% filter(is.na(LifetimeSA)==F) %>%
  transmute(ID,LifetimeSA,
            "wTIB" = (5*WD_TIB+2*WE_TIB)/7,
            "wTST" = 60*(5*WD_TST+2*WE_TST)/7,
            "wTWT" = 60*(5*WD_TWT+2*WE_TWT)/7,
            "wSE" = (5*WD_SE+2*WE_SE)/7,
            "wSQ" = (5*WD_SQ+2*WE_SQ)/7,
            MSFsc2,MSFsccat,SocialJetlag,absSocialJetlag,brisc_total,
            isi_cat2n,ddnsi_catn,
            Age,Sex,Race,Ethnicity,cesd_total,gad_total,CESDcat,GADcat,isi_total,ddnsi_total,
            "Alcohol" = subs_alc_1, "Cannabis" = subs_can_1)
# Run iterative linear models
outcomes <- c("wTIB","wTST","wTWT","wSE","wSQ","SocialJetlag","absSocialJetlag","MSFsc2","brisc_total")
covariates <- c("LifetimeSA","LifetimeSA+Sex+Race+Ethnicity+CESDcat+GADcat+Alcohol+Cannabis")
modelnames <- c("Unadjusted","Adjusted")
samodels.lm <- tibble()
for(idx in 1:length(outcomes)){
   for(jdx in 1:length(covariates)){
      samodels.lm <- bind_rows(
         samodels.lm,
         tidy(lm(as.formula(paste(outcomes[idx],"~",covariates[jdx])),lifetimesa.data),conf.int=T) %>%
            filter(term == "LifetimeSA") %>%
            transmute("Outcome" = outcomes[idx],
                      "Model" = modelnames[jdx], 
                      estimate = estimate, 
                      "Lower" = (conf.low), 
                      "Upper" = (conf.high), 
                      "P" = round(p.value,4))
      )
   }
}
# Run iterative robust Poisson models
outcomes <- c("isi_cat2n","ddnsi_catn")
samodels.glm <- tibble()
for(idx in 1:length(outcomes)){
   for(jdx in 1:length(covariates)){
      samodel.temp <- glm(as.formula(paste(outcomes[idx],"~",covariates[jdx])),lifetimesa.data,family="poisson")
      samodel.sandwich <- sandwich(samodel.temp)
      samodel.coef <- coeftest(samodel.temp,vcov. = samodel.sandwich) %>% tidy(conf.int=T) %>% 
         filter(term=="LifetimeSA") %>% 
         transmute("Outcome" = outcomes[idx],
                   "Model" = modelnames[jdx], 
                   estimate = exp(estimate), 
                   "Lower" = exp(conf.low), 
                   "Upper" = exp(conf.high), 
                   "P" = round(p.value,4))
      samodels.glm <- bind_rows(samodels.glm,samodel.coef)
   }
}
# Combine linear and Poisson models
samodels <- bind_rows(samodels.lm,samodels.glm) %>%
  mutate("95% CI" = paste("[",round(Lower,2),", ",round(Upper,2),"]",sep=""))
# Table 4
# write_csv(samodels, "c:/users/atala/sync/research/projects/to publish/answers study/phase I survey/suicidality/table4.csv")

## Part 2: Variable selection for Recent SI
recentsa.data <- answers %>% filter(is.na(RecentSA)==F) %>%
   dplyr::select("RecentSA","WD_TIB","WD_TST","WD_TWT","WD_SE","WD_SQ",
                 "WE_TIB","WE_TST","WE_TWT","WE_SE","WE_SQ",
                 "MSFsc2","SocialJetlag","absSocialJetlag",
                 "brisc_total","isi_total","ddnsi_total")
recentsa.step <- stepAIC(glm(RecentSA ~ ., recentsa.data,family="binomial"),
                         trace=F,direction = "both")
recentsa.m <- glm(RecentSA ~ WD_SQ+SocialJetlag,family="poisson",recentsa.data)
recentsa.m.final <- coeftest(recentsa.m,vcov. = sandwich(recentsa.m)) %>% tidy(conf.int=T) %>%
   filter(term!="(Intercept)") %>% 
   transmute("Predictor" = term, 
             "PR" = round(exp(estimate),2),
             "Lower" = round(exp(conf.low),2),
             "Upper" = round(exp(conf.high),2),
             "95% CI" = paste("[",Lower,", ",Upper,"]",sep=""),
             "pvalue" = round(p.value,4))

