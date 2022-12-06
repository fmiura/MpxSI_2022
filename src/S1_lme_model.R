############################################################
#Aim: to check the significant effects of reporting reliability on the onset interval (mixed error-component model)
#Final edit: 2 December 2022
#Editor: Fumi Miura
############################################################
###Procedure
#0. Package 
#1. NL data 
#2. mixed effect model
#3. Summary of results
############################################################

###0. Package -----
library(tidyverse)
library(nlme)
library(modelsummary)

###1. NL data -----
Anonym_All_data_score_indexed <- read_csv("Anonym_All_data_score_indexed.csv")
Anonym_All_data_score_indexed$GGD_index <- as.character(Anonym_All_data_score_indexed$GGD_index)

###2. mixed effect model -----
#fixed effects only
nl.lm <- lm(formula = SI ~ Conf_pair + Conf_sym,
            data = Anonym_All_data_score_indexed) 
#summary(nl.lm); round(coefficients(nl.lm),4)

#fixed effects + random effect
nl.lme <- lme(SI ~ Conf_pair + Conf_sym,
              random = ~ 1|GGD_index, 
              data = Anonym_All_data_score_indexed) 
#summary(nl.lme)$tTable; round(summary(nl.lme)$tTable,4)

###3. Summary of results -----
regs <- list()
regs[['Model A']] <- nl.lme
regs[['Model B']] <- nl.lm
msummary(regs, statistic = 'conf.int', conf_level = .95, fmt = '%.1f', 'tableS2.rtf')
#msummary(regs, fmt = '%.1f', stars = TRUE, 'tableS2_v2.rtf') #with stars; + p < 0.1, * p < 0.05, ** p < 0.01, *** p < 0.001