##################################################################################################
##################################################################################################
###
### Code to run Example 1 of data-driven simulations, for assessing the impact of omitting an  
### important prognostic factor (see the manuscript for details)
###
### Code by: Marie-Eve Beauchamp                                            
### Last update: March 19, 2024                                              
###
##################################################################################################
##################################################################################################

#rm(list=ls())

## Set the path to the directory for Example 1 (need to be changed by users)
setwd("C:/.../Example1/")

## Install the R package PermAlgo if needed
#install.packages("PermAlgo")
library(PermAlgo)

## Source the script containing the functions needed for the simulations and to produce results
source('./Code_ex1_functions.R')


##################################################################################################
# DATA PREPARATION AND PRE-ANALYSIS BEFORE SIMULATIONS
##################################################################################################

#-----------------
# Data preparation
#-----------------

## Colon dataset from the 'survival' package 

  # Documentation available at:
  # https://cran.r-project.org/web/packages/survival/survival.pdf

library(survival)
summary(colon)

## Selection of relevant variables 

colon2 <- colon[, c('id', 'rx', 'sex','age', 'obstruct', 'perfor', 'adhere', 'status', 'differ', 
                    'extent', 'surg', 'node4', 'time', 'etype')]

## Selection of observations for the event type "death"

table(colon2$etype)
colon2 <- colon2[colon2$etype == 2, ]
dim(colon2)
length(unique(colon2$id))

## Select complete cases 

apply(is.na(colon2), 2, sum) # 23 patients with NA for variable 'differ' are excluded
colon3 <- colon2[complete.cases(colon2[, c('differ')]), ]

## Creation of dummy variables for categorical variables 

table(colon3$rx) # reference level = 'Obs'
colon3$rx_Lev <- ifelse(colon3$rx == "Lev", 1, 0)
colon3$rx_Lev5FU <- ifelse(colon3$rx == "Lev+5FU", 1, 0)

table(colon3$differ) # reference level = 1
colon3$differ_mod <- ifelse(colon3$differ == 2, 1, 0)
colon3$differ_poor <- ifelse(colon3$differ == 3, 1, 0)

table(colon3$extent) # reference level = 1
colon3$extent_muscle <- ifelse(colon3$extent == 2, 1, 0)
colon3$extent_serosa <- ifelse(colon3$extent == 3, 1, 0)
colon3$extent_cs <- ifelse(colon3$extent == 4, 1, 0)

## Descriptive statistics for the final dataset used in the simulations

# Number of patients
length(unique(colon3$id))

# Number of death events
table(colon3$status)      

# Proportion of patients with colon obstruction (exposure) 
table(colon3$obstruct)     
(sum(colon3$obstruct == 1) / length(unique(colon3$id))) * 100

#-----------------------------------------------------------------------------------
# Logistic regression to estimate the association of each covariate with obstruction  
#-----------------------------------------------------------------------------------

lr.colon <- glm(obstruct ~ rx_Lev + rx_Lev5FU + sex + age + perfor + adhere + differ_mod + 
                  differ_poor + extent_muscle + extent_serosa + extent_cs + surg + node4,
                family = binomial, data = colon3)

# Results for Web Table 1 of the manuscript 
round(exp(cbind(lr.colon$coefficients, confint(lr.colon))), digits = 2)

#--------------------------------------------------------------------------------------
# Unadjusted odds ratios (ORs) between pairs of variables to evaluate their correlation
#--------------------------------------------------------------------------------------

# Create a variable measuring 10 years of increase in 'age'
colon3$age10 <- colon3$age / 10

# Select relevant variables for which ORs are calculated 
colon.selVar <- colon3[, c('obstruct', 'rx_Lev', 'rx_Lev5FU', 'sex', 'age10', 'perfor', 'adhere', 
                           'differ_mod', 'differ_poor', 'extent_muscle', 'extent_serosa', 
                           'extent_cs', 'surg', 'node4')]

# Create a matrix of ORs, with dependent and independent variables in rows and columns respectively
ORmat <- matrix(NA, nrow = ncol(colon.selVar), ncol = ncol(colon.selVar))
rownames(ORmat) <- colnames(colon.selVar)
colnames(ORmat) <- colnames(colon.selVar)

# Calculate ORs for each pair of variables
for (i in 1:ncol(colon.selVar)){
  for (j in 1:ncol(colon.selVar)){
    if (j < i){  # do not calculate OR for a variable with itself
      if (colnames(colon.selVar)[i] == 'age10'){  
        ORmat[i, j] <- NA  # 'age10' is not binary so it cannot be the dependent variable
      } else {
        ORmat[i, j] <- exp(glm(colon.selVar[, i] ~ colon.selVar[, j], family = binomial)$coef[2])
      }
    }
  }
}

# Replace ORs by NA for pairs of dummy variables from the same categorical variable 
ORmat['rx_Lev5FU', 'rx_Lev'] <- NA 
ORmat['differ_poor', 'differ_mod'] <- NA
ORmat['extent_serosa', 'extent_muscle'] <- NA; ORmat['extent_cs', 'extent_muscle'] <- NA;
ORmat['extent_cs', 'extent_serosa'] <- NA

# Results for Web Table 2 of the manuscript
round(ORmat, digits = 2)  

#------------------------------------------------------------------------------------------------
# Cox model for the adjusted associations of obstruction (exposure) and covariates with mortality 
# hazard
# [Step 2 (Initial analyses not corrected for the imperfection) of the proposed approach to
# data-driven simulations]
#------------------------------------------------------------------------------------------------

cox.colon <- coxph(Surv(time, status) ~ obstruct + rx_Lev + rx_Lev5FU + sex + age + perfor +  
                                        adhere + differ_mod + differ_poor + extent_muscle + 
                                        extent_serosa + extent_cs + surg + node4,
                   data = colon3)

# Results for Web Table 1 of the manuscript
round(summary(cox.colon)$conf.int, digits = 2)[, c(1, 3:4)]

# Save the coefficients for adjustment variables (i.e. excluding 'obstruct') for future use in the
# simulations
true.adj.coef <- cox.colon$coef[-1]


##################################################################################################
# SIMULATIONS
# [Steps 4-6 of the proposed approach to data-driven simulations]
##################################################################################################

# NOTE: Steps 4-6 are performed by the function 'simulations.ex1.fct'. For details see the script 
# 'Code_ex1_functions.R', which includes the code of 'simulations.ex1.fct' and the definition of  
# the arguments of the function.

## Run the simulations separately for each simulation scenario

options(warn = 1)  # To request warnings to be printed as they occur

# NOTE: this warning message below occurs for some simulation repetitions. However, it does not 
# affect results for the variable of interest (obstruction exposure) because the warning concerns
# the variables 'extent_muscle', 'extent_serosa', 'extent_cs'.
    #Warning in coxph.fit(X, Y, istrat, offset, init, control, weights = weights,  :
    #  Loglik converged before variable  10,11,12 ; coefficient may be infinite. 

# Scenario 1
set.seed(939001)  # Set seed to ensure reproducibility
simulations.ex1.fct(n.reps = 1000, colon.dat = colon3, HR.obstruct = 1.0, OR.obstruct.stage = 1.2,
                    HR.stage = 4.0, HRs.adj = exp(true.adj.coef))

# Scenario 2
set.seed(939002)  # Set seed to ensure reproducibility
simulations.ex1.fct(n.reps = 1000, colon.dat = colon3, HR.obstruct = 1.3, OR.obstruct.stage = 1.2,
                    HR.stage = 4.0, HRs.adj = exp(true.adj.coef))

# Scenario 3
set.seed(939003)  # Set seed to ensure reproducibility
simulations.ex1.fct(n.reps = 1000, colon.dat = colon3, HR.obstruct = 1.5, OR.obstruct.stage = 1.2,
                    HR.stage = 4.0, HRs.adj = exp(true.adj.coef))

# Scenario 4
set.seed(939004)  # Set seed to ensure reproducibility
simulations.ex1.fct(n.reps = 1000, colon.dat = colon3, HR.obstruct = 2.0, OR.obstruct.stage = 1.2,
                    HR.stage = 4.0, HRs.adj = exp(true.adj.coef))

# Scenario 5
set.seed(939005)  # Set seed to ensure reproducibility
simulations.ex1.fct(n.reps = 1000, colon.dat = colon3, HR.obstruct = 1.3, OR.obstruct.stage = 1.0,
                    HR.stage = 4.0, HRs.adj = exp(true.adj.coef))

# Scenario 6
set.seed(939006)  # Set seed to ensure reproducibility
simulations.ex1.fct(n.reps = 1000, colon.dat = colon3, HR.obstruct = 1.0, OR.obstruct.stage = 2.0,
                    HR.stage = 4.0, HRs.adj = exp(true.adj.coef))

# Scenario 7
set.seed(939007)  # Set seed to ensure reproducibility
simulations.ex1.fct(n.reps = 1000, colon.dat = colon3, HR.obstruct = 1.3, OR.obstruct.stage = 2.0,
                    HR.stage = 4.0, HRs.adj = exp(true.adj.coef))

# Scenario 8
set.seed(939008)  # Set seed to ensure reproducibility
simulations.ex1.fct(n.reps = 1000, colon.dat = colon3, HR.obstruct = 1.3, OR.obstruct.stage = 1.0,
                    HR.stage = 1.0, HRs.adj = exp(true.adj.coef))


##################################################################################################
# RESULTS: Table 1 of the manuscript
# [Step 7 (Summarizing results) of the proposed approach to data-driven simulations]
##################################################################################################

## List of all result files to be loaded

sc.toLoad <- c("./SimResults_sc1.RData",
               "./SimResults_sc2.RData",
               "./SimResults_sc3.RData",
               "./SimResults_sc4.RData",
               "./SimResults_sc5.RData",
               "./SimResults_sc6.RData",
               "./SimResults_sc7.RData",
               "./SimResults_sc8.RData")

## Create Table 1 by calculating performance measures for obstruction exposure for each scenario 
## in a loop

Table1 <- NULL

# The loop is executed for each scenario
for (i in 1:length(sc.toLoad)) {
  
  # Load results for the current scenario
  load(sc.toLoad[i])

  # Extract obstruction estimates for all simulation repetitions, which are at the position [1, 1] of  
  # the result matrices saved in the lists 'cox.wStage.coef' and 'cox.woStage.coef' respectively for 
  # the Cox model WITH stage (oracle data) and withOUT stage (imperfect data)
  cox.wStage.est.obstruct <- sapply(1:length(cox.wStage.coef), 
                                     function(i) cox.wStage.coef[[i]][1, 1])
  cox.woStage.est.obstruct <- sapply(1:length(cox.woStage.coef), 
                                      function(i) cox.woStage.coef[[i]][1, 1])

  # Extract standard errors for obstruction estimates for all simulation repetitions, which are at   
  # the position [1, 2] of the result matrices saved in the lists 'cox.wStage.coef' and 
  # 'cox.woStage.coef'
  cox.wStage.SE.obstruct <- sapply(1:length(cox.wStage.coef), 
                                     function(i) cox.wStage.coef[[i]][1, 2])
  cox.woStage.SE.obstruct <- sapply(1:length(cox.woStage.coef), 
                                      function(i) cox.woStage.coef[[i]][1, 2])  
  
  # Extract p-values for obstruction estimates for all simulation repetitions, which are at the  
  # position [1, 3] of the result matrices saved in the lists 'cox.wStage.coef' and 'cox.woStage.coef'
  cox.wStage.pval.obstruct <- sapply(1:length(cox.wStage.coef), 
                                   function(i) cox.wStage.coef[[i]][1, 3])
  cox.woStage.pval.obstruct <- sapply(1:length(cox.woStage.coef), 
                                    function(i) cox.woStage.coef[[i]][1, 3]) 

  # Calculate the performance measures for each of the two models estimated at step 6. Definitions 
  # of abbreviations: 
      # wStage, wS: Cox model WITH stage (oracle data)
      # woStage, woS: Cox model withOUT stage (imperfect data)
  cox.wStage.resu <- data.frame(wS.bias = bias.CIsign.fct(est = cox.wStage.est.obstruct, 
                                                           true.beta = log(HR.obstruct)),
                                wS.relBias = relBias.fct(est = cox.wStage.est.obstruct, 
                                                         true.beta = log(HR.obstruct)),
                                wS.empSE = round(sd(cox.wStage.est.obstruct), digits = 3),
                                wS.RMSE = rmse.fct(est = cox.wStage.est.obstruct, 
                                                   true.beta = log(HR.obstruct)),
                                wS.coverage = coverage.fct(est = cox.wStage.est.obstruct, 
                                                           se = cox.wStage.SE.obstruct,
                                                           true.beta = log(HR.obstruct)))

  cox.woStage.resu <- data.frame(woS.bias = bias.CIsign.fct(est = cox.woStage.est.obstruct, 
                                                            true.beta = log(HR.obstruct)),
                                 woS.relBias = relBias.fct(est = cox.woStage.est.obstruct, 
                                                          true.beta = log(HR.obstruct)),
                                 woS.empSE = round(sd(cox.woStage.est.obstruct), digits = 3),
                                 woS.RMSE = rmse.fct(est = cox.woStage.est.obstruct, 
                                                     true.beta = log(HR.obstruct)),
                                 woS.coverage = coverage.fct(est = cox.woStage.est.obstruct, 
                                                             se = cox.woStage.SE.obstruct, 
                                                             true.beta = log(HR.obstruct)))

  # Calculate either type I error rate or power, depending on settings of the current scenario 
  if (log(HR.obstruct) == log(1.0)){ 
    cox.wStage.resu$wS.typeI <- power.fct(pval = cox.wStage.pval.obstruct)
    cox.wStage.resu$wS.power <- ""     
    cox.woStage.resu$woS.typeI <- power.fct(pval = cox.woStage.pval.obstruct)
    cox.woStage.resu$wS.power <- ""   
  } else {
    cox.wStage.resu$wS.typeI <- ""
    cox.wStage.resu$wS.power <- power.fct(pval = cox.wStage.pval.obstruct)     
    cox.woStage.resu$woS.typeI <- ""
    cox.woStage.resu$wS.power <- power.fct(pval = cox.woStage.pval.obstruct)       
  }  
  
  # Add the input for the current scenario to Table 1
  info <- data.frame(Sc = i, 
                     HR.expo = HR.obstruct, 
                     logHR.expo = noquote(paste0("[", round(log(HR.obstruct), digits = 3), "]")),
                     OR.expo.stage = OR.obstruct.stage,
                     HR.stage = HR.stage)
  Table1 <- rbind(Table1, cbind(info, cox.wStage.resu, cox.woStage.resu))
} 

print(Table1, right = FALSE, row.names = FALSE) 

save(Table1, file = "./Table1.RData")


##################################################################################################
# VERIFICATION OF SIMULATION RESULTS
##################################################################################################

# It is essential to verify that the simulation results saved do not contain errors. This   
# verification could identify errors in the code. Here we are focusing on results for obstruction
# (exposure). Specifically, we are looking at results to confirm that:
  # - there are no missing values for estimates, standard errors, and p-values saved
  # - estimate values are within a reasonable range 
  # - all standard error values are greater than 0
  # - all p-values are between 0 and 1
  # - the number of events must be 441 for all simulated samples

## List of all result files to be loaded

sc.toLoad <- c("./SimResults_sc1.RData",
               "./SimResults_sc2.RData",
               "./SimResults_sc3.RData",
               "./SimResults_sc4.RData",
               "./SimResults_sc5.RData",
               "./SimResults_sc6.RData",
               "./SimResults_sc7.RData",
               "./SimResults_sc8.RData")

## Summary statistics for estimates, standard errors, and p-values for the variable 'obstruct' are 
## printed for each model estimated, for each scenario

for (i in 1:length(sc.toLoad)){
  load(sc.toLoad[i])

  # Extract obstruction estimates for all simulation repetitions, which are at the position [1, 1] of  
  # the result matrices saved in the lists 'cox.wStage.coef' and 'cox.woStage.coef' respectively for 
  # the Cox model WITH stage (oracle data) and withOUT stage (imperfect data)
  cox.wStage.est.obstruct <- sapply(1:length(cox.wStage.coef), 
                                     function(i) cox.wStage.coef[[i]][1, 1])
  cox.woStage.est.obstruct <- sapply(1:length(cox.woStage.coef), 
                                      function(i) cox.woStage.coef[[i]][1, 1])

  # Extract standard errors for obstruction estimates for all simulation repetitions, which are at   
  # the position [1, 2] of the result matrices saved in the lists 'cox.wStage.coef' and 
  # 'cox.woStage.coef'
  cox.wStage.SE.obstruct <- sapply(1:length(cox.wStage.coef), 
                                     function(i) cox.wStage.coef[[i]][1, 2])
  cox.woStage.SE.obstruct <- sapply(1:length(cox.woStage.coef), 
                                      function(i) cox.woStage.coef[[i]][1, 2])  
  
  # Extract p-values for obstruction estimates for all simulation repetitions, which are at the  
  # position [1, 3] of the result matrices saved in the lists 'cox.wStage.coef' and 'cox.woStage.coef'
  cox.wStage.pval.obstruct <- sapply(1:length(cox.wStage.coef), 
                                   function(i) cox.wStage.coef[[i]][1, 3])
  cox.woStage.pval.obstruct <- sapply(1:length(cox.woStage.coef), 
                                    function(i) cox.woStage.coef[[i]][1, 3]) 
  
  cat('\n\n', '*******  SCENARIO ', sc.no, '*******\n')

  cat('\n*** Obstruction estimates (true value = ', round(log(HR.obstruct), digits = 4), '):\n\n',
      ' cox.wStage.est.obstruct', '\n', sep = '')
  print(summary(cox.wStage.est.obstruct))
  cat('\n', 'cox.woStage.est.obstruct', '\n')
  print(summary(cox.woStage.est.obstruct))
  
  cat('\n*** Obstruction SEs:\n\n', 'cox.wStage.SE.obstruct', '\n')
  print(summary(cox.wStage.SE.obstruct))
  cat('\n', 'cox.woStage.SE.obstruct', '\n')
  print(summary(cox.woStage.SE.obstruct))
  
  cat('\n*** Obstruction p-values:\n\n', 'cox.wStage.pval.obstruct', '\n')
  print(summary(cox.wStage.pval.obstruct))
  cat('\n', 'cox.woStage.pval.obstruct', '\n')
  print(summary(cox.woStage.pval.obstruct))
  
  cat('\n*** Number of events (must always be 441):\n\n', 'n.events', '\n')
  print(summary(n.events))
}

