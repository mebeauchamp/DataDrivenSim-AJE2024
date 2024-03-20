##################################################################################################
##################################################################################################
###
### Code to run Example 2 of data-driven simulations, for assessing the impact of imprecise timing
### of an event on its association with a time-varying exposure (see the manuscript for details) 
###
### Code by: Marie-Eve Beauchamp                                            
### Last update: March 19, 2024                                             
###
##################################################################################################
##################################################################################################

#rm(list=ls())

## Set the path to the directory for Example 2 (need to be changed by users)
setwd("C:/.../Example2/")

## Source the script containing the functions needed for the simulations and to produce results
source('./Code_ex2_functions.R')

## Source the script containing the function 'permalgorithm.realdat', which is the adaptation of 
  # the permutational algorithm to input data including time-dependent covariates for which the 
  # values are known up to different follow-Up times across subjects
source('./Code_permalgorithm.realdat.R')


##################################################################################################
# DATA PREPARATION AND PRE-ANALYSIS BEFORE SIMULATIONS
##################################################################################################

#--------------------------
# Load the original dataset 
#--------------------------

load(file = "./dataOriginal_ex2.RData")

  # The object 'data.ori.evEnd' is the (synthetic) original dataset with each event time imputed  
  # at the end of the corresponding interval between the clinic visits immediately before and after  
  # the true unknown event time.

  # Variable description:
    # id: patient identifier
    # event: event status (1=event, 0=censoring)
    # start: start time of the interval  
    # stop: stop time of the interval 
    # anyUseLast14: time-varying binary indicator of benzodiazepine use in last 14 days (exposure)  
    # sex: patient's sex  
    # age: patient's age

head(data.ori.evEnd)
tail(data.ori.evEnd)

# Number of patients  
(N <- length(unique(data.ori.evEnd$id)))

# Number of events
(n.events.true <- sum(data.ori.evEnd$event))
(n.events.true / N) * 100

# Object 'visits' is a matrix with each row representing all visit times for one patient in 
# 'data.ori.evEnd'. First column indicates the patient's id. Columns 2-14 are the visit times (up 
# to 13 visits, with NA assigned to remaining columns after the last visit for a patient). 

dim(visits)
visits[1:15, ]
table(visits[, 'visit1']) # First visit is at time 0

# Average time between the visits for a patient (in months)
mean(t(apply(visits[, 2:14], 1, diff, lag = 1)) / 30, na.rm = TRUE) 

#---------------------------------------------------------------------------------------------
# Impute event times at the mid-point of between-visits intervals in which the events occurred 
#---------------------------------------------------------------------------------------------

# Create a list in which each item is the data for one patient in the original dataset
data.ori.evEnd.list <- split(data.ori.evEnd, data.ori.evEnd$id)

# Apply the function 'impute.evEndToMid.fct' to each patient, i.e. to each item of the list of 
# patients 'data.ori.evEnd.list' and using the visit times for the corresponding patient (note 
# that the 1st column of the object visits is the patients' identification number and is then 
# excluded when passing the visit times as an argument)
data.ori.evMid.tmp <- lapply(1:N, function(i) 
                                    impute.evEndToMid.fct(dat.evEnd.subj = data.ori.evEnd.list[[i]],
                                                          visits.subj = visits[i, -1]))  

# Convert the list of results (each item is for one patient) to a data frame
data.ori.evMid <- do.call("rbind", data.ori.evMid.tmp)

#-------------------------------------------------------------------------------------------------
# Identify the interval from the last visit before and first visit after each true (unknown) event 
# time in the original dataset
#-------------------------------------------------------------------------------------------------

# Extract the identification number and event time for patients with an event
ev.id.time <- data.ori.evMid[data.ori.evMid$event == 1, c('id', 'stop')]

# Extract all visit times for patients with an event 
visits.ev <- visits[visits[, 'id'] %in% ev.id.time[, 'id'], ]

# Create a matrix in which the 1st and 2nd column is respectively the visit time before and after
# a true event time, by applying the function 'interval.ev.fct' to each event time
interval.ev <- do.call('rbind', lapply(1:nrow(ev.id.time), function(i) 
                                       interval.ev.fct(ev.time = ev.id.time[i, 'stop'], 
                                                       visits.subj = visits.ev[i, ])))

#---------------------------------------------------------------------------------------------
# Analyze the original dataset with each (unknown) event time imputed at either: 
#   (i) the mid-point of the between-visits interval in which the event occurred, or 
#   (ii) the end of this interval, i.e. at the visit when the event was first reported 
# [Step 2 (Initial analyses not corrected for the imperfection) of the proposed approach to
# data-driven simulations]
#---------------------------------------------------------------------------------------------

## Descriptive statistics on the length of between-visits intervals (in days)

round(mean(interval.ev[, 2] - interval.ev[, 1]), digits = 1)
quantile(interval.ev[, 2] - interval.ev[, 1])

## Cox proportional hazards (PH) model with event times imputed at (i) the mid-point of the   
## between-visits interval in which the event occurred

library(survival)
cox.ori.evMid <- summary(coxph(Surv(start, stop, event) ~ anyUseLast14 + sex + age, 
                               data = data.ori.evMid))

# Hazard ratio (HR) for recent benzodiazepine exposure, adjusted for age and sex
round(cox.ori.evMid$conf.int[1, c(1, 3:4)], digits = 2)

## Cox PH model with event times impute at (ii) the end of the between-visits interval, i.e.    
## at the visit when the event was first reported

cox.ori.evEnd <- summary(coxph(Surv(start, stop, event) ~ anyUseLast14 + sex + age, 
                                data = data.ori.evEnd))

# HR for recent benzodiazepine exposure, adjusted for age and sex
round(cox.ori.evEnd$conf.int[1, c(1, 3:4)], digits = 2)


##################################################################################################
# SIMULATIONS  
# [Steps 4-6 of the proposed approach to data-driven simulations]
##################################################################################################

# NOTE: Steps 4-6 are performed by the function 'simulations.ex2.fct'. For details see the script 
# 'Code_ex2_functions.R', which includes the code of 'simulations.ex2.fct' and the definition of 
# the function arguments.

## Run the simulations separately for each simulation scenario

options(warn = 1)  # To request warnings to be printed as they occur

## Scenario: HR for exposure of 1.0
set.seed(860010)  # Set seed to ensure reproducibility
simulations.ex2.fct(n.reps = 1000, dat.ori.evEnd = data.ori.evEnd, visits = visits,
                    interval.ev = interval.ev, beta.expo = log(1.0), 
                    betas.covar = cox.ori.evMid$coef[2:3, 1])

## Scenario: HR for exposure of 1.5
set.seed(860015)  # Set seed to ensure reproducibility
simulations.ex2.fct(n.reps = 1000, dat.ori.evEnd = data.ori.evEnd, visits = visits,
                    interval.ev = interval.ev, beta.expo = log(1.5), 
                    betas.covar = cox.ori.evMid$coef[2:3, 1])

## Scenario: HR for exposure of 2.0
set.seed(860020)  # Set seed to ensure reproducibility
simulations.ex2.fct(n.reps = 1000, dat.ori.evEnd = data.ori.evEnd, visits = visits,
                    interval.ev = interval.ev, beta.expo = log(2.0), 
                    betas.covar = cox.ori.evMid$coef[2:3, 1])

## Scenario: HR for exposure of 2.5
set.seed(860025)  # Set seed to ensure reproducibility
simulations.ex2.fct(n.reps = 1000, dat.ori.evEnd = data.ori.evEnd, visits = visits,
                    interval.ev = interval.ev, beta.expo = log(2.5), 
                    betas.covar = cox.ori.evMid$coef[2:3, 1])


##################################################################################################
# RESULTS
# [Step 7 (Summarizing results) of the proposed approach to data-driven simulations]
##################################################################################################

## List of all result files to be loaded

sc.toLoad <- c("./SimResults_sc1.RData",
               "./SimResults_sc2.RData",
               "./SimResults_sc3.RData",
               "./SimResults_sc4.RData")

## Create Table 2 and Web Table 3 by calculating the performance measures for exposure for each 
## scenario in a loop

Table2 <- NULL
WebTable3 <- NULL

# The loop is executed for each scenario
for (i in 1:length(sc.toLoad)) {
  
  # Load results for the current scenario
  load(sc.toLoad[i])
  
  # Extract exposure estimates for all simulation repetitions, which are at the position [1, 1] of   
  # the result matrices saved in the lists 'cox.gen.est.expo', 'cox.mid.est.expo' and 
  # 'cox.end.est.expo', respectively for Cox PH models estimated on simulated data where event 
  # times were as generated, or imputed at the mid-point or end-point of the intervals between the  
  # clinic visits immediately before and after the generated event times
  cox.gen.est.expo <- sapply(1:length(cox.gen), function(i) cox.gen[[i]][1, 1])
  cox.mid.est.expo <- sapply(1:length(cox.mid), function(i) cox.mid[[i]][1, 1])
  cox.end.est.expo <- sapply(1:length(cox.end), function(i) cox.end[[i]][1, 1])
  
  # Extract standard errors for exposure estimates for all simulation repetitions, which are at 
  # the position [1, 2] of the result matrices saved in the lists 'cox.gen.est.expo', 
  # 'cox.mid.est.expo' and 'cox.end.est.expo'
  cox.gen.SE.expo <- sapply(1:length(cox.gen), function(i) cox.gen[[i]][1, 2])
  cox.mid.SE.expo <- sapply(1:length(cox.mid), function(i) cox.mid[[i]][1, 2])
  cox.end.SE.expo <- sapply(1:length(cox.end), function(i) cox.end[[i]][1, 2])
  
  # Calculate the performance measures for each of the three Cox PH models estimated at step 6 
  Model1.resu <- data.frame(M1.bias = bias.CIsign.fct(est = cox.gen.est.expo, 
                                                      true.beta = beta.expo),
                            M1.relBias = relBias.fct(est = cox.gen.est.expo, 
                                                     true.beta = beta.expo),
                            M1.empSE = round(sd(cox.gen.est.expo), digits = 3),
                            M1.RMSE = rmse.fct(est = cox.gen.est.expo, 
                                           true.beta = beta.expo),
                            M1.coverage = coverage.fct(est = cox.gen.est.expo, 
                                                       SE = cox.gen.SE.expo,
                                                       true.beta = beta.expo))
  
  Model2.resu <- data.frame(M2.bias = bias.CIsign.fct(est = cox.mid.est.expo, 
                                                      true.beta = beta.expo),
                            M2.relBias = relBias.fct(est = cox.mid.est.expo, 
                                                     true.beta = beta.expo),
                            M2.empSE = round(sd(cox.mid.est.expo), digits = 3),
                            M2.RMSE = rmse.fct(est = cox.mid.est.expo, 
                                               true.beta = beta.expo),
                            M2.coverage = coverage.fct(est = cox.mid.est.expo, 
                                                       SE = cox.mid.SE.expo,
                                                       true.beta = beta.expo))
  
  Model3.resu <- data.frame(M3.bias = bias.CIsign.fct(est = cox.end.est.expo, 
                                                  true.beta = beta.expo),
                            M3.relBias = relBias.fct(est = cox.end.est.expo, 
                                                     true.beta = beta.expo),
                            M3.empSE = round(sd(cox.end.est.expo), digits = 3),
                            M3.RMSE = rmse.fct(est = cox.end.est.expo, 
                                               true.beta = beta.expo),
                            M3.coverage = coverage.fct(est = cox.end.est.expo, 
                                                       SE = cox.end.SE.expo,
                                                       true.beta = beta.expo))
  
  # Ratios and comparison of performance measures between Model 2 and Model 3. Definitions of 
  # statistics calculated:
      # M3oM2.bias: ratio of bias Model 3 over Model 2  
      # M3oM2.RMSE: ratio of RMSE Model 3 over Model 2  
      # M2closerM3: % repetitions Model 2 closer to truth than Model 3
  if (beta.expo == log(1.0)){  # No ratios calculated when true HR for exposure is 1.0  
    comp.Model2.Model3 <- data.frame(M3oM2.bias = 'N/A',
                                     M3oM2.RMSE = round(Model3.resu$M3.RMSE / 
                                                          Model2.resu$M2.RMSE, digits = 2),
                                     M2closerM3 = closer.fct(estA = cox.mid.est.expo, 
                                                             estB = cox.end.est.expo, 
                                                             true.beta = beta.expo))
  } else {
    M2.bias <- bias.fct(est = cox.mid.est.expo, true.beta = beta.expo)
    M3.bias <- bias.fct(est = cox.end.est.expo, true.beta = beta.expo)
    comp.Model2.Model3 <- data.frame(M3oM2.bias = round(M3.bias / M2.bias, digits = 2),
                                     M3oM2.RMSE = round(Model3.resu$M3.RMSE / 
                                                          Model2.resu$M2.RMSE, digits = 2),
                                     M2closerM3 = closer.fct(estA = cox.mid.est.expo, 
                                                             estB = cox.end.est.expo, 
                                                             true.beta = beta.expo))
  }
  
  # Add the input for the current scenario to Table 2 and Web Table 3
  info <- data.frame(Scenario = i, 
                     TrueHR = exp(beta.expo), 
                     logHR = noquote(paste0("[", round(beta.expo, digits = 3), "]")))
  Table2 <- rbind(Table2, cbind(info, Model2.resu, Model3.resu, comp.Model2.Model3))
  WebTable3 <- rbind(WebTable3, cbind(info, Model1.resu))
} 

print(Table2, right = FALSE, row.names = FALSE) 
print(WebTable3, right = FALSE, row.names = FALSE) 

save(Table2, file = "./Table2.RData")
save(WebTable3, file = "./WebTable3.RData")

## Additional results presented

# Percentage of simulation repetitions for which Model 2 exposure estimate in scenario 1 was higher
# than step 2 initial estimate with event times imputed at mid-point of intervals (HR=1.47)  

load(sc.toLoad[1])
cox.mid.est.expo <- sapply(1:length(cox.mid), 
                           function(i) cox.mid[[i]][1, 1])
(sum(cox.mid.est.expo > 0.38667) / n.reps) * 100  # Note exp(0.38667) = 1.47

# Percentage of simulation repetitions for which Model 3 exposure CI includes HR=1.0 in scenario 2 
# (i.e. a p-value > 0.05)

load(sc.toLoad[2])
cox.end.pval.expo <- sapply(1:length(cox.end), 
                           function(i) cox.end[[i]][1, 3])
(sum(cox.end.pval.expo > 0.05) / n.reps) * 100  


##################################################################################################
# VERIFICATION OF SIMULATION RESULTS
##################################################################################################

# It is essential to verify that the simulation results saved do not contain errors. This   
  # verification could identify errors in the code. Here we are focusing on results for exposure.
  # Specifically, we are looking at results to confirm that:
  # - there are no missing values for estimates and standard errors saved
  # - estimate values are within a reasonable range 
  # - all standard error values are greater than 0
  # - the number of events must be 285 for all simulated samples

## List of all result files to be loaded

sc.toLoad <- c("./SimResults_sc1.RData",
               "./SimResults_sc2.RData",
               "./SimResults_sc3.RData",
               "./SimResults_sc4.RData")

## Summary statistics for estimates and standard errors for exposure are printed for each model
## estimated, for each scenario

for (i in 1:length(sc.toLoad)) {
  load(sc.toLoad[i])
  
  # Extract exposure estimates for all simulation repetitions, which are at the position [1, 1] of   
  # the result matrices saved in the lists 'cox.gen.est.expo', 'cox.mid.est.expo' and 
  # 'cox.end.est.expo', respectively for Cox PH models estimated on simulated data with event times
  # as they were generated, or imputed at the mid-point or end-point of the intervals between the  
  # clinic visits immediately before and after the generated event times
  cox.gen.est.expo <- sapply(1:length(cox.gen), 
                             function(i) cox.gen[[i]][1, 1])
  cox.mid.est.expo <- sapply(1:length(cox.mid), 
                             function(i) cox.mid[[i]][1, 1])
  cox.end.est.expo <- sapply(1:length(cox.end), 
                             function(i) cox.end[[i]][1, 1])
  
  # Extract standard errors for exposure estimates for all simulation repetitions, which are at 
  # the position [1, 2] of the result matrices saved in the lists 'cox.gen.est.expo', 
  # 'cox.mid.est.expo' and 'cox.end.est.expo'
  cox.gen.SE.expo <- sapply(1:length(cox.gen), 
                            function(i) cox.gen[[i]][1, 2])
  cox.mid.SE.expo <- sapply(1:length(cox.mid), 
                            function(i) cox.mid[[i]][1, 2])
  cox.end.SE.expo <- sapply(1:length(cox.end), 
                            function(i) cox.end[[i]][1, 2])
  
  cat('\n\n', '***** SCENARIO: TRUE HR EXPOSURE =', exp(beta.expo), '*****\n')
  
  cat('\nExposure estimates (true value = ', round(beta.expo, digits = 4), '):\n', 
      'cox.gen.est.expo', '\n', sep = '')
  print(summary(cox.gen.est.expo))
  cat('\n', 'cox.mid.est.expo', '\n')
  print(summary(cox.mid.est.expo))
  cat('\n', 'cox.end.est.expo', '\n')
  print(summary(cox.end.est.expo))
  
  cat('\nExposure SEs:\n', 'cox.gen.SE.expo', '\n')
  print(summary(cox.gen.SE.expo))
  cat('\n', 'cox.mid.SE.expo', '\n')
  print(summary(cox.mid.SE.expo))
  cat('\n', 'cox.end.SE.expo', '\n')
  print(summary(cox.end.SE.expo))
  
  cat('\nNumber of events (must always be 285):\n', 'n.events.gen', '\n')
  print(summary(n.events.gen))
  cat('\n', 'n.events.mid', '\n')
  print(summary(n.events.mid))
  cat('\n', 'n.events.end', '\n')
  print(summary(n.events.end))
}


##################################################################################################
# ADDITIONAL SENSITIVITY ANALYSES  
# Quantitative bias analysis (QBA)-like sensitivity analyses of the original dataset
##################################################################################################

options(warn = 1)  # To request warnings to be printed as they occur

## Scenario: QBA

set.seed(2860010)  # Set seed to ensure reproducibility
simulations.ex2.QBA.fct(n.reps = 1000, dat.ori.evEnd = data.ori.evEnd, interval.ev = interval.ev)

## Results

# Load results 
load(c("./SimResults_QBA.RData"))

# Extract exposure estimates for all simulation repetitions, which are at the position [1, 1] of   
# the result matrices saved in the list 'cox.gen.est.expo'
cox.gen.est.expo <- sapply(1:length(cox.gen), function(i) cox.gen[[i]][1, 1])

# Extract standard errors for exposure estimates for all simulation repetitions, which are at 
# the position [1, 2] of the result matrices saved in the list 'cox.gen.est.expo'
cox.gen.SE.expo <- sapply(1:length(cox.gen), function(i) cox.gen[[i]][1, 2])

# Mean of the 1,000 log(HR) estimates and the corresponding hazard ratio
round(mean(cox.gen.est.expo), digits = 3)
round(exp(mean(cox.gen.est.expo)), digits = 2)

# Empirical standard error of log(HR) estimates
round(sd(cox.gen.est.expo), digits = 3)

## Verification of simulation results
summary(n.events.gen)
summary(cox.gen.est.expo)

