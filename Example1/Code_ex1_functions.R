##################################################################################################
##################################################################################################
###
### Functions used in the script 'Code_ex1_simulations.R', which runs Example 1 of data-driven   
### simulations, for assessing the impact of omitting an important prognostic factor (see the
### manuscript for details)
###
### Code by: Marie-Eve Beauchamp                                            
### Last update: March 20, 2024                                              
###
##################################################################################################
##################################################################################################

library(survival)
library(PermAlgo)

##################################################################################################
# simulations.ex1.fct:                                                   
#   Function to generate the oracle and imperfect datasets, analyze these datasets, and save     
#   results for all repetitions of a given simulation scenario for Example 1.                        
##################################################################################################

## Arguments:
  # n.reps: number of repetitions per scenario.
  # colon.dat: final colon dataset prepared for the simulations, including the columns for the  
  #            time and event status, obstruction exposure and adjustment covariates, with these  
  #            exact names respectively: 'time', 'status', 'obstruct', 'rx_Lev', 'rx_Lev5FU',  
  #            'sex', 'age', 'perfor', 'adhere', 'differ_mod', 'differ_poor', 'extent_muscle', 
  #            'extent_serosa', 'extent_cs', 'surg', 'node4'. See the script 
  #            'Code_ex1_simulations.R' for the preparation of the final colon dataset. 
  # HR.obstruct: hazard ratio value for 'obstruct' (exposure) in the true Cox model for data  
  #              generation, for the current simulation scenario.
  # OR.obstruct.stage: odds ratio value for 'obstruct' in the logistic model to generate values 
  #                    for 'stage', for the current simulation scenario.
  # HR.stage: hazard ratio value for 'stage' in the true Cox model for data generation, for the  
  #           current simulation scenario.
  # HRs.adj: vector of hazard ratio values for the true Cox model for data generation, respectively   
  #          for the adjustment variables: 'rx_Lev', 'rx_Lev5FU', 'sex', 'age', 'perfor',     
  #          'adhere', 'differ_mod', 'differ_poor', 'extent_muscle', 'extent_serosa', 'extent_cs',    
  #          'surg', and 'node4' (with these exact names for the items in HRs.adj, i.e. that names   
  #          must match the names in 'colon.dat').

simulations.ex1.fct <- function(n.reps, colon.dat, HR.obstruct, OR.obstruct.stage, HR.stage, 
                                HRs.adj){

  #--------------------
  # Objects to be saved 
  #--------------------

  # Number of events in the generated oracle datasets
  n.events <- rep(NA, times = n.reps) 
  
  # Results from Cox models with and without 'stage' for the analysis of simulated datasets   
  cox.wStage.coef <- list()
  cox.woStage.coef <- list()

  #---------------------------------------------------------------------------------------------
  # Generation of oracle datasets 
  # [Step 4 (Oracle dataset generation) of the proposed approach to data-driven simulations]
  #---------------------------------------------------------------------------------------------   
  
  ## Extract observed values of colon obstruction exposure and all measured covariates from the  
  ## real-world colon data
  ## [Sub-step 4a]  
  
    # Note: this simply implies to use the corresponding columns from the object 'colon.dat'  
    # passed to the function 'simulations.ex1.fct'
  
  ## Extract the observed real-world outcomes (441 event and 465 censoring times) from the  
  ## real-world colon data
  ## [Sub-step 4b]
  
    # Note: event and censoring times (objects 'eventtimes' and 'censotimes') are prepared below
    # according to the format required by the 'permalgorithm' function from the 'PermAlgo' package,
    # which is used below to assign each outcome to a patient with the permutational algorithm. 
    # This function requires the arguments 'eventRandom' and 'censorRandom', both as vectors with a
    # length equal to the number of patients in the dataset. These arguments represent respectively
    # the potential event and censoring times. The minimum value for each pair of event and 
    # censoring times from the two vectors is then assigned to a patient (see details at 
    # https://cran.r-project.org/web/packages/PermAlgo/PermAlgo.pdf). 
    # For our simulation study, we want to ensure that the observed 441 event and 465 censoring   
    # times from the real-world data are the outcomes assigned. Therefore, for each truly observed   
    # event time, we indicate this value in the vector 'eventtimes' and we assign a mock greater     
    # value at the same position in the vector 'censotimes' (e.g., greater than the maximum    
    # follow-up time, i.e. 'max(colon.dat$time) + 1'). Consequently, the observed event time will  
    # be selected for that pair as the outcome value to assign to a patient. And vice versa for  
    # each truly observed censoring time in 'censotimes', a mock greater event time is assigned
    # in 'eventtimes' for that pair. 

  eventtimes <- ifelse(colon.dat$status==1, colon.dat$time, max(colon.dat$time) + 1)
  censotimes <- ifelse(colon.dat$status==0, colon.dat$time, max(colon.dat$time) + 1)  
  
  ## Loop repeated for each of the 'n.reps' repetitions

  cat('Simulation repetition:\n')
  for (k in 1:n.reps){
    cat('  k=', k, '\n')

    ## Generate a binary indicator of higher cancer stage for each patient from the binomial 
    ## distribution
    ## [Sub-step 4c]

    logitp <- log(0.2) + (log(OR.obstruct.stage) * colon.dat$obstruct) + 
      (log(2.5) * colon.dat$node4) + (log(1.5) * colon.dat$extent_muscle) + 
      (log(2.0) * colon.dat$extent_serosa) + (log(2.5) * colon.dat$extent_cs) 
    p.stage <- exp(logitp) / (1 + exp(logitp))  # Prob. of higher cancer stage for each patient
    #mean(p.stage)
    colon.dat$stage <- rbinom(nrow(colon.dat), 1, p.stage) 

    ## Preparation of exposure and covariates data for the permutational algorithm

      # Note: the 'permalgorithm' function in the 'PermAlgo' package allows to assign event times  
      # conditional on time-dependent and/or time-invariant covariates. This function requires as
      # an argument a matrix of covariate values ('Xmat') in the counting process format where    
      # each line represents the exposure/covariate values for a patient for one time unit. 
      # Exposure/covariate values must be provided for each time unit up to the maximum follow-up    
      # time in the data. Given that all covariates are time-invariant in the current simulation  
      # study, the lines of 'Xmat' for a patient are the repetition of her/his vector of covariate  
      # values.

    # Select from the real-world colon dataset the relevant columns, i.e. 'obstruct' and covariates  
    # for adjustment, and add the generated 'stage'
    Xmat.colon.dat <- colon.dat
    Xmat.colon.dat <- Xmat.colon.dat[, c('obstruct', names(HRs.adj), 'stage')]
    
    # Repeat each line (exposure and covariate values for a patient) a number of times equal to
    # the maximum follow-up in the colon data
    Xmat.colon.dat <- Xmat.colon.dat[rep(1:nrow(Xmat.colon.dat), each = max(colon.dat$time)), ]
      #dim(colon.dat)                      
      #dim(Xmat.colon.dat)  
      #nrow(colon.dat) * max(colon.dat$time)  # Must match the number of rows of Xmat.colon.dat
    
    # Transform to a matrix (requirement for the 'permalgorithm' function)
    Xmat.colon.dat <- as.matrix(Xmat.colon.dat)

    ## Assign each event and censoring times to a patient (see details of the function at 
    ## https://cran.r-project.org/web/packages/PermAlgo/PermAlgo.pdf)
    ## [Sub-step 4d]
    
    dat.oracle <- permalgorithm(numSubjects = nrow(colon.dat), 
                             maxTime = max(colon.dat$time), 
                             Xmat = Xmat.colon.dat,  XmatNames = colnames(Xmat.colon.dat),
                             eventRandom = eventtimes, 
                             censorRandom = censotimes,
                             betas = log(c(HR.obstruct, HRs.adj, HR.stage)))
    
      # Note: the output of the 'permalgorithm' function is a time-dependent dataset, with time
      # intervals indicated by columns 'Start' and 'Stop', and the outcome indicated by the column 
      # 'Event'. Note that the new 'Id' values in the resulting dataset does not match 'id' in 
      # 'colon.dat'.

    # Select only the last line of each patient because the data are actually time-invariant
    dat.oracle.lastLine <- dat.oracle[!duplicated(dat.oracle$Id, fromLast = T), 
                                      -which(colnames(dat.oracle) == 'Start')]
    
      # Note: the column 'Stop' corresponds to the event or censoring time for the patient. 
      # Column 'Start' was removed to avoid confusion.
    
    # Save the number of events assigned in this generated datasets for verification (must always 
    # be 441)
    n.events[k] <- sum(dat.oracle.lastLine$Event)
    
    #------------------------------------------------------------------------------------------------
    # Creation of the imperfect dataset by removing the column 'stage'
    # [Step 5 (Imperfect dataset generation)] of the proposed approach to data-driven simulations]
    #------------------------------------------------------------------------------------------------
    
    dat.imperfect.lastLine <- dat.oracle.lastLine[, -which(colnames(dat.oracle.lastLine) == 'stage')]
    
    #--------------------------------------------------------------------------------------------
    # [Step 6 (Analyses) of the proposed approach to data-driven simulations]
    #--------------------------------------------------------------------------------------------
   
    cox.wStage.coef[[k]] <- summary(coxph(Surv(Stop, Event) ~ obstruct + rx_Lev + rx_Lev5FU + 
                                            sex + age + perfor + adhere + differ_mod + 
                                            differ_poor + extent_muscle + extent_serosa + 
                                            extent_cs + surg + node4 + stage, 
                                          data = dat.oracle.lastLine))$coef[, c(1, 3, 5)]
      
    cox.woStage.coef[[k]] <- summary(coxph(Surv(Stop, Event) ~ obstruct + rx_Lev + rx_Lev5FU +
                                             sex + age + perfor + adhere + differ_mod + 
                                             differ_poor + extent_muscle + extent_serosa + 
                                             extent_cs + surg + node4, 
                                           data = dat.imperfect.lastLine))$coef[, c(1, 3, 5)]
  }

  #---------------
  # Saving results
  #--------------- 

  # Scenario number 
  if (HR.obstruct == 1.0 & OR.obstruct.stage == 1.2 & HR.stage == 4.0){
    sc.no <- 1
  } else if (HR.obstruct == 1.3 & OR.obstruct.stage == 1.2 & HR.stage == 4.0){
    sc.no <- 2
  } else if (HR.obstruct == 1.5 & OR.obstruct.stage == 1.2 & HR.stage == 4.0){
    sc.no <- 3
  } else if (HR.obstruct == 2.0 & OR.obstruct.stage == 1.2 & HR.stage == 4.0){
    sc.no <- 4
  } else if (HR.obstruct == 1.3 & OR.obstruct.stage == 1.0 & HR.stage == 4.0){
    sc.no <- 5
  } else if (HR.obstruct == 1.0 & OR.obstruct.stage == 2.0 & HR.stage == 4.0){
    sc.no <- 6
  } else if (HR.obstruct == 1.3 & OR.obstruct.stage == 2.0 & HR.stage == 4.0){
    sc.no <- 7
  } else if (HR.obstruct == 1.3 & OR.obstruct.stage == 1.0 & HR.stage == 1.0){
    sc.no <- 8
  }
  
  # Scenario name
  sc.name <- paste0('HRobstruct', gsub("\\.", "", HR.obstruct), 
                    ' ORobstruct', gsub("\\.", "", OR.obstruct.stage),
                    ' HRstage', gsub("\\.", "", HR.stage))

  save(n.reps, HR.obstruct, OR.obstruct.stage, HR.stage,  sc.name, sc.no, n.events, 
       cox.wStage.coef, cox.woStage.coef, 
       file = paste0('./SimResults_sc', sc.no, '.RData'))
}


##################################################################################################
# bias.fct, bias.CIsign.fct, relBias.fct, rmse.fct, coverage.fct, power.fct:                    
#   Functions to calculate respectively:                                      
#     - Bias
#     - Bias and indicate if bias 95% confidence interval (CI) excludes 0 with symbol *
#     - Relative bias
#     - Root mean squared error (RMSE)                                        
#     - Coverage rate of the 95% CI                          
#     - Power to detect a significant effect*
#   for log(HR) estimates for a variable for one simulation scenario.
#
#   * Note the power is equivalent to type I error rate for scenarios 1 and 6, with true HR of 1.0  
#     for obstruction exposure.
##################################################################################################

## Arguments:
  # est: vector of estimates for a variable from all simulation repetitions for a simulation 
  #      scenario.
  # true.beta: true parameter value for that variable for that simulation scenario.
  # se: vector of standard errors of the estimates for that variable from all repetitions for that  
  #     simulation scenario.
  # pval: vector of p-values of the estimates for that variable from all repetitions for that 
  #       simulation scenario.

bias.fct <- function(est, true.beta) {
  return(mean(est) - true.beta)
}

bias.CIsign.fct <- function(est, true.beta) {
  bias <- bias.fct(est, true.beta)
  biasCI <- bias + c(-1, 1) * qnorm(0.975) * sqrt(var(est) / length(est))
  biasCIsign <- !(biasCI[1] < 0 & biasCI[2] > 0)  
  if (biasCIsign == TRUE){
    return(noquote(paste0(round(bias, digits = 3), "*")))  # Symbol * indicates CI excludes 0   
  } else {
    return(round(bias, digits = 3))
  }
}

relBias.fct <- function(est, true.beta) {
  if (true.beta != 0){
    return(round((bias.fct(est, true.beta) / true.beta) * 100, digits = 1))
  } else {
    return('N/A')   
  }
}

rmse.fct <- function(est, true.beta){
  round(sqrt((mean(est) - true.beta) ^ 2 + var(est)), digits = 3)
}

coverage.fct <- function(est, se, true.beta){
  lowerBound <- est - qnorm(0.975) * se 
  upperBound <- est + qnorm(0.975) * se  
  inCI <- ifelse(true.beta > lowerBound & true.beta < upperBound, 1, 0)
  return(round(sum(inCI) / length(est), digits = 3))
}

power.fct <- function(pval){
  round(mean(pval < 0.05) * 100, digits = 1)
}
