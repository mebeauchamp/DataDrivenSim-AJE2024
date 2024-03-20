##################################################################################################
##################################################################################################
###
### Functions used in the script 'Code_ex2_simulations.R', which runs Example 2 of data-driven
### simulations, for assessing the impact of imprecise timing of an event on its association 
### with a time-varying exposure (see the manuscript for details)
###
### Code by: Marie-Eve Beauchamp                                            
### Last update: March 19, 2024                                              
###
##################################################################################################
##################################################################################################

library(survival)

##################################################################################################
# impute.evEndToMid.fct:   
#   Function to change the imputation of the event time for a patient from the end-point to 
#   mid-point of the interval between the last clinic visit before and the first visit after the
#   true event time; follow-up of a patient without an event is unchanged.                                                 
##################################################################################################

## Arguments:
  # dat.evEnd.subj: data for a subject, with the event (if applicable) imputed at the end-point of    
  #                 the relevant between-visits interval.
  # visits.subj: vector of all clinic visit times for the current subject.  

impute.evEndToMid.fct <- function(dat.evEnd.subj, visits.subj){

  # Remove NA values after the last visit (if any)
  visits.subj <- visits.subj[!is.na(visits.subj)]  
  
  # If current subject does not have the event: no changes to follow-up 
  if (tail(dat.evEnd.subj$event, n = 1) == 0){
    return(dat.evEnd.subj)  
    
  # If current subject has the event
  } else {
    
    # Impute event time at mid-point between visits before and after the true event time
    imp.evt <- ceiling((max(visits.subj[visits.subj < max(dat.evEnd.subj$stop)]) + 
                        min(visits.subj[visits.subj >= max(dat.evEnd.subj$stop)])) / 2)

    # Cut follow-up at the new imputed event time
    dat.subj.evMid <- dat.evEnd.subj[dat.evEnd.subj$stop <= imp.evt, ]
    
    # Assign an event at the last line of follow-up
    dat.subj.evMid$event[length(dat.subj.evMid$event)] <- 1
      
    return(dat.subj.evMid)
  }
} 


##################################################################################################
# interval.ev.fct:                                                                       
#   Function to identify the last clinic visit before and the first visit after an event time. 
##################################################################################################

## Arguments:
  # ev.time: an event time.
  # visits.subj: vector of all clinic visit times for the subject who had the current event time.

interval.ev.fct <- function(ev.time, visits.subj){
  interval <- c(max(visits.subj[visits.subj < ev.time], na.rm = T), 
                min(visits.subj[visits.subj >= ev.time], na.rm = T))
  return(interval)
}


##################################################################################################
# simulations.ex2.fct:                                                   
#   Function to generate the oracle and imperfect datasets, analyze these datasets, and save     
#   results for all repetitions of a given simulation scenario for Example 2.                         
##################################################################################################

## Arguments:
  # n.reps: number of simulation repetitions per scenario.
  # dat.ori.evEnd: original dataset with each event times imputed at the end of the interval
  #                between the last clinic visit before and the first visit after the true  
  #                (unknown) event time. The dataset must include the following columns (with the 
  #                names):
  #                - 'start' and 'stop' indicating the beginning and end the time interval of the  
  #                  current observation, 
  #                - 'event' indicating the the event status, 
  #                - 'anyUseLast14', 'sex' and 'age' for the time-varying exposure metric (any  
  #                   benzodiazepine in the last 14 days) and the covariates.
  # visits: matrix of all visits times for all patients, with the first column indicating the 
  #         patient's id.
  # interval.ev: matrix in which the 1st and 2nd column is respectively the visit time before and 
  #              after each true event time.
  # beta.expo: parameter value for 'anyUseLast14' exposure variable in true Cox proportional  
  #            hazards (PH) model for the data generation.
  # betas.covar: vector of parameter values respectively for 'sex' and 'age' in true Cox PH model  
  #              for the data generation.

simulations.ex2.fct <- function(n.reps, dat.ori.evEnd, visits, interval.ev, beta.expo, 
                                betas.covar){
  
  #-----------------
  # Data preparation
  #-----------------
  
  # Number of patients (N) and number of events in the original dataset
  N <- length(unique(dat.ori.evEnd$id))
  n.events.true <- sum(dat.ori.evEnd$event) 
  
  # Create a list in which each item is the original data for one patient (used later in the code)
  dat.ori.evEnd.list <- split(dat.ori.evEnd, dat.ori.evEnd$id)
  
  #--------------------
  # Objects to be saved
  #--------------------
  
  # To monitor to which patients the generated event times were assigned in each generated dataset
  ev.assigned.ids <- matrix(NA, nrow = n.reps, ncol = n.events.true)  

  # Number of events in the generated datasets with the event times either as generated or imputed
  # at the middle or end of the between-visits interval   
  n.events.gen <- rep(NA, n.reps)
  n.events.mid <- rep(NA, n.reps)
  n.events.end <- rep(NA, n.reps)
  
  # Results from Cox PH models estimated on the generated datasets with the event times either as
  # generated or imputed at the middle or end of the between-visits intervals
  cox.gen <- list()
  cox.mid <- list()
  cox.end <- list()

  #---------------------------------------------------------------------------------------------
  # Generation of oracle datasets 
  # [Step 4 (Oracle dataset generation) of the proposed approach to data-driven simulations]
  #---------------------------------------------------------------------------------------------   
  
  ## Use observed values of the time-varying exposure metric and covariates of the 1,250   
  ## individuals from the original dataset
  ## [Sub-step 4a]  

    # Note: this simply implies using the corresponding columns from the object 'dat.ori.evEnd'   
    # passed to the function 'simulations.ex2.fct'
  
  ## Use the times of the 285 visits when the observed events were reported and each corresponding 
  ## previous visit time from original data
  ## [Sub-step 4b]
  
    # Note: this simply implies using the object 'interval.ev' passed to the function 
    # 'simulations.ex2.fct'  

  ## Loop repeated for each of the 'n.reps' repetitions
  
  cat('Simulation repetition:\n')
  for (k in 1:n.reps){
    cat('  k=', k, '\n')

    ## Generate true event times over the interval between the visits just before and just after 
    ## the event times from original data
    ## [Sub-step 4c]    

    ev.gen <- sort(ceiling(runif(n.events.true, interval.ev[, 1], interval.ev[, 2])))
    
    ## Assign each event and censoring times to a patient, and prepare the generated dataset with
    ## event times as generated, where:
    ##  - patients assigned an event finish their follow-up at their assigned event time,
    ##  - patients not assigned an event are censored at the time their follow-up ended in the 
    ##    original dataset 
    ## [Sub-steps 4d-4e]
    
    dat.gen.gen <- permalgorithm.realdat(data = dat.ori.evEnd[, c('id', 'start', 'stop', 
                                                                  'anyUseLast14', 'sex', 'age')], 
                                         id = 'id', start = 'start', stop = 'stop', 
                                         covariates = c('anyUseLast14', 'sex', 'age'), 
                                         eventTimes = ev.gen, betas = c(beta.expo, betas.covar))
    dat.gen.gen$event <- dat.gen.gen$Event.NEW

    # Save the number of events in in this generated datasets for verification (must be 285 for all 
    # simulation repetitions)
    n.events.gen[k] <- sum(dat.gen.gen$event)
    
    #------------------------------------------------------------------------------------------------
    # Creation of the two versions of the generated imperfect datasets with:
    #  1) End-point imputation of event times, i.e. each event is imputed at the first clinic visit, 
    #     of the corresponding patient to whom this event was assigned, after the true event time;  
    #  2) Mid-point imputation of event times, i.e. each event time was imputed at the mid-point 
    #     between the visit times before and after the true event time for the corresponding 
    #     individual to whom the event was assigned.
    # [Step 5 (Imperfect dataset generation) of the proposed approach to data-driven simulations]
    #------------------------------------------------------------------------------------------------    

    ## End-point imputation of the generated dataset
    
    # Ids of subjects assigned an event and the event times
    ev.assigned.ids[k, ] <- dat.gen.gen$id[dat.gen.gen$event == 1]
    ev.assigned.times <- dat.gen.gen$stop[dat.gen.gen$event == 1]

    # Apply the function 'prep.datGen.evEnd.fct' to each patient and combine the resulting data 
    # for each patient with the 'rbind' function
    dat.gen.end <- do.call('rbind', lapply(1:N, function(i) 
                            prep.datGen.evEnd.fct(dat.ori.evEnd.subj = dat.ori.evEnd.list[[i]],
                                                  genEv.ids = ev.assigned.ids[k, ], 
                                                  genEv.times = ev.assigned.times, 
                                                  visits.subj = visits[i, -1]))) 
                                                    # 1st column in visits[i, ] is the patient's id,
                                                    # so it's excluded

    # Save the number of events in the resulting dataset (must be 285 for all simulation repetitions) 
    n.events.end[k] <- sum(dat.gen.end$event) 
    
    ## Mid-point imputation of the generated dataset
    
    # Create a list in which each item is the data for one patient in the dataset 'dat.gen.end'
    dat.gen.end.list <- split(dat.gen.end, dat.gen.end$id)
    
    # Apply the function 'impute.evEndToMid.fct' to each patient and combine the resulting data 
    # for each patient with the 'rbind' function
    dat.gen.mid <- do.call("rbind", lapply(1:N, function(i) 
                            impute.evEndToMid.fct(dat.evEnd.subj = dat.gen.end.list[[i]],
                                                  visits.subj = visits[i, -1]))) 
                                                    # 1st value in visits[i, ] is the patient's id
    
    # Save the number of events in the resulting dataset (must be 285 for all simulation repetitions)  
    n.events.mid[k] <- sum(dat.gen.mid$event) 
    
    #---------------------------------------------------------------------------------------------
    # Estimation of Cox PH models on the generated datasets with event times as generated or  
    # imputed at either the end-point or mid-point of intervals between the clinic visit   
    # immediately before and after the generated event times
    # [Step 6 (Analyses) of the proposed approach to data-driven simulations]
    #---------------------------------------------------------------------------------------------
      
    cox.gen[[k]] <- summary(coxph(Surv(start, stop, event) ~ anyUseLast14 + sex + age, 
                                  data = dat.gen.gen))$coef[,c(1,3,5)]
    cox.mid[[k]] <- summary(coxph(Surv(start, stop, event) ~ anyUseLast14 + sex + age, 
                                  data = dat.gen.mid))$coef[,c(1,3,5)]
    cox.end[[k]] <- summary(coxph(Surv(start, stop, event) ~ anyUseLast14 + sex + age, 
                                  data = dat.gen.end))$coef[,c(1,3,5)]
  }
 
  #---------------
  # Saving results
  #--------------- 
  
  # Scenario number 
  if (beta.expo == log(1)){
    sc.no <- 1
  } else if (beta.expo == log(1.5)){
    sc.no <- 2
  } else if (beta.expo == log(2)){
    sc.no <- 3
  } else if (beta.expo == log(2.5)){
    sc.no <- 4
  }

  save(n.reps, beta.expo, betas.covar, n.events.gen, n.events.end, n.events.mid, ev.assigned.ids,
       cox.gen, cox.mid, cox.end,
       file = paste0('./SimResults_sc', sc.no, 'samples.RData'))
}


##################################################################################################
# prep.datGen.evEnd.fct:
#   Function to prepare the generated data for a subject according to his/her assigned status:
#    a) if an event was assigned to the subject, the generated event time is imputed at the first 
#       clinic visit after the generated event time,
#    b) if no event was assigned, keep the full follow-up available in the original dataset for 
#       this subject. 
##################################################################################################

## Arguments:
  # dat.ori.evEnd.subj: data for a subject from the original dataset, with event times imputed at
  #                     first visit after the event time, i.e. at the end of relevant between- 
  #                     visits interval 
  # genEv.ids: vector of identification numbers of all subjects to whom an event was assigned.
  # genEv.times: vector of all generated event times.
  # visits.subj: vector of all clinic visit times of the current subject.

prep.datGen.evEnd.fct <- function(dat.ori.evEnd.subj, genEv.ids, genEv.times, visits.subj){
  
  # If current subject was assigned an event
  if (dat.ori.evEnd.subj$id[1] %in% genEv.ids){
    
    # Identify the assigned event time and the first visit after that event time
    evt.subj <- genEv.times[dat.ori.evEnd.subj$id[1] == genEv.ids]
    vis.after <- min(visits.subj[visits.subj >= evt.subj], na.rm = T)
    
    # Cut follow-up at the first visit after the generated event time
    dat.ori.evEnd.subj <- dat.ori.evEnd.subj[dat.ori.evEnd.subj$stop <= vis.after, ]
    
    # Assign an event at the last line of follow-up for thi subject
    dat.ori.evEnd.subj$event <- c(rep(0, times = (nrow(dat.ori.evEnd.subj) - 1)), 1)
    return(dat.ori.evEnd.subj)
    
  # If current subject was not assigned an event
  } else {
    
    # Assign no event to all lines of follow-up and return the full follow-up
    dat.ori.evEnd.subj$event <- 0
    return(dat.ori.evEnd.subj) 
  }
}


##################################################################################################
# bias.fct, bias.CIsign.fct, relBias.fct, rmse.fct, coverage.fct:                    
#   Functions to calculate respectively:                                      
#     - Bias
#     - Bias and indicate if bias 95% confidence interval (CI) excludes 0 with symbol *
#     - Relative bias
#     - Root mean squared error (RMSE)                                        
#     - Coverage rate of the 95% confidence interval 
#   for log(HR) estimates for a variable for one simulation scenario.                           
##################################################################################################

## Arguments:
  # est: vector of estimates for a variable from all simulation repetitions for a simulation 
  #      scenario.
  # true.beta: true parameter value for that variable for that simulation scenario.
  # SE: vector of standard errors of the estimates from all simulation repetitions for a  
  #     simulation scenario.

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

rmse.fct <- function(est, true.beta) {
  round(sqrt((mean(est) - true.beta) ^ 2 + var(est)), digits = 3)
}

coverage.fct <- function(est, SE, true.beta){
  lowerBound <- est - qnorm(0.975) * SE 
  upperBound <- est + qnorm(0.975) * SE  
  inCI <- ifelse(true.beta > lowerBound & true.beta < upperBound, 1, 0)
  return(round(sum(inCI) / length(est), digits = 3))
}


##################################################################################################
# closer.fct:                    
#   Function to calculate the percentage of simulation repetitions, for a simulation scenario, for  
#   which the estimates for a variable from a model A are closer to the truth parameter value  
#   than estimates from a model B.                           
##################################################################################################

## Arguments:
  # estA: vector of estimates for a variable for model A from all simulation repetitions for a 
  #       simulation scenario.
  # estB: vector of estimates for a variable for model B from all simulation repetitions for a 
  #       simulation scenario.
  # true.beta: true parameter value for that variable for the scenario.

closer.fct <- function(estA, estB, true.beta){
  round(sum(abs(estA - true.beta) < abs(estB - true.beta)) / length(estA) * 100, digits = 1)
}


##################################################################################################
# simulations.ex2.QBA.fct:
#   Function to run the additional quantitative bias analysis (QBA)-like sensitivity analyses of 
#   original data for Example 2.
##################################################################################################

## Arguments:
  # n.reps: number of simulation repetitions per scenario.
  # dat.ori.evEnd: original dataset with each event times imputed at the end of the interval
  #                between the last clinic visit before and the first visit after the true  
  #                (unknown) event time. The dataset must include the following columns (with the 
  #                names):
  #                - 'start' and 'stop' indicating the beginning and end the time interval of the  
  #                  current observation, 
  #                - 'event' indicating the the event status, 
  #                - 'anyUseLast14', 'sex' and 'age' for the time-varying exposure metric (any  
  #                   benzodiazepine in the last 14 days) and the covariates.
  # interval.ev: matrix in which the 1st and 2nd column is respectively the visit time before and 
  #              after each true event time.

simulations.ex2.QBA.fct <- function(n.reps, dat.ori.evEnd, interval.ev){
  
  #--------------------
  # Objects to be saved
  #--------------------
  
  # Save the number of events in the generated datasets 
  n.events.gen <- rep(NA, n.reps)
  
  # Results from the Cox PH models estimated on the generated datasets  
  cox.gen <- list()
  
  ## Loop repeated for each of the 'n.reps' repetitions
  
  cat('Simulation repetition:\n')
  for (k in 1:n.reps){
    cat('  k=', k, '\n')

    #---------------------------------------------------------------------------------------------
    # Generation of datasets 
    #---------------------------------------------------------------------------------------------   
    
    # Randomly assigned each event time over the interval between the visit just before and visit   
    # just after the event time in the original dataset
    ev.gen <- ceiling(runif(nrow(interval.ev), interval.ev[,1], interval.ev[,2]))

    # Cut the follow-up of subjects with an event at the randomly assigned new event times
    dat.gen.gen <- do.call('rbind', by(data = dat.ori.evEnd, dat.ori.evEnd$id, 
                                       cutFUPev.QBA.fct, eventTimes = ev.gen, 
                                       ev.assigned.ids = dat.ori.evEnd$id[dat.ori.evEnd$event==1]))
    
    # Number of events in the resulting dataset (must be 285 for all simulation repetitions)
    n.events.gen[k] <- sum(dat.gen.gen$Event.NEW)
    
    #---------------------------------------------------------------------------------------------
    # Estimation of Cox OH models on the generated datasets with event times as generated 
    #---------------------------------------------------------------------------------------------
    
    cox.gen[[k]] <- summary(coxph(Surv(start, stop, Event.NEW) ~ anyUseLast14 + sex + age, 
                                  data = dat.gen.gen))$coef[,c(1,3,5)]
  }
  
  #---------------
  # Saving results
  #--------------- 
  
  save(n.reps, n.events.gen, cox.gen, 
       file=paste0("./SimResults_QBA.RData"))
}


##################################################################################################
# cutFUPev.QBA.fct:
#   Function to cut the follow-up at the assigned event time in the additional quantitative bias
#   analysis (QBA)-like sensitivity analyses 
##################################################################################################

## Arguments:
  # curSubj: current subject on whom the function is applied
  # ev.assigned.ids: id of all subjects being assigned an event

cutFUPev.QBA.fct <- function(curSubj, eventTimes, ev.assigned.ids){

  # If 'curSubj' was assigned an event
  if (curSubj$id[1] %in% ev.assigned.ids){
    
    curSubj <- curSubj[1:(eventTimes[curSubj$id[1] == ev.assigned.ids]), ]
    curSubj[, 'Event.NEW'] <- c(rep(0, times=(nrow(curSubj) - 1)), 1)
    return(curSubj)
    
  # If 'curSubj' was not assigned an event
  } else {
    curSubj[, 'Event.NEW'] <- 0
    return(curSubj) 
  }
}
