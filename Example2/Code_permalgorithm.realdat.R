##################################################################################################
##################################################################################################
###
### Documentation and code for the function 'permalgorithm.realdat', which implements the   
### adaptation of the permutational algorithm to input data consisting of time-dependent      
### covariates with values known up to different follow-Up times across subjects
###
### Code by: Marie-Eve Beauchamp                                            
### Last update: March 15, 2024                                              
###
##################################################################################################
##################################################################################################

# NOTE: The version below of the function permalgorithm.realdat was used in the simulations for
  # Example 2 of data-driven simulations (see manuscript for details). For updates of the function
  # 'permalgorithm.realdat', see the version that will be implemented in the next version of the
  # PermAlgo package. The function 'permalgorithm.realdat' is complementary to the function
  # 'permalgorithm' included in the 'PermAlgo' R package, version 1.2. See the documentation below
  # for further details.


###################################################################################################
## DOCUMENTATION FOR THE FUNCTION 'permalgorithm.realdat'
###################################################################################################

#--------------------------------------------------------------------------------------------------
# permalgorithm.realdat  Generates Event Times Conditional On Time-Dependent  Covariates For Which 
#                        Values Are Known Up To Different Follow-Up Times Across Subjects, As In  
#                        Most Real-World Datasets
#--------------------------------------------------------------------------------------------------

### Description ###

  # The function 'permalgorithm.realdat' generates a dataset in which event times are conditional
  # on a user-specified list of covariates, some or all of which could be time-dependent. This 
  # function allows the use of covariate values taken from a real-world dataset with time-dependent
  # covariate values known up to different time points for different subjects. Note that time-
  # dependent covariate values are often unknown after the event or censoring time of a subject in 
  # a real-world data. In contrast, the function 'permalgorithm' in the 'PermAlgo' package   
  # requires that the matrix of covariate values (argument 'Xmat') passed to the function includes   
  # covariate values for each subject up to the maximum length of follow-up (argument 'maxTime').

### Usage ###

  # permalgorithm.realdat(data, id, start, stop, covariates, eventTimes, betas)

### Arguments ###

  # data:   is a data frame in the counting process format where every line represents a time
  #         interval which can either correspond to 1 or several time units during which all
  #         covariate values for a given subject remain constant. The data frame must include 
  #         columns indicating: the subject identifier, the beginning and end of time intervals,     
  #         and the covariates on which the hazard of an event depends. The length of follow-up may  
  #         vary across subjects but all subjects must start their follow-up at the same time (i.e.  
  #         delayed entry is not possible).
  #
  # id:     is a character string representing the name of the column in 'data' identifying the 
  #         subjects, e.g. "id".
  #
  # start:  is a character string representing the name of the column in 'data' identifying the   
  #         beginning of time intervals for the counting process format, e.g. "start".
  #
  # stop:   is a character string representing the name of the column in 'data' identifying the end   
  #         of the time interval for the counting process format, e.g. "stop".
  #
  # covariates: is a vector of character strings representing the names of covariates in 'data' on   
  #         which the hazard of an event depends, e.g. c("dose", "sex", "age"). Covariates cannot be  
  #         factors; instead, users need to code them with binary indicators. 
  #
  # eventTimes: is the vector of all event times to be assigned. Its length corresponds to the 
  #         number of events in the generated dataset and must be smaller or equal to the number of 
  #         subjects in 'data '. Event times must be smaller or equal to the maximum follow-up time 
  #         of subjects.
  #
  # betas:  is a vector of regression coefficients (log hazard) that represent the magnitude of the
  #         relationship between each covariate and the hazard of an event. The length of 'betas' 
  #         must correspond to the length of the vector 'covariates'.

### Details ###

  # The algorithm implemented in the function 'permalgorithm.realdat' is an adaptation of the 
  # permutational algorithm described in Sylvestre and Abrahamowicz (2008) and implemented in the
  # 'permalgorithm' function of the 'PermAlgo' package. This adaption allows using covariate  
  # values from a real-world dataset with time-dependent covariate values known up to different   
  # time points for different subjects. 

  # The current adaptation performs the matching of the v event times specified in 'eventTimes' 
  # with v vectors of covariates values. The matching is based on a permutation probability law
  # derived from the partial likelihood of Cox's Proportional Hazards (PH) model. The number of 
  # events and individual event times are fixed by the user. They can either be i) extracted from a
  # real-word dataset (e.g., the same dataset as for the argument 'data'), or ii) generated from a 
  # user-defined distribution. 

  # Event times in 'eventTimes' are ordered in increasing time order and a assigned sequentially. 
  # Each event time is randomly matched, with weighted sampling, to one vector of covariate values
  # (corresponding to one subject) among subjects available in the corresponding risk set. The risk
  # set for an event time is defined as the set of subjects: a) for whom the covariate values are 
  # known up to this event time, and b) who have not been selected yet for any earlier event

  # If there are no subjects available in the risk set of an event time to be assigned, then the 
  # function stops and returns an error message. This situation may occur if the number of events  
  # is too large with respect to the known follow-up of subjects (person-time in the cohort). 

  # After the v event times from the argument 'eventTimes' are assigned to v subjects, the other     
  # n-v subjects in 'data' are censored at the end of their known follow-up in 'data'. Note that     
  # this differs from the function 'permalgorithm', for which n observed times (event and censoring   
  # times) are randomly assigned to the n subjects in the dataset generated. 

  # Factor variables are not allowed in 'data'. Instead, users need to code them with binary 
  # indicators.

### Value ###

  # The function 'permalgorithm.realdat' returns a data frame including the same columns as in   
  # 'data' and the following two additional columns:

    # Event.NEW: indicator of generated event. 'Event.NEW'=1 when an event occurs and 0 otherwise.

    # Fup.NEW:   individual follow-up time for the corresponding subject in the generated data.

  # The columns indicating the beginning and end of time intervals in the output data frame have  
  # the same names as in 'data' (i.e. arguments 'start' and 'stop'). However, the end of time   
  # intervals are adjusted to event times assigned. All other columns in the 'data' argument 
  # remain unchanged in the output data frame. Note that if the column with event status from 
  # the original dataset is included in 'data', then this column will be included in the 
  # data output. Therefore, to avoid confusion with 'Event.NEW', it may be advisable to exclude 
  # the original event status from 'data'.

  # The number of subjects in the generated dataset is identical to 'data'.

### References ###

  # This algorithm is a variation of the permutational algorithm described and validated in  
  # Sylvestre and Abrahamowicz (2008), and originally presented by Abrahamowicz, MacKenzie and 
  # Esdaile (1996) and MacKenzie and Abrahamowicz (2002).

  # The current version of the permutational algorithm is a flexible tool to generate event times
  # that follow a user-specified or real-world distribution and that are conditional on user-
  # specified covariates. The function 'permalgorithm.realdat ' is especially useful when: 1) at 
  # least one of the covariate is time-dependent so that conventional inversion methods are 
  # difficult to implement, and 2) the time-dependent covariate values are known up to different 
  # time points across subjects, which prevents using the function 'permalgorithm ' from the 
  # 'PermAlgo' package.

  # Please reference the manuscript by Sylvestre and Abrahamowicz (2008), cited below, if this 
  # function is used in any published material.

  # Sylvestre M.-P., Abrahamowicz M. (2008) Comparison of algorithms to generate event times 
  # conditional on time-dependent covariates. Statistics in Medicine 27(14):2618–34.
  # 
  # Abrahamowicz M., MacKenzie T., Esdaile J.M. (1996) Time-dependent hazard ratio: modelling and
  # hypothesis testing with application in lupus nephritis. JASA 91:1432–9.
  # 
  # MacKenzie T., Abrahamowicz M. (2002) Marginal and hazard ratio specific random data generation:
  # Applications to semi-parametric bootstrapping. Statistics and Computing 12(3):245–252.

### Examples ###

  # # The two examples below are based on the original dataset used for Example 2 of data-driven  
  # # simulations presented in the manuscript.
  # 
  # # Set the path to the directory for Example 2 (need to be changed by users)
  # setwd("C:/.../Code/Example_2/")
  # 
  # # Load the original dataset ('data.ori.evEnd')
  # load(file = "./dataOriginal_ex2.RData")
  # 
  # # Description of variables in 'data.ori.evEnd':
  #   # id: patient identifier
  #   # event: event status
  #   # start: start time of the interval
  #   # stop: stop time of the interval
  #   # anyUseLast14: time-dependent binary indicator of any benzodiazepine use in last 14 days 
  #   #               (exposure metric)
  #   # sex: patient's sex
  #   # age: patient's age at baseline
  # 
  # head(data.ori.evEnd)
  # length(unique(data.ori.evEnd$id)) # 1250 subjects
  # sum(data.ori.evEnd$event == 1)    # 285 events
  # 
  # # Maximum follow-up time varies across subjects in 'data.ori.evEnd'
  # summary(by(data.ori.evEnd$stop, data.ori.evEnd$id, max))
  # 
  # ## EXAMPLE 1 - Generating a new dataset using the same event times and covariate values as in  
  # ## the original dataset 'data.ori.evEnd'
  # 
  # # Extract event times from the original dataset
  # true.evt <- data.ori.evEnd$stop[data.ori.evEnd$event == 1]
  # 
  # dat.new1 <- permalgorithm.realdat(data = data.ori.evEnd, id = "id", start = "start", stop = "stop",
  #                                   covariates = c("anyUseLast14", "sex", "age"),
  #                                   eventTimes = true.evt,
  #                                   betas = log(c(1.21, 1.10, 1.07)))
  # 
  # ## EXAMPLE 2 - Generating a new dataset using the covariate values from the original dataset 
  # ## but newly generated event times
  # 
  # # Generate 100 event times during the study follow-up
  # new.evt <- ceiling(runif(n = 100, min = 0, max = max(data.ori.evEnd$stop)))
  # 
  # dat.new2 <- permalgorithm.realdat(data = data.ori.evEnd, id = "id", start = "start", stop = "stop",
  #                                   covariates = c("anyUseLast14", "sex", "age"),
  #                                   eventTimes = new.evt,
  #                                   betas = log(c(1.21, 1.10, 1.07)))


###################################################################################################
## CODE FOR THE MAIN FUNCTION: permalgorithm.realdat
###################################################################################################

permalgorithm.realdat <- function(data, id, start, stop, covariates, eventTimes, betas){
  
  ### Verification that inputs are appropriate ###
  
  if (length(eventTimes) > length(unique(data[, id]))){
    stop("The number of events to assign must be smaller or equal to the number of subjects.")
  }
  if (max(eventTimes) > max(data[, stop])){
    stop("Some values in eventTimes are higher than the maximum follow-up time of subjects.") 
  }
  if (min(eventTimes) <= max(by(data[, start], data[, id], min))){
    stop("Some values in eventTimes are earlier than the beginning of follow-up of some subjects.") 
  }

  eventTimes <- sort(eventTimes)
  ids <- unique(data[, id])
  
  ### Calculate the numerator of sampling probability for subjects in risk set of each event time ###
    # numRS: 1 line per event time (increasing time order), 1 column per subject in 'data'.
    # NA assigned to subjects not in the risk set of an event time.

  if (max(data[, stop] - data[, start]) == 1){
    numRS <- do.call('cbind', by(data, data[, id], 
                                 numRS.1linePerUnit.fct, stop, covariates, eventTimes, betas, 
                                 simplify = FALSE))
  } else {
    numRS <- do.call('cbind', by(data, data[, id], 
                                 numRS.intervals.fct, start, stop, covariates, eventTimes, betas, 
                                 simplify = FALSE))
  }
  colnames(numRS) <- paste0('id.', ids)

  ### Assign event times ###
  
  ev.assigned.ids <- rep(NA, times = length(eventTimes)) # Store id being assigned an event
  
  # For each event time, in increasing time order
  for (i in 1:length(eventTimes)){

    if (sum(!is.na(numRS[i, ])) == 0){
      stop("Events assigned (in increasing time order) up to the ", i-1, "th out of ", 
      length(eventTimes), " events. 
  Cannot assign the remaining events because no subjects are available in their risk sets.")
    }
    
    # Sampling probabilities of being assigned the current event time (see equation (2) in 
    # Sylvestre and Abrahamowicz (2008))
    probPA <- numRS[i, ] / sum(numRS[i, ], na.rm = TRUE) 
    
    # Randomly sample 1 subject available in current risk set 
    if (sum(!is.na(numRS[i, ])) == 1){    
      ev.assigned.ids[i] <- ids[!is.na(probPA)]
    } else {
      ev.assigned.ids[i] <- sample(ids[!is.na(probPA)], size = 1, prob = probPA[!is.na(probPA)])
    }
    
    # Remove the subject selected from risk sets for event times to be subsequently assigned
    if (i < length(eventTimes)){
      numRS[(i+1):length(eventTimes), which(ids == ev.assigned.ids[i])] <- NA
    }
  }

  ### Construct the generated dataset ###
  
  # Notes: Subjects assigned an event finish their follow-up at the assigned event times.
  #        For those not assigned an event, their follow-up ends at same time as in 'data'.
  
  data.GEN <- do.call('rbind', by(data, data[, id], cutFUPev.fct, id, start, stop, eventTimes, 
                                  ev.assigned.ids))
  return(data.GEN)
}


###################################################################################################
## CODE FOR INTERNAL FUNCTIONS 
###################################################################################################

#-----------------------------------------------------------------------------------------------
# Function evaluating if each item in the vector 'p' is included in the interval 'int' delimited  
# by int[1] and int[2] 
#-----------------------------------------------------------------------------------------------

in.interval.fct <- function(int, p){
  p > int[1] & p <= int[2]
}

#------------------------------------------------------------------------------------------------
# Functions to calculate the numerators of sampling probabilities for subjects in the risk set of   
# each event time (NA assigned to subjects not in a risk set)  
#------------------------------------------------------------------------------------------------

### Version of the function for 'data' including time intervals longer than 1 time unit 
  # (function generalized to follow-up starting at other value than 0)

numRS.intervals.fct <- function(curSubj, start, stop, covariates, eventTimes, betas){
  # Arguments:
    # curSubj: current subject on whom the function is applied
  
  # Indicates if (TRUE/FALSE) each event time is included in follow-up of 'curSubj' 
  InRS <- eventTimes <= max(curSubj[, stop])
    # Note: the line above will need to be changed if delayed entry or gaps time intervals are allowed
  
  # Identifies the line in 'curSubj' for each event time during follow-up of 'curSubj' 
  tmp <- t(apply(curSubj[, c(start, stop)],  1, in.interval.fct, p = eventTimes[InRS]))
  tmp2 <- tmp * 1 * 1:nrow(curSubj)
  tmp2[tmp2 == 0] <- NA
  evLineNo.curSubj <- tmp2[!is.na(tmp2)]
  
  # Numerators of sampling probabilities for each event time in 'eventTimes'; NA indicated if 
  # 'curSubj' is not in the risk set for that event time  
  num <- rep(NA, times = length(InRS))
  num[InRS] <- exp(as.matrix(curSubj[evLineNo.curSubj, covariates]) %*%
                     as.matrix(betas, nrow = length(betas)))
    # Note: ties in event times are repeated above
  return(num)
}

### Version of the function for 'data' with the format of 1 line per time unit 
  # (function generalized to follow-up starting at other value than 0)

numRS.1linePerUnit.fct <- function(curSubj, stop, covariates, eventTimes, betas){
  # Arguments:
    # curSubj: current subject on whom the function is applied
  
  # Indicates if (TRUE/FALSE) each event time is included in follow-up of 'curSubj'
  InRS <- eventTimes <= max(curSubj[, stop])
    # Note: the line above will need to be changed if delayed entry or gaps time intervals are allowed
  
  # Numerators of sampling probabilities for each event time in 'eventTimes'; NA indicated if 
  # 'curSubj' is not in the risk set for that event time 
  num <- rep(NA, times = length(InRS))
  num[InRS] <- exp(as.matrix(curSubj[eventTimes[InRS] - (curSubj[1, stop] - 1), covariates]) %*%
                     as.matrix(betas, nrow = length(betas)))
    # Bote: "+ curSubj[1, stop] - 1" accounts for the fact that 1st stop value can be different from 0.
    # Note: ties in event times are repeated above.
  return(num)
}

#----------------------------------------------------------
# Function to cut the follow-up at the assigned event times
#----------------------------------------------------------

cutFUPev.fct <- function(curSubj, id, start, stop, eventTimes, ev.assigned.ids){
  # Arguments:
    # curSubj: current subject on whom the function is applied
    # ev.assigned.ids: id of all subjects being assigned an event
  
  # If 'curSubj' is assigned an event
  if (curSubj[1, id] %in% ev.assigned.ids){
    
    # Identifies the line of 'curSujb' in which falls the assigned event time
    tmp <- t(apply(curSubj[, c(start, stop)],  1, 
                   in.interval.fct, p = eventTimes[curSubj[1, id] == ev.assigned.ids]))
    evLine <- max(tmp * 1 * 1:nrow(curSubj))
    
    curSubj <- curSubj[1:evLine, ]
    curSubj[nrow(curSubj), stop] <- eventTimes[curSubj[1, id] == ev.assigned.ids]
    curSubj[, 'Event.NEW'] <- c(rep(0, times=(evLine - 1)), 1)
    curSubj[, 'Fup.NEW'] <- curSubj[nrow(curSubj), stop] 
    return(curSubj)
    
  # If 'curSubj' is not assigned an event
  } else {
    curSubj[, 'Event.NEW'] <- 0
    curSubj[, 'Fup.NEW'] <- curSubj[nrow(curSubj), stop]
    return(curSubj) 
  }
}
