SAM_uncertainty_had<-function(fit, n = 1000, print_screen = FALSE) {
  
  if (!is(fit, "sam")) stop("fit has to be class \"sam\"")
  
  ### index for fishing mortality ages
  idxF <- fit$conf$keyLogFsta[1, ] + dim(stk)[1] + 1
  idxF <- idxF[idxF != 0] ### remove 0s
  
  ### index for F variances (usually some ages are bound)
  idxNVar <- fit$conf$keyVarLogN 
  
  ### get ages used for calculating fbar
  bAges <- fit$conf$fbarRange
  bAges <- do.call(':', as.list(bAges))
  
  ### index for stock numbers ages
  #idxN <- 1:ncol(natural.mortality)
  idxN <- seq(min(fit$data$minAgePerFleet), max(fit$data$maxAgePerFleet))
  ### index for observation variances
  idxObs <- fit$conf$keyVarObs # starts at 0
  
  ##Resample estimated values to get N, F and q 
  
  ### calculate standard deviations of model parameters
  . <- capture.output(sds <- TMB::sdreport(obj = fit$obj, 
                                           par.fixed = fit$opt$par,
                                           getJointPrecision = TRUE))
  if (isTRUE(print_screen)) cat(paste0(., sep = "\n"))
  
  ### extract values for parameters
  est <- c(sds$par.fixed, sds$par.random)
  ### estimate covariance?
  cov <- solve(sds$jointPrecision)
  
  ### create random values based on estimation and covariance
  sim.states <- mvrnorm(n, est, cov)
  ### matrix, columns are values, rows are requested samples
  table(colnames(sim.states))
  ### contains, among others, logF, logN...
  
  ### combine SAM estimate and random samples
  #dat <- rbind(est, sim.states)
  dat <- sim.states
  
  ### ---------------------------------------------------------------------- ###
  ### stock ####
  ### ---------------------------------------------------------------------- ###
  
  ### stock characteristics
  min_age <- min(fit$data$minAgePerFleet[fit$data$fleetTypes == 0])
  max_age <- max(fit$data$maxAgePerFleet[fit$data$fleetTypes == 0])
  years <- fit$data$years
  ### FLQuant template for stock
  stk_template1 <- FLQuant(dimnames = list(age = min_age:max_age, year = years,
                                          iter = 1:n))
  stk_template2 <- FLQuant(dimnames = list(age = min_age:(max_age-1), year = years,
                                           iter = 1:n))
  ### numbers at age
  stock.n <- stk_template1
  stock.n[] <- exp(t(dat[, colnames(dat) == "logN"]))
  
  ### F at age
  harvest <- stk_template2
  harvest[] <- exp(t(dat[, colnames(dat) == "logF"]))
  
  #age 7 and age 8 are the same
  tmp<-stk_template1
  tmp[ac(0:7),]<-harvest
  tmp[ac(8),]<-harvest[ac(7),]
  harvest<-tmp
  
  ### ---------------------------------------------------------------------- ###
  ### surveys ####
  ### ---------------------------------------------------------------------- ###
  
  ### survey specs
  idx_surveys <- which(fit$data$fleetTypes > 0) ### which observation are surveys
  ### age range of surveys
  survey_ages <- lapply(seq_along(idx_surveys), function(x) {
    seq(fit$data$minAgePerFleet[idx_surveys][x],
        fit$data$maxAgePerFleet[idx_surveys][x])
  })
  ### index for estimated parameters
  sum(colnames(dat) == "logFpar") ### there are 9 parameters
  ### I assume: 5 for Q1, 4 for Q3, as they have this many ages...
  #survey_ages_idx <- split(seq(length(unlist(survey_ages))), 
  #                         rep(seq(survey_ages), sapply(survey_ages, length))) 
  survey_ages_idx<-list(c(1,2,3,4,4),c(5,6,7,8,9,9))
    
  ### get catchability at age (time-invariant) samples
  catchability <- lapply(seq_along(idx_surveys), function(x) {
    
    ### create FLQuant template
    tmp <- FLQuant(dimnames = list(age = survey_ages[[x]],
                                   year = "all", iter = 1:n))
    ### fill with catchability values
    tmp[] <- exp(t(dat[, colnames(dat) == "logFpar"][,survey_ages_idx[[x]]]))
    
    return(tmp)
    
  })
  
  ### ---------------------------------------------------------------------- ###
  ### standard deviation - catch ####
  ### ---------------------------------------------------------------------- ###
  ### time-invariant
  
  ### template
  catch_sd <- FLQuant(dimnames = list(age = dimnames(stock.n)$age, year = "all",
                                      iter = 1:n))
  
  ### index for catch sd (some ages are linked)
  catch_sd_idx <- idxObs[1, ][idxObs[1, ] > -1] + 1
  
  ### extract values
  catch_sd[] <- exp(t(dat[, colnames(dat) == "logSdLogObs"][, catch_sd_idx]))
  
  ### ---------------------------------------------------------------------- ###
  ### standard deviation - surveys ####
  ### ---------------------------------------------------------------------- ###
  ### time-invariant
  
  ### get catchability at age (time-invariant) samples
  survey_sd <- lapply(seq_along(idx_surveys), function(x) {
    
    ### create FLQuant template
    tmp <- FLQuant(dimnames = list(age = survey_ages[[x]],
                                   year = "all", iter = 1:n))
    ### index for sd (some ages are linked)
    idx <- idxObs[idx_surveys[x], ]
    idx <- idx[idx > -1] + 1
    
    ### fill with catchability values
    tmp[] <- exp(t(dat[, colnames(dat) == "logSdLogObs"][, idx]))
    
    return(tmp)
    
  })
  
  return(list(stock.n = stock.n, harvest = harvest,
              survey_catchability = catchability, catch_sd = catch_sd,
              survey_sd = survey_sd))
  
}
