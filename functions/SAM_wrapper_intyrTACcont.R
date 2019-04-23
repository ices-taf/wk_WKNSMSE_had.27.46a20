### ------------------------------------------------------------------------ ###
### stock assessment: SAM wrapper ####
### ------------------------------------------------------------------------ ###
### wrapper for calling SAM
### this makes used of the SAM wrapper in FLfse,
### which in turn calls SAM from the package stockassessment
SAM_wrapper_intyrTACcont <- function(stk, idx, tracking,
                                     genArgs, ### contains ay (assessment year)
                                     forecast = FALSE,
                                     fwd_trgt = "fsq", ### what to target in forecast
                                     fwd_yrs = 1, ### number of years to add
                                     fwd_yrs_average = -3:0, ### years used for averages
                                     fwd_yrs_rec_start = NULL, ### recruitment 
                                     fwd_yrs_sel = -3:-1, ### selectivity
                                     fwd_yrs_lf_remove = -2:-1,
                                     fwd_splitLD = TRUE,
                                     parallel = FALSE,
                                     conf = NULL, ### SAM configuration
                                     par_ini = NULL, ### initial parameters
                                     track_ini = FALSE, ### store ini for next year
                                     ...){
  
  ### get additional arguments
  args <- list(...)
  
  ### get current (assessment) year
  ay <- genArgs$ay
  
  ### check if initial parameter values for SAM exist from last year's fit
  ### and reuse if possible
  ### this overrides generic initial parameters, if they exist in par_ini
  ### (they are only used in first simulation year)
  if (isTRUE(track_ini) & !is.null(attr(tracking@units, "par_ini"))) {
    
    par_ini <- attr(tracking@units, "par_ini")
    
  }
  
  ### fit SAM to provided data
  fit <- FLR_SAM(stk = stk, idx = idx, conf = conf, par_ini = par_ini,
                 DoParallel = parallel, ...)
  
  ### store parameter values and provide them as initial values next year
  ### store in tracking object, this is the only object which does not get 
  ### overriden next year
  ### weird workaround: save as attribute of "unit" slot of tracking,
  ###                   otherwise the attributes will be lost later...
  if (isTRUE(track_ini)) {
    
    attr(tracking@units, "par_ini") <- sam_getpar(fit)
    
  }
  
  ### convert into FLStock
  stk0 <- SAM2FLStock(object = fit, stk = stk) 
  
  ### perform forecast to get SSB ay+1
  if (isTRUE(forecast)) {
    
    ### check how to do forecast
    ### currently, can only do F status quo
    if (!all(fwd_trgt %in% c("fsq","TAC","intyrTACcont"))) {
      stop("only fsq and TAC supported in forecast")
    }
    
    ### years for average values
    ave.years <- range(stk0)[["maxyear"]] + fwd_yrs_average
    ### years for sampling of recruitment years
    rec.years <- seq(from = fwd_yrs_rec_start, to = range(stk0)[["maxyear"]] - 1)
    ### years where selectivity is not used for mean in forecast
    overwriteSelYears <- range(stk0)[["maxyear"]] + fwd_yrs_sel
    
    lst_yr <- range(stk0)[["maxyear"]]
    
    ### extend stk0
    stk0 <- window(stk0, end = range(stk0)[["maxyear"]] + fwd_yrs)
    
    ### modify fwd_yrs in case last data year does not include catch
    if (all(is.na(catch(stk0)[, ac(lst_yr)]))) fwd_yrs <- fwd_yrs + 1
    
    ### forecast years
    yrs <- seq(to = dims(stk0)$maxyear, length.out = fwd_yrs)
    
    ### coerce fit into list if only 1 iteration
    if (is(fit, "sam")) {
      fit <- list(fit)
      class(fit) <- "sam_list"
    }
    
    ### template for forecast targets
    fscale <- ifelse(fwd_trgt == "fsq", 1, NA)
    catchval <- ifelse(fwd_trgt == "TAC", -1, NA)
    # if(all(is.na(catchval))){
    #    catchval <- ifelse(fwd_trgt == "intyrTACcont", -2, NA)
    #  }
    ### recycle target if neccessary
    if (fwd_yrs > length(fwd_trgt)) {
      fscale <- c(fscale, tail(fscale, 1))
      catchval <- c(catchval, tail(catchval, 1))
    }
    ### get recent TAC
    if (genArgs$iy == ay) {
      ### in first year of simulation, use value from OM saved earlier in ay
      TAC_last <- tracking["metric.is", ac(ay)]
    } else {
      ### in following years, use TAC advised the year before
      TAC_last <- tracking["metric.is", ac(ay - 1)]
    }
    
    ### do forecast for all iterations
    fc <- foreach(fit_i = fit, iter_i = seq_along(fit),
                  .errorhandling = "pass") %do% {
                    
                    ### overwrite landing fraction with last year, if requested
                    if (!is.null(fwd_yrs_lf_remove)) {
                      ### index for years to remove/overwrite
                      idx_remove <- nrow(fit_i$data$landFrac) + args$fwd_yrs_lf_remove
                      ### overwrite
                      fit_i$data$landFrac[idx_remove, ] <- fit_i$data$landFrac[rep(nrow(fit_i$data$landFrac), length(idx_remove)), ]
                    }
                    
                    ### define forecast targets for current iteration
                    fscale_i <- fscale
                    ### load TAC as catch target
                    catchval_i <- ifelse(catchval == -1, c(TAC_last[,,,,, iter_i]), catchval)
                    
                    if (fwd_trgt == "intyrTACcont"){
                      # use fsq unless it exceeds TAC
                      # fsq
                      fc_fsq <-stockassessment::forecast(fit = fit_i, 
                                                         fscale = 1, 
                                                         catchval = NA,
                                                         ave.years = ave.years,
                                                         rec.years = rec.years,
                                                         overwriteSelYears = overwriteSelYears,
                                                         splitLD = fwd_splitLD)
                      
                      if(attr(fc_fsq, "tab")[,"catch:median"]>TAC_last[,,,,, iter_i]){
                        # use intermediate year TAC constraint
                        catchval_i<--1
                        fscale_i<-NA
                        if (fwd_yrs > length(fwd_trgt)) {
                          fscale_i <- c(fscale_i, tail(fscale_i, 1))
                          catchval_i <- c(catchval_i, tail(catchval_i, 1))
                        }
                        fscale_i <- fscale_i
                        ### load TAC as catch target
                        catchval_i <- ifelse(catchval_i == -1, c(TAC_last[,,,,, iter_i]), catchval_i)
                      }else{ 
                        #use Fsq
                        catchval_i<-NA
                        fscale_i<-1
                        if (fwd_yrs > length(fwd_trgt)) {
                          fscale_i <- c(fscale_i, tail(fscale_i, 1))
                          catchval_i <- c(catchval_i, tail(catchval_i, 1))
                        }
                        fscale_i <- fscale_i
                        ### load TAC as catch target
                        catchval_i <- ifelse(catchval_i == -1, c(TAC_last[,,,,, iter_i]), catchval_i)
                      }
                    }else{
                      warning("forecast assumption is not intyrTACcont in all forecast years")
                      
                    }
                    
                    ### run forecast
                    fc_i <- stockassessment::forecast(fit = fit_i, 
                                                      fscale = fscale_i, 
                                                      catchval = catchval_i,
                                                      ave.years = ave.years,
                                                      rec.years = rec.years,
                                                      overwriteSelYears = overwriteSelYears,
                                                      splitLD = fwd_splitLD)
                    
                    ### get numbers at age for all forecast years
                    numbers <- lapply(seq(fwd_yrs), function(x) {
                      ### index for numbers at age
                      idx <- seq(length(fit_i$conf$keyLogFsta[1, ]))
                      ### get simulated numbers
                      n <- exp(fc_i[[x]]$sim[, idx])
                      ### median
                      apply(n, 2, median)
                    })
                    numbers <- do.call(cbind, numbers)
                    
                    return(list(stock.n = numbers))
                    
                  }
    ### if forecast failed for a iteration, the list element will for this
    ### iteration will be an error message
    
    ### get numbers
    fc_stock.n <- lapply(fc, "[[", "stock.n")
    
    ### insert stock numbers
    for (iter_i in seq(dim(stk0)[6])) {
      ### do not insert numbers if forecast failed
      if (!is.numeric(fc_stock.n[[iter_i]])) next()
      stock.n(stk0)[, ac(yrs),,,, iter_i] <- fc_stock.n[[iter_i]]
    }
    
    ### extend stock characteristics required for calculation of SSB, 
    ### weights, etc.
    
    ### find years to fill (do not fill years, if there is already data inside)
    yrs_fill <- setdiff(yrs, lst_yr)
    
    stock.wt(stk0)[, ac(yrs_fill)] <- yearMeans(stock.wt(stk0)[, ac(ave.years)])
    m(stk0)[, ac(yrs_fill)] <- yearMeans(m(stk0)[, ac(ave.years)])
    mat(stk0)[, ac(yrs_fill)] <- yearMeans(mat(stk0)[, ac(ave.years)])
    
    harvest.spwn(stk0)[, ac(yrs_fill)] <- yearMeans(harvest.spwn(stk0)[, 
                                                                       ac(ave.years)])
    m.spwn(stk0)[, ac(yrs_fill)] <- yearMeans(m.spwn(stk0)[, ac(ave.years)])
    harvest(stk0)[, ac(yrs_fill)] <- yearMeans(harvest(stk0)[, ac(ave.years)])
    
    ### PLEASE NOTE:
    ### SSB value slightly different from SSB value generated from SAM:
    ### SAM calculates SSB per simulation and gives median
    ### here: calculate median of numbers at age and calculate SSB from
    ###       median numbers
    
  }
  
  ### save convergence for all iterations
  tracking["conv.est", ac(ay)] <- sapply(fit, function(x) {
    if (isTRUE(is(x, "sam"))) {
      return(x$opt$convergence)
    } else {
      return(2)
    }
  })
  ### add perceived F and SSB
  ### done in mp()
  #tracking["F.est", ac(ay)] <- fbar(stk0)[, ac(ay - 1)]
  #tracking["B.est", ac(ay)] <- tail(ssb(stk0))
  
  ### save model fit (list) as attribute in stk0
  attr(stk0, "fit") <- fit
  
  ### return assessed stock, tracking & model fit 
  ### (model fit required later for TAC calculation)
  return(list(stk = stk0, tracking = tracking))
  
}

### ------------------------------------------------------------------------ ###