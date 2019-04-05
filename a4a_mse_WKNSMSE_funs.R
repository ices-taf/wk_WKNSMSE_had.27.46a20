### ------------------------------------------------------------------------ ###
### functions for running a4a standard MSE at WKNSMSE ####
### ------------------------------------------------------------------------ ###



### ------------------------------------------------------------------------ ###
### calculate survey index/indices from stock ####
### ------------------------------------------------------------------------ ###
#' calculate survey index/indices from FLStock
#'
#' This function calculates survey indices from the numbers at age of an 
#' FLStock object
#'
#' @param stk Object of class \linkS4class{FLStock} with stock and fishery data.
#' @param idx Object of class \linkS4class{FLIndices} or \linkS4class{FLIndex}.
#' @param use_q If \code{TRUE} index numbers are calculated by multiplying
#'   numbers at age in the stock with the catchability.
#' @param use_time If \code{TRUE} observed numbers in the survey are corrected
#'   for fishing and natural mortality.
#'
#' @return An object of class \code{FLIndex} or \code{FLIndex} with the survey
#'   index stored in the \code{index} slot.
#'
#' @export

setGeneric("calc_survey", function(stk, idx, use_q = TRUE, use_time = TRUE) {
  standardGeneric("calc_survey")
})

### stk = FLStock, idx = FLIndices
#' @rdname calc_survey
setMethod(f = "calc_survey",
          signature = signature(stk = "FLStock", idx = "FLIndices"),
          definition = function(stk, idx, use_q = TRUE, use_time = TRUE) {
            
            ### apply function to every element of FLIndices
            lapply(X = idx, FUN = calc_survey_ind, stk = stk, use_q = use_q, 
                   use_time = use_time)
            
          })
### stk = FLStock, idx = FLIndex
#' @rdname calc_survey
setMethod(f = "calc_survey",
          signature = signature(stk = "FLStock", idx = "FLIndex"),
          definition = function(stk, idx, use_q = TRUE, use_time = TRUE) {
            
            calc_survey_ind(stk = stk, idx = idx, use_q = use_q, use_time = use_time)
            
          })



calc_survey_ind <- function(stk, idx, 
                            use_q = TRUE, ### catchability
                            use_time = TRUE ### timing of survey
) {
  
  ### find ranges for years, ages & iters
  ages <- intersect(dimnames(index(idx))$age, dimnames(stock.n(stk))$age)
  years <- intersect(dimnames(index(idx))$year, dimnames(stock.n(stk))$year)
  iter <- intersect(dimnames(index(idx))$iter, dimnames(stock.n(stk))$iter)
  
  ### timing of survey
  if (isTRUE(use_time)) {
    ### use mean of fishing period
    time <- mean(range(idx)[c("startf", "endf")])
  } else {
    ### otherwise assume beginning of year
    time <- 0
  }
  
  ### extract stock numbers for requested/available dimensions
  index.n <- stock.n(stk)[ac(ages), ac(years),,,, ac(iter)]
  ### get Z = M & F
  Z <- m(stk)[ac(ages), ac(years),,,, ac(iter)] +
    harvest(stk)[ac(ages), ac(years),,,, ac(iter)]
  
  ### estimate stock numbers at time of survey
  index.n <- index.n * exp(-time * Z)
  
  ### add catchability, if requested
  if (isTRUE(use_q)) {
    index.n <- index.n * index.q(idx)[ac(ages), ac(years),,,, ac(iter)]
  }
  
  ### insert values into index
  index(idx)[ac(ages), ac(years),,,, ac(iter)] <- index.n
  
  return(idx)
  
}

### ------------------------------------------------------------------------ ###
### observations wrapper ####
### ------------------------------------------------------------------------ ###
oem_WKNSMSE <- function(stk, 
                        deviances, 
                        observations, 
                        genArgs, 
                        tracking, 
                        catch_timing = -1, ### catch timing relative to ay
                        idx_timing = -1, ### index timing relative to ay
                        use_catch_residuals = FALSE, ### use residuals for
                        use_idx_residuals = FALSE,   ### observations
                        use_stk_oem = FALSE, ### biological parameters, wts etc
                        dd_M = NULL,
                        ...) {
  #browser()
  ### current (assessment) year
  ay <- genArgs$ay
  
  
  ### Density-dependent M
  # Calculate 3-year means of M from the OM on key-run years
  # to simulate the SAM process
  if (!is.null(dd_M) & (ay %% 3 == 1)) {
    m(observations$stk)[, ac(ay:(ay+2))] <- yearMeans(m(stk)[, ac((ay-3):(ay-1))])
  }
  
  ### create object for observed stock
  if (!isTRUE(use_stk_oem)) {
    
    ### use OM as template
    stk0 <- stk
    
  } else {
    ### otherwise use observed stock provided in observations object
    ### this can include different biological data, e.g. weights at age
    stk0 <- observations$stk
    ### update fishery data
    catch.n(stk0) <- catch.n(stk)
    catch(stk0) <- catch(stk)
    discards.n(stk0) <- discards.n(stk)
    discards(stk0) <- discards(stk)
    landings.n(stk0) <- landings.n(stk)
    landings(stk0) <- landings(stk)
    
  }
  
  ### add uncertainty to catch
  if (isTRUE(use_catch_residuals)) {
    
    ### implement for catch at age
    catch.n(stk0) <- catch.n(stk) * deviances$stk$catch.dev
    
    ### split catch into discards and landings, based on landing fraction
    landings.n(stk0) <- catch.n(stk0) * (landings.n(stk) / catch.n(stk))
    discards.n(stk0) <- catch.n(stk0) * (1 - landings.n(stk) / catch.n(stk))
    
    ### update total catch/discards/landings
    catch(stk0) <- computeCatch(stk0)
    landings(stk0) <- computeLandings(stk0)
    discards(stk0) <- computeDiscards(stk0)
    
  }
  
  ### cut off years
  ### workaround for NS cod: survey until intermediate year, but catch stops
  ### 1 year earlier
  ### slots such as natural mortality, maturity need to be kept, otherwise
  ### SAM will fall over
  if (any(idx_timing > catch_timing)) {
    
    ### keep stock until last survey data year
    stk0 <- window(stk0, end = ay + max(idx_timing))
    ### find years to remove
    yrs_remove <- (ay + catch_timing + 1):ay
    ### remove catch data
    catch(stk0)[, ac(yrs_remove)] <- NA
    catch.n(stk0)[, ac(yrs_remove)] <- NA
    catch.wt(stk0)[, ac(yrs_remove)] <- NA
    landings(stk0)[, ac(yrs_remove)] <- NA
    landings.n(stk0)[, ac(yrs_remove)] <- NA
    landings.wt(stk0)[, ac(yrs_remove)] <- NA
    discards(stk0)[, ac(yrs_remove)] <- NA
    discards.n(stk0)[, ac(yrs_remove)] <- NA
    discards.wt(stk0)[, ac(yrs_remove)] <- NA
    
  } else {
    
    stk0 <- window(stk0, end = ay + catch_timing)
    
  }
  
  ### calculate index values
  observations$idx <- calc_survey(stk = stk, idx = observations$idx)
  
  ### observed survey
  idx0 <- observations$idx 
  
  ### add uncertainty to observed indices
  if (isTRUE(use_idx_residuals)) {
    
    idx0 <- lapply(seq_along(idx0), function(idx_i) {
      idx_tmp <- idx0[[idx_i]]
      index(idx_tmp) <- index(idx_tmp) * deviances$idx[[idx_i]]
      return(idx_tmp)
    })
    
  }
  
  ### timing of survey
  ### 0: intermediate/assessment year
  ### <0: fewer years available & vice versa
  if (length(idx_timing) < length(idx0)) { 
    idx_timing <- rep(idx_timing, length(idx0))
  }
  ### restrict years for indices based on timing
  idx0 <- lapply(seq_along(idx0), function(x) {
    window(idx0[[x]], end = ay + idx_timing[x])
  })
  idx0 <- FLIndices(idx0) ### restore class
  names(idx0) <- names(observations$idx)
  
  ### return observations
  return(list(stk = stk0, idx = idx0, observations = observations, 
              tracking = tracking))
  
}


### ------------------------------------------------------------------------ ###
### stock assessment: SAM wrapper ####
### ------------------------------------------------------------------------ ###
### wrapper for calling SAM
### this makes used of the SAM wrapper in FLfse,
### which in turn calls SAM from the package stockassessment
SAM_wrapper <- function(stk, idx, tracking,
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
    catchval <- ifelse(fwd_trgt == "intyrTACcont", -2, NA)
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
                    if (catchval == -2){
                      # use fsq unless it exceeds TAC
                      fscale<-1
                      catchval<-NA
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
                        catchval<--1
                        fscale<-NA
                      }
                    }
                    
                    fscale_i <- fscale
                    ### load TAC as catch target
                    catchval_i <- ifelse(catchval == -1, c(TAC_last[,,,,, iter_i]), catchval)
                    
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
### phcr: parameterize HCR ####
### ------------------------------------------------------------------------ ###
phcr_WKNSMSE <- function(Btrigger = NULL, Ftrgt = NULL, Bpa = NULL, Fpa = NULL,
                         Blim = NULL, tracking, ...) {
  
  ### coerce existing values into FLPar
  hcrpars <- list("Btrigger" = Btrigger, "Ftrgt" = Ftrgt, "Bpa" = Bpa, 
                  "Fpa" = Fpa, "Blim" = Blim)
  hcrpars <- hcrpars[!sapply(hcrpars, is.null)]
  hcrpars <- do.call(FLPar, hcrpars)
  
  ### if more iterations provided than neccessary, subset
  if (dims(hcrpars)$iter > dims(tracking)$iter) {
    hcrpars <- hcrpars[, dimnames(tracking)$iter]
  }
  
  ### return as list
  ### keep tracking unchanged
  return(list(hcrpars = hcrpars, tracking = tracking))
  
}

### ------------------------------------------------------------------------ ###
### hcr: ICES advice rule and modificatinons ####
### ------------------------------------------------------------------------ ###

### target Ftrgt
### A: if SSB < Btrigger, target reduced: Ftrgt * SSB/Btrigger
### B & C: different behaviour if SSB < Blim
hcr_WKNSME <- function(stk, genArgs, hcrpars, tracking, 
                       option = "A", ### WKNSMSE options
                       ...) {
  
  ### target F
  ### reduce F if SSB below Btrigger
  
  ### get current (assessment) year
  ay <- genArgs$ay
  ### number of iterations
  it <- dim(stk)[6]
  
  ### get reference values and extend iteration dimension
  Ftrgt <- propagate(FLPar(hcrpars["Ftrgt"]), it)
  Btrigger <- propagate(FLPar(hcrpars["Btrigger"]), it)
  
  ### get Blim
  ### if Blim not specified, set to zero
  if ("Blim" %in% dimnames(hcrpars)$params) {
    
    Blim <- propagate(FLPar(hcrpars["Blim"]), it)
    
  } else {
    
    Blim <- propagate(FLPar(0), it)
    
  } 
  
  ### SSB status (last year only) relative to Btrigger
  status_Btrigger <- tail(ssb(stk), 1) / Btrigger 
  ### SSB status (last year only) relative to Blim
  status_Blim <- tail(ssb(stk), 1) / Blim
  
  ### positions (iterations) where SSB is below Btrigger
  pos_Btrigger <- which(status_Btrigger < 1)
  ### below Blim
  pos_Blim <- which(status_Blim < 1)
  
  ### default ICES HCR (option A):
  ### if SSB<Btrigger => F = Ftrgt * SSB/Btrigger
  if (option == "A") {
    
    mult <- ifelse(status_Btrigger < 1, status_Btrigger, 1)
    
  } else if (option == "B") {
    ### option B:
    ### if SSB<Btrigger => F = Ftrgt * SSB/Btrigger
    ### if SSB<Blim => F = Ftrgt * 0.25
    ### SSB < Btrigger
    mult <- ifelse(status_Btrigger < 1, status_Btrigger, 1)
    ### set to 0.25 if SSB < Blim
    if (length(pos_Blim) > 0) {
      mult[,,,,, pos_Blim] <- 0.25
    }
    
  } else if (option == "C") {
    ### option C:
    ### if SSB<Btrigger => F = Ftrgt * SSB/Btrigger
    ### if SSB<Blim => F = max(Ftrgt * 0.25, Ftrgt * SSB/Btrigger)
    ### SSB < Btrigger
    mult <- ifelse(status_Btrigger < 1, status_Btrigger, 1)
    ### if SSB < Blim, use maximum of SSB/Btrigger or 0.25
    ### i.e. limit ratio to 0.25
    if (length(pos_Blim) > 0) {
      mult[,,,,, pos_Blim] <- c(ifelse(mult[,,,,, pos_Blim] < 0.25, 0.25, 
                                       mult[,,,,, pos_Blim]))
    }
    
  } else {
    
    ### unknown options
    mult <- 1
    
  }
  
  ### new target
  Ftrgt <- Ftrgt * mult
  
  ### create ctrl object
  ctrl <- getCtrl(values = Ftrgt, quantity = "f", years = ay + 1, it = it)
  
  ### save in tracking
  ### done later
  # tracking["advice", ac(ay)] <- ctrl@trgtArray[ac(ay + 1), "val", ]
  
  return(list(ctrl = ctrl, tracking = tracking))
  
}



### ------------------------------------------------------------------------ ###
### management implementation: short term forecast with SAM ####
### ------------------------------------------------------------------------ ###
### short term forecast with SAM
### including TAC constraint
is_WKNSMSE <- function(stk, tracking, ctrl,
                       genArgs, ### contains ay (assessment year)
                       TAC_constraint = c(FALSE, TRUE),
                       upper = Inf, lower = -Inf, Btrigger_cond = FALSE,
                       ### short term forecast
                       fwd_trgt = c("fsq", "hcr"), ### target in forecast
                       fwd_yrs = 2, ### number of years to add
                       fwd_yrs_average = -3:0, ### years used for averages
                       fwd_yrs_rec_start = NULL, ### recruitment 
                       fwd_yrs_sel = -3:-1, ### selectivity
                       fwd_yrs_lf_remove = -2:-1,
                       fwd_splitLD = TRUE, 
                       ### banking and borrowing
                       BB = FALSE, ### banking and borrowing
                       ### check stock status before applying BB
                       BB_check_hcr = FALSE, ### check status before forecast
                       BB_check_fc = FALSE, ### check status after forecast
                       BB_rho, ### definition of BB
                       ### reference points
                       hcrpars = list(),
                       ...) {
  
  ### get current (assessment) year
  ay <- genArgs$ay
  ### number of iterations
  it <- dim(stk)[6]
  
  ### retrieve SAM model fit (list)
  fit <- attr(stk, "fit")
  
  ### check class of model fit(s)
  if (!class(fit) %in% c("sam", "sam_list")) 
    stop("attr(stk0, \"fit\") has to be class sam or sam_list")
  
  ### if single fit, turn into list
  if (is(fit, "sam")) fit <- list(fit)
  
  ### if conditional banking & borrowing applied, extend forecast for one more
  ### year to check SSB in year after advice year
  ### for this additional forecast assume Fsq as target
  ### i.e. target F from HCR twice, in analogy to intermediate year assumption
  if (isTRUE(BB) & isTRUE(BB_check_fc)) {
    
    ### duplicate last target value
    fwd_trgt <- c(fwd_trgt, tail(fwd_trgt, 1))
    
  }
  
  ### get recent TAC
  if (genArgs$iy == ay) {
    ### in first year of simulation, use value from OM saved earlier in ay
    TAC_last <- tracking["metric.is", ac(ay)]
  } else {
    ### in following years, use TAC advised the year before
    TAC_last <- tracking["metric.is", ac(ay - 1)]
  }
  
  ### go through all model fits
  fc <- foreach(fit_i = fit, iter_i = seq_along(fit), 
                .errorhandling = "pass") %do% {
                  
                  ### overwrite landing fraction with last year, if requested
                  if (!is.null(fwd_yrs_lf_remove)) {
                    ### index for years to remove/overwrite
                    idx_remove <- nrow(fit_i$data$landFrac) + fwd_yrs_lf_remove
                    ### overwrite
                    fit_i$data$landFrac[idx_remove, ] <- fit_i$data$landFrac[rep(nrow(fit_i$data$landFrac), length(idx_remove)), ]
                  }
                  
                  ### check how to do forecast
                  ### can handle F status quo, F target from ctrl object and TAC
                  ### scaled F
                  fscale <- ifelse(fwd_trgt == "fsq", 1, NA)
                  ### target F values
                  fval <- ifelse(fwd_trgt == "hcr", ctrl@trgtArray[, "val", iter_i], NA)
                  ### target catch values
                  catchval <- ifelse(fwd_trgt == "TAC", c(TAC_last[,,,,, iter_i]), NA)
                  
                  ### years for average values
                  ave.years <- max(fit_i$data$years) + fwd_yrs_average
                  ### years for sampling of recruitment years
                  if (is.null(fwd_yrs_rec_start)) {
                    rec.years <- fit_i$data$years ### use all years, if not defined
                  } else {
                    rec.years <- seq(from = fwd_yrs_rec_start, max(fit_i$data$years))
                  }
                  
                  ### years where selectivity is not used for mean in forecast
                  overwriteSelYears <- max(fit_i$data$years) + fwd_yrs_sel
                  
                  ### forecast 
                  fc_i <- stockassessment::forecast(fit = fit_i, 
                                                    fscale = fscale, 
                                                    fval = fval, 
                                                    catchval = catchval,
                                                    ave.years = ave.years,
                                                    rec.years = rec.years,
                                                    overwriteSelYears = overwriteSelYears,
                                                    splitLD = fwd_splitLD)
                  
                  ### return forecast table
                  return(attr(fc_i, "tab"))
                  
                }
  ### if forecast fails, error message returned
  ### replace error message with NA
  ### extract catch target
  catch_target <- sapply(fc, function(x) {
    if (is(x, "error")) {
      return(NA)
    } else {
      return(x[ac(ay + 1), "catch:median"])
    }
  })
  
  ### get reference points
  hcrpars <- hcrpars[!sapply(hcrpars, is.null)]
  hcrpars <- do.call(FLPar, hcrpars)
  ### if more iterations provided than neccessary, subset
  if (dims(hcrpars)$iter > dims(tracking)$iter) {
    hcrpars <- hcrpars[, dimnames(tracking)$iter]
  } else if (isTRUE(dim(hcrpars)[2] < it)) {
    hcrpars <- propagate(hcrpars, it)
  }
  
  ### ---------------------------------------------------------------------- ###
  ### TAC constraint ####
  if (isTRUE(TAC_constraint)) {
    
    ### target year
    yr_target <- ctrl@target$year
    
    ### get previous target from tracking
    catch_prev <- tracking["metric.is", ac(yr_target - 1), drop = TRUE]
    
    ### change in advice, % of last advice
    change <- (catch_target / catch_prev) * 100
    
    ### limit changes
    changes_new <- change
    ### upper limit
    changes_new <- ifelse(changes_new > upper, upper, changes_new) 
    ### lower limit
    changes_new <- ifelse(changes_new < lower, lower, changes_new) 
    
    ### find positions which exceed limits
    pos <- which(changes_new != change)
    
    ### conditional constraint based on SSB>=Btrigger?
    if (isTRUE(Btrigger_cond)) {
      
      ### get SSB in TAC year from forecast
      SSB_TACyr <- sapply(fc, function(x) {
        if (is(x, "error")) {
          return(NA)
        } else {
          return(x[ac(ay + 1), "ssb:median"])
        }
      })
      
      ### iterations where SSB is at or above Btrigger at start of TAC year
      pos_Btrigger <- which(SSB_TACyr >= c(hcrpars["Btrigger"]))
      ### only apply TAC constraint if both
      ### - TAC change exceeds limit
      ### - stock at or above Btrigger
      pos <- intersect(pos, pos_Btrigger)
      
    }
    
    ### modify advice
    catch_target[pos] <- catch_prev[pos] * changes_new[pos]/100
    
  }
  
  ### ---------------------------------------------------------------------- ###
  ### banking and borrowing ####
  
  if (isTRUE(BB)) {
    
    ### get current rho
    BB_rho_i <- tail(rep(BB_rho, length.out = (ay - genArgs$y0)), 1)
    
    ### get catch borrowed last year
    BB_return <- tracking["BB_borrow", ac(ay - 1)]
    ### assume nothing borrowed if NA
    BB_return <- ifelse(!is.na(BB_return), BB_return, 0)
    
    ### get catch banked last year
    BB_bank_use <- tracking["BB_bank", ac(ay - 1)]
    BB_bank_use <- ifelse(!is.na(BB_bank_use), BB_bank_use, 0)
    
    ### bank for next year
    if (BB_rho_i < 0) {
      BB_bank <- catch_target * abs(BB_rho_i)
    } else {
      BB_bank <- rep(0, it)
    }
    
    ### borrow from next year
    if (BB_rho_i > 0) {
      BB_borrow <- catch_target * abs(BB_rho_i)
    } else {
      BB_borrow <- rep(0, it)
    }
    
    ### conditional banking and borrowing?
    ### first condition: for HCR option A (stability option D2)
    ### apply BB only if HCR option A1 (not A2) is applied
    ### i.e. stop BB if SSB is below Btrigger in TAC year
    
    ### find iterations where SSB is below Btriggerin TAC year
    if (isTRUE(BB_check_hcr)) {
      pos_hcr <- which(c(tail(ssb(stk), 1)) < c(hcrpars["Btrigger"]))
    } else {
      pos_hcr <- integer(0)
    }
    
    ### second condition: for HCR options B & C (stability option E2)
    ### stop BB if either
    ### - if SSB is below Bpa AND F above Fpa in TAC year
    ### - if SSB is below Bpa in TAC year and year after
    ### if TAC restricted by TAC constraint, additional forecast required
    ### to estimate stock status when only TAC is fished
    pos_fc <- integer(0)
    if (any(c(BB_bank, BB_borrow) > 0) & isTRUE(BB_check_fc)) {
      
      ### if TAC constraint activated, do a forecast
      ### check if TAC constraint used
      if (isTRUE(TAC_constraint)) {
        ### check where TAC constraint is implemented
        pos_constr <- which(changes_new != change)
        if (isTRUE(length(pos_constr) > 0)) {
          
          ### go through model fits of requested iterations
          fc_new <- foreach(fit_i = fit[pos_constr], 
                            iter_i = seq_along(fit)[pos_constr],
                            .errorhandling = "pass") %do% {
                              
                              ### overwrite landing fraction with last year, if requested
                              if (!is.null(fwd_yrs_lf_remove)) {
                                ### index for years to remove/overwrite
                                idx_remove <- nrow(fit_i$data$landFrac) + fwd_yrs_lf_remove
                                ### overwrite
                                fit_i$data$landFrac[idx_remove, ] <- fit_i$data$landFrac[rep(nrow(fit_i$data$landFrac), length(idx_remove)), ]
                              }
                              
                              ### check how to do forecast
                              ### target recently advised TAC
                              catchval <- ifelse(fwd_trgt == "TAC", c(TAC_last[,,,,, iter_i]), NA)
                              ### scaled F
                              fscale <- ifelse(fwd_trgt == "fsq", 1, NA)
                              ### target new TAC after implementation of TAC constraint TAC year
                              catchval[which(fwd_trgt == "hcr")[1]] <- catch_target[iter_i]
                              ### target fsq (F that corresponds to TAC) in following years
                              fscale[which(fwd_trgt == "hcr")[-1]] <- 1
                              
                              ### years for average values
                              ave.years <- max(fit_i$data$years) + fwd_yrs_average
                              ### years for sampling of recruitment years
                              if (is.null(fwd_yrs_rec_start)) {
                                rec.years <- fit_i$data$years ### use all years, if not defined
                              } else {
                                rec.years <- seq(from = fwd_yrs_rec_start, max(fit_i$data$years))
                              }
                              
                              ### years where selectivity is not used for mean in forecast
                              overwriteSelYears <- max(fit_i$data$years) + fwd_yrs_sel
                              
                              ### forecast 
                              fc_i <- stockassessment::forecast(fit = fit_i, 
                                                                fscale = fscale, 
                                                                catchval = catchval,
                                                                ave.years = ave.years,
                                                                rec.years = rec.years,
                                                                overwriteSelYears = overwriteSelYears,
                                                                splitLD = fwd_splitLD)
                              
                              ### return forecast table
                              return(attr(fc_i, "tab"))
                              
                            }
          
          ### overwrite updated forecasts
          fc[pos_constr] <- fc_new
          
          ### get SSB in TAC year
          SSB_TACyr <- sapply(fc, function(x) {
            if (is(x, "error")) {
              return(NA)
            } else {
              return(x[ac(ay + 1), "ssb:median"])
            }
          })
          ### get SSB in year after TAC year
          SSB_TACyr1 <- sapply(fc, function(x) {
            if (is(x, "error")) {
              return(NA)
            } else {
              return(x[ac(ay + 2), "ssb:median"])
            }
          })
          ### get F in TAC year
          F_TACyr <- sapply(fc, function(x) {
            if (is(x, "error")) {
              return(NA)
            } else {
              return(x[ac(ay + 1), "fbar:median"])
            }
          })
          
          ### check if SSB in TAC year below Bpa AND F above Fpa
          pos_fc1 <- which(SSB_TACyr < c(hcrpars["Bpa"]) & 
                             F_TACyr > c(hcrpars["Fpa"]))
          ### check if SSB below Bpa in TAC year and year after
          pos_fc2 <- which(SSB_TACyr < c(hcrpars["Bpa"]) & 
                             SSB_TACyr1 < c(hcrpars["Bpa"]))
          ### combine both conditions
          pos_fc <- union(pos_fc1, pos_fc2)
          
        }
        
      }
      
    }
    
    ### get position where BB is stopped
    pos_stop <- union(pos_hcr, pos_fc)
    ### stop banking and borrowing
    ### (but still pay back/use from last year)
    BB_bank[pos_stop] <- 0
    BB_borrow[pos_stop] <- 0
    
    ### correct target catch later in iem module
    # catch_target <- catch_target - c(BB_return) +
    #   c(BB_bank_use) - c(BB_bank) + c(BB_borrow)
    
  } else {
    ### if B&B not applied, store 0 
    BB_return <- BB_bank_use <- BB_bank <- BB_borrow <- 0
    
  }
  
  ### save B&B transfers in  tracking
  tracking["BB_return", ac(ay)] <- BB_return
  tracking["BB_bank_use", ac(ay)] <- BB_bank_use
  tracking["BB_bank", ac(ay)] <- BB_bank
  tracking["BB_borrow", ac(ay)] <- BB_borrow
  
  ### create ctrl object
  ctrl <- getCtrl(values = catch_target, quantity = "catch", 
                  years = ctrl@target$year, it = it)
  
  ### return catch target and tracking
  return(list(ctrl = ctrl, tracking = tracking))
  
}

### ------------------------------------------------------------------------ ###
### implementation error: banking and borrowing ####
### ------------------------------------------------------------------------ ###
### so far only banking and borrowing (B&B) implemented

iem_WKNSMSE <- function(tracking, ctrl,
                        genArgs, ### contains ay (assessment year)
                        BB = FALSE, ### apply banking and borrowing
                        ...) {
  
  if (isTRUE(BB)) {
    
    ### get current (assessment) year
    ay <- genArgs$ay
    
    ### retrieve banking and borrowing values
    BB_return <- tracking["BB_return", ac(ay)]
    BB_bank_use <- tracking["BB_bank_use", ac(ay)]
    BB_bank <- tracking["BB_bank", ac(ay)]
    BB_borrow <- tracking["BB_borrow", ac(ay)]
    
    ### replace NAs with 0
    BB_return <- ifelse(!is.na(BB_return), BB_return, 0)
    BB_bank_use <- ifelse(!is.na(BB_bank_use), BB_bank_use, 0)
    BB_bank <- ifelse(!is.na(BB_bank), BB_bank, 0)
    BB_borrow <- ifelse(!is.na(BB_borrow), BB_borrow, 0)
    
    ### get catch target(s)
    catch_target <- ctrl@trgtArray[, "val", ]
    
    ### implement B&B
    catch_target <- catch_target - c(BB_return) + c(BB_bank_use) - 
      c(BB_bank) + c(BB_borrow)
    
    ### save in ctrl object
    ctrl@trgtArray[, "val", ] <- catch_target
    
  }
  
  ### return catch target and tracking
  return(list(ctrl = ctrl, tracking = tracking))
  
}

### ------------------------------------------------------------------------ ###
### forward projection ####
### ------------------------------------------------------------------------ ###
### including process error on stock.n
### implemented by simply multiplying numbers at age from fwd with noise factor

fwd_WKNSMSE <- function(stk, ctrl,
                        sr, ### stock recruitment model
                        sr.residuals, ### recruitment residuals
                        sr.residuals.mult = TRUE, ### are res multiplicative?
                        maxF = 2, ### maximum allowed Fbar
                        proc_res = NULL, ### process error noise,
                        dd_M = NULL, relation = NULL, ### density-dependent M
                        ...) {
  
  ### calculate density-dependent natural mortality if required
  if (!is.null(dd_M)) {
    
    ### overwrite m in the target year before projecting forward
    m(stk)[, ac(ctrl@target[, "year"])] <- calculate_ddM(stk, ctrl@target[, "year"], relation = relation)
    
  }
  
  ### project forward with FLash::fwd
  stk[] <- fwd(object = stk, control = ctrl, sr = sr, 
               sr.residuals = sr.residuals, 
               sr.residuals.mult = sr.residuals.mult,
               maxF = maxF)
  
  ### add process error noise, if supplied
  if (!is.null(proc_res)) {
    
    ### projected years
    yrs_new <- seq(from = ctrl@target[, "year"], to = range(stk)[["maxyear"]])
    
    ### workaround to get residuals
    ### they are saved in the "fitted" slot of sr...
    if (!isTRUE(proc_res == "fitted")) {
      
      stop("survival process error inacessible")
      
    } else {
      
      ### implement process error
      stock.n(stk)[, ac(yrs_new)] <- stock.n(stk)[, ac(yrs_new)] *
        fitted(sr)[, ac(yrs_new)]
      ### update stock biomass
      stock(stk) <- computeStock(stk)
      
    }
    
  }
  
  ### return stock
  return(list(object = stk))
  
}


### ------------------------------------------------------------------------ ###
### extract uncertainty from SAM object ####
### ------------------------------------------------------------------------ ###
### function for creating iterations based on estimation of uncertainty in SAM
### note: SAM works on a log scale and all reported parameters are also on a 
###       scale, even standard deviations.
###       Therefore, the values returned from this function are exponentiated

SAM_uncertainty <- function(fit, n = 1000, print_screen = FALSE,
                            idx_cov = FALSE, ### get idx covariance
                            catch_est = FALSE ### get SAM catch estimate
) {
  
  if (!is(fit, "sam")) stop("fit has to be class \"sam\"")
  
  ### index for fishing mortality ages
  idxF <- fit$conf$keyLogFsta[1, ] + 1# +dim(stk)[1]
  idxF <- idxF[idxF != 0] ### remove 0s
  
  ### index for F variances (usually some ages are bound)
  #idxNVar <- fit$conf$keyVarLogN 
  
  ### get ages used for calculating fbar
  #bAges <- fit$conf$fbarRange
  #bAges <- do.call(':', as.list(bAges))
  
  ### index for stock numbers ages
  #idxN <- 1:ncol(natural.mortality)
  #idxN <- seq(min(fit$data$minAgePerFleet), max(fit$data$maxAgePerFleet))
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
  ### get covariance matrix of all model parameters
  cov <- solve(sds$jointPrecision)
  
  ### create random values based on estimation and covariance
  sim.states <- mvrnorm(n, est, cov) ### REPLACE with rmvnorm()
  ### matrix, columns are values, rows are requested samples
  table(colnames(sim.states))
  ### contains, among others, logF, logN...
  
  ### combine SAM estimate and random samples
  #dat <- rbind(est, sim.states)
  dat <- sim.states
  
  ### workaround if only one iteration/replicate
  if (!is.matrix(dat)) dat <- t(as.matrix(dat))
  
  ### ---------------------------------------------------------------------- ###
  ### stock ####
  ### ---------------------------------------------------------------------- ###
  
  ### stock characteristics
  min_age <- min(fit$data$minAgePerFleet[fit$data$fleetTypes == 0])
  max_age <- max(fit$data$maxAgePerFleet[fit$data$fleetTypes == 0])
  years <- fit$data$years
  ### FLQuant template for stock
  stk_template <- FLQuant(dimnames = list(age = min_age:max_age, year = years,
                                          iter = 1:n))
  
  ### numbers at age
  stock.n <- stk_template
  stock.n[] <- exp(t(dat[, colnames(dat) == "logN"]))
  
  ### F at age
  harvest <- stk_template
  ### insert values for estimated ages
  harvest[unique(idxF)] <- exp(t(dat[, colnames(dat) == "logF"]))
  ### duplicate ages, if some of them are bound
  harvest[] <- harvest[idxF]
  
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
  idx_LogFpar <- lapply(idx_surveys, function(x) {
    tmp <- fit$conf$keyLogFpar[x, ] + 1
    tmp <- tmp[tmp > 0]
  })
  
  sum(colnames(dat) == "logFpar") ### there are 9 parameters for cod
  ### 5 for Q1 (ages 1-5), 4 for Q3 (ages 1-4).
  survey_ages_idx <- split(seq(length(unlist(survey_ages))), 
                           rep(seq(survey_ages), sapply(survey_ages, length))) 
  
  ### get catchability at age (time-invariant) samples
  catchability <- lapply(seq_along(idx_surveys), function(x) {
    
    ### create FLQuant template
    tmp <- FLQuant(dimnames = list(age = survey_ages[[x]],
                                   year = "all", iter = 1:n))
    ### fill with catchability values
    tmp[] <- 
      exp(t(dat[, colnames(dat) == "logFpar", drop = FALSE][, idx_LogFpar[[x]]]))
    
    return(tmp)
    
  })
  
  ### ---------------------------------------------------------------------- ###
  ### standard deviation - catch ####
  ### ---------------------------------------------------------------------- ###
  ### time-invariant
  
  ### template
  catch_sd <- FLQuant(dimnames = list(age = dimnames(stock.n)$age, year = "all",
                                      iter = 1:n))
  # sum(colnames(dat) == "logSdLogObs") 
  
  ### index for catch sd (some ages are linked)
  catch_sd_idx <- idxObs[1, ][idxObs[1, ] > -1] + 1
  
  ### extract values
  catch_sd[] <- 
    exp(t(dat[, colnames(dat) == "logSdLogObs", drop = FALSE][, catch_sd_idx]))
  ### logSdLogObs is the log of the SD of the log observations
  ### exponentiate values here to get SD (and not logSD)
  
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
    idx_sd <- idxObs[idx_surveys[x], ]
    idx_sd <- idx_sd[idx_sd > -1] + 1
    
    ### fill with catchability values
    tmp[] <- exp(t(dat[, colnames(dat) == "logSdLogObs", drop = FALSE][, idx_sd]))
    
    return(tmp)
    
  })
  
  ### ---------------------------------------------------------------------- ###
  ### surveys - get covariance ####
  ### ---------------------------------------------------------------------- ###
  
  if (isTRUE(idx_cov)) {
    . <- capture.output(survey_cov <- foreach(iter_i = 1:n) %do% {
      
      fit$obj$fn(sim.states[iter_i, 1:length(sds$par.fixed)])
      cov <- fit$obj$report()$obsCov
      
      return(cov[idx_surveys])
      
    })
    
  } else {
    
    idx_cov <- NULL
    
  }
  
  ### ---------------------------------------------------------------------- ###
  ### stock.n uncertainty ####
  ### ---------------------------------------------------------------------- ###
  ### process error
  
  ### get index for ages
  idx_SdLogN <- fit$conf$keyVarLogN + 1
  
  ### FLQuant template
  SdLogN <- catch_sd %=% NA_integer_
  
  ### fill
  SdLogN[] <- 
    exp(t(dat[, colnames(dat) == "logSdLogN", drop = FALSE][, idx_SdLogN]))
  
  ### ---------------------------------------------------------------------- ###
  ### 
  ### ---------------------------------------------------------------------- ###
  ### not used so far: 
  ### process error:
  ### logSdLogN (2 values per sim)
  ### logSdLogFsta (1)
  
  ### ---------------------------------------------------------------------- ###
  ### get catch estimates ####
  ### ---------------------------------------------------------------------- ###
  
  if (isTRUE(catch_est)) {
    
    ### catch fleet index/indices
    catch_fleets <- which(fit$data$fleetTypes == 0)
    catch_desc <- fit$data$aux
    
    
    . <- capture.output(res_n <- foreach(iter_i = 1:n, .combine = c
    ) %do% {
      
      fit$obj$fn(sim.states[iter_i, 1:length(sds$par.fixed)])
      tmp <- cbind(fit$data$aux, est = fit$obj$report()$predObs)
      tmp <- tmp[tmp[, 2] == catch_fleets, ]
      tmp <- stats::aggregate(est ~ age + year, tmp,
                              FUN = sum)
      tmp_full <- expand.grid(age = unique(tmp$age),
                              year = unique(tmp$year))
      tmp <- merge(x = tmp, y = tmp_full, all = TRUE)
      tmp <- tmp[order(tmp$year, tmp$age), ]
      return(tmp$est)
      
    })
    catch_n <- FLQuant(NA, 
                       dimnames = list(year = rownames(fit$data$catchMeanWeight), 
                                       age = colnames(fit$data$catchMeanWeight), 
                                       iter = seq(n)))
    catch_n[] <- exp(res_n)
    
    ### apply catch multiplier, if used
    if (fit$conf$noScaledYears > 0) {
      
      ### get catch multiplier dimensions
      ages_mult <- fit$conf$minAge:fit$conf$maxAge
      yrs_mult <- fit$conf$keyScaledYears
      ### get multiplier positions in created data set
      pos_mult <- which(colnames(dat) == "logScale")
      ### create catch multiplier FLQuant
      catch_mult_data <- FLQuant(NA, 
                                 dimnames = list(year = fit$conf$keyScaledYears,
                                                 age = fit$conf$minAge:fit$conf$maxAge,
                                                 iter = seq(nrow(dat))))
      catch_mult_data[] <- rep(t(dat[, pos_mult]), each = length(ages_mult))
      ### model values in log scale, exponentiate
      catch_mult_data <- exp(catch_mult_data)
      ### for simpler calculations, use catch_n as template
      catch_mult <- catch_n %=% 1
      catch_mult[ac(ages_mult), ac(yrs_mult)] <- catch_mult_data
      
      ### correct catch.n
      catch_n <- catch_n * catch_mult
      
    }
    
  } else {
    
    catch_n <- NULL
    
  }
  
  return(list(stock.n = stock.n, harvest = harvest,
              survey_catchability = catchability, catch_sd = catch_sd,
              survey_sd = survey_sd, survey_cov = survey_cov, 
              proc_error = SdLogN, catch_n = catch_n))
  
}

### ------------------------------------------------------------------------ ###
### additional functions for SAM ####
### ------------------------------------------------------------------------ ###

### extract parameters of model fit
sam_getpar <- function(fit) {
  
  ### check class
  if (is(fit, "sam") | is(fit, "list")) {
    fit <- list(fit)
  } else if (is(fit, "sam_list")) {
    ### do nothing
  } else {
    stop("unknown fit object")
  }
  
  ### go through fits and extract parameters
  pars <- foreach(fit_i = fit, .errorhandling = "pass") %do% {
    
    p <- fit_i$pl
    p$missing <- NULL
    attr(p, "what") <- NULL
    
    return(p)
    
  }
  
  ### return list or params
  if (length(pars) == 1) {
    return(pars[[1]])
  } else {
    return(pars)
  }
  
}


### ------------------------------------------------------------------------ ###
### run fmle for FLSR in parallel ####
### ------------------------------------------------------------------------ ###

fmle_parallel <- function(sr, cl, seed = NULL) {
  
  ### split into parts
  it_parts <- split(1:dim(sr)[6], cut(1:dim(sr)[6], length(cl)))
  names(it_parts) <- NULL
  
  ### split sr into junks
  sr_chunks <- lapply(it_parts, function(x) {
    FLCore::iter(sr, x)
  })
  
  ### set RNG seed
  if (!is.null(seed)) {
    clusterEvalQ(cl = cl, expr = {set.seed(seed)})
  }
  
  ### fit model to each junk
  fits <- foreach(sr_chunk = sr_chunks, .packages = "FLCore") %dopar% {
    
    fmle(sr_chunk)
    
  }
  
  ### combine into single object
  ### slots rec, ssb, covar, logerrror, model, logl, gr, distribution, inital,
  ### name, desc, range
  ### do not need to be changed
  logLik(sr) <- logLik(fits[[1]])
  sr@logLik[seq(dims(sr)$iter)] <- c(sapply(fits, logLik))
  
  fitted(sr) <- do.call(FLCore::combine, lapply(fits, "fitted"))
  
  residuals(sr) <- do.call(FLCore::combine, lapply(fits, "residuals"))
  
  vcov(sr) <- array(data = unlist(lapply(fits, vcov)), 
                    dim = c(dim(fits[[1]]@vcov)[1], 
                            dim(fits[[1]]@vcov)[2], 
                            dims(sr)$iter), 
                    dimnames = list(dimnames(fits[[1]]@vcov)[[1]], 
                                    dimnames(fits[[1]]@vcov)[[2]], 
                                    iter = ac(dimnames(sr)$iter)))
  
  params(sr) <- propagate(params(sr), dims(sr)$iter)
  params(sr)[] <- unlist(lapply(fits, function(x) {
    split(x@params@.Data, slice.index(x@params, MARGIN = 2))
  }))
  
  hessian(sr) <- array(data = unlist(lapply(fits, hessian)), 
                       dim = c(dim(fits[[1]]@hessian)[1], 
                               dim(fits[[1]]@hessian)[2], 
                               dims(sr)$iter), 
                       dimnames = list(dimnames(fits[[1]]@hessian)[[1]], 
                                       dimnames(fits[[1]]@hessian)[[2]], 
                                       iter = ac(dimnames(sr)$iter)))
  
  details(sr) <- details(fits[[1]])
  
  return(sr)
  
}

### ------------------------------------------------------------------------ ###
### inter-annual variability ####
### ------------------------------------------------------------------------ ###
#' calculate inter-annual variability of FLQuant
#'
#' This function calculates survey indices from the numbers at age of an 
#' FLStock object
#'
#' @param object Object of class \linkS4class{FLQuant} with values.
#' @param period Select every n-th year, e.g. biennial (optional).
#' @param from,to Optional year range for analysis.
#' @param summary_per_iter Function for summarising per iter. Defaults to mean.
#' @param summary Function for summarising over iter. Defaults to mean.
#' @return An object of class \code{FLQuant} with inter-annual variability.
#'
#' @export
#' 
setGeneric("iav", function(object, period, from, to, summary_per_iter, 
                           summary_year, summary_all) {
  standardGeneric("iav")
})

### object = FLQuant
#' @rdname iav
setMethod(f = "iav",
          signature = signature(object = "FLQuant"),
          definition = function(object, 
                                period, ### periodicity, e.g. use every 2nd value 
                                from, to,### year range
                                summary_per_iter, ### summarise values per iteration
                                summary_year,
                                summary_all) {
            
            ### subset years
            if (!missing(from)) object <- window(object, start = from)
            if (!missing(to)) object <- window(object, end = from)
            
            ### get years in object
            yrs <- dimnames(object)$year
            
            ### select every n-th value, if requested
            if (!missing(period)) {
              yrs <- yrs[seq(from = 1, to = length(yrs), by = period)]
            }
            
            ### reference years
            yrs_ref <- yrs[-length(yrs)]
            ### years to compare
            yrs_comp <- yrs[-1]
            
            ### calculate variation (absolute values, ignore pos/neg)
            res <- abs(1 - object[, yrs_comp] / object[, yrs_ref])
            
            ### replace Inf with NA (compared to 0 catch)
            res <- ifelse(is.finite(res), res, NA)
            
            ### summarise per iteration
            if (!missing(summary_per_iter)) {
              res <- apply(res, 6, summary_per_iter, na.rm = TRUE)
            }
            
            ### summarise per year
            if (!missing(summary_year)) {
              res <- apply(res, 1:5, summary_year, na.rm = TRUE)
            }
            
            ### summarise over everything
            if (!missing(summary_all)) {
              
              res <- summary_all(res, na.rm = TRUE)
              
            }
            
            return(res)
            
          })

### ------------------------------------------------------------------------ ###
### density-dependent natural mortality ####
### ------------------------------------------------------------------------ ###
### calculate density-dependent M from M2 relationships and stock.n 

calculate_ddM <- function(stk,
                          yr,
                          relation) {
  
  ### read in M2 relationships
  M2 <- relation
  M2$pM2 <- M2$Nprey <- NA
  M2 <- array(rep(unlist(M2), dim(stk)[6]), dim=c(nrow(M2), ncol(M2), dim(stk)[6]))
  dimnames(M2)[[2]] <- c("age","predator","intercept","logbPred","logbPrey","nPred","nPrey","pM2") 
  
  M_new <- m(stk)[, ac(yr)] %=% 0
  
  for (iter_i in seq(dim(stk)[6])) {
    
    ### extract number of predator and prey cod-at-age from stock.n in the specified year (000s)
    M2[,"nPrey", iter_i] <- stock.n(stk)[M2[,"age",iter_i], ac(yr),,,,iter_i]
    M2[is.na(M2[,"nPred", iter_i]) ,"nPred", iter_i] <- stock.n(stk)[M2[!is.na(M2[,"predator",iter_i]), "predator", iter_i], ac(yr),,,,iter_i]
    
    ### calculate partial M2s from predator and prey abundances
    M2[, "pM2", iter_i] <- M2[, "intercept", iter_i] + M2[, "logbPrey", iter_i] * log(M2[, "nPrey", iter_i])
    M2[!is.na(M2[, "logbPred", iter_i]), "pM2", iter_i] <- M2[!is.na(M2[, "logbPred", iter_i]), "pM2", iter_i] + M2[!is.na(M2[, "logbPred", iter_i]), "logbPred", iter_i] * log(M2[!is.na(M2[, "logbPred", iter_i]), "nPred", iter_i])
    M2[, "pM2", iter_i] <- exp(M2[, "pM2", iter_i])
    
    ### sum M2s for each prey age class
    M2age <- aggregate(M2[,"pM2", iter_i], by = list(age = M2[,"age",iter_i]), sum)
    
    ### overwrite m in the specified year
    M_new[M2age$age,,,,,iter_i] <- M2age$x
    
  }
  ### and add 0.2 for non-predation mortality
  M_new <- M_new + 0.2
  
  return(M_new)
  
}