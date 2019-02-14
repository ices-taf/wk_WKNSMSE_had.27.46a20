### management implementation: short term forecast with SAM ####
### ------------------------------------------------------------------------ ###
### short term forecast with SAM
### including TAC constraint
is_WKNSMSE_JarworskiGrowth <- function(stk, tracking, ctrl,
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
                       BB_conditional = FALSE, ### check status before BB
                       BB_rho, ### definition of BB
                       ### reference points
                       hcrpars = list(),
                       ...) {
  
  JarworskiGrowthFrcst<-function(dat.wt,is.landings=F){
    
    #ages
    ages<-as.numeric(colnames(dat.wt))
    years<-as.numeric(rownames(dat.wt))
    plus.group<-max(as.numeric(ages))
    
    # change zeros to NA
    dat.wt[dat.wt==0]<-NA
    
    # extend st.wt object and create new data frame for model results
    dat.wt.ext<-rbind(dat.wt,as.data.frame(matrix(NA,ncol=length(ages),nrow=3,
                                                  dimnames=list((c(years[length(years)]+1):(years[length(years)]+3)),as.character(ages)))))
    # results matrix
    dat.wt.mod<-dat.wt.ext
    
    
    # make cohort key
    cohort<-dat.wt.ext*NA
    cohort$"0"<-as.numeric(row.names(dat.wt.ext))
    
    # add cohort year to cell
    co.idx<-function(x,ages){
      x2<-rep(x["0"],times=length(ages))-as.numeric(ages[])
      return(x2)
    }
    
    cohort<-apply(cohort,1,co.idx,ages)
    
    # bit of rearrangement needed
    cohort<-(t(cohort))
    colnames(cohort)<-as.character(ages)
    
    
    # calculate cohort weights
    
    #fill wts 
    start.yr<-(max(as.numeric(row.names(dat.wt.ext)))-2)-max(as.numeric(ages))
    
    # which years
    yr.lst<-as.list(c(as.numeric(start.yr:max(as.numeric(row.names(dat.wt))))))
    
    # model growth per cohort
    new.wts<-lapply(yr.lst,function(x,ages.=ages,dat.wt.in=dat.wt.ext,cohort.=cohort){
      
      tom<-data.frame(ages=ages.)  
      tom$wts<-NA
      tom$wts[1:length(dat.wt.in[cohort. == x])]<-dat.wt.in[cohort. == x]
      
      # only use ages 0:8 for regression
      if (sum(!(is.na(tom$wts)))>=3){
        params<-lm(wts~ages,tom[tom$ages %in% ages,])
        
        tom$wts[is.na(tom$wts)]<-params[[1]][1]+ params[[1]][2]*tom$ages[is.na(tom$wts)]
      }
      
      return(tom$wts)
    })
    
    
    new.wts<-matrix(unlist(new.wts),ncol=length(ages),nrow=length(yr.lst),byrow=T)
    
    
    # put wts into results matrix
    for (x in yr.lst){
      dat.wt.mod[cohort == x]<-new.wts[yr.lst == x][1:length(dat.wt.mod[cohort == x])]
    }
    
    # for age 0 - only keep actual observed weights
    dat.wt.mod[1:dim(dat.wt)[1],"0"]<-dat.wt[,c("0")] # only really needed when doing landings - discards and stock should have observed weights
    
    
    # fill with 3 year average for cohorts with less than 3 data points in the last 3 years
    ave<-colMeans(dat.wt.mod[c((dim(dat.wt.mod)[1]-5):(dim(dat.wt.mod)[1]-3)),],na.rm=T)
    
    # fill ages 0, 1 and 2 with 3 year average
    dat.wt.mod[c((dim(dat.wt.mod)[1]-2):dim(dat.wt.mod)[1]),c("0")]<-ave["0"]
    dat.wt.mod[c((dim(dat.wt.mod)[1]-2):dim(dat.wt.mod)[1]),c("1")]<-ave["1"]
    dat.wt.mod[c((dim(dat.wt.mod)[1]-2):dim(dat.wt.mod)[1]),c("2")]<-ave["2"]  
    
    
    if (isTRUE(is.landings)){
      # for landings we fill age 3 with 3 year average
      dat.wt.mod[c((dim(dat.wt.mod)[1]-2):dim(dat.wt.mod)[1]),c("3")]<-ave["3"]  
      
      # age 4 - use 3 year average but only fill last 2 years
      dat.wt.mod[c((dim(dat.wt.mod)[1]-1):dim(dat.wt.mod)[1]),c("4")]<-ave["4"]  
      
      #age 5 - use last 3 year average but only fill last year
      dat.wt.mod[c((dim(dat.wt.mod)[1]):dim(dat.wt.mod)[1]),c("5")]<-ave["5"] 
      
    }else{
      # age 3 - use 3 year average but only fill last 2 years
      dat.wt.mod[c((dim(dat.wt.mod)[1]-1):dim(dat.wt.mod)[1]),c("3")]<-ave[4]  
      
      #age 4 - use last 3 year average but only fill last year
      dat.wt.mod[c((dim(dat.wt.mod)[1]):dim(dat.wt.mod)[1]),c("4")]<-ave[5]  
      
      
    }
    
    #then fill any others which don't have an average
    idx<-is.na(dat.wt.mod[c((dim(dat.wt.mod)[1]-2):(dim(dat.wt.mod)[1])),])
    if(sum(idx)>0){ 
      idx2<-which(is.na(dat.wt.mod[c((dim(dat.wt.mod)[1]-2):(dim(dat.wt.mod)[1])),]),arr.ind=T)
      dat.wt.mod[idx2]<-ave[idx2[,2]]
    }
    
    
    
    # round to 3dp
    idx<-(nrow(dat.wt.mod)-2):nrow(dat.wt.mod)
    frcst.wts<-round(dat.wt.mod[idx,],digits=3)
    
    # fill NA in age 0 with 0 for landings
    if(isTRUE(is.landings)){
      frcst.wts[,c("0")]<-0.000
    }
    return(frcst.wts)
  }
  
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
  if (isTRUE(BB) & isTRUE(BB_conditional)) {
    
    ### duplicate last target value
    fwd_trgt <- c(fwd_trgt, tail(fwd_trgt, 1))
    
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
                  ### currently, can only do F status quo and F target from ctrl object
                  fscale <- ifelse(fwd_trgt == "fsq", 1, NA) ### scaled Fsq
                  ### target F values
                  fval <- ifelse(fwd_trgt == "hcr", ctrl@trgtArray[, "val", iter_i], NA) 
                  
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
                  

		  ### update mean weights with Jarworski growth model values
                  st.wt<-fit_i$data$stockMeanWeight[ac(min(fit_i$data$years):(ay-1)),]
                  ld.wt<-fit_i$data$landMeanWeight[ac(min(fit_i$data$years):(ay-1)),]
                  ds.wt<-fit_i$data$disMeanWeight[ac(min(fit_i$data$years):(ay-1)),]
                  ct.wt<-fit_i$data$catchMeanWeight[ac(min(fit_i$data$years):(ay-1)),]
  
                  ### start function
                  st.wt.fcst<-JarworskiGrowthFrcst(dat.wt=st.wt,is.landings=F)
                  ld.wt.fcst<-JarworskiGrowthFrcst(dat.wt=ld.wt,is.landings=T)
                  ds.wt.fcst<-JarworskiGrowthFrcst(dat.wt=ds.wt,is.landings=F)
                  ct.wt.fcst<-JarworskiGrowthFrcst(dat.wt=ct.wt,is.landings=F)
                    
                  fit_i$data$stockMeanWeight<-rbind(st.wt,st.wt.fcst)
                  fit_i$data$landMeanWeight<-rbind(ld.wt,ld.wt.fcst)
                  fit_i$data$disMeanWeight<-rbind(ds.wt,ds.wt.fcst)
                  fit_i$data$catchMeanWeight<-rbind(ct.wt,ct.wt.fcst)

                  ### forecast 
                  fc <- stockassessment::forecast(fit = fit_i, fscale = fscale, fval = fval,
                                                  ave.years = ave.years,
                                                  rec.years = rec.years,
                                                  overwriteSelYears = overwriteSelYears,
                                                  splitLD = fwd_splitLD)
                  
                  ### return forecast table
                  return(attr(fc, "tab"))#[ac(ay + 1), "catch:median"])
                  
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
  hcrpars <- FLPar(hcrpars)
  if (dim(hcrpars)[2] == 1) hcrpars <- propagate(hcrpars, it)
  
  ### ---------------------------------------------------------------------- ###
  ### TAC constraint ####
  if (isTRUE(TAC_constraint)) {
    
    ### target year
    yr_target <- ctrl@target$year
    
    ### get previous target from tracking
    catch_prev <- tracking["metric.is", ac(yr_target - 1), drop = TRUE]
    
    ### change in advice in %
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
    BB_rho_i <- tail(rep(BB_rho, length.out = (ay - genArgs$y0 + 1)), 1)
    
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
    if (isTRUE(BB_conditional)) {
      
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
      pos_neg1 <- which(SSB_TACyr < c(hcrpars["Bpa"]) & 
                          F_TACyr > c(hcrpars["Fpa"]))
      ### check if SSB below Bpa in TAC year and year after
      pos_neg2 <- which(SSB_TACyr < c(hcrpars["Bpa"]) & 
                          SSB_TACyr1 < c(hcrpars["Bpa"]))
      ### combine both conditions
      pos_neg <- union(pos_neg1, pos_neg2)
      
      ### if status evaluated as not precautionary
      ### stop banking and borrowing
      ### (but still pay back/use from last year)
      BB_bank[pos_neg] <- 0
      BB_borrow[pos_neg] <- 0
      
    }
    
    # ### correct target catch later in iem module
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