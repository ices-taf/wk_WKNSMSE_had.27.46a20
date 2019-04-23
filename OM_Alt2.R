### ------------------------------------------------------------------------ ###
### Alternative OM2 for had.27.46a20 ####
### ------------------------------------------------------------------------ ###
### base on SAM assessment
# fix regularity of recruitment spikes

rm(list=ls())
#windows()
# set directories
#setwd("N:\\Stock assessment\\Haddock WKNSMSE 2018\\WK_WKNSMSE_had.27.46a20")

# ### check versions of required R packages
if (packageVersion("FLCore") < "2.6.11.9001") 
  stop("please update FLCore")
if (packageVersion("FLfse") < "0.0.0.9003") 
  stop("please update FLfse")
if (packageVersion("stockassessment") < "0.8.1") 
  stop("please update stockassessment")
if (packageVersion("mse") < "0.9.1") 
  stop("please update mse")

### load packages
library(FLfse)
library(stockassessment)
library(ggplotFL)
library(FLAssess)
library(mse)
### load files from package mse for easier debugging
#devtools::load_all("../mse/")
library(FLash)
library(tidyr)
library(dplyr)

source("a4a_mse_WKNSMSE_funs.R")

dir.create(path = "input/had4", recursive = TRUE)
dir.create(path = "output/runs/had4", recursive = TRUE)

### source the scripts from functions folder
# also sources TSA functions "NS haddock functions 20190116.r" for TSA uncertainty calcs
invisible(lapply(list.files(path = "functions/", pattern = "*.R$", 
                            full.names = TRUE), source))


### ------------------------------------------------------------------------ ###
### 1. Set simulation specifications ####
### ------------------------------------------------------------------------ ###

### number of iterations
n <- 1000 #1000
### number of years
n_years <- 20 #30

### last data year
yr_data <- 2018

omName<-"Alt2"
stkname<-"NoSh Haddock"

verbose<-F

# set some useful index values
eqSim_yrs<-(yr_data-10):(yr_data-1)
rec.period<-2000:2017
hist_yrs<-1972:2017
ages<-0:8
sim_yrs<-(yr_data):(yr_data+n_years)


### ------------------------------------------------------------------------ ###
### 2. WGNSSK 2018 results ####
### ------------------------------------------------------------------------ ###
# recreate SAM fit from WGNSSK
# use input data provided in FLfse

# start time series in 1972 as done with TSA
had4_stk<-window(had4_stk,start=1972)
had4_idx[[1]]<-window(had4_idx[[1]],start=1983)

# check input data is the same as used in TSA - it is!

# load in from WGNSSK 2018 results
load("input/had4/had2746a20_FLStock objects_2018.Rdata")
# catch includes the 4 components. 
# convert discards to include bms and ibc components

# update numbers
# use numbers to fidn weighted mean wt at age
# recalc "discards" yield
x.pg<-window(x.pg,end=2018)
x.bms.pg<-window(x.bms.pg,end=2018)
x.ibc.pg<-window(x.ibc.pg,end=2018)
disN<-x.pg@discards.n+x.bms.pg@landings.n+x.ibc.pg@landings.n
disW<-(x.pg@discards.wt*x.pg@discards.n+
         x.bms.pg@landings.wt*x.bms.pg@landings.n+
         x.ibc.pg@landings.wt*x.ibc.pg@landings.n)/disN
disW[is.na(disW)]<-0

#update had4_stk
had4_stk@discards.n<-disN
had4_stk@discards.wt<-disW
had4_stk@discards <- computeDiscards(had4_stk)


### set units
units(had4_stk)[1:17] <- as.list(c(rep(c("t", "1000", "kg"), 4),
                                   "NA", "NA", "f", "NA", "NA"))

### save for later comparison
stk_orig <- had4_stk

# rename had4_stk
stk<-had4_stk

rm(list=c("x.pg","x.bms.pg","x.ibc.pg","disN","disW"))

### ------------------------------------------------------------------------ ###
### 3. SAM fit ####
### ------------------------------------------------------------------------ ###

#update conf file to optimised SAM
had4_conf_sam$fracMixF<-0
had4_conf_sam$fracMixN<-0
had4_conf_sam$fracMixObs<-c(0,0,0)

# fit SAM
fit <- FLR_SAM(stk = had4_stk, idx = had4_idx, conf = had4_conf_sam)

### fit SAM as it is done during the MSE simulation
fit_est <- FLR_SAM(stk = had4_stk, idx = had4_idx, 
                   conf = had4_conf_sam, 
                   newtonsteps = 0, rel.tol = 0.001)
### extract model parameters and use them in the simulation as starting values
sam_initial <- sam_getpar(fit_est)

### ------------------------------------------------------------------------ ###
### 4. create FLStock ####
### ------------------------------------------------------------------------ ###
### create template with 1 iteration

### the results (stock numbers & harvest) are used from the real WGNSSK fit
stk <- SAM2FLStock(object = fit, stk = had4_stk)

### set units
units(stk)[1:17] <- as.list(c(rep(c("t", "1000", "kg"), 4),
                              "", "", "f", "", ""))

### save for later comparison
stk_orig <- stk

### ------------------------------------------------------------------------ ###
### 5. add uncertainty ####
### ------------------------------------------------------------------------ ###
### use variance-covariance


### add iteration dimension
stk <- FLCore::propagate(stk, n)
dim(stk)

### add uncertainty estimated by SAM as iterations
set.seed(1)
uncertainty <- SAM_uncertainty(fit = fit, n = n, print_screen = FALSE, 
                               idx_cov = TRUE, catch_est = TRUE)
### add noise to stock
stock.n(stk)[] <- uncertainty$stock.n
stock(stk)[] <- computeStock(stk)
### add noise to F
harvest(stk)[] <- uncertainty$harvest

### catch noise added later

### maximum observed F
max(fbar(stk))
# 1.652771 in year 1991 
max(harvest(stk))
# 2.270628 in year 1972 for age 7 

### get estimated catch numbers
catch_n <- uncertainty$catch_n

### ------------------------------------------------------------------------ ###
### 6. extend stock for MSE simulation ####
### ------------------------------------------------------------------------ ###

### stock weights, M etc in 2018 based on a three year average to enable 
### calculation of SSB. Though these will get overwritten for haddock as we use cohort growth
### although SAM estimates F in 2018, this is not reported or taken forward into
### forcasts by the WG

stk_stf2018 <- stf(window(stk, end = 2018), n_years)
### use all available data
stk_stf <- stk_stf2018


### ------------------------------------------------------------------------ ###
### 7. biological data for OM ####
### ------------------------------------------------------------------------ ###

### Resample weights, maturity and natural mortality from the last 10 years 
### (2008-2017)
### set up an array with one resampled year for each projection year 
### (including intermediate year) and replicate
### use the same resampled year for all biological parameters
### this is the approach used in eqsim for Northern shelf haddock

set.seed(2)

### use last five data years to sample biological parameters
sample_yrs <- eqSim_yrs
### get year position of sample years
sample_yrs_pos <- which(dimnames(stk_stf)$year %in% sample_yrs)

sel_sample_yrs <- eqSim_yrs[6:10]
### get year position of sample years
sel_sample_yrs_pos <- which(dimnames(stk_stf)$year %in% sel_sample_yrs)

### create samples for biological data (weights, etc.)
### the historical biological parameters are identical for all iterations
### and consequently do not need to be treated individually
### (but keep age structure)
### create vector with resampled years
bio_samples <- sample(x = sample_yrs_pos, 
                      size = (n_years + 1) * n, replace = TRUE)
### do the same for selectivity
sel_samples <- sample(x = sel_sample_yrs_pos, 
                      size = (n_years + 1) * n, replace = TRUE)
### years to be populated
bio_yrs <- which(dimnames(stk_stf)$year %in% 2018:dims(stk_stf)$maxyear)


### insert values
catch.wt(stk_stf)[, bio_yrs] <- c(catch.wt(stk)[, bio_samples,,,, 1])
stock.wt(stk_stf)[, bio_yrs] <- c(stock.wt(stk)[, bio_samples,,,, 1])
landings.wt(stk_stf)[, bio_yrs] <- c(landings.wt(stk)[, bio_samples,,,, 1])
discards.wt(stk_stf)[, bio_yrs] <- c(discards.wt(stk)[, bio_samples,,,, 1])
m(stk_stf)[, bio_yrs] <- c(m(stk)[, bio_samples,,,, 1])
mat(stk_stf)[, bio_yrs] <- c(mat(stk)[, bio_samples,,,, 1])

### use different samples for selectivity
harvest(stk_stf)[, bio_yrs] <- c(harvest(stk)[, sel_samples,,,, 1])

### ------------------------------------------------------------------------ ###
### 8. stock recruitment ####
### ------------------------------------------------------------------------ ###
### fit hockey-stick model
### get residuals from smoothed residuals

### use only data from 2000 and later
sr <- as.FLSR(window(stk_stf, start = min(rec.period)), model = "segreg")
### fit model individually to each iteration and suppress output to screen
suppressWarnings(. <- capture.output(sr <- fmle(sr)))
sr_fit<-sr
### run in parallel
# library(doParallel)
# cl <- makeCluster(10)
# registerDoParallel(cl)
# sr <- fmle_parallel(sr, cl)

# ### run again for failed iterations
# pos_error <- which(is.na(params(sr)["a"]))
# sr_corrected <- fmle(FLCore::iter(sr, pos_error))
# sr[,,,,, pos_error] <- sr_corrected[]
# params(sr)[, pos_error] <- params(sr_corrected)
 
### generate residuals for MSE
### years with missing residuals
yrs_res <- dimnames(sr)$year[which(is.na(iterMeans(rec(sr))))]
spike.yrs<-c(2005,2009,2014)

res_first <- foreach(iter_i = seq(dim(sr)[6]), .packages = "FLCore", 
                     #.errorhandling = "pass") %dopar% {
                     .errorhandling = "pass") %do% {  
                       set.seed(iter_i)
                       
                       ### get residuals for current iteration
                       res_i <- c(FLCore::iter(residuals(sr), iter_i))
                       res_i <- res_i[!is.na(res_i)]
                       #exclude spike yrs
                       res_i <- res_i[-which(dimnames(residuals(sr))$year %in% spike.yrs)]
                       
                       ### calculate kernel density of residuals
                       density <- density(x = res_i)
                       ### sample residuals
                       mu <- sample(x = res_i, size = length(yrs_res), replace = TRUE)
                       ### "smooth", i.e. sample from density distribution
                       res_first <- rnorm(n = length(yrs_res), mean = mu, sd = density$bw)
                       
                       return(res_first)
                       
                     }


# second stage - model spikes
mean.int<-mean(diff(spike.yrs))
n.spikes<-ceiling(diff(range(sim_yrs))/mean.int)

res_mod<-foreach(iter_i = seq(dim(sr)[6]), .packages = "FLCore", 
                 .errorhandling = "pass") %do% {
                   
                   set.seed(iter_i^2)
                   
                   # calc spike yrs
                   x<-runif(n=n.spikes, min = 0, max = 1)
                   s<-0.5 # spikes are fairly regular
                   y<-round(mean.int*(s*x+1-(s/2)))
                   future.spike.yrs<-cumsum(c(max(spike.yrs),y))
                   future.spike.yrs<-future.spike.yrs[future.spike.yrs>2017]
                   
                   
                   ### get residuals for spike yrs in current iteration
                   res_i <- c(FLCore::iter(residuals(sr), iter_i))
                   res_i <- res_i[!is.na(res_i)]
                   res_i <- res_i[dimnames(residuals(sr))$year %in% spike.yrs]
                   
                   # get baseline residuals and add spike residuals in for future spike yrs
                   res_orig<-res_first[[iter_i]]
                   idx<-which(yrs_res %in% future.spike.yrs)
                   
                   ### sample residuals
                   mu <- sample(x = res_i, size = length(future.spike.yrs), replace = TRUE)
                   ### "smooth", i.e. sample from density distribution
                   res_new <- res_orig
                   res_new[idx]<-mu
                   
                   res_mod<-list(res_new,future.spike.yrs)
                   return(res_mod)
                   
                 }

res_new<-lapply(res_mod,function(x){return(x[[1]])})
future.spike.yrs<-lapply(res_mod,function(x){return(x[[2]])})

summary(exp(unlist(res_new)))
### insert into model
residuals(sr)[, ac(yrs_res)] <- unlist(res_new)
### exponeniate residuals to get factor
residuals(sr) <- exp(residuals(sr))
sr_res <- residuals(sr)

if (isTRUE(verbose)) plot(sr_res)

if(verbose){
  windows()
  plot(sr_res)
  savePlot(paste0("output/OM_plots/OM_",omName,"_",stkname,"_recruitment - residuals"),type="png",res=600)
  
  i_samp <- sample(seq(dim(sr)[6]), min(c(n,10)), replace=FALSE)
  
  ggplot(data=iter(sr_res,i_samp), aes(year, data)) + geom_line(aes(group=iter, colour=factor(iter))) + 
    facet_wrap(~iter, scales="free", nrow=5) + theme(legend.position = "none") 
  ggsave(paste0("output/OM_plots/OM_",omName,"_",stkname,"_recruitment - worm plot.png"), height=7, width=11,dpi=300)
  
  windows()
  layout(cbind(c(1,2,3),c(4,5,6)))
  tmp<-unlist(lapply(future.spike.yrs,function(x){return(x[1])}))
  hist(tmp,breaks=sim_yrs,col="grey80",main="Timing of first spike",xlab="")
  tmp<-unlist(lapply(future.spike.yrs,function(x){return(x[2])}))
  hist(tmp,breaks=sim_yrs,col="grey80",main="Timing of second spike",xlab="")
  tmp<-unlist(lapply(future.spike.yrs,function(x){return(x[3])}))
  hist(tmp,breaks=sim_yrs,col="grey80",main="Timing of third spike",xlab="")
  tmp<-unlist(lapply(future.spike.yrs,function(x){return(x[4])}))
  hist(tmp,breaks=sim_yrs,col="grey80",main="Timing of fourth spike",xlab="")
  tmp<-unlist(lapply(future.spike.yrs,function(x){return(x[5])}))
  tmp[tmp>max(sim_yrs)]<-NA
  hist(tmp,breaks=sim_yrs,col="grey80",main="Timing of fifth spike",xlab="")
  
  savePlot(paste0("output/OM_plots/OM_",omName,"_",stkname,"_recruitment - spike timing"),type="png",res=600)
  
}

### ------------------------------------------------------------------------ ###
### 9. process noise ####
### ------------------------------------------------------------------------ ###
### create FLQuant with process noise
### this will be added to the values obtained from fwd() in the MSE

### create noise for process error
set.seed(3)
proc_res <- stock.n(stk_stf) %=% 0 ### template FLQuant
proc_res[] <- stats::rnorm(n = length(proc_res), mean = 0, 
                           sd = uncertainty$proc_error)

### the proc_res values are on a normale scale,
### exponentiate to get log-normal 
proc_res <- exp(proc_res)
### proc_res is a factor by which the numbers at age are multiplied

### for historical period, numbers already include process error from SAM
### -> remove deviation
proc_res[, dimnames(proc_res)$year <= 2017] <- 1

### remove deviation for first age class (recruits)
proc_res[1, ] <- 1

### try saving in stock recruitment model ... 
### this gets passed on to the projection module
fitted(sr) <- proc_res

### ------------------------------------------------------------------------ ###
### 10. stf for 2018: assume catch advice is taken ####
### ------------------------------------------------------------------------ ###
c2018 <- 48990
ctrl <- fwdControl(data.frame(year = 2018, quantity = "catch", 
                              val = c2018))

### project forward for intermediate year (2018)
stk_int <- stk_stf
stk_int[] <- fwd(stk_stf, ctrl = ctrl, sr = sr, sr.residuals = sr_res,
                 sr.residuals.mult = TRUE, maxF = 5)[]
### add process noise
stock.n(stk_int) <- stock.n(stk_int) * proc_res
stock(stk_int)[] <- computeStock(stk_int)

### create stock for MSE simulation
stk_fwd <- stk_stf
### insert values for 2018
stk_fwd[, ac(2018)] <- stk_int[, ac(2018)]
### insert stock number for 2019 in order to calculate SSB at beginning of 
### 2019
stock.n(stk_fwd)[, ac(2019)] <- stock.n(stk_int)[, ac(2019)]
stock(stk_fwd)[, ac(2019)] <- computeStock(stk_fwd[, ac(2019)])

#all.equal(window(stk_fwd, end = 2018), window(stk_stf, end = 2018))


### ------------------------------------------------------------------------ ###
### 11. biological data for OEM ####
### ------------------------------------------------------------------------ ###

### base on OM stock
stk_oem <- stk_fwd

### projection years
proj_yrs <- 2018:range(stk_oem)[["maxyear"]]

### use means of sampled values for projection period
catch.wt(stk_oem)[, ac(proj_yrs)] <- 
  yearMeans(catch.wt(stk_oem)[, ac(sample_yrs)])
landings.wt(stk_oem)[, ac(proj_yrs)] <- 
  yearMeans(landings.wt(stk_oem)[, ac(sample_yrs)])
discards.wt(stk_oem)[, ac(proj_yrs)] <- 
  yearMeans(discards.wt(stk_oem)[, ac(sample_yrs)])
stock.wt(stk_oem)[, ac(proj_yrs)] <- 
  yearMeans(stock.wt(stk_oem)[, ac(sample_yrs)])
m(stk_oem)[, ac(proj_yrs)] <- yearMeans(m(stk_oem)[, ac(sample_yrs)])
mat(stk_oem)[, ac(proj_yrs)] <- yearMeans(mat(stk_oem)[, ac(sample_yrs)])

### remove stock assessment results
stock.n(stk_oem)[] <- stock(stk_oem)[] <- harvest(stk_oem)[] <- NA


### ------------------------------------------------------------------------ ###
### 12. indices ####
### ------------------------------------------------------------------------ ###
### use real FLIndices object as template (included in FLfse)
idx <- had4_idx
### extend for simulation period
idx <- window(idx, end = yr_data + n_years)
### add iterations
idx <- lapply(idx, propagate, n)

### insert catchability
for (idx_i in seq_along(idx)) {
  
  ### set catchability for projection
  index.q(idx[[idx_i]])[] <- uncertainty$survey_catchability[[idx_i]]
  
}
### create copy of index with original values
idx_raw <- lapply(idx ,index)
### calculate index values
idx <- calc_survey(stk = stk_fwd, idx = idx)

### create deviances for indices
### first, get template
idx_dev <- lapply(idx, index)
### create random noise based on sd
set.seed(4)
for (idx_i in seq_along(idx_dev)) {
  ### insert sd
  idx_dev[[idx_i]][] <- uncertainty$survey_sd[[idx_i]]
  ### noise
  idx_dev[[idx_i]][] <- stats::rnorm(n = length(idx_dev[[idx_i]]),
                                     mean = 0, sd = idx_dev[[idx_i]])
  ### exponentiate to get from normal to log-normal scale
  idx_dev[[idx_i]] <- exp(idx_dev[[idx_i]])
}


### modify residuals for historical period so that index values passed to 
### stock assessment are the ones observed in reality
### IBTS Q1, values up to 2018
idx_dev$NSQ1[, dimnames(idx_dev$NSQ1)$year <= 2018] <- 
  idx_raw$NSQ1[, dimnames(idx_raw$NSQ1)$year <= 2018] /
  index(idx$NSQ1)[, dimnames(idx$NSQ1@index)$year <= 2018]
### IBTS Q3, values up to 2017
idx_dev$NSQ3[, dimnames(idx_dev$NSQ3)$year <= 2017] <- 
  idx_raw$NSQ3[, dimnames(idx_raw$NSQ3)$year <= 2017] /
  index(idx$NSQ3)[, dimnames(idx$NSQ3@index)$year <= 2017]



### check survey
# idx0 <- calc_survey(stk = stk_fwd, idx = idx)
# idx0 <- lapply(seq_along(idx0), function(idx_i) {
#   idx_tmp <- idx0[[idx_i]]
#   index(idx_tmp) <- index(idx_tmp) * idx_dev[[idx_i]]
#   return(idx_tmp)
# })
# plot(index(idx0[[2]]))


### ------------------------------------------------------------------------ ###
### 13. catch noise ####
### ------------------------------------------------------------------------ ###
### take estimates from sam: uncertainty$catch_sd is "logSdLogObs"
### assume catch observed by SAM in projection is log-normally distributed
### around operating model catch

### create noise for catch
set.seed(5)
catch_res <- catch.n(stk_fwd) %=% 0 ### template FLQuant
catch_res[] <- stats::rnorm(n = length(catch_res), mean = 0, 
                            sd = uncertainty$catch_sd)
### the catch_res values are on a normale scale,
### exponentiate to get log-normal 
catch_res <- exp(catch_res)
### catch_res is a factor by which the numbers at age are multiplied

### for historical period, pass on real observed catch
### -> remove deviation
catch_res[, dimnames(catch_res)$year <= 2017] <- 1



### ------------------------------------------------------------------------ ###
### save OM ####
### ------------------------------------------------------------------------ ###

### path
### path
input_path <- paste0("input/had4/", n, "_", n_years, "/")
dir.create(input_path)
### stock
saveRDS(stk_fwd, file = paste0(input_path,omName, "_stk.rds"))
### stock recruitment
saveRDS(sr, file = paste0(input_path, omName, "_sr.rds"))
### recruitment residuals
saveRDS(sr_res, file = paste0(input_path, omName, "_sr_res.rds"))
### surveys
saveRDS(idx, file = paste0(input_path, omName, "_idx.rds"))
saveRDS(idx_dev, file = paste0(input_path, omName, "_idx_dev.rds"))
### catch noise
saveRDS(catch_res, file = paste0(input_path, omName, "_catch_res.rds"))
### process error
saveRDS(proc_res, file = paste0(input_path, omName, "_proc_res.rds"))
### observed stock
saveRDS(stk_oem, file = paste0(input_path, omName, "_stk_oem.rds"))
### sam initial parameters
saveRDS(sam_initial, file = paste0(input_path, omName, "_sam_initial.rds"))
### sam configuration
saveRDS(had4_conf_sam, file = paste0(input_path, omName, "_had4_conf_sam.rds"))
### catch numbers
saveRDS(catch_n, file = paste0(input_path, omName, "_catch_n.rds"))


#save.image(file = paste0(input_path, "image.RData"))
save.image(file=paste0("input/had4/image_OM",omName,"_",n,".RData"))
load(file=paste0("input/had4/image_OM",omName,"_",n,".RData"))

### ------------------------------------------------------------------------ ###
### prepare objects for new a4a standard mse package ####
### ------------------------------------------------------------------------ ###
### https://github.com/flr/mse

if(0){
### reference points
refpts_mse <- list(Btrigger = 132000,
                   Blim=94000,
                   #  Flim=0.38,
                   Ftrgt = 0.194,
                   Fpa = 0.274,
                   Bpa = 132000)

### some specifications for short term forecast with SAM
had4_stf_def <- list(fwd_yrs_average = -3:-1, #3 year average
                    fwd_yrs_rec_start = 2000,
                    fwd_yrs_sel = -1, #same as last year
                    fwd_yrs_lf_remove = -2:-1,
                    fwd_splitLD = FALSE)

### some arguments (passed to mp())
genArgs <- list(fy = dims(stk_fwd)$maxyear, ### final simulation year
                y0 = dims(stk_fwd)$minyear, ### first data year
                iy = yr_data, ### first simulation (intermediate) year
                nsqy = 3, ### not used, but has to provided
                nblocks = 2, ### block for parallel processing
                seed = 1 ### random number seed before starting MSE
)

### operating model
om <- FLom(stock = stk_fwd, ### stock 
           sr = sr, ### stock recruitment and precompiled residuals
           projection = mseCtrl(method = fwd_WKNSMSE, 
                                args = list(maxF = 2,
                                            ### process noise on stock.n
                                            proc_res = "fitted"
                                ))
)

### observation (error) model
oem <- FLoem(method = oem_WKNSMSE,
             observations = list(stk = stk_oem, idx = idx), 
             deviances = list(stk = FLQuants(catch.dev = catch_res), 
                              idx = idx_dev),
             args = list(idx_timing = c(0, -1),
                         catch_timing = -1,
                         use_catch_residuals = TRUE, 
                         use_idx_residuals = TRUE,
                         use_stk_oem = TRUE))
### implementation error model (banking and borrowing)
# iem <- FLiem(method = iem_WKNSMSE, 
#              args = list(BB = TRUE))

### default management
ctrl_obj <- mpCtrl(list(
  ctrl.est = mseCtrl(method = SAM_wrapper,
                     args = c(### short term forecast specifications
                       forecast = TRUE, 
                       fwd_trgt = "fsq", fwd_yrs = 1, 
                       had4_stf_def,
                       ### speeding SAM up
                       newtonsteps = 0, rel.tol = 0.001,
                       par_ini = list(sam_initial),
                       track_ini = TRUE, ### store ini for next year
                       ### SAM model specifications
                       conf = list(had4_conf_sam),
                       parallel = FALSE ### TESTING ONLY
                     )),
  ctrl.phcr = mseCtrl(method = phcr_WKNSMSE,
                      args = refpts_mse),
  ctrl.hcr = mseCtrl(method = hcr_WKNSME, args = list(option = "A")),
  ctrl.is = mseCtrl(method = is_WKNSMSE_JarworskiGrowth_intyrTACcont, 
                    args = c(hcrpars = list(refpts_mse),
                             ### for short term forecast
                             fwd_trgt = list(c("fsq", "hcr")), fwd_yrs = 2,
                             had4_stf_def#,
                             ### TAC constraint
                             #TAC_constraint = TRUE,
                             #lower = -Inf, upper = Inf,
                             #Btrigger_cond = FALSE,
                             ### banking and borrowing 
                             #BB = TRUE,
                             #BB_check_hcr = FALSE,
                             #BB_check_fc = TRUE,
                             #BB_rho = list(c(-0.1, 0.1))
                    ))#,
  #ctrl.tm = NULL
))
### additional tracking metrics
tracking_add <- c("BB_return", "BB_bank_use", "BB_bank", "BB_borrow")

### save mse objects
input <- list(om = om, oem = oem, ctrl.mp = ctrl_obj,
              genArgs = genArgs, tracking = tracking_add)
saveRDS(object = input, 
        file = paste0(input_path,"MPbase_OM",omName,"_",n,".rds"))

}

### ------------------------------------------------------------------------ ###
### Validation plots ####
### ------------------------------------------------------------------------ ###

if(F){

layout(1)
# calc spike yrs
spike.yrs<-c(2005,2009,2014)
mean.int<-mean(diff(spike.yrs))
n.spikes<-ceiling(diff(range(sim_yrs))/mean.int)
future.spike.yrs<-NULL

s <- vector(mode="list",length=10)
s<-as.list(seq(0.1,1,0.1))

spk.yrs<-sapply(s,function(sval){
  s.yrs<-NULL
  
  for (i in 1:100){
    x<-runif(n=n.spikes, min = 0, max = 1)
    #s<-0.5 # spikes are fairly regular
    y<-round(mean.int*(sval*x+1-(sval/2)))
    future.spike.yrs<-cumsum(c(max(spike.yrs),y))
    future.spike.yrs<-future.spike.yrs[future.spike.yrs>2017]
    s.yrs<-c(s.yrs,future.spike.yrs)
  } 
  return(s.yrs)
})
windows()
layout(1)
for (ind in 1:length(s)){
  if(ind==1){
    plot(jitter(spk.yrs[[ind]]),jitter(rep(s[[ind]],length(spk.yrs[[ind]]))),ylim=c(0.1,1),pch=1,ylab="s",xlab="Future spike year")
    #plot first iter#
    idx<-c(1:1000)[diff(spk.yrs[[ind]])<0]
    idx<-1:(idx[1])
    points(spk.yrs[[ind]][idx],rep(s[[ind]],length(spk.yrs[[ind]][idx])),col="green",pch=16)
    
  }else{
    points(jitter(spk.yrs[[ind]]),jitter(rep(s[[ind]],length(spk.yrs[[ind]]))))
    
    #plot first iter#
    idx<-c(1:1000)[diff(spk.yrs[[ind]])<0]
    idx<-1:(idx[1])
    points(spk.yrs[[ind]][idx],rep(s[[ind]],length(spk.yrs[[ind]][idx])),col="green",pch=16)
  }
}
savePlot(filename=paste0("output/OM_plots/OM_",omName,"_",stkname,"_spike yr variability"),type="png",res=600)

}