### ------------------------------------------------------------------------ ###
### Baseline OM for had.27.46a20 ####
### ------------------------------------------------------------------------ ###
### base on TSA assessment

rm(list=ls())

# set directories
setwd("N:\\Stock assessment\\Haddock WKNSMSE 2018\\WK_WKNSMSE_had.27.46a20")

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
# library(mse)
### load files from package mse for easier debugging
#devtools::load_all("../mse/")
library(FLash)
library(tidyr)
library(dplyr)

source("a4a_mse_WKNSMSE_funs.R")

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

omName<-"Alt1"
stkname<-"NoSh Haddock"

doPlot<-T

# set some useful index values
eqSim_yrs<-(yr_data-10):(yr_data-1)
rec.period<-2000:2017
hist_yrs<-1972:2017
ages<-0:8
sim_yrs<-(yr_data):(yr_data+n_years)

### ------------------------------------------------------------------------ ###
### 2. TSA  fit and WGNSSK 2018 results ####
### ------------------------------------------------------------------------ ###
### use input data provided in input folder (FLR data)

### recreate in the WGNSSK 2018 had assessment
load("input/haddock fit 2018 04 27.RData")
tsa_fit<-haddock_fit

# load in from WGNSSK 2018 results
load("input/had2746a20_FLStock objects_2018.Rdata")

stk <- x.pg # catch includes the 4 components. Yields are n*wt
stk.bms<-x.bms.pg
stk.ibc<-x.ibc.pg

summary(stk)

# load in indices
had_idx <- readFLIndices(file.path("input/FLR data/nor_had_ef.txt"))
had_idx[[1]]@index<-had_idx[[1]]@index*100
had_idx[[2]]@index<-had_idx[[2]]@index*100
had_idx[[1]]<-window(had_idx[[1]],start=1983)
names(had_idx)<-c("NSQ1","NSQ3")


# convert discards to include bms and ibc components
# think of this as wanted (landings) and unwanted(dis, ibc, bms) catch
# This is what happens in TSA.
# BMS is less than 0.5% of total catch (2016 onwards)
# IBC has reduced in importance over time and has been less than 1% of catch since 2004 and less than 0.12% in last 5 years
#100*(stk.bms@landings/stk@catch)
#100*(stk.ibc@landings/stk@catch)

# update numbers
# use numbers to fidn weighted mean wt at age
# recalc "discards" yield
disN<-stk@discards.n+stk.bms@landings.n+stk.ibc@landings.n
disW<-(stk@discards.wt*stk@discards.n+
         stk.bms@landings.wt*stk.bms@landings.n+
         stk.ibc@landings.wt*stk.ibc@landings.n)/disN
disW[is.na(disW)]<-0

stk@discards.n<-disN
stk@discards.wt<-disW
discards(stk) <- computeDiscards(stk)

### set units
units(stk)[1:17] <- as.list(c(rep(c("tonnes", "thousands", "kg"), 3),
                              c("tonnes", "thousands * 10000", "kg"),
                              
                              "NA", "NA", "f", "NA", "NA"))

### use input data provided in FLfse but overwrite with TSA results

# edit had_stk
had_stk<-window(stk,end=2018)
had_stk@stock[]<-NA
had_stk@stock.n[]<-NA
had_stk@harvest[]<-NA
had_stk@stock.wt[,"2018"]<-had_stk@stock.wt[,"2017"]
had_stk@m[,"2018"]<-had_stk@m[,"2017"]
had_stk@mat[,"2018"]<-had_stk@mat[,"2017"]
units(had_stk@stock.n)<-"thousands"

# then we need to overwrite had4_stk with values from had_stk becaue otherwise SAM crashes
had4_stk<-window(had4_stk,start=1972)

had4_stk@catch<-had_stk@catch
had4_stk@catch.n<-had_stk@catch.n
had4_stk@catch.wt<-had_stk@catch.wt
had4_stk@landings<-had_stk@landings
had4_stk@landings.n<-had_stk@landings.n
had4_stk@landings.wt<-had_stk@landings.wt
had4_stk@discards<-had_stk@discards
had4_stk@discards.n<-had_stk@discards.n
had4_stk@discards.wt<-had_stk@discards.wt
had4_stk@stock<-had_stk@stock
had4_stk@stock.n<-had_stk@stock.n
had4_stk@stock.wt<-had_stk@stock.wt
had4_stk@m<-had_stk@m
had4_stk@mat<-had_stk@mat
had4_stk@harvest<-had_stk@harvest

### save for later comparison
stk_orig <- had4_stk

# rename had4_stk
stk<-had4_stk


rm(list=c("x.pg","x.bms.pg","x.ibc.pg","haddock_fit","disN","disW"))

### ------------------------------------------------------------------------ ###
### 2a. SAM fit ####
### ------------------------------------------------------------------------ ###

### fit SAM to get starting values


# get conf file
had_conf_sam<-had4_conf_sam

fit <- FLR_SAM(stk = stk, idx = had_idx, conf = had_conf_sam)


is(fit)
fit

### extract model parameters and use them in the simulation as starting values
sam_initial <- sam_getpar(fit)
sam_initial$logScale <- numeric(0)


#check survey calc
sam_uncertainty <- SAM_uncertainty_had(fit = fit, n = n, print_screen = FALSE)


### ------------------------------------------------------------------------ ###
### 3. add uncertainty ####
### ------------------------------------------------------------------------ ###
### ideal approach: use variance-covariance

# first approach - just use variance.  No covariance
#rnorm(n,mean=,sd=)

### add iteration dimension
stk <- propagate(stk, n)
dim(stk)

# Get N and F uncertainty from TSA fit
wk <- simulate_stock(tsa_fit,n_sim=n,seed = 43395)

### add noise to stock
uncertainty<-vector("list",5)

names(uncertainty)<-c("stock.n","harvest","survey_catchability","survey_sd","catch_sd")

flq_template <- FLQuant(dimnames = list(age = range(stk)["min"]:range(stk)["max"], 
                                        year = range(stk)["minyear"]:range(stk)["maxyear"],iter = 1:n))

uncertainty$stock.n<-uncertainty$harvest<-catch_sd<-flq_template

tmp<-sapply(wk,function(x){return(t(x$N_sim))},simplify="array")
uncertainty$stock.n[]<-tmp[, ac(range(stk)["minyear"]:range(stk)["maxyear"]),]*10000

tmp<-sapply(wk,function(x){return(t(x$F_sim))},simplify="array")
uncertainty$harvest[]<-tmp[, ac(range(stk)["minyear"]:range(stk)["maxyear"]),]

stock.n(stk)[] <-uncertainty$stock.n
stock(stk)[] <- computeStock(stk)
### add noise to F
harvest(stk)[] <- uncertainty$harvest

rm(list=c("wk","tmp","flq_template"))

# things to check:
# compare to orig N and F
# is N negative at anypoint
# check cohorts decrease over time

range(stock.n(stk))
range(harvest(stk))


### ------------------------------------------------------------------------ ###
### 6. extend stock for MSE simulation ####
### ------------------------------------------------------------------------ ###

### stock weights, M etc in 2018 based on a three year average to enable 
### calculation of SSB. Though these will get overwritten for haddock as we use cohort growth
### although SAM/TSA estimates F in 2018, this is not reported or taken forward into
### forcasts by the WG
stk_stf2018 <- stf(window(stk, end = 2018), n_years)
### use all available data
stk_stf <- stk_stf2018


dim(stk_stf)
range(stk_stf)

### ------------------------------------------------------------------------ ###
### biological data for OM ####
### ------------------------------------------------------------------------ ###
# Baseline uses EqSim assumptions for biological parameters and  S-R relationship
# EqSim draws biological parameters from last 10 years

# Randomly sample last 10 years (EqSIm assumption)
# draw random samples
set.seed(1000001)
samp_years<-sample(eqSim_yrs,n,replace=T)
sel_samp_years<-sample(eqSim_yrs,n,replace=T)


### add random draws to iteration dimension
stf_template <- FLQuant(dimnames = list(age = range(stk_stf)["min"]:range(stk_stf)["max"], 
                                        year = sim_yrs,iter = 1:n))

### mortality
m_stf <- stf_template
m_stf[] <- (m(stk_stf)[, ac(samp_years),,,,1])
m(stk_stf)[, ac(sim_yrs)]<-m_stf

### maturity
mat_stf <- stf_template
mat_stf[] <- (mat(stk_stf)[, ac(samp_years),,,,1])
mat(stk_stf)[, ac(sim_yrs)]<-mat_stf

# weights at age - stock
stock.wt_stf <- stf_template
stock.wt_stf[] <- (stock.wt(stk_stf)[, ac(samp_years),,,,1])
stock.wt(stk_stf)[, ac(sim_yrs)]<-stock.wt_stf

# weights at age - catch
catch.wt_stf <- stf_template
catch.wt_stf[] <- (catch.wt(stk_stf)[, ac(samp_years),,,,1])
catch.wt(stk_stf)[, ac(sim_yrs)]<-catch.wt_stf

# weights at age - landings
landings.wt_stf <- stf_template
landings.wt_stf[] <- (landings.wt(stk_stf)[, ac(samp_years),,,,1])
landings.wt(stk_stf)[, ac(sim_yrs)]<-landings.wt_stf

# weights at age - discards
discards.wt_stf <- stf_template
discards.wt_stf[] <- (discards.wt(stk_stf)[, ac(samp_years),,,,1])
discards.wt(stk_stf)[, ac(sim_yrs)]<-discards.wt_stf

# selectivity
harvest_stf <- stf_template
harvest_stf[] <- (harvest(stk_stf)[, ac(sel_samp_years),,,,1])
harvest(stk_stf)[, ac(sim_yrs)]<-harvest_stf


rm(list=c("discards.wt_stf","landings.wt_stf","catch.wt_stf","stock.wt_stf","harvest_stf","m_stf","mat_stf",
          "stf_template"))

### ------------------------------------------------------------------------ ###
### 7. stock recruitment ####
### ------------------------------------------------------------------------ ###

# start with EqSim assumptions:
# years 2000 onwards

### fit hockey-stick model
### get residuals from smoothed residuals

### use only data from 2000
sr <- as.FLSR(window(stk_stf, start = min(rec.period)), model = "segreg")
### fit model individually to each iteration and suppress output to screen
suppressWarnings(. <- capture.output(sr <- fmle(sr,fixed=list(b=94000))))
### run in parallel
# library(doParallel)
# cl <- makeCluster(10)
# registerDoParallel(cl)
# sr <- fmle_parallel(sr, cl)

#do problem replicate and transfer results

pos_error <- which(is.na(params(sr)["a"]))
if(length(pos_error)){
sr_corrected <- fmle(FLCore::iter(sr, pos_error))
sr[,,,,, pos_error] <- sr_corrected[]
params(sr)[, pos_error] <- params(sr_corrected)
}

### check breakpoints
summary(params(sr)["b"])

### plot model and data
if(doPlot){
  windows()
  srP<-iter(sr,1)
  plot(srP)
  savePlot(paste0("output/had4/OM_plots/OM_",omName,"_",stkname,"_recruitment - model fit diagnostics"),type="png",res=600)
  
  as.data.frame(FLQuants(fitted = sr@fitted, rec = sr@rec, SSB = sr@ssb)) %>%
    mutate(age = NULL,
           year = ifelse(qname == "SSB", year, year)) %>%
    tidyr::spread(key = qname, value = data) %>%
    ggplot() +
    geom_point(aes(x = SSB, y = rec, group = iter), 
               alpha = 0.5, colour = "grey", shape = 1) +
    geom_line(aes(x = SSB, y = fitted, group = iter)) +
    theme_bw() + xlim(0, NA) + ylim(0, NA)
  ggsave(paste0("output/had4/OM_plots/OM_",omName,"_",stkname,"_recruitment - model fit.png"), height=4.5, width=6,dpi=600)
}

### Check extent of autocorrelation
# Not significant, so no need to account for it in this OM
#acf(window(stock.n(stk_orig)[1], start = 2000))


### Check method proposed for generating recruitment compares with past recruitment estimates
test <- as.data.frame(FLQuants(fitted = sr@fitted, rec = sr@rec, SSB = sr@ssb))
test <- mutate(test, age = NULL, year = ifelse(qname == "SSB", year, year))
test <- tidyr::spread(test, key = qname, value = data)
test <- test[complete.cases(test),]
test$res <- rep(NA, nrow(test))

# Generate residuals for future recruitments
foreach(iter_i = seq(dim(sr)[6]), .packages = "FLCore", 
        .errorhandling = "pass") %do% {
          
          set.seed(iter_i^2)
          
          ### get residuals for current iteration
          res_i <- c(FLCore::iter(residuals(sr), iter_i))
          res_i <- res_i[!is.na(res_i)]
          
          ### calculate kernel density of residuals
          density <- density(x = res_i)
          ### sample residuals
          mu <- sample(x = res_i, size = length(res_i), replace = TRUE)
          ### "smooth", i.e. sample from density distribution
          test$res[test$iter==iter_i] <- rnorm(n = length(res_i), mean = mu, sd = density$bw)
          
        }

# Generate future recruitments from past SSBs and generated residuals
test$future <- test$fitted * exp(test$res)

# 10 randomly selected iters for plotting
# should probably increase the number later
i_samp <- sample(seq(dim(sr)[6]), min(c(n,20)), replace=FALSE)

if(doPlot){
  
  # Plot past and future stock recruit pairs for selected iters
  ggplot(test[is.element(test$iter, i_samp),]) +
    geom_point(aes(x = SSB, y = rec), 
               alpha = 0.5, colour = "red", shape = 19) +
    geom_point(aes(x = SSB, y = future), 
               alpha = 0.5, colour = "black", shape = 19) +
    geom_line(aes(x = SSB, y = fitted)) +
    facet_wrap(~iter) +
    theme_bw() + xlim(0, NA) + ylim(0, NA)
  ggsave(paste0("output/had4/OM_plots/OM_",omName,"_",stkname,"_recruitment - diagnostics iter sr-pairs.png"), height=4.5, width=6,dpi=600)
  
  # Empirical cumulative distributions for the same iters
  ggplot(test[is.element(test$iter, i_samp),]) +
    stat_ecdf(aes(rec), geom = "step", colour = "red") +
    stat_ecdf(aes(future), geom = "step", colour = "black") +
    facet_wrap(~iter) +
    theme_bw() + xlim(0, NA) + ylim(0, NA)
  ggsave(paste0("output/had4/OM_plots/OM_",omName,"_",stkname,"_recruitment - diagnostics iter ECD.png"), height=4.5, width=6,dpi=600)
  
  # Combine previous two plots over iters
  # Stock recruit pairs
  ggplot(test[is.element(test$iter, i_samp),]) +
    geom_point(aes(x = SSB, y = rec), 
               alpha = 0.5, colour = "red", shape = 19) +
    geom_point(aes(x = SSB, y = future), 
               alpha = 0.5, colour = "black", shape = 19) +
    theme_bw() + xlim(0, NA) + ylim(0, NA)
  ggsave(paste0("output/had4/OM_plots/OM_",omName,"_",stkname,"_recruitment - diagnostics sr-pairs.png"), height=4.5, width=6,dpi=600)
  
  # Empirical cumulative distribution
  ggplot(test[is.element(test$iter, i_samp),]) +
    stat_ecdf(aes(rec), geom = "step", colour = "red") +
    stat_ecdf(aes(future), geom = "step", colour = "black") +
    theme_bw() + xlim(0, NA) + ylim(0, NA)
  ggsave(paste0("output/had4/OM_plots/OM_",omName,"_",stkname,"_recruitment - diagnostics ECD.png"), height=4.5, width=6,dpi=600)
  
  
  #rm(test, i_samp)
}

### generate residuals for MSE
### years with missing residuals
yrs_res <- ac(sim_yrs) #colnames(rec(sr))[which(is.na(iterMeans(rec(sr))))]
spike.yrs<-c(2005,2009,2014)

### go through iterations and create residuals
### use kernel density to create smooth distribution of residuals
### and sample from this distribution
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

if(doPlot){
  windows()
  plot(sr_res)
  savePlot(paste0("output/had4/OM_plots/OM_",omName,"_",stkname,"_recruitment - residuals"),type="png",res=600)

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
  
  savePlot(paste0("output/had4/OM_plots/OM_",omName,"_",stkname,"_recruitment - spike timing"),type="png",res=600)
  
  # plot some individual lines
  i_samp <- sample(seq(dim(sr)[6]), min(c(n,10)), replace=FALSE)
  
  ggplot(data=iter(sr_res,i_samp), aes(year, data)) + geom_line(aes(group=iter, colour=factor(iter))) + 
    facet_wrap(~iter, scales="free", nrow=5) + theme(legend.position = "none") 
  
  ggsave(paste0("output/had4/OM_plots/OM_",omName,"_",stkname,"_recruitment - worm plot.png"), height=4.5, width=6,dpi=600)
  
  
  }



### ------------------------------------------------------------------------ ###
### stf for 2018: assume catch advice is taken ####
### ------------------------------------------------------------------------ ###
c2018 <- 48990 # advised catch for 2018
ctrl <- fwdControl(data.frame(year = 2018, quantity = "catch", 
                              val = c2018))

### project forward for intermediate year (2018)
stk_int <- stk_stf
stk_int[] <- fwd(stk_stf, ctrl = ctrl, sr = sr, sr.residuals = sr_res,
                 sr.residuals.mult = T, maxF = 5)[]

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
### biological data for OEM ####
### ------------------------------------------------------------------------ ###

### base on OM stock
stk_oem <- stk_fwd

### use means of sampled values for projection period
catch.wt(stk_oem)[, ac(sim_yrs)] <- 
  yearMeans(catch.wt(stk_oem)[, ac(eqSim_yrs)])
landings.wt(stk_oem)[, ac(sim_yrs)] <- 
  yearMeans(landings.wt(stk_oem)[, ac(eqSim_yrs)])
discards.wt(stk_oem)[, ac(sim_yrs)] <- 
  yearMeans(discards.wt(stk_oem)[, ac(eqSim_yrs)])
stock.wt(stk_oem)[, ac(sim_yrs)] <- 
  yearMeans(stock.wt(stk_oem)[, ac(eqSim_yrs)])
m(stk_oem)[, ac(sim_yrs)] <- yearMeans(m(stk_oem)[, ac(eqSim_yrs)])
mat(stk_oem)[, ac(sim_yrs)] <- yearMeans(mat(stk_oem)[, ac(eqSim_yrs)])

### remove stock assessment results
stock.n(stk_oem)[] <- stock(stk_oem)[] <- harvest(stk_oem)[] <- NA


### ------------------------------------------------------------------------ ###
### indices TSA ####
### ------------------------------------------------------------------------ ###
### use real FLIndices object as template (included in FLfse)

idx<-had_idx

### extend for simulation period
idx <- window(idx, end = max(sim_yrs))
### add iterations
idx <- lapply(idx, propagate, n)

## get catchabilities and survey noise
#make templates
flq_template1 <- FLQuant(dimnames = list(age = range(idx[[1]])["min"]:range(idx[[1]])["max"], 
                                         iter = 1:n))
flq_template2 <- FLQuant(dimnames = list(age = range(idx[[2]])["min"]:range(idx[[2]])["max"], 
                                         iter = 1:n))
uncertainty$survey_sd<-uncertainty$survey_catchability<-list("NSQ1"=flq_template1,"NSQ3"=flq_template2)

wk_NSQ1 <- simulate_survey(tsa_fit, "NSQ1", seed = 72,n_sim=n)
wk_NSQ3 <- simulate_survey(tsa_fit, "NSQ3", seed = 55,n_sim=n)

uncertainty$survey_sd[[1]][]<-t(wk_NSQ1$cv_at_age)
uncertainty$survey_sd[[2]][]<-t(wk_NSQ3$cv_at_age)
uncertainty$survey_catchability[[1]][]<-t(wk_NSQ1$selection)/1000
uncertainty$survey_catchability[[2]][]<-t(wk_NSQ3$selection)/1000


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
set.seed(1000002)
for (idx_i in seq_along(idx_dev)) {
  ### insert sd
  idx_dev[[idx_i]][] <- uncertainty$survey_sd[[idx_i]]
  ### noise
  idx_dev[[idx_i]][] <- stats::rnorm(n = length(idx_dev[[idx_i]]),
                                     mean = 1, sd = idx_dev[[idx_i]])
} #TSA is cv not sd

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

hadTSA_idx_dev<-idx_dev
hadTSA_idx<-idx

### ------------------------------------------------------------------------ ###
### indices SAM ####
### ------------------------------------------------------------------------ ###
### use real FLIndices object as template (included in FLfse)
idx <- had_idx
### extend for simulation period
idx <- window(idx, end = yr_data + n_years)
### add iterations
idx <- lapply(idx, propagate, n)

### insert catchability
for (idx_i in seq_along(idx)) {
  
  ### set catchability for projection
  index.q(idx[[idx_i]])[] <- sam_uncertainty$survey_catchability[[idx_i]]
  
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
  idx_dev[[idx_i]][] <- sam_uncertainty$survey_sd[[idx_i]]
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

hadSAM_idx_dev<-idx_dev
hadSAM_idx<-idx

### ------------------------------------------------------------------------ ###
### catch noise ####
### ------------------------------------------------------------------------ ###
### take estimates from sam: uncertainty$catch_sd is "logSdLogObs"
### assume catch observed by SAM in projection is log-normally distributed
### around operating model catch

## get catchabilities and survey noise
#make templates
flq_template <- FLQuant(dimnames = list(age = range(stk)["min"]:range(stk)["max"], 
                                         iter = 1:n))

uncertainty$catch_sd<-list("landings"=flq_template,"discards"=flq_template)

wk_cvmult <- list(
  landings = c(41, 14, 11, 12, 15, 18, 30, 31) / 11, 
  discards = c(48, 40, 24, 36, 43, 57) / 24
)

wk_cv <- simulate_catch_cv(tsa_fit, seed = 33, cvmult = wk_cvmult,n_sim=n)

tmp<-cbind("age 0"=0,wk_cv$landings)
uncertainty$catch_sd[[1]][]<-t(tmp)
tmp<-cbind(wk_cv$discards,"age 6"=0,"age 7"=0,"age 8"=0)
uncertainty$catch_sd[[2]][]<-t(tmp)

### create noise for landings
set.seed(1000003)
#first need to make up some landings.n and discards for 2018
#assume landings.n in 2018 are same split as 2017
landings.n(stk_fwd)[,"2018"]<-catch.n(stk_fwd)[,"2018"]*(landings.n(stk_fwd)[,"2017"]/catch.n(stk_fwd)[,"2017"])
discards.n(stk_fwd)[,"2018"]<-catch.n(stk_fwd)[,"2018"]*(1 - landings.n(stk_fwd)[,"2017"]/catch.n(stk_fwd)[,"2017"])

lan_res <- landings.n(stk_fwd) %=% 0 ### template FLQuant
lan_res[] <- stats::rnorm(n = length(lan_res), mean = 1, 
                            sd = uncertainty$catch_sd[[1]])

lan1<-lan_res*landings.n(stk_fwd)

### create noise for discards
dis_res <- discards.n(stk_fwd) %=% 0 ### template FLQuant
dis_res[] <- stats::rnorm(n = length(dis_res), mean = 1, 
                          sd = uncertainty$catch_sd[[2]])

dis1<-dis_res*discards.n(stk_fwd)

#combine to calc catch CV
catch1<-lan1+dis1
tmp1<-catch.n(stk_fwd)
tmp1[is.na(tmp1)]<-1
catch_res<-catch1/tmp1

#put landings.n and discards.n in 2018 back to NA
landings.n(stk_fwd)[,"2018"]<-NA
discards.n(stk_fwd)[,"2018"]<-NA


# makes historic = 1
catch_res[, dimnames(catch_res)$year <= 2017] <- 1


### ------------------------------------------------------------------------ ###
### save OM ####
### ------------------------------------------------------------------------ ###

# choose index
idx<-hadTSA_idx
idx_dev<-hadTSA_idx_dev

### path
input_path <- paste0("input/had4/", n, "_", n_years, "/")
dir.create(input_path)

### stock
saveRDS(stk_fwd, file = paste0(input_path, omName,"_stk.rds"))
### stock recruitment
saveRDS(sr, file = paste0(input_path, omName,"_sr.rds"))
### recruitment residuals
saveRDS(sr_res, file = paste0(input_path, omName,"_sr_res.rds"))
### surveys
saveRDS(idx, file = paste0(input_path, omName,"_idx.rds"))
saveRDS(idx_dev, file = paste0(input_path, omName,"_idx_dev.rds"))
### catch noise
saveRDS(catch_res, file = paste0(input_path, omName,"_catch_res.rds"))
### process error
#saveRDS(proc_res, file = paste0(input_path, "proc_res.rds"))
### observed stock
saveRDS(stk_oem, file = paste0(input_path, omName,"_stk_oem.rds"))

### sam initial parameters
saveRDS(sam_initial, file = paste0(input_path, omName,"_sam_initial.rds"))
### sam configuration
saveRDS(had_conf_sam, file = paste0(input_path, omName,"_had_conf_sam"))


save.image(file=paste0("input/had4/image_OM",omName,"_",n,".RData"))
#load(file=paste0("input/had4/image_OM",omName,"_",n,".RData"))

# 
# ### reference points
# refpts_mse <- list(Btrigger = 132000,
#                    Ftrgt = 0.194,
#                    Fpa = 0.274,
#                    Bpa = 94000)
# 
# ### some specifications for short term forecast with SAM
# had_stf_def <- list(fwd_yrs_average = -3:0, #3 year average
#                     fwd_yrs_rec_start = 2000,
#                     fwd_yrs_sel = -1, #same as last year
#                     fwd_yrs_lf_remove = -2:-1,
#                     fwd_splitLD = FALSE)
# 
# ### some arguments (passed to mp())
# genArgs <- list(fy = dims(stk_fwd)$maxyear, ### final simulation year
#                 y0 = dims(stk_fwd)$minyear, ### first data year
#                 iy = yr_data, ### first simulation (intermediate) year
#                 nsqy = 3, ### not used, but has to provided
#                 nblocks = 1, ### block for parallel processing
#                 seed = 1 ### random number seed before starting MSE
# )
# 
# ### operating model
# om <- FLom(stock = stk_fwd, ### stock 
#            sr = sr, ### stock recruitment and precompiled residuals
#            projection = mseCtrl(method = fwd_WKNSMSE, 
#                                 args = list(maxF = 2,
#                                             ### process noise on stock.n
#                                             proc_res = "fitted"
#                                 ))
# )
# 
# ### observation (error) model
# oem <- FLoem(method = oem_WKNSMSE,
#              observations = list(stk = stk_oem, idx = idx), 
#              deviances = list(stk = FLQuants(catch.dev = catch_res), 
#                               idx = idx_dev),
#              args = list(idx_timing = c(0, -1),
#                          catch_timing = -1,
#                          use_catch_residuals = TRUE, 
#                          use_idx_residuals = TRUE,
#                          use_stk_oem = TRUE))
# ### implementation error model (banking and borrowing)
# # iem <- FLiem(method = iem_WKNSMSE, 
# #              args = list(BB = TRUE))
# 
# ### default management
# ctrl_obj <- mpCtrl(list(
#   ctrl.est = mseCtrl(method = SAM_wrapper,
#                      args = c(### short term forecast specifications
#                        forecast = TRUE, 
#                        fwd_trgt = "fsq", fwd_yrs = 1, 
#                        had_stf_def,
#                        ### speeding SAM up
#                        newtonsteps = 0, rel.tol = 0.001,
#                        par_ini = list(sam_initial),
#                        track_ini = TRUE, ### store ini for next year
#                        ### SAM model specifications
#                        conf = list(had_conf_sam),
#                        parallel = FALSE ### TESTING ONLY
#                      )),
#   ctrl.phcr = mseCtrl(method = phcr_WKNSMSE,
#                       args = refpts_mse),
#   ctrl.hcr = mseCtrl(method = hcr_WKNSME, args = list(option = "A")),
#   ctrl.is = mseCtrl(method = is_WKNSMSE, 
#                     args = c(hcrpars = list(refpts_mse),
#                              ### for short term forecast
#                              fwd_trgt = c("fsq", "hcr"), fwd_yrs = 2,
#                              had_stf_def#,
#                              ### TAC constraint
#                              #TAC_constraint = TRUE,
#                              #lower = -Inf, upper = Inf,
#                              #Btrigger_cond = FALSE,
#                              ### banking and borrowing 
#                              #BB = TRUE,
#                              #BB_conditional = TRUE,
#                              #BB_rho = list(c(-0.1, 0.1))
#                     ))#,
#   #ctrl.tm = NULL
# ))
# ### additional tracking metrics
# tracking_add <- c("BB_return", "BB_bank_use", "BB_bank", "BB_borrow")
# 
# 
# ### save mse objects
# input <- list(om = om, oem = oem, ctrl.mp = ctrl_obj,
#               genArgs = genArgs, tracking = tracking_add)
# saveRDS(object = input, 
#         file = paste0(input_path, "base_run.rds"))
# # input <- readRDS(paste0(input_path, "/base_run.rds"))
# 
# ### ------------------------------------------------------------------------ ###
# ### run MSE ####
# ### ------------------------------------------------------------------------ ###
# 
# 
# ### run MSE
# ### WARNING: takes a while...
# ### check normal execution
# res1 <- mp(om = input$om,
#            oem = input$oem,
#            #iem = iem,
#            ctrl.mp = input$ctrl.mp,
#            genArgs = input$genArgs,
#            tracking = input$tracking)

