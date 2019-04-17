### ------------------------------------------------------------------------ ###
### Baseline OM for had.27.46a20 ####
### ------------------------------------------------------------------------ ###
### base on SAM assessment

rm(list=ls())

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

omName<-"Baseline"
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
# think of this as wanted (landings) and unwanted(dis, ibc, bms) catch
# This is what happens in TSA.
# BMS is less than 0.5% of total catch (2016 onwards)
# IBC has reduced in importance over time and has been less than 1% of catch since 2004 and less than 0.12% in last 5 years
#100*(stk.bms@landings/stk@catch)
#100*(stk.ibc@landings/stk@catch)

if(verbose){
  ggplot(data=FLQuants(catch=catch(x.pg), landings=landings(x.pg), discards=discards(x.pg),
                       BMS=landings(x.bms.pg),IBC=landings(x.ibc.pg)), aes(year, data)) + 
    geom_line(aes(group=qname,colour=factor(qname)))+ theme(legend.position = c(0.8,0.7))+labs(color="Yield",y="Yield (t)") 
  ggsave(paste0("output/OM_plots/OM_",omName,"_",stkname,"_data - catch components.png"), height=4.5, width=6,dpi=600)
  
  pct.bms<-100*(x.bms.pg@landings/x.pg@catch)
  pct.ibc<-100*(x.ibc.pg@landings/x.pg@catch)
  ggplot(data=FLQuants(BMS=pct.bms,IBC=pct.ibc),aes(year,data))+
    geom_line(aes(group=qname,colour=factor(qname)))+ theme(legend.position = c(0.85,0.85))+labs(color="% of catch",y="Yield (t)") 
  ggsave(paste0("output/OM_plots/OM_",omName,"_",stkname,"_data - catch comp percent of catch 1.png"), height=4.5, width=6,dpi=600)
  
  ggplot(data=FLQuants(BMS=pct.bms[,ac(rec.period)],IBC=pct.ibc[,ac(rec.period)]),aes(year,data))+
    geom_line(aes(group=qname,colour=factor(qname)))+ theme(legend.position = c(0.85,0.85))+labs(color="% of catch",y="Yield (t)") 
  ggsave(paste0("output/OM_plots/OM_",omName,"_",stkname,"_data - catch comp percent of catch 2.png"), height=4.5, width=6,dpi=600)
}

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

if (isTRUE(verbose)) {
  is(fit)
  fit
  plot(fit)
  savePlot(paste0("output/OM_plots/OM_",omName,"_",stkname,"_SAM fit - stock summary"),res=600,type="png")
  
  plot(c(1972,2018),c(0,700000),type="n",ylab="SSB (t)",xlab="")
  ssbplot(fit,add=T)
    title("SSB")
  lines(hist_yrs,rep(132000,length(hist_yrs)),lty=2) # Bpa=Btrigger
  lines(hist_yrs,rep(94000,length(hist_yrs)),lty=3) #Blim
legend("topright",inset=0.02,legend=c("Bpa=Btrigger","Blim"),col="black",lty=c(2,3),cex=0.8)
savePlot(paste0("output/OM_plots/OM_",omName,"_",stkname,"_SAM fit - SSB"),res=600,type="png")

  #add Fmsy, Fpa, Flim
  # add labels and title
  plot(c(1972,2018),c(0,1.5),type="n",ylab="Fbar (2-4)",xlab="")
  fbarplot(fit,add=T,partial=F)
  title("Fbar")
  lines(hist_yrs,rep(0.194,length(hist_yrs)),lty=4) # Fmsy
  lines(hist_yrs,rep(0.274,length(hist_yrs)),lty=2) # Fpa
  lines(hist_yrs,rep(0.384,length(hist_yrs)),lty=3) # Flim
  legend("topright",inset=0.02,legend=c("Flim","Fpa","Fmsy"),col="black",lty=c(3,2,4),cex=0.8)
  savePlot(paste0("output/OM_plots/OM_",omName,"_",stkname,"_SAM fit - Fbar"),res=600,type="png")
  
  plot(c(1972,2018),c(0,5000000),type="n",ylab="TSB (t)",xlab="")
  tsbplot(fit,add=T)
  title("TSB")
  savePlot(paste0("output/OM_plots/OM_",omName,"_",stkname,"_SAM fit - TSB"),res=600,type="png")

plot(c(1972,2018),c(0,600000),type="n",ylab="Catch (t)",xlab="")
  catchplot(fit,add=T)
  title("Catch")
  savePlot(paste0("output/OM_plots/OM_",omName,"_",stkname,"_SAM fit - Catch"),res=600,type="png")
  
  #add title
  recplot(fit)
  plot(c(1972,2018),c(0,104000000),type="n",ylab="Recruitment (1000)",xlab="")
  recplot(fit,add=T)
  title("Recruitment")
  savePlot(paste0("output/OM_plots/OM_",omName,"_",stkname,"_SAM fit - Recruitment"),res=600,type="png")
  
  fitplot(fit,fleet=1) #catch
  savePlot(paste0("output/OM_plots/OM_",omName,"_",stkname,"_SAM fit - fit to catch"),res=600,type="png")
  
  fitplot(fit,fleet=2) #NSQ1
  savePlot(paste0("output/OM_plots/OM_",omName,"_",stkname,"_SAM fit - fit to NSQ1"),res=600,type="png")
  
  fitplot(fit,fleet=3) #NSQ3
  savePlot(paste0("output/OM_plots/OM_",omName,"_",stkname,"_SAM fit - fit to NSQ3"),res=600,type="png")
  
  res <- residuals(fit)
  plot(res)
  savePlot(paste0("output/OM_plots/OM_",omName,"_",stkname,"_SAM fit - residuals"),res=600,type="png")
  
  resp <- procres(fit)
  plot(resp)
  savePlot(paste0("output/OM_plots/OM_",omName,"_",stkname,"_SAM fit - process error residuals"),res=600,type="png")
  
  #lo <- leaveout(fit)
  #plot(lo)
}

### ------------------------------------------------------------------------ ###
### 4. create FLStock ####
### ------------------------------------------------------------------------ ###
### create template with 1 iteration

### the results (stock numbers & harvest) are used from the real WGNSSK fit
stk <- SAM2FLStock(object = fit, stk = had4_stk)
if (isTRUE(verbose)) summary(stk)

### set units
units(stk)[1:17] <- as.list(c(rep(c("t", "1000", "kg"), 4),
                              "", "", "f", "", ""))

if (isTRUE(verbose)) plot(stk)

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

if (isTRUE(verbose)) plot(stk, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))

### maximum observed F
max(fbar(stk))
# 1.652771 in year 1991 
max(harvest(stk))
# 2.270628 in year 1972 for age 7 

### get estimated catch numbers
catch_n <- uncertainty$catch_n

#plots
if(verbose){
  ggplot(data=stk@stock.n, aes(year, data)) + geom_line(aes(group=age, colour=factor(age))) + 
    geom_flquantiles(probs = c(0.05,0.5,0.95),aes(fill=factor(age)),alpha=0.8)+
    facet_wrap(~age, scales="free_y", nrow=3) + theme(legend.position = "none") + ggtitle("stock.n")
  ggsave(paste0("output/OM_plots/OM_",omName,"_",stkname,"_parameter uncertainty - stock n.png"), height=4.5, width=6,dpi=600)
  
  ggplot(data=stk@harvest, aes(year, data)) + geom_line(aes(group=age, colour=factor(age))) + 
    geom_flquantiles(probs = c(0.05,0.5,0.95),aes(fill=factor(age)),alpha=0.8)+
    facet_wrap(~age, scales="free_y", nrow=3) + theme(legend.position = "none") + ggtitle("F")
  ggsave(paste0("output/OM_plots/OM_",omName,"_",stkname,"_parameter uncertainty - harvest.png"), height=4.5, width=6,dpi=500)
  
  plot(stk, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
  ggsave(paste0("output/OM_plots/OM_",omName,"_",stkname,"_parameter uncertainty - stock summary.png"), height=4.5, width=6,dpi=600)
}

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

if (isTRUE(verbose)) plot(stk_stf)


if(verbose){
  #plot
  plot(stk_stf@m)
  ggsave(paste0("output/OM_plots/OM_",omName,"_",stkname,"_biological parameters - mortality.png"), height=4.5, width=6,dpi=600)
  plot(stk_stf@mat)
  ggsave(paste0("output/OM_plots/OM_",omName,"_",stkname,"_biological parameters - maturity.png"), height=4.5, width=6,dpi=600)
  plot(stk_stf@stock.wt)
  ggsave(paste0("output/OM_plots/OM_",omName,"_",stkname,"_biological parameters - stockwt.png"), height=4.5, width=6,dpi=600)
  plot(stk_stf@landings.wt)
  ggsave(paste0("output/OM_plots/OM_",omName,"_",stkname,"_biological parameters - landingswt.png"), height=4.5, width=6,dpi=600)
  plot(stk_stf@discards.wt)
  ggsave(paste0("output/OM_plots/OM_",omName,"_",stkname,"_biological parameters - discardswt.png"), height=4.5, width=6,dpi=600)
  plot(stk_stf@catch.wt)
  ggsave(paste0("output/OM_plots/OM_",omName,"_",stkname,"_biological parameters - catchwt.png"), height=4.5, width=6,dpi=600)
  plot(stk_stf@harvest)
  ggsave(paste0("output/OM_plots/OM_",omName,"_",stkname,"_biological parameters - selectivity.png"), height=4.5, width=6,dpi=600)
}

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

if (isTRUE(verbose)) {
  
  # plot(sr)
  sr_orig <- as.FLSR(window(stk_orig, start = min(rec.period)), model = "segreg")
  suppressWarnings(. <- capture.output(sr_orig <- fmle(sr_orig)))
  #srP<-iter(sr,1)
  plot(sr_orig)
  savePlot(paste0("output/OM_plots/OM_",omName,"_",stkname,"_recruitment - model fit diagnostics"),type="png",res=600)
  
  ### check breakpoints
  summary(params(sr_fit)["b"])
  
  ### plot model and data
  as.data.frame(FLQuants(fitted = sr_fit@fitted, rec = sr_fit@rec, SSB = sr_fit@ssb)) %>%
    mutate(age = NULL,
           year = ifelse(qname == "SSB", year+1, year)) %>%
    tidyr::spread(key = qname, value = data) %>%
    ggplot() +
    geom_point(aes(x = SSB, y = rec, group = iter), 
               alpha = 0.5, colour = "grey", shape = 1) +
    geom_line(aes(x = SSB, y = fitted, group = iter)) +
    theme_bw() + xlim(0, NA) + ylim(0, NA)
  ggsave(paste0("output/OM_plots/OM_",omName,"_",stkname,"_recruitment - model fit.png"), height=4.5, width=6,dpi=600)
  
  ### Check extent of autocorrelation
  # Not significant, so no need to account for it in this OM
  acf(window(stock.n(stk_orig)[1], start = min(rec.period)))
  
  ### Check method proposed for generating recruitment compares with past recruitment estimates
  test <- as.data.frame(FLQuants(fitted = sr_fit@fitted, rec = sr_fit@rec, SSB = sr_fit@ssb))
  test <- mutate(test, age = NULL, year = ifelse(qname == "SSB", year+1, year))
  test <- tidyr::spread(test, key = qname, value = data)
  test <- test[complete.cases(test),]
  test$res <- rep(NA, nrow(test))
  
  # Generate residuals for future recruitments
  foreach(iter_i = seq(dim(sr)[6]), .packages = "FLCore", 
          .errorhandling = "pass") %do% {
            
            set.seed(iter_i^2)
            
            ### get residuals for current iteration
            res_i <- c(FLCore::iter(residuals(sr_fit), iter_i))
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
  i_samp <- sample(seq(dim(sr_fit)[6]), 20, replace=FALSE)
  
  # Plot past and future stock recruit pairs for selected iters
  ggplot(test[is.element(test$iter, i_samp),]) +
    geom_point(aes(x = SSB, y = rec), 
               alpha = 0.5, colour = "red", shape = 19) +
    geom_point(aes(x = SSB, y = future), 
               alpha = 0.5, colour = "black", shape = 19) +
    geom_line(aes(x = SSB, y = fitted)) +
    facet_wrap(~iter) +
    theme_bw() + xlim(0, NA) + ylim(0, NA)
  ggsave(paste0("output/OM_plots/OM_",omName,"_",stkname,"_recruitment - diagnostics iter sr-pairs.png"), height=4.5, width=6,dpi=600)
  
  
  # Empirical cumulative distributions for the same iters
  ggplot(test[is.element(test$iter, i_samp),]) +
    stat_ecdf(aes(rec), geom = "step", colour = "red") +
    stat_ecdf(aes(future), geom = "step", colour = "black") +
    facet_wrap(~iter) +
    theme_bw() + xlim(0, NA) + ylim(0, NA)
  ggsave(paste0("output/OM_plots/OM_",omName,"_",stkname,"_recruitment - diagnostics iter ECD.png"), height=4.5, width=6,dpi=600)
  
  # Combine previous two plots over iters
  # Stock recruit pairs
  i_samp <- sample(seq(dim(sr_fit)[6]), 1000, replace=FALSE)
  
  ggplot(test[is.element(test$iter, i_samp),]) +
    geom_point(aes(x = SSB, y = rec), 
               alpha = 0.5, colour = "red", shape = 19) +
    geom_point(aes(x = SSB, y = future), 
               alpha = 0.5, colour = "black", shape = 19) +
    theme_bw() + xlim(0, NA) + ylim(0, NA)
  ggsave(paste0("output/OM_plots/OM_",omName,"_",stkname,"_recruitment - diagnostics sr-pairs.png"), height=4.5, width=6,dpi=600)
  
  # Empirical cumulative distribution
  ggplot(test[is.element(test$iter, i_samp),]) +
    stat_ecdf(aes(rec), geom = "step", colour = "red") +
    stat_ecdf(aes(future), geom = "step", colour = "black") +
    theme_bw() + xlim(0, NA) + ylim(0, NA)
  ggsave(paste0("output/OM_plots/OM_",omName,"_",stkname,"_recruitment - diagnostics ECD.png"), height=4.5, width=6,dpi=600)
  
  rm(test, i_samp)
}

### generate residuals for MSE
### years with missing residuals
yrs_res <- dimnames(sr)$year[which(is.na(iterMeans(rec(sr))))]
#yrs_res <- colnames(rec(sr))[which(is.na(iterMeans(rec(sr))))]

### go through iterations and create residuals
### use kernel density to create smooth distribution of residuals
### and sample from this distribution
res_new <- foreach(iter_i = seq(dim(sr)[6]), .packages = "FLCore", 
                   .errorhandling = "pass") %do% {
                     
                     set.seed(iter_i)
                     
                     ### get residuals for current iteration
                     res_i <- c(FLCore::iter(residuals(sr), iter_i))
                     res_i <- res_i[!is.na(res_i)]
                     
                     ### calculate kernel density of residuals
                     density <- density(x = res_i)
                     ### sample residuals
                     mu <- sample(x = res_i, size = length(yrs_res), replace = TRUE)
                     ### "smooth", i.e. sample from density distribution
                     res_new <- rnorm(n = length(yrs_res), mean = mu, sd = density$bw)
                     
                     return(res_new)
                     
                   }
summary(exp(unlist(res_new)))
### insert into model
residuals(sr)[, yrs_res] <- unlist(res_new)
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

if (isTRUE(verbose)){
  plot(proc_res)
  ggsave(paste0("output/OM_plots/OM_",omName,"_",stkname,"_process error.png"), height=7, width=11,dpi=300)
}


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

if(verbose){
  plot(idx_dev[[1]])+ggtitle("NSQ1")
  ggsave(paste0("output/OM_plots/OM_",omName,"_",stkname,"_survey - NSQ1 deviations.png"), height=4.5, width=6,dpi=600)
  
  plot(idx_dev[[2]])+ggtitle("NSQ3")
  ggsave(paste0("output/OM_plots/OM_",omName,"_",stkname,"_survey - NSQ3 deviations.png"), height=4.5, width=6,dpi=600)
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


if (isTRUE(verbose)) {
  
  ### compare simulated to original survey(s)
  as.data.frame(FLQuants(had4_q1 = index(had4_idx$NSQ1), 
                         had4_q3 = index(had4_idx$NSQ3),
                         sim_q1 = (index(idx$NSQ1)),
                         sim_q3 = (index(idx$NSQ3))
  )) %>%
    mutate(survey = ifelse(grepl(x = qname, pattern = "*_q1$"), "Q1", "Q3"),
           source = ifelse(grepl(x = qname, pattern = "^sim*"), "sim", "data")) %>%
    filter(year <= 2018) %>%
    ggplot(aes(x = year, y = data, colour = source)) +
    facet_grid(paste("age", age) ~ paste("IBTS", survey), scales = "free_y") +
    stat_summary(fun.y = quantile, fun.args = 0.25, geom = "line",
                 alpha = 0.5) +
    stat_summary(fun.y = quantile, fun.args = 0.75, geom = "line",
                 alpha = 0.5) +
    stat_summary(fun.y = median, geom = "line") +
    theme_bw()
  ggsave(paste0("output/OM_plots/OM_",omName,"_",stkname,"_survey_sim_survey.png"), height=4.5, width=6,dpi=600)
  
  plot(idx_dev[[1]])+ggtitle("NSQ1")
  ggsave(paste0("output/OM_plots/OM_",omName,"_",stkname,"_survey - NSQ1 deviations with historical corrections.png"), height=4.5, width=6,dpi=600)
  
  plot(idx_dev[[2]])+ggtitle("NSQ3")
  ggsave(paste0("output/OM_plots/OM_",omName,"_",stkname,"_survey - NSQ3 deviations with historical corrections.png"), height=4.5, width=6,dpi=600)
  
}

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

if (isTRUE(verbose)) {
  plot(catch_res, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
  ggsave(paste0("output/OM_plots/OM_",omName,"_",stkname,"_catch - deviations.png"), height=4.5, width=6,dpi=600)
}


### ------------------------------------------------------------------------ ###
### save OM ####
### ------------------------------------------------------------------------ ###

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
                       fwd_trgt = "intyrTACcont", fwd_yrs = 1, 
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
                             fwd_trgt = list(c("intyrTACcont", "hcr")), fwd_yrs = 2,
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
  library(pracma)
  
  windows()
  # biological params
  
  #### check for trends 
  # First check for trends in last 10 years
  
  samp_period<-(yr_data-9):(yr_data)
  
  ggplot(data=stk_orig@stock.wt[,ac(samp_period)], aes(year, data)) + geom_line(aes(group=age, colour=factor(age))) + 
    facet_wrap(~age, scales="free_y", nrow=3) + theme(legend.position = "none") + ggtitle("stock.wt")
  ggsave(paste0("output/OM_plots/OM_",omName,"_",stkname,"_biological parameters - trend check - stock wt.png"), height=4.5, width=6,dpi=600)
  
  ggplot(data=stk_orig@landings.wt[,ac(samp_period)], aes(year, data)) + geom_line(aes(group=age, colour=factor(age))) + 
    facet_wrap(~age, scales="free_y", nrow=3) + theme(legend.position = "none") + ggtitle("landings.wt")
  ggsave(paste0("output/OM_plots/OM_",omName,"_",stkname,"_biological parameters - trend check - landings wt.png"), height=4.5, width=6,dpi=600)
  
  ggplot(data=stk_orig@discards.wt[,ac(samp_period)], aes(year, data)) + geom_line(aes(group=age, colour=factor(age))) + 
    facet_wrap(~age, scales="free_y", nrow=3) + theme(legend.position = "none") + ggtitle("discards.wt")
  ggsave(paste0("output/OM_plots/OM_",omName,"_",stkname,"_biological parameters - trend check - discards wt.png"), height=4.5, width=6,dpi=600)
  
  ggplot(data=stk_orig@catch.wt[,ac(samp_period)], aes(year, data)) + geom_line(aes(group=age, colour=factor(age))) + 
    facet_wrap(~age, scales="free_y", nrow=3) + theme(legend.position = "none") + ggtitle("catch.wt")
  ggsave(paste0("output/OM_plots/OM_",omName,"_",stkname,"_biological parameters - trend check - catch wt.png"), height=4.5, width=6,dpi=600)
  
  ggplot(data=stk_orig@m[,ac(samp_period)], aes(year, data)) + geom_line(aes(group=age, colour=factor(age))) + 
    facet_wrap(~age, scales="free_y", nrow=3) + theme(legend.position = "none") + ggtitle("mortality")
  ggsave(paste0("output/OM_plots/OM_",omName,"_",stkname,"_biological parameters - trend check - mortality.png"), height=4.5, width=6,dpi=600)
  
  ggplot(data=stk_orig@harvest[,ac(samp_period)], aes(year, data)) + geom_line(aes(group=age, colour=factor(age))) + 
    facet_wrap(~age, scales="free_y", nrow=3) + theme(legend.position = "none") + ggtitle("selectivity")
  ggsave(paste0("output/OM_plots/OM_",omName,"_",stkname,"_biological parameters - trend check - selectivity.png"), height=4.5, width=6,dpi=600)
  
  
  #### check for autocorrelation 
  
  for (var in c("stock.wt","catch.wt","landings.wt","discards.wt","m","harvest")){
    varName<-var
    ages<-range(stk_orig)["min"]:range(stk_orig)["max"]
    
    layout(matrix(1:length(ages),nrow=3))
    for (i in 1:length(ages)){
      dat<-eval(parse(text=paste0(varName,"(stk_orig[ac(ages[i]),,,,,1])")))[,ac(eqSim_yrs)]
      tmp<-detrend(matrix(dat))
      # tmp[tmp==0]<-NA
      acf(tmp,main=paste0(varName, " ",ages[i]))
      savePlot(paste0("output/OM_plots/OM_",omName,"_",stkname,"_biological parameters - AR check - ",gsub("[.]","",varName),".png"), type="png",res=600)
    }
  }
  
  
  #### compare distirubtions
  layout(mat=matrix(c(1,2)))
  
  for (var in c("stock.wt","catch.wt","landings.wt","discards.wt","m","mat","harvest")){
    varName<-var
    dat<-eval(parse(text=paste0(varName,"(stk_stf)")))
    rng<-c(min(floor(range(dat))),max(ceiling(range(dat))))
    
    h1<-hist(dat[,ac(eqSim_yrs)],plot=F,breaks=seq(rng[1],rng[2],0.1))
    h2<-hist(dat[,ac(sim_yrs)],plot=F,breaks=seq(rng[1],rng[2],0.1))
    barplot(h1$counts,names.arg=h1$mids,main=paste0(varName," ",min(eqSim_yrs),"-",max(eqSim_yrs)))
    barplot(h2$counts,names.arg=h2$mids,main=paste0(varName," ",min(sim_yrs),"-",max(sim_yrs)))
    savePlot(paste0("output/OM_plots/OM_",omName,"_",stkname,"_biological parameters dist - ",gsub("[.]","",varName),".png"), type="png",res=600)
  }
  
  
  #### check selelctivities
  layout(1)
  sel<-matrix(harvest(stk_stf)[, ac(eqSim_yrs),,,,1],ncol=length(eqSim_yrs),dimnames = list(ages,eqSim_yrs))
  sel2<-matrix(harvest(stk_stf)[, ac(c(eqSim_yrs-10,eqSim_yrs)),,,,1],ncol=2*length(eqSim_yrs),dimnames = list(ages,c(eqSim_yrs-10,eqSim_yrs)))
  sel<-sel/colSums(sel)
  sel2<-sel2/colSums(sel2)
  mSel10<-apply(sel,1,mean)
  mSel3<-apply(sel[,ac(2015:2017)],1,mean)
  mSel5<-apply(sel[,ac(2013:2017)],1,mean)
  mSelp5<-apply(sel[,ac(2008:2012)],1,mean)
  mSel20<-apply(sel2,1,mean)
  
 
  # how do years compare to means
  plot(ages,mSel10,type="l",main="Selection curves 2008-2017",ylab="Selectivity",xlab="age",col="grey20",ylim=c(0,0.5),
       lwd=3)
  cols<-rainbow(length(eqSim_yrs))
  for (i in 1:length(eqSim_yrs)){
    lines(ages,sel[,ac(eqSim_yrs[i])],col=cols[i])
  }
  lines(ages,mSel5,col="grey60",lwd=3)
  lines(ages,mSel3,col="grey80",lwd=3)
  lines(ages,mSelp5,col="grey40",lwd=3)
  lines(ages,mSelp10,col="grey20",lwd=3)
  lines(ages,mSel20,col="black",lwd=3)

  legend("topleft",inset=0.02,legend=c(eqSim_yrs,"20 yr mean","10 yr mean","5 yr mean (08-12)","5 yr mean (13-17)","3 yr mean"),
         col=c(cols,"black","grey20","grey40","grey60","grey80"),lty=1,
         lwd=c(rep(1,length(eqSim_yrs)),3,3,3,3,3),cex=0.7)
  savePlot(paste0("output/OM_plots/OM_",omName,"_",stkname,"_biological parameters - selection curves"),type="png")
  
  #are any years stat sig diff to mean?
  sigP<-NULL
  x<-mSelp5 #2008-2012
  y<-sel[,ac(2008:2012)]
    for (i in 1:ncol(y)){
    print(colnames(y)[i])
    p<-(t.test(x=x,y=y[,i],paired=T)$p.value)
    print(p)
    if(p<0.05){
      sigP<-rbind(sigP,colnames(y)[i])
    }
  }
  
  x<-mSel5 #2013-2017
  y<-sel[,ac(2013:2017)]
  
  for (i in 1:ncol(y)){
    print(colnames(y)[i])
    p<-(t.test(x=x,y=y[,i],paired=T)$p.value)
    print(p)
    if(p<0.05){
      sigP<-rbind(sigP,colnames(y)[i])
    }
  }
  
  # when does LYC contribute to stock
  tmp<-iterMedians(stock.n(stk_stf)*stock.wt(stk_stf))[,ac(2000:2017)]
  a<-matrix(rep(t(c(colSums(tmp,na.rm=T))),9),ncol=9)
  biomass<-FLQuant(0,dim=c(9,length(2000:2017)),dimnames=list(age=ac(0:8),year=ac(2000:2017)))
  biomass[]<-t(a)
  tmp1<-tmp/biomass
  
  layout(1:length(sigP))
  for (i in 1:length(sigP)){
    plot(ages,c(tmp1[,ac(sigP[i])]),type="o",pch=16,main=sigP[i],ylab="proportion of stock biomass")
  }
  savePlot(paste0("output/OM_plots/OM_",omName,"_",stkname,"_biological parameters - selection curves LYC check"),type="png")
  
  
  
  # when does LYC contribute to stock
  tmp<-iterMedians(stock.n(stk_stf)*stock.wt(stk_stf))[,ac(2000:2017)]
  a<-matrix(rep(t(c(colSums(tmp,na.rm=T))),9),ncol=9)
  biomass<-FLQuant(0,dim=c(9,length(2000:2017)),dimnames=list(age=ac(0:8),year=ac(2000:2017)))
  biomass[]<-t(a)
  tmp1<-tmp/biomass
  
  layout(1:length(sigP))
  for (i in 1:length(sigP)){
  plot(ages,c(tmp1[,ac(sigP[i])]),pch=16,main=sigP[i])
  }
  
  
  ### recruitment
  
  # test AR
  
  acf(detrend(matrix(window(rec(stk_orig), start = min(rec.period)))),main="Recruitment")
  savePlot(paste0("output/OM_plots/OM_",omName,"_",stkname,"_recruitment -  acf"),type="png",res=600)


  ### Cohort growth 

  # Checks:
  # Are large year classes different to other cohorts?
  
  #plot cohort growth
  #stk_coh is stock.wt
  
  stk_coh<-FLCohort(window(stock.wt(stk_fwd),end=max(hist_yrs)))
  
  dat<-FLQuant(stk_coh,dimnames=list(age=ages,year=1964:2017,iter=1:1000))
  
  yr<-sample(x=hist_yrs[1:(length(hist_yrs)-5)],size=10,replace=F)
  ggplot(data=dat[,ac(yr)], aes(age, data)) + geom_line(aes(group=year, colour=factor(year))) + 
    facet_wrap(~year, scales="free_y", nrow=3) + theme(legend.position = "none") + ggtitle("cohort")
  
  savePlot(filename=paste0("output/OM_plots/OM_",omName,"_",stkname,"_cohort growth"),type="png",res=600)
  
  #Are LYC years sig diff from mean?
  LYC<-c(1974,1979,1999,2005,2009,2014)
  years<-as.numeric(dimnames(stk_coh)$cohort)
  tmp<-array(stk_coh,dim=c(9,54,1000),dimnames=list(ages,dimnames(stk_coh)$cohort,1:1000))
  mcoh<-apply(tmp[,ac(years[-which(years %in% LYC)]),],c(1),median,na.rm=T)
  
  
  #Are recent spike years sig diff from mean?
  windows()
  LYC<-c(2005,2009,2014)
  #plot all cohorts
  plot(ages,mcoh,ylim=range(tmp,na.rm=T),ylab="stock wt",col="red", main="stock.wt by age by cohort",type="n")
  for (i in years[-which(years %in% LYC)]){
    lines(ages,tmp[,ac(i),1],col="red",lty=2)
  }
  #for (i in 1:3){
  lines(ages,tmp[,ac(LYC[1]),1],col="blue")
  lines(ages,tmp[,ac(LYC[2]),1],col="blue")
  lines(ages,tmp[,ac(LYC[3]),1],col="blue")
  
  legend("topleft",inset=0.02,legend=c("recent large years","other years (2000+)"),col=c("blue","red"),lty=c(1,2),cex=0.7)
  #  }
  savePlot(filename=paste0("output/OM_plots/OM_",omName,"_",stkname,"_cohort growth recent LYC check"),type="png",res=600)
  
  
  
  # growth rate vs time and rec size
  tmp<-array(stk_coh,dim=c(9,54,1000),dimnames=list(ages,dimnames(stk_coh)$cohort,1:1000))
  tmp<-tmp[,ac(unique(c(eqSim_yrs-8,eqSim_yrs))),]
  tmp_med<-apply(tmp,c(1,2),median,na.rm=T) # median of 1000 iters
  #cut down to full years
  years<-colnames(tmp_med)
  #tmp_med<-tmp_med[,ac(years[years<2013])]
  g<-matrix(NA,nrow=ncol(tmp_med))
  for (i in 1:ncol(tmp_med)){
    data<-data.frame(age=ages,wt=tmp_med[,i])
    g[i]<-coef(lm(wt~age,data=data))["age"]
  }
  
  ##
  plot(years,g,xlab="years",ylab="slope",pch=16,cex=1.2,main="Growth rate by cohort")
  points("2005",g[which(years == 2005)],pch=16,col="blue",cex=1.2)
  points("2009",g[which(years == 2009)],pch=16,col="blue",cex=1.2)
  points("2014",g[which(years == 2014)],pch=16,col="blue",cex=1.2)
  legend("topleft",inset=0.02,legend=c("recent large years","other years (2000+)"),col=c("blue","black"),pch=16,cex=0.7)
  savePlot(filename=paste0("output/OM_plots/OM_",omName,"_",stkname,"_cohort growth rate recent LYC check"),type="png",res=600)
  

  ### OEM check
 
  # Check:
  # comparison of 10 yr mean with actual values
  windows()
  param.list<-c("stock.wt","catch.wt","landings.wt","discards.wt","mat","m")
  col.list<-rainbow(length(eqSim_yrs))
  for (par.nm in param.list){
    
    dat<-eval(parse(text=paste0("stk_oem@",par.nm)))
    plot(ages,c(dat[,ac(2018),,,,1]),col="black",lwd=3,type="n",ylim=range(dat),ylab=gsub("[.]"," ",par.nm))
    # plot each eqSim year
    for (i in 1:length(eqSim_yrs)){
      lines(ages,c(dat[,ac(eqSim_yrs[i]),,,,1]),col=col.list[i],lty=2)
    }
    # plot OEM value
    lines(ages,c(dat[,ac(2018),,,,1]),col="black",lwd=2)
    legend("topleft",inset=0.02,legend=c(eqSim_yrs,"OEM"),col=c(col.list,"black"),lty=c(rep(2,length(eqSim_yrs)),1))
    savePlot(paste0("output/OM_plots/OM_",omName,"_",stkname,"_biological parameters - oem check - ",gsub("[.]"," ",par.nm),".png"), type="png",res=600)
  }
  

}
