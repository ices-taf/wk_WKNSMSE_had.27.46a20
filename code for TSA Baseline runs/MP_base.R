
rm(list=ls())

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

setwd(paste("/home/coleh/WKNSMSE/wk_WKNSMSE_had.27.46a20", sep=""))


### ------------------------------------------------------------------------ ###
### load OM ####
### ------------------------------------------------------------------------ ###

#use Baseline OM with 10 iters

# OM ver
omName<-"Baseline"
#omName<-"Alt1"
#set n
n<-999
#n<-1000

load(file=paste0("input/had4/image_OM",omName,"_",n,".RData"))

source("a4a_mse_WKNSMSE_funs.R")
invisible(lapply(list.files(path = "functions/", pattern = "*.R$", 
                            full.names = TRUE), source))
# MP run name
MPrunName<-"BaseIntYrTACcont"
#
#test_iter<-66
#stk_fwd<-iter(stk_fwd,test_iter)
#stk_oem<-iter(stk_oem,test_iter)
#sr<-iter(sr,test_iter)
#sr_res<-iter(sr_res,test_iter)
#catch_res<-iter(catch_res,test_iter)
#idx[[1]]<-iter(idx[[1]],test_iter)
#idx[[2]]<-iter(idx[[2]],test_iter)
#idx_dev[[1]]<-iter(idx_dev[[1]],test_iter)
#idx_dev[[2]]<-iter(idx_dev[[2]],test_iter)
#
#n<-1
#iters<-1


### ------------------------------------------------------------------------ ###
### Set up for MP ####
### ------------------------------------------------------------------------ ###

### reference points
refpts_mse <- list(Btrigger = 132000,
                   Blim=94000,
                 #  Flim=0.38,
                   Ftrgt = 0.194,
                   Fpa = 0.274,
                   Bpa = 132000)

### some specifications for short term forecast with SAM
had_stf_def <- list(fwd_yrs_average = -3:-1, #3 year average
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
                                            proc_res = NULL
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
                       fwd_trgt = "TAC", fwd_yrs = 1, 
                       had_stf_def,
                       ### speeding SAM up
                       newtonsteps = 0, rel.tol = 0.001,
                       par_ini = list(sam_initial),
                       track_ini = TRUE, ### store ini for next year
                       ### SAM model specifications
                       conf = list(had_conf_sam),
                       parallel = FALSE ### TESTING ONLY
                     )),
  ctrl.phcr = mseCtrl(method = phcr_WKNSMSE,
                      args = refpts_mse),
  ctrl.hcr = mseCtrl(method = hcr_WKNSME, args = list(option = "A")),
  ctrl.is = mseCtrl(method = is_WKNSMSE_JarworskiGrowth_intyrTACcont, 
                    args = c(hcrpars = list(refpts_mse),
                             ### for short term forecast
                             fwd_trgt = c("TAC", "hcr"), fwd_yrs = 2,
                             had_stf_def#,
                             ### TAC constraint
                             #TAC_constraint = TRUE,
                             #lower = -Inf, upper = Inf,  # needs to be in 80 and 125 % catch_target/catch prev*100
                             #Btrigger_cond = FALSE,
                             ### banking and borrowing option D and E
                             #BB = FALSE,
                             #BB_conditional = FALSE,
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
        file = paste0("input/had4/MP_base/MP",MPrunName,"_OM",omName,"_",n,".rds"))

### ------------------------------------------------------------------------ ###
### run MSE ####
### ------------------------------------------------------------------------ ###

#
#### run MSE
#### WARNING: takes a while...
#### check normal execution
# res1 <- mp(om = input$om,
#           oem = input$oem,
#           #iem = iem,
#           ctrl.mp = input$ctrl.mp,
#           genArgs = input$genArgs,
#           tracking = input$tracking)

# saveRDS(object = res1, 
#        file = paste0("output/had4/MP_base/MP",MPrunName,"_OM",omName,"_",n,"_results.rds"))
