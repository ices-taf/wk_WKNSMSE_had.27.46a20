### ------------------------------------------------------------------------ ###
### R script to run WKNSMSE had MSE on HPC ####
### ------------------------------------------------------------------------ ###
### This is designed to be called by a job submission script
### run_mse.qsub for systems using PBS and the qsub commands
### run_mse.bsub for system using LSF and the bsub commands


### ------------------------------------------------------------------------ ###
### load arguments from job script ####
### ------------------------------------------------------------------------ ###

### load arguments
args <- commandArgs(TRUE)
print("arguments passed on to this script:")
print(args)

### evaluate arguments, if they are passed to R:
if (length(args) > 0) {
  
  ### extract arguments
  for (i in seq_along(args)) eval(parse(text = args[[i]]))
  
  ### parallelisation environment
  if (!exists("par_env")) par_env <- 1
  if (!exists("n_workers")) n_workers <- 1
  
} else {
  
  stop("no argument passed to R")
  
}

### ------------------------------------------------------------------------ ###
### set up environment ####
### ------------------------------------------------------------------------ ###

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

### load additional functions
setwd(paste("/home/coleh/WKNSMSE/wk_WKNSMSE_had.27.46a20", sep=""))
source("a4a_mse_WKNSMSE_funs.R")
invisible(lapply(list.files(path = "functions/", pattern = "*.R$", 
                            full.names = TRUE), source))
path_out <- paste0("output/had4/runs/", iters, "_", years)
### ------------------------------------------------------------------------ ###
### setup parallel environment ####
### ------------------------------------------------------------------------ ###
### par_env=1 -> MPI (Rmpi, DoMPI)
### par_env=2 -> DoParallel

### ------------------------------------------------------------------------ ###
### set HCR parameters 
  ### ------------------------------------------------------------------------ ###

    if (HCRoption %in% 1:6) {
    #options 1:12
hcr_vals <- expand.grid(
  Ftrgt = seq(from = 0.1, to = 0.3, by = 0.1),
  Btrigger = seq(from = 120000, to = 150000, by = 10000))

# option 13:24
hcr_vals<-rbind(hcr_vals,
                expand.grid(Ftrgt = seq(from = 0.28, to = 0.30, by = 0.01),
                            Btrigger = 180000),
                expand.grid(Ftrgt = seq(from = 0.27, to = 0.29, by = 0.01),
                            Btrigger = 170000),
                expand.grid(Ftrgt = seq(from = 0.26, to = 0.28, by = 0.01),
                            Btrigger = 160000),
                expand.grid(Ftrgt = seq(from = 0.25, to = 0.27, by = 0.01),
                            Btrigger = 150000)
)

hcr_vals<-rbind(hcr_vals,expand.grid(
  Ftrgt = seq(from = 0.2, to = 0.3, by = 0.1),
  Btrigger = seq(from = 160000, to = 180000, by = 10000)))
  
  hcr_vals<-rbind(hcr_vals,c(0.29,160000),c(0.28,150000),c(0.27,180000))
  
  # extra for AD
  hcr_vals<-rbind(hcr_vals,c(0.25,140000),c(0.26,140000),c(0.24,140000))
  
  # extra for BE
  hcr_vals<-rbind(hcr_vals,c(0.29,160000),c(0.28,150000),c(0.27,140000),c(0.26,140000),c(0.25,140000),c(0.24,140000))
  
  #extra for C
  hcr_vals<-rbind(hcr_vals,c(0.29,160000),c(0.28,150000))
}

  if (par_env == 1) {
    
    library(doMPI)
    cl <- startMPIcluster()
    registerDoMPI(cl)
    cl_length <- cl$workerCount
    
  } else if (par_env == 2) {
    
    library(doParallel)
    cl <- makeCluster(n_workers)
    registerDoParallel(cl)
    cl_length <- length(cl)
    
  }
  
  
#if(FALSE){  
  ### load packages and functions into workers
  . <- foreach(i = seq(cl_length)) %dopar% {
    #devtools::load_all("../mse/")
    library(mse)
    library(FLash)
    library(FLfse)
    library(stockassessment)
    library(foreach)
    library(doRNG)
    setwd(paste("/home/coleh/WKNSMSE/wk_WKNSMSE_had.27.46a20", sep=""))
    
    source("a4a_mse_WKNSMSE_funs.R")
    invisible(lapply(list.files(path = "functions/", pattern = "*.R$", 
                                full.names = TRUE), source))
  }
  
  ### set random seed for reproducibility
  library(doRNG)
  registerDoRNG(123)
 # }
  ### ------------------------------------------------------------------------ ###
  ### load data for MSE ####
  ### ------------------------------------------------------------------------ ###
  
  ### data path
  path_data <- paste0("input/had4/")
  
  omName<-"Baseline"
  MPrunName<-"Base"
  n<-iters
  
  ### load input objects
  input<-readRDS(file = paste0(path_data,"/MP_base/MP",MPrunName,"_OM",omName,"_",n,".rds"))
  
  ### modify input for running in parallel
  input$genArgs$nblocks <- nblocks
  


  ### ------------------------------------------------------------------------ ###
  ### set HCR option: A, B, C
  if (exists("HCRoption")) {
    
    input$ctrl.mp$ctrl.hcr@args$option <- switch(HCRoption, 
                                                 "1" = "A", 
                                                 "2" = "B", 
                                                 "3" = "C",
                                                 "4" = "A",
                                                 "5" = "B",
                                                 "6" = "C")
    
    cat(paste0("\nSetting custom HCR option: HCRoption = ", HCRoption, 
               " => HCR ", input$ctrl.mp$ctrl.hcr@args$option, "\n\n"))
    
  } else {
    
    cat(paste0("\nUsing default HCR option: HCR ", 
               input$ctrl.mp$ctrl.hcr@args$option, "\n\n"))
    HCRoption <- 0
  }
  
input_bckp <- input
for (HCR_comb in 37:42) {
  input <- input_bckp
  
  ### implement
  if (exists("HCR_comb")) {
    
    ### set Btrigger
    Btrigger <- hcr_vals[HCR_comb, "Btrigger"]
    input$ctrl.mp$ctrl.phcr@args$Btrigger <- Btrigger
    input$ctrl.mp$ctrl.is@args$hcrpars$Btrigger <- Btrigger
    
    ### set Ftrgt
    Ftrgt <- hcr_vals[HCR_comb, "Ftrgt"]
    input$ctrl.mp$ctrl.phcr@args$Ftrgt <- Ftrgt
    input$ctrl.mp$ctrl.is@args$hcrpars$Ftrgt <- Ftrgt
    
    cat(paste0("\nSetting custom Btrigger/Ftrgt values.\n",
               "Using HCR_comb = ", HCR_comb, "\n",
               "Ftrgt = ", Ftrgt, "\n",
               "Btrigger = ", Btrigger, "\n\n"))
    
  } else {
    input$ctrl.mp$ctrl.phcr@args$Ftrgt <- Ftrgt
    input$ctrl.mp$ctrl.is@args$hcrpars$Ftrgt <- Ftrgt
    input$ctrl.mp$ctrl.phcr@args$Btrigger <- Btrigger
    input$ctrl.mp$ctrl.is@args$hcrpars$Btrigger <- Btrigger
    
    cat(paste0("\nUsing default Btrigger/Ftrgt values.\n",
               "Ftrgt = ", input$ctrl.mp$ctrl.phcr@args$Ftrgt, "\n",
               "Btrigger = ", input$ctrl.mp$ctrl.phcr@args$Btrigger, "\n\n"))
    
  }
  
  ### ------------------------------------------------------------------------ ###
  ## TAC constraint
  input$ctrl.mp$ctrl.is@args$TAC_constraint <- FALSE
  ### check conditions
  ### either manually requested or as part of HCR options 4-6 
  if (exists("TAC_constraint")) {
    
    if (isTRUE(as.logical(TAC_constraint))) {
      input$ctrl.mp$ctrl.is@args$TAC_constraint <- TRUE
    }
  }
  if (HCRoption %in% 4:6) {
    input$ctrl.mp$ctrl.is@args$TAC_constraint <- TRUE
  }
  ### implement
  if (isTRUE(input$ctrl.mp$ctrl.is@args$TAC_constraint)) {
    
    input$ctrl.mp$ctrl.is@args$lower <- 80
    input$ctrl.mp$ctrl.is@args$upper <- 125
    input$ctrl.mp$ctrl.is@args$Btrigger_cond <- TRUE
    
    cat(paste0("\nImplementing TAC constraint.\n\n"))
    
  } else {
    
    cat(paste0("\nTAC constraint NOT implemented.\n\n"))
    
  }
  
  
  ### ------------------------------------------------------------------------ ###
  ### banking & borrowing
  ### banking & borrowing
  input$ctrl.mp$ctrl.is@args$BB <- FALSE
  input$iem <- NULL
  
  ### check conditions
  ### either manually requested or as part of HCR options 4-6 
  if (exists("BB")) {
    if (isTRUE(as.logical(BB))) {
      
      input$iem <- FLiem(method = iem_WKNSMSE, args = list(BB = TRUE))
      input$ctrl.mp$ctrl.is@args$BB <- TRUE
      input$ctrl.mp$ctrl.is@args$BB_check_hcr <- TRUE
      input$ctrl.mp$ctrl.is@args$BB_check_fc <- TRUE
      input$ctrl.mp$ctrl.is@args$BB_rho <- c(-0.1, 0.1)
      
    }
    
  }
  if (HCRoption %in% 4:6) {
    
    input$iem <- FLiem(method = iem_WKNSMSE, args = list(BB = TRUE))
    input$ctrl.mp$ctrl.is@args$BB <- TRUE
    input$ctrl.mp$ctrl.is@args$BB_rho <- c(-0.1, 0.1)
    input$ctrl.mp$ctrl.is@args$BB_check_hcr <- FALSE
    input$ctrl.mp$ctrl.is@args$BB_check_fc <- FALSE
    
    if (HCRoption %in% 4) {
      
      input$ctrl.mp$ctrl.is@args$BB_check_hcr <- TRUE
      
    } else if (HCRoption %in% 5:6) {
      
      input$ctrl.mp$ctrl.is@args$BB_check_fc <- TRUE
      
    }
    
  }
  
  if (!is.null(input$iem)) {
    
    cat(paste0("\nImplementing banking and borrowing.\n\n"))
    
  } else {
    
    cat(paste0("\nBanking and borrowing NOT implemented.\n\n"))
    
  }
  
  
  ### ------------------------------------------------------------------------ ###
  ### run MSE ####
  ### ------------------------------------------------------------------------ ###
  #debugonce(mse:::goFish)
  print(Sys.time())
  ### run MSE
  res1 <- mp(om = input$om,
             oem = input$oem,
             iem = input$iem,
             ctrl.mp = input$ctrl.mp,
             genArgs = input$genArgs,
             tracking = input$tracking)
  
  print(Sys.time())
  
  ### save results
   path_out <- paste0("output/had4/runs/",omName,"/", n, "_", years)
 
  dir.create(path = path_out, recursive = TRUE)
  file_out <- paste0("HCR-", input$ctrl.mp$ctrl.hcr@args$option,
                     "_Ftrgt-", input$ctrl.mp$ctrl.phcr@args$Ftrgt,
                     "_Btrigger-", input$ctrl.mp$ctrl.phcr@args$Btrigger,
                     "_TACconstr-", input$ctrl.mp$ctrl.is@args$TAC_constraint,
                     "_BB-", input$ctrl.mp$ctrl.is@args$BB
  )
  
  saveRDS(object = res1, paste0(path_out, "/", file_out, ".rds"))
}
  ### ------------------------------------------------------------------------ ###
  ### combine and plot ####
  ### ------------------------------------------------------------------------ ###
  if(FALSE){output <- readRDS(paste0(path_out, "/", file_out, ".rds"))
  # ### get stock before simulation
  stk <- input$om@stock
  # ### add simulated data
  stk[, dimnames(res1@stock)$year] <- res1@stock
  # ### save
  saveRDS(object = stk, file = paste0(path_out,"/",file_out,"_base_full_stk.rds"))
  
  # ### plot
  plot(stk, probs = c(0.05, 0.25, 0.5, 0.75, 0.95)) + 
    xlab("year") + geom_vline(xintercept = 2018.5) +
    geom_hline(data = data.frame(qname = "SSB", data = 94000), #blim
               aes(yintercept = data), linetype = "dashed") +
    geom_hline(data = data.frame(qname = "SSB", data = 132000), #bpa
               aes(yintercept = data), linetype = "solid") +
    geom_hline(data = data.frame(qname = "F", data = 0.38), # flim
               aes(yintercept = data), linetype = "dashed") +
    geom_hline(data = data.frame(qname = "F", data = 0.274), #fpa
               aes(yintercept = data), linetype = "solid") +
    theme_bw()
  ggsave(filename = paste0(path_out,"/",file_out,"_base_full_stk.png"), 
         width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
  
  }

                                      ### ------------------------------------------------------------------------ ###
                                      ### terminate ####
                                      ### ------------------------------------------------------------------------ ###
                                      
                                      ### close R
                                      # mpi.finalize()
                                      ### mpi.finalize() or mpi.quit() hang...
                                      ### -> kill R, the MPI processes stop afterwards
                                      
                                      ### try killing current job...
                                      # if (par_env == 1 & exists("kill")) {
                                      #   system("bkill $LSB_JOBID")
                                      # }
                                      
                                      quit(save = "no")
                                      
                                      
                                      
                                      