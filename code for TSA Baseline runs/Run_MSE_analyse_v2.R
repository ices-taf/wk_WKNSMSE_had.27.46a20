### ------------------------------------------------------------------------ ###
### process results ####
### ------------------------------------------------------------------------ ###
library(FLCore)
library(ggplot2)
library(tidyr)
library(cowplot)
library(dplyr)

library(doParallel)
cl <- makeCluster(1)
registerDoParallel(cl)

setwd(paste("/home/coleh/WKNSMSE/wk_WKNSMSE_had.27.46a20", sep=""))
### load additional functions
source("a4a_mse_WKNSMSE_funs.R")


. <- foreach(i = seq(length(cl))) %dopar% {
  #devtools::load_all("../mse/")
library(FLCore)
library(ggplot2)
library(tidyr)
library(cowplot)
library(dplyr)
library(ggplotFL)

setwd(paste("~/WKNSMSE/wk_WKNSMSE_had.27.46a20", sep=""))
### load additional functions
source("a4a_mse_WKNSMSE_funs.R")

}


### reference points
refpts <- list(msyBtrigger = 132000,
               Blim=94000,
               Flim=0.38,
               Fmsy = 0.194,
               Fpa = 0.274,
               Bpa = 132000)

omName<-"Baseline"
#omName<-"Alt1"
n<-999 # number of iters
#n<-1000
### ------------------------------------------------------------------------ ###
### get files ####
### ------------------------------------------------------------------------ ###

path_res <- paste0("output/had4/runs/",omName,"/",n,"_20/")
path_out <- paste0("output/had4/runs/",omName,"/",n,"_20/plots/")
path_mp<-"input/had4/MP_base/"
mp_file<-paste0("MPBase_OM",omName,"_",n,".rds")

files_res <- data.frame(file = list.files(path_res, pattern = "*.rds"), 
                        stringsAsFactors = FALSE)
files_res <- files_res[!grepl(x = files_res$file, pattern = "*stats*"),, 
                       drop = FALSE]
files_res <- files_res[!grepl(x = files_res$file, pattern = "F0.rds"),, 
                       drop = FALSE]


files_res$OM <- omName 
files_res$Ftrgt <- as.numeric(
  gsub(x = regmatches(x = files_res$file, 
                      m = regexpr(text = files_res$file, 
                                  pattern = "Ftrgt-0.[0-9]{1,}")),
       pattern = "Ftrgt-", replacement = ""))
files_res$Btrigger <- as.numeric(
  gsub(x = regmatches(x = files_res$file, 
                      m = regexpr(text = files_res$file,
                                  pattern = "Btrigger-[0-9]{1,}[e+]{0,}[0-9]{0,}")),
       pattern = "Btrigger-", replacement = ""))
files_res$HCR <- gsub(x = regmatches(x = files_res$file, 
                                     m = regexpr(text = files_res$file, 
                                                 pattern = "HCR-[A-Z]{1,}")),
                      pattern = "HCR-", replacement = "")
files_res$TACconstr <- as.logical(
  gsub(x = regmatches(x = files_res$file, 
                      m = regexpr(text = files_res$file, 
                                  pattern = "TACconstr-TRUE|TACconstr-FALSE")),
       pattern = "TACconstr-", replacement = ""))
files_res$BB <-  as.logical(
  gsub(x = regmatches(x = files_res$file,
                      m = regexpr(text = files_res$file, 
                                  pattern = "BB-TRUE|BB-FALSE")),
       pattern = "BB-", replacement = ""))


stats <- readRDS(paste0("output/had4/res_plots/","stats.rds"))
stats <- read.csv(paste0("output/had4/res_plots/", "stats.csv"))
stats_new <- merge(stats, files_res, all = TRUE)

# set Blim depending on OM
stats_new$Blim <-94000

### keep only new files
if(omName %in% "Baseline"){
stats_new <- stats_new[!stats_new$file %in% stats$file, ]
}else{
stats_new <- stats_new[is.na(stats_new$F_maxed), ]
}
res_list <- foreach(i = seq(nrow(stats_new))) %dopar% {
  readRDS(paste0(path_res, stats_new$file[i]))
}

### ------------------------------------------------------------------------ ###
### calculate non convergence and remove iterations for stats calcs ####
### ------------------------------------------------------------------------ ###

### SAM convergence
stats_new$conv_failed <- foreach(x = res_list, .packages = "FLCore",
                                 .combine = "c") %dopar% {
                                   sum(x@tracking["conv.est", ac(2018:2037)] != 0)
                                 }
all(stats_new$conv_failed == 0)

pos <- which(stats_new$conv_failed != 0)

#remove iterations that failed to converge
if(length(pos)>0){
  for (i in 1:length(pos)){
    
    #fl<- #readRDS(paste0(path_res, as.character(tmp$file[i])))
    pos_iter<-NULL
    for (k in 1:n){
      if(sum(res_list[[pos[i]]]@tracking["conv.est",,,,,k],na.rm=T)>0){
        pos_iter<-c(pos_iter,k)
      }
    }
    stk2<-iter(res_list[[pos[i]]]@stock,c(1:999)[-pos_iter])
    res_list[[pos[i]]]@stock<-stk2
  }
}


### ------------------------------------------------------------------------ ###
### calculate summary statistics ####
### ------------------------------------------------------------------------ ###
### calculate for short- (year 1-5), medium- (year 6-10) and 
### long-term (year 11-20)
### risk 1: proportion of stock below Blim, average over iterations and period
### risk 3: maximum of annual proportions
### catch: mean catch in period (mean over years and iterations)
### iav: inter-annual variation of catch, average over years and iterations


### catch
stats_new$catch_long <- foreach(x = res_list, .packages = "FLCore",
                                .combine = "c") %dopar% {
                                  mean(window(catch(x@stock), start = 2029))
                                }
stats_new$catch_short <- foreach(x = res_list, .packages = "FLCore",
                                 .combine = "c") %dopar% {
                                   mean(window(catch(x@stock), start = 2019, end = 2023))
                                 }
stats_new$catch_medium <- foreach(x = res_list, .packages = "FLCore",
                                  .combine = "c") %dopar% {
                                    mean(window(catch(x@stock), start = 2024, end = 2028))
                                  }
### catch median
stats_new$catch_median_long <- foreach(x = res_list, .packages = "FLCore",
                                       .combine = "c") %dopar% {
                                         median(window(catch(x@stock), start = 2029))
                                       }
stats_new$catch_median_short <- foreach(x = res_list, .packages = "FLCore",
                                        .combine = "c") %dopar% {
                                          median(window(catch(x@stock), start = 2019, end = 2023))
                                        }
stats_new$catch_median_medium <- foreach(x = res_list, .packages = "FLCore",
                                         .combine = "c") %dopar% {
                                           median(window(catch(x@stock), start = 2024, end = 2028))
                                         }
### risks
stats_new$risk1_full <- foreach(x = res_list, 
                                Blim = stats_new$Blim[seq(nrow(stats_new))],
                                .packages = "FLCore",
                                .combine = "c") %dopar% {
                                  mean(window(ssb(x@stock), start = 2019) < Blim)
                                }
stats_new$risk1_long <- foreach(x = res_list, 
                                Blim = stats_new$Blim[seq(nrow(stats_new))],
                                .packages = "FLCore",
                                .combine = "c") %dopar% {
                                  mean(window(ssb(x@stock), start = 2029) < Blim)
                                }
stats_new$risk1_short <- foreach(x = res_list, 
                                 Blim = stats_new$Blim[seq(nrow(stats_new))],
                                 .packages = "FLCore",
                                 .combine = "c") %dopar% {
                                   mean(window(ssb(x@stock), start = 2019, end = 2023) < Blim)
                                 }
stats_new$risk1_medium <- foreach(x = res_list, 
                                  Blim = stats_new$Blim[seq(nrow(stats_new))],
                                  .packages = "FLCore",
                                  .combine = "c") %dopar% {
                                    mean(window(ssb(x@stock), start = 2024, end = 2028) < Blim)
                                  }
stats_new$risk3_long <- foreach(x = res_list, 
                                Blim = stats_new$Blim[seq(nrow(stats_new))],
                                .packages = "FLCore",
                                .combine = "c") %dopar% {
                                  max(iterMeans(window(ssb(x@stock), start = 2029) < Blim))
                                }
stats_new$risk3_short <- foreach(x = res_list, 
                                 Blim = stats_new$Blim[seq(nrow(stats_new))],
                                 .packages = "FLCore",
                                 .combine = "c") %dopar% {
                                   max(iterMeans(window(ssb(x@stock), start = 2019, end = 2023) < Blim))
                                 }
stats_new$risk3_medium <- foreach(x = res_list, 
                                  Blim = stats_new$Blim[seq(nrow(stats_new))],
                                  .packages = "FLCore",
                                  .combine = "c") %dopar% {
                                    max(iterMeans(window(ssb(x@stock), start = 2024, end = 2028) < Blim))
                                  }
### inter-annual variation of catch
stats_new$iav_long <- foreach(x = res_list, .packages = "FLCore",
                              .combine = "c") %dopar% {
                                iav(object = catch(window(stock(x), start = 2028)), summary_all = median)
                              }
stats_new$iav_short <- foreach(x = res_list, .packages = "FLCore",
                               .combine = "c") %dopar% {
                                 iav(object = catch(window(stock(x), start = 2018, end = 2023)), 
                                     summary_all = median)
                               }
stats_new$iav_medium <- foreach(x = res_list, .packages = "FLCore",
                                .combine = "c") %dopar% {
                                  iav(object = catch(window(stock(x), start = 2023, end = 2028)), 
                                      summary_all = mean)
                                }
### SSB
stats_new$ssb_median_long <- foreach(x = res_list, .packages = "FLCore",
                                     .combine = "c") %dopar% {
                                       median(window(ssb(x@stock), start = 2029))
                                     }
stats_new$ssb_median_short <- foreach(x = res_list, .packages = "FLCore",
                                      .combine = "c") %dopar% {
                                        median(window(ssb(x@stock), start = 2019, end = 2023))
                                      }
stats_new$ssb_median_medium <- foreach(x = res_list, .packages = "FLCore",
                                       .combine = "c") %dopar% {
                                         median(window(ssb(x@stock), start = 2024, end = 2028))
                                       }
### time to recovery
MSYBtrigger <- refpts$msyBtrigger
stats_new$recovery_proportion <- foreach(x = res_list, .packages = "FLCore",
                                         .combine = "c") %dopar% {
                                           mean(apply(window(ssb(x@stock), start = 2019) >= MSYBtrigger, 6, max))
                                         }
stats_new$recovery_time <- foreach(x = res_list, .packages = "FLCore",
                                   .combine = "c") %dopar% {
                                     median(apply(window(ssb(x@stock), start = 2019)@.Data >= MSYBtrigger, 6, 
                                                  function(x) {
                                                    if (any(x)) {
                                                      which(x)[1]
                                                    } else {
                                                      Inf
                                                    }
                                                  }))
                                   }
### fbar
stats_new$F_median_long <- foreach(x = res_list, .packages = "FLCore",
                                   .combine = "c") %dopar% {
                                     median(window(fbar(x@stock), start = 2029))
                                   }

stats_new$F_median_medium <- foreach(x = res_list, .packages = "FLCore",
                                     .combine = "c") %dopar% {
                                       median(window(fbar(x@stock),  start = 2024, end = 2028))
                                     }
stats_new$F_median_short<- foreach(x = res_list, .packages = "FLCore",
                                   .combine = "c") %dopar% {
                                     median(window(fbar(x@stock),start = 2019, end = 2023) )
                                   }


### check F maxed (2)
stats_new$F_maxed <- foreach(x = res_list, .packages = "FLCore",
                             .combine = "c") %dopar% {
                               sum(window(fbar(stock(x)), start = 2019) >= 2)
                             }

stats <- rbind(stats, stats_new)
stats <- stats[order(stats$file), ]
saveRDS(object = stats, file = paste0("output/had4/res_plots/", "stats.rds"))
write.csv(x = stats, file = paste0("output/had4/res_plots","/stats.csv"), row.names = FALSE)

### ------------------------------------------------------------------------ ###
### plot functions ####
### ------------------------------------------------------------------------ ###


### plot function
grid <- function(dat, HCR = "A",
                 time = c("long", "short", "medium"),
                 add_risk1 = FALSE, highlight_max = FALSE) {
  
  ### catch
  dat$catch <- dat[, paste0("catch_median_", time)]
  dat$risk <- dat[, paste0("risk3_", time)]
  dat$iav <- dat[, paste0("iav_", time)]
  dat$ssb <- dat[, paste0("ssb_median_", time)]
  dat$risk1 <- dat[, paste0("risk1_", time)]
  dat$Btrigger <- dat$Btrigger / 1000
  ### find yield maximum
  dat_max <- dat %>% 
    filter(risk <= 0.05) %>%
    filter(catch == max(catch)) %>%
    select(Ftrgt, Btrigger)
  # p1 <- ggplot(data = dat, 
  #               aes(x = Btrigger, y = Ftrgt, fill = catch)) +
  #   geom_raster() +
  #   scale_fill_gradient(paste0(time, "-term\ncatch (median)"), low = "red", 
  #                       high = "green") +
  #   geom_text(aes(label = round(catch), colour = risk <= 0.05),
  #             size = 2) +
  #   scale_colour_manual("risk <= 0.05", 
  #                       values = c("FALSE" = "red", "TRUE" = "black")) +
  #   theme_bw()
  p1 <- ggplot() +
    geom_raster(data = dat %>% 
                  filter(risk <= 0.05) %>%
                  filter(catch >= 0.95 * max(catch)),
                aes(x = Btrigger, y = Ftrgt, fill = catch)) +
    scale_fill_gradient(paste0("yield maximum\narea [t]"), low = "red",
                        high = "green") +
    geom_text(data = dat, 
              aes(x = Btrigger, y = Ftrgt, 
                  label = round(catch), colour = risk <= 0.05),
              size = 2) +
    scale_colour_manual("risk <= 0.05", 
                        values = c("FALSE" = "red", "TRUE" = "black")) +
    theme_bw() +
    facet_wrap(~ paste0("median ", time, "-term catch [t]")) +
    scale_x_continuous(breaks = c(seq(from = 110, to = 190, by = 10))) +
    labs(x = expression(B[trigger]~"[1000t]"),
         y = expression(F[trgt]))
  ### risk
  p2 <- ggplot(data = dat, 
               aes(x = Btrigger, y = Ftrgt)) +
    geom_raster(aes(fill = risk), alpha = 0.75) +
    scale_fill_gradient(paste0(time, "-term\nrisk 3"), 
                        low = "green", high = "red") +
    geom_text(aes(label = round(risk, 3), colour = risk <= 0.05),
              size = 2) +
    scale_colour_manual("risk <= 0.05", 
                        values = c("FALSE" = "red", "TRUE" = "black")) +
    theme_bw() +
    facet_wrap(~ paste0(time, "-term risk 3")) +
    scale_x_continuous(breaks = c(seq(from = 110, to = 190, by = 10))) +
    labs(x = expression(B[trigger]~"[1000t]"),
         y = expression(F[trgt]))
  ### iav
  p3 <- ggplot(data = dat, 
               aes(x = Btrigger, y = Ftrgt)) +
    geom_raster(aes(fill = iav)) +
    geom_text(aes(label = round(iav, 3), colour = risk <= 0.05),
              size = 2) +
    scale_colour_manual("risk <= 0.05",
                        values = c("FALSE" = "red", "TRUE" = "black")) +
    scale_fill_gradient(paste0("median\n", time,
                               "-term\ninter-annual\ncatch variability"),
                        low = "green", high = "red") +
    theme_bw() +
    facet_wrap(~ paste0("median ", time, 
                        "-term inter-annual\ catch variability")) +
    scale_x_continuous(breaks = c(seq(from = 110, to = 190, by = 10))) +
    labs(x = expression(B[trigger]~"[1000t]"),
         y = expression(F[trgt]))
  ### SSB
  p4 <- ggplot(data = dat, 
               aes(x = Btrigger, y = Ftrgt)) +
    geom_raster(aes(fill = ssb), alpha = 0.5) +
    geom_text(aes(label = round(ssb), colour = risk <= 0.05),
              size = 2) +
    scale_colour_manual("risk <= 0.05",
                        values = c("FALSE" = "red", "TRUE" = "black")) +
    scale_fill_gradient(paste0("median\n", time,
                               "-term\nSSB [t]"),
                        low = "red", high = "green") +
    theme_bw() +
    facet_wrap(~ paste0("median ", time, 
                        "-term SSB [t]")) +
    scale_x_continuous(breaks = c(seq(from = 110, to = 190, by = 10))) +
    labs(x = expression(B[trigger]~"[1000t]"),
         y = expression(F[trgt]))
  ### risk1
  p5 <- ggplot(data = dat, 
               aes(x = Btrigger, y = Ftrgt)) +
    geom_raster(aes(fill = risk1), alpha = 0.75) +
    scale_fill_gradient(paste0(time, "-term\nrisk 1"), 
                        low = "green", high = "red") +
    geom_text(aes(label = round(risk1, 3), colour = risk1 <= 0.05),
              size = 2) +
    scale_colour_manual("risk <= 0.05", 
                        values = c("FALSE" = "red", "TRUE" = "black")) +
    theme_bw() +
    facet_wrap(~ paste0(time, "-term risk 1")) +
    scale_x_continuous(breaks = c(seq(from = 110, to = 190, by = 10))) +
    labs(x = expression(B[trigger]~"[1000t]"),
         y = expression(F[trgt]))
  ### highlight maximum
  if (isTRUE(highlight_max)) {
    p_add <- geom_tile(data = dat_max, aes(x = Btrigger, y = Ftrgt),
                       width = 10, height = 0.01, linetype = "solid",
                       alpha = 0, colour = "black", size = 0.3)
    p1 <- p1 + p_add
    p2 <- p2 + p_add
    p3 <- p3 + p_add
    p4 <- p4 + p_add
    p5 <- p5 + p_add
  }
  
  if (isTRUE(add_risk1)) {
    plot_grid(p1, p2, p3, p5, nrow = 2, ncol = 2, align = "hv")
  } else {
    plot_grid(p1, p2, p3, p4, nrow = 2, ncol = 2, align = "hv")
  }
  
  
}



### ------------------------------------------------------------------------ ###
### plot grid ####
### ------------------------------------------------------------------------ ###

### A
grid(dat = stats %>%
       filter(HCR == "A" & BB == FALSE & TACconstr == FALSE &
                Ftrgt %in% round(seq(0, 1, 0.01), 2) & OM == "Baseline"), 
     HCR = "A", time = "long", add_risk1 = FALSE, highlight_max = TRUE)
ggsave(filename = paste0(path_out,"grid_A_long.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
grid(dat = stats %>%
       filter(HCR == "A" & BB == FALSE & TACconstr == FALSE &
                Ftrgt %in% round(seq(0, 1, 0.01), 2) & OM == "Baseline"),HCR = "A", time = "medium")
ggsave(filename = paste0(path_out,"grid_A_medium.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
grid(dat = stats %>%
       filter(HCR == "A" & BB == FALSE & TACconstr == FALSE &
                Ftrgt %in% round(seq(0, 1, 0.01), 2) & OM == "Baseline"),HCR = "A", time = "short")
ggsave(filename = paste0(path_out,"/grid_A_short.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")

### B
grid(dat = stats %>%
       filter(HCR == "B" & BB == FALSE & TACconstr == FALSE &
                Ftrgt %in% round(seq(0, 1, 0.01), 2) & OM == "Baseline"), 
     HCR = "B", time = "long", add_risk1 = FALSE, highlight_max = TRUE)
ggsave(filename = paste0(path_out,"/grid_B_long.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
grid(dat = stats %>%
       filter(HCR == "B" & BB == FALSE & TACconstr == FALSE &
                Ftrgt %in% round(seq(0, 1, 0.01), 2) & OM == "Baseline"), 
     HCR = "B", time = "medium")
ggsave(filename = paste0(path_out,"/grid_B_medium.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
grid(dat = stats %>%
       filter(HCR == "B" & BB == FALSE & TACconstr == FALSE &
                Ftrgt %in% round(seq(0, 1, 0.01), 2) & OM == "Baseline"), 
     HCR = "B", time = "short")
ggsave(filename = paste0(path_out,"/grid_B_short.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")

### C
grid(dat = stats %>%
       filter(HCR == "C" & BB == FALSE & TACconstr == FALSE &
                Ftrgt %in% round(seq(0, 1, 0.01), 2) & OM == "Baseline"), 
     HCR = "C", time = "long", add_risk1 = FALSE, highlight_max = TRUE)
ggsave(filename = paste0(path_out,"/grid_C_long.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
grid(dat = stats %>%
       filter(HCR == "C" & BB == FALSE & TACconstr == FALSE &
                Ftrgt %in% round(seq(0, 1, 0.01), 2) & OM == "Baseline"), 
     HCR = "C", time = "medium")
ggsave(filename = paste0(path_out,"/grid_C_medium.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
grid(dat = stats %>%
       filter(HCR == "C" & BB == FALSE & TACconstr == FALSE &
                Ftrgt %in% round(seq(0, 1, 0.01), 2) & OM == "Baseline"), 
     HCR = "C", time = "short")
ggsave(filename = paste0(path_out,"/grid_C_short.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")

### AD
grid(dat = stats %>%
       filter(HCR == "A" & BB == TRUE & TACconstr == TRUE &
                Ftrgt %in% round(seq(0, 1, 0.01), 2) & OM == "Baseline"), 
     HCR = "A", time = "long", add_risk1 = FALSE, highlight_max = TRUE)
ggsave(filename = paste0(path_out,"/grid_AD_long.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
grid(dat = stats %>%
       filter(HCR == "A" & BB == TRUE & TACconstr == TRUE &
                Ftrgt %in% round(seq(0, 1, 0.01), 2) & OM == "Baseline"), 
     HCR = "A", time = "medium")
ggsave(filename = paste0(path_out,"/grid_AD_medium.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
grid(dat = stats %>%
       filter(HCR == "A" & BB == TRUE & TACconstr == TRUE &
                Ftrgt %in% round(seq(0, 1, 0.01), 2) & OM == "Baseline"), 
     HCR = "A", time = "short")
ggsave(filename = paste0(path_out,"/grid_AD_short.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")


### BE
grid(dat = stats %>%
       filter(HCR == "B" & BB == TRUE & TACconstr == TRUE &
                Ftrgt %in% round(seq(0, 1, 0.01), 2) & OM == "Baseline"), 
     HCR = "B", time = "long", add_risk1 = FALSE, highlight_max = TRUE)
ggsave(filename = paste0(path_out,"/grid_BE_long.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
grid(dat = stats %>%
       filter(HCR == "B" & BB == TRUE & TACconstr == TRUE &
                Ftrgt %in% round(seq(0, 1, 0.01), 2) & OM == "Baseline"), 
     HCR = "B", time = "medium")
ggsave(filename = paste0(path_out,"/grid_BE_medium.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
grid(dat = stats %>%
       filter(HCR == "B" & BB == TRUE & TACconstr == TRUE &
                Ftrgt %in% round(seq(0, 1, 0.01), 2) & OM == "Baseline"), 
     HCR = "B", time = "short")
ggsave(filename = paste0(path_out,"/grid_BE_short.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")


### CE
grid(dat = stats %>%
       filter(HCR == "C" & BB == TRUE & TACconstr == TRUE &
                Ftrgt %in% round(seq(0, 1, 0.01), 2) & OM == "Baseline"), 
     HCR = "C", time = "long", add_risk1 = FALSE, highlight_max = TRUE)
ggsave(filename = paste0(path_out,"/grid_CE_long.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
grid(dat = stats %>%
       filter(HCR == "C" & BB == TRUE & TACconstr == TRUE &
                Ftrgt %in% round(seq(0, 1, 0.01), 2) & OM == "Baseline"), 
     HCR = "C", time = "medium")
ggsave(filename = paste0(path_out,"/grid_CE_medium.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
grid(dat = stats %>%
       filter(HCR == "C" & BB == TRUE & TACconstr == TRUE &
                Ftrgt %in% round(seq(0, 1, 0.01), 2) & OM == "Baseline"), 
     HCR = "C", time = "short")
ggsave(filename = paste0(path_out,"/grid_CE_short.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")



### ------------------------------------------------------------------------ ###
### plot runs that reached Fmax ####
### ------------------------------------------------------------------------ ###

### #plot up Fmax

pos <- which(stats$F_maxed != 0)
tmp<-stats[pos, ]

if(length(pos)>0){
  for (i in 1:length(pos)){
    file.nm<-tmp$file[i]
    fl<-readRDS(paste0(path_res, tmp$file[i]))
    
    p1<-plot(stock(fl),main=file.nm)
    
    pos_iter <- which(apply(fbar(stock(fl)), c(1, 6), max) >= 2)
    for(j in 1:length(pos_iter)){
      if(j==1){
        p2<-plot(stock(fl)[,,,,, pos_iter[j]],main=paste0(file.nm," iter:",pos_iter[j]))
        
      }else{
        p3<-plot(stock(fl)[,,,,, pos_iter[j]],main=paste0(file.nm," iter:",pos_iter[j]))
        
      }
      
    }
    if(length(pos_iter)==1){
      plot_grid(p1, p2, nrow = 2, ncol = 2, align = "hv")
    }
    if(length(pos_iter)==2){
      plot_grid(p1, p2, p3, nrow = 2, ncol = 2, align = "hv")
    }
    ggsave(filename = paste0(path_out,gsub(".rds","",file.nm),"_Fmax.png"), 
           width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
  }
}

### ------------------------------------------------------------------------ ###
### extrapolate 5% risk line ####
### ------------------------------------------------------------------------ ###

# df_risks <- data.frame(Btrigger = c(seq(from = 110000, to = 190000, 
#                                         by = 10000)),
#                        Ftrgt = c(0.355, 0.355, 0.365, 0.375, 0.375, 0.385,
#                                  0.395, 0.405, 0.425))
# 
# lm_risks <- lm(formula = Ftrgt ~ Btrigger, data = tail(df_risks, 5))
# plot(Ftrgt ~ Btrigger, data = df_risks, 
#      xlim = c(110000, 300000), ylim = c(0.3, 0.6))
# abline(lm_risks)

