### ------------------------------------------------------------------------ ###
### process results ####
### ------------------------------------------------------------------------ ###
library(FLCore)
library(ggplot2)
library(tidyr)
library(cowplot)
library(dplyr)

library(doParallel)
cl <- makeCluster(30)
registerDoParallel(cl)

### load additional functions
source("a4a_mse_WKNSMSE_funs.R")

### ------------------------------------------------------------------------ ###
### get files ####
### ------------------------------------------------------------------------ ###

path_res <- "output/runs/had4/1000_20/"
path_res_stats<-"output/runs/had4/"
path_plot<-"output/runs/had4/1000_20/plots/"

# add Alt1
#path_res <- "output/runs/had4/999_20/Alt1/"
#path_res_stats<-"output/runs/had4/"

files_res <- data.frame(file = list.files(path_res, pattern = "*.rds"), 
                        stringsAsFactors = FALSE)
files_res <- files_res[!grepl(x = files_res$file, pattern = "*stats*"),, 
                       drop = FALSE]
files_res <- files_res[!grepl(x = files_res$file, pattern = "F0.rds"),, 
                       drop = FALSE]

# files_rename <- files_res$file
# files_rename <- files_rename[grep(x = files_rename, pattern = "^HCR*")]
# file.rename(from = paste0("output/runs/cod4/1000_20/", files_rename),
#             to = paste0("output/runs/cod4/1000_20/cod4_", files_rename))
files_res$OM <- sapply(strsplit(x = files_res$file, split = "\\_HCR"), "[[", 1)
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


#stats <- readRDS(paste0(path_res_stats, "stats.rds"))
stats<-read.csv(paste0(path_res_stats,"stats.csv"))
stats_new <- merge(stats, files_res, all = TRUE)
### set Blim depending on OM
stats_new$Blim <- sapply(stats_new$OM, function(x) {
  switch(x, "Baseline" = 94000,"Alt1" = 94000, "Alt2" = 94000,
         "intyrTACcont" = 94000)})
### keep only new files
stats_new <- stats_new[!stats_new$file %in% stats$file, ]

res_list <- foreach(i = seq(nrow(stats_new)), .packages = "mse") %dopar% {
  tmp <- readRDS(paste0(path_res, stats_new$file[i]))
  tmp@oem <- FLoem()
  tmp
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


# first remove failed iters
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
    for (k in 1:dim(res_list[[pos[i]]]@tracking)[6]){
      if(sum(res_list[[pos[i]]]@tracking["conv.est",,,,,k],na.rm=T)>0){
        pos_iter<-c(pos_iter,k)
      }
    }
    stk2<-iter(res_list[[pos[i]]]@stock,c(1:999)[-pos_iter])
    res_list[[pos[i]]]@stock<-stk2
  }
}


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
                                      summary_all = median)
                                }
### inter-annual variation of TAC
stats_new$iavTAC_long <- foreach(x = res_list, .packages = "FLCore",
                                 .combine = "c") %dopar% {
                                   iav(object = window(x@tracking["metric.is"], start = 2028, end = 2037), 
                                       summary_all = median)
                                 }
stats_new$iavTAC_short <- foreach(x = res_list, .packages = "FLCore",
                                  .combine = "c") %dopar% {
                                    iav(object = window(x@tracking["metric.is"], start = 2018, end = 2022), 
                                        summary_all = median)
                                  }
stats_new$iavTAC_medium <- foreach(x = res_list, .packages = "FLCore",
                                   .combine = "c") %dopar% {
                                     iav(object = window(x@tracking["metric.is"], start = 2023, end = 2027), 
                                         summary_all = median)
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
### fbar
stats_new$fbar_median_long <- foreach(x = res_list, .packages = "FLCore",
                                      .combine = "c") %dopar% {
                                        median(window(fbar(x@stock), start = 2029))
                                      }
stats_new$fbar_median_short <- foreach(x = res_list, .packages = "FLCore",
                                       .combine = "c") %dopar% {
                                         median(window(fbar(x@stock), start = 2019, end = 2023))
                                       }
stats_new$fbar_median_medium <- foreach(x = res_list, .packages = "FLCore",
                                        .combine = "c") %dopar% {
                                          median(window(fbar(x@stock), start = 2024, end = 2028))
                                        }
### time to recovery
MSYBtrigger <- 132000
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

### check F maxed (2)
stats_new$F_maxed <- foreach(x = res_list, .packages = "FLCore",
                             .combine = "c") %dopar% {
                               sum(window(fbar(stock(x)), start = 2019) >= 2)
                             }

### F maxed in 12 scenarios, 
### always only once in one iteration

### proportion where MP is on slope
stats_new$slope_long <- foreach(x = res_list, Ftrgt = stats_new$Ftrgt,
                                .packages = "FLCore", .combine = "c") %dopar% {
                                  mean(c(window(x@tracking["metric.hcr"], start = 2028) < 
                                           (Ftrgt * (1 - 1e-16))), 
                                       na.rm = TRUE)
                                }
stats_new$slope_medium <- foreach(x = res_list, Ftrgt = stats_new$Ftrgt,
                                  .packages = "FLCore", .combine = "c") %dopar% {
                                    mean(c(window(x@tracking["metric.hcr"], start = 2023, end = 2027) < 
                                             (Ftrgt * (1 - 1e-16))), 
                                         na.rm = TRUE)
                                  }
stats_new$slope_short <- foreach(x = res_list, Ftrgt = stats_new$Ftrgt,
                                 .packages = "FLCore", .combine = "c") %dopar% {
                                   mean(c(window(x@tracking["metric.hcr"], start = 2018, end = 2022) < 
                                            (Ftrgt * (1 - 1e-16))), 
                                        na.rm = TRUE)
                                 }

stats <- rbind(stats, stats_new)
stats <- stats[order(stats$file), ]
saveRDS(object = stats, file = paste0(path_res_stats, "stats.rds"))
write.csv(x = stats, file = paste0(path_res_stats,"stats.csv"), row.names = FALSE)


#plot Fmax
pos <- which(stats$F_maxed != 0)
if(length(pos)>0){
  
  stats[pos, ]
  
  for (i in (1:length(pos))[-2]){
    res_i <- readRDS(paste0(path_res, stats$file[pos[i]]))
    
    plot(stock(res_i))
    ggsave(filename = paste0(path_plot,"/Fmax/","Fmax_",gsub(".rds","",stats$file[pos[i]]),".png"), 
           width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
    pos_iter <- which(apply(fbar(stock(res_i)), c(1, 6), max) >= 2)
    for (j in 1:length(pos_iter)){
      plot(stock(res_i)[,,,,, pos_iter[j]])
      ggsave(filename = paste0(path_plot,"/Fmax/","Fmax_iter_",pos_iter[j],"_",gsub(".rds","",stats$file[pos[i]]),".png"), 
             width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
        }
  }
}
### ------------------------------------------------------------------------ ###
### plot ####
### ------------------------------------------------------------------------ ###


### plot function
grid <- function(dat, HCR = "A",
                 time = c("long", "short", "medium"),
                 add_risk1 = FALSE, highlight_max = FALSE,
                 add_slope = FALSE) {
  
  ### catch
  dat$catch <- dat[, paste0("catch_median_", time)]
  dat$risk <- dat[, paste0("risk3_", time)]
  dat$iav <- dat[, paste0("iav_", time)]
  dat$ssb <- dat[, paste0("ssb_median_", time)]
  dat$risk1 <- dat[, paste0("risk1_", time)]
  dat$Btrigger <- dat$Btrigger / 1000
  dat$slope <- dat[, paste0("slope_", time)]
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
    scale_x_continuous(breaks = sort(unique(dat$Btrigger))) +
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
    scale_x_continuous(breaks = sort(unique(dat$Btrigger))) +
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
    scale_x_continuous(breaks = sort(unique(dat$Btrigger))) +
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
    scale_x_continuous(breaks = sort(unique(dat$Btrigger))) +
    labs(x = expression(B[trigger]~"[1000t]"),
         y = expression(F[trgt]))
  ### risk1
  p5 <- ggplot(data = dat, 
               aes(x = Btrigger, y = Ftrgt)) +
    geom_raster(aes(fill = risk1), alpha = 0.75) +
    scale_fill_gradient(paste0(time, "-term\nrisk 1"), 
                        low = "green", high = "red") +
    geom_text(aes(label = round(risk1, 3), colour = risk <= 0.05),
              size = 2) +
    scale_colour_manual("risk <= 0.05", 
                        values = c("FALSE" = "red", "TRUE" = "black")) +
    theme_bw() +
    facet_wrap(~ paste0(time, "-term risk 1")) +
    scale_x_continuous(breaks = sort(unique(dat$Btrigger))) +
    labs(x = expression(B[trigger]~"[1000t]"),
         y = expression(F[trgt]))
  p6 <- ggplot(data = dat, 
               aes(x = Btrigger, y = Ftrgt)) +
    geom_raster(aes(fill = slope), alpha = 0.75) +
    scale_fill_gradient(paste0(time, "-term\nproportion on\nHCR slope"), 
                        low = "green", high = "red") +
    geom_text(aes(label = round(slope, 3), colour = risk <= 0.05),
              size = 2) +
    scale_colour_manual("risk <= 0.05", 
                        values = c("FALSE" = "red", "TRUE" = "black")) +
    theme_bw() +
    facet_wrap(~ paste0(time, "-term proportion on HCR slope")) +
    scale_x_continuous(breaks = sort(unique(dat$Btrigger))) +
    labs(x = expression(B[trigger]~"[1000t]"),
         y = expression(F[trgt]))
  ### list with plots
  ps <- list(p1, p2, p3, p4, p5, p6)[c(rep(TRUE, 4), add_risk1, add_slope)]
  ### highlight maximum
  if (isTRUE(highlight_max)) {
    p_add <- geom_tile(data = dat_max, aes(x = Btrigger, y = Ftrgt),
                       width = 10, height = 0.01, linetype = "solid",
                       alpha = 0, colour = "black", size = 0.3)
    ps <- lapply(ps, function(x) {x + p_add})
  }
  ps$align = "hv"
  do.call(plot_grid, ps)
  
}

### A
grid(dat = stats %>%
       filter(HCR == "A" & BB == FALSE & TACconstr == FALSE &
                Ftrgt %in% c(0.194,round(seq(0, 1, 0.01), 2)) & OM == "Baseline"), 
     HCR = "A", time = "long", add_risk1 = FALSE, highlight_max = TRUE)
ggsave(filename = "output/runs/had4/1000_20/plots/grid/grid_A_long.png", 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
grid(dat = stats %>%
       filter(HCR == "A" & BB == FALSE & TACconstr == FALSE &
                Ftrgt %in% c(0.194,round(seq(0, 1, 0.01), 2)) & OM == "Baseline"), 
     HCR = "A", time = "medium")
ggsave(filename = "output/runs/had4/1000_20/plots/grid/grid_A_medium.png", 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
grid(dat = stats %>%
       filter(HCR == "A" & BB == FALSE & TACconstr == FALSE &
                Ftrgt %in% c(0.194,round(seq(0, 1, 0.01), 2)) & OM == "Baseline"), 
     HCR = "A", time = "short")
ggsave(filename = "output/runs/had4/1000_20/plots/grid/grid_A_short.png", 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")

### B
grid(dat = stats %>%
       filter(HCR == "B" & BB == FALSE & TACconstr == FALSE &
                Ftrgt %in% c(0.194,round(seq(0, 1, 0.01), 2)) & OM == "Baseline"), 
     HCR = "B", time = "long", add_risk1 = FALSE, highlight_max = TRUE)
ggsave(filename = "output/runs/had4/1000_20/plots/grid/grid_B_long.png", 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
grid(dat = stats %>%
       filter(HCR == "B" & BB == FALSE & TACconstr == FALSE &
                Ftrgt %in% c(0.194,round(seq(0, 1, 0.01), 2)) & OM == "Baseline"), 
     HCR = "B", time = "medium")
ggsave(filename = "output/runs/had4/1000_20/plots/grid/grid_B_medium.png", 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
grid(dat = stats %>%
       filter(HCR == "B" & BB == FALSE & TACconstr == FALSE &
                Ftrgt %in% c(0.194,round(seq(0, 1, 0.01), 2)) & OM == "Baseline"), 
     HCR = "B", time = "short")
ggsave(filename = "output/runs/had4/1000_20/plots/grid/grid_B_short.png", 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")

### C
grid(dat = stats %>%
       filter(HCR == "C" & BB == FALSE & TACconstr == FALSE &
                Ftrgt %in% c(0.194,round(seq(0, 1, 0.01), 2)) & OM == "Baseline"), 
     HCR = "C", time = "long", add_risk1 = FALSE, highlight_max = TRUE)
ggsave(filename = "output/runs/had4/1000_20/Baseline/plots/grid/grid_C_long.png", 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
grid(dat = stats %>%
       filter(HCR == "C" & BB == FALSE & TACconstr == FALSE &
                Ftrgt %in% c(0.194,round(seq(0, 1, 0.01), 2)) & OM == "Baseline"), 
     HCR = "C", time = "medium")
ggsave(filename = "output/runs/had4/1000_20/plots/grid/grid_C_medium.png", 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
grid(dat = stats %>%
       filter(HCR == "C" & BB == FALSE & TACconstr == FALSE &
                Ftrgt %in% c(0.194,round(seq(0, 1, 0.01), 2)) & OM == "Baseline"), 
     HCR = "C", time = "short")
ggsave(filename = "output/runs/had4/1000_20/plots/grid/grid_C_short.png", 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")

### AD
grid(dat = stats %>%
       filter(HCR == "A" & BB == TRUE & TACconstr == TRUE &
                Ftrgt %in% c(0.194,round(seq(0, 1, 0.01), 2)) & OM == "Baseline"), 
     HCR = "A", time = "long", add_risk1 = FALSE, highlight_max = TRUE)
ggsave(filename = "output/runs/had4/1000_20/plots/grid/grid_AD_long.png", 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
grid(dat = stats %>%
       filter(HCR == "A" & BB == TRUE & TACconstr == TRUE &
                Ftrgt %in% c(0.194,round(seq(0, 1, 0.01), 2)) & OM == "Baseline"), 
     HCR = "A", time = "medium")
ggsave(filename = "output/runs/had4/1000_20/plots/grid/grid_AD_medium.png", 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
grid(dat = stats %>%
       filter(HCR == "A" & BB == TRUE & TACconstr == TRUE &
                Ftrgt %in% c(0.194,round(seq(0, 1, 0.01), 2)) & OM == "Baseline"), 
     HCR = "A", time = "short")
ggsave(filename = "output/runs/had4/1000_20/plots/grid/grid_AD_short.png", 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")


### BE
grid(dat = stats %>%
       filter(HCR == "B" & BB == TRUE & TACconstr == TRUE &
                Ftrgt %in% c(0.194,round(seq(0, 1, 0.01), 2)) & OM == "Baseline"), 
     HCR = "B", time = "long", add_risk1 = FALSE, highlight_max = TRUE)
ggsave(filename = "output/runs/had4/1000_20/plots/grid/grid_BE_long.png", 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
grid(dat = stats %>%
       filter(HCR == "B" & BB == TRUE & TACconstr == TRUE &
                Ftrgt %in% c(0.194,round(seq(0, 1, 0.01), 2)) & OM == "Baseline"), 
     HCR = "B", time = "medium")
ggsave(filename = "output/runs/had4/1000_20/plots/grid/grid_BE_medium.png", 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
grid(dat = stats %>%
       filter(HCR == "B" & BB == TRUE & TACconstr == TRUE &
                Ftrgt %in% c(0.194,round(seq(0, 1, 0.01), 2)) & OM == "Baseline"), 
     HCR = "B", time = "short")
ggsave(filename = "output/runs/had4/1000_20/plots/grid/grid_BE_short.png", 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")


### CE
grid(dat = stats %>%
       filter(HCR == "C" & BB == TRUE & TACconstr == TRUE &
                Ftrgt %in% c(0.194,round(seq(0, 1, 0.01), 2)) & OM == "Baseline"), 
     HCR = "C", time = "long", add_risk1 = FALSE, highlight_max = TRUE)
ggsave(filename = "output/runs/had4/1000_20/plots/grid/grid_CE_long.png", 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
grid(dat = stats %>%
       filter(HCR == "C" & BB == TRUE & TACconstr == TRUE &
                Ftrgt %in% c(0.194,round(seq(0, 1, 0.01), 2)) & OM == "Baseline"), 
     HCR = "C", time = "medium")
ggsave(filename = "output/runs/had4/1000_20/plots/grid/grid_CE_medium.png", 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
grid(dat = stats %>%
       filter(HCR == "C" & BB == TRUE & TACconstr == TRUE &
                Ftrgt %in% c(0.194,round(seq(0, 1, 0.01), 2)) & OM == "Baseline"), 
     HCR = "C", time = "short")
ggsave(filename = "output/runs/had4/1000_20/plots/grid/grid_CE_short.png", 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")


### ------------------------------------------------------------------------ ###
### plot some runs ####
### ------------------------------------------------------------------------ ###
### function for plotting
plot_stk <- function(stats, OM_ = "Baseline", HCR_ = "A", Ftrgt_ = 0.28,
                     Btrigger_ = 180000, TACconstr_ = FALSE, BB_ = FALSE,
                     probs = c(0.05, 0.25, 0.5, 0.75, 0.95),
                     path_out = "output/runs/had4/1000_20/",
                     path_in = "input/had4/1000_20/", file_in = "MPbase_OMBaseline_1000.rds",
                     path_res = "output/runs/had4/1000_20/plots/stock_plots/",
                     save_plot = TRUE,
                     get_history = TRUE, overwrite_catch_history = FALSE,
                     yr_start = 2018.5,
                     Blim = 107000, MSYBtrigger = 150000,
                     Flim = 0.54, Fmsy = 0.31,
                     plot_iter = FALSE, iters_plot = 0,
                     width = 30, height = 20) {
  
  ### get scenario
  stats_i <- stats %>% filter(OM == OM_ & HCR == HCR_ & Ftrgt == Ftrgt_ & 
                                Btrigger == Btrigger_, TACconstr == TACconstr_ &
                                BB == BB_)
  ### load MSE results
  res_i <- readRDS(paste0(path_out, stats_i$file))
  
  ### load stock history
  if (isTRUE(get_history)) {
    hist_i <- readRDS(paste0(path_in, file_in))
    stk_i <- hist_i$om@stock
    stk_i[, dimnames(res_i@stock)$year] <- res_i@stock
    ### overwrite catch history
    if (isTRUE(overwrite_catch_history)) {
      catch_hist <- readRDS(paste0(path_in, OM_,"_catch_n.rds"))
      catch.n(stk_i)[, dimnames(catch_hist)$year] <- catch_hist
      catch(stk_i) <- computeCatch(stk_i)
    }
  }
  
  nonConv<-c(sum(res_i@tracking["conv.est", ac(2018:2037)] != 0))
  nIter<-dim(res_i@stock)[6]
  if (nonConv>0){
    pos_iter<-NULL 
    for (k in 1:nIter){
      if(sum(res_i@tracking["conv.est",,,,,k],na.rm=T)>0){
        pos_iter<-c(pos_iter,k)
      }
    }
    stk2<-iter(stk_i,c(1:nIter)[-pos_iter])
    stk_i<-stk2
  }
  
  ### plot
  if (!isTRUE(plot_iter)) {
    ### plot percentiles
    p <- plot(stk_i, probs = probs, iter = iters_plot) +
      xlab("year") + geom_vline(xintercept = yr_start) +
      geom_hline(data = data.frame(qname = "SSB", data = Blim),
                 aes(yintercept = data), linetype = "dashed") +
      geom_hline(data = data.frame(qname = "SSB", data = MSYBtrigger),
                 aes(yintercept = data), linetype = "solid") +
      geom_hline(data = data.frame(qname = "F", data = Flim),
                 aes(yintercept = data), linetype = "dashed") +
      geom_hline(data = data.frame(qname = "F", data = Fmsy),
                 aes(yintercept = data), linetype = "solid") +
      theme_bw() + ylim(0, NA) + theme(legend.position = 0) + 
      facet_wrap(~ qname, ncol = 1, strip.position = "right", scales = "free_y",
                 labeller = as_labeller(c(
                   "Rec" = "Rec [1000]",
                   "SSB" = "SSB [t]",
                   "Catch" = "Catch [t]",
                   "F" = paste0("F (ages ", paste(range(stk_i)["minfbar"], 
                                                  range(stk_i)["maxfbar"], sep = "-"),
                                ")"))))
  } else {
    ### plot iterations
    stk_df <- as.data.frame(FLQuants(`Rec [1000]` = rec(stk_i),
                                     `SSB [t]` = ssb(stk_i),
                                     `Catch [t]` = catch(stk_i),
                                     `F` = fbar(stk_i)))
    p <- ggplot(data = stk_df, 
                aes(x = year, y = data, group = iter)) +
      geom_line(alpha = 0.025) +
      facet_wrap(~ qname, ncol = 1, strip.position = "right",
                 scale = "free_y") +
      theme_bw() +
      ylim(c(0, NA)) + labs(y = "") + 
      geom_vline(xintercept = yr_start) +
      geom_hline(data = data.frame(qname = "SSB [t]", data = Blim),
                 aes(yintercept = data), linetype = "dashed") +
      geom_hline(data = data.frame(qname = "SSB [t]", data = MSYBtrigger),
                 aes(yintercept = data), linetype = "solid") +
      geom_hline(data = data.frame(qname = "F", data = Flim),
                 aes(yintercept = data), linetype = "dashed") +
      geom_hline(data = data.frame(qname = "F", data = Fmsy),
                 aes(yintercept = data), linetype = "solid")
    
  }
  print(p)
  ### save plot
  if (isTRUE(save_plot)) {
    filename <- paste0(path_res, gsub(x = stats_i$file, pattern = ".rds",
                                      replacement = ""), 
                       ifelse(isTRUE(iters_plot == 0), "", "_iters"),
                       ".png")
    ggsave(filename = filename,
           width = width, height = height, units = "cm", dpi = 300, 
           type = "cairo")
    
  }
  
}

### base OM
### current HCR
plot_stk(stats = stats, OM_ = "Baseline", HCR_ = "A", Ftrgt_ = 0.194, 
         Btrigger_ = 132000, TACconstr_ = TRUE, BB_ = TRUE, 
         overwrite_catch_history = TRUE, iters_plot = 7:12)
### A, B, C, AD, BE, CE
combs <- data.frame(name = c( "A*", "A", "B", "C", "AD", "BE", "CE"),
                    HCR = c( "A", "A", "B", "C", "A", "B", "C"),
                    BB = c(rep(FALSE, 4), TRUE, TRUE, TRUE),
                    TACconstr = c(rep(FALSE, 4), TRUE, TRUE, TRUE),
                    Btrigger = c(132000, 180000, 190000, 180000, 180000,
                                 170000, 160000),
                    Ftrgt = c(0.194, 0.28, 0.29, 0.28, 0.28, 0.27, 0.26))
combs$OM<-"Baseline"

lapply(X = split(combs, seq(nrow(combs))), FUN = function(i) {
  plot_stk(stats = stats, OM_ = i$OM, HCR_ = i$HCR, Ftrgt_ = i$Ftrgt, 
           Btrigger_ = i$Btrigger, TACconstr_ = i$TACconstr, BB_ = i$BB, 
           overwrite_catch_history = TRUE,file_in = paste0("MPbase_OM",i$OM,"_1000.rds"))
})
### plot with first five iterations
lapply(X = split(combs, seq(nrow(combs))), FUN = function(i) {
  plot_stk(stats = stats, OM_ = i$OM, HCR_ = i$HCR, Ftrgt_ = i$Ftrgt, 
           Btrigger_ = i$Btrigger, TACconstr_ = i$TACconstr, BB_ = i$BB, 
           overwrite_catch_history = TRUE, iters_plot = 7:12,file_in = paste0("MPbase_OM",i$OM,"_1000.rds"))
})

### add plot for A*D
plot_stk(stats = stats, OM_ = "Baseline", HCR_ = "A", Ftrgt_ = 0.194, 
         Btrigger_ = 132000, TACconstr_ = TRUE, BB_ = TRUE, 
         overwrite_catch_history = TRUE)
plot_stk(stats = stats, OM_ = "Baseline", HCR_ = "A", Ftrgt_ = 0.194, 
         Btrigger_ = 132000, TACconstr_ = TRUE, BB_ = TRUE, 
         overwrite_catch_history = TRUE, iters_plot = 7:12)

# Alt 1 OM
### A, B, C, AD, BE, CE
combs <- data.frame(name = c( "A*", "A", "B", "C", "AD", "BE", "CE"),
                    HCR = c( "A", "A", "B", "C", "A", "B", "C"),
                    BB = c(rep(FALSE, 4), TRUE, TRUE, TRUE),
                    TACconstr = c(rep(FALSE, 4), TRUE, TRUE, TRUE),
                    Btrigger = c(132000, 180000, 190000, 180000, 180000,
                                 170000, 160000),
                    Ftrgt = c(0.194, 0.28, 0.29, 0.28, 0.28, 0.27, 0.26))
combs$OM<-"Alt1"

lapply(X = split(combs, seq(nrow(combs))), FUN = function(i) {
  plot_stk(stats = stats, OM_ = i$OM, HCR_ = i$HCR, Ftrgt_ = i$Ftrgt, 
           Btrigger_ = i$Btrigger, TACconstr_ = i$TACconstr, BB_ = i$BB, 
           overwrite_catch_history = F,file_in = paste0("MPbase_OM",i$OM,"_999.rds"),
           path_in = "input/had4/999_20/")
})
### plot with first five iterations
lapply(X = split(combs, seq(nrow(combs))), FUN = function(i) {
  plot_stk(stats = stats, OM_ = i$OM, HCR_ = i$HCR, Ftrgt_ = i$Ftrgt, 
           Btrigger_ = i$Btrigger, TACconstr_ = i$TACconstr, BB_ = i$BB, 
           overwrite_catch_history = F, iters_plot = 7:12,file_in = paste0("MPbase_OM",i$OM,"_999.rds"),
           path_in = "input/had4/999_20/")
})


# Alt2 OM
### A, B, C, AD, BE, CE
combs <- data.frame(name = c( "A*", "A", "B", "C", "AD", "BE", "CE"),
                    HCR = c( "A", "A", "B", "C", "A", "B", "C"),
                    BB = c(rep(FALSE, 4), TRUE, TRUE, TRUE),
                    TACconstr = c(rep(FALSE, 4), TRUE, TRUE, TRUE),
                    Btrigger = c(132000, 180000, 190000, 180000, 180000,
                                 170000, 160000),
                    Ftrgt = c(0.194, 0.28, 0.29, 0.28, 0.28, 0.27, 0.26))
combs$OM<-"Alt2"

lapply(X = split(combs, seq(nrow(combs))), FUN = function(i) {
  plot_stk(stats = stats, OM_ = i$OM, HCR_ = i$HCR, Ftrgt_ = i$Ftrgt, 
           Btrigger_ = i$Btrigger, TACconstr_ = i$TACconstr, BB_ = i$BB, 
           overwrite_catch_history = T)
})
### plot with first five iterations
lapply(X = split(combs, seq(nrow(combs))), FUN = function(i) {
  plot_stk(stats = stats, OM_ = i$OM, HCR_ = i$HCR, Ftrgt_ = i$Ftrgt, 
           Btrigger_ = i$Btrigger, TACconstr_ = i$TACconstr, BB_ = i$BB, 
           overwrite_catch_history = T, iters_plot = 7:12)
})

# ### ------------------------------------------------------------------------ ###
# ### F=0 plot ####
# ### ------------------------------------------------------------------------ ###
# 
# stkF0 <- readRDS(file = "input/cod4/10000_100/data_F0.RData")$om@stock
# stkF0_res <- readRDS("output/runs/cod4/F0_10000_100.rds")@stock
# 
# catch_n <- readRDS("input/cod4/10000_100/catch_n.rds")
# catch.n(stkF0)[dimnames(catch_n)$age, dimnames(catch_n)$year] <- 
#   catch_n
# catch(stkF0) <- computeCatch(stkF0)
# 
# stkF0[, ac(2018:2118)] <- stkF0_res[, ac(2018:2118)]
# plot(stkF0, probs = c(0.05, 0.25, 0.5, 0.75, 0.95)) +
#   xlab("year") + geom_vline(xintercept = 2018.5) +
#   geom_hline(data = data.frame(qname = "SSB", data = 107000),
#              aes(yintercept = data), linetype = "dashed") +
#   geom_hline(data = data.frame(qname = "SSB", data = 150000),
#              aes(yintercept = data), linetype = "solid") +
#   geom_hline(data = data.frame(qname = "F", data = 0.54),
#              aes(yintercept = data), linetype = "dashed") +
#   geom_hline(data = data.frame(qname = "F", data = 0.31),
#              aes(yintercept = data), linetype = "solid") +
#   theme_bw() #+
# # geom_blank(data = as.data.frame(FLQuants(`Rec` = rec(stkF0),
# #                                          `SSB` = ssb(stkF0),
# #                                          `Catch` = catch(stkF0),
# #                                          `F` = fbar(stkF0))), 
# #            aes(x = year, y = data, group = iter))
# ggsave(filename = paste0("output/runs/cod4/1000_20/plots/stock_plots/", 
#                          "stk_F0_10000iters.png"),
#        width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
# ### iterations
# stkF0_df <- as.data.frame(FLQuants(`Rec [1000]` = rec(stkF0),
#                                    `SSB [t]` = ssb(stkF0),
#                                    `Catch [t]` = catch(stkF0),
#                                    `F` = fbar(stkF0)))
# ggplot(data = stkF0_df[stkF0_df$iter %in% 1:1000, ], 
#        aes(x = year, y = data, group = iter)) +
#   geom_line(alpha = 0.025) +
#   facet_wrap(~ qname, ncol = 1, strip.position = "right",
#              scale = "free_y") +
#   theme_bw() +
#   ylim(c(0, NA)) + labs(y = "") + 
#   geom_vline(xintercept = 2018.5) +
#   geom_hline(data = data.frame(qname = "SSB [t]", data = 107000),
#              aes(yintercept = data), linetype = "dashed") +
#   geom_hline(data = data.frame(qname = "SSB [t]", data = 150000),
#              aes(yintercept = data), linetype = "solid") +
#   geom_hline(data = data.frame(qname = "F", data = 0.54),
#              aes(yintercept = data), linetype = "dashed") +
#   geom_hline(data = data.frame(qname = "F", data = 0.31),
#              aes(yintercept = data), linetype = "solid")
# ggsave(filename = paste0("output/runs/cod4/1000_20/plots/stock_plots/", 
#                          "stk_F0_10000iters_iters.png"),
#        width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")


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


### ------------------------------------------------------------------------ ###
### summary plots: compare HCR options ####
### ------------------------------------------------------------------------ ###

### select maximum yield combinations
combs <- data.frame(name = c( "A*", "A", "B", "C", "AD", "BE", "CE"),
                    OM = c("Baseline"),
                    HCR = c( "A", "A", "B", "C", "A", "B", "C"),
                    BB = c(rep(FALSE, 4), TRUE, TRUE, TRUE),
                    TACconstr = c(rep(FALSE, 4), TRUE, TRUE, TRUE),
                    Btrigger = c(132000, 180000, 190000, 180000, 180000,
                                 170000, 160000),
                    Ftrgt = c(0.194, 0.28, 0.29, 0.28, 0.28, 0.27, 0.26),
                    scenario=0)

combs <- merge(combs, stats, all.x = TRUE)
combs2 <- gather(data = combs, key = "key", value = "value",
                 catch_median_long, risk3_long, iav_long,
                 ssb_median_long, recovery_proportion, recovery_time)
combs2$name <- factor(combs2$name, levels = c("A*", "A", "B", "C", "AD", "BE", "CE"))

ggplot(data = combs2, 
       mapping = aes(x = name, y = value, group = name)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ key, scales = "free_y") +
  theme_bw()

### load entire distribution for stats
stats_full <- function(data) {
  combs_full <- foreach(i = split(data, seq(nrow(data))), 
                        .packages = "FLCore", .combine = rbind) %do% {
                          
                          stk_i <- readRDS(paste0("output/runs/had4/1000_20/", i$file))
                          MSYBtrigger <- 132000
                          Blim <- 94000 #ifelse(!i$OM == "cod4_alt2", 107000, 110000)
                          
                          nonConv<-c(sum(stk_i@tracking["conv.est", ac(2018:2037)] != 0))
                          nIter<-dim(stk_i@stock)[6]
                          if (nonConv>0){
                            pos_iter<-NULL 
                            for (k in 1:nIter){
                              if(sum(stk_i@tracking["conv.est",,,,,k],na.rm=T)>0){
                                pos_iter<-c(pos_iter,k)
                              }
                            }
                            stk2<-iter(stk_i@stock,c(1:nIter)[-pos_iter])
                            stk_i@stock<-stk2
                          }
                          
                          res <- rbind(
                            data.frame(name = i$name, scenario = i$scenario,
                                       key = "catch_long",
                                       value = c(window(catch(stk_i@stock), start = 2029))),
                            data.frame(name = i$name, scenario = i$scenario,
                                       key = "catch_medium",
                                       value = c(window(catch(stk_i@stock), start = 2019, end = 2023))),
                            data.frame(name = i$name, scenario = i$scenario,
                                       key = "catch_short",
                                       value = c(window(catch(stk_i@stock), start = 2024, end = 2028))),
                            data.frame(name = i$name, scenario = i$scenario,
                                       key = "risk1_long",
                                       value = mean(window(ssb(stk_i@stock), start = 2029) < Blim)),
                            data.frame(name = i$name, scenario = i$scenario,
                                       key = "risk1_medium",
                                       value = mean(window(ssb(stk_i@stock), 
                                                           start = 2024, end = 2028) < Blim)),
                            data.frame(name = i$name, scenario = i$scenario,
                                       key = "risk1_short",
                                       value = mean(window(ssb(stk_i@stock), 
                                                           start = 2019, end = 2023) < Blim)),
                            data.frame(name = i$name, scenario = i$scenario,
                                       key = "risk3_long",
                                       value = max(iterMeans(window(ssb(stk_i@stock), 
                                                                    start = 2029) < Blim))),
                            data.frame(name = i$name, scenario = i$scenario,
                                       key = "risk3_medium",
                                       value = max(iterMeans(window(ssb(stk_i@stock), 
                                                                    start = 2024, end = 2028) < Blim))),
                            data.frame(name = i$name, scenario = i$scenario,
                                       key = "risk3_short",
                                       value = max(iterMeans(window(ssb(stk_i@stock), 
                                                                    start = 2019, end = 2023) < Blim))),
                            data.frame(name = i$name, scenario = i$scenario,
                                       key = "iav_long",
                                       value = c(iav(object = catch(window(stock(stk_i), start = 2028))))),
                            data.frame(name = i$name, scenario = i$scenario,
                                       key = "iav_medium",
                                       value = c(iav(object = catch(window(stock(stk_i), 
                                                                           start = 2023, end = 2028))))),
                            data.frame(name = i$name, scenario = i$scenario,
                                       key = "iav_short",
                                       value = c(iav(object = catch(window(stock(stk_i), 
                                                                           start = 2018, end = 2023))))),
                            data.frame(name = i$name, scenario = i$scenario,
                                       key = "ssb_long",
                                       value = c(window(ssb(stk_i@stock), start = 2029))),
                            data.frame(name = i$name, scenario = i$scenario,
                                       key = "ssb_medium",
                                       value = c(window(ssb(stk_i@stock), start = 2024, end = 2028))),
                            data.frame(name = i$name, scenario = i$scenario,
                                       key = "ssb_short",
                                       value = c(window(ssb(stk_i@stock), start = 2019, end = 2023))),
                            data.frame(name = i$name, scenario = i$scenario,
                                       key = "recovery_proportion",
                                       value = mean(apply(window(ssb(stk_i@stock), 
                                                                 start = 2019) >= MSYBtrigger, 6, max))),
                            data.frame(name = i$name, scenario = i$scenario,
                                       key = "recovery_time",
                                       value = c(apply(window(ssb(stk_i@stock), 
                                                              start = 2019)@.Data >= MSYBtrigger, 6,
                                                       function(x) {
                                                         if (any(x)) {which(x)[1]} else {Inf}})))
                          )
                          if (i$HCR == "F0") {
                            res$value[res$key %in% c("catch_long", "catch_medium", "catch_short",
                                                     "iav_long", "iav_medium", "iav_short")] <- 0
                          }
                          res <- merge(res, i[, c("name", "OM", "HCR", "BB", "TACconstr", "Btrigger",
                                                  "Ftrgt")])
                          return(res)
                        }
  combs_full$name <- factor(combs_full$name, 
                            levels = c("A*", "A", "B", "C", "AD", "BE", "CE"))
  return(combs_full)
}

### base OM
combs_base <- stats_full(data = combs)
ggplot(data = combs_base, 
       mapping = aes(x = name, y = value, group = name)) +
  #geom_bar(stat = "identity") +
  geom_boxplot() + 
  facet_wrap(~ key, scales = "free_y") +
  theme_bw() +
  ylim(0, NA)

### get median for option A*
combs_base <- left_join(combs_base, 
                        combs_base %>%
                          group_by(key, OM, name) %>%
                          summarise(value_median = median(value)) %>%
                          filter(name == "A*") %>%
                          select(-name))
p_catch_long <- ggplot(data = combs_base[combs_base$key == "catch_long", ], 
                       mapping = aes(x = name, y = value)) +
  geom_boxplot() + theme_bw() + ylim(0, NA) +
  geom_hline(aes(yintercept = value_median), colour = "red", alpha = 0.5) +
  labs(x = "", y = "long-term catch [t]") +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE, digits = 1))
p_catch_medium <- ggplot(data = combs_base[combs_base$key == "catch_medium", ], 
                         mapping = aes(x = name, y = value)) +
  geom_boxplot() + theme_bw() + ylim(0, NA) +
  geom_hline(aes(yintercept = value_median), colour = "red", alpha = 0.5) +
  labs(x = "", y = "medium-term catch [t]") +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE, digits = 1))
p_catch_short <- ggplot(data = combs_base[combs_base$key == "catch_short", ], 
                        mapping = aes(x = name, y = value)) +
  geom_boxplot() + theme_bw() + ylim(0, NA) +
  geom_hline(aes(yintercept = value_median), colour = "red", alpha = 0.5) +
  labs(x = "", y = "short-term catch [t]") +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE, digits = 1))
p_risk1_long <- ggplot(data = combs_base[combs_base$key == "risk1_long", ], 
                       mapping = aes(x = name, y = value)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  geom_blank(data = combs_base[combs_base$key == "risk3_long", ]) +
  geom_hline(aes(yintercept = value_median), colour = "red", alpha = 0.5) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "long-term risk 1")
p_risk1_medium <- ggplot(data = combs_base[combs_base$key == "risk1_medium", ], 
                         mapping = aes(x = name, y = value)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") + 
  geom_blank(data = combs_base[combs_base$key == "risk3_medium", ]) +
  geom_hline(aes(yintercept = value_median), colour = "red", alpha = 0.5) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "medium-term risk 1")
p_risk1_short <- ggplot(data = combs_base[combs_base$key == "risk1_short", ], 
                        mapping = aes(x = name, y = value)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") + 
  geom_blank(data = combs_base[combs_base$key == "risk3_short", ]) +
  geom_hline(aes(yintercept = value_median), colour = "red", alpha = 0.5) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "short-term risk 1")
p_risk3_long <- ggplot(data = combs_base[combs_base$key == "risk3_long", ], 
                       mapping = aes(x = name, y = value)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") + 
  geom_blank(data = combs_base[combs_base$key == "risk1_long", ]) +
  geom_hline(aes(yintercept = value_median), colour = "red", alpha = 0.5) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "long-term risk 3")
p_risk3_medium <- ggplot(data = combs_base[combs_base$key == "risk3_medium", ], 
                         mapping = aes(x = name, y = value)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") + 
  geom_blank(data = combs_base[combs_base$key == "risk1_medium", ]) +
  geom_hline(aes(yintercept = value_median), colour = "red", alpha = 0.5) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "medium-term risk 3")
p_risk3_short <- ggplot(data = combs_base[combs_base$key == "risk3_short", ], 
                        mapping = aes(x = name, y = value)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") + 
  geom_blank(data = combs_base[combs_base$key == "risk1_short", ]) +
  geom_hline(aes(yintercept = value_median), colour = "red", alpha = 0.5) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "short-term risk 3")
p_iav_long <- ggplot(data = combs_base[combs_base$key == "iav_long", ], 
                     mapping = aes(x = name, y = value)) +
  geom_boxplot() + theme_bw() + ylim(0, NA) +
  geom_hline(aes(yintercept = value_median), colour = "red", alpha = 0.5) +
  coord_cartesian(ylim = c(0, 1)) + 
  labs(x = "", y = "long-term inter-annual catch variability")
p_iav_medium <- ggplot(data = combs_base[combs_base$key == "iav_medium", ], 
                       mapping = aes(x = name, y = value)) +
  geom_boxplot() + theme_bw() + ylim(0, NA) +
  geom_hline(aes(yintercept = value_median), colour = "red", alpha = 0.5) +
  coord_cartesian(ylim = c(0, 1)) + 
  labs(x = "", y = "medium-term inter-annual catch variability")
p_iav_short <- ggplot(data = combs_base[combs_base$key == "iav_short", ], 
                      mapping = aes(x = name, y = value)) +
  geom_boxplot() + theme_bw() + ylim(0, NA) +
  geom_hline(aes(yintercept = value_median), colour = "red", alpha = 0.5) +
  coord_cartesian(ylim = c(0, 1.2)) + 
  labs(x = "", y = "short-term inter-annual catch variability")
p_ssb_long <- ggplot(data = combs_base[combs_base$key == "ssb_long", ], 
                     mapping = aes(x = name, y = value)) +
  geom_boxplot() + theme_bw() + ylim(0, NA) +
  geom_hline(aes(yintercept = value_median), colour = "red", alpha = 0.5) +
  coord_cartesian(ylim = c(0, 1.2e+06)) +
  labs(x = "", y = "long-term SSB [t]") +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE, digits = 2))
p_ssb_medium <- ggplot(data = combs_base[combs_base$key == "ssb_medium", ], 
                       mapping = aes(x = name, y = value)) +
  geom_boxplot() + theme_bw() + ylim(0, NA) +
  geom_hline(aes(yintercept = value_median), colour = "red", alpha = 0.5) +
  coord_cartesian(ylim = c(0, 0.9e+06)) +
  labs(x = "", y = "medium-term SSB [t]") +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE, digits = 2),
                     limits = c(0, NA))
p_ssb_short <- ggplot(data = combs_base[combs_base$key == "ssb_short", ], 
                      mapping = aes(x = name, y = value)) +
  geom_boxplot() + theme_bw() + ylim(0, NA) +
  geom_hline(aes(yintercept = value_median), colour = "red", alpha = 0.5) +
  coord_cartesian(ylim = c(0, 0.65e+06)) +
  labs(x = "", y = "short-term SSB [t]") +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE, digits = 2),
                     limits = c(0, NA))
p_recovery_proportion <- 
  ggplot(data = combs_base[combs_base$key == "recovery_proportion", ], 
         mapping = aes(x = name, y = value)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") + 
  geom_hline(aes(yintercept = value_median), colour = "red", alpha = 0.5) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "recovery proportion")
p_recovery_time <- 
  ggplot(data = combs_base[combs_base$key == "recovery_time", ], 
         mapping = aes(x = name, y = value)) +
  geom_boxplot() + theme_bw() + ylim(0, NA) +
  geom_hline(aes(yintercept = value_median), colour = "red", alpha = 0.5) +
  labs(x = "", y = "recovery time [years]")

plot_grid(p_catch_long, p_risk1_long, p_risk3_long, p_iav_long, p_ssb_long,
          align = "hv")
ggsave(filename = paste0("output/runs/had4/1000_20/plots/baseOM_stats/", 
                         "summary_baseOM_long.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
plot_grid(p_catch_medium, p_risk1_medium, p_risk3_medium, p_iav_medium, 
          p_ssb_medium,
          align = "hv")
ggsave(filename = paste0("output/runs/had4/1000_20/plots/baseOM_stats/", 
                         "summary_baseOM_medium.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
plot_grid(p_catch_short, p_risk1_short, p_risk3_short, p_iav_short, p_ssb_short,
          align = "hv")
ggsave(filename = paste0("output/runs/had4/1000_20/plots/baseOM_stats/", 
                         "summary_baseOM_short.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
plot_grid(p_recovery_proportion, p_recovery_time,
          align = "hv")
ggsave(filename = paste0("output/runs/had4/1000_20/plots/baseOM_stats/", 
                         "summary_baseOM_recovery.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")

### ------------------------------------------------------------------------ ###
### summary plots: base OM, additional scenarios around maximum yield ####
### ------------------------------------------------------------------------ ###

### select maximum yield combinations
### add: 0.9 & 1.1 * Ftrgt
###      Fmsylower, Fmsyupper
combs <- data.frame(name = rep(c("A", "B", "C", "AD", "BE", "CE"), each = 5),
                    HCR = rep(c("A", "B", "C", "A", "B", "C"), each = 5),
                    BB = rep(c(rep(FALSE, 3), rep(TRUE, 3)), each = 5),
                    TACconstr = rep(c(rep(FALSE, 3), rep(TRUE, 3)), each = 5),
                    Btrigger = rep(c(180000, 190000, 180000, 180000, 170000,
                                     160000), each = 5),
                    Ftrgt = c(0.28 * c(0.9, 1, 1.1), 0.167, 0.194,
                              0.29 * c(0.9, 1, 1.1), 0.167, 0.194,
                              0.28 * c(0.9, 1, 1.1), 0.167, 0.194,
                              0.28 * c(0.9, 1, 1.1), 0.167, 0.194,
                              0.27 * c(0.9, 1, 1.1), 0.167, 0.194,
                              0.26 * c(0.9, 1, 1.1), 0.167, 0.194),
                    scenario = c("0.9*Ftrgt", "Ftrgt", "1.1*Ftrgt",
                                 "Fmsylower", "Fmsyupper"),
                    OM = "Baseline")
combs<-rbind(combs,data.frame(name = rep(c("A", "AD"), each = 2),
                              HCR = rep(c("A"), each = 4),
                              BB = rep(c(rep(FALSE, 1), rep(TRUE, 1)), each = 2),
                              TACconstr = rep(c(rep(FALSE, 1), rep(TRUE, 1)), each = 2),
                              Btrigger = c(180000*c(1.5,2),
                                           180000*c(1.5,2)),
                              Ftrgt = rep(c(0.28, 0.28), each = 2), 
                              scenario = c("1.5*Btrigger", "2*Btrigger"),
                              OM = "Baseline"))
combs <- merge(combs, stats)
combs_dat <- stats_full(data = combs)
combs_dat$scenario <- factor(combs_dat$scenario, 
                             levels = c("Fmsylower", "0.9*Ftrgt", "Ftrgt", 
                                        "1.1*Ftrgt", "Fmsyupper","1.5*Btrigger","2*Btrigger"))
ggplot(data = combs_dat, 
       mapping = aes(x = name, y = value, group = interaction(scenario, name), 
                     colour = scenario)) +
  geom_boxplot() + 
  facet_wrap(~ key, scales = "free_y") +
  theme_bw() +
  ylim(0, NA)

p_catch_long <- ggplot(data = combs_dat[combs_dat$key == "catch_long", ], 
                       mapping = aes(x = name, y = value, colour = scenario)) +
  geom_boxplot(show.legend = FALSE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "long-term catch [t]")
p_catch_medium <- ggplot(data = combs_dat[combs_dat$key == "catch_medium", ], 
                         mapping = aes(x = name, y = value, colour = scenario)) +
  geom_boxplot(show.legend = FALSE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) + 
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "medium-term catch [t]")
p_catch_short <- ggplot(data = combs_dat[combs_dat$key == "catch_short", ], 
                        mapping = aes(x = name, y = value, colour = scenario)) +
  geom_boxplot(show.legend = FALSE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) + 
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "short-term catch [t]")
p_risk1_long <- ggplot(data = combs_dat[combs_dat$key == "risk1_long", ], 
                       mapping = aes(x = name, y = value, fill = scenario,
                                     colour = scenario)) +
  geom_bar(show.legend = FALSE, stat = "identity", 
           position = position_dodge2(preserve = "single")) + 
  geom_blank(data = combs_dat[combs_dat$key == "risk3_long", ]) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "long-term risk 1")
p_risk1_medium <- ggplot(data = combs_dat[combs_dat$key == "risk1_medium", ], 
                         mapping = aes(x = name, y = value, fill = scenario,
                                       colour = scenario)) +
  geom_bar(show.legend = FALSE, stat = "identity", 
           position = position_dodge2(preserve = "single")) + 
  geom_blank(data = combs_dat[combs_dat$key == "risk3_medium", ]) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "medium-term risk 1")
p_risk1_short <- ggplot(data = combs_dat[combs_dat$key == "risk1_short", ], 
                        mapping = aes(x = name, y = value, fill = scenario,
                                      colour = scenario)) +
  geom_bar(show.legend = FALSE, stat = "identity", 
           position = position_dodge2(preserve = "single")) + 
  geom_blank(data = combs_dat[combs_dat$key == "risk3_short", ]) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "short-term risk 1")
p_risk3_long <- ggplot(data = combs_dat[combs_dat$key == "risk3_long", ], 
                       mapping = aes(x = name, y = value, fill = scenario,
                                     colour = scenario)) +
  geom_bar(show.legend = FALSE, stat = "identity", 
           position = position_dodge2(preserve = "single")) + 
  geom_blank(data = combs_dat[combs_dat$key == "risk1_long", ]) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "long-term risk 3")
p_risk3_medium <- ggplot(data = combs_dat[combs_dat$key == "risk3_medium", ], 
                         mapping = aes(x = name, y = value, fill = scenario,
                                       colour = scenario)) +
  geom_bar(show.legend = FALSE, stat = "identity", 
           position = position_dodge2(preserve = "single")) + 
  geom_blank(data = combs_dat[combs_dat$key == "risk1_medium", ]) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "medium-term risk 3")
p_risk3_short <- ggplot(data = combs_dat[combs_dat$key == "risk3_short", ], 
                        mapping = aes(x = name, y = value, fill = scenario,
                                      colour = scenario)) +
  geom_bar(show.legend = FALSE, stat = "identity", 
           position = position_dodge2(preserve = "single")) + 
  geom_blank(data = combs_dat[combs_dat$key == "risk1_short", ]) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "short-term risk 3")
p_iav_long <- ggplot(data = combs_dat[combs_dat$key == "iav_long", ], 
                     mapping = aes(x = name, y = value, colour = scenario)) +
  geom_boxplot(show.legend = FALSE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) +
  coord_cartesian(ylim = c(0, 1.3)) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "long-term inter-annual catch variability")
p_iav_medium <- ggplot(data = combs_dat[combs_dat$key == "iav_medium", ], 
                       mapping = aes(x = name, y = value, colour = scenario)) +
  geom_boxplot(show.legend = FALSE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) +
  coord_cartesian(ylim = c(0, 1.3)) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "medium-term inter-annual catch variability")
p_iav_short <- ggplot(data = combs_dat[combs_dat$key == "iav_short", ], 
                      mapping = aes(x = name, y = value, colour = scenario)) +
  geom_boxplot(show.legend = FALSE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) +
  coord_cartesian(ylim = c(0, 1.5)) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "short-term inter-annual catch variability")
p_ssb_long <- ggplot(data = combs_dat[combs_dat$key == "ssb_long", ], 
                     mapping = aes(x = name, y = value, colour = scenario)) +
  geom_boxplot(show.legend = TRUE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) +
  coord_cartesian(ylim = c(0, 5e+5)) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "long-term SSB [t]") +
  theme(legend.direction = "horizontal")
p_ssb_medium <- ggplot(data = combs_dat[combs_dat$key == "ssb_medium", ], 
                       mapping = aes(x = name, y = value, colour = scenario)) +
  geom_boxplot(show.legend = TRUE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) + 
  coord_cartesian(ylim = c(0, 5e+5)) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "medium-term SSB [t]") +
  theme(legend.direction = "horizontal")
p_ssb_short <- ggplot(data = combs_dat[combs_dat$key == "ssb_short", ], 
                      mapping = aes(x = name, y = value, colour = scenario)) +
  geom_boxplot(show.legend = TRUE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) +
  coord_cartesian(ylim = c(0, 4.5e+5)) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "short-term SSB [t]") +
  theme(legend.direction = "horizontal") +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE, digits = 1),
                     limits = c(0, NA))
p_recovery_proportion <- 
  ggplot(data = combs_dat[combs_dat$key == "recovery_proportion", ], 
         mapping = aes(x = name, y = value, fill = scenario, colour = scenario)) +
  geom_bar(show.legend = FALSE, stat = "identity",
           position = position_dodge2(preserve = "single")) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "recovery proportion")
p_recovery_time <- 
  ggplot(data = combs_dat[combs_dat$key == "recovery_time", ], 
         mapping = aes(x = name, y = value, colour = scenario)) +
  geom_boxplot(show.legend = TRUE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "recovery time [years]")

plot_grid(plot_grid(p_catch_long, p_risk1_long, p_risk3_long, p_iav_long,
                    p_ssb_long + theme(legend.position = "none"),
                    align = "hv"),
          get_legend(p_ssb_long), nrow = 2, rel_heights = c(1, 0.1))
ggsave(filename = paste0("output/runs/had4/1000_20/plots/baseOM_stats_combs/", 
                         "baseOM_combs_long.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
plot_grid(plot_grid(p_catch_medium, p_risk1_medium, p_risk3_medium, p_iav_medium,
                    p_ssb_medium + theme(legend.position = "none"),
                    align = "hv"),
          get_legend(p_ssb_medium), nrow = 2, rel_heights = c(1, 0.1))
ggsave(filename = paste0("output/runs/had4/1000_20/plots/baseOM_stats_combs/", 
                         "baseOM_combs_medium.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
plot_grid(plot_grid(p_catch_short, p_risk1_short, p_risk3_short, p_iav_short,
                    p_ssb_short + theme(legend.position = "none"),
                    align = "hv"),
          get_legend(p_ssb_short), nrow = 2, rel_heights = c(1, 0.1))
ggsave(filename = paste0("output/runs/had4/1000_20/plots/baseOM_stats_combs/", 
                         "baseOM_combs_short.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
plot_grid(plot_grid(p_recovery_proportion, 
                    p_recovery_time + theme(legend.position = "none"),
                    align = "hv"),
          get_legend(p_recovery_time), ncol = 2, rel_widths = c(0.5, 0.1))
ggsave(filename = paste0("output/runs/had4/1000_20/plots/baseOM_stats_combs/", 
                         "baseOM_combs_recovery.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")


### ------------------------------------------------------------------------ ###
### summary plots: compare alternative OMs ####
### ------------------------------------------------------------------------ ###

### alternative OMs
### select maximum yield combinations
combs_alt <- data.frame(name = c("A*", "A", "B", "C", "AD", "BE", "CE"),
                        HCR = c("A", "A", "B", "C", "A", "B", "C"),
                        BB = c(rep(FALSE, 4), rep(TRUE, 3)),
                        TACconstr = c(rep(FALSE, 4), rep(TRUE, 3)),
                        Btrigger = c(132000, 180000, 190000, 180000, 180000,
                                     170000, 160000),
                        Ftrgt = c(0.194, 0.28, 0.29, 0.28, 0.28, 0.27, 0.26),
                        scenario = 0)
combs_alt <- rbind(cbind(combs_alt, OM = "Baseline"),
                   cbind(combs_alt, OM = "Alt1"),
                   cbind(combs_alt, OM = "Alt2"))
combs_alt <- merge(combs_alt, stats)
combs_alt <- stats_full(data = combs_alt)
ggplot(data = combs_alt, 
       mapping = aes(x = name, y = value, group = interaction(OM, name), 
                     colour = OM)) +
  geom_boxplot() + 
  facet_wrap(~ key, scales = "free_y") +
  theme_bw() +
  ylim(0, NA)

p_catch_long <- ggplot(data = combs_alt[combs_alt$key == "catch_long", ], 
                       mapping = aes(x = name, y = value, colour = OM)) +
  geom_boxplot(show.legend = FALSE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "long-term catch [t]") +
  coord_cartesian(ylim = c(0, 2e+05)) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE, digits = 1),
                     limits = c(0, NA))
p_catch_medium <- ggplot(data = combs_alt[combs_alt$key == "catch_medium", ], 
                         mapping = aes(x = name, y = value, colour = OM)) +
  geom_boxplot(show.legend = FALSE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) + 
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "medium-term catch [t]") +
  coord_cartesian(ylim = c(0, 2e+05)) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE, digits = 2),
                     limits = c(0, NA))
p_catch_short <- ggplot(data = combs_alt[combs_alt$key == "catch_short", ], 
                        mapping = aes(x = name, y = value, colour = OM)) +
  geom_boxplot(show.legend = FALSE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) + 
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "short-term catch [t]") +
  coord_cartesian(ylim = c(0, 2.1e+05)) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE, digits = 2),
                     limits = c(0, NA))
p_risk1_long <- ggplot(data = combs_alt[combs_alt$key == "risk1_long", ], 
                       mapping = aes(x = name, y = value, fill = OM, 
                                     colour = OM)) +
  geom_bar(show.legend = FALSE, stat = "identity", 
           position = position_dodge2(preserve = "single")) + 
  geom_blank(data = combs_alt[combs_alt$key == "risk3_long", ]) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "long-term risk 1") +
  scale_y_continuous(breaks = c(0, 0.025, 0.05, 0.075, 0.1))
p_risk1_medium <- ggplot(data = combs_alt[combs_alt$key == "risk1_medium", ], 
                         mapping = aes(x = name, y = value, fill = OM,
                                       colour = OM)) +
  geom_bar(show.legend = FALSE, stat = "identity", 
           position = position_dodge2(preserve = "single")) + 
  geom_blank(data = combs_alt[combs_alt$key == "risk3_medium", ]) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "medium-term risk 1") +
  scale_y_continuous(breaks = c(0, 0.025, 0.05, 0.075, 0.1))
p_risk1_short <- ggplot(data = combs_alt[combs_alt$key == "risk1_short", ], 
                        mapping = aes(x = name, y = value, fill = OM,
                                      colour = OM)) +
  geom_bar(show.legend = FALSE, stat = "identity", 
           position = position_dodge2(preserve = "single")) + 
  geom_blank(data = combs_alt[combs_alt$key == "risk3_short", ]) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "short-term risk 1")
p_risk3_long <- ggplot(data = combs_alt[combs_alt$key == "risk3_long", ], 
                       mapping = aes(x = name, y = value, fill = OM,
                                     colour = OM)) +
  geom_bar(show.legend = FALSE, stat = "identity", 
           position = position_dodge2(preserve = "single")) + 
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "long-term risk 3") +
  scale_y_continuous(breaks = c(0, 0.025, 0.05, 0.075, 0.1))
p_risk3_medium <- ggplot(data = combs_alt[combs_alt$key == "risk3_medium", ], 
                         mapping = aes(x = name, y = value, fill = OM,
                                       colour = OM)) +
  geom_bar(show.legend = FALSE, stat = "identity", 
           position = position_dodge2(preserve = "single")) + 
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "medium-term risk 3") +
  scale_y_continuous(breaks = c(0, 0.025, 0.05, 0.075, 0.1))
p_risk3_short <- ggplot(data = combs_alt[combs_alt$key == "risk3_short", ], 
                        mapping = aes(x = name, y = value, fill = OM,
                                      colour = OM)) +
  geom_bar(show.legend = FALSE, stat = "identity", 
           position = position_dodge2(preserve = "single")) + 
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "short-term risk 3")
p_iav_long <- ggplot(data = combs_alt[combs_alt$key == "iav_long", ], 
                     mapping = aes(x = name, y = value, colour = OM)) +
  geom_boxplot(show.legend = FALSE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "long-term inter-annual catch variability")
p_iav_medium <- ggplot(data = combs_alt[combs_alt$key == "iav_medium", ], 
                       mapping = aes(x = name, y = value, colour = OM)) +
  geom_boxplot(show.legend = FALSE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "medium-term inter-annual catch variability")
p_iav_short <- ggplot(data = combs_alt[combs_alt$key == "iav_short", ], 
                      mapping = aes(x = name, y = value, colour = OM)) +
  geom_boxplot(show.legend = FALSE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) +
  coord_cartesian(ylim = c(0, 1.5)) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "short-term inter-annual catch variability")
p_ssb_long <- ggplot(data = combs_alt[combs_alt$key == "ssb_long", ], 
                     mapping = aes(x = name, y = value, colour = OM)) +
  geom_boxplot(show.legend = TRUE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "long-term SSB [t]") +
  theme(legend.direction = "horizontal") +
  coord_cartesian(ylim = c(0, 2.3e+06))
p_ssb_medium <- ggplot(data = combs_alt[combs_alt$key == "ssb_medium", ], 
                       mapping = aes(x = name, y = value, colour = OM)) +
  geom_boxplot(show.legend = TRUE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) + 
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "medium-term SSB [t]") +
  theme(legend.direction = "horizontal") +
  coord_cartesian(ylim = c(0, 1.8e+06))
p_ssb_short <- ggplot(data = combs_alt[combs_alt$key == "ssb_short", ], 
                      mapping = aes(x = name, y = value, colour = OM)) +
  geom_boxplot(show.legend = TRUE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "short-term SSB [t]") +
  theme(legend.direction = "horizontal") +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE, digits = 1),
                     limits = c(0, NA)) +
  coord_cartesian(ylim = c(0, 8e+05))
p_recovery_proportion <- 
  ggplot(data = combs_alt[combs_alt$key == "recovery_proportion", ], 
         mapping = aes(x = name, y = value, fill = OM, colour = OM)) +
  geom_bar(show.legend = FALSE, stat = "identity",
           position = position_dodge2(preserve = "single")) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "recovery proportion")
p_recovery_time <- 
  ggplot(data = combs_alt[combs_alt$key == "recovery_time", ], 
         mapping = aes(x = name, y = value, colour = OM)) +
  geom_boxplot(show.legend = TRUE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "recovery time [years]")

plot_grid(plot_grid(p_catch_long, p_risk1_long, p_risk3_long, p_iav_long,
                    p_ssb_long + theme(legend.position = "none"),
                    align = "hv"),
          get_legend(p_ssb_long), nrow = 2, rel_heights = c(1, 0.1))
ggsave(filename = paste0("output/runs/had4/1000_20/plots/altOMs_stats/", 
                         "summary_altOM_long.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
plot_grid(plot_grid(p_catch_medium, p_risk1_medium, p_risk3_medium, p_iav_medium,
                    p_ssb_medium + theme(legend.position = "none"),
                    align = "hv"),
          get_legend(p_ssb_medium), nrow = 2, rel_heights = c(1, 0.1))
ggsave(filename = paste0("output/runs/had4/1000_20/plots/altOMs_stats/", 
                         "summary_altOM_medium.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
plot_grid(plot_grid(p_catch_short, p_risk1_short, p_risk3_short, p_iav_short,
                    p_ssb_short + theme(legend.position = "none"),
                    align = "hv"),
          get_legend(p_ssb_short), nrow = 2, rel_heights = c(1, 0.1))
ggsave(filename = paste0("output/runs/had4/1000_20/plots/altOMs_stats/", 
                         "summary_altOM_short.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
plot_grid(plot_grid(p_recovery_proportion, 
                    p_recovery_time + theme(legend.position = "none"),
                    align = "hv"),
          get_legend(p_recovery_time), ncol = 2, rel_widths = c(0.5, 0.1))
ggsave(filename = paste0("output/runs/had4/1000_20/plots/altOMs_stats/", 
                         "summary_altOM_recovery.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")


### ------------------------------------------------------------------------ ###
### compare base OM with OM_alt3 (density dependent M) ####
### ------------------------------------------------------------------------ ###

OM_0_3_stats <- stats %>% 
  filter(OM %in% c("cod4", "cod4_alt3") &
           HCR == "A" & Btrigger == 170000 & Ftrgt == 0.38 &
           BB == FALSE)
OMs <- list(
  base = readRDS(paste0("output/runs/cod4/1000_20/",
                        "cod4_HCR-A_Ftrgt-0.38_Btrigger-170000_TACconstr-FALSE",
                        "_BB-FALSE.rds"))@stock,
  ddM = readRDS(paste0("output/runs/cod4/1000_20/",
                       "cod4_alt3_HCR-A_Ftrgt-0.38_Btrigger-170000_TACconstr",
                       "-FALSE_BB-FALSE.rds"))@stock)
input_OM <- list(
  base = readRDS("input/cod4/1000_20/base_run.rds")$om@stock,
  ddM = readRDS("input/cod4/1000_20/base_run.rds")$om@stock
)
OMs_df <- lapply(seq(input_OM), function(x) {
  M <- m(input_OM[[x]])
  M[, dimnames(m(OMs[[x]]))$year] <- m(OMs[[x]])
  SSB <- ssb(input_OM[[x]])
  SSB[, dimnames(m(OMs[[x]]))$year] <- ssb(OMs[[x]])
  qnts <- FLQuants("M at age 1" = M[1,], "M at age 2" = M[2,], 
                   "M at age 3" = M[3,], "M at age 4" = M[4,], 
                   "M at age 5" = M[5,], "M at age 6" = M[6,],
                   "SSB" = SSB)
  cbind(as.data.frame(qnts), OM = names(OMs)[x])
})
OMs_df <- do.call(rbind, OMs_df)
OMs_df <- OMs_df %>% group_by(age, year, qname, OM) %>%
  summarise(X0.05 = quantile(data, probs = 0.05),
            X0.25 = quantile(data, probs = 0.25),
            X0.50 = quantile(data, probs = 0.50),
            X0.75 = quantile(data, probs = 0.75),
            X0.95 = quantile(data, probs = 0.95))
p1 <- ggplot(data = OMs_df %>% filter(year >= 2000 & age == 1),
             aes(x = year, y = X0.50, fill = OM)) +
  geom_ribbon(aes(ymin = X0.05, ymax = X0.95), alpha = 0.3,show.legend = FALSE) +
  geom_ribbon(aes(ymin = X0.25, ymax = X0.75), alpha = 0.6, show.legend = FALSE) +
  geom_line(aes(colour = OM), show.legend = FALSE) +
  facet_grid(qname ~ OM, scales = "free_y") +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank())
p2 <- ggplot(data = OMs_df %>% filter(year >= 2000 & age == 2),
             aes(x = year, y = X0.50, fill = OM)) +
  geom_ribbon(aes(ymin = X0.05, ymax = X0.95), alpha = 0.3,show.legend = FALSE) +
  geom_ribbon(aes(ymin = X0.25, ymax = X0.75), alpha = 0.6, show.legend = FALSE) +
  geom_line(aes(colour = OM), show.legend = FALSE) +
  facet_grid(qname ~ OM, scales = "free_y") +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        strip.text.x = element_blank())
p3 <- ggplot(data = OMs_df %>% filter(year >= 2000 & age == 3),
             aes(x = year, y = X0.50, fill = OM)) +
  geom_ribbon(aes(ymin = X0.05, ymax = X0.95), alpha = 0.3,show.legend = FALSE) +
  geom_ribbon(aes(ymin = X0.25, ymax = X0.75), alpha = 0.6, show.legend = FALSE) +
  geom_line(aes(colour = OM), show.legend = FALSE) +
  facet_grid(qname ~ OM, scales = "free_y") +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        strip.text.x = element_blank())
p4 <- ggplot(data = OMs_df %>% filter(year >= 2000 & age == 4),
             aes(x = year, y = X0.50, fill = OM)) +
  geom_line(aes(colour = OM), show.legend = FALSE) +
  facet_grid(qname ~ OM, scales = "free_y") +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        strip.text.x = element_blank()) + ylim(c(0.199, 0.201))
p5 <- ggplot(data = OMs_df %>% filter(year >= 2000 & age == 5),
             aes(x = year, y = X0.50, fill = OM)) +
  geom_line(aes(colour = OM), show.legend = FALSE) +
  facet_grid(qname ~ OM, scales = "free_y") +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        strip.text.x = element_blank()) + ylim(c(0.199, 0.201))
p6 <- ggplot(data = OMs_df %>% filter(year >= 2000 & age == 6),
             aes(x = year, y = X0.50, fill = OM)) +
  geom_line(aes(colour = OM), show.legend = FALSE) +
  facet_grid(qname ~ OM, scales = "free_y") +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        strip.text.x = element_blank()) + ylim(c(0.199, 0.201))
p7 <- ggplot(data = OMs_df %>% filter(year >= 2000 & qname == "SSB"),
             aes(x = year, y = X0.50, fill = OM)) +
  geom_ribbon(aes(ymin = X0.05, ymax = X0.95), alpha = 0.3,show.legend = FALSE) +
  geom_ribbon(aes(ymin = X0.25, ymax = X0.75), alpha = 0.6, show.legend = FALSE) +
  geom_line(aes(colour = OM), show.legend = FALSE) +
  facet_grid(qname ~ OM, scales = "free_y") +
  theme_bw() +
  theme(axis.title.y = element_blank(), strip.text.x = element_blank())
plot_grid(p1, p2, p3, p4, p5, p6, p7, ncol = 1, align = "vh")
ggsave(filename = paste0("output/runs/cod4/1000_20/plots/altOMs_stats/", 
                         "baseOM_vs_ddM.png"), 
       width = 20, height = 20, units = "cm", dpi = 300, type = "cairo")


### ------------------------------------------------------------------------ ###
### compare OM SSB/F with MP SSB/F ####
### ------------------------------------------------------------------------ ###

df <- foreach(OM = c("Baseline", "Alt1", "Alt2"),
              .combine = rbind) %do% {
                res <- readRDS(paste0("output/runs/had4/1000_20/", OM, "_HCR-A_Ftrgt-0.28",
                                      "_Btrigger-180000_TACconstr-FALSE_BB-FALSE.rds"))
                if(OM %in% "Alt1"){
                  res_input <- readRDS(paste0("input/had4/999_20/MPbase_OM",OM,"_999.rds"))$om@stock
                  }else{
                  res_input <- readRDS(paste0("input/had4/1000_20/MPbase_OM",OM,"_1000.rds"))$om@stock
                }
                tmp <- FLQuant(NA, dimnames = list(
                  metric = c("F.om", "SSB.om", "F.est", "SSB.est", "F", "SSB"), 
                  year = dimnames(res_input)$year,
                  iter = dimnames(res_input)$iter))
                tmp["F.est", ac(an(dimnames(res@tracking)$year) - 1)] <- 
                  res@tracking["F.est"]
                tmp["SSB.est", ac(an(dimnames(res@tracking)$year) - 1)] <- 
                  res@tracking["B.est"]
                tmp["SSB.om", dimnames(res_input)$year] <- ssb(res_input)
                tmp["F.om", dimnames(res_input)$year] <- fbar(res_input)
                tmp["SSB.om", dimnames(res@stock)$year] <- ssb(res@stock)
                tmp["F.om", dimnames(res@stock)$year] <- fbar(res@stock)
                tmp["SSB"] <- tmp["SSB.est"] / tmp["SSB.om"]
                tmp["F"] <- tmp["F.est"] / tmp["F.om"]
                tmp_out <- cbind(as.data.frame(tmp), OM = OM)
                tmp_out <- tmp_out[!is.na(tmp_out$data) & is.finite(tmp_out$data), ]
                tmp_out
                
              }
df <- df %>% group_by(metric, year, OM) %>%
  summarise(X0.05 = quantile(data, probs = 0.05),
            X0.25 = quantile(data, probs = 0.25),
            X0.50 = quantile(data, probs = 0.50),
            X0.75 = quantile(data, probs = 0.75),
            X0.95 = quantile(data, probs = 0.95))
ggplot(data = df %>% filter(metric %in% c("F", "SSB")),
       aes(x = year, y = X0.50)) +
  geom_ribbon(aes(ymin = X0.05, ymax = X0.95), alpha = 0.3, fill = "#F8766D",
              show.legend = FALSE) +
  geom_ribbon(aes(ymin = X0.25, ymax = X0.75), alpha = 0.6, fill = "#F8766D",
              show.legend = FALSE) +
  geom_line(show.legend = FALSE) +
  facet_grid(OM ~ metric) +
  theme_bw() +
  geom_hline(yintercept = 1, alpha = 0.5) +
  labs(y = "MP value / OM value")
ggsave(filename = paste0("output/runs/had4/1000_20/plots/altOMs_stats/", 
                         "MP_vs_OM.png"), 
       width = 20, height = 20, units = "cm", dpi = 300, type = "cairo")


### ------------------------------------------------------------------------ ###
### OM_alt4: very low recruitment & stock below Blim ####
### ------------------------------------------------------------------------ ###

A <- readRDS("output/runs/cod4/1000_30/cod4_alt4_HCR-A_Ftrgt-0.38_Btrigger-170000_TACconstr-FALSE_BB-FALSE.rds")
plot(A@stock)
B <- readRDS("output/runs/cod4/1000_30/cod4_alt4_HCR-B_Ftrgt-0.38_Btrigger-160000_TACconstr-FALSE_BB-FALSE.rds")
plot(B@stock)
C <- readRDS("output/runs/cod4/1000_30/cod4_alt4_HCR-C_Ftrgt-0.38_Btrigger-170000_TACconstr-FALSE_BB-FALSE.rds")
AD <- readRDS("output/runs/cod4/1000_30/cod4_alt4_HCR-A_Ftrgt-0.4_Btrigger-190000_TACconstr-TRUE_BB-TRUE.rds")
BE <- AD <- readRDS("output/runs/cod4/1000_30/cod4_alt4_HCR-B_Ftrgt-0.36_Btrigger-130000_TACconstr-TRUE_BB-TRUE.rds")
CE <- readRDS("output/runs/cod4/1000_30/cod4_alt4_HCR-C_Ftrgt-0.36_Btrigger-140000_TACconstr-TRUE_BB-TRUE.rds")
OM_alt4 <- FLStocks(A = A@stock, B = B@stock, C = C@stock,
                    AD = AD@stock, BE = BE@stock, CE = CE@stock)
plot(OM_alt4[1:3])
ggsave(filename = paste0("output/runs/cod4/1000_20/plots/altOMs_stats/", 
                         "OM_alt4_ABC_stock.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")


### find minimum SSB
sapply(OM_alt4, function(x) which.min(iterMedians(ssb(x))))
dimnames(ssb(OM_alt4$A))$year[which.min(iterMedians(ssb(OM_alt4$A)))]
### 2034
df_4 <- as.data.frame(FLQuants(lapply(window(OM_alt4[1:3], 
                                             start = 2034, end = 2034),
                                      ssb)))

ggplot(data = df_4,
       aes(x = qname, y = data, group = qname, colour = qname)) +
  geom_boxplot() +
  theme_bw() +
  labs(y = "SSB [t]", x = "HCR option") +
  scale_color_discrete("")
ggsave(filename = paste0("output/runs/cod4/1000_20/plots/altOMs_stats/", 
                         "OM_alt4_ABC_lowest_SSB.png"), 
       width = 20, height = 20, units = "cm", dpi = 300, type = "cairo")
### recovery time
ssb_4 <- lapply(window(OM_alt4[1:3], start = 2034),
                ssb)
ssb_4 <- lapply(ssb_4, function(xi) {
  FLQuant(xi >= 107000)
})
ssb_4 <- foreach(i = seq_along(ssb_4), .combine = "rbind") %do% {
  data.frame(
    data = apply(ssb_4[[i]]@.Data, 6, function(x) {
      c(if (any(as.logical(x))) {which(as.logical(x))[1]} else {Inf})
    }),
    HCR = names(ssb_4)[i])
}
ggplot(data = ssb_4,
       aes(x = HCR, y = data, colour = HCR)) +
  geom_boxplot() +
  theme_bw() +
  labs(y = "recovery time to Blim [years]")
ggsave(filename = paste0("output/runs/cod4/1000_20/plots/altOMs_stats/", 
                         "OM_alt4_ABC_recovery.png"), 
       width = 20, height = 20, units = "cm", dpi = 300, type = "cairo")


### ------------------------------------------------------------------------ ###
### base OM optimised options + TAC constraint, no BB ####
### ------------------------------------------------------------------------ ###

### select maximum yield combinations
combs <- data.frame(name = rep(c("A", "B", "C", "AD", "BE", "CE"), each = 2),
                    HCR = rep(c("A", "B", "C", "A", "B", "C"), each = 2),
                    BB = c(rep(FALSE, 6), rep(c(TRUE, FALSE), 3)),
                    TACconstr = c(rep(c(FALSE, TRUE), 3), 
                                  rep(c(TRUE, TRUE), 3)),
                    Btrigger = rep(c(170000, 160000, 170000, 190000, 130000,
                                     140000), each = 2),
                    Ftrgt = rep(c(0.38, 0.38, 0.38,
                                  0.40, 0.36, 0.36), each = 2),
                    scenario = c("default", "TAC constraint"),
                    OM = "cod4")
combs <- merge(combs, stats)
combs_dat <- stats_full(data = combs)
combs_dat$scenario <- factor(combs_dat$scenario, 
                             levels = c("default", "TAC constraint"))
ggplot(data = combs_dat, 
       mapping = aes(x = name, y = value, group = interaction(scenario, name), 
                     colour = scenario)) +
  geom_boxplot() + 
  facet_wrap(~ key, scales = "free_y") +
  theme_bw() +
  ylim(0, NA)

p_catch_long <- ggplot(data = combs_dat[combs_dat$key == "catch_long", ], 
                       mapping = aes(x = name, y = value, colour = scenario)) +
  geom_boxplot(show.legend = FALSE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "long-term catch [t]")
p_catch_medium <- ggplot(data = combs_dat[combs_dat$key == "catch_medium", ], 
                         mapping = aes(x = name, y = value, colour = scenario)) +
  geom_boxplot(show.legend = FALSE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) + 
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "medium-term catch [t]")
p_catch_short <- ggplot(data = combs_dat[combs_dat$key == "catch_short", ], 
                        mapping = aes(x = name, y = value, colour = scenario)) +
  geom_boxplot(show.legend = FALSE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) + 
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "short-term catch [t]")
p_risk1_long <- ggplot(data = combs_dat[combs_dat$key == "risk1_long", ], 
                       mapping = aes(x = name, y = value, fill = scenario,
                                     colour = scenario)) +
  geom_bar(show.legend = FALSE, stat = "identity", 
           position = position_dodge2(preserve = "single")) + 
  geom_blank(data = combs_dat[combs_dat$key == "risk3_long", ]) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "long-term risk 1")
p_risk1_medium <- ggplot(data = combs_dat[combs_dat$key == "risk1_medium", ], 
                         mapping = aes(x = name, y = value, fill = scenario,
                                       colour = scenario)) +
  geom_bar(show.legend = FALSE, stat = "identity", 
           position = position_dodge2(preserve = "single")) + 
  geom_blank(data = combs_dat[combs_dat$key == "risk3_medium", ]) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "medium-term risk 1")
p_risk1_short <- ggplot(data = combs_dat[combs_dat$key == "risk1_short", ], 
                        mapping = aes(x = name, y = value, fill = scenario,
                                      colour = scenario)) +
  geom_bar(show.legend = FALSE, stat = "identity", 
           position = position_dodge2(preserve = "single")) + 
  geom_blank(data = combs_dat[combs_dat$key == "risk3_short", ]) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "short-term risk 1")
p_risk3_long <- ggplot(data = combs_dat[combs_dat$key == "risk3_long", ], 
                       mapping = aes(x = name, y = value, fill = scenario,
                                     colour = scenario)) +
  geom_bar(show.legend = FALSE, stat = "identity", 
           position = position_dodge2(preserve = "single")) + 
  geom_blank(data = combs_dat[combs_dat$key == "risk1_long", ]) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "long-term risk 3")
p_risk3_medium <- ggplot(data = combs_dat[combs_dat$key == "risk3_medium", ], 
                         mapping = aes(x = name, y = value, fill = scenario,
                                       colour = scenario)) +
  geom_bar(show.legend = FALSE, stat = "identity", 
           position = position_dodge2(preserve = "single")) + 
  geom_blank(data = combs_dat[combs_dat$key == "risk1_medium", ]) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "medium-term risk 3")
p_risk3_short <- ggplot(data = combs_dat[combs_dat$key == "risk3_short", ], 
                        mapping = aes(x = name, y = value, fill = scenario,
                                      colour = scenario)) +
  geom_bar(show.legend = FALSE, stat = "identity", 
           position = position_dodge2(preserve = "single")) + 
  geom_blank(data = combs_dat[combs_dat$key == "risk1_short", ]) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "short-term risk 3")
p_iav_long <- ggplot(data = combs_dat[combs_dat$key == "iav_long", ], 
                     mapping = aes(x = name, y = value, colour = scenario)) +
  geom_boxplot(show.legend = FALSE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "long-term inter-annual catch variability")
p_iav_medium <- ggplot(data = combs_dat[combs_dat$key == "iav_medium", ], 
                       mapping = aes(x = name, y = value, colour = scenario)) +
  geom_boxplot(show.legend = FALSE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "medium-term inter-annual catch variability")
p_iav_short <- ggplot(data = combs_dat[combs_dat$key == "iav_short", ], 
                      mapping = aes(x = name, y = value, colour = scenario)) +
  geom_boxplot(show.legend = FALSE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) +
  coord_cartesian(ylim = c(0, 1.5)) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "short-term inter-annual catch variability")
p_ssb_long <- ggplot(data = combs_dat[combs_dat$key == "ssb_long", ], 
                     mapping = aes(x = name, y = value, colour = scenario)) +
  geom_boxplot(show.legend = TRUE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) +
  coord_cartesian(ylim = c(0, 3e+5)) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "long-term SSB [t]") +
  theme(legend.direction = "horizontal")
p_ssb_medium <- ggplot(data = combs_dat[combs_dat$key == "ssb_medium", ], 
                       mapping = aes(x = name, y = value, colour = scenario)) +
  geom_boxplot(show.legend = TRUE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) + 
  coord_cartesian(ylim = c(0, 3e+5)) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "medium-term SSB [t]") +
  theme(legend.direction = "horizontal")
p_ssb_short <- ggplot(data = combs_dat[combs_dat$key == "ssb_short", ], 
                      mapping = aes(x = name, y = value, colour = scenario)) +
  geom_boxplot(show.legend = TRUE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) +
  coord_cartesian(ylim = c(0, 3e+5)) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "short-term SSB [t]") +
  theme(legend.direction = "horizontal") +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE, digits = 1),
                     limits = c(0, NA))
p_recovery_proportion <- 
  ggplot(data = combs_dat[combs_dat$key == "recovery_proportion", ], 
         mapping = aes(x = name, y = value, fill = scenario, colour = scenario)) +
  geom_bar(show.legend = FALSE, stat = "identity",
           position = position_dodge2(preserve = "single")) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "recovery proportion")
p_recovery_time <- 
  ggplot(data = combs_dat[combs_dat$key == "recovery_time", ], 
         mapping = aes(x = name, y = value, colour = scenario)) +
  geom_boxplot(show.legend = TRUE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "recovery time [years]")

plot_grid(plot_grid(p_catch_long, p_risk1_long, p_risk3_long, p_iav_long,
                    p_ssb_long + theme(legend.position = "none"),
                    align = "hv"),
          get_legend(p_ssb_long), nrow = 2, rel_heights = c(1, 0.1))
ggsave(filename = paste0("output/runs/cod4/1000_20/plots/baseOM_TAC_constraint/", 
                         "baseOM_combs_long.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
plot_grid(plot_grid(p_catch_medium, p_risk1_medium, p_risk3_medium, p_iav_medium,
                    p_ssb_medium + theme(legend.position = "none"),
                    align = "hv"),
          get_legend(p_ssb_medium), nrow = 2, rel_heights = c(1, 0.1))
ggsave(filename = paste0("output/runs/cod4/1000_20/plots/baseOM_TAC_constraint/", 
                         "baseOM_combs_medium.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
plot_grid(plot_grid(p_catch_short, p_risk1_short, p_risk3_short, p_iav_short,
                    p_ssb_short + theme(legend.position = "none"),
                    align = "hv"),
          get_legend(p_ssb_short), nrow = 2, rel_heights = c(1, 0.1))
ggsave(filename = paste0("output/runs/cod4/1000_20/plots/baseOM_TAC_constraint/", 
                         "baseOM_combs_short.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
plot_grid(plot_grid(p_recovery_proportion, 
                    p_recovery_time + theme(legend.position = "none"),
                    align = "hv"),
          get_legend(p_recovery_time), ncol = 2, rel_widths = c(0.5, 0.1))
ggsave(filename = paste0("output/runs/cod4/1000_20/plots/baseOM_TAC_constraint/", 
                         "baseOM_combs_recovery.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")


### ------------------------------------------------------------------------ ###
### plot risk over time ####
### ------------------------------------------------------------------------ ###

### get optimized A
stkA_file <- stats %>% filter(OM == "Baseline" & Ftrgt == 0.28 & Btrigger == 180000 &
                                TACconstr == FALSE & BB == FALSE & HCR == "A")
### get simulated SSB
ssbA_new <- ssb(readRDS(paste0("output/runs/had4/1000_20/", 
                               stkA_file$file))@stock)
### get historical SSB
ssbA <- ssb(readRDS("input/had4/1000_20/MPbase_OMBaseline_1000.rds")$om@stock)
### combine
ssbA[, dimnames(ssbA_new)$year] <- ssbA_new
### calculate annual risk
riskA <- apply((ssbA < stkA_file$Blim), 2, mean)
### plot
ggplot(data = as.data.frame(window(riskA, start = 2018)), 
       aes(x = year , y = data)) +
  geom_line() +
  theme_bw() +
  labs(x = "year", y = "p(SSB<Blim)") +
  geom_vline(data = data.frame(x = c(2018.5, 2023.5, 2028.5)),
             aes(xintercept = x), linetype = "dashed")
ggsave(filename = paste0("output/runs/had4/1000_20/plots/stock_plots/", 
                         "Baseline_risk_A_optimized.png"), 
       width = 15, height = 10, units = "cm", dpi = 300, type = "cairo")



#Alt1
### get optimized A
stkA_file <- stats %>% filter(OM == "Alt1" & Ftrgt == 0.28 & Btrigger == 180000 &
                                TACconstr == FALSE & BB == FALSE & HCR == "A")
### get simulated SSB
ssbA_new <- ssb(readRDS(paste0("output/runs/had4/1000_20/", 
                               stkA_file$file))@stock)
### get historical SSB
ssbA <- ssb(readRDS("input/had4/999_20/MPbase_OMAlt1_999.rds")$om@stock)
### combine
ssbA[, dimnames(ssbA_new)$year] <- ssbA_new
### calculate annual risk
riskA <- apply((ssbA < stkA_file$Blim), 2, mean)
### plot
ggplot(data = as.data.frame(window(riskA, start = 2018)), 
       aes(x = year , y = data)) +
  geom_line() +
  theme_bw() +
  labs(x = "year", y = "p(SSB<Blim)") +
  geom_vline(data = data.frame(x = c(2018.5, 2023.5, 2028.5)),
             aes(xintercept = x), linetype = "dashed")
ggsave(filename = paste0("output/runs/had4/1000_20/plots/stock_plots/", 
                         "Alt1_risk_A_optimized.png"), 
       width = 15, height = 10, units = "cm", dpi = 300, type = "cairo")

### get optimized A
stkA_file <- stats %>% filter(OM == "Alt2" & Ftrgt == 0.28 & Btrigger == 180000 &
                                TACconstr == FALSE & BB == FALSE & HCR == "A")
### get simulated SSB
ssbA_new <- ssb(readRDS(paste0("output/runs/had4/1000_20/", 
                               stkA_file$file))@stock)
### get historical SSB
ssbA <- ssb(readRDS("input/had4/1000_20/MPbase_OMAlt2_1000.rds")$om@stock)
### combine
ssbA[, dimnames(ssbA_new)$year] <- ssbA_new
### calculate annual risk
riskA <- apply((ssbA < stkA_file$Blim), 2, mean)
### plot
ggplot(data = as.data.frame(window(riskA, start = 2018)), 
       aes(x = year , y = data)) +
  geom_line() +
  theme_bw() +
  labs(x = "year", y = "p(SSB<Blim)") +
  geom_vline(data = data.frame(x = c(2018.5, 2023.5, 2028.5)),
             aes(xintercept = x), linetype = "dashed")
ggsave(filename = paste0("output/runs/had4/1000_20/plots/stock_plots/", 
                         "Alt2_risk_A_optimized.png"), 
       width = 15, height = 10, units = "cm", dpi = 300, type = "cairo")