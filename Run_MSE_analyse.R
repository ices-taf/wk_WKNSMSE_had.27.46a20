### ------------------------------------------------------------------------ ###
### process results ####
### ------------------------------------------------------------------------ ###
library(FLCore)
library(ggplot2)
library(tidyr)
library(cowplot)
library(dplyr)

library(doParallel)
cl <- makeCluster(2)
registerDoParallel(cl)

setwd(paste("/home/coleh/WKNSMSE/wk_WKNSMSE_had.27.46a20", sep=""))
### load additional functions
source("a4a_mse_WKNSMSE_funs.R")

### reference points
refpts <- list(msyBtrigger = 132000,
                   Blim=94000,
                   Flim=0.38,
                   Fmsy = 0.194,
                   Fpa = 0.274,
                   Bpa = 132000)

### ------------------------------------------------------------------------ ###
### get files ####
### ------------------------------------------------------------------------ ###

path_res <- "output/had4/runs/Baseline/999_20/"
files_res <- data.frame(file = list.files(path_res, pattern = "*.rds"), 
                        stringsAsFactors = FALSE)
files_res <- files_res[files_res$file != "stats.rds",, drop = FALSE]

files_res$Ftrgt <- as.numeric(lapply(lapply(strsplit(files_res$file, 
                                                     split = "_", fixed = TRUE),
                                            "[[", 2), gsub,
                                     pattern = "Ftrgt-", replacement = ""))
files_res$Btrigger <- as.numeric(lapply(lapply(strsplit(files_res$file, 
                                                        split = "_", fixed = TRUE),
                                               "[[", 3), gsub,
                                        pattern = "Btrigger-", replacement = ""))
files_res$HCR <- sapply(lapply(strsplit(files_res$file, 
                                        split = "_", fixed = TRUE),
                               "[[", 1), gsub,
                        pattern = "HCR-", replacement = "")
files_res$TACconstr <- as.logical(sapply(lapply(strsplit(files_res$file, 
                                                         split = "_", fixed = TRUE),
                                                "[[", 4), gsub,
                                         pattern = "TACconstr-", replacement = ""))
files_res$BB <- as.logical(
  sapply(lapply(lapply(strsplit(files_res$file, split = "_", fixed = TRUE), 
                       "[[", 5), gsub, pattern = "BB-", replacement = ""), 
         gsub, pattern = ".rds", replacement = ""))

#stats <- readRDS(paste0(path_res, "stats.rds"))
stats <- read.csv(paste0(path_res, "stats.csv"))
stats_new <- merge(stats, files_res, all = TRUE)
stats_new <- stats_new[!stats_new$file %in% stats$file, ]

res_list <- vector(mode = "list", length = nrow(stats_new))
res_list <- foreach(i = seq(nrow(stats_new))) %do% { #%dopar% {
  readRDS(paste0(path_res, stats_new$file[i]))
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
stats_new$risk1_full <- foreach(x = res_list, .packages = "FLCore",
                                .combine = "c") %dopar% {
                                  mean(window(ssb(x@stock), start = 2019) < 107000)
                                }
stats_new$risk1_long <- foreach(x = res_list, .packages = "FLCore",
                                .combine = "c") %dopar% {
                                  mean(window(ssb(x@stock), start = 2029) < refpts$Blim)
                                }
stats_new$risk1_short <- foreach(x = res_list, .packages = "FLCore",
                                 .combine = "c") %dopar% {
                                   mean(window(ssb(x@stock), start = 2019, end = 2023) < refpts$Blim)
                                 }
stats_new$risk1_medium <- foreach(x = res_list, .packages = "FLCore",
                                  .combine = "c") %dopar% {
                                    mean(window(ssb(x@stock), start = 2024, end = 2028) < refpts$Blim)
                                  }
stats_new$risk3_long <- foreach(x = res_list, .packages = "FLCore",
                                .combine = "c") %dopar% {
                                  max(iterMeans(window(ssb(x@stock), start = 2029) < refpts$Blim))
                                }
stats_new$risk3_short <- foreach(x = res_list, .packages = "FLCore",
                                 .combine = "c") %dopar% {
                                   max(iterMeans(window(ssb(x@stock), start = 2019, end = 2023) < refpts$Blim))
                                 }
stats_new$risk3_medium <- foreach(x = res_list, .packages = "FLCore",
                                  .combine = "c") %dopar% {
                                    max(iterMeans(window(ssb(x@stock), start = 2024, end = 2028) < refpts$Blim))
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
### SAM convergence
stats_new$conv_failed <- foreach(x = res_list, .packages = "FLCore",
                                 .combine = "c") %dopar% {
                                   sum(x@tracking["conv.est", ac(2018:2037)] != 0)
                                 }
all(stats_new$conv_failed == 0)
### check F maxed (2)
stats_new$F_maxed <- foreach(x = res_list, .packages = "FLCore",
                             .combine = "c") %dopar% {
                               sum(window(fbar(stock(x)), start = 2019) >= 2)
                             }
pos <- which(stats_new$F_maxed != 0)
stats_new[pos, ]
if(length(pos)>0){
plot(stock(res_list[[pos]]))
pos_iter <- which(apply(fbar(stock(res_list[[pos]])), c(1, 6), max) >= 2)
plot(stock(res_list[[pos]])[,,,,, pos_iter])
}
### F maxed in THREE scenarios, 
### in each of them in ONE iteration and ONE time only
### HCR B, Btrigger = 110000, Ftrgt = 0.5, iteration 29
### HCR AD, Btrigger = 170000, Ftrgt = 0.5, iteration 442
### HCR AD, Btrigger = 190000, Ftrgt = 0.5, iteration 442
### HCR BE, Btrigger = 190000, Ftrgt = 0.5, iteration 442
### HCR CE, Btrigger = 190000, Ftrgt = 0.5, iteration 442

stats <- rbind(stats, stats_new)
saveRDS(object = stats, file = paste0(path_res, "stats.rds"))
write.csv(x = stats, file = paste0(path_res, "stats.csv"), row.names = FALSE)

### ------------------------------------------------------------------------ ###
### plot ####
### ------------------------------------------------------------------------ ###


### plot function
grid <- function(dat, HCR = "A",
                 time = c("long", "short", "medium"),
                 add_risk1 = FALSE) {
  
  ### catch
  dat$catch <- dat[, paste0("catch_median_", time)]
  dat$risk <- dat[, paste0("risk3_", time)]
  dat$iav <- dat[, paste0("iav_", time)]
  dat$risk1 <- dat[, paste0("risk1_", time)]
  p1 <- ggplot(data = dat, 
               aes(x = Btrigger, y = Ftrgt, fill = catch)) +
    geom_raster() +
    scale_fill_gradient(paste0(time, "-term\ncatch (median)"), low = "red", 
                        high = "green") +
    geom_text(aes(label = round(catch), colour = risk <= 0.05),
              size = 2) +
    scale_colour_manual("risk <= 0.05", 
                        values = c("FALSE" = "red", "TRUE" = "black")) +
    theme_bw()
  ### risk
  p2 <- ggplot(data = dat, 
               aes(x = Btrigger, y = Ftrgt, fill = risk)) +
    geom_raster(alpha = 0.75) +
    scale_fill_gradient(paste0(time, "-term\nrisk 3"), 
                        low = "green", high = "red") +
    geom_text(aes(label = round(risk, 3), colour = risk <= 0.05),
              size = 2) +
    scale_colour_manual("risk <= 0.05", 
                        values = c("FALSE" = "red", "TRUE" = "black")) +
    theme_bw()
  ### iav
  p3 <- ggplot(data = dat, 
               aes(x = Btrigger, y = Ftrgt, fill = iav)) +
    geom_raster() +
    geom_text(aes(label = round(iav, 3), colour = risk <= 0.05),
              size = 2) +
    scale_colour_manual("risk <= 0.05",
                        values = c("FALSE" = "red", "TRUE" = "black")) +
    scale_fill_gradient(paste0(time,"-term\ninter-annual\nvariability of catch",
                               "\nmedian"),
                        low = "green", high = "red") +
    theme_bw()
  ### risk1
  p4 <- ggplot(data = dat, 
               aes(x = Btrigger, y = Ftrgt, fill = risk1)) +
    geom_raster(alpha = 0.75) +
    scale_fill_gradient(paste0(time, "-term\nrisk 1"), 
                        low = "green", high = "red") +
    geom_text(aes(label = round(risk1, 3), colour = risk1 <= 0.05),
              size = 2) +
    scale_colour_manual("risk <= 0.05", 
                        values = c("FALSE" = "red", "TRUE" = "black")) +
    theme_bw()
  
  if (isTRUE(add_risk1)) {
    plot_grid(p1, p2, p3, p4, nrow = 2, ncol = 2, align = "hv")
  } else {
    plot_grid(p1, p2, p3, nrow = 2, ncol = 2, align = "hv")
  }
  
  
}

### A
grid(dat = stats %>%
       filter(HCR == "A" & BB == FALSE & TACconstr == FALSE) %>%
       filter(Ftrgt %in% round(seq(0, 1, 0.01), 2)), 
     HCR = "A", time = "long", add_risk1 = FALSE)
ggsave(filename = "output/had4/runs/Baseline/999_20/grid_A_long.png", 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
grid(dat = stats %>%
       filter(HCR == "A" & BB == FALSE & TACconstr == FALSE) %>%
       filter(Ftrgt %in% round(seq(0, 1, 0.01), 2)), 
     HCR = "A", time = "medium")
ggsave(filename = "output/had4/runs/Baseline/999_20/grid_A_medium.png", 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
grid(dat = stats %>%
       filter(HCR == "A" & BB == FALSE & TACconstr == FALSE) %>%
       filter(Ftrgt %in% round(seq(0, 1, 0.01), 2)), 
     HCR = "A", time = "short")
ggsave(filename = "output/had4/runs/Baseline/999_20/grid_A_short.png", 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")

### B
grid(dat = stats %>%
       filter(HCR == "B" & BB == FALSE & TACconstr == FALSE) %>%
       filter(Ftrgt %in% round(seq(0, 1, 0.01), 2)), 
     HCR = "B", time = "long", add_risk1 = FALSE)
ggsave(filename = "output/had4/runs/1000_20/grid_B_long.png", 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
grid(dat = stats %>%
       filter(HCR == "B" & BB == FALSE & TACconstr == FALSE) %>%
       filter(Ftrgt %in% round(seq(0, 1, 0.01), 2)), 
     HCR = "B", time = "medium")
ggsave(filename = "output/had4/runs/1000_20/grid_B_medium.png", 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
grid(dat = stats %>%
       filter(HCR == "B" & BB == FALSE & TACconstr == FALSE) %>%
       filter(Ftrgt %in% round(seq(0, 1, 0.01), 2)), 
     HCR = "B", time = "short")
ggsave(filename = "output/had4/runs/1000_20/grid_B_short.png", 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")

### C
grid(dat = stats %>%
       filter(HCR == "C" & BB == FALSE & TACconstr == FALSE) %>%
       filter(Ftrgt %in% round(seq(0, 1, 0.01), 2)), 
     HCR = "C", time = "long", add_risk1 = FALSE)
ggsave(filename = "output/had4/runs/1000_20/grid_C_long.png", 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
grid(dat = stats %>%
       filter(HCR == "C" & BB == FALSE & TACconstr == FALSE) %>%
       filter(Ftrgt %in% round(seq(0, 1, 0.01), 2)), 
     HCR = "C", time = "medium")
ggsave(filename = "output/had4/runs/1000_20/grid_C_medium.png", 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
grid(dat = stats %>%
       filter(HCR == "C" & BB == FALSE & TACconstr == FALSE) %>%
       filter(Ftrgt %in% round(seq(0, 1, 0.01), 2)), 
     HCR = "C", time = "short")
ggsave(filename = "output/had4/runs/1000_20/grid_C_short.png", 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")

### AD
grid(dat = stats %>%
       filter(HCR == "A" & BB == TRUE & TACconstr == TRUE) %>%
       filter(Ftrgt %in% round(seq(0, 1, 0.01), 2)), 
     HCR = "A", time = "long", add_risk1 = FALSE)
ggsave(filename = "output/runs/cod4/1000_20/grid_AD_long.png", 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
grid(dat = stats %>%
       filter(HCR == "A" & BB == TRUE & TACconstr == TRUE) %>%
       filter(Ftrgt %in% round(seq(0, 1, 0.01), 2)), 
     HCR = "A", time = "medium")
ggsave(filename = "output/had4/runs/1000_20/grid_AD_medium.png", 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
grid(dat = stats %>%
       filter(HCR == "A" & BB == TRUE & TACconstr == TRUE) %>%
       filter(Ftrgt %in% round(seq(0, 1, 0.01), 2)), 
     HCR = "A", time = "short")
ggsave(filename = "output/had4/runs/1000_20/grid_AD_short.png", 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")


### BE
grid(dat = stats %>%
       filter(HCR == "B" & BB == TRUE & TACconstr == TRUE) %>%
       filter(Ftrgt %in% round(seq(0, 1, 0.01), 2)), 
     HCR = "B", time = "long", add_risk1 = FALSE)
ggsave(filename = "output/had4/runs/1000_20/grid_BE_long.png", 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
grid(dat = stats %>%
       filter(HCR == "B" & BB == TRUE & TACconstr == TRUE) %>%
       filter(Ftrgt %in% round(seq(0, 1, 0.01), 2)), 
     HCR = "B", time = "medium")
ggsave(filename = "output/had4/runs/1000_20/grid_BE_medium.png", 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
grid(dat = stats %>%
       filter(HCR == "B" & BB == TRUE & TACconstr == TRUE) %>%
       filter(Ftrgt %in% round(seq(0, 1, 0.01), 2)), 
     HCR = "B", time = "short")
ggsave(filename = "output/had4/runs/1000_20/grid_BE_short.png", 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")


### CE
grid(dat = stats %>%
       filter(HCR == "C" & BB == TRUE & TACconstr == TRUE) %>%
       filter(Ftrgt %in% round(seq(0, 1, 0.01), 2)), 
     HCR = "C", time = "long", add_risk1 = FALSE)
ggsave(filename = "output/had4/runs/1000_20/grid_CE_long.png", 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
grid(dat = stats %>%
       filter(HCR == "C" & BB == TRUE & TACconstr == TRUE) %>%
       filter(Ftrgt %in% round(seq(0, 1, 0.01), 2)), 
     HCR = "C", time = "medium")
ggsave(filename = "output/had4/runs/1000_20/grid_CE_medium.png", 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
grid(dat = stats %>%
       filter(HCR == "C" & BB == TRUE & TACconstr == TRUE) %>%
       filter(Ftrgt %in% round(seq(0, 1, 0.01), 2)), 
     HCR = "C", time = "short")
ggsave(filename = "output/had4/runs/1000_20/grid_CE_short.png", 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")



### ------------------------------------------------------------------------ ###
### check runs ####
### ------------------------------------------------------------------------ ###

if(FALSE){
pos <- which(stats$Ftrgt == 0.1 & stats$Btrigger == 110000 & 
               stats$HCR == "A" & stats$TACconstr == FALSE)
stats[pos, ]
plot(FLStocks(A = res_list[[1]]@stock, stability = res_list[[162]]@stock))

### load historical values
input <- readRDS("input/cod4/1000_20/base_run.rds")
### load estimated catch numbers
catch_n <- readRDS("input/cod4/1000_20/catch_n.rds")
catch.n(input$om@stock)[dimnames(catch_n)$age, dimnames(catch_n)$year] <- 
  catch_n
catch(input$om@stock) <- computeCatch(input$om@stock)
### plot "optimum" A
### Ftrgt = 0.37 & Btrigger = 150,000
#posA <- which(stats$Ftrgt == 0.37 & stats$Btrigger == 150000 & stats$HCR == "A")
stkA <- input$om@stock
stkA_res <- readRDS(paste0("output/had4/runs/1000_20/HCR-A_Ftrgt-0.37_Btrigger",
                           "-150000_TACconstr-FALSE_BB-FALSE.rds"))@stock
stkA[, ac(2018:2038)] <- stkA_res[, ac(2018:2038)]
plot(stkA, probs = c(0.05, 0.25, 0.5, 0.75, 0.95)) +
  xlab("year") + geom_vline(xintercept = 2018.5) +
  geom_hline(data = data.frame(qname = "SSB", data = 107000),
             aes(yintercept = data), linetype = "dashed") +
  geom_hline(data = data.frame(qname = "SSB", data = 150000),
             aes(yintercept = data), linetype = "solid") +
  geom_hline(data = data.frame(qname = "F", data = 0.54),
             aes(yintercept = data), linetype = "dashed") +
  geom_hline(data = data.frame(qname = "F", data = 0.31),
             aes(yintercept = data), linetype = "solid") +
  theme_bw() +
  geom_blank(data = as.data.frame(FLQuants(`Rec` = rec(stkA),
                                           `SSB` = ssb(stkA),
                                           `Catch` = catch(stkA),
                                           `F` = fbar(stkA))), 
             aes(x = year, y = data, group = iter))
ggsave(filename = paste0("output/had4/runs/1000_20/", 
                         "stk_Ftrgt=0.37_Btrigger=150000_HCR=A.png"),
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
### iterations
stkA_df <- as.data.frame(FLQuants(`Rec [1000]` = rec(stkA),
                                  `SSB [t]` = ssb(stkA),
                                  `Catch [t]` = catch(stkA),
                                  `F` = fbar(stkA)))
ggplot(data = stkA_df[stkA_df$iter %in% 1:1000, ], 
       aes(x = year, y = data, group = iter)) +
  geom_line(alpha = 0.025) +
  facet_wrap(~ qname, ncol = 1, strip.position = "right",
             scale = "free_y") +
  theme_bw() +
  ylim(c(0, NA)) + labs(y = "") + 
  geom_vline(xintercept = 2018.5) +
  geom_hline(data = data.frame(qname = "SSB [t]", data = 107000),
             aes(yintercept = data), linetype = "dashed") +
  geom_hline(data = data.frame(qname = "SSB [t]", data = 150000),
             aes(yintercept = data), linetype = "solid") +
  geom_hline(data = data.frame(qname = "F", data = 0.54),
             aes(yintercept = data), linetype = "dashed") +
  geom_hline(data = data.frame(qname = "F", data = 0.31),
             aes(yintercept = data), linetype = "solid")
ggsave(filename = paste0("output/had4/runs/1000_20/", 
                         "stk_Ftrgt=0.37_Btrigger=150000_HCR=A_iter.png"),
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")

}

### ------------------------------------------------------------------------ ###
### F=0 plot ####
### ------------------------------------------------------------------------ ###
# 
# stkF0 <- readRDS(file = "input/had4/10000_100/data_F0.RData")$om@stock
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
# ggsave(filename = paste0("output/runs/cod4/1000_20/", 
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
# ggsave(filename = paste0("output/runs/cod4/1000_20/", 
#                          "stk_F0_10000iters_iters.png"),
#        width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
# 

### ------------------------------------------------------------------------ ###
### extrapolate 5% risk line ####
### ------------------------------------------------------------------------ ###

df_risks <- data.frame(Btrigger = c(seq(from = 120000, to = 190000, 
                                        by = 10000)),
                       Ftrgt = c(0.228, 0.24, 0.248, 0.2575, 0.27, 0.28,
                                 0.29, 0.30))

lm_risks <- lm(formula = Ftrgt ~ Btrigger, data = df_risks[1:4,])
plot(Ftrgt ~ Btrigger, data = df_risks, 
     xlim = c(110000, 200000), ylim = c(0.2, 0.3))
abline(lm_risks)