### ------------------------------------------------------------------------ ###
### process results ####
### ------------------------------------------------------------------------ ###
library(FLCore)
library(ggplot2)
library(tidyr)
library(cowplot)
library(dplyr)
library(ggplotFL)

setwd(paste("~/WKNSMSE/wk_WKNSMSE_had.27.46a20", sep=""))
### load additional functions
source("a4a_mse_WKNSMSE_funs.R")

library(doParallel)
cl <- makeCluster(1)
registerDoParallel(cl)


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

omName<-c("Baseline","Alt1")
n<-c(999,1000) # number of iters

### ------------------------------------------------------------------------ ###
### get files ####
### ------------------------------------------------------------------------ ###

path_res <- c(paste0("output/had4/runs/",omName[1],"/",n[1],"_20/"),paste0("output/had4/runs/",omName[2],"/",n[2],"_20/"))
path_out <- paste0("output/had4/res_plots/plots/")
path_mp<-"input/had4/MP_base/"
mp_file<-c(paste0("MPBase_OM",omName[1],"_",n[1],".rds"),paste0("MPBase_OM",omName[2],"_",n[2],".rds"))


stats <- readRDS(paste0("output/had4/res_plots/", "stats.rds"))
stats <- read.csv(paste0("output/had4/res_plots/", "stats.csv"))

# remove intyrTACcont
stats_orig<-stats
idx<-grep("intyrTACcont",stats$file)
stats<-stats[-idx,]

### ------------------------------------------------------------------------ ###
### plot functions ####
### ------------------------------------------------------------------------ ###

### ------------------------------------------------------------------------ ###
### plot some runs ####
### ------------------------------------------------------------------------ ###

### function for plotting
plot_stk <- function(stats, OM_ = "Baseline", HCR_ = "A", Ftrgt_ = 0.28,
                     Btrigger_ = 160000, TACconstr_ = FALSE, BB_ = FALSE,
                     probs = c(0.05, 0.25, 0.5, 0.75, 0.95),
                     path_out = "output/had4/runs/Baseline/999_20/",
                     path_in = "input/had4/MP_base/", file_in = "MPBase_OMBaseline_999.rds",
                     path_res = "output/had4/runs/Baseline/999_20/plots/",
                     path_catch=paste0("input/had4/999_20/Baseline_stk.rds"),
                     save_plot = TRUE,
                     get_history = TRUE, overwrite_catch_history = FALSE,
                     yr_start = 2018.5,
                     Blim = 94000, MSYBtrigger = 132000,
                     Flim = 0.38, Fmsy = 0.194,
                     plot_iter = FALSE, iters_plot=0,
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
      #catch_hist <- readRDS(paste0(path_in, "catch_n.rds"))
      catch_hist <- catch.n(readRDS(path_catch))
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
    filename <- paste0(path_res,OM_, gsub(x = stats_i$file, pattern = ".rds",
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
plot_stk(stats = stats, OM_ = "Baseline", HCR_ = "A", Ftrgt_ = 0.28, 
         Btrigger_ = 160000, TACconstr_ = FALSE, BB_ = FALSE, 
         overwrite_catch_history = F, Blim=refpts$Blim,
         MSYBtrigger=refpts$msyBtrigger, Flim=refpts$Flim, Fmsy=refpts$Fmsy,
         path_out=path_res, path_in=path_mp,path_res=path_out,file_in=mp_file,
         path_catch=paste0("input/had4/999_20/Baseline_stk.rds"))

### A, B, C, AD, BE, CE
combs <- data.frame(name = c("A*", "A", "B", "C", "AD", "BE", "CE"),
                    HCR = c("A", "A", "B", "C", "A", "B", "C"),
                    BB = c(rep(FALSE, 4), TRUE, TRUE, TRUE),
                    TACconstr = c(rep(FALSE, 4), TRUE, TRUE, TRUE),
                    Btrigger = c(132000, 160000, 160000,160000, 160000, 140000, 150000),
                    Ftrgt = c(0.194, 0.28, 0.29, 0.28, 0.27, 0.27, 0.26))
combs <- rbind(cbind(combs, OM = "Baseline"),
               cbind(combs, OM = "Alt1"))
combs <- combs[-c(8), ]

path_catch<-c(paste0("input/had4/999_20/Baseline_stk.rds"),paste0("input/had4/1000_20/Alt1_stk.rds"))
lapply(X = split(combs, seq(nrow(combs))), FUN = function(i) {
j<-ifelse(i$OM %in% "Baseline",1,2)
  plot_stk(stats = stats, OM_ = i$OM, HCR_ = i$HCR, Ftrgt_ = i$Ftrgt, 
           Btrigger_ = i$Btrigger, TACconstr_ = i$TACconstr, BB_ = i$BB, 
           overwrite_catch_history = F,
           Blim=refpts$Blim,
           MSYBtrigger=refpts$msyBtrigger, Flim=refpts$Flim, Fmsy=refpts$Fmsy,
           path_out=path_res[j], path_in=path_mp,path_res=path_out,file_in=mp_file[j],
           path_catch=path_catch[j])
})
### plot with first five iterations
lapply(X = split(combs, seq(nrow(combs))), FUN = function(i) {
j<-ifelse(i$OM %in% "Baseline",1,2)
  plot_stk(stats = stats, OM_ = i$OM, HCR_ = i$HCR, Ftrgt_ = i$Ftrgt, 
           Btrigger_ = i$Btrigger, TACconstr_ = i$TACconstr, BB_ = i$BB, 
           overwrite_catch_history = F, iters_plot = 2:6,
           Blim=refpts$Blim,
           MSYBtrigger=refpts$msyBtrigger, Flim=refpts$Flim, Fmsy=refpts$Fmsy,
           path_out=path_res[j], path_in=path_mp,path_res=path_out,file_in=mp_file[j],
           path_catch=path_catch[j])
})


 #path_in = switch(i$OM,"Baseline"=paste0("input/", i$OM, "/999_20/"),Alt1=paste0("input/", i$OM, "/1000_20/")))


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


### ------------------------------------------------------------------------ ###
### summary plots: compare HCR options ####
### ------------------------------------------------------------------------ ###

### select maximum yield combinations
combs <- data.frame(name = c("A*", "A", "B", "C", "AD", "BE", "CE"),
                    OM=c("Baseline"),
                    HCR = c("A", "A", "B", "C", "A", "B", "C"),
                    BB = c(rep(FALSE, 4), TRUE, TRUE, TRUE),
                    TACconstr = c(rep(FALSE, 4), TRUE, TRUE, TRUE),
                    Btrigger = c(132000, 160000, 160000,160000, 160000, 140000, 150000),
                    Ftrgt = c(0.194, 0.28, 0.29, 0.28, 0.27, 0.27, 0.26),
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
                        .packages = "FLCore", .combine = rbind) %dopar% {
                          
                          if(i$OM %in% "Baseline"){
                          stk_i <- readRDS(paste0("output/had4/runs/Baseline/999_20/", i$file))
                          }
                       if(i$OM %in% "Alt1"){
                          stk_i <- readRDS(paste0("output/had4/runs/Alt1/1000_20/", i$file))
                          }
                          
                          MSYBtrigger <-  132000 
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
ggsave(filename = paste0(path_out, 
                         "summary_baseOM_long.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
plot_grid(p_catch_medium, p_risk1_medium, p_risk3_medium, p_iav_medium, 
          p_ssb_medium,
          align = "hv")
ggsave(filename = paste0(path_out, 
                         "summary_baseOM_medium.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
plot_grid(p_catch_short, p_risk1_short, p_risk3_short, p_iav_short, p_ssb_short,
          align = "hv")
ggsave(filename = paste0(path_out, 
                         "summary_baseOM_short.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
plot_grid(p_recovery_proportion, p_recovery_time,
          align = "hv")
ggsave(filename = paste0(path_out,
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
                    Btrigger = rep(c(160000, 160000, 160000, 160000, 140000,
                                     150000), each = 5),
                    Ftrgt = c(0.28 * c(0.9, 1, 1.1), 0.167, 0.194,
                              0.29 * c(0.9, 1, 1.1), 0.167, 0.194,
                              0.28 * c(0.9, 1, 1.1), 0.167, 0.194,
                              0.27 * c(0.9, 1, 1.1), 0.167, 0.194,
                              0.27 * c(0.9, 1, 1.1), 0.167, 0.194,
                              0.26 * c(0.9, 1, 1.1), 0.167, 0.194), 
                    scenario = c("0.9*Ftrgt", "Ftrgt", "1.1*Ftrgt",
                                 "Fmsylower", "Fmsyupper"),
                    OM = "Baseline")
combs<-rbind(combs,data.frame(name = rep(c("A", "AD"), each = 2),
           HCR = rep(c("A"), each = 4),
           BB = rep(c(rep(FALSE, 1), rep(TRUE, 1)), each = 2),
           TACconstr = rep(c(rep(FALSE, 1), rep(TRUE, 1)), each = 2),
           Btrigger = c(160000*c(1.5,2),
                        160000*c(1.5,2)),
           Ftrgt = rep(c(0.28, 0.27), each = 2), 
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
ggsave(filename = paste0(path_out, 
                         "baseOM_combs_long.png"), 
       width = 30, height = 20, units = "cm", dpi = 250, type = "cairo")
plot_grid(plot_grid(p_catch_medium, p_risk1_medium, p_risk3_medium, p_iav_medium,
                    p_ssb_medium + theme(legend.position = "none"),
                    align = "hv"),
          get_legend(p_ssb_medium), nrow = 2, rel_heights = c(1, 0.1))
ggsave(filename = paste0(path_out,#"output/runs/cod4/1000_20/plots/baseOM_stats_combs/", 
                         "baseOM_combs_medium.png"), 
       width = 30, height = 20, units = "cm", dpi = 200, type = "cairo")
plot_grid(plot_grid(p_catch_short, p_risk1_short, p_risk3_short, p_iav_short,
                    p_ssb_short + theme(legend.position = "none"),
                    align = "hv"),
          get_legend(p_ssb_short), nrow = 2, rel_heights = c(1, 0.1))
ggsave(filename = paste0(path_out,#"output/runs/cod4/1000_20/plots/baseOM_stats_combs/", 
                         "baseOM_combs_short.png"), 
       width = 30, height = 20, units = "cm", dpi = 200, type = "cairo")
plot_grid(plot_grid(p_recovery_proportion, 
                    p_recovery_time + theme(legend.position = "none"),
                    align = "hv"),
          get_legend(p_recovery_time), ncol = 2, rel_widths = c(0.5, 0.1))
ggsave(filename = paste0(path_out,#"output/runs/cod4/1000_20/plots/baseOM_stats_combs/", 
                         "baseOM_combs_recovery.png"), 
       width = 30, height = 20, units = "cm", dpi = 200, type = "cairo")


### ------------------------------------------------------------------------ ###
### summary plots: compare alternative OMs ####
### ------------------------------------------------------------------------ ###

### alternative OMs
### select maximum yield combinations
combs_alt <- data.frame(name = c("A*", "A", "B", "C", "AD", "BE", "CE"),
                        HCR = c("A", "A", "B", "C", "A", "B", "C"),
                        BB = c(rep(FALSE, 4), TRUE, TRUE, TRUE),
                        TACconstr = c(rep(FALSE, 4), TRUE, TRUE, TRUE),
                        Btrigger = c(132000, 160000, 160000, 160000, 160000, 140000,
                                     150000),
                        Ftrgt = c(0.194, 0.28, 0.29, 0.28, 0.27, 0.27, 0.26),
                        scenario = 0)
combs_alt <- rbind(cbind(combs_alt, OM = "Baseline"),
                   cbind(combs_alt, OM = "Alt1"))
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
ggsave(filename = paste0(path_out, 
                         "summary_altOM_long.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
plot_grid(plot_grid(p_catch_medium, p_risk1_medium, p_risk3_medium, p_iav_medium,
                    p_ssb_medium + theme(legend.position = "none"),
                    align = "hv"),
          get_legend(p_ssb_medium), nrow = 2, rel_heights = c(1, 0.1))
ggsave(filename = paste0(path_out, 
                         "summary_altOM_medium.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
plot_grid(plot_grid(p_catch_short, p_risk1_short, p_risk3_short, p_iav_short,
                    p_ssb_short + theme(legend.position = "none"),
                    align = "hv"),
          get_legend(p_ssb_short), nrow = 2, rel_heights = c(1, 0.1))
ggsave(filename = paste0(path_out, 
                         "summary_altOM_short.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
plot_grid(plot_grid(p_recovery_proportion, 
                    p_recovery_time + theme(legend.position = "none"),
                    align = "hv"),
          get_legend(p_recovery_time), ncol = 2, rel_widths = c(0.5, 0.1))
ggsave(filename = paste0(path_out, 
                         "summary_altOM_recovery.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")

#########################################################################################
#########################################################################################

### get optimized A

stkA_file <- stats %>% filter(OM == "Baseline" & Ftrgt == 0.28 & Btrigger == 160000 &
                   TACconstr == FALSE & BB == FALSE & HCR == "A")

### get simulated SSB
ssbA_new <- ssb(readRDS(paste0("output/had4/runs/Baseline/999_20/", 
                               stkA_file$file))@stock)

### get historical SSB
ssbA <- ssb(readRDS(paste0(path_mp,mp_file[1]))$om@stock)

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

ggsave(filename = paste0(path_out, 
                         "risk_A_optimized.png"), 
       width = 15, height = 10, units = "cm", dpi = 300, type = "cairo")


### get optimized A

stkA_file <- stats %>% filter(OM == "Alt1" & Ftrgt == 0.28 & Btrigger == 160000 &
                   TACconstr == FALSE & BB == FALSE & HCR == "A")

### get simulated SSB
ssbA_new <- ssb(readRDS(paste0("output/had4/runs/Alt1/1000_20/", 
                               stkA_file$file))@stock)

### get historical SSB
ssbA <- ssb(readRDS(paste0(path_mp,mp_file[2]))$om@stock)

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

ggsave(filename = paste0(path_out, 
                         "risk_A_Alt1optimized.png"), 
       width = 15, height = 10, units = "cm", dpi = 300, type = "cairo")
       
ssbA<-apply(ssbA,2,median)
ggplot(data = as.data.frame(window(ssbA, start = 2018)), 
       aes(x = year , y = data)) +
  geom_line() +
  theme_bw() +
  labs(x = "year", y = "SSB") +
  geom_vline(data = data.frame(x = c(2018.5, 2023.5, 2028.5)),
             aes(xintercept = x), linetype = "dashed")

ggsave(filename = paste0(path_out, 
                         "ssb_A_optimized.png"), 
       width = 15, height = 10, units = "cm", dpi = 300, type = "cairo")


stkA_file <- stats %>% filter(OM == "Baseline" & Ftrgt == 0.28 & Btrigger == 160000 &
                   TACconstr == FALSE & BB == FALSE & HCR == "A")

### get simulated SSB
fbarA_new <- fbar(readRDS(paste0("output/had4/runs/Baseline/999_20/", 
                               stkA_file$file))@stock)

### get historical SSB
fbarA <- fbar(readRDS(paste0(path_mp,mp_file[1]))$om@stock)

### combine
fbarA[, dimnames(fbarA_new)$year] <- fbarA_new
fbarA<-apply(fbarA,2,median)

### plot
ggplot(data = as.data.frame(window(fbarA, start = 2018)), 
       aes(x = year , y = data)) +
  geom_line() +
  theme_bw() +
  labs(x = "year", y = "Fbar") +
  geom_vline(data = data.frame(x = c(2018.5, 2023.5, 2028.5)),
             aes(xintercept = x), linetype = "dashed")

ggsave(filename = paste0(path_out, 
                         "fbar_A_optimized.png"), 
       width = 15, height = 10, units = "cm", dpi = 300, type = "cairo")


### get simulated SSB
catA_new <- catch(readRDS(paste0("output/had4/runs/Baseline/999_20/", 
                               stkA_file$file))@stock)

### get historical SSB
catA <- catch(readRDS(paste0(path_mp,mp_file[1]))$om@stock)

### combine
catA[, dimnames(catA_new)$year] <- catA_new
catA<-apply(catA,2,median)

### plot
ggplot(data = as.data.frame(window(catA, start = 2018)), 
       aes(x = year , y = data)) +
  geom_line() +
  theme_bw() +
  labs(x = "year", y = "Catch") +
  geom_vline(data = data.frame(x = c(2018.5, 2023.5, 2028.5)),
             aes(xintercept = x), linetype = "dashed")

ggsave(filename = paste0(path_out, 
                         "catch_A_optimized.png"), 
       width = 15, height = 10, units = "cm", dpi = 300, type = "cairo")
       
       ### get simulated SSB
methcrA_new <- (readRDS(paste0("output/had4/runs/Baseline/999_20/", 
                               stkA_file$file))@tracking)["metric.hcr",]

tmp <- FLQuant(NA, dimnames = list(
    metric = c("metric.hcr"), 
    year = dimnames(methcrA_new)$year,
    iter = dimnames(methcrA_new)$iter))
  tmp["metric.hcr", ac(an(dimnames(methcrA_new)$year) )] <- methcrA_new 

    
methcrA<-apply(tmp,2,median,na.rm=T)

### plot
ggplot(data = as.data.frame(window(methcrA, start = 2018)), 
       aes(x = year , y = data)) +
  geom_line() +
  theme_bw() +
  labs(x = "year", y = "metric.hcr") +
  geom_vline(data = data.frame(x = c(2018.5, 2023.5, 2028.5)),
             aes(xintercept = x), linetype = "dashed")

ggsave(filename = paste0(path_out, 
                         "catch_A_optimized.png"), 
       width = 15, height = 10, units = "cm", dpi = 300, type = "cairo")
### ------------------------------------------------------------------------ ###

### compare OM SSB/F with MP SSB/F ####

### ------------------------------------------------------------------------ ###



df <- foreach(OM = c("Baseline", "Alt1"),
              .combine = rbind) %do% {
              j<-ifelse(OM %in% "Baseline",1,2)
  res <- readRDS(paste0(path_res[j],"HCR-A_Ftrgt-0.28",
                        "_Btrigger-160000_TACconstr-FALSE_BB-FALSE.rds"))

  res_input <- readRDS(paste0(path_mp,mp_file[j]))$om@stock
  
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
ggsave(filename = paste0(path_out, 
                         "MP_vs_OM.png"), 
       width = 20, height = 20, units = "cm", dpi = 300, type = "cairo")




###### look at tracking for higher risk in medium term 
# ICV is also high so what is driving it?


df <- foreach(OM = c("Baseline", "Alt1"),
              .combine = rbind) %do% {
              j<-ifelse(OM %in% "Baseline",1,2)
  res <- readRDS(paste0(path_res[j],"HCR-A_Ftrgt-0.28",
                        "_Btrigger-160000_TACconstr-FALSE_BB-FALSE.rds"))

res_input <- readRDS(paste0(path_mp,mp_file[j]))$om@stock
  
  tmp <- FLQuant(NA, dimnames = list(
    metric = c("metric.hcr", "metric.is","B.est"), 
    year = dimnames(res_input)$year,
    iter = dimnames(res_input)$iter))
  tmp["metric.hcr", ac(an(dimnames(res@tracking)$year) )] <- 
    res@tracking["metric.hcr"]
  tmp["metric.is", ac(an(dimnames(res@tracking)$year) )] <- 
    res@tracking["metric.is"]
      tmp["B.est", ac(an(dimnames(res@tracking)$year) )] <- 
    res@tracking["B.est"]
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


ggplot(data = df %>% filter(metric %in% c("metric.is")),
       aes(x = year, y = X0.50)) +
  geom_ribbon(aes(ymin = X0.05, ymax = X0.95), alpha = 0.3, fill = "#F8766D",
              show.legend = FALSE) +
  geom_ribbon(aes(ymin = X0.25, ymax = X0.75), alpha = 0.6, fill = "#F8766D",
              show.legend = FALSE) +
  geom_line(show.legend = FALSE) +
  facet_grid(OM ~ metric) +
  theme_bw()+
    labs(y = "Advised TAC")
  
ggplot(data = df %>% filter(metric %in% c("metric.hcr")),
       aes(x = year, y = X0.50)) +
  geom_ribbon(aes(ymin = X0.05, ymax = X0.95), alpha = 0.3, fill = "#F8766D",
              show.legend = FALSE) +
  geom_ribbon(aes(ymin = X0.25, ymax = X0.75), alpha = 0.6, fill = "#F8766D",
              show.legend = FALSE) +
  geom_line(show.legend = FALSE) +
  facet_grid(OM ~ metric) +
  theme_bw() +
    geom_hline(yintercept = 0.28, alpha = 0.5)+
  labs(y = "Ftrgt (hcr)")
    ggsave(filename = paste0(path_out, 
                         "Ftrgt_hcr.png"), 
       width = 20, height = 20, units = "cm", dpi = 300, type = "cairo")

  
ggplot(data = df %>% filter(metric %in% c("B.est")),
       aes(x = year, y = X0.50)) +
  geom_ribbon(aes(ymin = X0.05, ymax = X0.95), alpha = 0.3, fill = "#F8766D",
              show.legend = FALSE) +
  geom_ribbon(aes(ymin = X0.25, ymax = X0.75), alpha = 0.6, fill = "#F8766D",
              show.legend = FALSE) +
  geom_line(show.legend = FALSE) +
  facet_grid(OM ~ metric) +
  theme_bw() +
  geom_hline(yintercept = 160000, alpha = 0.5)+
  labs(y = "B.est")
  ggsave(filename = paste0(path_out, 
                         "B_est.png"), 
       width = 20, height = 20, units = "cm", dpi = 300, type = "cairo")

##########################################################################################
  
       #### look at why intyr changes things?
  
idx<-grep("intyrTACcont",stats_orig$file)

stats_int<-stats_orig[idx,]

plot_stk(stats = stats_int, OM_ = "Baseline", HCR_ = "A", Ftrgt_ = 0.28, 
         Btrigger_ = 160000, TACconstr_ = FALSE, BB_ = FALSE, 
         overwrite_catch_history = F, Blim=refpts$Blim,
         MSYBtrigger=refpts$msyBtrigger, Flim=refpts$Flim, Fmsy=refpts$Fmsy,
         path_out=path_res[1], path_in=path_mp[1],path_res=path_out[1],file_in=mp_file[1],
         path_catch=paste0("input/had4/999_20/Baseline_stk.rds"))

plot_stk(stats = stats_int, OM_ = "Baseline", HCR_ = "A", Ftrgt_ = 0.28, 
         Btrigger_ = 160000, TACconstr_ = FALSE, BB_ = FALSE, 
         overwrite_catch_history = F, Blim=refpts$Blim,iters_plot=2:6,
         MSYBtrigger=refpts$msyBtrigger, Flim=refpts$Flim, Fmsy=refpts$Fmsy,
         path_out=path_res[1], path_in=path_mp[1],path_res=path_out[1],file_in=mp_file[1],
         path_catch=paste0("input/had4/999_20/Baseline_stk.rds"))

stats_int<-rbind(stats_int,filter(stats,OM == "Baseline" & Ftrgt == 0.28 & Btrigger == 160000 &
                   TACconstr == FALSE & BB == FALSE & HCR == "A"))
                   
### select maximum yield combinations
combs <- data.frame(name = c("A_intyrTAC", "A"),
                    OM=c("Baseline"),
                    HCR = c("A", "A"),
                    BB = c(rep(FALSE,2)),
                    TACconstr = c(rep(FALSE, 2)),
                    Btrigger = c(160000, 160000),
                    Ftrgt = c( 0.28,0.28),
                    scenario=0)

combs <- merge(combs, stats_int, all.x = TRUE)[c(1,4),]
combs2 <- gather(data = combs, key = "key", value = "value",
                 catch_median_long, risk3_long, iav_long,
                 ssb_median_long, recovery_proportion, recovery_time)
combs2$name <- factor(combs2$name, levels = c("A_intyrTAC", "A"))

ggplot(data = combs2, 
       mapping = aes(x = name, y = value, group = name)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ key, scales = "free_y") +
  theme_bw()


combs_base <- stats_full(data = combs)
levels(combs_base$name)<-c(levels(combs_base$name),"A_intyrTAC")
combs_base$name[is.na(combs_base$name)]<-"A_intyrTAC"
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
                          filter(name == "A") %>%
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
ggsave(filename = paste0(path_out, 
                         "summary_AvsintyrTAC_long.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
plot_grid(p_catch_medium, p_risk1_medium, p_risk3_medium, p_iav_medium, 
          p_ssb_medium,
          align = "hv")
ggsave(filename = paste0(path_out, 
                         "summary_AvsintyrTAC_medium.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
plot_grid(p_catch_short, p_risk1_short, p_risk3_short, p_iav_short, p_ssb_short,
          align = "hv")
ggsave(filename = paste0(path_out, 
                         "summary_AvsintyrTAC_short.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
plot_grid(p_recovery_proportion, p_recovery_time,
          align = "hv")
ggsave(filename = paste0(path_out,
                         "summary_AvsintyrTAC_recovery.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")