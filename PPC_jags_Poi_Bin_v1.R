# PPC of zero counts for JAGS
# 

# Bring in the results from the JAGS run:

rm(list = ls())
library(tidyverse)
library(posterior)
library(bayesplot)
library(ggplot2)

data.path <- "data\\Formatted Annual Data Combined v3\\"
FILES <- list.files(path = data.path,
                    pattern = "Formatted_combined_v3.csv")

stan.fixed.year.run.date <- "2026-09-11"
fixed.year.stan.info <- readRDS(paste0("RData\\GWCalfCount_Stan_nb_fixedyear_",
                                       stan.fixed.year.run.date, ".info"))

ppc.out <- list()
i <- 15
for(i in 1:length(FILES)){
  
  year <- as.numeric(str_split(FILES[i], " Formatted_combined_v3.csv")[[1]][1])
  if (year < 2026){

    fixed.year.stan.info$pre.stan.data %>%
      filter(calendar_year == year) %>%
      mutate(count_per_hr = count_obs/effort_hours) %>%
      group_by(calendar_year) %>%
      mutate(cumu_effort = cumsum(effort_hours)) %>%
      ungroup()  -> data.y
    
    stan.data.count.obs <- data.y %>%
      filter(effort_hours > 0) %>% #-> tmp
      select(count_obs) %>% pull()
    
    out.file.name <- paste0("RData\\calf_estimates_v3",
                            "_Mv1_", year, ".rds") 
    
    jm.out <- readRDS(file = out.file.name)
    
    # Pull out jags input data and filter out effort == 0
    jags.data.df <- data.frame(effort = jm.out$jags.data$effort,
                               week = jm.out$jags.data$week,
                               count.obs = jm.out$jags.data$count.obs) %>%
      dplyr::filter(effort > 0)
    
    # Error checking for data - jags vs stan
    # if (abs(length(jags.data.df$count.obs) - length(stan.data.count.obs)) > 0)
    #   print(paste0("data lengths are not the same for year ", year))
    # 
    # if (mean(jags.data.df$count.obs == 0) != mean(stan.data.count.obs == 0))
    #   print(paste0("means of zero counts are not the same for year ", year))
    
    # Find all weeks with > 0 effort:
    unique.weeks <- unique(jags.data.df$week)
    
    # map each row's week to a column index in the subsetted lambda
    week_col <- match(jags.data.df$week, unique.weeks)
    
    lambda <- jm.out$jm$sims.list$lambda[, unique.weeks] # posterior samples of lambda
    p.obs <- jm.out$jm$sims.list$p.obs  # posterior samples of p.obs
    obs_prop0 <- mean(jags.data.df$count.obs == 0)  # find count == 0
    n.draws <- nrow(lambda)    # # of posterior samples
    prop0_rep <- numeric(n.draws)   # make sure it is numeric
    n.int <- nrow(jags.data.df)     # find # of non-zero effort shifts 
    y_rep <- matrix(NA_integer_, nrow = n.draws, ncol = n.int)
    
    # for each draw (i.e., posterior sample)
    for (s in 1:n.draws){
      
      lam <- lambda[s, week_col] # Find appropriate lambda for all weeks
      pcorr <- pmin(pmax(p.obs[s] * jags.data.df$effort/3, 0), 1) # Contain pcorr to be [0,1]
      true_rep <- rpois(n.int, lam) # First find the True^{rep}
      y_rep[s,] <- rbinom(n.int, size = true_rep, prob = pcorr) # simulate y_rep
      
    }
    prop0_rep <- rowMeans(y_rep == 0)  # proportion of zeros
    bayes_p <- mean(prop0_rep >= obs_prop0)  # Bayesian p value
    
    ppc.out[[as.character(year)]] <- list(y = jags.data.df$count.obs,
                                          yrep = y_rep,
                                          obs_prop0 = obs_prop0,
                                          bayes_p = bayes_p,
                                          prop0_rep = prop0_rep,
                                          year = year)
    
  }
}

# align on a common number of draws (in case years differ), then cbind by interval
# This (different # reps) shouldn't happen if all datasets were ran at the same time
S <- min(sapply(ppc.out, function(z) nrow(z$yrep)))
y_all    <- unlist(lapply(ppc.out, function(z) z$y), use.names = FALSE)
yrep_all <- do.call(cbind, 
                    lapply(ppc.out, function(z) z$yrep[1:S, , drop = FALSE]))

ppc_bars(y = y_all, 
         yrep = yrep_all[sample(S, min(500, S)), ]) +
  ggplot2::coord_cartesian(xlim = c(-0.5, 6.5))   # counts run 0..~7; trim empty tail

# pooled zero-proportion p-value, comparable to the NB one (compute NB on same years)
obs_prop0_all  <- mean(y_all == 0)
prop0_rep_all  <- rowMeans(yrep_all == 0)
bayes_p_pooled <- mean(prop0_rep_all >= obs_prop0_all)

# Poisson pooled zero-proportions. Rename them to be compared with NB output below
poi_prop0 <- prop0_rep_all               # rowMeans(yrep_all == 0)
poi_obs   <- obs_prop0_all               # mean(y_all == 0)

# NB fixed-year: compute the SAME statistic on the SAME 1994-2025 intervals.
# nb_yrep = draws x intervals matrix of NB replicated counts, restricted to year < 2026
# nb_prop0 <- rowMeans(nb_yrep == 0)
# nb_obs   <- mean(nb_y == 0)              # should equal poi_obs if same intervals
# Bring in the .info file, which contains various objects that were used for 
# running the Stan code
stan.fixed.year.run.date <- "2026-09-11"
fixed.year.stan.info <- readRDS(paste0("RData//GWCalfCount_Stan_nb_fixedyear_",
                                       stan.fixed.year.run.date, ".info"))

# Find out how many entries are for 2026
fixed.year.stan.info$pre.stan.data %>%
  filter(calendar_year == 2026,
         effort_hours > 0) %>%
  mutate(count_per_hr = count_obs/effort_hours) %>%
  group_by(calendar_year) %>%
  mutate(cumu_effort = cumsum(effort_hours)) %>%
  ungroup()  -> data.2026

# Bring in the PPC output for the same model (fixedyear)
fixedyear.ppc <- readRDS("RData//PPC_out_nb_fixedyear.rds")

# Remove the 2026 part from the data and PPC output:
y.rep.1994to2025 <- fixedyear.ppc$y.rep[,1:(dim(fixedyear.ppc$y.rep)[2] - dim(data.2026)[1])]
y.obs.1994to2025 <- fixedyear.ppc$y.obs[1:(dim(fixedyear.ppc$y.rep)[2] - dim(data.2026)[1])]

# Calculate proportion of zeros for simulated and observed data:
nb_prop0_1994to2025 <- rowMeans(y.rep.1994to2025 == 0)
nb_obs_1994to2025 <- mean(y.obs.1994to2025 == 0)

# poi_obs and nb_obs_1994to2025 should be identical

## FOR STAN, OFFSHORE AND INSHORE DATA WERE USED WHEREAS
## FOR JAGS, ONLY INSHORE WAS USED. THIS IS FIXED AND JAGS WERE RUN AGAIN FOR ALL YEARS 2026-09-18 

# Combine poisson and negative binomial results for plotting
reps <- rbind(
  data.frame(model = "Poisson-Binomial", prop0 = poi_prop0),
  data.frame(model = "Negative Binomial", prop0 = nb_prop0_1994to2025)
)

# observed data for plotting and Bayes p values
obs <- data.frame(
  model = c("Poisson-Binomial", "Negative Binomial"),
  obs   = c(poi_obs, nb_obs_1994to2025),
  bayes_p = c(mean(poi_prop0 >= poi_obs), mean(nb_prop0_1994to2025 >= nb_obs_1994to2025))
)


p.PPC <- ggplot(reps, aes(prop0)) +
  geom_histogram(bins = 40) +
  geom_vline(data = obs, aes(xintercept = obs), linewidth = 1) +
  #facet_wrap(~model, scales = "free_x")
  facet_grid(rows = vars(model)) +
  labs(x = "Proportion of zeros", y = "Count") +
  theme_bw()

PPC.compared <- list(y.reps = reps,
                     y.obs = obs,
                     plot = p.PPC)

saveRDS(PPC.compared,
        file = "RData//PPC_compared.rds")

# ggsave(filename = "figures//PPC_compared.png",
#        plot = p.PPC,
#        dpi = 600, device = "png", bg = "transparent")
