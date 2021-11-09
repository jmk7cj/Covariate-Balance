# ----------------------------------------------------------------------------------------------------------------- #  
# Author: Joseph Kush (jkush1@jhu.edu) 
# 
# Title: Covariate Balance for Observational Scale-ups: A Comparison of Matching and Weighting
#
# Date: 10/25/2021 (under construction)
#
# Purpose: Master .R file to conduct Monte Carlo simulations comparing propensity score 
#          matching and weighting methods for achieving covariate balance in observational 
#          studies. 
#
# NOTE: This file is actively under construction
# ----------------------------------------------------------------------------------------------------------------- #  


# ----------------------------------------------------------------------------------------------------------------- #  
# Load packages, set working directory
# ----------------------------------------------------------------------------------------------------------------- #  
rm(list = ls())
gc()
closeAllConnections()
library("MASS")
library("MatchIt")
library("WeightIt")
library("cobalt")
library("survey")
library("data.table")
library("ggplot2")
library("doParallel")
library("doSNOW")

setwd("/Users/my folder")
# ----------------------------------------------------------------------------------------------------------------- #  


# ----------------------------------------------------------------------------------------------------------------- #  
# Simulation factors to vary
# ----------------------------------------------------------------------------------------------------------------- #  
iterations <- 10 # for example 
sample_size <- c(500) 
treatment_effect <- c(0.5) 
treatment_prevalence <- c(0.2, 0.4, 0.6, 0.8) 
n_covs <- c(10)
smd_covs <- c(0.1, 0.2, 0.3, 0.5)
smd_y <- c(0.2)
r_square <- c(1)

# Design matrix
design_matrix <- data.frame(expand.grid(iteration = 1:iterations,
                                        sample_size = sample_size,
                                        treatment_effect = treatment_effect,
                                        treatment_prevalence = treatment_prevalence,
                                        n_covs = n_covs,
                                        smd_covs = smd_covs,
                                        smd_y = smd_y,
                                        r_square = r_square))
design_matrix$i <- 1:nrow(design_matrix)


# prepare machine for parallel processing
processors <- makeCluster(detectCores()[1]-1) 
registerDoSNOW(processors)
pb <- txtProgressBar(max=nrow(design_matrix), style=3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
# ----------------------------------------------------------------------------------------------------------------- #  



# ----------------------------------------------------------------------------------------------------------------- #  
# Function to generate data
# ----------------------------------------------------------------------------------------------------------------- #  
generate_data <- function(sample_size, treatment_prevalence, n_covs, smd_covs, smd_y, treatment_effect, r_square) {
  
  # generate binary treatment variable according to desired treatment prevalence
  treat <- c(rep(1, times = (sample_size * treatment_prevalence)),
             rep(0, times = (sample_size * (1-treatment_prevalence))))
  
  # generate covariates, values determined by treatment status
  for(i in 1:n_covs) {
    assign(paste("x",i, sep=""), ifelse(treat == 1,
                                        rnorm(sample_size, 0, 1),
                                        rnorm(sample_size, smd_covs, 1)),
           envir = .GlobalEnv)
  }
  
  d <- data.frame(cbind(do.call(cbind, lapply(paste0("x", 1:n_covs), get), envir = .GlobalEnv),
                        treat))
  
  # generate potential outcomes
  y <- NULL
  for(i in 1:n_covs) {
    x <- mapply(`*`, d[,i], smd_y)
    y <- data.frame(cbind(y, x))
  }
  
  y0 <- rowSums(y) + treatment_effect*0 + rnorm(nrow(y), 0, r_square)
  y1 <- rowSums(y) + treatment_effect*1 + rnorm(nrow(y), 0, r_square)
  
  d$y0 <- y0
  d$y1 <- y1
  
  d$y_observed <- ifelse(d$treat == 1, d$y1, d$y0)
  
  colnames(d)[1:n_covs] <- paste0("x", seq_len(n_covs))
  
  return(d)
}
# ----------------------------------------------------------------------------------------------------------------- #  


# ----------------------------------------------------------------------------------------------------------------- #  
# Function to estimate propensity scores
# ----------------------------------------------------------------------------------------------------------------- #  
estimate_ps <- function(d) {
  
  prop_treat <- mean(d$treat)
  pop_ATE <- mean(d$y1 - d$y0)
  pop_ATT <- by(d$y1, d$treat, mean)[[2]] - by(d$y0, d$treat, mean)[[2]]
  n_covs <- which(colnames(d) == "treat") - 1
  ps_form <- reformulate(paste0("x", seq_len(n_covs)), "treat")
  
  avg_unbal_smd <- mean(abs(bal.tab(list(treat = d$treat, covs = d[,1:n_covs], s.d.denom = "pooled"))$Balance[,2]))  
  
  match_wr1 <- matchit(ps_form, data=d, link="linear.logit", method="nearest", 
                       estimand="ATT", replace=T, ratio=1)
  
  match_wr3 <- matchit(ps_form, data=d, link="linear.logit", method="nearest", 
                       estimand="ATT", replace=T, ratio=3)
  
  match_wr5 <- matchit(ps_form, data=d, link="linear.logit", method="nearest", 
                       estimand="ATT", replace=T, ratio=5)
  
  weight <- weightit(ps_form, data=d, method = "ps", estimand = "ATT")
  
  # Assess covariate balance 
  avg_mwr1_bal_smd <- mean(abs(bal.tab(match_wr1)$Balance[,3]))
  avg_mwr3_bal_smd <- mean(abs(bal.tab(match_wr3)$Balance[,3]))
  avg_mwr5_bal_smd <- mean(abs(bal.tab(match_wr5)$Balance[,3]))
  avg_wbal_smd <- mean(abs(bal.tab(weight)$Balance[,3]))
  
  # Retain matched units only  
  matchwr1_data <- match.data(match_wr1)
  matchwr3_data <- match.data(match_wr3)
  matchwr5_data <- match.data(match_wr5)
  # none for weighting
  
  # Use weights as matching with replacement
  est_att_mwr1 <- svydesign(ids=~1, weights =~ weights, data = matchwr1_data)
  est_att_mwr3 <- svydesign(ids=~1, weights =~ weights, data = matchwr3_data)
  est_att_mwr5 <- svydesign(ids=~1, weights =~ weights, data = matchwr5_data)
  est_att_weight <- svydesign(~1, weights = weight$weights, data = d)
  
  # Estimate ATT
  est_ATT_mwr1 <- svyglm(reformulate(c("treat", paste0("x", seq_len(n_covs))), response = "y_observed"), est_att_mwr1, family=gaussian())
  est_ATT_mwr3 <- svyglm(reformulate(c("treat", paste0("x", seq_len(n_covs))), response = "y_observed"), est_att_mwr3, family=gaussian())
  est_ATT_mwr5 <- svyglm(reformulate(c("treat", paste0("x", seq_len(n_covs))), response = "y_observed"), est_att_mwr5, family=gaussian())
  est_ATT_weight <- svyglm(reformulate(c("treat", paste0("x", seq_len(n_covs))), response = "y_observed"), est_att_weight, family=gaussian())
  
  att_est_mwr1 <- est_ATT_mwr1$coefficients["treat"]
  att_est_mwr3 <- est_ATT_mwr3$coefficients["treat"]
  att_est_mwr5 <- est_ATT_mwr5$coefficients["treat"]
  att_est_weight <- est_ATT_weight$coefficients["treat"]
  
  return(list(prop_treat, pop_ATE, pop_ATT, avg_unbal_smd, 
              avg_mwr1_bal_smd, avg_mwr3_bal_smd, 
              avg_mwr5_bal_smd, avg_wbal_smd,
              att_est_mwr1, att_est_mwr3, 
              att_est_mwr5, att_est_weight))
}
# ----------------------------------------------------------------------------------------------------------------- #  



# ----------------------------------------------------------------------------------------------------------------- #  
# Conduct simulations using dopar
# ----------------------------------------------------------------------------------------------------------------- #  
results.list <- foreach(i = 1:nrow(design_matrix), .options.snow=opts) %dopar% {
  
  library("MatchIt")
  library("WeightIt")
  library("cobalt")
  library("survey")
  
  d <- generate_data(sample_size = as.numeric(design_matrix[i,2]),
                     treatment_effect = as.numeric(design_matrix[i,3]),
                     treatment_prevalence = as.numeric(design_matrix[i,4]),
                     n_covs = as.numeric(design_matrix[i,5]),
                     smd_covs = as.numeric(design_matrix[i,6]),
                     smd_y = as.numeric(design_matrix[i,7]),
                     r_square = as.numeric(design_matrix[i,8]))
  
  estimate_ps(d)
  
}

results.df <- data.frame(matrix(unlist(results.list), nrow = nrow(design_matrix), byrow = T),
                      stringsAsFactors=FALSE)
names(results.df) <- c("prop_treat", "pop_ATE", "pop_ATT", "avg_unbal_smd", 
                       "avg_mwr1_bal_smd", "avg_mwr3_bal_smd", 
                       "avg_mwr5_bal_smd", "avg_wbal_smd",
                       "att_est_mwr1", "att_est_mwr3", 
                       "att_est_mwr5", "att_est_weight")
head(results.df)
write.csv(results.df, "results_list_10_25_2021.csv", row.names = F)
# ----------------------------------------------------------------------------------------------------------------- #  



# ----------------------------------------------------------------------------------------------------------------- #  
# Bind results list together into a data.frame 
# ----------------------------------------------------------------------------------------------------------------- #  
results <- data.table::rbindlist(results.list, idcol = "i")

# Merge with design matrix
results <- merge(design_matrix, results, by = "i")

colnames(results)[(ncol(design_matrix)+1):ncol(results)] <- c("prop_treat", "pop_ATE", "pop_ATT","avg_unbal_smd", 
                                                              "avg_mwr1_bal_smd", "avg_mwr3_bal_smd", 
                                                              "avg_mwr5_bal_smd", "avg_wbal_smd",
                                                              "att_est_mwr1", "att_est_mwr3", 
                                                              "att_est_mwr5", "att_est_weight") #, 

results$mwr1_bal_smd_reduction <- (results$avg_unbal_smd - results$avg_mwr1_bal_smd) / results$avg_unbal_smd
results$mwr3_bal_smd_reduction <- (results$avg_unbal_smd - results$avg_mwr3_bal_smd) / results$avg_unbal_smd
results$mwr5_bal_smd_reduction <- (results$avg_unbal_smd - results$avg_mwr5_bal_smd) / results$avg_unbal_smd
results$wbal_smd_reduction <- (results$avg_unbal_smd - results$avg_wbal_smd) / results$avg_unbal_smd

results$bias_att_mwr1 <- (results$att_est_mwr1 - results$pop_ATT) / results$pop_ATT
results$bias_att_mwr3 <- (results$att_est_mwr3 - results$pop_ATT) / results$pop_ATT
results$bias_att_mwr5 <- (results$att_est_mwr5 - results$pop_ATT) / results$pop_ATT
results$bias_att_weight <- (results$att_est_weight - results$pop_ATT) / results$pop_ATT

ag_results <- aggregate(. ~ sample_size + treatment_effect + treatment_prevalence + n_covs + smd_covs + smd_y + r_square,
                        data = results[,-1], FUN = mean)

ag_results <- ag_results[rep(seq_len(nrow(ag_results)), each=4), ]
ag_results$method <- rep(c("match_wr1", "match_wr3", "match_wr5", "weight"), each=1)

ag_results$avg_bal <- ifelse(ag_results$method == "match_wr1", ag_results$avg_mwr1_bal_smd, 
                             ifelse(ag_results$method == "match_wr3", ag_results$avg_mwr3_bal_smd, 
                                    ifelse(ag_results$method == "match_wr5", ag_results$avg_mwr5_bal_smd, 
                                           ifelse(ag_results$method == "weight", ag_results$avg_wbal_smd, NA))))

ag_results$smd_reduction <- ifelse(ag_results$method == "match_wr1", ag_results$mwr1_bal_smd_reduction,
                                   ifelse(ag_results$method == "match_wr3", ag_results$mwr3_bal_smd_reduction,
                                          ifelse(ag_results$method == "match_wr5", ag_results$mwr5_bal_smd_reduction,
                                                 ifelse(ag_results$method == "weight", ag_results$wbal_smd_reduction, NA))))

ag_results$att_est <- ifelse(ag_results$method == "match_wr1", ag_results$att_est_mwr1,
                             ifelse(ag_results$method == "match_wr3", ag_results$att_est_mwr3,
                                    ifelse(ag_results$method == "match_wr5", ag_results$att_est_mwr5,
                                           ifelse(ag_results$method == "weight", ag_results$att_est_weight, NA))))

ag_results$bias_att <- ifelse(ag_results$method == "match_wr1", ag_results$bias_att_mwr1,
                              ifelse(ag_results$method == "match_wr3", ag_results$bias_att_mwr3,
                                     ifelse(ag_results$method == "match_wr5", ag_results$bias_att_mwr5,
                                            ifelse(ag_results$method == "weight", ag_results$bias_att_weight, NA))))

ag_results <- subset(ag_results, select = -c(avg_mwr1_bal_smd, 
                                             avg_mwr3_bal_smd, 
                                             avg_mwr5_bal_smd, 
                                             avg_wbal_smd,
                                             mwr1_bal_smd_reduction, 
                                             mwr3_bal_smd_reduction, 
                                             mwr5_bal_smd_reduction, 
                                             wbal_smd_reduction, 
                                             att_est_mwr1,
                                             att_est_mwr3,
                                             att_est_mwr5,
                                             att_est_weight,
                                             bias_att_mwr1, 
                                             bias_att_mwr3, 
                                             bias_att_mwr5, 
                                             bias_att_weight))

head(ag_results)
# ----------------------------------------------------------------------------------------------------------------- #  




# ----------------------------------------------------------------------------------------------------------------- #  
# Plot balance 
# ----------------------------------------------------------------------------------------------------------------- #  
label_method <- c("1:1 Matching with Replacement",
                  "3:1 Matching with Replacement",
                  "5:1 Matching with Replacement",
                  "IPTW")

breaks_method <- c("match_wr1",
                   "match_wr3",
                   "match_wr5",
                   "weight")

label_smd_covs <- c("0.1" = "Covariate Imbalance SMD = 0.1",
                    "0.2" = "Covariate Imbalance SMD = 0.2",
                    "0.3" = "Covariate Imbalance SMD = 0.3",
                    "0.5" = "Covariate Imbalance SMD = 0.5")


plot <- ggplot(data = ag_results, aes(y = (smd_reduction)*100, 
                                      x = as.factor(treatment_prevalence), 
                                      group = as.factor(method))) +
  geom_line(aes(linetype=as.factor(method), color=as.factor(method))) +
  facet_wrap(~ smd_covs, labeller = labeller(.multi_line = T, smd_covs = label_smd_covs)) +
  geom_point(aes(shape=as.factor(method), color=as.factor(method)), size=2) +
  theme_bw() +
  labs(x = "Treatment Prevalence") +
  labs(y = "Percentage Reduction in Covariate Imbalance") + 
  scale_colour_discrete(name = "Propensity Score Method",
                        breaks = breaks_method, labels = label_method) +
  scale_shape_discrete(name = "Propensity Score Method",
                       breaks = breaks_method, labels = label_method) + 
  scale_linetype_discrete(name = "Propensity Score Method",
                          breaks = breaks_method, labels = label_method)
plot
# ----------------------------------------------------------------------------------------------------------------- #  
