# ----------------------------------------------------------------------------------------------------------------- #  
# Author: Joseph Kush (jkush1@jhu.edu) 
# 
# Title: Covariate Balance for Observational Scale-ups: A Comparison of Matching and Weighting
#
# Date: 11/1/2021 (under construction)
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
variance_ratio <- c(0.8)
knn <- c(1, 3, 5)
psmw_method <- c("matching", "weighting")

# Design matrix
design_matrix <- data.frame(expand.grid(iteration = 1:iterations,
                                        sample_size = sample_size,
                                        treatment_effect = treatment_effect,
                                        treatment_prevalence = treatment_prevalence,
                                        n_covs = n_covs,
                                        smd_covs = smd_covs,
                                        smd_y = smd_y,
                                        r_square = r_square, 
                                        variance_ratio = variance_ratio, 
                                        knn = knn, 
                                        psmw_method = psmw_method))

design_matrix <- design_matrix[!(design_matrix$knn > 1 & design_matrix$psmw_method == "weighting"), ]
#design_matrix$knn <- ifelse(design_matrix$psmw_method == "weighting", NA, design_matrix$knn)
design_matrix$knn <- ifelse(design_matrix$psmw_method == "weighting", -99, design_matrix$knn)
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
generate_data <- function(sample_size, treatment_prevalence, n_covs, smd_covs, smd_y, treatment_effect, r_square, variance_ratio) {
  
  # generate binary treatment variable according to desired treatment prevalence
  treat <- c(rep(1, times = (sample_size * treatment_prevalence)),
             rep(0, times = (sample_size * (1-treatment_prevalence))))
  
  # generate covariates, values determined by treatment status
  for(i in 1:n_covs) {
    assign(paste("x",i, sep=""), ifelse(treat == 1,
                                        rnorm(sample_size, 0, sqrt(variance_ratio)),
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
estimate_ps <- function(d, r, psmw_method) {
  
  prop_treat <- mean(d$treat)
  pop_ATE <- mean(d$y1 - d$y0)
  pop_ATT <- by(d$y1, d$treat, mean)[[2]] - by(d$y0, d$treat, mean)[[2]]
  n_covs <- which(colnames(d) == "treat") - 1
  ps_form <- reformulate(paste0("x", seq_len(n_covs)), "treat")
  
  if(psmw_method == "matching") {
    ps_model <- matchit(ps_form, data=d, link="linear.logit", method="nearest", 
                        estimand="ATT", replace=T, ratio=r)
  }
  
  if(psmw_method == "weighting") {
    ps_model <- weightit(ps_form, data=d, method = "ps", estimand = "ATT")
  }
  
  
  # average SMD across all covariates at baseline
  avg_smd_baseline <- mean(abs(bal.tab(ps_model, stats = c("mean.diffs", "variance.ratios"), un=T)$Balance[-c(1),"Diff.Un"]))

  # average variance ratio across all covariates at baseline
  avg_varratio_baseline <- mean(abs(bal.tab(ps_model, stats = c("mean.diffs", "variance.ratios"), un=T)$Balance[-c(1),"V.Ratio.Un"]))

  # average SMD across all covariates after PSM/W
  avg_smd_ps <- mean(abs(bal.tab(ps_model, stats = c("mean.diffs", "variance.ratios"), un=T)$Balance[-c(1),"Diff.Adj"]))

  # average variance ratio across all covariates after PSM/W
  avg_varratio_ps <- mean(abs(bal.tab(ps_model, stats = c("mean.diffs", "variance.ratios"), un=T)$Balance[-c(1),"V.Ratio.Adj"]))

  # average reduction in SMD across all covariates after PSM/W
  avg_smd_reduction <- (mean((avg_smd_ps - avg_smd_baseline) / avg_smd_baseline))*-100
  
  # average reduction in variance ratios across all covariates after PSM/W
  avg_varratio_reduction <- (mean((log(avg_varratio_ps) - log(avg_varratio_baseline)) / log(avg_varratio_baseline)))*-100
  # avg_varratio_reduction <- (mean(abs(((abs(log(bal.tab(ps_model, stats = c("mean.diffs", "variance.ratios"), un=T)$Balance[-c(1),"V.Ratio.Adj"])) - 
  #                                        abs(log(bal.tab(ps_model, stats = c("mean.diffs", "variance.ratios"), un=T)$Balance[-c(1),"V.Ratio.Un"]))) / 
  #                                       abs(log(bal.tab(ps_model, stats = c("mean.diffs", "variance.ratios"), un=T)$Balance[-c(1),"V.Ratio.Un"]))))))*-100
  
  # Retain matched units only  
  if(psmw_method == "matching") {
    match_data <- match.data(ps_model)
  }
  
  # Use weights as matching with replacement
  if(psmw_method == "matching") { 
    est_att <- svydesign(ids=~1, weights =~ weights, data = match_data)
  }
  
  if(psmw_method == "weighting") {
    est_att <- svydesign(~1, weights = ps_model$weights, data = d)
  }
                    
  # Estimate ATT
  est_ATT <- svyglm(reformulate(c("treat", paste0("x", seq_len(n_covs))), response = "y_observed"), est_att, family=gaussian())
  
  att_est <- est_ATT$coefficients["treat"]
  
  return(list(prop_treat, pop_ATE, pop_ATT, 
              avg_smd_baseline, avg_varratio_baseline, 
              avg_smd_ps, avg_varratio_ps, 
              avg_smd_reduction, avg_varratio_reduction, 
              att_est))
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
                     r_square = as.numeric(design_matrix[i,8]), 
                     variance_ratio = as.numeric(design_matrix[i,9]))
  
  estimate_ps(d, r = as.numeric(design_matrix[i,10]), psmw_method = design_matrix[i,11]) 
  
}

results.df <- data.frame(matrix(unlist(results.list), nrow = nrow(design_matrix), byrow = T), 
                             stringsAsFactors=FALSE)
names(results.df) <- c("prop_treat", "pop_ATE", "pop_ATT", 
                       "avg_smd_baseline", "avg_varratio_baseline", 
                        "avg_smd_ps", "avg_varratio_ps", 
                        "avg_smd_reduction", "avg_varratio_reduction", 
                        "att_est") 
head(results.df)
#write.csv(results.df, "results_list_11_1_2021.csv", row.names = F)
# ----------------------------------------------------------------------------------------------------------------- #  






# ----------------------------------------------------------------------------------------------------------------- #  
# Bind results list together into a data.frame 
# ----------------------------------------------------------------------------------------------------------------- #  
results <- data.table::rbindlist(results.list, idcol = "i")

# Merge with design matrix
results <- merge(design_matrix, results, by = "i")

colnames(results)[(ncol(design_matrix)+1):ncol(results)] <- c("prop_treat", "pop_ATE", "pop_ATT", 
                                                              "avg_smd_baseline", "avg_varratio_baseline", 
                                                              "avg_smd_ps", "avg_varratio_ps", 
                                                              "avg_smd_reduction", "avg_varratio_reduction", 
                                                              "att_est") 
# results across all iterations 
head(results)
table(results$psmw_method)

# aggregated results (across unique simulation cells)
ag_results <- aggregate(. ~ sample_size + treatment_effect + treatment_prevalence + n_covs + 
                          smd_covs + smd_y + r_square + variance_ratio + knn + psmw_method,
                        data = results[,-1], FUN = mean)

ag_results$treatment_prevalence <- as.factor(ag_results$treatment_prevalence)
ag_results$psmw_method <- as.factor(ag_results$psmw_method)
ag_results$method <- as.factor(ifelse(ag_results$knn == 1, "psm_1", 
                                      ifelse(ag_results$knn == 3, "psm_3", 
                                             ifelse(ag_results$knn == 5, "psm_5", "weighting"))))
ag_results$bias_att <- ((ag_results$att_est - ag_results$pop_ATT) / ag_results$pop_ATT)*100

head(ag_results)
# ----------------------------------------------------------------------------------------------------------------- #  




# ----------------------------------------------------------------------------------------------------------------- #  
# Plot balance with respect to standardized mean differences 
# ----------------------------------------------------------------------------------------------------------------- #  
label_method <- c("1:1 matching with replacement",
                  "3:1 matching with replacement",
                  "5:1 matching with replacement",
                  "ATT weighting")

breaks_method <- c("psm_1",
                   "psm_3",
                   "psm_5",
                   "weighting")

label_smd_covs <- c("0.1" = "Covariate imbalance SMD = 0.1",
                    "0.2" = "Covariate imbalance SMD = 0.2",
                    "0.3" = "Covariate imbalance SMD = 0.3",
                    "0.5" = "Covariate imbalance SMD = 0.5")


plot <- ggplot(data = ag_results, aes(y = avg_smd_reduction, 
                                      x = treatment_prevalence, 
                                      group = method)) +  
  geom_line(aes(linetype = method, 
                color = method)) +  
  geom_point(aes(shape = method, 
                 color = method), size = 2) + 
  facet_wrap(~ smd_covs, 
             labeller = labeller(.multi_line = T, smd_covs = label_smd_covs)) + 
  labs(x = "Treatment Prevalence") +
  labs(y = "Percentage Reduction in Covariate Imbalance") + 
  scale_colour_discrete(name = "Propensity Score Method",
                        breaks = breaks_method, labels = label_method) +
  scale_shape_discrete(name = "Propensity Score Method",
                       breaks = breaks_method, labels = label_method) + 
  scale_linetype_discrete(name = "Propensity Score Method",
                          breaks = breaks_method, labels = label_method) + 
  theme_bw()

plot
# ----------------------------------------------------------------------------------------------------------------- # 





# ----------------------------------------------------------------------------------------------------------------- #  
# Plot bias in treatment effect estimates
# ----------------------------------------------------------------------------------------------------------------- #  
plot <- ggplot(data = ag_results, aes(y = bias_att, 
                                      x = treatment_prevalence, 
                                      group = method)) +  
  geom_line(aes(linetype = method, 
                color = method)) +  
  geom_point(aes(shape = method, 
                 color = method), size = 2) + 
  facet_wrap(~ smd_covs, 
             labeller = labeller(.multi_line = T, smd_covs = label_smd_covs)) + 
  labs(x = "Treatment Prevalence") +
  labs(y = "Percentage Bias in ATT Estimate") + 
  scale_colour_discrete(name = "Propensity Score Method",
                        breaks = breaks_method, labels = label_method) +
  scale_shape_discrete(name = "Propensity Score Method",
                       breaks = breaks_method, labels = label_method) + 
  scale_linetype_discrete(name = "Propensity Score Method",
                          breaks = breaks_method, labels = label_method) + 
  theme_bw()

plot
# ----------------------------------------------------------------------------------------------------------------- # 
