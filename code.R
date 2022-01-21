# ----------------------------------------------------------------------------------------------------------------- #  
# Author: Joseph Kush (jkush1@jhu.edu) 
# 
# Title: Covariate Balance for Observational Scale-ups: A Comparison of Matching and Weighting
#
# Date: 1/21/2022 
#
# Purpose: Master .R file to conduct Monte Carlo simulations comparing propensity score 
#          matching and weighting methods for achieving covariate balance in observational 
#          studies. 
# ----------------------------------------------------------------------------------------------------------------- #  


# ----------------------------------------------------------------------------------------------------------------- #  
# Load packages, set working directory
# ----------------------------------------------------------------------------------------------------------------- #  
rm(list = ls())
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
library("tidyverse")

setwd("/Users/joekush/Desktop/myfolder")

iterations <- 10 # for example 
sample_size <- c(250, 500)
n_covs <- c(10, 20)
treatment_effect <- c(0.5)
treatment_prevalence <- c(0.2, 0.4, 0.6, 0.8)
smd_covs <- c(0.2, 0.4, 0.6)
r_square <- c(0.3)
knn <- c(1, 3, 5)
psmw_method <- c("matching", "weighting")

# Design matrix
design_matrix <- data.frame(expand.grid(iteration = 1:iterations,
                                        sample_size = sample_size,
                                        treatment_effect = treatment_effect,
                                        treatment_prevalence = treatment_prevalence,
                                        n_covs = n_covs,
                                        smd_covs = smd_covs,
                                        r_square = r_square,
                                        knn = knn,
                                        psmw_method = psmw_method))

design_matrix <- design_matrix[!(design_matrix$knn > 1 & design_matrix$psmw_method == "weighting"), ]
design_matrix$knn <- ifelse(design_matrix$psmw_method == "weighting", -99, design_matrix$knn)
design_matrix$i <- 1:nrow(design_matrix)

processors <- makeCluster(detectCores()[1]-1)
registerDoSNOW(processors)
pb <- txtProgressBar(max=nrow(design_matrix), style=3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
# ----------------------------------------------------------------------------------------------------------------- #  



# ----------------------------------------------------------------------------------------------------------------- #  
# Function to generate data
# ----------------------------------------------------------------------------------------------------------------- #
generate_data <- function(sample_size, treatment_prevalence, n_covs, smd_covs, r_square, treatment_effect) {
  
  for(i in 1:n_covs) {
    assign(paste("x",i, sep=""), rnorm(sample_size, 0, 1), envir = .GlobalEnv)
  }
  
  d <- data.frame(cbind(do.call(cbind, lapply(paste0("x", 1:n_covs), get), envir = .GlobalEnv)))
  
  # Function to output the % of sample receiving treatment (treat=1)
  if(n_covs == 10) {
    determine_intercept <- function(int) {
      logit_treat <- int + smd_covs*d[,1] + smd_covs*d[,2] + smd_covs*d[,3] + smd_covs*d[,4] + smd_covs*d[,5] +
        smd_covs*d[,6] + smd_covs*d[,7] + smd_covs*d[,4]*d[,5] + smd_covs*d[,1]*d[,1] + smd_covs*d[,7]*d[,7]
      
      # convert logit odds into probability
      prob_treat <- exp(logit_treat)/(1 + exp(logit_treat))
      return(mean(prob_treat))
    }
  }
  
  if(n_covs == 20) {
    determine_intercept <- function(int) {
      logit_treat <- int + smd_covs*d[,1] + smd_covs*d[,2] + smd_covs*d[,3] + smd_covs*d[,4] + smd_covs*d[,5] +
        smd_covs*d[,6] + smd_covs*d[,7] + smd_covs*d[,8] + smd_covs*d[,9] + smd_covs*d[,10] +
        smd_covs*d[,11] + smd_covs*d[,12] + smd_covs*d[,13] + smd_covs*d[,14] +
        smd_covs*d[,9]*d[,10] + smd_covs*d[,7]*d[,8] + smd_covs*d[,1]*d[,1] + smd_covs*d[,10]*d[,10] +
        smd_covs*d[,3]*d[,3] + smd_covs*d[,7]*d[,7]
      
      # convert logit odds into probability
      prob_treat <- exp(logit_treat)/(1 + exp(logit_treat))
      return(mean(prob_treat))
    }
  }
  
  # Bisection approach to determine correct intercept value corresponding to a desired treatment prevalence
  bisection <- function(treatment_prevalence, lower=-10, upper=10, maxError=0.0001, maxIter=1000) {
    
    a <- lower
    b <- upper
    iter <- 0
    
    fa <- determine_intercept(int = a)
    fb <- determine_intercept(int = b)
    
    c <- (a+b) / 2
    
    diff_a <- abs(fa - treatment_prevalence)
    diff_b <- abs(fb - treatment_prevalence)
    
    while(iter <= maxIter) {
      c <- (a+b) / 2
      fc <- determine_intercept(int = c)
      
      if(abs(fc - treatment_prevalence) < maxError) {
        break
      }
      
      fa <- determine_intercept(int = a)
      fb <- determine_intercept(int = b)
      
      diff_a <- fa - treatment_prevalence
      diff_c <- fc - treatment_prevalence
      
      if(sign(diff_c) == sign(diff_a)) {
        a <- c
      } else {
        b <- c
      }
      
      iter <- iter+1
    }
    
    if(iter >= maxIter) {
      print(c)
      stop("Search failed. Increase number of iterations")
    }
    
    return(c)
    
  }
  
  
  # Now we can determine the correct intercept value corresponding to a desired treatment exposure level
  if(n_covs == 10) {
    int <- bisection(treatment_prevalence, lower=-10, upper=10, maxError=0.0001, maxIter=1000)
    logit_treat <- int + smd_covs*d[,1] + smd_covs*d[,2] + smd_covs*d[,3] + smd_covs*d[,4] + smd_covs*d[,5] +
      smd_covs*d[,6] + smd_covs*d[,7] + smd_covs*d[,4]*d[,5] + smd_covs*d[,1]*d[,1] + smd_covs*d[,7]*d[,7]
    
    # convert logit odds into probability
    prob_treat <- exp(logit_treat)/(1 + exp(logit_treat))
    d$prob_treat <- prob_treat
    
    # generate binary treatment indicator
    treat <- rbinom(nrow(d), 1, prob_treat)
    d$treat <- treat
    
    
    # create potential outcomes 
      b <- 0.195
      y0 <- rnorm(nrow(d), (0 + treatment_effect*0 + b*d[,1] + b*d[,2] + b*d[,3] + b*d[,4] +
                              b*d[,5] + b*d[,8] + b*d[,9] + b*d[,2]*d[,4] + b*d[,3]*d[,5] + b*d[,1]*d[,1]), 1)
      
      y1 <- rnorm(nrow(d), (0 + treatment_effect*1 + b*d[,1] + b*d[,2] + b*d[,3] + b*d[,4] + 
                              b*d[,5] + b*d[,8] + b*d[,9] + b*d[,2]*d[,4] + b*d[,3]*d[,5] + b*d[,1]*d[,1]), 1)
                              
  }
  
  
  if(n_covs == 20) {
    int <- bisection(treatment_prevalence, lower=-10, upper=10, maxError=0.0001, maxIter=1000)
    logit_treat <- int + smd_covs*d[,1] + smd_covs*d[,2] + smd_covs*d[,3] + smd_covs*d[,4] + smd_covs*d[,5] +
      smd_covs*d[,6] + smd_covs*d[,7] + smd_covs*d[,8] + smd_covs*d[,9] + smd_covs*d[,10] +
      smd_covs*d[,11] + smd_covs*d[,12] + smd_covs*d[,13] + smd_covs*d[,14] +
      smd_covs*d[,9]*d[,10] + smd_covs*d[,7]*d[,8] + smd_covs*d[,1]*d[,1] + smd_covs*d[,10]*d[,10] +
      smd_covs*d[,3]*d[,3] + smd_covs*d[,7]*d[,7]
    
    # convert logit odds into probability
    prob_treat <- exp(logit_treat)/(1 + exp(logit_treat))
    d$prob_treat <- prob_treat
    
    # generate binary treatment indicator
    treat <- rbinom(nrow(d), 1, prob_treat)
    d$treat <- treat
    
    # create potential outcomes 
      b <- 0.185
      y0 <- rnorm(nrow(d), (0 + treatment_effect*0 +
                              b*d[,1] + b*d[,2] + b*d[,3] + b*d[,4] + b*d[,5] +
                              b*d[,6] + b*d[,7] + b*d[,8] + b*d[,9] + b*d[,10] +
                              b*d[,15] + b*d[,16] + b*d[,17] + b*d[,18] +
                              b*d[,3]*d[,5] + b*d[,1]*d[,1]), 1)
      
      y1 <- rnorm(nrow(d), (0 + treatment_effect*1 +
                              b*d[,1] + b*d[,2] + b*d[,3] + b*d[,4] + b*d[,5] +
                              b*d[,6] + b*d[,7] + b*d[,8] + b*d[,9] + b*d[,10] +
                              b*d[,15] + b*d[,16] + b*d[,17] + b*d[,18] +
                              b*d[,3]*d[,5] + b*d[,1]*d[,1]), 1)
   }
  
  d$y0 <- y0
  d$y1 <- y1
  d$logit_treat <- logit_treat
  d$prob_treat <- prob_treat
  
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
  n_covs <- which(colnames(d) == "treat") - 2
  ps_form <- reformulate(paste0("x", seq_len(n_covs)), "treat")
  
  if(psmw_method == "matching") {
    ps_model <- matchit(ps_form, data=d, link="linear.logit", method="nearest",
                        estimand="ATT", replace=T, ratio=r)
  }
  
  if(psmw_method == "weighting") {
    ps_model <- weightit(ps_form, data=d, method = "ps", estimand = "ATT")
  }
  
  # average SMD across all covariates at baseline
  avg_smd_baseline <- mean(abs(bal.tab(ps_model, stats = c("mean.diffs", "variance.ratios"), 
                                       un=T, s.d.denom = "treated")$Balance[-c(1),"Diff.Un"]))
  
  # average variance ratio across all covariates at baseline
  avg_varratio_baseline <- mean(abs(bal.tab(ps_model, stats = c("mean.diffs", "variance.ratios"), 
                                            un=T, s.d.denom = "treated")$Balance[-c(1),"V.Ratio.Un"]))
  
  # average SMD across all covariates after PSM/W
  avg_smd_ps <- mean(abs(bal.tab(ps_model, stats = c("mean.diffs", "variance.ratios"), 
                                 un=T, s.d.denom = "treated")$Balance[-c(1),"Diff.Adj"]))
  
  # average variance ratio across all covariates after PSM/W
  avg_varratio_ps <- mean(abs(bal.tab(ps_model, stats = c("mean.diffs", "variance.ratios"),
                                      un=T, s.d.denom = "treated")$Balance[-c(1),"V.Ratio.Adj"]))
  
  # average reduction in SMD across all covariates after PSM/W
  avg_smd_reduction <- ((avg_smd_ps - avg_smd_baseline) / avg_smd_baseline)*-100
  
  # average reduction in variance ratios across all covariates after PSM/W
  avg_varratio_reduction <- ((avg_varratio_ps - avg_varratio_baseline) / avg_varratio_baseline)*-100
  
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
  
  # Estimate ATT (unadjusted)
  est_ATT_unadj <- svyglm(y_observed ~ treat, est_att, family=gaussian())
  att_est_unadj <- est_ATT_unadj$coefficients["treat"]
  att_se_unadj <- SE(est_ATT_unadj)["treat"]
  
  # Estimate ATT (adjusted via linear combination of predictors)
  est_ATT_linear <- svyglm(reformulate(c("treat", paste0("x", seq_len(n_covs))), "y_observed"), est_att, family=gaussian())
  att_est_linear <- est_ATT_linear$coefficients["treat"]
  att_se_linear <- SE(est_ATT_linear)["treat"]
  
  return(list(prop_treat, pop_ATE, pop_ATT,
              avg_smd_baseline, avg_varratio_baseline,
              avg_smd_ps, avg_varratio_ps,
              avg_smd_reduction, avg_varratio_reduction,
              att_est_unadj, att_se_unadj, 
              att_est_linear, att_se_linear))
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
                     r_square = as.numeric(design_matrix[i,7]))
  
  estimate_ps(d, r = as.numeric(design_matrix[i,8]), psmw_method = design_matrix[i,9])
  
}

results.df <- data.frame(matrix(unlist(results.list), nrow = nrow(design_matrix), byrow = T),
                         stringsAsFactors=FALSE)
names(results.df) <- c("prop_treat", "pop_ATE", "pop_ATT",
                       "avg_smd_baseline", "avg_varratio_baseline",
                       "avg_smd_ps", "avg_varratio_ps",
                       "avg_smd_reduction", "avg_varratio_reduction",
                       "att_est_unadj", "att_se_unadj",
                       "att_est_linear", "att_se_linear")
head(results.df)
write.csv(results.df, "results_df.csv", row.names = F)
# ----------------------------------------------------------------------------------------------------------------- # 
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
                                                               "att_est_unadj", "att_se_unadj", 
                                                              "att_est_linear", "att_se_linear")

# calculate variance of ATT estimates (variance across replications - need later for MSE)
results[seq(from = 1, to = nrow(results), by = iterations), "var_att_est_unadj"] <-aggregate(
  att_est_unadj ~ sample_size + treatment_effect + treatment_prevalence + n_covs + 
    smd_covs + r_square + knn + psmw_method, data = results, FUN = var)[,"att_est_unadj"]

results[seq(from = 1, to = nrow(results), by = iterations), "var_att_est_linear"] <-aggregate(
  att_est_linear ~ sample_size + treatment_effect + treatment_prevalence + n_covs + 
    smd_covs + r_square + knn + psmw_method, data = results, FUN = var)[,"att_est_linear"]

results <- results %>% 
  group_by(sample_size, treatment_effect, treatment_prevalence, n_covs, smd_covs, r_square, knn, psmw_method) %>% 
  fill(var_att_est_unadj) %>% 
  fill(var_att_est_unadj, .direction = "down")

results <- results %>% 
  group_by(sample_size, treatment_effect, treatment_prevalence, n_covs, smd_covs, r_square, knn, psmw_method) %>% 
  fill(var_att_est_linear) %>% 
  fill(var_att_est_linear, .direction = "down")

results <- data.frame(results)

# results across all iterations
head(results)

# aggregated results (across unique simulation cells)
ag_results <- aggregate(. ~ sample_size + treatment_effect + treatment_prevalence + n_covs +
                          smd_covs + r_square + knn + psmw_method,
                        data = results[,-1], FUN = mean)

ag_results$treatment_prevalence <- as.factor(ag_results$treatment_prevalence)
ag_results$psmw_method <- as.factor(ag_results$psmw_method)
ag_results$r_square <- as.factor(ag_results$r_square)
ag_results$method <- as.factor(ifelse(ag_results$knn == 1, "psm_1",
                                      ifelse(ag_results$knn == 3, "psm_3",
                                             ifelse(ag_results$knn == 5, "psm_5", "weighting"))))

ag_results$bias_att_unadj <- ((ag_results$att_est_unadj - ag_results$pop_ATT) / ag_results$pop_ATT)*100
ag_results$bias_att_linear <- ((ag_results$att_est_linear - ag_results$pop_ATT) / ag_results$pop_ATT)*100
ag_results$mse_unadj <- (ag_results$bias_att_unadj * ag_results$bias_att_unadj) + ag_results$var_att_est_unadj
ag_results$mse_linear <- (ag_results$bias_att_linear * ag_results$bias_att_linear) + ag_results$var_att_est_linear

head(ag_results)
# ----------------------------------------------------------------------------------------------------------------- #  





# ----------------------------------------------------------------------------------------------------------------- #  
# Plot results (create labels for graphs)
# ----------------------------------------------------------------------------------------------------------------- #  
label_method <- c("1:1 matching with replacement",
                  "3:1 matching with replacement",
                  "5:1 matching with replacement",
                  "ATT weighting")

breaks_method <- c("psm_1",
                   "psm_3",
                   "psm_5",
                   "weighting")

label_smd_covs <- c("0.2" = "Std. Mean Difference = 0.2",
                    "0.4" = "Std. Mean Difference = 0.4",
                    "0.6" = "Std. Mean Difference = 0.6")

label_n_covs <- c("10" = "Covariates = 10",
                  "20" = "Covariates = 20")

label_sample_size <- c("250" = "Sample Size = 250",
                       "500" = "Sample Size = 500")

label_r_square <- c("0.1" = "R-square = 0.1", 
                    "0.3" = "R-square = 0.3",
                    "0.5" = "R-square = 0.5")

label_treatment_effect <- c("0.2" = "Treatment Effect = 0.2", 
                            "0.5" = "Treatment Effect = 0.5")
# ----------------------------------------------------------------------------------------------------------------- #  





# ----------------------------------------------------------------------------------------------------------------- #  
# Reduction in Covariate Imbalance
# ----------------------------------------------------------------------------------------------------------------- #
p1 <- ggplot(data = subset(ag_results, subset = sample_size == 500 & n_covs == 10), 
       aes(y = avg_smd_ps, x = treatment_prevalence, group = method)) +  
  geom_line(aes(linetype = method, color = method)) +  
  geom_point(aes(shape = method, color = method), size = 2) +
  facet_wrap(~ smd_covs, labeller = labeller(.multi_line = T, smd_covs = label_smd_covs)) +
  labs(x = "Treatment Prevalence") +
  labs(y = "Standardized Mean Difference") +
  scale_colour_discrete(name = "Propensity Score Method",
                        breaks = breaks_method, labels = label_method) +
  scale_shape_discrete(name = "Propensity Score Method",
                       breaks = breaks_method, labels = label_method) +
  scale_linetype_discrete(name = "Propensity Score Method",
                          breaks = breaks_method, labels = label_method) +
  theme_bw() + ggtitle("Covariate Imbalance") + theme(plot.title = element_text(hjust = 0.5)) + 
  theme(text = element_text(size = 14))
p1
# ----------------------------------------------------------------------------------------------------------------- #




# ----------------------------------------------------------------------------------------------------------------- #  
# Bias in unadjusted ATT Estimate
# ----------------------------------------------------------------------------------------------------------------- #
p2 <- ggplot(data = subset(ag_results, subset = sample_size == 500 & n_covs == 10), 
       aes(y = bias_att_unadj, x = treatment_prevalence, group = method)) +  
  geom_line(aes(linetype = method, color = method)) +  
  geom_point(aes(shape = method, color = method), size = 2) +
  facet_wrap(~ smd_covs, labeller = labeller(.multi_line = T, smd_covs = label_smd_covs)) +
  labs(x = "Treatment Prevalence") +
  labs(y = "% Bias") +
  scale_colour_discrete(name = "Propensity Score Method",
                        breaks = breaks_method, labels = label_method) +
  scale_shape_discrete(name = "Propensity Score Method",
                       breaks = breaks_method, labels = label_method) +
  scale_linetype_discrete(name = "Propensity Score Method",
                          breaks = breaks_method, labels = label_method) +
  theme_bw() + ggtitle("Bias of Unadjusted ATT Estimates") + theme(plot.title = element_text(hjust = 0.5)) + 
  theme(text = element_text(size = 14))
p2
# ----------------------------------------------------------------------------------------------------------------- #




# ----------------------------------------------------------------------------------------------------------------- #  
# MSE unadjusted
# ----------------------------------------------------------------------------------------------------------------- #
p3 <- ggplot(data = subset(ag_results, subset = sample_size == 500 & n_covs == 10), 
       aes(y = mse_unadj, x = treatment_prevalence, group = method)) +  
  geom_line(aes(linetype = method, color = method)) +  
  geom_point(aes(shape = method, color = method), size = 2) +
  facet_wrap(~ smd_covs, labeller = labeller(.multi_line = T, smd_covs = label_smd_covs)) +
  labs(x = "Treatment Prevalence") +
  labs(y = "MSE") +
  scale_colour_discrete(name = "Propensity Score Method",
                        breaks = breaks_method, labels = label_method) +
  scale_shape_discrete(name = "Propensity Score Method",
                       breaks = breaks_method, labels = label_method) +
  scale_linetype_discrete(name = "Propensity Score Method",
                          breaks = breaks_method, labels = label_method) +
  theme_bw() + ggtitle("MSE of Unadjusted ATT Estimates") + theme(plot.title = element_text(hjust = 0.5)) + 
  theme(text = element_text(size = 14))
p3
# ----------------------------------------------------------------------------------------------------------------- #



# ----------------------------------------------------------------------------------------------------------------- #
# Combine all plots into single graph
ggarrange(p1,p2,p3, ncol = 1, nrow=3, common.legend = T, legend = "right")
# ----------------------------------------------------------------------------------------------------------------- #
# END
# ----------------------------------------------------------------------------------------------------------------- #
