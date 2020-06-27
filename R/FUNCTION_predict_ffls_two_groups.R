#from prll_load_02_12


#two-group function: predict ffls in unique of one of two biological groups

#' @title predict_ffls_two_groups
#' @description function to predict FFLs unique to one of two biological groups
#'
#' @param mirna_expr_g1 Group 1's miRNA expression data (dataframe: miRNAs x samples, see \code{sample_mirna_expr_g1})
#' @param mrna_expr_g1 Group 1's mRNA expression data (dataframe: mRNAs x samples, see \code{sample_mrna_expr_g1})
#' @param mirna_expr_g2 Group 2's miRNA expression data (dataframe: miRNAs x samples, see \code{sample_mirna_expr_g2})
#' @param mrna_expr_g2 Group 2's miRNA expression data (dataframe: miRNAs x samples, see \code{sample_mirna_expr_g2})
#' @param ffl_type FFL type (character: "miRNA" or "TF")
#' @param candidate_ffls candidate FFLs (dataframe: candidate FFLs x 10, see \code{sample_candidate_ffls_mirna})
#' @param first_row first row of candidate_ffls dataframe to include in analyses (integer: default = 1)
#' @param last_row last row of candidate_ffls dataframe to include in analyses (integer: default = nrow(candidate_ffls))
#' @param num_bootstrap_samples number of bootstrap samples (integer: default = 1000)
#' @param num_permutations number of permutations (integer: default = 1000)
#' @param p_value_adjust_method method to adjust p-values for multiple testing (character: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "locfdr", or "none")
#' @param seed random seed (integer: default = 12345)
#' @param parallel whether to run analyses in parallel (boolean: default = TRUE)
#'
#' @return candidate FFLs with delta-P(FFL) values, p-values, and coefficient estimates (dataframe: number of candidate FFLs x 85 for miRNA-FFLs; number of candidate FFLs x 86 for TF-FFLs; see sample output)
#' @export


predict_ffls_two_groups <- function(mirna_expr_g1, mrna_expr_g1, mirna_expr_g2, mrna_expr_g2,
                                    ffl_type = c("miRNA", "TF"),
                                    candidate_ffls, first_row = 1, last_row = nrow(candidate_ffls),
                                    num_bootstrap_samples = 1000, num_permutations = 1000,
                                    p_value_adjust_method = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "locfdr", "none"),
                                    seed = 12345, parallel = TRUE){
  #check for errors
  if(ncol(mirna_expr_g1) == 0) stop("mirna_expr_g1 is empty (0 rows)")
  if(ncol(mrna_expr_g1) == 0) stop("mrna_expr_g1 is empty (0 rows)")
  if(ncol(mirna_expr_g2) == 0) stop("mirna_expr_g2 is empty (0 rows)")
  if(ncol(mrna_expr_g2) == 0) stop("mrna_expr_g2 is empty (0 rows)")
  if(nrow(candidate_ffls) == 0) stop("candidate_ffls is empty (0 rows)")
  if((first_row > nrow(candidate_ffls)) | (last_row > nrow(candidate_ffls))) stop("first_row and/or last_row larger than number of rows in candidate_ffls")
  if(!ncol(mirna_expr_g1) == ncol(mrna_expr_g1)) stop("mirna_expr_g1 and mrna_expr_g1 do not have the same number of samples (rows)")
  if(!ncol(mirna_expr_g2) == ncol(mrna_expr_g2)) stop("mirna_expr_g2 and mrna_expr_g2 do not have the same number of samples (rows)")

  #functions used in foreach
  #mediation function (mirna-ffls)
  mediation_ffl <- function(mirna, tf, targetgene, ffl_type = c("miRNA", "TF")){
    ###1. set variables
    if(ffl_type == "miRNA"){X <- mirna; M <- tf; Y <- targetgene}
    if(ffl_type == "TF"){X <- tf; M <- mirna; Y <- targetgene}

    ###2. fit models; store coefficients
    #model1: Y ~ X
    model1 <- lm(Y ~ X)
    alpha1 <- summary(model1)$coefficients["X", ]
    #model2: M ~ X
    model2 <- lm(M ~ X)
    beta1 <- summary(model2)$coefficients["X", ]
    #model3: Y ~ X + M
    model3 <- lm(Y ~ X + M)
    gamma1 <- summary(model3)$coefficients["X", ]
    gamma2 <- summary(model3)$coefficients["M", ]

    ###3. determine if candidate ffl meets mediation conditions
    #mirna-ffl conditions
    if(ffl_type == "miRNA"){
      alpha1_crit <- alpha1["Estimate"] < 0
      beta1_crit <- beta1["Estimate"] < 0
      gamma1_crit <- (gamma1["Estimate"] < 0) & (abs(gamma1["Estimate"]) < abs(alpha1["Estimate"]))
      gamma2_crit <- gamma2["Estimate"] > 0
    }
    #tf-ffl conditions
    if(ffl_type == "TF"){
      alpha1_crit <- alpha1["Estimate"] > 0
      beta1_crit <- beta1["Estimate"] > 0
      gamma1_crit <- (gamma1["Estimate"] > 0) & (abs(gamma1["Estimate"]) > abs(alpha1["Estimate"]))
      gamma2_crit <- gamma2["Estimate"] < 0
    }

    #determine if candidate ffl meets mediation conditions
    meet_conditions <- alpha1_crit & beta1_crit & gamma1_crit & gamma2_crit

    return(meet_conditions)
  }
  #bootstrap function
  bootstrap_stat <- function(data, indices, ffl_type){
    #allows boot to select sample
    d <- data[indices, ]
    #determine if each bootstrap sample meets mediation analysis conditions
    meet_conditions <- mediation_ffl(mirna=d$mirna, tf=d$tf, targetgene=d$targetgene, ffl_type=ffl_type)
    return(meet_conditions)
  }

  #start parallelization
  if(parallel){num_cores <- detectCores()-1}
  else{num_cores <- 1}
  cl <- makeCluster(num_cores, outfile = "")
  registerDoParallel(cl)

  #iterate through each candidate ffl
  candidate_ffls <- candidate_ffls[first_row:last_row, ] #subset candidate ffls
  predicted_ffls <- foreach(i = 1:nrow(candidate_ffls), .combine = rbind) %dopar% {
    ######step1. identify candidate ffl: extract mirna, tf, and targetgene expression levels
    row <- candidate_ffls[i, ] #row is the candidate ffl
    #group1
    mirna_g1 <- t(mirna_expr_g1[row[["mirna"]], ])
    tf_g1 <- t(mrna_expr_g1[row[["tf"]], ])
    targetgene_g1 <- t(mrna_expr_g1[row[["targetgene"]], ])
    #group2
    mirna_g2 <- t(mirna_expr_g2[row[["mirna"]], ])
    tf_g2 <- t(mrna_expr_g2[row[["tf"]], ])
    targetgene_g2 <- t(mrna_expr_g2[row[["targetgene"]], ])

    ######step2. calculate regression models, store coefficients
    ###2a. set variables
    if(ffl_type == "miRNA"){
      X_g1 <- mirna_g1
      M_g1 <- tf_g1
      Y_g1 <- targetgene_g1
      X_g2 <- mirna_g2
      M_g2 <- tf_g2
      Y_g2 <- targetgene_g2
    }
    if(ffl_type == "TF"){
      X_g1 <- tf_g1
      M_g1 <- mirna_g1
      Y_g1 <- targetgene_g1
      X_g2 <- tf_g2
      M_g2 <- mirna_g2
      Y_g2 <- targetgene_g2
    }

    ###2b. fit models --- group1
    #model1: Y ~ X (targetgene ~ mirna)
    model1 <- lm(Y_g1 ~ X_g1)
    alpha0 <- summary(model1)$coefficients["(Intercept)", ]
    alpha1 <- summary(model1)$coefficients["X_g1", ]
    #model2: M ~ X (tf ~ mirna)
    model2 <- lm(M_g1 ~ X_g1)
    beta0 <- summary(model2)$coefficients["(Intercept)", ]
    beta1 <- summary(model2)$coefficients["X_g1", ]
    #model3: Y ~ X + M (targetgene ~ mirna + tf)
    model3 <- lm(Y_g1 ~ X_g1 + M_g1)
    gamma0 <- summary(model3)$coefficients["(Intercept)", ]
    gamma1 <- summary(model3)$coefficients["X_g1", ]
    gamma2 <- summary(model3)$coefficients["M_g1", ]
    ###2c. store information --- group1
    #coefficients
    #alpha0
    row[["alpha0_g1_estimate"]] <- alpha0["Estimate"]
    row[["alpha0_g1_se"]] <- alpha0["Std. Error"]
    row[["alpha0_g1_t"]] <- alpha0["t value"]
    row[["alpha0_g1_p_value"]] <- alpha0["Pr(>|t|)"]
    #alpha1
    row[["alpha1_g1_estimate"]] <- alpha1["Estimate"]
    row[["alpha1_g1_se"]] <- alpha1["Std. Error"]
    row[["alpha1_g1_t"]] <- alpha1["t value"]
    row[["alpha1_g1_p_value"]] <- alpha1["Pr(>|t|)"]
    #beta0
    row[["beta0_g1_estimate"]] <- beta0["Estimate"]
    row[["beta0_g1_se"]] <- beta0["Std. Error"]
    row[["beta0_g1_t"]] <- beta0["t value"]
    row[["beta0_g1_p_value"]] <- beta0["Pr(>|t|)"]
    #beta1
    row[["beta1_g1_estimate"]] <- beta1["Estimate"]
    row[["beta1_g1_se"]] <- beta1["Std. Error"]
    row[["beta1_g1_t"]] <- beta1["t value"]
    row[["beta1_g1_p_value"]] <- beta1["Pr(>|t|)"]
    #gamma0
    row[["gamma0_g1_estimate"]] <- gamma0["Estimate"]
    row[["gamma0_g1_se"]] <- gamma0["Std. Error"]
    row[["gamma0_g1_t"]] <- gamma0["t value"]
    row[["gamma0_g1_p_value"]] <- gamma0["Pr(>|t|)"]
    #gamma1
    row[["gamma1_g1_estimate"]] <- gamma1["Estimate"]
    row[["gamma1_g1_se"]] <- gamma1["Std. Error"]
    row[["gamma1_g1_t"]] <- gamma1["t value"]
    row[["gamma1_g1_p_value"]] <- gamma1["Pr(>|t|)"]
    #gamma2
    row[["gamma2_g1_estimate"]] <- gamma2["Estimate"]
    row[["gamma2_g1_se"]] <- gamma2["Std. Error"]
    row[["gamma2_g1_t"]] <- gamma2["t value"]
    row[["gamma2_g1_p_value"]] <- gamma2["Pr(>|t|)"]
    #means
    row[["meanX_g1"]] <- mean(X_g1)
    row[["meanM_g1"]] <- mean(M_g1)
    #variances
    row[["varX_g1"]] <- var(X_g1)*(nrow(X_g1)-1)/nrow(X_g1)
    row[["varM_g1"]] <- var(M_g1)*(nrow(M_g1)-1)/nrow(M_g1)
    row[["covXM_g1"]] <- cov(X_g1, M_g1)*(nrow(X_g1)-1)/nrow(X_g1)
    row[["var_epsilon_g1"]] <- summary(model3)$sigma^2

    ###2b. fit models --- group2
    #model1: Y ~ X (targetgene ~ mirna)
    model1 <- lm(Y_g2 ~ X_g2)
    alpha0 <- summary(model1)$coefficients["(Intercept)", ]
    alpha1 <- summary(model1)$coefficients["X_g2", ]
    #model2: M ~ X (tf ~ mirna)
    model2 <- lm(M_g2 ~ X_g2)
    beta0 <- summary(model2)$coefficients["(Intercept)", ]
    beta1 <- summary(model2)$coefficients["X_g2", ]
    #model3: Y ~ X + M (targetgene ~ mirna + tf)
    model3 <- lm(Y_g2 ~ X_g2 + M_g2)
    gamma0 <- summary(model3)$coefficients["(Intercept)", ]
    gamma1 <- summary(model3)$coefficients["X_g2", ]
    gamma2 <- summary(model3)$coefficients["M_g2", ]
    ###2c. store information --- group2
    #coefficients
    #alpha0
    row[["alpha0_g2_estimate"]] <- alpha0["Estimate"]
    row[["alpha0_g2_se"]] <- alpha0["Std. Error"]
    row[["alpha0_g2_t"]] <- alpha0["t value"]
    row[["alpha0_g2_p_value"]] <- alpha0["Pr(>|t|)"]
    #alpha1
    row[["alpha1_g2_estimate"]] <- alpha1["Estimate"]
    row[["alpha1_g2_se"]] <- alpha1["Std. Error"]
    row[["alpha1_g2_t"]] <- alpha1["t value"]
    row[["alpha1_g2_p_value"]] <- alpha1["Pr(>|t|)"]
    #beta0
    row[["beta0_g2_estimate"]] <- beta0["Estimate"]
    row[["beta0_g2_se"]] <- beta0["Std. Error"]
    row[["beta0_g2_t"]] <- beta0["t value"]
    row[["beta0_g2_p_value"]] <- beta0["Pr(>|t|)"]
    #beta1
    row[["beta1_g2_estimate"]] <- beta1["Estimate"]
    row[["beta1_g2_se"]] <- beta1["Std. Error"]
    row[["beta1_g2_t"]] <- beta1["t value"]
    row[["beta1_g2_p_value"]] <- beta1["Pr(>|t|)"]
    #gamma0
    row[["gamma0_g2_estimate"]] <- gamma0["Estimate"]
    row[["gamma0_g2_se"]] <- gamma0["Std. Error"]
    row[["gamma0_g2_t"]] <- gamma0["t value"]
    row[["gamma0_g2_p_value"]] <- gamma0["Pr(>|t|)"]
    #gamma1
    row[["gamma1_g2_estimate"]] <- gamma1["Estimate"]
    row[["gamma1_g2_se"]] <- gamma1["Std. Error"]
    row[["gamma1_g2_t"]] <- gamma1["t value"]
    row[["gamma1_g2_p_value"]] <- gamma1["Pr(>|t|)"]
    #gamma2
    row[["gamma2_g2_estimate"]] <- gamma2["Estimate"]
    row[["gamma2_g2_se"]] <- gamma2["Std. Error"]
    row[["gamma2_g2_t"]] <- gamma2["t value"]
    row[["gamma2_g2_p_value"]] <- gamma2["Pr(>|t|)"]
    #means
    row[["meanX_g2"]] <- mean(X_g2)
    row[["meanM_g2"]] <- mean(M_g2)
    #variances
    row[["varX_g2"]] <- var(X_g2)*(nrow(X_g2)-1)/nrow(X_g2)
    row[["varM_g2"]] <- var(M_g2)*(nrow(M_g2)-1)/nrow(M_g2)
    row[["covXM_g2"]] <- cov(X_g2, M_g2)*(nrow(X_g2)-1)/nrow(X_g2)
    row[["var_epsilon_g2"]] <- summary(model3)$sigma^2

    #####step3. calculate delta-p(FFL)
    #create expression df (expr_data)
    expr_data_g1 <- data.frame(mirna_g1, tf_g1, targetgene_g1) #group1
    colnames(expr_data_g1) <- c("mirna", "tf", "targetgene")
    expr_data_g2 <- data.frame(mirna_g2, tf_g2, targetgene_g2) #group2
    colnames(expr_data_g2) <- c("mirna", "tf", "targetgene")
    #bootstrap to calculate delta-p(FFL)
    set.seed(seed)
    boostrap_results_g1 <- boot::boot(data = expr_data_g1, statistic = bootstrap_stat, R = num_bootstrap_samples, ffl_type = ffl_type, parallel = "multicore")
    boostrap_results_g2 <- boot::boot(data = expr_data_g2, statistic = bootstrap_stat, R = num_bootstrap_samples, ffl_type = ffl_type, parallel = "multicore")
    obs_delta_pffl <- mean(boostrap_results_g1[["t"]]) - mean(boostrap_results_g2[["t"]]) #calculate delta-p(FFL)
    row[["delta_pFFL_(group1-group2)"]] <- obs_delta_pffl
    row[["pFFL_group1"]] <- mean(boostrap_results_g1[["t"]])
    row[["pFFL_group2"]] <- mean(boostrap_results_g2[["t"]])

    #####step4. calculate p-value
    #create empty vector to store results
    null_delta_pffls <- rep(NA, num_permutations)
    #conduct permutations
    set.seed(seed)
    for(perm in 1:num_permutations){
      ###a. permute data --- shuffle group assignment
      expr_data_combined <- rbind(expr_data_g1, expr_data_g2)
      row_nums_perm_g1 <- sample(nrow(expr_data_combined), nrow(expr_data_g1)) #sample rows to be in group 1
      all_rows <- seq(1, nrow(expr_data_combined))
      row_nums_perm_g2 <- all_rows[!(all_rows %in% row_nums_perm_g1)] #remaining rows that are not sampled are in group 2
      expr_data_perm_g1 <- expr_data_combined[row_nums_perm_g1, ]
      expr_data_perm_g2 <- expr_data_combined[row_nums_perm_g2, ]
      ###b. calculate delta-p(FFL) of permuted data
      #bootstrap
      boostrap_results_perm_g1 <- boot::boot(data = expr_data_perm_g1, statistic = bootstrap_stat, R = num_bootstrap_samples, ffl_type = ffl_type, parallel = "multicore")
      boostrap_results_perm_g2 <- boot::boot(data = expr_data_perm_g2, statistic = bootstrap_stat, R = num_bootstrap_samples, ffl_type = ffl_type, parallel = "multicore")
      #calculate null delta-p(FFL) values
      null_delta_pffls[perm] <- mean(boostrap_results_perm_g1[["t"]]) - mean(boostrap_results_perm_g2[["t"]])
    }

    #calculate p-value/percentile of observed delta-p(FFL) among null delta-p(FFL)s
    p_val <- mean(abs(null_delta_pffls) >= abs(obs_delta_pffl)) #two-sided test since delta-p(FFL) exists in [-1, 1]
    row[["p_value"]] <- p_val

    #print progress
    if(p_val < 0.05){print(paste0("*****candidate FFL ", i, " / ", nrow(candidate_ffls), " ----- ", "delta-p(FFL) = ", obs_delta_pffl, "; p-value = ",  p_val, sep = ""))}
    else{print(paste0("candidate FFL ", i, " / ", nrow(candidate_ffls), " ----- ", "delta-p(FFL) = ", obs_delta_pffl, "; p-value = ",  p_val, sep = ""))}

    return(row)
  }

  #end parallelization
  stopCluster(cl)

  #assemble results
  #add mirna, tf, targetgene names
  predicted_ffls$targetgene_name <- mapIds(org.Hs.eg.db, keys = predicted_ffls$targetgene, column = c("SYMBOL"), keytype = "ENSEMBL")
  predicted_ffls$tf_name <- mapIds(org.Hs.eg.db, keys = predicted_ffls$tf, column = c("SYMBOL"), keytype = "ENSEMBL")
  predicted_ffls$mirna_name <- miRNA_AccessionToName(predicted_ffls$mirna, targetVersion = "v22")$TargetName
  #adjust for multiple testing
  if(!p_value_adjust_method == "locfdr"){
    predicted_ffls$p_value_adj <- p.adjust(predicted_ffls$p_value, method = p_value_adjust_method)
  }
  if(p_value_adjust_method == "locfdr"){
    results_vec <- fdrtool(predicted_ffls$p_value, statistic = "pvalue", plot = FALSE, color.figure = FALSE, verbose = FALSE)
    predicted_ffls$p_value_adj <- results_vec["lfdr"][[1]]
  }
  #order rows based on adjusted p-value
  predicted_ffls[order(predicted_ffls$p_value_adj), ]
  #order columns
  if(ffl_type == "miRNA"){
    col_order <- c("mirna_name", "tf_name", "targetgene_name", "pFFL_group1", "pFFL_group2", "delta_pFFL_(group1-group2)", "p_value_adj", "p_value",
                   "mirna", "tf", "targetgene", "TARGETSCAN", "MIRTARBASE", "MIRDB",
                   "MIRANDA", "TRRUST", "ENCODE",
                   "alpha0_g1_estimate", "alpha0_g1_se", "alpha0_g1_t", "alpha0_g1_p_value",
                   "alpha1_g1_estimate", "alpha1_g1_se", "alpha1_g1_t", "alpha1_g1_p_value",
                   "beta0_g1_estimate", "beta0_g1_se", "beta0_g1_t", "beta0_g1_p_value",
                   "beta1_g1_estimate", "beta1_g1_se", "beta1_g1_t", "beta1_g1_p_value",
                   "gamma0_g1_estimate", "gamma0_g1_se", "gamma0_g1_t", "gamma0_g1_p_value",
                   "gamma1_g1_estimate", "gamma1_g1_se", "gamma1_g1_t", "gamma1_g1_p_value",
                   "gamma2_g1_estimate", "gamma2_g1_se", "gamma2_g1_t", "gamma2_g1_p_value",
                   "meanX_g1", "meanM_g1", "varX_g1", "varM_g1", "covXM_g1", "var_epsilon_g1",
                   "alpha0_g2_estimate", "alpha0_g2_se", "alpha0_g2_t", "alpha0_g2_p_value",
                   "alpha1_g2_estimate", "alpha1_g2_se", "alpha1_g2_t", "alpha1_g2_p_value",
                   "beta0_g2_estimate", "beta0_g2_se", "beta0_g2_t", "beta0_g2_p_value",
                   "beta1_g2_estimate", "beta1_g2_se", "beta1_g2_t", "beta1_g2_p_value",
                   "gamma0_g2_estimate", "gamma0_g2_se", "gamma0_g2_t", "gamma0_g2_p_value",
                   "gamma1_g2_estimate", "gamma1_g2_se", "gamma1_g2_t", "gamma1_g2_p_value",
                   "gamma2_g2_estimate", "gamma2_g2_se", "gamma2_g2_t", "gamma2_g2_p_value",
                   "meanX_g2", "meanM_g2", "varX_g2", "varM_g2", "covXM_g2", "var_epsilon_g2")
    predicted_ffls <- predicted_ffls[col_order]
  }
  if(ffl_type == "TF"){
    col_order <- c("tf_name", "mirna_name", "targetgene_name", "pFFL_group1", "pFFL_group2", "delta_pFFL_(group1-group2)", "p_value_adj", "p_value",
                   "tf", "mirna", "targetgene", "TRANSMIR", "TARGETSCAN", "MIRTARBASE", "MIRDB",
                   "MIRANDA", "TRRUST", "ENCODE",
                   "alpha0_g1_estimate", "alpha0_g1_se", "alpha0_g1_t", "alpha0_g1_p_value",
                   "alpha1_g1_estimate", "alpha1_g1_se", "alpha1_g1_t", "alpha1_g1_p_value",
                   "beta0_g1_estimate", "beta0_g1_se", "beta0_g1_t", "beta0_g1_p_value",
                   "beta1_g1_estimate", "beta1_g1_se", "beta1_g1_t", "beta1_g1_p_value",
                   "gamma0_g1_estimate", "gamma0_g1_se", "gamma0_g1_t", "gamma0_g1_p_value",
                   "gamma1_g1_estimate", "gamma1_g1_se", "gamma1_g1_t", "gamma1_g1_p_value",
                   "gamma2_g1_estimate", "gamma2_g1_se", "gamma2_g1_t", "gamma2_g1_p_value",
                   "meanX_g1", "meanM_g1", "varX_g1", "varM_g1", "covXM_g1", "var_epsilon_g1",
                   "alpha0_g2_estimate", "alpha0_g2_se", "alpha0_g2_t", "alpha0_g2_p_value",
                   "alpha1_g2_estimate", "alpha1_g2_se", "alpha1_g2_t", "alpha1_g2_p_value",
                   "beta0_g2_estimate", "beta0_g2_se", "beta0_g2_t", "beta0_g2_p_value",
                   "beta1_g2_estimate", "beta1_g2_se", "beta1_g2_t", "beta1_g2_p_value",
                   "gamma0_g2_estimate", "gamma0_g2_se", "gamma0_g2_t", "gamma0_g2_p_value",
                   "gamma1_g2_estimate", "gamma1_g2_se", "gamma1_g2_t", "gamma1_g2_p_value",
                   "gamma2_g2_estimate", "gamma2_g2_se", "gamma2_g2_t", "gamma2_g2_p_value",
                   "meanX_g2", "meanM_g2", "varX_g2", "varM_g2", "covXM_g2", "var_epsilon_g2")
    predicted_ffls <- predicted_ffls[col_order]
  }

  #return predictions
  return(predicted_ffls)
}

