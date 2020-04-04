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
#' @param first_row first row of candidate FFL dataframe to include in analyses (integer: default = 1)
#' @param last_row last row of candidate FFL dataframe to include in analyses (integer: default = nrow(candidate_ffls))
#' @param num_bootstrap_samples number of bootstrap samples (integer: default = 1000)
#' @param num_permutations number of permutations (integer: default = 1000)
#' @param p_value_adjust_method method to adjust p-values for multiple testing (character: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", or "none")
#' @param seed random seed (integer: default = 12345)
#'
#' @return candidate FFLs with delta-pFFLs and p-values (dataframe: candidate FFLs x 18, see sample output)
#' @export


predict_ffls_two_groups <- function(mirna_expr_g1, mrna_expr_g1,
                                    mirna_expr_g2, mrna_expr_g2,
                                    ffl_type = c("miRNA", "TF"),
                                    candidate_ffls, first_row = 1, last_row = nrow(candidate_ffls),
                                    num_bootstrap_samples = 1000, num_permutations = 1000,
                                    p_value_adjust_method = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"),
                                    seed = 12345){

  #functions used in foreach
  #####mediation function (mirna-ffls)
  mediation_ffl <- function(mirna, tf, targetgene, ffl_type = c("miRNA", "TF")){
    ###mirna-ffl
    if(ffl_type == "miRNA"){
      #model1: M ~ X (tf ~ mirna)
      model1 <- lm(tf ~ mirna)
      alpha1 <- summary(model1)$coefficients["mirna", ]
      alpha1_crit <- alpha1["Estimate"] < 0
      #model2: Y ~ X (targetgene ~ mirna)
      model2 <- lm(targetgene ~ mirna)
      beta1 <- summary(model2)$coefficients["mirna", ]
      beta1_crit <- beta1["Estimate"] < 0
      #model3: Y ~ X + M (targetgene ~ mirna + tf)
      model3 <- lm(targetgene ~ mirna + tf)
      gamma1 <- summary(model3)$coefficients["mirna", ]
      gamma2 <- summary(model3)$coefficients["tf", ]
      gamma1_crit <- (gamma1["Estimate"] < 0) & (abs(gamma1["Estimate"]) < abs(beta1["Estimate"]))
      gamma2_crit <- gamma2["Estimate"] > 0
    }

    ###tf-ffls
    if(ffl_type == "TF"){
      #model1: M ~ X (mirna ~ tf)
      model1 <- lm(mirna ~ tf)
      alpha1 <- summary(model1)$coefficients["tf", ]
      alpha1_crit <- alpha1["Estimate"] > 0
      #model2: Y ~ X (targetgene ~ tf)
      model2 <- lm(targetgene ~ tf)
      beta1 <- summary(model2)$coefficients["tf", ]
      beta1_crit <- beta1["Estimate"] > 0
      #model3: Y ~ X + M (gene ~ tf + mirna)
      model3 <- lm(targetgene ~ tf + mirna)
      gamma1 <- summary(model3)$coefficients["tf", ]
      gamma2 <- summary(model3)$coefficients["mirna", ]
      gamma1_crit <- (gamma1["Estimate"] > 0) & (abs(gamma1["Estimate"]) > abs(beta1["Estimate"]))
      gamma2_crit <- gamma2["Estimate"] < 0
    }

    #determine if candidate ffl meets mediation conditions
    meet_conditions <- alpha1_crit & beta1_crit & gamma1_crit & gamma2_crit
    return(meet_conditions)
  }
  #####bootstrap function
  bootstrap_stat <- function(data, indices, ffl_type){
    if(ffl_type == "miRNA"){
      d <- data[indices, ] #allows boot to select sample
      #model1: M ~ X (tf ~ mirna)
      model1 <- lm(tf ~ mirna, data = d)
      alpha1 <- summary(model1)$coefficients["mirna", ]
      alpha1_crit <- alpha1["Estimate"] < 0
      #model2: Y ~ X (targetgene ~ mirna)
      model2 <- lm(targetgene ~ mirna, data = d)
      beta1 <- summary(model2)$coefficients["mirna", ]
      beta1_crit <- beta1["Estimate"] < 0
      #model3: Y ~ X + M (targetgene ~ mirna + tf)
      model3 <- lm(targetgene ~ mirna + tf, data = d)
      gamma1 <- summary(model3)$coefficients["mirna", ]
      gamma2 <- summary(model3)$coefficients["tf", ]
      gamma1_crit <- (gamma1["Estimate"] < 0) & (abs(gamma1["Estimate"]) < abs(beta1["Estimate"]))
      gamma2_crit <- gamma2["Estimate"] > 0
    }
    if(ffl_type == "TF"){
      d <- data[indices, ] #allows boot to select sample
      #model1: M ~ X (mirna ~ tf)
      model1 <- lm(mirna ~ tf, data = d)
      alpha1 <- summary(model1)$coefficients["tf", ]
      alpha1_crit <- alpha1["Estimate"] > 0
      #model2: Y ~ X (targetgene ~ tf)
      model2 <- lm(targetgene ~ tf, data = d)
      beta1 <- summary(model2)$coefficients["tf", ]
      beta1_crit <- beta1["Estimate"] > 0
      #model3: Y ~ X + M (targetgene ~ tf + mirna)
      model3 <- lm(targetgene ~ tf + mirna, data = d)
      gamma1 <- summary(model3)$coefficients["tf", ]
      gamma2 <- summary(model3)$coefficients["mirna", ]
      gamma1_crit <- (gamma1["Estimate"] > 0) & (abs(gamma1["Estimate"]) > abs(beta1["Estimate"]))
      gamma2_crit <- gamma2["Estimate"] < 0
    }
    #determine if candidate ffl meets mediation model conditions
    meet_conditions <- alpha1_crit & beta1_crit & gamma1_crit & gamma2_crit
    return(meet_conditions)
  }

  #start parallelization
  num_cores <- detectCores()
  cl <- makeCluster(num_cores, outfile = "")
  registerDoParallel(cl)

  set.seed(seed)
  predicted_ffls <- foreach(i = first_row:last_row, .combine = rbind) %dopar% {
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

    #####step3. calculate delta-p(FFL)
    #create expression df (expr_data)
    expr_data_g1 <- data.frame(mirna_g1, tf_g1, targetgene_g1) #group1
    colnames(expr_data_g1) <- c("mirna", "tf", "targetgene")
    expr_data_g2 <- data.frame(mirna_g2, tf_g2, targetgene_g2) #group2
    colnames(expr_data_g2) <- c("mirna", "tf", "targetgene")

    #bootstrap to calculate delta-p(FFL)
    boostrap_results_g1 <- boot::boot(data = expr_data_g1, statistic = bootstrap_stat, R = num_bootstrap_samples, ffl_type = ffl_type, parallel = "multicore")
    boostrap_results_g2 <- boot::boot(data = expr_data_g2, statistic = bootstrap_stat, R = num_bootstrap_samples, ffl_type = ffl_type, parallel = "multicore")
    obs_delta_pffl <- mean(boostrap_results_g1[["t"]]) - mean(boostrap_results_g2[["t"]])
    row[["pFFL_group1"]] <- mean(boostrap_results_g1[["t"]])
    row[["pFFL_group2"]] <- mean(boostrap_results_g2[["t"]])
    row[["delta_pFFL_(group1-group2)"]] <- obs_delta_pffl

    #####step4. calculate p-value
    #create empty vectors to store results
    perm_mirna_vector <- rep(NA, ceiling(num_permutations/4))
    perm_tf_vector <- rep(NA, ceiling(num_permutations/4))
    perm_targetgene_vector <- rep(NA, ceiling(num_permutations/4))
    perm_all_vector <- rep(NA, ceiling(num_permutations/4))

    #conduct permutations
    for(perm in 1:ceiling(num_permutations/4)){
      ###a. permute data -- shuffle 1) mirna; 2) tf; 3) targetgene; 4) mirna and tf
      perm_mirna <- list(data.frame(mirna = sample(expr_data_g1$mirna), tf = expr_data_g1$tf, targetgene = expr_data_g1$targetgene),
                         data.frame(mirna = sample(expr_data_g2$mirna), tf = expr_data_g2$tf, targetgene = expr_data_g2$targetgene))
      perm_tf <- list(data.frame(mirna = expr_data_g1$mirna, tf = sample(expr_data_g1$tf), targetgene = expr_data_g1$targetgene),
                      data.frame(mirna = expr_data_g2$mirna, tf = sample(expr_data_g2$tf), targetgene = expr_data_g2$targetgene))
      perm_targetgene <- list(data.frame(mirna = expr_data_g1$mirna, tf = expr_data_g1$tf, targetgene = sample(expr_data_g1$targetgene)),
                              data.frame(mirna = expr_data_g2$mirna, tf = expr_data_g2$tf, targetgene = sample(expr_data_g2$targetgene)))
      perm_all <- list(data.frame(mirna = sample(expr_data_g1$mirna), tf = sample(expr_data_g1$tf), targetgene = expr_data_g1$targetgene),
                       data.frame(mirna = sample(expr_data_g2$mirna), tf = sample(expr_data_g2$tf), targetgene = expr_data_g2$targetgene))

      ###b. calculate delta-p(FFL) of permuted data -- through bootstrapping
      boostrap_perm_mirna <- list(boot::boot(data = perm_mirna[[1]], statistic = bootstrap_stat, R = num_bootstrap_samples, ffl_type = ffl_type, parallel = "multicore"),
                                  boot::boot(data = perm_mirna[[2]], statistic = bootstrap_stat, R = num_bootstrap_samples, ffl_type = ffl_type, parallel = "multicore"))
      boostrap_perm_tf <- list(boot::boot(data = perm_tf[[1]], statistic = bootstrap_stat, R = num_bootstrap_samples, ffl_type = ffl_type, parallel = "multicore"),
                               boot::boot(data = perm_tf[[2]], statistic = bootstrap_stat, R = num_bootstrap_samples, ffl_type = ffl_type, parallel = "multicore"))
      boostrap_perm_targetgene <- list(boot::boot(data = perm_targetgene[[1]], statistic = bootstrap_stat, R = num_bootstrap_samples, ffl_type = ffl_type, parallel = "multicore"),
                                       boot::boot(data = perm_targetgene[[2]], statistic = bootstrap_stat, R = num_bootstrap_samples, ffl_type = ffl_type, parallel = "multicore"))
      boostrap_perm_all <- list(boot::boot(data = perm_all[[1]], statistic = bootstrap_stat, R = num_bootstrap_samples, ffl_type = ffl_type, parallel = "multicore"),
                                boot::boot(data = perm_all[[2]], statistic = bootstrap_stat, R = num_bootstrap_samples, ffl_type = ffl_type, parallel = "multicore"))
      #store results
      perm_mirna_vector[perm] <- mean(boostrap_perm_mirna[[1]][["t"]]) - mean(boostrap_perm_mirna[[2]][["t"]])
      perm_tf_vector[perm] <- mean(boostrap_perm_tf[[1]][["t"]]) - mean(boostrap_perm_tf[[2]][["t"]])
      perm_targetgene_vector[perm] <- mean(boostrap_perm_targetgene[[1]][["t"]]) - mean(boostrap_perm_targetgene[[2]][["t"]])
      perm_all_vector[perm] <- mean(boostrap_perm_all[[1]][["t"]]) - mean(boostrap_perm_all[[2]][["t"]])
    }

    #calculate p-value/percentile of observed p(FFL) among null p(FFL)s (2-sided permutation p-value)
    null_delta_pffls <- c(perm_mirna_vector, perm_tf_vector, perm_targetgene_vector, perm_all_vector)
    p_val <- mean(abs(null_delta_pffls) >= abs(obs_delta_pffl))
    row[["p_value"]] <- p_val

    #print results
    if(p_val < 0.05){print(paste0("*****candidate FFL ", i, " / ", last_row, " ----- ", "delta-p(FFL) = ", obs_delta_pffl, "; p-value = ",  p_val, sep = ""))}
    else{print(paste0("candidate FFL ", i, " / ", last_row, " ----- ", "delta-p(FFL) = ", obs_delta_pffl, "; p-value = ",  p_val, sep = ""))}

    return(row)

  }
  #end parallelization
  stopCluster(cl)

  #add mirna, tf, targetgene names
  if(dim(predicted_ffls)[1] > 0){
    predicted_ffls$targetgene_name <- mapIds(org.Hs.eg.db, keys = predicted_ffls$targetgene, column = c("SYMBOL"), keytype = "ENSEMBL")
    predicted_ffls$tf_name <- mapIds(org.Hs.eg.db, keys = predicted_ffls$tf, column = c("SYMBOL"), keytype = "ENSEMBL")
    predicted_ffls$mirna_name <- miRNA_AccessionToName(predicted_ffls$mirna, targetVersion = "v22")$TargetName
    predicted_ffls$p_value_adj <- p.adjust(predicted_ffls$p_value, method = p_value_adjust_method)
    predicted_ffls[order(predicted_ffls$p_value_adj), ] #order rows based on adjusted p-value
    #order columns
    if(ffl_type == "miRNA"){
      col_order <- c("mirna_name", "tf_name", "targetgene_name", "pFFL_group1", "pFFL_group2", "delta_pFFL_(group1-group2)", "p_value_adj", "p_value",
                     "mirna", "tf", "targetgene", "TARGETSCAN", "MIRTARBASE", "MIRDB",
                     "MIRANDA", "TRRUST", "ENCODE")
      predicted_ffls <- predicted_ffls[col_order] #order columns
    }
    if(ffl_type == "TF"){
      col_order <- c("tf_name", "mirna_name", "targetgene_name", "pFFL_group1", "pFFL_group2", "delta_pFFL_(group1-group2)", "p_value_adj", "p_value",
                     "tf", "mirna", "targetgene", "TRANSMIR", "TARGETSCAN", "MIRTARBASE", "MIRDB",
                     "MIRANDA", "TRRUST", "ENCODE")
      predicted_ffls <- predicted_ffls[col_order] #order columns
    }
  }
  #return prediction
  return(predicted_ffls)
}

