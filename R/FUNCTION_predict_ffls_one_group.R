#from prll_load_02_12


globalVariables(c("i"))

#' @import boot
#' @import stats
#' @import parallel
#' @import foreach
#' @import iterators
#' @import doParallel
#' @import org.Hs.eg.db
#' @import miRBaseConverter
#' @import AnnotationDbi

#one-group function: predict ffls in one biological group

#' @title predict_ffls_one_group
#' @description function to predict FFLs in one biological group
#'
#' @param mirna_expr miRNA expression data (dataframe: miRNAs x samples, see \code{sample_mirna_expr})
#' @param mrna_expr mRNA expression data (dataframe: mRNAs x samples, see \code{sample_mrna_expr})
#' @param ffl_type FFL type (character: "miRNA" or "TF")
#' @param candidate_ffls candidate FFLs (dataframe: candidate FFLs x 10, see \code{sample_candidate_ffls_mirna})
#' @param first_row first row of candidate FFL dataframe to include in analyses (integer: default = 1)
#' @param last_row last row of candidate FFL dataframe to include in analyses (integer: default = nrow(candidate_ffls))
#' @param num_bootstrap_samples number of bootstrap samples (integer: default = 1000)
#' @param num_permutations number of permutations (integer: default = 1000)
#' @param p_value_adjust_method method to adjust p-values for multiple testing (character: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", or "none")
#' @param seed random seed (integer: default = 12345)
#'
#' @return candidate FFLs with pFFLs and p-values (dataframe: candidate FFLs x 16, see sample output)
#' @export


predict_ffls_one_group <- function(mirna_expr, mrna_expr,
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
    mirna <- t(mirna_expr[row[["mirna"]], ])
    tf <- t(mrna_expr[row[["tf"]], ])
    targetgene <- t(mrna_expr[row[["targetgene"]], ])

    #####step3. calculate p(FFL)
    #create expression df (expr_data)
    expr_data <- data.frame(mirna, tf, targetgene)
    colnames(expr_data) <- c("mirna", "tf", "targetgene")
    #bootstrap to calculate p(FFL)
    boostrap_results <- boot::boot(data = expr_data, statistic = bootstrap_stat, R = num_bootstrap_samples, ffl_type = ffl_type, parallel = "multicore")
    obs_pffl <- mean(boostrap_results[["t"]])
    row[["pFFL"]] <- obs_pffl

    #####step4. calculate p-value
    #create empty vectors to store results
    perm_mirna_vector <- rep(NA, ceiling(num_permutations/4))
    perm_tf_vector <- rep(NA, ceiling(num_permutations/4))
    perm_targetgene_vector <- rep(NA, ceiling(num_permutations/4))
    perm_all_vector <- rep(NA, ceiling(num_permutations/4))
    #conduct permutations
    for(perm in 1:ceiling(num_permutations/4)){
      ###a. permute data -- shuffle 1) mirna; 2) tf; 3) targetgene; 4) mirna and tf
      perm_mirna <- data.frame(mirna = sample(expr_data$mirna), tf = expr_data$tf, targetgene = expr_data$targetgene)
      perm_tf <- data.frame(mirna = expr_data$mirna, tf = sample(expr_data$tf), targetgene = expr_data$targetgene)
      perm_targetgene <- data.frame(mirna = expr_data$mirna, tf = expr_data$tf, targetgene = sample(expr_data$targetgene))
      perm_all <- data.frame(mirna = sample(expr_data$mirna), tf = sample(expr_data$tf), targetgene = expr_data$targetgene)
      ###b. calculate p(FFL) of permuted data -- through bootstrapping
      boostrap_perm_mirna <- boot::boot(data = perm_mirna, statistic = bootstrap_stat, R = num_bootstrap_samples, ffl_type = ffl_type, parallel = "multicore")
      boostrap_perm_tf <- boot::boot(data = perm_tf, statistic = bootstrap_stat, R = num_bootstrap_samples, ffl_type = ffl_type, parallel = "multicore")
      boostrap_perm_targetgene <- boot::boot(data = perm_targetgene, statistic = bootstrap_stat, R = num_bootstrap_samples, ffl_type = ffl_type, parallel = "multicore")
      boostrap_perm_all <- boot::boot(data = perm_all, statistic = bootstrap_stat, R = num_bootstrap_samples, ffl_type = ffl_type, parallel = "multicore")
      #store results
      perm_mirna_vector[perm] <- mean(boostrap_perm_mirna[["t"]])
      perm_tf_vector[perm] <- mean(boostrap_perm_tf[["t"]])
      perm_targetgene_vector[perm] <- mean(boostrap_perm_targetgene[["t"]])
      perm_all_vector[perm] <- mean(boostrap_perm_all[["t"]])
    }

    #calculate p-value/percentile of observed p(FFL) among null p(FFL)s
    null_pffls <- c(perm_mirna_vector, perm_tf_vector, perm_targetgene_vector, perm_all_vector)
    p_val <- mean(null_pffls >= obs_pffl)
    row[["p_value"]] <- p_val

    #print results
    if(p_val < 0.05){print(paste0("*****candidate FFL ", i, " / ", last_row, " ----- ", "p(FFL) = ", obs_pffl, "; p-value = ",  p_val, sep = ""))}
    else{print(paste0("candidate FFL ", i, " / ", last_row, " ----- ", "p(FFL) = ", obs_pffl, "; p-value = ",  p_val, sep = ""))}

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
      col_order <- c("mirna_name", "tf_name", "targetgene_name", "pFFL", "p_value_adj", "p_value",
                     "mirna", "tf", "targetgene", "TARGETSCAN", "MIRTARBASE", "MIRDB",
                     "MIRANDA", "TRRUST", "ENCODE")
      predicted_ffls <- predicted_ffls[col_order] #order columns
    }
    if(ffl_type == "TF"){
      col_order <- c("tf_name", "mirna_name", "targetgene_name", "pFFL", "p_value_adj", "p_value",
                     "tf", "mirna", "targetgene", "TRANSMIR", "TARGETSCAN", "MIRTARBASE", "MIRDB",
                     "MIRANDA", "TRRUST", "ENCODE")
      predicted_ffls <- predicted_ffls[col_order] #order columns
    }
  }
  #return prediction
  return(predicted_ffls)
}
