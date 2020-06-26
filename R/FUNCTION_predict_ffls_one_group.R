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
#' @import fdrtool

#one-group function: predict ffls in one biological group

#' @title predict_ffls_one_group
#' @description function to predict FFLs in one biological group
#'
#' @param mirna_expr miRNA expression data (dataframe: miRNAs x samples, see \code{sample_mirna_expr})
#' @param mrna_expr mRNA expression data (dataframe: mRNAs x samples, see \code{sample_mrna_expr})
#' @param ffl_type FFL type (character: "miRNA" or "TF")
#' @param candidate_ffls candidate FFLs (dataframe: candidate FFLs x 10, see \code{sample_candidate_ffls_mirna})
#' @param first_row first row of candidate_ffls dataframe to include in analyses (integer: default = 1)
#' @param last_row last row of candidate_ffls dataframe to include in analyses (integer: default = nrow(candidate_ffls))
#' @param num_bootstrap_samples number of bootstrap samples (integer: default = 1000)
#' @param num_permutations number of permutations (integer: default = 1000)
#' @param p_value_adjust_method method to adjust p-values for multiple testing (character: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "locfdr", or "none")
#' @param seed random seed (integer: default = 12345)
#'
#' @return candidate FFLs with P(FFL) values, p-values, and coefficient estimates (dataframe: # candidate FFLs x 49 for miRNA-FFLs; # candidate FFLs x 50 for TF-FFLs; see sample output)
#' @export


predict_ffls_one_group <- function(mirna_expr, mrna_expr,
                                   ffl_type = c("miRNA", "TF"),
                                   candidate_ffls, first_row = 1, last_row = nrow(candidate_ffls),
                                   num_bootstrap_samples = 1000, num_permutations = 1000,
                                   p_value_adjust_method = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "locfdr", "none"),
                                   seed = 12345){

  #check for errors
  if(ncol(mirna_expr) == 0) stop("mirna_expr is empty (0 rows)")
  if(ncol(mrna_expr) == 0) stop("mrna_expr is empty (0 rows)")
  if(nrow(candidate_ffls) == 0) stop("candidate_ffls is empty (0 rows)")
  if((first_row > nrow(candidate_ffls)) | (last_row > nrow(candidate_ffls))) stop("first_row and/or last_row larger than number of rows in candidate_ffls")
  if(!ncol(mirna_expr) == ncol(mrna_expr)) stop("mirna_expr and mrna_expr do not have the same number of samples (rows)")

  #functions used in foreach: mediation_ffl and bootstrap_stat
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
    meet_conditions <- mediation_ffl(mirna = d$mirna, tf = d$tf, targetgene = d$targetgene, ffl_type = ffl_type)
    return(meet_conditions)
  }

  #start parallelization
  num_cores <- detectCores()-1
  cl <- makeCluster(num_cores, outfile = "")
  registerDoParallel(cl)

  #iterate through each candidate ffl
  candidate_ffls <- candidate_ffls[first_row:last_row, ] #subset candidate ffls
  predicted_ffls <- foreach(i = 1:nrow(candidate_ffls), .combine = rbind) %dopar% {
    ######step1. identify candidate ffl: extract mirna, tf, and targetgene expression levels
    row <- candidate_ffls[i, ] #row is the candidate ffl
    mirna <- t(mirna_expr[row[["mirna"]], ])
    tf <- t(mrna_expr[row[["tf"]], ])
    targetgene <- t(mrna_expr[row[["targetgene"]], ])

    ######step2. calculate regression models, store coefficients
    ###2a. set variables
    if(ffl_type == "miRNA"){X <- mirna; M <- tf; Y <- targetgene}
    if(ffl_type == "TF"){X <- tf; M <- mirna; Y <- targetgene}
    ###2b. fit models
    #model1: Y ~ X
    model1 <- lm(Y ~ X)
    alpha0 <- summary(model1)$coefficients["(Intercept)", ]
    alpha1 <- summary(model1)$coefficients["X", ]
    #model2: M ~ X
    model2 <- lm(M ~ X)
    beta0 <- summary(model2)$coefficients["(Intercept)", ]
    beta1 <- summary(model2)$coefficients["X", ]
    #model3: Y ~ X + M
    model3 <- lm(Y ~ X + M)
    gamma0 <- summary(model3)$coefficients["(Intercept)", ]
    gamma1 <- summary(model3)$coefficients["X", ]
    gamma2 <- summary(model3)$coefficients["M", ]
    ###2c. store information
    #coefficients
    #alpha0
    row[["alpha0_estimate"]] <- alpha0["Estimate"]
    row[["alpha0_se"]] <- alpha0["Std. Error"]
    row[["alpha0_t"]] <- alpha0["t value"]
    row[["alpha0_p_value"]] <- alpha0["Pr(>|t|)"]
    #alpha1
    row[["alpha1_estimate"]] <- alpha1["Estimate"]
    row[["alpha1_se"]] <- alpha1["Std. Error"]
    row[["alpha1_t"]] <- alpha1["t value"]
    row[["alpha1_p_value"]] <- alpha1["Pr(>|t|)"]
    #beta0
    row[["beta0_estimate"]] <- beta0["Estimate"]
    row[["beta0_se"]] <- beta0["Std. Error"]
    row[["beta0_t"]] <- beta0["t value"]
    row[["beta0_p_value"]] <- beta0["Pr(>|t|)"]
    #beta1
    row[["beta1_estimate"]] <- beta1["Estimate"]
    row[["beta1_se"]] <- beta1["Std. Error"]
    row[["beta1_t"]] <- beta1["t value"]
    row[["beta1_p_value"]] <- beta1["Pr(>|t|)"]
    #gamma0
    row[["gamma0_estimate"]] <- gamma0["Estimate"]
    row[["gamma0_se"]] <- gamma0["Std. Error"]
    row[["gamma0_t"]] <- gamma0["t value"]
    row[["gamma0_p_value"]] <- gamma0["Pr(>|t|)"]
    #gamma1
    row[["gamma1_estimate"]] <- gamma1["Estimate"]
    row[["gamma1_se"]] <- gamma1["Std. Error"]
    row[["gamma1_t"]] <- gamma1["t value"]
    row[["gamma1_p_value"]] <- gamma1["Pr(>|t|)"]
    #gamma2
    row[["gamma2_estimate"]] <- gamma2["Estimate"]
    row[["gamma2_se"]] <- gamma2["Std. Error"]
    row[["gamma2_t"]] <- gamma2["t value"]
    row[["gamma2_p_value"]] <- gamma2["Pr(>|t|)"]
    #means
    row[["meanX"]] <- mean(X)
    row[["meanM"]] <- mean(M)
    #variances
    row[["varX"]] <- var(X)*(nrow(X)-1)/nrow(X)
    row[["varM"]] <- var(M)*(nrow(M)-1)/nrow(M)
    row[["covXM"]] <- cov(X, M)*(nrow(X)-1)/nrow(X)
    row[["var_epsilon"]] <- summary(model3)$sigma^2

    #####step3. calculate p(FFL)
    #create expression df (expr_data)
    expr_data <- data.frame(mirna, tf, targetgene)
    colnames(expr_data) <- c("mirna", "tf", "targetgene")
    #bootstrap to calculate p(FFL)
    set.seed(seed)
    boostrap_results <- boot::boot(data = expr_data, statistic = bootstrap_stat, R = num_bootstrap_samples, ffl_type = ffl_type, parallel = "multicore")
    obs_pffl <- mean(boostrap_results[["t"]]) #calculate p(FFL)
    row[["pFFL"]] <- obs_pffl

    #####step4. calculate p-value
    #create empty vectors to store results
    perm_mirna_vector <- rep(NA, ceiling(num_permutations/4))
    perm_tf_vector <- rep(NA, ceiling(num_permutations/4))
    perm_targetgene_vector <- rep(NA, ceiling(num_permutations/4))
    perm_all_vector <- rep(NA, ceiling(num_permutations/4))
    #conduct permutations
    set.seed(seed)
    for(perm in 1:ceiling(num_permutations/4)){
      ###a. permute data -- shuffle 1) mirna; 2) tf; 3) targetgene; 4) mirna and tf
      perm_mirna <- data.frame(mirna = sample(expr_data$mirna), tf = expr_data$tf, targetgene = expr_data$targetgene)
      perm_tf <- data.frame(mirna = expr_data$mirna, tf = sample(expr_data$tf), targetgene = expr_data$targetgene)
      perm_targetgene <- data.frame(mirna = expr_data$mirna, tf = expr_data$tf, targetgene = sample(expr_data$targetgene))
      perm_all <- data.frame(mirna = sample(expr_data$mirna), tf = sample(expr_data$tf), targetgene = expr_data$targetgene)
      ###b. calculate p(FFL) of permuted data -- through bootstrapping
      #bootstrap
      boostrap_perm_mirna <- boot::boot(data = perm_mirna, statistic = bootstrap_stat, R = num_bootstrap_samples, ffl_type = ffl_type, parallel = "multicore")
      boostrap_perm_tf <- boot::boot(data = perm_tf, statistic = bootstrap_stat, R = num_bootstrap_samples, ffl_type = ffl_type, parallel = "multicore")
      boostrap_perm_targetgene <- boot::boot(data = perm_targetgene, statistic = bootstrap_stat, R = num_bootstrap_samples, ffl_type = ffl_type, parallel = "multicore")
      boostrap_perm_all <- boot::boot(data = perm_all, statistic = bootstrap_stat, R = num_bootstrap_samples, ffl_type = ffl_type, parallel = "multicore")
      ##calculate null p(FFL) values
      perm_mirna_vector[perm] <- mean(boostrap_perm_mirna[["t"]])
      perm_tf_vector[perm] <- mean(boostrap_perm_tf[["t"]])
      perm_targetgene_vector[perm] <- mean(boostrap_perm_targetgene[["t"]])
      perm_all_vector[perm] <- mean(boostrap_perm_all[["t"]])
    }

    #calculate p-value/percentile of observed p(FFL) among null p(FFL)s
    null_pffls <- c(perm_mirna_vector, perm_tf_vector, perm_targetgene_vector, perm_all_vector)
    p_val <- mean(null_pffls >= obs_pffl) #one-sided test since p(FFL) exists in [0, 1]
    row[["p_value"]] <- p_val

    #print progress
    if(p_val < 0.05){print(paste0("*****candidate FFL ", i, " / ", nrow(candidate_ffls), " ----- ", "p(FFL) = ", obs_pffl, "; p-value = ",  p_val, sep = ""))}
    else{print(paste0("candidate FFL ", i, " / ", nrow(candidate_ffls), " ----- ", "p(FFL) = ", obs_pffl, "; p-value = ",  p_val, sep = ""))}

    return(row)
  }

  #end parallelization
  stopCluster(cl)

  #####assemble results
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
    col_order <- c("mirna_name", "tf_name", "targetgene_name", "pFFL", "p_value_adj", "p_value",
                   "mirna", "tf", "targetgene", "TARGETSCAN", "MIRTARBASE", "MIRDB",
                   "MIRANDA", "TRRUST", "ENCODE",
                   "alpha0_estimate", "alpha0_se", "alpha0_t", "alpha0_p_value",
                   "alpha1_estimate", "alpha1_se", "alpha1_t", "alpha1_p_value",
                   "beta0_estimate", "beta0_se", "beta0_t", "beta0_p_value",
                   "beta1_estimate", "beta1_se", "beta1_t", "beta1_p_value",
                   "gamma0_estimate", "gamma0_se", "gamma0_t", "gamma0_p_value",
                   "gamma1_estimate", "gamma1_se", "gamma1_t", "gamma1_p_value",
                   "gamma2_estimate", "gamma2_se", "gamma2_t", "gamma2_p_value",
                   "meanX", "meanM", "varX", "varM", "covXM", "var_epsilon")
    predicted_ffls <- predicted_ffls[col_order]
  }
  if(ffl_type == "TF"){
    col_order <- c("tf_name", "mirna_name", "targetgene_name", "pFFL", "p_value_adj", "p_value",
                   "tf", "mirna", "targetgene", "TRANSMIR", "TARGETSCAN", "MIRTARBASE", "MIRDB",
                   "MIRANDA", "TRRUST", "ENCODE",
                   "alpha0_estimate", "alpha0_se", "alpha0_t", "alpha0_p_value",
                   "alpha1_estimate", "alpha1_se", "alpha1_t", "alpha1_p_value",
                   "beta0_estimate", "beta0_se", "beta0_t", "beta0_p_value",
                   "beta1_estimate", "beta1_se", "beta1_t", "beta1_p_value",
                   "gamma0_estimate", "gamma0_se", "gamma0_t", "gamma0_p_value",
                   "gamma1_estimate", "gamma1_se", "gamma1_t", "gamma1_p_value",
                   "gamma2_estimate", "gamma2_se", "gamma2_t", "gamma2_p_value",
                   "meanX", "meanM", "varX", "varM", "covXM", "var_epsilon")
    predicted_ffls <- predicted_ffls[col_order]
  }

  #return predictions
  return(predicted_ffls)
}
