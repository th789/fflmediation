

globalVariables(c("names_mirna_db", "names_tf_db", "names_targetgene_db",
                  "mirna_targetgene_db", "tf_targetgene_db", "i"))

#' @importFrom parallel detectCores
#' @importFrom boot boot
#' @importFrom stats lm
#' @import foreach
#' @import iterators
#' @import doParallel
#' @import snow
#' @import org.Hs.eg.db
#' @import miRBaseConverter
#' @import AnnotationDbi



# functions used ----------------------------------------------------------

#####mediation function (mirna-ffls)
#' @title mediation_ffl
#' @description Determines whether a set of miRNA, TF, and target gene expression data meet the conditions of the mediation model
#' @param mirna Vector of miRNA expression data
#' @param tf Vector of mRNA expression data
#' @param targetgene Vector target gene expression data
#' @param ffl_type Character ("miRNA" or "TF") indicating the FFL type (miRNA-FFL or TF-FFL)
#' @param alpha Significance level of coefficients in the mediation model's linear equations
#' @return Boolean indicating whether the set of miRNA, TF, and target gene expression data meet the conditions of the mediation model
mediation_ffl <- function(mirna, tf, targetgene, ffl_type = c("miRNA", "TF"), alpha){
  ###mirna-ffl
  if(ffl_type == "miRNA"){
    #model1: M ~ X (tf ~ mirna)
    model1 <- lm(tf ~ mirna)
    alpha1 <- summary(model1)$coefficients["mirna", ]
    alpha1_crit <- (alpha1["Pr(>|t|)"] < alpha) & (alpha1["Estimate"] < 0)
    #model2: Y ~ X (targetgene ~ mirna)
    model2 <- lm(targetgene ~ mirna)
    beta1 <- summary(model2)$coefficients["mirna", ]
    beta1_crit <- (beta1["Pr(>|t|)"] < alpha) & (beta1["Estimate"] < 0)
    #model3: Y ~ X + M (targetgene ~ mirna + tf)
    model3 <- lm(targetgene ~ mirna + tf)
    gamma1 <- summary(model3)$coefficients["mirna", ]
    gamma2 <- summary(model3)$coefficients["tf", ]
    gamma1_crit_partialmed <- (gamma1["Pr(>|t|)"] < alpha) & (gamma1["Estimate"] < 0) & (abs(gamma1["Estimate"]) < abs(beta1["Estimate"]))
    gamma2_crit <- (gamma2["Pr(>|t|)"] < alpha) & (gamma2["Estimate"] > 0)
    #determine if candidate ffl has full or partial mediation
    partialmed_crit <- alpha1_crit & beta1_crit & gamma1_crit_partialmed & gamma2_crit
  }

  ###tf-ffls
  if(ffl_type == "TF"){
    #model1: M ~ X (mirna ~ tf)
    model1 <- lm(mirna ~ tf)
    alpha1 <- summary(model1)$coefficients["tf", ]
    alpha1_crit <- (alpha1["Pr(>|t|)"]) < alpha & (alpha1["Estimate"] > 0)
    #model2: Y ~ X (targetgene ~ tf)
    model2 <- lm(targetgene ~ tf)
    beta1 <- summary(model2)$coefficients["tf", ]
    beta1_crit_partialmed <- (beta1["Pr(>|t|)"] < alpha) & (beta1["Estimate"] > 0)
    #model3: Y ~ X + M (gene ~ tf + mirna)
    model3 <- lm(targetgene ~ tf + mirna)
    gamma1 <- summary(model3)$coefficients["tf", ]
    gamma2 <- summary(model3)$coefficients["mirna", ]
    gamma1_crit_partialmed <- (gamma1["Pr(>|t|)"] < alpha) & (gamma1["Estimate"] > 0) & (abs(gamma1["Estimate"]) < abs(beta1["Estimate"]))
    gamma2_crit <- (gamma2["Pr(>|t|)"] < alpha) & (gamma2["Estimate"] < 0)
    #determine if candidate ffl has full or partial mediation
    partialmed_crit <- alpha1_crit & beta1_crit_partialmed & gamma1_crit_partialmed & gamma2_crit
  }

  #return whether loop meets all conditions
  return(partialmed_crit)
}




# full function -----------------------------------------------------------

#' @title predict_ffls
#' @description Predict FFLs using mediation model approach
#' @param mirna_expr Dataframe (miRNA x samples) of miRNA expression data with miRNAs in accession number format
#' @param mrna_expr Dataframe (mRNA x samples) of mRNA expression data with genes in EnsemblID format
#' @param ffl_type Character ("miRNA" or "TF") indicating the FFL type (miRNA-FFL or TF-FFL)
#' @param num_bootstrap_samples Number of bootstrap samples when calculating p(FFL) (default is 1000)
#' @param alpha_mediation Significance level of coefficients in the mediation model's linear equations (default is 0.05)
#' @param num_permutations Number of permutations in permutation test (default is 1000)
#' @param alpha_permutation_test Significance level of FFLs from permutation test (default is 0.05)
#' @param seed random seed for bootstrapping and permutation test (default is 12345)
#' @return Dataframes of candidate FFLs that meet mediation model conditions and are statistically significant)
#' @export

predict_mirna_ffls <- function(mirna_expr, mrna_expr, ffl_type,
                               num_bootstrap_samples, alpha_mediation,
                               num_permutations, alpha_permutation_test,
                               seed){
  #set seed
  set.seed(seed)

  #####0. subset databases
  #identify mirnas, tfs, and targetgenes in both expression data and databases
  mirnas <- rownames(mirna_expr)[rownames(mirna_expr) %in% names_mirna_db]
  tfs <- rownames(mrna_expr)[rownames(mrna_expr) %in% names_tf_db]
  targetgenes <- rownames(mrna_expr)[rownames(mrna_expr) %in% names_targetgene_db]
  #subset databases to include only mirnas, tfs, and targetgenes in expression data
  mirna_targetgene_db_small <- mirna_targetgene_db[mirna_targetgene_db$mirna %in% mirnas & mirna_targetgene_db$targetgene %in% c(tfs, targetgenes), ]
  mirna_tf_db_small <- mirna_targetgene_db_small[mirna_targetgene_db_small$targetgene %in% tfs, ]
  rownames(mirna_tf_db_small) = seq(length = nrow(mirna_tf_db_small)) #reset row numbers
  colnames(mirna_tf_db_small)[colnames(mirna_tf_db_small) == "targetgene"] <- "tf" #rename "targetgene" column to "tf"
  tf_targetgene_db_small <- tf_targetgene_db[tf_targetgene_db$tf %in% tfs & tf_targetgene_db$targetgene %in% c(tfs, targetgenes), ]

  #create function to apply to each row in step2
  mediation_ffl_row <- function(row){
    mirna <- t(mirna_expr[row["mirna"], ])
    tf <- t(mrna_expr[row["tf"], ])
    targetgene <- t(mrna_expr[row["targetgene"], ])
    ffl_meet_conditions <- mediation_ffl(mirna = mirna, tf = tf, targetgene = targetgene, ffl_type = ffl_type, alpha = alpha_mediation)
    return(ffl_meet_conditions)}

  #boostrap statistic function (mirna-ffls)
  bootstrap_stat_mirna <- function(data, indices, mediation_type) {
    d <- data[indices, ] # allows boot to select sample
    #model1: M ~ X (tf ~ mirna)
    model1 <- lm(tf ~ mirna, data = d)
    alpha1 <- summary(model1)$coefficients["mirna", ]
    alpha1_crit <- (alpha1["Pr(>|t|)"] < alpha_mediation) & (alpha1["Estimate"] < 0)
    #model2: Y ~ X (targetgene ~ mirna)
    model2 <- lm(targetgene ~ mirna, data = d)
    beta1 <- summary(model2)$coefficients["mirna", ]
    beta1_crit <- (beta1["Pr(>|t|)"] < alpha_mediation) & (beta1["Estimate"] < 0)
    #model3: Y ~ X + M (targetgene ~ mirna + tf)
    model3 <- lm(targetgene ~ mirna + tf, data = d)
    gamma1 <- summary(model3)$coefficients["mirna", ]
    gamma2 <- summary(model3)$coefficients["tf", ]
    gamma1_crit_fullmed <- gamma1["Pr(>|t|)"] > alpha_mediation
    gamma1_crit_partialmed <- (gamma1["Pr(>|t|)"] < alpha_mediation) & (gamma1["Estimate"] < 0) & (abs(gamma1["Estimate"]) < abs(beta1["Estimate"]))
    gamma2_crit <- (gamma2["Pr(>|t|)"] < alpha_mediation) & (gamma2["Estimate"] > 0)
    #determine if candidate ffl has full or partial mediation
    fullmed_crit <- alpha1_crit & beta1_crit & gamma1_crit_fullmed & gamma2_crit
    partialmed_crit <- alpha1_crit & beta1_crit & gamma1_crit_partialmed & gamma2_crit
    if(mediation_type == "full_mediation"){return(fullmed_crit)}
    if(mediation_type == "partial_mediation"){return(partialmed_crit)}
  }

  #function to apply to each row
  mirna_tf_pair_row_fn <- function(i){
    #####1. assemble candidate ffls
    mirna <- mirna_tf_db_small[i, "mirna"] #mirna in ffl
    tf <- mirna_tf_db_small[i, "tf"] #tf in ffl
    #find targetgenes that are shared by mirna and tf
    mirna_targetgene <- mirna_targetgene_db_small[mirna_targetgene_db_small$mirna == mirna, ] #mirna's targetgenes
    tf_targetgene <- tf_targetgene_db_small[tf_targetgene_db_small$tf == tf, ] #tf's targetgenes
    ffls_candidate <- merge(mirna_targetgene, tf_targetgene, by = "targetgene") #candidate ffls: mirna, tf, shared targetgenes

    if(dim(ffls_candidate)[1] > 0){
      #####2. determine whether candidate ffls meet mediation model conditions
      #apply mediation_ffl_row function to each row
      ffls_candidate$partial_med <- apply(ffls_candidate, 1, mediation_ffl_row)
      #keep candidate ffls that show partial mediation
      ffls_partialmed <- ffls_candidate[ffls_candidate$partial_med, ]

      #####3&4. calculate p(FFL) and p-value for full mediation (linear)
      if(dim(ffls_partialmed)[1] > 0){
        #####3. calculate p(FFL) for candidate ffls that meet mediation model conditions
        ffls_partialmed$p_FFL <- NA #add p_FFL column (to store p(FFL) values)
        ffls_partialmed$p_value <- NA #add p_FFL column (to store p-values)
        #iterate through ffls_fullmed rows
        for(k in 1:nrow(ffls_partialmed)){
          #create expression df (expr_data)
          targetgene <- ffls_partialmed[k, "targetgene"]
          expr_data <- as.data.frame(t(mirna_expr[mirna, ]))
          expr_data$tf <- t(mrna_expr[tf, ])
          expr_data$targetgene <- t(mrna_expr[targetgene, ])
          colnames(expr_data) <- c("mirna", "tf", "targetgene")
          #bootstrap to calculate p(FFL)
          boostrap_results <- boot::boot(data = expr_data, statistic = bootstrap_stat_mirna, R = num_bootstrap_samples,
                                         mediation_type = "partial_mediation", parallel = "multicore")
          #add p(FFL) value to column
          obs_pffl <- mean(boostrap_results[[2]])
          ffls_partialmed[k, "p_FFL"] <- obs_pffl

          #####4. calculate p-value for candidate ffls that meet mediation model conditions
          #create empty vectors to store results
          perm_mirna_vector <- rep(NA, ceiling(num_permutations/4))
          perm_tf_vector <- rep(NA, ceiling(num_permutations/4))
          perm_targetgene_vector <- rep(NA, ceiling(num_permutations/4))
          perm_all_vector <- rep(NA, ceiling(num_permutations/4))
          ###permutation test
          for(perm in 1:ceiling(num_permutations/4)){
            ###a. permute data -- shuffle 1) mirna; 2) tf; 3) targetgene; 4) mirna and tf
            perm_mirna <- data.frame(mirna = sample(expr_data$mirna), tf = expr_data$tf, targetgene = expr_data$targetgene)
            colnames(perm_mirna) <- c("mirna", "tf", "targetgene")
            perm_tf <- data.frame(mirna = expr_data$mirna, tf = sample(expr_data$tf), targetgene = expr_data$targetgene)
            colnames(perm_tf) <- c("mirna", "tf", "targetgene")
            perm_targetgene <- data.frame(mirna = expr_data$mirna, tf = expr_data$tf, targetgene = sample(expr_data$targetgene))
            colnames(perm_targetgene) <- c("mirna", "tf", "targetgene")
            perm_all <- data.frame(mirna = sample(expr_data$mirna), tf = sample(expr_data$tf), targetgene = expr_data$targetgene)
            colnames(perm_all) <- c("mirna", "tf", "targetgene")
            ###b. calculate p(FFL) of permuted data -- through bootstrapping
            boostrap_perm_mirna <- boot::boot(data = perm_mirna, statistic = bootstrap_stat_mirna, R = num_bootstrap_samples,
                                              mediation_type = "partial_mediation", parallel = "multicore")
            boostrap_perm_tf <- boot::boot(data = perm_tf, statistic = bootstrap_stat_mirna, R = num_bootstrap_samples,
                                           mediation_type = "partial_mediation", parallel = "multicore")
            boostrap_perm_targetgene <- boot::boot(data = perm_targetgene, statistic = bootstrap_stat_mirna, R = num_bootstrap_samples,
                                                   mediation_type = "partial_mediation", parallel = "multicore")
            boostrap_perm_all <- boot::boot(data = perm_all, statistic = bootstrap_stat_mirna, R = num_bootstrap_samples,
                                            mediation_type = "partial_mediation", parallel = "multicore")
            #store results
            perm_mirna_vector[perm] <- mean(boostrap_perm_mirna[[2]])
            perm_tf_vector[perm] <- mean(boostrap_perm_tf[[2]])
            perm_targetgene_vector[perm] <- mean(boostrap_perm_targetgene[[2]])
            perm_all_vector[perm] <- mean(boostrap_perm_all[[2]])
          }
          #get p-value/percentile of observed p(FFL) among null p(FFL)s
          perm_pffls <- c(perm_mirna_vector, perm_tf_vector, perm_targetgene_vector, perm_all_vector)
          p_val <- mean(perm_pffls > obs_pffl)
          ffls_partialmed[k, "p_value"] <- p_val
        }

        #return hits
        ffls_hits <- ffls_partialmed[ffls_partialmed$p_value < alpha_permutation_test, ]
        return(ffls_hits)
      }
    }
  }

  #Start a cluster
  num_cores <- detectCores()
  cl <- makeCluster(num_cores, type='SOCK')
  registerDoParallel(cl)

  # Run the loop in parallel
  n <- nrow(mirna_tf_db_small)
  ffl_hits <- foreach(i = 1:n, .combine=rbind) %dopar% {mirna_tf_pair_row_fn(i)}

  #Stop the cluster
  stopCluster(cl)

  return(ffl_hits)
}
