globalVariables(c("mirna_targetgene_db", "tf_mirna_db", "tf_targetgene_db",
                  "names_mirna_db", "names_tf_db", "names_targetgene_db",
                  "mirna", "targetgene",
                  "lm", "mapIds"))
#step1
#' @import dplyr
#for left_join, transmute
#' @import tidyr
#for replace_na
#' @import org.Hs.eg.db
#for mapIds
#' @import miRBaseConverter
#for miRNA_AccessionToName


#step2
#' @importFrom stats lm


# 1. generate list of candidate ffls --------------------------------------

#' @title step1_candidate_ffls
#' @description Generate candidate miRNA-FFLs and TF-FFLs
#' @param mirna_expr Dataframe (miRNA x samples) of miRNA expression data with miRNAs in accession number format
#' @param mrna_expr Dataframe (mRNA x samples) of mRNA expression data with genes in EnsemblID format
#' @param ffl_type Character ("miRNA" or "TF") indicating the FFL type (miRNA-FFL or TF-FFL)
#' @return Vector of two dataframes (miRNA-FFLs and TF-FFLs)
#####step1_candidate_ffls
step1_candidate_ffls <- function(mirna_expr, mrna_expr, ffl_type = c("miRNA", "TF")){
  #####extract mirnas, tfs, and targetgenes from expression data
  names_mirna <- rownames(mirna_expr)[rownames(mirna_expr) %in% names_mirna_db]
  names_tf <- rownames(mrna_expr)[rownames(mrna_expr) %in% names_tf_db]
  names_targetgene <- rownames(mrna_expr)[rownames(mrna_expr) %in% names_targetgene_db]

  #####candidate mirna-ffls
  #mirna-tf arm
  if(ffl_type == "miRNA"){
    ffls_mirna <- mirna_targetgene_db[mirna_targetgene_db$mirna %in% names_mirna & mirna_targetgene_db$targetgene %in% names_tf, ]
    names(ffls_mirna)[names(ffls_mirna) == "targetgene"] <- "tf" #rename "targetgene" col to "tf" since all targetgenes are tfs (based on the line above)
    #tf-targetgene arm
    #!!!!!subset of data -- need to efficiently merge both dfs (ffls_mirna & tf_targetgene_db)
    ffls_mirna <- merge(ffls_mirna, tf_targetgene_db[1:1000, ], by = "tf")
    ffls_mirna <- ffls_mirna[ffls_mirna$targetgene %in% names_targetgene, ] #keep only triplets where targetgene is in the expr. data
    #mirna-targetgene arm
    ffls_mirna <- ffls_mirna %>%
      left_join(mirna_targetgene_db %>% transmute(mirna, targetgene, closed_loop = "yes"), by = c("mirna", "targetgene")) %>%
      replace_na(list(closed_loop = "no")) #check whether mirna-targetgene pair of triplet is in mirna_targetgene_db
    #keep candidate ffls
    ffls_mirna <- ffls_mirna[ffls_mirna$closed_loop == "yes", ]
    #add mirna_name, tf_symbol, targetgene_symbol columns
    mirna_name_df <- miRNA_AccessionToName(ffls_mirna$mirna, targetVersion = "v22") #mirna
    ffls_mirna$mirna_name <- mirna_name_df$TargetName
    ffls_mirna$tf_symbol <- mapIds(org.Hs.eg.db, keys = ffls_mirna$tf, column = c("SYMBOL"), keytype = "ENSEMBL") #tf
    ffls_mirna$targetgene_symbol <- mapIds(org.Hs.eg.db, keys = ffls_mirna$targetgene, column = c("SYMBOL"), keytype = "ENSEMBL") #targetgene
    #re-order columns, drop "closed_loop" column
    ffls_mirna <- ffls_mirna[c("mirna", "tf", "targetgene", "mirna_name", "tf_symbol", "targetgene_symbol", "TARGETSCAN", "MIRTARBASE", "MIRDB", "MIRANDA", "TRRUST", "ENCODE")]
    #####return candidate ffls
    print(paste0(dim(ffls_mirna)[1], " candidate miRNA-FFLs"))
    return(ffls_mirna)
  }

  #####candidate tf-ffls
  #tf-mirna arm
  if(ffl_type == "TF"){
    ffls_tf <- tf_mirna_db[tf_mirna_db$tf %in% names_tf & tf_mirna_db$mirna %in% names_mirna, ]
    #tf-targetgene arm
    #!!!!!subset of data -- need to efficiently merge both dfs (ffls_tf & tf_targetgene_db)
    ffls_tf <- merge(ffls_tf, tf_targetgene_db[7000:10000, ], by = "tf")
    ffls_tf <- ffls_tf[ffls_tf$targetgene %in% names_targetgene, ] #keep only triplets where targetgene is in the expr. data
    #mirna-targetgene arm
    ffls_tf <- ffls_tf %>%
      left_join(mirna_targetgene_db %>% transmute(mirna, targetgene, closed_loop = "yes"), by = c("mirna", "targetgene")) %>%
      replace_na(list(closed_loop = "no")) #check whether mirna-targetgene pair of triplet is in mirna_targetgene_db
    #keep candidate ffls
    ffls_tf <- ffls_tf[ffls_tf$closed_loop == "yes", ]
    #add mirna_name, tf_symbol, targetgene_symbol columns
    mirna_name_df <- miRNA_AccessionToName(ffls_tf$mirna, targetVersion = "v22") #mirna
    ffls_tf$mirna_name <- mirna_name_df$TargetName
    ffls_tf$tf_symbol <- mapIds(org.Hs.eg.db, keys = ffls_tf$tf, column = c("SYMBOL"), keytype = "ENSEMBL") #tf
    ffls_tf$targetgene_symbol <- mapIds(org.Hs.eg.db, keys = ffls_tf$targetgene, column = c("SYMBOL"), keytype = "ENSEMBL") #targetgene
    #re-order columns, drop "closed_loop" column
    ffls_tf <- ffls_tf[c("tf", "mirna", "targetgene", "tf_symbol", "mirna_name", "targetgene_symbol", "TRANSMIR", "TRRUST", "ENCODE")]
    #####return candidate ffls
    print(paste0(dim(ffls_tf)[1], " candidate TF-FFLs"))
    return(ffls_tf)
  }
}
#####fin


# 2. identify FFLs that meet mediation model conditions -------------------------------

#' @title mediation_ffl
#' @description Determines whether a set of miRNA, TF, and target gene expression data meet the conditions of the mediation model
#' @param mirna Vector of miRNA expression data
#' @param tf Vector of mRNA expression data
#' @param targetgene Vector target gene expression data
#' @param ffl_type Character ("miRNA" or "TF") indicating the FFL type (miRNA-FFL or TF-FFL)
#' @param alpha Significance level of coefficients in the mediation model's linear equations
#' @return Boolean indicating whether the set of miRNA, TF, and target gene expression data meet the conditions of the mediation model

#####mediation model function: does a set of mirna, tf, targetgene data meet mediation model conditions?
mediation_ffl <- function(mirna, tf, targetgene, ffl_type = c("miRNA", "TF"), alpha){
  ###mirna-ffl
  if(ffl_type == "miRNA"){
    #model1: tf ~ mirna
    model1 <- lm(tf ~ mirna)
    alpha1 <- summary(model1)$coefficients["mirna", ]
    model1_crit <- alpha1["Estimate"] < 0 & alpha1["Pr(>|t|)"] < alpha
    #model2: targetgene ~ mirna
    model2 <- lm(targetgene ~ mirna)
    beta1 <- summary(model2)$coefficients["mirna", ]
    model2_crit <- beta1["Estimate"] < 0 & beta1["Pr(>|t|)"] < alpha
    #model3: gene ~ mirna + tf
    model3 <- lm(targetgene ~ mirna + tf)
    gamma1 <- summary(model3)$coefficients["mirna", ]
    gamma2 <- summary(model3)$coefficients["tf", ]
    model3_crit <- (gamma2["Estimate"] > 0 & gamma2["Pr(>|t|)"] < alpha) &
      (gamma1["Estimate"] < 0 & gamma1["Pr(>|t|)"] < alpha & (abs(gamma1["Estimate"]) < abs(alpha1["Estimate"])))
  }
  ###tf-ffls
  if(ffl_type == "TF"){
    #model1: mirna ~ tf
    model1 <- lm(mirna ~ tf)
    alpha1 <- summary(model1)$coefficients["tf", ]
    model1_crit <- alpha1["Estimate"] > 0 & alpha1["Pr(>|t|)"] < alpha
    #model2: targetgene ~ tf
    model2 <- lm(targetgene ~ tf)
    beta1 <- summary(model2)$coefficients["tf", ]
    model2_crit <- beta1["Estimate"] > 0 & beta1["Pr(>|t|)"] < alpha
    #model3: gene ~ tf + mirna
    model3 <- lm(targetgene ~ tf + mirna)
    gamma1 <- summary(model3)$coefficients["tf", ]
    gamma2 <- summary(model3)$coefficients["mirna", ]
    model3_crit <- (gamma2["Estimate"] < 0 & gamma2["Pr(>|t|)"] < alpha) &
      (gamma1["Estimate"] > 0 & gamma1["Pr(>|t|)"] < alpha & (abs(gamma1["Estimate"]) < abs(alpha1["Estimate"])))
  }
  #return whether loop meets all conditions
  return(model1_crit & model2_crit & model3_crit)
}
#####fin


#' @title step2_mediation
#' @description Apply mediation model to each candidate FFL
#' @param mirna_expr Dataframe (miRNA x samples) of miRNA expression data with miRNAs in accession number format
#' @param mrna_expr Dataframe (mRNA x samples) of mRNA expression data with genes in EnsemblID format
#' @param candidate_ffls Dataframe of candidate ffls (output from \code{\link{step1_candidate_ffls}})
#' @param ffl_type Character ("miRNA" or "TF") indicating the FFL type (miRNA-FFL or TF-FFL)
#' @param alpha Significance level of coefficients in the mediation model's linear equations
#' @return Dataframe of FFLs among candidate FFLs that meet the mediation model's conditions
#####step2_mediation
step2_mediation <- function(mirna_expr, mrna_expr,
                            candidate_ffls, ffl_type = c("miRNA", "TF"),
                            alpha){
  #####function for each row: does ffl meet mediation model conditions?
  mediation_ffl_row <- function(row){
    mirna <- t(mirna_expr[row["mirna"], ])
    tf <- t(mrna_expr[row["tf"], ])
    targetgene <- t(mrna_expr[row["targetgene"], ])
    ffl_meet_conditions <- mediation_ffl(mirna = mirna, tf = tf, targetgene = targetgene, ffl_type = ffl_type, alpha = alpha)
    return(ffl_meet_conditions)}

  #####see if each ffl meets mediation model criteria
  candidate_ffls$mediation_analysis <- apply(candidate_ffls, 1, mediation_ffl_row)
  #return candidate ffls that meet criteria
  ffls_mediation <- candidate_ffls[candidate_ffls$mediation_analysis, ]
  print(paste0(dim(ffls_mediation)[1], "/", dim(candidate_ffls)[1], " candidate ", ffl_type, "-FFLs meet mediation model conditions"))
  return(ffls_mediation)
}
#####fin


# 3. calculate p(FFL) through bootstrapping -------------------------------

#' @title step3_pffl
#' @description Calculate p(FFL) for FFLs that meet mediation model conditions through bootstrapping
#' @param mirna_expr Dataframe (miRNA x samples) of miRNA expression data with miRNAs in accession number format
#' @param mrna_expr Dataframe (mRNA x samples) of mRNA expression data with genes in EnsemblID forma
#' @param ffls Dataframe of FFLs that meet mediation model conditions (output of step2_mediation)
#' @param ffl_type Character ("miRNA" or "TF") indicating the FFL type (miRNA-FFL or TF-FFL
#' @param num_bootstrap_samples Number of bootstrap samples when calculating p(FFL) (default is 1000)
#' @param seed Random seed
#' @param alpha Significance level of coefficients in the mediation model's linear equations
#' @return \code{ffls} dataframe with added column of p(FFL) values

#####step3_pffl
step3_pffl <- function(mirna_expr, mrna_expr, ffls, ffl_type = c("miRNA", "TF"),
                       num_bootstrap_samples, seed, alpha){
  #function to apply to each row
  step3_bootstrap <- function(row){
    #vector to store result of each bootstrap sample
    bootstrap_results <- rep(NA, num_bootstrap_samples)
    #bootstrapping
    for(i in 1:num_bootstrap_samples){
      #expression data
      mirna <- t(mirna_expr[row["mirna"], ])
      tf <- t(mrna_expr[row["tf"], ])
      targetgene <- t(mrna_expr[row["targetgene"], ])
      #bootstrap sample
      sampleIDs <- colnames(mirna_expr)
      sampleIDs_boot <- sample(x = sampleIDs, size = length(sampleIDs), replace = TRUE)
      mirna_boot <- mirna[sampleIDs_boot, ]
      tf_boot <- tf[sampleIDs_boot, ]
      targetgene_boot <- targetgene[sampleIDs_boot, ]
      #apply mediation model
      bootstrap_results[i] <- mediation_ffl(mirna = mirna_boot, tf = tf_boot, targetgene = targetgene_boot, ffl_type = ffl_type, alpha = alpha)
    }
    #return p(FFL)
    return(mean(bootstrap_results))
  }
  #apply function to every row
  set.seed(seed)
  ffls$p_ffl <- apply(ffls, 1, step3_bootstrap)
  return(ffls)
}
#####fin



# 4. calculate statistical significance through permutation test ----------

#' @title step4_permutation_test
#' @description Calculate p-value associated with p(FFL) for FFLs that meet mediation model conditions through a permutation test
#' @param mirna_expr Dataframe (miRNA x samples) of miRNA expression data with miRNAs in accession number format
#' @param mrna_expr Dataframe (mRNA x samples) of mRNA expression data with genes in EnsemblID forma
#' @param ffls Dataframe of FFLs that meet mediation model conditions (output of step2_mediation)
#' @param ffl_type Character ("miRNA" or "TF") indicating the FFL type (miRNA-FFL or TF-FFL
#' @param num_permutations Number of permutations in permutation test (default is 1000)
#' @param num_bootstrap_samples Number of bootstrap samples when calculating p(FFL) (default is 1000)
#' @param alpha Significance level of coefficients in the mediation model's linear equations
#' @param seed Random seed
#' @return \code{ffls} dataframe with added column of p-values associated with p(FFL) values

#####step4_permutation_test
step4_permutation_test <- function(mirna_expr, mrna_expr, ffls, ffl_type = c("miRNA", "TF"),
                                   num_permutations, num_bootstrap_samples, alpha, seed){
  #####step4_bootstrap function: calculate p(FFL) of permuted data (used in step4_ffl_pval function)
  step4_bootstrap <- function(expr_df, num_bootstrap_samples, ffl_type = c("miRNA", "TF"), alpha){
    #set.seed(seed) #set seed in step4_permutation_test function
    #vector to store result of each bootstrap sample
    bootstrap_results <- rep(NA, num_bootstrap_samples)
    #bootstrapping
    for(i in 1:num_bootstrap_samples){
      #bootstrap sample
      rows <- 1:nrow(expr_df)
      rows_boot <- sample(x = rows, size = length(rows), replace = TRUE)
      expr_boot <- expr_df[rows_boot, ]
      #apply mediation model, store results
      bootstrap_results[i] <- mediation_ffl(mirna = expr_boot$mirna, tf = expr_boot$tf, targetgene = expr_boot$targetgene, ffl_type = ffl_type, alpha = alpha)
    }
    #return p(FFL)
    return(mean(bootstrap_results))
  }

  #####step4_ffl_pval function: conduct permutation test for each ffl
  step4_ffl_pval <- function(row){
    #make expr (df with mirna, tf, and targetgene values for a row's ffl)
    mirna <- t(mirna_expr[row["mirna"], ])
    tf <- t(mrna_expr[row["tf"], ])
    targetgene <- t(mrna_expr[row["targetgene"], ])
    expr <- data.frame(mirna = mirna, tf = tf, targetgene = targetgene)
    colnames(expr) <- c("mirna", "tf", "targetgene")

    #empty vectors to store p(FFL) values of permuted data
    pffl_perm_mirna_vector <- rep(NA, ceiling(num_permutations/4))
    pffl_perm_tf_vector <- rep(NA, ceiling(num_permutations/4))
    pffl_perm_targetgene_vector <- rep(NA, ceiling(num_permutations/4))
    pffl_perm_all_vector <- rep(NA, ceiling(num_permutations/4))

    ###permutation test
    for(i in 1:ceiling(num_permutations/4)){
      ###1. permute data -- shuffle 1) mirna; 2) tf; 3) targetgene; 4) mirna and tf
      expr_perm_mirna <- data.frame(mirna = sample(expr$mirna), tf = expr$tf, targetgene = expr$targetgene)
      expr_perm_tf <- data.frame(mirna = expr$mirna, tf = sample(expr$tf), targetgene = expr$targetgene)
      expr_perm_targetgene <- data.frame(mirna = expr$mirna, tf = expr$tf, targetgene = sample(expr$targetgene))
      expr_perm_all <- data.frame(mirna = sample(expr$mirna), tf = sample(expr$tf), targetgene = expr$targetgene)
      ###2. calculate p(FFL) of permuted data -- through bootstrapping
      pffl_perm_mirna <- step4_bootstrap(expr_perm_mirna, num_bootstrap_samples, ffl_type = ffl_type, alpha = alpha)
      pffl_perm_tf <- step4_bootstrap(expr_perm_tf, num_bootstrap_samples, ffl_type = ffl_type, alpha = alpha)
      pffl_perm_targetgene <- step4_bootstrap(expr_perm_targetgene, num_bootstrap_samples, ffl_type = ffl_type, alpha = alpha)
      pffl_perm_all <- step4_bootstrap(expr_perm_all, num_bootstrap_samples, ffl_type = ffl_type, alpha = alpha)
      #store results
      pffl_perm_mirna_vector[i] <- pffl_perm_mirna
      pffl_perm_tf_vector[i] <- pffl_perm_tf
      pffl_perm_targetgene_vector[i] <- pffl_perm_targetgene
      pffl_perm_all_vector[i] <- pffl_perm_all
    }
    #get p-value/percentile of observed p(FFL) among null p(FFL)s
    perm_pffls <- c(pffl_perm_mirna_vector, pffl_perm_tf_vector, pffl_perm_targetgene_vector, pffl_perm_all_vector)
    obs_pffl <- row["p_ffl"][[1]]
    p_val <- mean(perm_pffls > obs_pffl)
    return(p_val)
  }

  #####apply function to every row of ffls df
  set.seed(seed)
  ffls$p_val <- apply(ffls, 1, step4_ffl_pval)
  #return ffls df with added column of p-values
  return(ffls)
}
#####fin
