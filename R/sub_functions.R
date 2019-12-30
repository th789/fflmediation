globalVariables(c("mirna_targetgene_db", "tf_mirna_db", "tf_targetgene_db",
                  "names_mirna_db", "names_tf_db", "names_targetgene_db",
                  "mirna", "targetgene",
                  "lm"))
#step1
#' @import dplyr
#for left_join, transmute
#' @import tidyr
#for replace_na

#step2
#' @importFrom stats lm


# 1. generate list of candidate ffls --------------------------------------

#' @title step1_candidate_ffls
#' @description Generate candidate miRNA-FFLs and TF-FFLs
#' @param mirna_expr Dataframe (miRNA x samples) of miRNA expression data with miRNAs in accession number format
#' @param mrna_expr Dataframe (mRNA x samples) of mRNA expression data with genes in EnsemblID format
#' @return Vector of two dataframes (miRNA-FFLs and TF-FFLs)

step1_candidate_ffls <- function(mirna_expr, mrna_expr){
  #####extract mirnas, tfs, and targetgenes from expression data
  names_mirna <- rownames(mirna_expr)[rownames(mirna_expr) %in% names_mirna_db]
  names_tf <- rownames(mrna_expr)[rownames(mrna_expr) %in% names_tf_db]
  names_targetgene <- rownames(mrna_expr)[rownames(mrna_expr) %in% names_targetgene_db]

  #####candidate mirna-ffls
  #mirna-tf arm
  ffls_mirna <- mirna_targetgene_db[mirna_targetgene_db$mirna %in% names_mirna & mirna_targetgene_db$targetgene %in% names_tf, ]
  names(ffls_mirna)[names(ffls_mirna) == "targetgene"] <- "tf" #rename "targetgene" col to "tf" since all targetgenes are tfs (based on the line above)
  #tf-targetgene arm
  #!!!!!subset of data -- need to efficiently merge both dfs (ffls_mirna & tf_targetgene_db)
  ffls_mirna <- merge(ffls_mirna, tf_targetgene_db[1:1000, ], by = "tf")
  ffls_mirna <- ffls_mirna[ffls_mirna$targetgene %in% names_targetgene, ] #keep only triplets where targetgene is in the expr. data
  #mirna-targetgene arm
  ffls_mirna <- ffls_mirna %>%
    left_join(mirna_targetgene_db %>% transmute(mirna, targetgene, closed_loop = "yes")) %>%
    replace_na(list(closed_loop = "no")) #check whether mirna-targetgene pair of triplet is in mirna_targetgene_db
  #keep candidate ffls
  ffls_mirna <- ffls_mirna[ffls_mirna$closed_loop == "yes", ]
  #re-order columns, drop "closed_loop" column
  ffls_mirna <- ffls_mirna[c("mirna", "tf", "targetgene", "TARGETSCAN", "MIRTARBASE", "MIRDB", "MIRANDA", "TRRUST", "ENCODE")]

  #####candidate tf-ffls
  #tf-mirna arm
  ffls_tf <- tf_mirna_db[tf_mirna_db$tf %in% names_tf & tf_mirna_db$mirna %in% names_mirna, ]
  #tf-targetgene arm
  #!!!!!subset of data -- need to efficiently merge both dfs (ffls_tf & tf_targetgene_db)
  ffls_tf <- merge(ffls_tf, tf_targetgene_db[7000:10000, ], by = "tf")
  ffls_tf <- ffls_tf[ffls_tf$targetgene %in% names_targetgene, ] #keep only triplets where targetgene is in the expr. data
  #mirna-targetgene arm
  ffls_tf <- ffls_tf %>%
    left_join(mirna_targetgene_db %>% transmute(mirna, targetgene, closed_loop = "yes")) %>%
    replace_na(list(closed_loop = "no")) #check whether mirna-targetgene pair of triplet is in mirna_targetgene_db
  #keep candidate ffls
  ffls_tf <- ffls_tf[ffls_tf$closed_loop == "yes", ]
  #re-order columns, drop "closed_loop" column
  ffls_tf <- ffls_tf[c("tf", "mirna", "targetgene", "TRANSMIR", "TRRUST", "ENCODE")]

  #####return candidate ffls
  print(paste0(dim(ffls_mirna)[1], " candidate miRNA-FFLs"))
  print(paste0(dim(ffls_tf)[1], " candidate TF-FFLs"))
  return(c("mirna_ffls" = ffls_mirna, "tf_ffls" = ffls_tf))
}



# 2. mediation model ------------------------------------------------------

#' @title mediation_ffl
#' @description Determines whether a set of miRNA, TF, and target gene expression data meet the requirements of the mediation model
#' @param mirna Vector of miRNA expression data
#' @param tf Vector of mRNA expression data
#' @param targetgene Vector target gene expression data
#' @param ffl_type Character ("miRNA" or "TF") indicating the FFL type (miRNA-FFL or TF-FFL)
#' @param alpha Significance level of coefficients in the mediation model's linear equations (default is 0.05)
#' @return Boolean indicating whether the set of miRNA, TF, and target gene expression data meet the requirements of the mediation model

#####mediation model function: does a set of mirna, tf, targetgene data meet mediation model requirements?
mediation_ffl <- function(mirna, tf, targetgene, ffl_type = c("miRNA", "TF"), alpha = 0.05){
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
#' @param alpha Significance level of coefficients in the mediation model's linear equations (default is 0.05)
#' @return Dataframe of FFLs among candidate FFLs that meet the mediation model's requirements
#####step2_mediation
step2_mediation <- function(mirna_expr, mrna_expr,
                            candidate_ffls, ffl_type = c("miRNA", "TF"),
                            alpha = 0.05){
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
  print(paste0(dim(ffls_mediation)[1], "/", dim(candidate_ffls)[1], " candidate ", ffl_type, "-FFLs meet mediation model requirements"))
  return(ffls_mediation)
}
#####fin


