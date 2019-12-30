#' @title predict_ffls
#' @description Predict FFLs using mediation model approach
#' @param mirna_expr Dataframe (miRNA x samples) of miRNA expression data with miRNAs in accession number format
#' @param mrna_expr Dataframe (mRNA x samples) of mRNA expression data with genes in EnsemblID format
#' @param ffl_type Character ("miRNA" or "TF") indicating the FFL type (miRNA-FFL or TF-FFL)
#' @param mediation_alpha Significance level of coefficients in the mediation model's linear equations (default is 0.05)
#' @param num_bootstrap_samples Number of bootstrap samples (default is 1000)
#' @param seed random seed for bootstrapping and permutation test (default is 12345)
#' @return !!!!!not sure yet: List of length X
#' @export

#####predict_ffls
predict_ffls <- function(mirna_expr,
                         mrna_expr,
                         ffl_type = c("miRNA", "TF"),
                         mediation_alpha = 0.05,
                         num_bootstrap_samples = 1000,
                         seed = 12345){
  #####1. generate list of candidate ffls
  print(paste0("step 1: generate list of candidate ", ffl_type, "-FFLs"))
  candidate_ffls <- step1_candidate_ffls(mirna_expr = mirna_expr,
                                         mrna_expr = mrna_expr,
                                         ffl_type = ffl_type)
  #error message
  if(dim(candidate_ffls)[1] == 0) stop(paste0("no candidate ", ffl_type, "-FFLs predicted by databases"))

  #####2. use mediation model to evaluate each candidate ffls
  print(paste0("step 2: identify candidate ", ffl_type, "-FFLs that meet mediation model conditions"))
  ffls_mediation <- step2_mediation(mirna_expr = mirna_expr,
                                    mrna_expr = mrna_expr,
                                    candidate_ffls = candidate_ffls,
                                    ffl_type = ffl_type,
                                    alpha = mediation_alpha)
  #error message
  if(dim(ffls_mediation)[1] == 0) stop(paste0("no candidate ", ffl_type, "-FFLs meet mediation model conditions"))

  #####3. for candidate ffls that meet mediation model conditions, calculate p(ffl) via bootstrapping
  print(paste0("step 3: calculate p(FFL) of the ", ffl_type, "-FFLs that meet mediation model conditions through bootstrapping"))
  ffls_mediation <- step3_pffl(mirna_expr = mirna_expr,
                               mrna_expr = mrna_expr,
                               ffls = ffls_mediation,
                               ffl_type = ffl_type,
                               num_bootstrap_samples = num_bootstrap_samples,
                               seed = seed,
                               alpha = mediation_alpha)

  #####4. for candidate ffls that meet mediation model conditions, calculate statistical significance via permutation test
  print(paste0("step 4: calculate statistical significance of the ", ffl_type, "-FFLs that meet mediation model conditions through permutation test"))

  ###!!!!!need to fix return, filler for now, need to add roxygen
  return(list("candidate_ffls" = candidate_ffls,
              "ffls_mediation" = ffls_mediation))
}
#####fin
