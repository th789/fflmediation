#' @title predict_ffls
#' @description Predict FFLs using mediation model approach
#' @param mirna_expr Dataframe (miRNA x samples) of miRNA expression data with miRNAs in accession number format
#' @param mrna_expr Dataframe (mRNA x samples) of mRNA expression data with genes in EnsemblID format
#' @param ffl_type Character ("miRNA" or "TF") indicating the FFL type (miRNA-FFL or TF-FFL)
#' @param mediation_alpha Significance level of coefficients in the mediation model's linear equations (default is 0.05)
#' @param num_bootstrap_samples Number of bootstrap samples when calculating p(FFL) (default is 1000)
#' @param num_permutations Number of permutations in permutation test (default is 1000)
#' @param permutation_test_alpha Significance level of FFLs from permutation test (default is 0.05)
#' @param seed random seed for bootstrapping and permutation test (default is 12345)
#' @return List (length = 3) containing three dataframes ("ffls_candidate" is a dataframe of candidate FFLs; "ffls_mediation" is a dataframe of candidate FFLs that meet mediation model conditions; "ffls_significant" is a dataframe of candidate FFLs that meet mediation model conditions and are statistically significant)
#' @export

#####predict_ffls function
predict_ffls <- function(mirna_expr,
                         mrna_expr,
                         ffl_type = c("miRNA", "TF"),
                         mediation_alpha = 0.05,
                         num_bootstrap_samples = 1000,
                         num_permutations = 1000,
                         permutation_test_alpha = 0.05,
                         seed = 12345){
  print("***analyses beginning***")
  #####1. generate list of candidate ffls
  print(paste0("step 1: generate list of candidate ", ffl_type, "-FFLs"))
  ffls_candidate <- step1_candidate_ffls(mirna_expr = mirna_expr,
                                         mrna_expr = mrna_expr,
                                         ffl_type = ffl_type)
  #error message
  if(dim(ffls_candidate)[1] == 0) stop(paste0("no candidate ", ffl_type, "-FFLs predicted by databases"))

  #####2. use mediation model to evaluate each candidate ffls
  print(paste0("step 2: identify candidate ", ffl_type, "-FFLs that meet mediation model conditions"))
  ffls_mediation <- step2_mediation(mirna_expr = mirna_expr,
                                    mrna_expr = mrna_expr,
                                    candidate_ffls = ffls_candidate,
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
  ffls_mediation <- step4_permutation_test(mirna_expr = mirna_expr,
                                           mrna_expr = mrna_expr,
                                           ffls = ffls_mediation,
                                           ffl_type = ffl_type,
                                           num_permutations = num_permutations,
                                           num_bootstrap_samples = num_bootstrap_samples,
                                           alpha = mediation_alpha,
                                           seed = seed)
  #subset ffls_mediation that are statistically significant
  ffls_significant <- ffls_mediation[ffls_mediation$p_val < permutation_test_alpha, ]

  #####5. return final results
  print("***analyses complete***")
  return(list("ffls_candidate" = ffls_candidate,
              "ffls_mediation" = ffls_mediation,
              "ffls_significant" = ffls_significant))
}
#####fin
