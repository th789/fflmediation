#' @title predict_ffls
#' @description Predict FFLs using mediation model approach
#' @param mirna_expr Dataframe (miRNA x samples) of miRNA expression data with miRNAs in accession number format
#' @param mrna_expr Dataframe (mRNA x samples) of mRNA expression data with genes in EnsemblID format
#' @param mediation_alpha Significance level of coefficients in the mediation model's linear equations (default is 0.05)
#' @return !!!!!not sure yet
#' @export !!!!!not sure yet

predict_ffls <- function(mirna_expr, mrna_expr, mediation_alpha = 0.05){
  #####1. generate list of candidate ffls
  print("step 1: generate list of candidate ffls")
  candidate_ffls <- step1_candidate_ffls(mirna_expr = mirna_expr,
                                         mrna_expr = mrna_expr)
  candidate_ffls_mirna <- candidate_ffls[["mirna_ffls"]]
  candidate_ffls_tf <- candidate_ffls[["tf_ffls"]]
  #error message
  if(dim(candidate_ffls_mirna)[1] == 0 & dim(candidate_ffls_tf)[1] == 0) stop("no candidate FFLs predicted by databases")

  #####2. use mediation model to evaluate each candidate ffls
  print("step 2: evaluate candidate ffls using mediation model")
  #mirna-ffls
  ffls_mirna <- data.frame() #empty dataframe (to create ffls_mirna variable) so that error message works
  if(dim(candidate_ffls_mirna)[1] > 0){
    ffls_mirna <- step2_mediation(mirna_expr = mirna_expr,
                                mrna_expr = mrna_expr,
                                candidate_ffls = candidate_ffls_mirna,
                                ffl_type = "miRNA",
                                alpha = mediation_alpha)}
  #tf-ffls
  ffls_tf <- data.frame() #empty dataframe (to create ffls_tf variable) so that error message works
  if(dim(candidate_ffls_tf)[1] > 0){
    ffls_tf <- step2_mediation(mirna_expr = mirna_expr,
                             mrna_expr = mrna_expr,
                             candidate_ffls = candidate_ffls_tf,
                             ffl_type = "TF",
                             alpha = mediation_alpha)}
  #error message
  if(dim(ffls_mirna)[1] == 0 & dim(ffls_tf)[1] == 0) stop("no candidate FFLs meet mediation model requirements")

  #####3. for candidate ffls that meet mediation model requirements, calculate p(ffl) via bootstrapping
  print("step 3: calculate p(FFL) of each ffl through bootstrapping")
  #!!! if df has at least one row

  #####4. for candidate ffls that meet mediation model requirements, calculate statistical significance via permutation test
  print("step 4: calculate statistical significance of each ffl through permutation test")

  ###!!!!!need to fix return, filler for now, need to add roxygen
  return(list("candidate_ffls_mirna" = candidate_ffls_mirna,
              "candidate_ffls_tf" = candidate_ffls_tf,
              "ffls_mirna" = ffls_mirna,
              "ffls_tf" = ffls_tf))
}
