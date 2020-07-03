## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = ""
)

## ----setup---------------------------------------------------------------
library(fflmediation)

## ---- include = TRUE-----------------------------------------------------
#sample miRNA expression data: miRNAs x samples (in this example: 3 x 50)
sample_mirna_expr

#sample mRNA expression data: mRNAs x samples (in this example 5 x 50)
sample_mrna_expr

## ---- include = TRUE-----------------------------------------------------
#candidate_ffls_mirna_all
candidate_ffls_mirna_all <- rbind(candidate_ffls_mirna_all_part1, candidate_ffls_mirna_all_part2)
head(candidate_ffls_mirna_all)

#candidate_ffls_tf_all
head(candidate_ffls_tf_all)

#sample candidate miRNA-FFLs for analyses (a subset of the rows in candidate_ffls_mirna_all): 3 x 9
sample_candidate_ffls_mirna

## ---- include = TRUE-----------------------------------------------------
#one-group method, miRNA-FFL example
sample_output_onegroup_mirna_ffls <- 
  predict_ffls_one_group(mirna_expr = sample_mirna_expr, 
                         mrna_expr = sample_mrna_expr,
                         ffl_type = "miRNA",
                         candidate_ffls = sample_candidate_ffls_mirna,
                         num_bootstrap_samples = 100, num_permutations = 100,
                         p_value_adjust_method = "fdr", 
                         seed = 12345, 
                         parallel = FALSE)

sample_output_onegroup_mirna_ffls

## ---- include = TRUE-----------------------------------------------------
#####group 1
#sample miRNA expression data: miRNAs x samples (in this example: 3 x 50)
sample_mirna_expr_g1

#sample mRNA expression data: mRNAs x samples (in this example 5 x 50)
sample_mrna_expr_g1

#####group 2
#sample miRNA expression data: miRNAs x samples (in this example: 3 x 50)
sample_mirna_expr_g2

#sample mRNA expression data: mRNAs x samples (in this example 5 x 50)
sample_mrna_expr_g2

## ---- include = TRUE-----------------------------------------------------
#two-group method, miRNA-FFL example
sample_output_twogroups_mirna_ffls <- 
  predict_ffls_two_groups(mirna_expr_g1 = sample_mirna_expr_g1, 
                          mrna_expr_g1 = sample_mrna_expr_g1,
                          mirna_expr_g2 = sample_mirna_expr_g2, 
                          mrna_expr_g2 = sample_mrna_expr_g2,
                          ffl_type = "miRNA",
                          candidate_ffls = sample_candidate_ffls_mirna,
                          num_bootstrap_samples = 100, 
                          num_permutations = 100,
                          p_value_adjust_method = "fdr", 
                          seed = 12345,
                          parallel = FALSE)
sample_output_twogroups_mirna_ffls

