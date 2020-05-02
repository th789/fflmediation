

#sample code for miRNA-FFLs (analogous for TF-FFLs)



# one-group analysis ------------------------------------------------------

#####step1. create a list of candidate miRNA-FFLs in format of candidate_ffls_mirna_all
#####can subset rows from candidate_ffls_mirna_all, or create your own dataframe
#sample_candidate_ffls_mirna

#####step2. gather miRNA and mRNA expression data in the format of sample_mirna_expr and sample_mrna_expr, respectively
#sample_mirna_expr
#sample_mrna_expr

#####step3. perform one-group analysis
# sample_output_onegroup_mirna_ffls <- predict_ffls_one_group(mirna_expr = sample_mirna_expr, mrna_expr = sample_mrna_expr,
#                                                              ffl_type = "miRNA",
#                                                              candidate_ffls = sample_candidate_ffls_mirna,
#                                                              first_row = 1, last_row = nrow(sample_candidate_ffls_mirna),
#                                                              num_bootstrap_samples = 100, num_permutations = 100,
#                                                              p_value_adjust_method = "fdr", seed = 12345)

#save(sample_output_onegroup_mirna_ffls, file = "~/Dropbox/research/fflmediation/data/sample_output_onegroup_mirna_ffls.rda")


# two-group analysis ------------------------------------------------------


#####step1. create a list of candidate miRNA-FFLs in format of candidate_ffls_mirna_all
#####can subset rows from candidate_ffls_mirna_all, or create your own dataframe
#sample_candidate_ffls_mirna

#####step2. gather miRNA and mRNA expression data in the format of sample_mirna_expr and sample_mrna_expr, respectively, for both groups
###group1
#sample_mirna_expr_g1
#sample_mrna_expr_g1
###group2
#sample_mirna_expr_g2
#sample_mrna_expr_g2

#####step3. perform two-group analysis
# sample_output_twogroups_mirna_ffls <- predict_ffls_two_groups(mirna_expr_g1 = sample_mirna_expr_g1, mrna_expr_g1 = sample_mrna_expr_g1,
#                                                               mirna_expr_g2 = sample_mirna_expr_g2, mrna_expr_g2 = sample_mrna_expr_g2,
#                                                               ffl_type = "miRNA",
#                                                               candidate_ffls = sample_candidate_ffls_mirna,
#                                                               first_row = 1, last_row = nrow(sample_candidate_ffls_mirna),
#                                                               num_bootstrap_samples = 100, num_permutations = 100,
#                                                               p_value_adjust_method = "fdr", seed = 12345)

#save(sample_output_twogroups_mirna_ffls, file = "~/Dropbox/research/fflmediation/data/sample_output_twogroups_mirna_ffls.rda")

