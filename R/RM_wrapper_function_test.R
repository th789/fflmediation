
#databases
#load(file = "/Users/than/Dropbox/research/work_draft/databases/tf_mirna_db.rda")
#load(file = "/Users/than/Dropbox/research/work_draft/databases/tf_targetgene_db.rda")
#load(file = "/Users/than/Dropbox/research/work_draft/databases/mirna_targetgene_db.rda")

#mirna, tf, targetgene lists
#load(file = "/Users/than/Dropbox/research/work_draft/databases/names_tf_db.rda")
#load(file = "/Users/than/Dropbox/research/work_draft/databases/names_mirna_db.rda")
#load(file = "/Users/than/Dropbox/research/work_draft/databases/names_targetgene_db.rda")

#expression data
load(file = "/Users/than/Dropbox/research/work_draft/expression_data_subgroups.rda")


t1114_results_mirna <- predict_ffls(mirna_expr = miRNA_t1114,
                                    mrna_expr = mRNA_t1114,
                                    ffl_type = "miRNA",
                                    mediation_alpha = 0.05,
                                    num_bootstrap_samples = 100,
                                    num_permutations = 100,
                                    permutation_test_alpha = 0.05,
                                    seed = 12345)

t1114_results_tf <- predict_ffls(mirna_expr = miRNA_t1114,
                                 mrna_expr = mRNA_t1114,
                                 ffl_type = "TF",
                                 mediation_alpha = 0.05,
                                 num_bootstrap_samples = 100,
                                 num_permutations = 100,
                                 permutation_test_alpha = 0.05,
                                 seed = 12345)
