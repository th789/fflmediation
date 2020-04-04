


# candidate ffls ----------------------------------------------------------



#all candidate mirna ffls --- assembled from databases
load(file = "/Users/than/Dropbox/research/work_draft/package_data/candidate_ffls_mirna_all.rda")
usethis::use_data(candidate_ffls_mirna_all, compress = "xz")

candidate_ffls_mirna_all_part1 <- candidate_ffls_mirna_all[1:12000000, ]
candidate_ffls_mirna_all_part2 <- candidate_ffls_mirna_all[12000000:25240973, ]
usethis::use_data(candidate_ffls_mirna_all_part1, compress = "xz")
usethis::use_data(candidate_ffls_mirna_all_part2, compress = "xz")


#all candidate tf ffls --- assembled from databases
load(file = "/Users/than/Dropbox/research/work_draft/package_data/candidate_ffls_tf_all.rda")
usethis::use_data(candidate_ffls_tf_all, compress = "xz")


#sample candidate mirna ffls --- for function demonstration
load(file = "/Users/than/Dropbox/research/work_draft/package_data/sample_candidate_ffls_mirna.rda")
usethis::use_data(sample_candidate_ffls_mirna, compress = "xz")


# one-group data ----------------------------------------------------------

#sample mirna_expr
load(file = "/Users/than/Dropbox/research/work_draft/package_data/sample_mirna_expr.rda")
usethis::use_data(sample_mirna_expr, compress = "xz")

#sample mrna_expr
load(file = "/Users/than/Dropbox/research/work_draft/package_data/sample_mrna_expr.rda")
usethis::use_data(sample_mrna_expr, compress = "xz")

#sample one-group output
load(file = "/Users/than/Dropbox/research/work_draft/package_data/sample_output_onegroup_mirna_ffls.rda")
usethis::use_data(sample_output_onegroup_mirna_ffls, compress = "xz")


# two-group data ----------------------------------------------------------

#sample mirna_expr_g1
load(file = "/Users/than/Dropbox/research/work_draft/package_data/sample_mirna_expr_g1.rda")
usethis::use_data(sample_mirna_expr_g1, compress = "xz")

#sample mrna_expr_g1
load(file = "/Users/than/Dropbox/research/work_draft/package_data/sample_mrna_expr_g1.rda")
usethis::use_data(sample_mrna_expr_g1, compress = "xz")

#sample mirna_expr_g2
load(file = "/Users/than/Dropbox/research/work_draft/package_data/sample_mirna_expr_g2.rda")
usethis::use_data(sample_mirna_expr_g2, compress = "xz")

#sample mrna_expr_g2
load(file = "/Users/than/Dropbox/research/work_draft/package_data/sample_mrna_expr_g2.rda")
usethis::use_data(sample_mrna_expr_g2, compress = "xz")

#sample mirna_expr
load(file = "/Users/than/Dropbox/research/work_draft/package_data/sample_output_twogroups_mirna_ffls.rda")
usethis::use_data(sample_output_twogroups_mirna_ffls, compress = "xz")

