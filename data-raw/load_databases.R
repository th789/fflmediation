

#tf-mirna interactions
load(file = "/Users/than/Dropbox/research/work_draft/databases/tf_mirna_db.rda")
usethis::use_data(tf_mirna_db, compress = "xz")

#tf-targetgene interactions
load(file = "/Users/than/Dropbox/research/work_draft/databases/tf_targetgene_db.rda")
usethis::use_data(tf_targetgene_db, compress = "xz")

#mirna-targetgene interactions
load(file = "/Users/than/Dropbox/research/work_draft/databases/mirna_targetgene_db.rda")
usethis::use_data(mirna_targetgene_db, compress = "xz")


