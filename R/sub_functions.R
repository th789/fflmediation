globalVariables(c("mirna_targetgene_db", "tf_mirna_db", "tf_targetgene_db",
                  "names_mirna_db", "names_tf_db", "names_targetgene_db",
                  "mirna", "targetgene"))

#' @import dplyr
# for left_join, transmute

#' @import tidyr
# for replace_na


# 1. generate list of candidate ffls --------------------------------------

#' @title candidate_ffls
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
