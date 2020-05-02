



# candidate ffls ----------------------------------------------------------


#all candidate mirna ffls --- assembled from databases
#' @title candidate_ffls_mirna_all_part1
#' @description First half of all possible candidate miRNA FFLs assembled from databases
#' @format A DataFrame (12,000,000 x 9) with the following columns:
#' \describe{
#' \item{mirna}{miRNA in candidate FFL [Accession Number]}
#' \item{tf}{TF in candidate FFL [Ensembl ID]}
#' \item{targetgene}{target gene in candidate FFL [Ensembl ID]}
#' \item{TARGETSCAN}{indicator of whether miRNA-target gene or miRNA-TF interaction is in TargetScan v7.2 database}
#' \item{MIRTARBASE}{indicator of whether miRNA-target gene or miRNA-TF interaction is in miRTarBase v7.0 database}
#' \item{MIRDB}{indicator of whether miRNA-target gene or miRNA-TF interaction is in miRDB v6.0 database}
#' \item{MIRANDA}{indicator of whether miRNA-target gene or miRNA-TF interaction is in miRanda v3.0 database}
#' \item{TRRUST}{indicator of whether TF-target gene or miRNA-TF interaction is in TRRUST v2.0 database}
#' \item{ENCODE}{indicator of whether TF-target gene or miRNA-TF interaction is in ENCODE database}
#' }
"candidate_ffls_mirna_all_part1"


#all candidate mirna ffls --- assembled from databases
#' @title candidate_ffls_mirna_all_part2
#' @description Second half of all possible candidate miRNA FFLs assembled from databases
#' @format A DataFrame (13240974 x 9) with the following columns:
#' \describe{
#' \item{mirna}{miRNA in candidate FFL [Accession Number]}
#' \item{tf}{TF in candidate FFL [Ensembl ID]}
#' \item{targetgene}{target gene in candidate FFL [Ensembl ID]}
#' \item{TARGETSCAN}{indicator of whether miRNA-target gene or miRNA-TF interaction is in TargetScan v7.2 database}
#' \item{MIRTARBASE}{indicator of whether miRNA-target gene or miRNA-TF interaction is in miRTarBase v7.0 database}
#' \item{MIRDB}{indicator of whether miRNA-target gene or miRNA-TF interaction is in miRDB v6.0 database}
#' \item{MIRANDA}{indicator of whether miRNA-target gene or miRNA-TF interaction is in miRanda v3.0 database}
#' \item{TRRUST}{indicator of whether TF-target gene or miRNA-TF interaction is in TRRUST v2.0 database}
#' \item{ENCODE}{indicator of whether TF-target gene or miRNA-TF interaction is in ENCODE database}
#' }
"candidate_ffls_mirna_all_part2"

#all candidate tf ffls --- assembled from databases
#' @title candidate_ffls_tf_all
#' @description All possible candidate TF FFLs assembled from databases
#' @format A DataFrame (178,872 x 10) with the following columns:
#' \describe{
#' \item{tf}{TF in candidate FFL [Ensembl ID]}
#' \item{mirna}{miRNA in candidate FFL [Accession Number]}
#' \item{targetgene}{target gene in candidate FFL [Ensembl ID]}
#' \item{TRANSMIR}{indicator of whether miRNA-TF interaction is in TransmiR v2.0 database}
#' \item{TARGETSCAN}{indicator of whether miRNA-target gene or miRNA-TF interaction is in TargetScan v7.2 database}
#' \item{MIRTARBASE}{indicator of whether miRNA-target gene or miRNA-TF interaction is in miRTarBase v7.0 database}
#' \item{MIRDB}{indicator of whether miRNA-target gene or miRNA-TF interaction is in miRDB v6.0 database}
#' \item{MIRANDA}{indicator of whether miRNA-target gene or miRNA-TF interaction is in miRanda v3.0 database}
#' \item{TRRUST}{indicator of whether TF-target gene or miRNA-TF interaction is in TRRUST v2.0 database}
#' \item{ENCODE}{indicator of whether TF-target gene or miRNA-TF interaction is in ENCODE database}
#' }
"candidate_ffls_tf_all"


#sample candidate mirna ffls --- input into functions
#' @title sample_candidate_ffls_mirna
#' @description Sample candidate miRNA FFLs serving as input for \code{\link{predict_ffls_one_group}} or \code{\link{predict_ffls_two_groups}}
#' @format A DataFrame (3 x 9) with the following columns:
#' \describe{
#' \item{mirna}{miRNA in candidate FFL [Accession Number]}
#' \item{tf}{TF in candidate FFL [Ensembl ID]}
#' \item{targetgene}{target gene in candidate FFL [Ensembl ID]}
#' \item{TARGETSCAN}{indicator of whether miRNA-target gene or miRNA-TF interaction is in TargetScan v7.2 database}
#' \item{MIRTARBASE}{indicator of whether miRNA-target gene or miRNA-TF interaction is in miRTarBase v7.0 database}
#' \item{MIRDB}{indicator of whether miRNA-target gene or miRNA-TF interaction is in miRDB v6.0 database}
#' \item{MIRANDA}{indicator of whether miRNA-target gene or miRNA-TF interaction is in miRanda v3.0 database}
#' \item{TRRUST}{indicator of whether TF-target gene or miRNA-TF interaction is in TRRUST v2.0 database}
#' \item{ENCODE}{indicator of whether TF-target gene or miRNA-TF interaction is in ENCODE database}
#' }
"sample_candidate_ffls_mirna"



# one-group sample data ---------------------------------------------------


#sample mirna_expr
#' @title sample_mirna_expr
#' @description Sample mirna_expr serving as input for \code{\link{predict_ffls_one_group}}
#' @format A DataFrame (3 x 50) representing expression data of 3 miRNAs for 50 samples
"sample_mirna_expr"


#sample mrna_expr
#' @title sample_mrna_expr
#' @description Sample mrna_expr serving as input for \code{\link{predict_ffls_one_group}}
#' @format A DataFrame (5 x 50) representing expression data of 5 mRNAs for 50 samples
"sample_mrna_expr"


#sample output for predict_ffls_one_group for miRNA-FFLs
#' @title sample_output_onegroup_mirna_ffls
#' @description Sample output of \code{\link{predict_ffls_one_group}}
#' @format A DataFrame (3 x 43) representing expression data of 5 mRNAs for 50 samples
#' \describe{
#' \item{mirna_name}{miRNA in candidate FFL [miRNA name]}
#' \item{tf_name}{TF in candidate FFL [SYMBOL]}
#' \item{targetgene_name}{target gene in candidate FFL [SYMBOL]}
#' \item{pFFL}{pFFL of candidate FFL}
#' \item{p_value_adj}{adjusted p-value of candidate FFL}
#' \item{p_value}{unadjusted p-value of candidate FFL}
#' \item{mirna}{miRNA in candidate FFL [Accession Number]}
#' \item{tf}{TF in candidate FFL [Ensembl ID]}
#' \item{targetgene}{target gene in candidate FFL [Ensembl ID]}
#' \item{TARGETSCAN}{indicator of whether miRNA-target gene or miRNA-TF interaction is in TargetScan v7.2 database}
#' \item{MIRTARBASE}{indicator of whether miRNA-target gene or miRNA-TF interaction is in miRTarBase v7.0 database}
#' \item{MIRDB}{indicator of whether miRNA-target gene or miRNA-TF interaction is in miRDB v6.0 database}
#' \item{MIRANDA}{indicator of whether miRNA-target gene or miRNA-TF interaction is in miRanda v3.0 database}
#' \item{TRRUST}{indicator of whether TF-target gene or miRNA-TF interaction is in TRRUST v2.0 database}
#' \item{ENCODE}{indicator of whether TF-target gene or miRNA-TF interaction is in ENCODE database}
#' \item{remaining columns}{estimate, standard error (se), t-statistic (t), and p-value of regression coefficients}
#' }
"sample_output_onegroup_mirna_ffls"




# two-group sample data ---------------------------------------------------


#sample mirna_expr_g1
#' @title sample_mirna_expr_g1
#' @description Sample mirna_expr_g1 serving as input for \code{\link{predict_ffls_two_groups}}
#' @format A DataFrame (3 x 50) representing expression data of 3 miRNAs for 50 samples in Group 1
"sample_mirna_expr_g1"

#sample mrna_expr_g1
#' @title sample_mrna_expr_g1
#' @description Sample mrna_expr_g1 serving as input for \code{\link{predict_ffls_two_groups}}
#' @format A DataFrame (5 x 50) representing expression data of 5 mRNAs for 50 samples in Group 1
"sample_mrna_expr_g1"

#sample mirna_expr_g2
#' @title sample_mirna_expr_g2
#' @description Sample mirna_expr_g2 serving as input for \code{\link{predict_ffls_two_groups}}
#' @format A DataFrame (3 x 50) representing expression data of 3 miRNAs for 50 samples in Group 2
"sample_mirna_expr_g2"

#sample mrna_expr_g2
#' @title sample_mrna_expr_g2
#' @description Sample mrna_expr_g2 serving as input for \code{\link{predict_ffls_two_groups}}
#' @format A DataFrame (5 x 50) representing expression data of 5 mRNAs for 50 samples in Group 2
"sample_mrna_expr_g2"

#sample output for predict_ffls_two_groups for miRNA-FFLs
#' @title sample_output_twogroups_mirna_ffls
#' @description Sample output of \code{\link{predict_ffls_two_groups}}
#' @format A DataFrame (3 x 45) representing expression data of 5 mRNAs for 50 samples
#' \describe{
#' \item{mirna_name}{miRNA in candidate FFL [miRNA name]}
#' \item{tf_name}{TF in candidate FFL [SYMBOL]}
#' \item{targetgene_name}{target gene in candidate FFL [SYMBOL]}
#' \item{pFFL_group1}{pFFL of candidate FFL in Group 1}
#' \item{pFFL_group2}{pFFL of candidate FFL in Group 2}
#' \item{delta_pFFL_(group1-group2)}{delta-pFFL of candidate FFL}
#' \item{p_value_adj}{adjusted p-value of candidate FFL}
#' \item{p_value}{unadjusted p-value of candidate FFL}
#' \item{mirna}{miRNA in candidate FFL [Accession Number]}
#' \item{tf}{TF in candidate FFL [Ensembl ID]}
#' \item{targetgene}{target gene in candidate FFL [Ensembl ID]}
#' \item{TARGETSCAN}{indicator of whether miRNA-target gene or miRNA-TF interaction is in TargetScan v7.2 database}
#' \item{MIRTARBASE}{indicator of whether miRNA-target gene or miRNA-TF interaction is in miRTarBase v7.0 database}
#' \item{MIRDB}{indicator of whether miRNA-target gene or miRNA-TF interaction is in miRDB v6.0 database}
#' \item{MIRANDA}{indicator of whether miRNA-target gene or miRNA-TF interaction is in miRanda v3.0 database}
#' \item{TRRUST}{indicator of whether TF-target gene or miRNA-TF interaction is in TRRUST v2.0 database}
#' \item{ENCODE}{indicator of whether TF-target gene or miRNA-TF interaction is in ENCODE database}
#' \item{remaining columns}{estimate, standard error (se), t-statistic (t), and p-value of regression coefficients for group1 and group2}
#' }
"sample_output_twogroups_mirna_ffls"

