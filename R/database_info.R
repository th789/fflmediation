
###miRNA-target gene databases
#' @title miRNA-Target Gene Interactions
#' @description A DataFrame of miRNA-target gene interactions
#' compiled from databases (TargetScan v7.2, miRTarBase v7.0,
#' miRDB v6.0, and miRanda v3.0).
#' @format A DataFrame (1757820 x 6) with the following columns:
#' \describe{
#' \item{mirna}{miRNA involved in the interaction}
#' \item{targetgene}{target gene involved in the interaction (gene targeted by miRNA)}
#' \item{TARGETSCAN}{indicator of whether interaction is in TargetScan v7.2 database}
#' \item{MIRTARBASE}{indicator of whether interaction is in miRTarBase v7.0 database}
#' \item{MIRDB}{indicator of whether interaction is in miRDB v6.0 database}
#' \item{MIRANDA}{indicator of whether interaction is in miRanda v3.0 database}
#' }
"mirna_targetgene_db"
