
###TF-miRNA interactions
#' @title Transcription Factor (TF)-miRNA Interactions
#' @description A DataFrame of TF-miRNA interactions
#' compiled from databases (TransmiR v2.0).
#' @format A DataFrame (117543 x 3) with the following columns:
#' \describe{
#' \item{tf}{TF involved in the interaction [Ensembl ID]}
#' \item{mirna}{miRNA involved in the interaction (miRNA targeted by TF) [Accession Number]}
#' \item{TRANSMIR}{indicator of whether interaction is in TransmiR v2.0 database}
#' }
"tf_mirna_db"


###TF-target gene interactions
#' @title Transcription Factor (TF)-Target Gene Interactions
#' @description A DataFrame of TF-target gene interactions
#' compiled from databases (ENCODE and TRRUST v2.0).
#' @format A DataFrame (1624664 x 4) with the following columns:
#' \describe{
#' \item{tf}{TF involved in the interaction [Ensembl ID]}
#' \item{targetgene}{target gene involved in the interaction (gene targeted by TF) [Ensembl ID]}
#' \item{TRRUST}{indicator of whether interaction is in TRRUST v2.0 database}
#' \item{ENCODE}{indicator of whether interaction is in ENCODE database}
#' }
"tf_targetgene_db"


###miRNA-target gene interactions
#' @title miRNA-Target Gene Interactions
#' @description A DataFrame of miRNA-target gene interactions
#' compiled from databases (TargetScan v7.2, miRTarBase v7.0,
#' miRDB v6.0, and miRanda v3.0).
#' @format A DataFrame (2477728 x 6) with the following columns:
#' \describe{
#' \item{mirna}{miRNA involved in the interaction [Accession Number]}
#' \item{targetgene}{target gene involved in the interaction (gene targeted by miRNA) [Ensembl ID]}
#' \item{TARGETSCAN}{indicator of whether interaction is in TargetScan v7.2 database}
#' \item{MIRTARBASE}{indicator of whether interaction is in miRTarBase v7.0 database}
#' \item{MIRDB}{indicator of whether interaction is in miRDB v6.0 database}
#' \item{MIRANDA}{indicator of whether interaction is in miRanda v3.0 database}
#' }
"mirna_targetgene_db"
