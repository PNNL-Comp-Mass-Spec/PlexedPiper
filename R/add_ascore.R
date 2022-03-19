#' Remap Sequence IDs in FASTA File
#'
#' The function remaps the IDs in the FASTA file from RefSeq to Gene symbols.
#' In case of multiple RefSeq sequences available, only the first longest is
#' retained.
#'
#' @param msnid (MSnID object) MS/MS ID data
#' @param ascore (data.frame) AScore results
#' @return (MSnID object) MS/MS ID data with added AScore
#'
#'
#' @importMethodsFrom MSnID $<-
#' @importFrom Biostrings readAAStringSet width writeXStringSet
#' @importFrom AnnotationHub AnnotationHub query
#' @importFrom AnnotationDbi select
#' @importFrom dplyr bind_cols filter inner_join slice pull ungroup select rename
#' @importFrom tools file_path_sans_ext file_ext
#' @export best_PTM_location_by_ascore

# best_PTM_location_by_ascore_old <- function(msnid, ascore){
#    #
#    ascore <- ascore %>%
#       group_by(Job, Scan, OriginalSequence, BestSequence) %>%
#       summarise(maxAScore = max(AScore))
# 
#    x <- ascore %>%
#       rename(Peptide = OriginalSequence) %>%
#       inner_join(psms(msnid)) %>%
#       ungroup() %>%
#       select(-Peptide) %>%
#       rename(Peptide = BestSequence)
# 
#    psms(msnid) <- x
# 
#    return(msnid)
# 
# }

best_PTM_location_by_ascore <- function(msnid, ascore){
   #
   ascore <- ascore %>%
      group_by(spectrumFile, Scan, OriginalSequence, BestSequence) %>%
      summarise(maxAScore = max(AScore)) %>%
      ungroup() %>%
      rename(peptide = OriginalSequence)

   x <- ascore %>%
      inner_join(psms(msnid), ., by=c("spectrumFile", "Scan", "peptide")) %>%
      rename(OriginalPeptide = peptide,
             peptide = BestSequence)

   psms(msnid) <- x

   return(msnid)

}



