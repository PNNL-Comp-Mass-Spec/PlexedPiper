#' @title Improve Phosphosite Localization
#'
#' @description Use AScore results to improve phosphosite localization.
#'
#' @param msnid (MSnID object) MS/MS ID data
#' @param ascore (data.frame) AScore results
#'
#' @return (MSnID object) MS/MS ID data with added AScore
#'
#' @importMethodsFrom MSnID $<-
#' @importFrom dplyr group_by summarise ungroup rename inner_join
#'
#' @export best_PTM_location_by_ascore
#'
#' @examples
#' \dontrun{
#' # Get AScore results
#' path_to_ascore <- system.file("extdata/phospho/ascore_output",
#'                               package = "PlexedPiperTestData")
#' ascore <- read_AScore_results(path_to_ascore)
#'
#' # MS/MS IDs
#' path_to_MSGF_results <- system.file("extdata/global/msgf_output",
#'                                     package = "PlexedPiperTestData")
#' msnid <- read_msgf_data(path_to_MSGF_results)
#'
#' # Improve phosphosite localization
#' msnid <- best_PTM_location_by_ascore(msnid, ascore)
#' }


best_PTM_location_by_ascore <- function(msnid, ascore){
   #
   ascore <- ascore %>%
      group_by(spectrumFile, Scan, OriginalSequence, BestSequence) %>%
      summarise(maxAScore = max(AScore)) %>%
      ungroup() %>%
      dplyr::rename(peptide = OriginalSequence)

   psms(msnid) <- inner_join(psms(msnid), ascore,
                             by = c("spectrumFile", "Scan", "peptide")) %>%
      dplyr::rename(OriginalPeptide = peptide,
                    peptide = BestSequence)

   return(msnid)
}

utils::globalVariables(c("spectrumFile", "Scan", "OriginalSequence",
                         "BestSequence", "AScore"))

