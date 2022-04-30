#' Filtering MSGF Data
#'
#' Filtering MSGF data by number of missed cleavages.
#'
#' @param msnid (MSnID object) collated MSGF output
#' @param missed_cleavages_threshold (numeric) Maximum acceptable number of
#'   missed cleavages. Default is 2.
#'
#' @return (MSnID object) filtered MSGF output
#'
#' @seealso
#' \code{\link[MSnID]{assess_missed_cleavages}}
#'
#' @importFrom MSnID MSnIDFilter apply_filter assess_missed_cleavages
#'
#' @export filter_msgf_data_missed_cleavages
#'
#' @examples
#' \dontrun{
#' path_to_MSGF_results <- system.file("extdata/global/msgf_output",
#'                                     package = "PlexedPiperTestData")
#' msnid <- read_msgf_data(path_to_MSGF_results)
#' msnid <- MSnID::correct_peak_selection(msnid)
#' show(msnid)
#' # At most 2 missed cleavages
#' msnid <-  filter_msgf_data_missed_cleavages(msnid, 2)
#' show(msnid)
#' }


filter_msgf_data_missed_cleavages <- function(msnid,
                                              missed_cleavages_threshold = 2)
{
  msnid <- assess_missed_cleavages(msnid)

  filtObj <- MSnIDFilter(msnid)
  filtObj$numMissCleavages <- list(comparison="<=",
                                   threshold=missed_cleavages_threshold)

  apply_filter(msnid, filtObj)
}

