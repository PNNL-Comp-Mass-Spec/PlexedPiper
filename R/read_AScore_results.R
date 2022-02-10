#' Read AScore results from folder
#'
#' This function reads AScore results from a local folder.
#'
#' @param path_to_AScore_results (character) path to folder containing
#'   a file ending in "_syn_ascore.txt".
#'
#' @return \code{\link[base]{data.frame}} of AScore results
#'
#' @importFrom dplyr rename
#'
#' @export read_AScore_results
#'
#' @examples
#' \dontrun{
#' # Get AScore results
#' path_to_AScore_results <- system.file("extdata/phospho/ascore_output",
#'                                       package = "PlexedPiperTestData")
#' ascore <- read_AScore_results(path_to_AScore_results)
#' }


read_AScore_results <- function(path_to_AScore_results) {

  ascore <- collate_files(path_to_AScore_results, "_syn_ascore.txt") %>%
    dplyr::rename(spectrumFile = Dataset)

  return(ascore)
}

utils::globalVariables("Dataset")

