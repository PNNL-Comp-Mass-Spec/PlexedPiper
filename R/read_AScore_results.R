#' Read AScore results from folder
#'
#' This function read AScore results from a local folder.
#'
#' @param path_to_AScore_results (character) Path to AScore results
#' @return ascore (data.frame) AScore results
#'
#' @importFrom dplyr rename
#' @export read_AScore_results
#'
#' @examples \dontrun{
#' ascore <- read_AScore_results("./data/ascore_outputs")
#' msnid <- best_PTM_location_by_ascore(msnid, ascore)
#' }


read_AScore_results <- function(path_to_AScore_results) {
  
  ascore <- collate_files(path_to_AScore_results, "_syn_ascore.txt")
  
  ascore <- ascore %>% rename(spectrumFile = Dataset)
  
  return(ascore)
  
}
