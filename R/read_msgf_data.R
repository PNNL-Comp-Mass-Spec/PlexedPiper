#' Reading MS-GF+ Results. Generic.
#'
#' Reading MS-GF+ output from a single directory.
#'
#' @param path_to_MSGF_results (character) path to directory with MSGF results
#'   for all datasets.
#' @param suffix (character) optional file suffix. Either `"_msgfplus_syn.txt"`,
#'   `"_msgfdb_syn.txt"`, or `"_syn.txt"`.
#' @param use_mzIdentML (logical) whether to read mzIdentML files into `psms`
#'   `data.table` slot of the `MSnID`.
#'
#' @md
#'
#' @return (MSnID) MSnID object
#'
#' @importFrom dplyr mutate
#' @importFrom MSnID MSnID convert_msgf_output_to_msnid
#' @importMethodsFrom MSnID psms<-
#'
#' @examples
#' \dontrun{
#' path_to_MSGF_results <- system.file("extdata/global/msgf_output",
#'                                     package = "PlexedPiperTestData")
#' msnid <- read_msgf_data(path_to_MSGF_results)
#' print(msnid)
#' head(MSnID::psms(msnid))
#' }


#' @export
read_msgf_data <- function(path_to_MSGF_results,
                           suffix = character(0),
                           use_mzIdentML = FALSE)
{
  if (use_mzIdentML) {
    mzid_files <- list.files(path_to_MSGF_results,
                             pattern = "mzid",
                             full.names = TRUE)

    if (length(mzid_files) == 0) {
      stop("mzid files not found.")
    }

    # Read mzIdentML files into psms slot
    # Note: adding a progress bar with llply drastically
    # increases memory usage and time
    suppressMessages(
      msnid <- MSnID::read_mzIDs(MSnID(), mzid_files, backend = "mzR")
    )
  } else {
    if (identical(suffix, character(0))) {
      for (pattern in c("_msgfplus_syn.txt", "_msgfdb_syn.txt", "_syn.txt")) {
        if (length(list.files(path_to_MSGF_results, pattern)) > 0) {
          suffix <- pattern
          break
        }
      }
    }

    if (identical(suffix, character(0)) |
        (length(list.files(path_to_MSGF_results, suffix)) == 0)) {
      stop("MS-GF+ results not found.")
    }

    x <- collate_files(path_to_MSGF_results, suffix)
    msnid <- convert_msgf_output_to_msnid(x)
  }
  return(msnid)
}



#' # (yet) non-exported helper function
#' #' @importFrom Biostrings AA_STANDARD
#' #' @importFrom purrr map_chr
#' #' @importFrom stringr str_replace_all
#' peptide_to_sequence <- function(pep){
#'    pep_no_flank <- sub(".\\.(.*)\\..", "\\1", pep)
#'    present_chars <- paste0(pep_no_flank, collapse = '') %>%
#'       strsplit(split='') %>%
#'       `[[`(1) %>%
#'       unique()
#'    other_chars <- setdiff(present_chars, AA_STANDARD)
#'    if(length(other_chars) > 0){
#'       message("Detected extra chararacters in the peptide sequences!")
#'       # erase other chars in TrimmedPeptide
#'       other_chars_pttrn <- other_chars %>%
#'          map_chr(~paste0("\\",.x)) %>%
#'          paste0(collapse='') %>%
#'          paste0("[",.,"]")
#'       pep_clean <- str_replace_all(pep_no_flank, other_chars_pttrn, "")
#'    }
#'    return(pep_clean)
#' }
