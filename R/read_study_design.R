#' Utilities for reading study design files.
#'
#' @description
#' Fetches study design results from local folder.
#' Checks that study design files are internally consistent.
#'
#' @description
#' * `read_study_design()`: returns a list of study design tables, accessible by $
#'
#' @param path_to_folder (string) path to folder containing study design files
#' @param dataPkgNumber (integer) data package number for DMS
#'
#' @importFrom readr read_tsv
#' @importFrom dplyr filter select rename %>%
#' @name read_study_design
#'
#' @examples
#' study_design <- read_study_design("data/study_design_folder")
#' 
#' fractions  <- study_design$fractions
#' samples    <- study_design$samples
#' references <- study_design$references
#' 


#' @export
#' @rdname read_study_design
# gets 3 study design files from local directory
read_study_design <- function(path_to_study_design) {
  
  
  ## fetch samples.txt
  pathToFile <- list.files(path=path_to_study_design,
                           pattern="^samples.txt$",
                           full.names=T)
  if(length(pathToFile) == 0) {
    stop("'samples.txt' not found.")
  }
  
  samples <- readr::read_tsv(pathToFile,
                             col_types=readr::cols(.default = "c"),
                             progress=FALSE)
  required_samples_columns <- c("PlexID",
                                "ReporterName",
                                "ReporterAlias",
                                "MeasurementName")
  if (!all(required_samples_columns %in% colnames(samples))) {
    message("\nRequired column(s) not found in the 'samples.txt' file: ",
            paste(required_samples_columns[!(required_samples_columns %in% colnames(samples))], collapse = ", "))
    stop("Incorrect column names or missing columns in the 'samples' study design table.")
  }
  
  ## fetch fractions.txt
  pathToFile <- list.files(path=path_to_study_design,
                           pattern="^fractions.txt$",
                           full.names=T)
  if (length(pathToFile) == 0){
    stop("'fractions.txt' not found.")
  }
  
  fractions <- readr::read_tsv(pathToFile,
                               col_types=readr::cols(.default = "c"),
                               progress=FALSE)
  required_fractions_columns <- c("PlexID",
                                  "Dataset")
  if (!all(required_fractions_columns %in% colnames(fractions))) {
    message("\nRequired column(s) not found in the 'fractions.txt' file: ", 
            paste(required_fractions_columns[!(required_fractions_columns %in% colnames(fractions))], collapse = ", "))
    stop("Incorrect column names or missing columns in the 'fractions' table.")
  }
  
  ## fetch references.txt
  pathToFile <- list.files(path=path_to_study_design,
                           pattern="^references.txt$",
                           full.names=T)
  if(length(pathToFile) > 0) {
    references <- readr::read_tsv(pathToFile,
                                  col_types=readr::cols(.default = "c"),
                                  progress=FALSE)
    required_references_columns <- c("PlexID",
                                     "Reference")
    if (!all(required_references_columns %in% colnames(references))) {
      message("\nRequired column(s) not found in the 'references.txt' file: ",
              paste(required_references_columns[!(required_references_columns %in% colnames(references))], collapse = ", "))
      stop("Incorrect column names or missing columns in the 'references' table.")
    }
  } else {
    warning("'references.txt' not found. It will be made automatically from `samples.txt`")
    references <- filter(samples, is.na(MeasurementName)) %>%
      select(PlexID, QuantBlock, ReporterAlias) %>%
      rename(Reference = ReporterAlias)
  }
  
  # Check for duplicates
  if (any(duplicated(fractions$Dataset))) {
    stop("Duplicate datasets in 'fractions.txt'")
  }
  if (any(duplicated(samples$MeasurementName[!is.na(samples$MeasurementName)]))) {
    stop("Duplicate sample names in 'samples.txt'")
  }
  if (!setequal(fractions$PlexID, samples$PlexID)) {
    stop("Plex IDs in 'fractions.txt' and 'samples.txt' do not match.")
  }
  
  study_design <- list(samples = samples,
                    fractions = fractions,
                    references = references)
  
  return(study_design)
}
