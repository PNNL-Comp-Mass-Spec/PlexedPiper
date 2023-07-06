#' Utilities for reading study design files.
#'
#' @description Fetches study design results from local folder. Checks that
#'   study design files are internally consistent.
#'
#' @md
#'
#' @param path_to_study_design (character) path to folder containing study
#'   design files.
#' @param prefix (character) optional prefix for study design tables.
#'
#' @return A list of study design tables, accessible by `$` or `[]`.
#'
#' @importFrom readr read_tsv
#' @importFrom dplyr filter select rename %>%
#' @importFrom data.table data.table
#'
#' @export read_study_design
#'
#' @examples
#' path_to_study_design <- system.file("extdata/study_design",
#'                                     package = "PlexedPiperTestData")
#' study_design <- read_study_design(path_to_study_design)
#'
#' fractions  <- study_design$fractions
#' samples    <- study_design$samples
#' references <- study_design$references
#'


# gets 3 study design files from local directory
read_study_design <- function(path_to_study_design,
                              prefix = character(0))
{
  # Columns that must be present in their respective tables
  required_cols <- list(
    "samples" = c("PlexID", "ReporterName", "ReporterAlias", "MeasurementName"),
    "fractions" = c("PlexID", "Dataset"),
    "references" = c("PlexID", "Reference")
  )

  study_design <- lapply(names(required_cols), function(name_i) {
    pttrn <- paste0(prefix, name_i, ".txt")
    pathToFile <- list.files(path = path_to_study_design,
                             pattern = paste0("^", pttrn, "$"),
                             full.names = TRUE)
    if(length(pathToFile) == 0) {
      if (name_i == "references") {
        warning(sprintf("'%s' not found. It will be made automatically from `%s`",
                        pttrn, sub("references", "samples", pttrn)))
        references <- filter(samples, is.na(MeasurementName)) %>%
          select(PlexID, QuantBlock, ReporterAlias) %>%
          rename(Reference = ReporterAlias)
      } else {
        stop(sprintf("'%s' not found.", pttrn))
      }
    } else if (length(pathToFile) > 1) {
      stop(sprintf("Multiple %s files discovered. Check prefix."))
    }

    tbl <- readr::read_tsv(pathToFile,
                           col_types = readr::cols(.default = "c"),
                           progress = FALSE)

    required_cols_i <- required_cols[[name_i]]
    col_in_tbl <- required_cols_i %in% colnames(tbl)

    if (!all(col_in_tbl)) {
      message(sprintf("\nRequired column(s) not found in the '%s' file: ",
                      pttrn),
              paste(required_cols_i[!col_in_tbl], collapse = ", "))
      stop(sprintf(
        "Incorrect column names or missing columns in the '%s' study design table.",
        name_i))
    }

    return(tbl)
  })

  names(study_design) <- names(required_cols)

  # Check for duplicates
  if (any(duplicated(study_design$fractions$Dataset))) {
    stop(sprintf("Duplicate Dataset entries in '%sfractions.txt'", prefix))
  }
  if (any(duplicated(study_design$samples$MeasurementName,
                     incomparables = NA))) {
    stop(sprintf("Duplicate MeasurementName entries in '%ssamples.txt'",
                 prefix))
  }
  if (!setequal(study_design$fractions$PlexID,
                study_design$samples$PlexID)) {
    stop(sprintf(
      "Plex IDs in '%sfractions.txt' and '%ssamples.txt' do not match.",
      prefix, prefix))
  }

  # Check for valid TMT layout
  # Set reporter_converter to prevent "no visible binding" note
  reporter_converter <- PlexedPiper::reporter_converter
  for (i in seq_along(reporter_converter)) {
    if (setequal(study_design$samples$ReporterName,
                 reporter_converter[[i]]$ReporterName)) {
      converter <- data.table(reporter_converter[[i]])
      break
    }
  }

  if (!exists("converter")) {
    stop(paste("No reporter ion converter tables match reporter ions in",
               "samples$ReporterName"))
  }

  return(study_design)
}

utils::globalVariables(c("MeasurementName", "PlexID",
                         "QuantBlock", "ReporterAlias"))

