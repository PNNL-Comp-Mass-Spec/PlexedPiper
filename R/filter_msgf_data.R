#' Filtering MSGF Data
#'
#' Filtering MSGF data. In this implementation the peptide level filter optimizes both ppm and
#' PepQValue thresholds to achieve maximum number of peptide identifications within
#' given FDR constrain.
#'
#' @param msnid (MSnID object) collated MSGF output
#' @param fdr.max (numeric) Maximum acceptable FDR rate. Default is 0.01.
#' @param level (character) Level at which to perform FDR filter
#' @param n.iter.grid (numeric) number of grid-distributed evaluation points. Default 500.
#' @param n.iter.nm (numeric) number of iterations for Nelder-Mead optimization algorithm. Default 100.
#' @param ... Take arguments from `filter_msgf_data`?
#'
#' @return (MSnID object) filtered MSGF output
#'
#' @importFrom MSnID MSnIDFilter MSnIDFilter optimize_filter mass_measurement_error apply_filter
#'
#' @examples \dontrun{
#' path_to_MSGF_results <- system.file("extdata/global/msgf_output", 
#'                         package = "PlexedPiperTestData")
#' msnid <- read_msgf_data(path_to_MSGF_results)
#' msnid <- MSnID::correct_peak_selection(msnid)
#' show(msnid)
#' msnid <- filter_msgf_data(msnid, "peptide", 0.01) # 1% FDR at peptide level
#' show(msnid)
#' path_to_FASTA <- system.file("extdata/Rattus_norvegicus_NCBI_RefSeq_2018-04-10.fasta.gz", 
#'                         package = "PlexedPiperTestData")
#' msnid <- compute_num_peptides_per_1000aa(msnid, path_to_FASTA)
#' msnid <- filter_msgf_data(msnid, "accession", 0.01) # 1% FDR at protein level
#' show(msnid)
#' }
#' @export
filter_msgf_data <- function(msnid,
                             level,
                             fdr.max=0.01,
                             n.iter.grid=500,
                             n.iter.nm=100){

  # Clean up on exit
  on.exit(rm(list = ls()))
  on.exit(gc(verbose = FALSE), add = TRUE)

  # Check input
  level <- match.arg(level, choices = c("peptide", "accession"))

  keep_cols <- c(level, "isDecoy") # columns to calculate FDR

  # Setup
  if (level == "peptide") {
    # Create columns for peptide filtering
    # Can not use data.table syntax if the msnid has been modified at all,
    # as it results in the "Invalid .internal.selfref" warning and
    # columns not being created.
    msnid$msmsScore <- -log10(msnid$PepQValue)
    msnid$absParentMassErrorPPM <- abs(mass_measurement_error(msnid))

    # Create MSnID of minimum size
    suppressMessages(msnid_small <- MSnID())
    # add filter criteria columns
    keep_cols <- c(keep_cols, "msmsScore", "absParentMassErrorPPM")
    msnid_small@psms <- unique(msnid@psms[, keep_cols, with = FALSE])

    # Create filter object
    filtObj <- MSnIDFilter(msnid_small)
    filtObj$absParentMassErrorPPM <- list(comparison = "<", threshold = 10)
    filtObj$msmsScore <- list(comparison = ">", threshold = 2)
    method <- "Nelder-Mead"
  } else {
    # Create MSnID of minimum size
    suppressMessages(msnid_small <- MSnID())
    keep_cols <- c(keep_cols, "peptides_per_1000aa") # add filter criteria column
    msnid_small@psms <- unique(msnid@psms[, keep_cols, with = FALSE])

    # Create filter object
    filtObj <- MSnIDFilter(msnid_small)
    filtObj$peptides_per_1000aa <- list(comparison = ">", threshold = 1)
    method <- "SANN"
  }

  # step 1
  filtObj.grid <- optimize_filter(filtObj,
                                  msnid_small,
                                  fdr.max=fdr.max,
                                  method="Grid",
                                  level=level,
                                  n.iter=n.iter.grid)
  # step 2
  filtObj.nm <- optimize_filter(filtObj.grid,
                                msnid_small,
                                fdr.max=fdr.max,
                                method=method,
                                level=level,
                                n.iter=n.iter.nm)
  return(apply_filter(msnid, filtObj.nm))
}


#' @export
#' @rdname filter_msgf_data
filter_msgf_data_peptide_level <- function(msnid, ...) {
  filter_msgf_data(msnid, level="peptide", ...)
}

#' @export
#' @rdname filter_msgf_data
filter_msgf_data_protein_level <- function(msnid, ...) {
  filter_msgf_data(msnid, level="accession", ...)
}
