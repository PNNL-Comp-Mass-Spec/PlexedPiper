#' Filtering msfragger Data
#'
#' Filtering msfragger data. In this implementation, the peptide-level filter
#' optimizes both ppm and one of Expectation or PeptideProphet Probability thresholds
#' to achieve maximum number of peptide identifications within a given FDR constraint.
#' The accession-level filter optimizes based on `peptides_per_1000aa`, so
#' \code{\link{compute_num_peptides_per_1000aa}} must be used first.
#'
#' @md
#'
#' @param msnid (MSnID object) collated msfragger output
#' @param fdr.max (numeric) Maximum acceptable FDR. Default is 0.01 (1%).
#' @param level (character) Level at which to perform FDR filter. The name of a
#'   column in `psms(msnid)`. Currently, only `"peptide"` or `"accession"` are
#'   supported. The added level `SiteID` makes sense only for PTM data and 
#'   first requires mapping of the modification site using
#'    `MSnID::map_mod_sites`.
#' @param filtering_criterion (character) One of "evalue" which is 
#'    expectation value or "pp_prob" - peptide prophet probability. Default is
#'    "pp_prob".
#' @param n.iter.grid (numeric) number of grid-distributed evaluation points.
#' @param n.iter.nm (numeric) number of iterations for Nelder-Mead optimization
#'   algorithm.
#' @param ... arguments passed to `filter_msfragger_data`.
#'
#' @return (MSnID object) filtered msfragger output
#'
#' @seealso
#' \code{\link[MSnID]{MSnIDFilter}}
#' \code{\link[MSnID]{optimize_filter}}
#' \code{\link[MSnID]{apply_filter}}
#'
#' @importFrom MSnID MSnIDFilter optimize_filter
#'   mass_measurement_error apply_filter


#' @export
filter_msfragger_data <- function(msnid,
                             level,
                             filtering_criterion = c("pp_prob","evalue"),
                             fdr.max=0.01,
                             n.iter.grid=500,
                             n.iter.nm=100){

  # Clean up on exit
  on.exit(rm(list = ls()))
  on.exit(gc(verbose = FALSE), add = TRUE)

  # Check input
  level <- match.arg(level, choices = c("peptide", "accession", "SiteID"))
  filtering_criterion <- match.arg(filtering_criterion)
  
  if(level == "SiteID" & !("SiteID" %in% names(msnid)))
     stop("Column 'SiteID' is not in the MSnID object. Please map the PTMs first.")

  keep_cols <- c(level, "isDecoy") # columns to calculate FDR

  # Create MSnID of minimum size
  suppressMessages(msnid_small <- MSnID())
  
  # Setup
  if (level == "accession") {

    # Add filter criteria column
    keep_cols <- c(keep_cols, "peptides_per_1000aa")
    msnid_small@psms <- unique(msnid@psms[, keep_cols, with = FALSE])

    # Create filter object
    filtObj <- MSnIDFilter(msnid_small)
    filtObj$peptides_per_1000aa <- list(comparison = ">", threshold = 1)
    method <- "SANN"
  } else {
     #Choose filter object probability value
     if (filtering_criterion == "evalue") {
        msnid$msmsScore <- -log10(msnid$Expectation)
     }
     if (filtering_criterion == "pp_prob") {
        msnid$msmsScore <- msnid$`PeptideProphet Probability`
     }
     # Create columns for peptide filtering
     # Can not use data.table syntax if the msnid has been modified at all,
     # as it results in the "Invalid .internal.selfref" warning and
     # columns not being created.
     msnid$absParentMassErrorPPM <- abs(mass_measurement_error(msnid))
     
     # Add filter criteria columns
     keep_cols <- c(keep_cols, "msmsScore", "absParentMassErrorPPM")
     msnid_small@psms <- unique(msnid@psms[, keep_cols, with = FALSE])
     
     # Create filter object
     filtObj <- MSnIDFilter(msnid_small)
     filtObj$absParentMassErrorPPM <- list(comparison = "<", threshold = 10)
     if (filtering_criterion == "evalue") {
        filtObj$msmsScore <- list(comparison = ">", threshold = 2)
     }
     if (filtering_criterion == "pp_prob") {
        filtObj$msmsScore <- list(comparison = ">", threshold = 0.99)
     }
     method <- "Nelder-Mead"
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
