#' Makes MSnSet from Quantitative Crosstab.
#'
#' Creates MSnSet object from quantitative crosstab and attached phenotype data.
#'
#' @param crosstab (matrix) with log2-transformed relative reporter ion intensities.
#' Row names are the names of the measured species.
#' Column names are the names of the samples.
#' @param samples (table) matches sample names to reporter ions
#' 
#' @return (MSnSet object) with log2-transformed relative reporter ion intensities as expression
#' data and the samples table as phenotype data.
#' 
#' @importFrom MSnbase MSnSet pData
#' @export create_msnset
#' @examples \dontrun{
#' #' # Prepare MS/MS IDs
#' path_to_MSGF_results <- system.file("extdata/global/msgf_output", 
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
#' show(msnid)package = "PlexedPiperTestData")
#' msnid <- read_msgf_data(path_to_MSGF_results)
#' msnid <- MSnID::correct_peak_selection(msnid)
#' show(msnid)
#' msnid <- filter_msgf_data_peptide_level(msnid, 0.01)
#' show(msnid)
#' path_to_FASTA <- system.file("extdata/Rattus_norvegicus_NCBI_RefSeq_2018-04-10.fasta.gz", 
#'                             package = "PlexedPiperTestData")
#' msnid <- compute_num_peptides_per_1000aa(msnid, path_to_FASTA)
#' msnid <- filter_msgf_data_protein_level(msnid, 0.01)
#' show(msnid)
#'
#' # Prepare table with reporter ion intensities
#' path_to_MASIC_results <- system.file("extdata/global/masic_output", 
#'                             package = "PlexedPiperTestData")
#' masic_data <- read_masic_data(path_to_MASIC_results, interference_score=TRUE)
#' masic_data <- filter_masic_data(masic_data, 0.5, 0)
#'
#' # Creating cross-tab
#' library(readr)
#' fractions <- read_tsv(system.file("extdata/study_design/fractions.txt", 
#'                             package = "PlexedPiperTestData"))
#' samples <- read_tsv(system.file("extdata/study_design/samples.txt", 
#'                             package = "PlexedPiperTestData"))
#' references <- read_tsv(system.file("extdata/study_design/references.txt", 
#'                             package = "PlexedPiperTestData"))
#' aggregation_level <- c("accession")
#'
#' crosstab <- create_crosstab(msnid, masic_data, aggregation_level, 
#'                             fractions, samples, references)
#' m <- create_msnset(crosstab, samples)
#' dim(m)
#' head(exprs(m))
#' }

create_msnset <- function(crosstab, samples) {
   
   m <- MSnSet(as.matrix(crosstab))
   
   p <- as.data.frame(samples)
   p <- p[!is.na(p$MeasurementName), ]
   rownames(p) <- p$MeasurementName
   p <- p[sampleNames(m),]

   pData(m) <- p
   
   return(m)
}
