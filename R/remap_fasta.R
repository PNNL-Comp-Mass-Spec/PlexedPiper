#' Remap Sequence IDs in FASTA File
#'
#' @description These functions remap the IDs in the FASTA file from RefSeq IDs
#'   or UniProt accessions to Gene symbols. If a single gene matches more than
#'   one protein sequence, only the longest sequence is retained.
#'
#' @md
#'
#' @param path_to_FASTA (character) path to FASTA file.
#' @param organism_name (character) scientific name of organism (e.g. `"Homo
#'   sapiens"`, `"Rattus norvegicus"`, `"Mus musculus"`, etc.). Not required if
#'   `conversion_table` is supplied.
#' @param conversion_table (data.frame) optional data frame with two columns:
#'   REFSEQ or UNIPROT and SYMBOL (in that order). If not provided, it will be
#'   created with \code{\link[MSnID]{fetch_conversion_table}}.
#'
#' @return Path to remapped FASTA file.
#'
#' @importMethodsFrom MSnID $<-
#' @importFrom MSnID fetch_conversion_table remap_fasta_entry_names
#' @importFrom Biostrings readAAStringSet
#'
#' @name remap_accessions_fasta
#'
#' @examples
#' \dontrun{
#' path_to_FASTA <- system.file(
#'   "extdata/Rattus_norvegicus_NCBI_RefSeq_2018-04-10.fasta.gz",
#'   package = "PlexedPiperTestData"
#' )
#' temp_work_dir <- tempdir() # can be set to "." or getwd(), if done carefully
#' file.copy(path_to_FASTA, temp_work_dir)
#' path_to_FASTA <- file.path(temp_work_dir, basename(path_to_FASTA))
#' library(Biostrings)
#' readAAStringSet(path_to_FASTA) # RefSeq IDs
#' path_to_new_FASTA <- remap_accessions_refseq_to_gene_fasta(
#'   path_to_FASTA, organism_name = "Rattus norvegicus"
#' )
#' readAAStringSet(path_to_new_FASTA) # gene symbols
#' }


#' @export
#' @rdname remap_accessions_fasta
remap_accessions_refseq_to_gene_fasta <- function(path_to_FASTA,
                                                  organism_name,
                                                  conversion_table)
{
  if(missing(conversion_table)){
    stopifnot(!missing(organism_name))
    message("Fetching REFSEQ to SYMBOL conversion table")
    suppressWarnings(
      conversion_table <- MSnID::fetch_conversion_table(
        organism_name, from = "REFSEQ", to = "SYMBOL"
      )
    )
  }
  # Names must be the same and in this order
  stopifnot(identical(names(conversion_table), c("REFSEQ", "SYMBOL")))

  # Create new FASTA file and output file path
  path_to_FASTA_remapped <- remap_fasta_entry_names(
    path_to_FASTA,
    conversion_table = conversion_table,
    extraction_pttrn = "^([A-Z]P_\\d+)"
  )

  return(path_to_FASTA_remapped)
}




#' @export
#' @rdname remap_accessions_fasta
# This uses the gene symbols from the FASTA entry names instead of
# creating a conversion table with fetch_conversion_table
remap_accessions_uniprot_to_gene_fasta <- function(path_to_FASTA)
{
  # UniProt accession extraction pattern
  extraction_pttrn <- "\\|([^|-]+)(-\\d+)?\\|"

  mySequences <- readAAStringSet(path_to_FASTA)

  # Subset to only sequences with known genes
  mySequences <- mySequences[grepl("GN=", names(mySequences))]

  # Create conversion table from FASTA entry names
  conversion_table <- data.frame(
    UNIPROT = gsub(paste0(".*", extraction_pttrn, ".*"),
                   "\\1", names(mySequences)),
    SYMBOL = gsub(".*GN=([^ ]*).*", "\\1", names(mySequences))
  )

  # Create new FASTA file and output file path
  path_to_FASTA_remapped <- remap_fasta_entry_names(
    path_to_FASTA,
    conversion_table = conversion_table,
    extraction_pttrn = extraction_pttrn
  )

  return(path_to_FASTA_remapped)
}

