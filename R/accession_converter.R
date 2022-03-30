#' Remap accessions
#'
#' Converting accessions from RefSeq to Gene. If \code{conversion_table} is not
#' supplied, the function leverages \code{\link[MSnID]{fetch_conversion_table}}.
#'
#' @param msnid (MSnID object) MS/MS ID data
#' @param organism_name (character) Scientific name of organism.
#' @param conversion_table (data.frame) Optional data.frame with two columns.
#'   One column should match the accessions from the `msnid` object (e.g.
#'   RefSeq). The other column is for the alternate annotation to map to (e.g.
#'   gene symbol).
#'
#' @return (MSnID object) MS/MS ID data with computed number of peptides per
#'   1000 aa. Added column name - "peptides_per_1000aa".
#'
#' @importMethodsFrom MSnID $<- accessions
#' @importFrom MSnID fetch_conversion_table
#' @importFrom dplyr bind_cols inner_join
#'
#' @name remap_accessions
#'
#' @examples
#' \dontrun{
#' path_to_MSGF_results <- system.file("extdata/global/msgf_output",
#'                                     package = "PlexedPiperTestData")
#' msnid <- read_msgf_data(path_to_MSGF_results)
#' show(msnid)
#' msnid <- remap_accessions_refseq_to_gene(msnid,
#'                                          organism_name="Rattus norvegicus")
#' show(msnid)
#' }

#' @export
#' @rdname remap_accessions
remap_accessions_refseq_to_gene <- function(msnid,
                                            organism_name,
                                            conversion_table)
{
   remap_accessions_to_gene(msnid, organism_name,
                            conversion_table, from = "REFSEQ")
}


#' @export
#' @rdname remap_accessions
remap_accessions_uniprot_to_gene <- function(msnid,
                                             organism_name,
                                             conversion_table)
{
   remap_accessions_to_gene(msnid, organism_name,
                            conversion_table, from = "UNIPROT")
}



# Generic remap function - not exported
remap_accessions_to_gene <- function(msnid, organism_name, conversion_table,
                                     from = c("REFSEQ", "UNIPROT"))
{
   # Check input
   from <- match.arg(from, choices = c("REFSEQ", "UNIPROT"))

   pttrn <- ifelse(from == "UNIPROT",
                   "^.*\\|(.+?)(-\\d+)?\\|.*", # UniProt pattern
                   ".*?(.P_\\d+)(\\.\\d+)?") # RefSeq pattern

   acc_full <- accessions(msnid)
   # this drops ^XXX and isoform number if present
   acc <- sub(pttrn, "\\1", acc_full)

   acc.x <- bind_cols(accessions = acc_full, from = acc)
   # Rename "from" to "REFSEQ" or "UNIPROT"
   names(acc.x)[names(acc.x) == "from"] <- from

   if (!missing(conversion_table)) {
      stopifnot(identical(names(conversion_table), c(from, "SYMBOL")))
      conversion_table <- filter(conversion_table,
                                 !is.na(conversion_table$SYMBOL))
      conv_map <- inner_join(acc.x, conversion_table, by = from)
   } else {
      message(sprintf("Fetching %s to SYMBOL conversion table", from))
      suppressWarnings(
         conv_ann <- fetch_conversion_table(organism_name,
                                            from = from, to = "SYMBOL")
      )
      conv_map <- inner_join(acc.x, conv_ann, by = from)
   }

   conv_vec <- conv_map$SYMBOL
   names(conv_vec) <- conv_map$accessions

   msnid$accession <- conv_vec[msnid$accession]

   # make sure decoy accessions start with XXX
   msnid$accession <- ifelse(msnid$isDecoy,
                             paste0("XXX_", msnid$accession),
                             msnid$accession)

   return(msnid)
}

utils::globalVariables("SYMBOL")

