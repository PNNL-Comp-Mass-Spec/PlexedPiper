#' Compute number of peptides per 1000 aa
#'
#' Computing number of peptides per 1000 amino acids.
#'
#' @param msnid (MSnID object) MS/MS ID data
#' @param path_to_FASTA (character) Path to FASTA file.
#'
#' @return (MSnID object) MS/MS ID data with computed number of peptides per
#'   1000 aa. Added column name - "peptides_per_1000aa".
#'
#' @importFrom dplyr mutate select distinct group_by summarise n inner_join
#' @importFrom MSnID psms
#' @importFrom Biostrings readAAStringSet width
#'
#' @export compute_num_peptides_per_1000aa
#'
#' @examples
#' \dontrun{
#' path_to_MSGF_results <- system.file("extdata/global/msgf_output",
#'                                     package = "PlexedPiperTestData")
#' msnid <- read_msgf_data(path_to_MSGF_results)
#' path_to_FASTA <- system.file(
#'   "extdata/Rattus_norvegicus_NCBI_RefSeq_2018-04-10.fasta.gz",
#'   package = "PlexedPiperTestData"
#' )
#' msnid <- compute_num_peptides_per_1000aa(msnid, path_to_FASTA)
#' hist(msnid$peptides_per_1000aa)
#' }


compute_num_peptides_per_1000aa <- function(msnid,
                                            path_to_FASTA){
   # fasta
   mySequences <- readAAStringSet(path_to_FASTA)
   # compute protein lengths
   prot_lengths <-
      data.frame(accession = sub("^(\\S+)\\s.*","\\1",names(mySequences)),
                 Length = width(mySequences),
                 stringsAsFactors = FALSE)
   # append with lengths of reverse sequences
   prot_lengths <- mutate(prot_lengths, accession = paste0("XXX_", accession)) %>%
      rbind(prot_lengths, .)

   # count peptides per 1000 aa
   pepN1000 <- psms(msnid) %>%
      dplyr::select(accession, peptide, isDecoy) %>%
      distinct() %>%
      group_by(accession) %>%
      summarise(pepN = n()) %>%
      inner_join(prot_lengths, by = "accession") %>%
      mutate(pep_per_1000 = 1000*pepN/Length)
   pepN1000_vec <- pepN1000$pep_per_1000
   names(pepN1000_vec) <- pepN1000$accession

   msnid$peptides_per_1000aa <- pepN1000_vec[msnid$accession]
   return(msnid)
}

utils::globalVariables(c("isDecoy", "pepN", "Length"))


