# Undeclared variables
utils::globalVariables(
  c(
    ".",
    "...",
    ".SD",
    ":=",
    "accession",
    "..aggregation_level",
    "d",
    "ANNOTATION",
    "AScore",
    "assess_missed_cleavages",
    "best_PTM_location_by_ascore",
    "BestSequence",
    "channel",
    "confident_score",
    "Dataset",
    "ENTREZID",
    "everything",
    "fetch_conversion_table",
    "flanking_sequence",
    "IDScan",
    "if_else",
    "ind",
    "intensity",
    "InterferenceScore",
    "isDecoy",
    "Job",
    "len",
    "Length",
    "link_min_vectorized",
    "maxAScore",
    "MeasurementName",
    "MSGFDB_SpecEValue",
    "noninferableProteins",
    "OriginalSequence",
    "pData<-",
    "pepN",
    "peptide",
    "Peptide",
    "peptide_score",
    "peptideSignature",
    "percent_coverage",
    "percentAACoverage",
    "PlexID",
    "protein_id",
    "ptm_id",
    "ptm_peptide",
    "QuantBlock",
    "quantile",
    "QuantScan",
    "read.delim",
    "redundant_ids",
    "redundantAccessions",
    "remap_accessions_uniprot_to_gene_fasta",
    "reporter_converter",
    "ReporterAlias",
    "ReporterIon",
    "s2n",
    "sampleNames",
    "Scan Filter Text",
    "Scan",
    "ScanFilterText",
    "ScanNumber",
    "sequenceWindow",
    "setTxtProgressBar",
    "SiteID",
    "Specie",
    "spectrumFile",
    "SYMBOL",
    "txtProgressBar"
  )
)

#' Remap accessions
#'
#' Converting accessions from RefSeq to Gene. If `conversion_table` is not
#' supplied, the function leverages `AnnotationHub()`
#'
#' @param msnid (MSnID object) MS/MS ID data
#' @param organism_name (string) Official organism name
#' @param conversion_table (data.frame) data frame with two columns
#' one should with named accessions and contain accessions from `msnid` object (e.g. RefSeq) the other is with alternative annotation to map to (e.g. gene symbol).
#' @return (MSnID object) MS/MS ID data with computed number of peptides per 1000 aa. Added column name - "peptides_per_1000aa".
#' @importMethodsFrom MSnID $<- accessions
#' @importFrom Biostrings readAAStringSet width
#' @importFrom AnnotationHub AnnotationHub query
#' @importFrom AnnotationDbi select
#' @importFrom dplyr bind_cols
#' @importFrom tools file_path_sans_ext file_ext
#'
#' @name remap_accessions
#'
#' @examples \dontrun{
#' path_to_MSGF_results <- system.file("extdata/global/msgf_output", package = "PlexedPiperTestData")
#' msnid <- read_msgf_data(path_to_MSGF_results)
#' show(msnid)
#' msnid <- remap_accessions_refseq_to_gene(msnid, organism_name="Rattus norvegicus")
#' show(msnid)
#' }

#' @export
#' @rdname remap_accessions
remap_accessions_refseq_to_gene <- function(msnid,
                                            organism_name,
                                            conversion_table){
   
   acc.d <- accessions(msnid)
   # by default RefSeq ID are with version digit(s)
   acc <- sub(".*?(.P_\\d+)(\\.\\d+)?", "\\1", acc.d) # this drops ^XXX
   
   acc.x <- bind_cols(accessions = acc.d, REFSEQ = acc)
   
   if(!missing(conversion_table)){
      conversion_table <- filter(conversion_table, 
                                 !is.na(conversion_table$SYMBOL))
      conv_map <- inner_join(acc.x, conversion_table, by = "REFSEQ")
   }
   
   else{
      ah <- AnnotationHub()
      orgs <- subset(ah, ah$rdataclass == "OrgDb")
      db <- query(orgs, organism_name)
      db <- db[[1]] # wierd command
      # columns(db)
      # keytypes(db)
      conv_ann <- AnnotationDbi::select(db,
                                       keys=acc,
                                       columns="SYMBOL",
                                       keytype="REFSEQ")
      conv_ann <- filter(conv_ann, !is.na(SYMBOL))
      conv_map <- inner_join(acc.x, conv_ann, by="REFSEQ")
   }

   # # add reverse IDs
   # conv_vec_rev <- paste0("XXX_",conv_vec)
   # names(conv_vec_rev) <- paste0("XXX_",names(conv_vec_rev))
   # conv_vec <- c(conv_vec, conv_vec_rev)
   
   conv_vec <- conv_map$SYMBOL
   names(conv_vec) <- conv_map$accessions
   
   msnid$accession <- conv_vec[msnid$accession]

   # make sure decoy accessions start with XXX
   msnid$accession <- ifelse(msnid$isDecoy,
                             paste0("XXX_", msnid$accession),
                             msnid$accession)

   return(msnid)
}










#' @export
#' @rdname remap_accessions
remap_accessions_uniprot_to_gene <- function(msnid,
                                            organism_name,
                                            conversion_table){
   acc_full <- accessions(msnid)
   # this drops ^XXX and isoform number if present
   acc <- sub("^.*\\|(.+?)(-\\d+)?\\|.*", "\\1", acc_full)
   
   acc.x <- bind_cols(accessions = acc_full, UNIPROT = acc)
   
   if(!missing(conversion_table)){
      conversion_table <- filter(conversion_table, 
                                 !is.na(conversion_table$SYMBOL))
      conv_map <- inner_join(acc.x, conversion_table, by = "UNIPROT")
   }
   
   else{ 
      ah <- AnnotationHub()
      orgs <- subset(ah, ah$rdataclass == "OrgDb")
      db <- query(orgs, organism_name)
      db <- db[[1]] # wierd command
      # columns(db)
      # keytypes(db)

      conv_ann <- AnnotationDbi::select(db,
                                        keys=acc,
                                        columns="SYMBOL",
                                        keytype="UNIPROT")
      conv_ann <- filter(conv_ann, !is.na(SYMBOL))
      conv_map <- inner_join(acc.x, conv_ann, by="UNIPROT")
   }

   # # add reverse IDs
   # conv_vec_rev <- paste0("XXX_",conv_vec)
   # names(conv_vec_rev) <- paste0("XXX_",names(conv_vec_rev))
   # conv_vec <- c(conv_vec, conv_vec_rev)
   
   conv_vec <- conv_map$SYMBOL
   names(conv_vec) <- conv_map$accessions

   msnid$accession <- conv_vec[msnid$accession]

   # make sure decoy accessions start with XXX
   msnid$accession <- ifelse(msnid$isDecoy,
                             paste0("XXX_", msnid$accession),
                             msnid$accession)

   return(msnid)
}

