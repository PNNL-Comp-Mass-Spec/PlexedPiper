#' Remap Sequence IDs in FASTA File
#'
#' The function remaps the IDs in the FASTA file from RefSeq to Gene symbols.
#' In case of multiple RefSeq sequences available, only the first longest is
#' retained.
#'
#' @param path_to_FASTA (MSnID object) MS/MS ID data
#' @param organism_name (string) Official organism name
#' @param conversion_table (data.frame) data frame with two columns
#' one should with named accessions and contain accessions from `msnid` object
#' (e.g. RefSeq) the other is with alternative annotation to map to
#' (e.g. gene symbol).
#' @return (MSnID object) MS/MS ID data with computed number of
#' peptides per 1000 aa. Added column name - "peptides_per_1000aa".
#' @importMethodsFrom MSnID $<-
#' @importFrom Biostrings readAAStringSet width writeXStringSet
#' @importFrom AnnotationHub AnnotationHub query
#' @importFrom AnnotationDbi select
#' @importFrom dplyr bind_cols filter inner_join slice pull arrange
#' @importFrom tools file_path_sans_ext file_ext
#' @importFrom tibble tibble
#' @importFrom tools file_path_sans_ext
#'
#' @name remap_accessions_refseq_to_gene_fasta
#'
#' @examples \dontrun{
#' #' path_to_FASTA <- system.file("extdata/Rattus_norvegicus_NCBI_RefSeq_2018-04-10.fasta.gz", 
#'                                 package = "PlexedPiperTestData")
#' temp_work_dir <- tempdir() # can be set to "." or getwd(), if done carefully
#' file.copy(path_to_FASTA, temp_work_dir)
#' path_to_FASTA <- file.path(temp_work_dir, basename(path_to_FASTA))
#' library(Biostrings)
#' readAAStringSet(path_to_FASTA) # refseq IDs
#' path_to_new_FASTA <- remap_accessions_refseq_to_gene_fasta(path_to_FASTA,"Rattus norvegicus")
#' readAAStringSet(path_to_new_FASTA) # gene IDs
#' }
#' @export
remap_accessions_refseq_to_gene_fasta <- function(path_to_FASTA,
                                                  organism_name,
                                                  conversion_table){
                                            

   is_compressed <- FALSE
   if(grepl("[.]gz$", path_to_FASTA)){
      is_compressed <- TRUE
   }else if(grepl("[.](bz2|xz|zip)$", path_to_FASTA)){
      stop("The only supported compression is gzip!")
   }

   mySequences <- readAAStringSet(path_to_FASTA)
   names(mySequences) <- sub("^(.P_\\d+)(\\.\\d+)?\\s.*","\\1",names(mySequences))
   prot_lengths <- data.frame(REFSEQ = names(mySequences),
                              Length = width(mySequences),
                              stringsAsFactors = FALSE)

   if(missing(conversion_table)){
      ah <- AnnotationHub()
      orgs <- subset(ah, ah$rdataclass == "OrgDb")
      db <- query(orgs, organism_name)
      db <- db[[1]]
      conversion_table <- AnnotationDbi::select(db,
                                       keys=names(mySequences),
                                       columns="SYMBOL",
                                       keytype="REFSEQ")
      conversion_table <- conversion_table[,c("REFSEQ","SYMBOL")]
   }
   # I don't know what column names are, but I require that refseq is first,
   # gene second
   colnames(conversion_table) <- c("REFSEQ","SYMBOL")
   # resolving mapping ambiguity by selecting longest RefSeq per gene
   conversion_table <- conversion_table %>%
      filter(!is.na(SYMBOL)) %>%
      inner_join(prot_lengths, by = "REFSEQ") %>%
      group_by(SYMBOL) %>%
      slice(which.max(Length))
   conversion_vec <- conversion_table %>% pull(2)
   names(conversion_vec) <- conversion_table %>% pull(1)
   mySequences <- mySequences[names(conversion_vec)]
   names(mySequences) <- conversion_vec[names(mySequences)]

   file_no_ext <- file_path_sans_ext(path_to_FASTA, compression=TRUE)
   ext <- sub(file_no_ext, "", path_to_FASTA, fixed=TRUE)

   path_to_FASTA_gene <- paste0(file_no_ext, '_gene', ext)

   writeXStringSet(mySequences, path_to_FASTA_gene, compress = is_compressed)
   return(path_to_FASTA_gene)
}



#' Remap Sequence IDs from Uniprot in FASTA File
#'
#' The function remaps the IDs in the FASTA file from Uniprot to Gene symbols.
#' In case of multiple Uniprot sequences available, only the first longest is
#' retained.
#'
#' @param path_to_FASTA (MSnID object) MS/MS ID data
#' 
#' @return (MSnID object) MS/MS ID data with computed number of
#' peptides per 1000 aa. Added column name - "peptides_per_1000aa".
#' 
#' @importFrom Biostrings readAAStringSet
#' @name remap_accessions_refseq_to_gene_fasta
#'
#' @export
remap_accessions_uniprot_to_gene_fasta <- function(path_to_FASTA){

   is_compressed <- FALSE
   if(grepl("[.]gz$", path_to_FASTA)){
      is_compressed <- TRUE
   }else if(grepl("[.](bz2|xz|zip)$", path_to_FASTA)){
      stop("The only supported compression is gzip!")
   }

   mySequences <- readAAStringSet(path_to_FASTA)
   # subsetting to only with known genes
   mySequences <- mySequences[grepl("GN=",names(mySequences))]
   genes <- sub(".*GN=([^ ]*).*","\\1",names(mySequences))

   # retain only the longest form per gene
   top_len <- tibble(genes, len=width(mySequences), ind=seq_along(genes)) %>%
      group_by(genes) %>%
      dplyr::slice(which.max(len)) %>%
      ungroup() %>%
      arrange(ind)

   mySequences <- mySequences[top_len$ind]
   names(mySequences) <- top_len$genes

   file_no_ext <- file_path_sans_ext(path_to_FASTA, compression=TRUE)
   ext <- sub(file_no_ext, "", path_to_FASTA)
   path_to_FASTA_gene <- paste0(file_no_ext, '_gene', ext)
   writeXStringSet(mySequences, path_to_FASTA_gene, compress = is_compressed)
   return(path_to_FASTA_gene)
}


