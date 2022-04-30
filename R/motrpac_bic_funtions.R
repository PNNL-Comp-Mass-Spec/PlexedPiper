#' @title Format Tables for BIC
#'
#' @description Assembles data in format compliant with BIC requirements.
#'
#' @details
#' The `ratio` and `rii` functions require columns "redundantAccessions",
#' "noninferableProteins" and "percentAACoverage" (created with
#' \code{\link[MSnID]{compute_accession_coverage}}) to be present in
#' `psms(msnid)`.
#'
#' * `make_rii_peptide_gl`: returns 'RII_peptide.txt' table (global)
#' * `make_results_ratio_gl`: returns 'results_ratio.txt' table (global)
#' * `make_rii_peptide_ph`: returns 'RII_peptide.txt' table (phospho)
#' * `make_results_ratio_ph`: returns 'results_ratio.txt' table (phospho)
#' * `assess_redundant_protein_matches`: appends proteins matched to multiple
#'   peptides. Creates the "redundantAccessions" column in `psms(msnid)`.
#' * `assess_noninferable_proteins`: appends proteins with identical peptide
#'   sets. Creates the "noninferableProteins" column in `psms(msnid)`.
#'
#' @md
#'
#' @param msnid (MSnID object) final filtered version of MSnID object
#' @param masic_data (object coercible to `data.table`) final filtered version
#'   of MASIC table
#' @param fractions (object coercible to `data.table`) study design table
#'   linking Dataset with PlexID
#' @param samples (object coercible to `data.table`) study design table linking
#'   sample names with TMT channels and PlexID
#' @param references (object coercible to `data.table`) study design table
#'   describing reference value calculation
#' @param annotation (character) format of `accessions(msnid)`. Either
#'   `"refseq"` or `"uniprot"` (case insensitive).
#' @param org_name (character) scientific name of organism (e.g. `"Homo
#'   sapiens"`, `"Rattus norvegicus"`, `"Mus musculus"`, etc.). Case sensitive.
#' @param sep (character) used to concatenate protein, SiteID, and peptide.
#' @param collapse (character) used to collapse proteins in
#'   `assess_redundant_protein_matches`
#'
#' @importFrom MSnID fetch_conversion_table
#' @importFrom dplyr select inner_join left_join mutate %>% case_when rename
#'   group_by summarize if_else arrange
#' @importFrom tibble rownames_to_column
#' @importFrom IRanges IRanges IRangesList reduce
#'
#' @name motrpac_bic_output
#'
#' @examples
#' \dontrun{
#' # Prepare MS/MS IDs ----
#' path_to_MSGF_results <- system.file("extdata/global/msgf_output",
#'                                     package = "PlexedPiperTestData")
#' msnid <- read_msgf_data(path_to_MSGF_results)
#' msnid <- MSnID::correct_peak_selection(msnid)
#' show(msnid)
#' msnid <- filter_msgf_data_peptide_level(msnid, 0.01)
#' show(msnid)
#' path_to_FASTA <- system.file(
#'   "extdata/Rattus_norvegicus_NCBI_RefSeq_2018-04-10.fasta.gz",
#'   package = "PlexedPiperTestData"
#' )
#' msnid <- compute_num_peptides_per_1000aa(msnid, path_to_FASTA)
#' msnid <- filter_msgf_data_protein_level(msnid, 0.01)
#' show(msnid)
#' msnid <- assess_redundant_protein_matches(msnid)
#' msnid <- assess_noninferable_proteins(msnid)
#' fst <- Biostrings::readAAStringSet(path_to_FASTA)
#' names(fst) <- gsub(" .*", "", names(fst)) # make names match accessions
#' msnid <- MSnID::compute_accession_coverage(msnid, fst)
#' head(psms(msnid))
#'
#' # Prepare table with reporter ion intensities ----
#' path_to_MASIC_results <- system.file("extdata/global/masic_output",
#'                                      package = "PlexedPiperTestData")
#' masic_data <- read_masic_data(path_to_MASIC_results,
#'                               interference_score = TRUE)
#' masic_data <- filter_masic_data(masic_data, 0.5, 0)
#'
#' # Read study design files ----
#' library(readr)
#' fractions <- read_tsv(system.file("extdata/study_design/fractions.txt",
#'                                   package = "PlexedPiperTestData"))
#' samples <- read_tsv(system.file("extdata/study_design/samples.txt",
#'                                 package = "PlexedPiperTestData"))
#' references <- read_tsv(system.file("extdata/study_design/references.txt",
#'                                    package = "PlexedPiperTestData"))
#'
#' # Create final tables ----
#' results_ratio <- make_results_ratio_gl(msnid, masic_data,
#'                                        fractions, samples, references,
#'                                        annotation = "RefSeq",
#'                                        org_name = "Rattus norvegicus")
#' head(results_ratio, 10)
#'
#' rii_peptide <- make_rii_peptide_gl(msnid, masic_data,
#'                                    fractions, samples, references,
#'                                    annotation = "RefSeq",
#'                                    org_name = "Rattus norvegicus")
#' head(rii_peptide, 10)
#'
#' # Clean-up cache
#' unlink(".Rcache", recursive = TRUE)
#' }


#' @export
#' @rdname motrpac_bic_output
make_rii_peptide_gl <- function(msnid, masic_data,
                                fractions, samples, references,
                                annotation, org_name = "Rattus norvegicus") {

  ## Make RII study design tables
  if (any(duplicated(samples$ReporterAlias))) {
    samples_rii <- samples %>%
      mutate(MeasurementName = if_else(is.na(MeasurementName),
                                       paste(ReporterAlias, PlexID, sep="_"),
                                       ReporterAlias))
  } else {
    samples_rii <- samples %>%
      mutate(MeasurementName = ReporterAlias)
  }

  references_rii <- references %>%
    mutate(Reference = 1)

  ## Create crosstab
  aggregation_level <- c("accession", "peptide")
  crosstab <- create_crosstab(msnid,
                              masic_data,
                              aggregation_level,
                              fractions, samples_rii, references_rii)
  crosstab <- 2^crosstab # undo log2

  crosstab <- as.data.frame(crosstab) %>%
    rownames_to_column("Specie")

  ## Fetch conversion table
  stopifnot(annotation %in% c("RefSeq", "UniProt"))
  if (annotation == "RefSeq") {
    rgx <- "(^.*)\\.\\d+"
    grp <- "\\1"
  } else if (annotation == "UniProt") {
    rgx <- "((sp|tr)\\|)?([^\\|]*)(.*)?"
    grp <- "\\3"
  }
  conv <- suppressWarnings(fetch_conversion_table(org_name,
                                                  from = toupper(annotation),
                                                  to = c("SYMBOL", "ENTREZID")))

  # Feature data
  feature_data <- crosstab %>%
    select(Specie) %>%
    mutate(protein_id = sub("(^.*)@(.*)", "\\1", Specie),
           sequence = sub("(^.*)@(.*)", "\\2", Specie),
           organism_name = org_name) %>%
    mutate(ANNOTATION = sub(rgx, grp, protein_id)) %>%
    left_join(conv, by=c("ANNOTATION" = toupper(annotation))) %>%
    select(-ANNOTATION) %>%
    rename(gene_symbol = SYMBOL,
           entrez_id = ENTREZID)

  ## Additional info from MS/MS
  ids <- psms(msnid) %>%
    select(accession, peptide, noninferableProteins, MSGFDB_SpecEValue) %>%
    rename(protein_id = accession,
           sequence = peptide,
           redundant_ids = noninferableProteins) %>%
    group_by(protein_id, sequence, redundant_ids) %>%
    summarize(peptide_score = min(MSGFDB_SpecEValue)) %>%
    mutate(is_contaminant = grepl("Contaminant", protein_id))

  # Final table
  rii_peptide <- feature_data %>%
    inner_join(ids, by=c("protein_id", "sequence")) %>%
    inner_join(crosstab, by="Specie") %>%
    select(-Specie)
}

utils::globalVariables(c("MeasurementName", "ReporterAlias", "PlexID", "Specie",
                         "protein_id", "ANNOTATION", "SYMBOL", "ENTREZID",
                         "accession", "peptide", "noninferableProteins",
                         "MSGFDB_SpecEValue", "redundant_ids"))




#' @export
#' @rdname motrpac_bic_output
make_results_ratio_gl <- function(msnid, masic_data,
                                  fractions, samples, references,
                                  annotation, org_name = "Rattus norvegicus") {
  ## Create crosstab
  aggregation_level <- c("accession")
  crosstab <- create_crosstab(msnid, masic_data, aggregation_level, fractions,
                              samples, references)

  crosstab <- as.data.frame(crosstab) %>%
    rownames_to_column("protein_id")

  ## Fetch conversation table
  stopifnot(annotation %in% c("RefSeq", "UniProt"))
  if (annotation == "RefSeq") {
    rgx <- "(^.*)\\.\\d+"
    grp <- "\\1"
  } else if (annotation == "UniProt") {
    rgx <- "((sp|tr)\\|)?([^\\|]*)(.*)?"
    grp <- "\\3"
  }
  conv <- suppressWarnings(fetch_conversion_table(org_name,
                                                  from = toupper(annotation),
                                                  to = c("SYMBOL", "ENTREZID")))

  ## Create feature data
  feature_data <- crosstab %>%
    select(protein_id) %>%
    mutate(organism_name = org_name) %>%
    mutate(ANNOTATION = sub(rgx, grp, protein_id)) %>%
    left_join(conv, by=c("ANNOTATION" = toupper(annotation))) %>%
    select(-ANNOTATION) %>%
    rename(gene_symbol = SYMBOL,
           entrez_id = ENTREZID)


  ## Additional info from MS/MS
  ids <- psms(msnid) %>%
    select(accession, peptide, noninferableProteins, percentAACoverage,
           MSGFDB_SpecEValue) %>%
    rename(protein_id = accession,
           sequence = peptide,
           redundant_ids = noninferableProteins,
           percent_coverage = percentAACoverage) %>%
    group_by(protein_id, sequence, redundant_ids, percent_coverage) %>%
    summarize(peptide_score = min(MSGFDB_SpecEValue)) %>%
    group_by(protein_id, redundant_ids, percent_coverage) %>%
    summarize(protein_score = min(peptide_score),
              num_peptides = n()) %>%
    mutate(is_contaminant = grepl("Contaminant", protein_id))

  results_ratio <- feature_data %>%
    inner_join(ids, by="protein_id") %>%
    inner_join(crosstab, by="protein_id")
}

utils::globalVariables(c("noninferableProteins", "percentAACoverage",
                         "percent_coverage"))



#' @export
#' @rdname motrpac_bic_output
make_rii_peptide_ph <- function(msnid, masic_data,
                                fractions, samples, references,
                                annotation, org_name = "Rattus norvegicus",
                                sep = "_") {
  ## Make RII study design tables
  if (any(duplicated(samples$ReporterAlias))) {
    samples_rii <- samples %>%
      mutate(MeasurementName = if_else(is.na(MeasurementName),
                                       paste(ReporterAlias, PlexID, sep="_"),
                                       ReporterAlias))
  } else {
    samples_rii <- samples %>%
      mutate(MeasurementName = ReporterAlias)
  }

  references_rii <- references %>%
    mutate(Reference = 1)

  ## Create crosstab
  aggregation_level <- c("accession", "peptide", "SiteID")
  crosstab <- create_crosstab(msnid,
                              masic_data,
                              aggregation_level,
                              fractions, samples_rii, references_rii)

  crosstab <- 2^crosstab # undo log2

  crosstab <- as.data.frame(crosstab) %>%
    rownames_to_column("Specie")

  ## Fetch conversation table
  stopifnot(annotation %in% c("RefSeq", "UniProt"))
  if (annotation == "RefSeq") {
    rgx <- "(^.*)\\.\\d+"
    grp <- "\\1"
  } else if (annotation == "UniProt") {
    rgx <- "((sp|tr)\\|)?([^\\|]*)(.*)?"
    grp <- "\\3"
  }
  conv <- suppressWarnings(fetch_conversion_table(org_name,
                                                  from = toupper(annotation),
                                                  to = c("SYMBOL", "ENTREZID")))


  ## Create RII peptide table
  feature_data <- crosstab %>%
    select(Specie) %>%
    mutate(protein_id = sub("(^.*)@(.*)@(.*)", "\\1", Specie),
           sequence = sub("(^.*)@(.*)@(.*)", "\\2", Specie),
           ptm_id = sub("(^.*)@(.*)@(.*)", "\\3", Specie)) %>%
    mutate(ptm_peptide = paste(ptm_id, sequence, sep=sep),
           organism_name = org_name) %>%
    mutate(ANNOTATION = sub(rgx, grp, protein_id)) %>%
    left_join(conv, by=c("ANNOTATION" = toupper(annotation))) %>%
    select(-ANNOTATION) %>%
    rename(gene_symbol = SYMBOL,
           entrez_id = ENTREZID)

  ## Additional info from MS/MS
  ids <- psms(msnid) %>%
    select(protein_id=accession,
           sequence=peptide,
           ptm_id=SiteID,
           flanking_sequence=sequenceWindow,
           redundant_ids=redundantAccessions,
           MSGFDB_SpecEValue,
           maxAScore) %>%
    group_by(protein_id, sequence, ptm_id, flanking_sequence, redundant_ids) %>%
    summarize(peptide_score = min(MSGFDB_SpecEValue),
              confident_score = max(maxAScore)) %>%
    mutate(confident_site = case_when(confident_score >= 17 ~ TRUE,
                                      confident_score < 17 ~ FALSE),
           is_contaminant = grepl("Contaminant", protein_id))

  rii_peptide <- feature_data %>%
    inner_join(ids, by = c("protein_id", "sequence", "ptm_id")) %>%
    mutate(ptm_id = gsub("-", sep, ptm_id),
           ptm_peptide = gsub("-", sep, ptm_peptide)) %>%
    inner_join(crosstab, by="Specie") %>%
    select(-Specie)
}

utils::globalVariables(c("ptm_peptide", "MeasurementName", "ReporterAlias",
                         "PlexID", "Specie", "ptm_id", "protein_id",
                         "ANNOTATION", "SYMBOL", "ENTREZID", "accession",
                         "peptide", "SiteID", "sequenceWindow",
                         "redundantAccessions", "MSGFDB_SpecEValue",
                         "maxAScore", "flanking_sequence", "redundant_ids"))




#' @export
#' @rdname motrpac_bic_output
make_results_ratio_ph <- function(msnid, masic_data,
                                  fractions, samples, references,
                                  annotation, org_name = "Rattus norvegicus",
                                  sep = "_") {

  aggregation_level <- c("accession", "SiteID")
  crosstab <- create_crosstab(msnid, masic_data, aggregation_level, fractions,
                              samples, references)
  crosstab <- as.data.frame(crosstab) %>%
    rownames_to_column("Specie")

  ## Fetch conversation table
  stopifnot(annotation %in% c("RefSeq", "UniProt"))
  if (annotation == "RefSeq") {
    rgx <- "(^.*)\\.\\d+"
    grp <- "\\1"
  } else if (annotation == "UniProt") {
    rgx <- "((sp|tr)\\|)?([^\\|]*)(.*)?"
    grp <- "\\3"
  }
  conv <- suppressWarnings(fetch_conversion_table(org_name,
                                                  from = toupper(annotation),
                                                  to = c("SYMBOL", "ENTREZID")))

  ## Create RII peptide table
  feature_data <- crosstab %>%
    select(Specie) %>%
    mutate(protein_id = sub("(^.*)@(.*)", "\\1", Specie),
           ptm_id = sub("(^.*)@(.*)", "\\2", Specie),
           organism_name = org_name) %>%
    mutate(ANNOTATION = sub(rgx, grp, protein_id)) %>%
    left_join(conv, by=c("ANNOTATION" = toupper(annotation))) %>%
    select(-ANNOTATION) %>%
    rename(gene_symbol = SYMBOL,
           entrez_id = ENTREZID)

  ## Additional info from MS/MS
  ids <- psms(msnid) %>%
    select(protein_id = accession,
           sequence = peptide,
           ptm_id = SiteID,
           redundant_ids = noninferableProteins,
           flanking_sequence = sequenceWindow,
           MSGFDB_SpecEValue,
           maxAScore) %>%
    # group at peptide level to calculate peptide score, confident score
    group_by(protein_id, sequence, ptm_id, flanking_sequence, redundant_ids) %>%
    summarize(peptide_score = min(MSGFDB_SpecEValue),
              confident_score = max(maxAScore)) %>%
    # regroup at siteID level and recalculate ptm score
    group_by(protein_id, ptm_id, flanking_sequence, redundant_ids) %>%
    summarize(ptm_score = min(peptide_score),
              confident_score = max(confident_score)) %>%
    mutate(confident_site = case_when(confident_score >= 17 ~ TRUE,
                                      confident_score < 17 ~ FALSE),
           is_contaminant = grepl("Contaminant", protein_id))

  results_ratio <- feature_data %>%
    inner_join(ids, by=c("protein_id", "ptm_id")) %>%
    mutate(ptm_id = gsub("-", sep, ptm_id)) %>%
    inner_join(crosstab, by="Specie") %>%
    select(-Specie)
}

utils::globalVariables(c("flanking_sequence", "redundant_ids", "peptide_score",
                         "confident_score", "noninferableProteins"))




#' @export
#' @rdname motrpac_bic_output
assess_redundant_protein_matches <- function(msnid, collapse="|") {

  res <- psms(msnid) %>%
    select(accession, peptide) %>%
    distinct() %>%
    group_by(peptide) %>%
    summarize(redundantAccessions = paste(accession, collapse=collapse))

  res <- left_join(psms(msnid), res, by="peptide")

  psms(msnid) <- res

  return(msnid)
}




#' @export
#' @rdname motrpac_bic_output
assess_noninferable_proteins <- function(msnid, collapse="|") {

  # assign each accession to its peptide signature
  res <- psms(msnid) %>%
    select(accession, peptide) %>%
    distinct() %>%
    group_by(accession) %>%
    arrange(peptide) %>%
    summarize(peptideSignature = paste(peptide, collapse=collapse))

  # group together accessions with identical peptide signature
  res <- res %>%
    group_by(peptideSignature) %>%
    summarize(noninferableProteins = paste(accession, collapse=collapse)) %>%
    left_join(res, by="peptideSignature") %>%
    select(-peptideSignature)

  res <- left_join(psms(msnid), res, by="accession")

  psms(msnid) <- res

  return(msnid)
}

utils::globalVariables(c("peptideSignature"))

