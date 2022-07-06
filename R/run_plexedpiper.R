#' @title Wrapper for PlexedPiper
#'
#' @description A wrapper function for PlexedPiper that performs all the steps
#'   necessary to go from MS-GF+ and MASIC data to Reporter Ion Intensity (RII)
#'   and Ratio Results tables used by the MoTrPAC Bioinformatics Center (BIC).
#'
#' @param msgf_output_folder (character) Path to MSGF+ results folder
#' @param fasta_file (character) Path to FASTA file
#' @param masic_output_folder (character) MASIC results folder
#' @param ascore_output_folder (character) AScore results folder
#' @param proteomics (character) Either "pr" - proteomics, "ph" -
#'   phosphoproteomics, "ub" - ubiquitinomics, or "ac" - acetylomics
#' @param study_design_folder (character) Folder containing the three study
#'   design tables: fractions.txt, samples.txt, and references.txt
#' @param prefix (character) optional prefix for study design tables.
#' @param species (character) Scientific name of species (e.g. "Rattus
#'   norvegicus", "Homo sapiens", etc.)
#' @param annotation (character) Source for annotations: either `RefSeq` or
#'   `UniProt` (case insensitive).
#' @param global_results (character) Only for PTM experiments. Ratio results
#'   from a global protein abundance experiment. If provided, prioritized
#'   inference will be performed. Otherwise, parsimonious inference is performed
#'   without a prior. See \code{\link[MSnID]{infer_parsimonious_accessions}} for
#'   details.
#' @param unique_only (logical) Whether to discard peptides that match multiple
#'   proteins in the parsimonious protein inference step. Default \code{FALSE}.
#' @param refine_prior (logical) if \code{FALSE} (default), peptides are allowed
#'   to match multiple proteins in the prior. That is, the greedy set cover
#'   algorithm is only applied to the set of proteins not in the prior. If
#'   \code{TRUE}, the algorithm is applied to the prior and non-prior sets
#'   separately before combining. See
#'   \code{\link[MSnID]{infer_parsimonious_accessions}} for more details.
#' @param write_results_to_file (logical) Whether to write the results to files.
#' @param output_folder (character) Output folder name to save results. If not
#'   provided, it will save it to the current directory.
#' @param file_prefix (character) Prefix for the file name outputs.
#' @param save_env (logical) Whether to save the R environment to the output
#'   folder.
#' @param return_results (logical) Whether to return the ratio and rii results.
#' @param verbose (logical) Whether to show messages describing steps.
#'
#' @return (list) If `return_results` is `TRUE`, it returns list with ratio and
#'   RII data frames.
#'
#'
#' @md
#'
#' @importFrom Biostrings readAAStringSet
#' @importFrom utils read.table write.table
#' @importFrom MSnID compute_accession_coverage correct_peak_selection
#'   extract_sequence_window infer_parsimonious_accessions map_mod_sites
#'
#' @examples \dontrun{
#' results <- run_plexedpiper(msgf_output_folder = "~/path/to/msgfplus/",
#'                            fasta_file  = "~/path/to/fasta/sequence.fasta",
#'                            masic_output_folder = "~/path/to/masic-results/",
#'                            ascore_output_folder = "~/path/to/ascore-results/",
#'                            proteomics = "ph",
#'                            study_design_folder = "~/path/to/study-design/",
#'                            species = "Rattus norvegicus",
#'                            annotation = "RefSeq",
#'                            global_results = "~/path/to/global/ratio.txt",
#'                            unique_only = FALSE,
#'                            refine_prior = FALSE,
#'                            output_folder = "~/path/to/pp-results/",
#'                            file_prefix = "msgfplus-pp-results",
#'                            return_results = TRUE,
#'                            verbose = TRUE)
#' }
#' @export


run_plexedpiper <- function(msgf_output_folder,
                            fasta_file,
                            masic_output_folder,
                            ascore_output_folder = NULL,
                            proteomics,
                            study_design_folder,
                            prefix = character(0),
                            species,
                            annotation,
                            file_prefix = NULL,
                            unique_only = FALSE,
                            global_results = NULL,
                            refine_prior = FALSE,
                            write_results_to_file = TRUE,
                            output_folder = NULL,
                            save_env = FALSE,
                            return_results = FALSE,
                            verbose = TRUE) {

  if( is.null(write_results_to_file) & is.null(return_results) ){
    stop("\nProvide either <write_results_to_file = TRUE> or <return_results = TRUE> or both. Both cannot be FALSE.")
  }

  if (verbose) {
    message("Running PlexedPiper (the MSGF+ pipeline wrapper), with the following parameters:")
    message("- Proteomics experiment:\"", proteomics, "\"")
    message("- Species: \"", species, "\"")
    message("- Annotation: \"", annotation, "\"")
  }

  if(is.null(file_prefix)){
    file_prefix <- paste0("MSGFPLUS_", toupper(proteomics))
    if(!is.null(global_results)){
      file_prefix <- paste0(file_prefix,"-ip")
    }
  }

  # Data loading
  message("- Fetch study design tables")

  study_design <- read_study_design(study_design_folder, prefix = prefix)
  msnid <- read_msgf_data(msgf_output_folder)
  if(!is.null(ascore_output_folder)) ascore <- read_AScore_results(ascore_output_folder)
  masic_data <- read_masic_data(masic_output_folder,
                                interference_score = TRUE)

  fst <- Biostrings::readAAStringSet(fasta_file)
  names(fst) <- sub("^(\\S*)\\s.*", "\\1", names(fst))

  if (verbose) message("- Filtering MS-GF+ results.")
  if (proteomics == "pr") {
    if (verbose) message("   + Correct for isotope selection error")
    msnid <- correct_peak_selection(msnid)
  }

  if (proteomics != "pr") {
    if (verbose) message("   + Select best PTM location by AScore")
    msnid <- best_PTM_location_by_ascore(msnid, ascore)

    if (verbose) message("   + Apply PTM filter")
    if (proteomics == "ph") {
      reg.expr <- "grepl(\"\\\\*\", peptide)"
    } else if (proteomics == "ac") {
      reg.expr <- "grepl(\"#\", peptide)"
    } else if (proteomics == "ub") {
      psms(msnid) <- psms(msnid) %>% mutate(peptide = gsub("@", "#", peptide))
      reg.expr <- "grepl(\"#\", peptide)"
    }
    msnid <- apply_filter(msnid, reg.expr)
  }

  if(verbose) message("   + Peptide-level FDR filter")
  msnid <- filter_msgf_data(msnid, level = "peptide", fdr.max = 0.01)

  if (proteomics == "pr") {
    if (verbose) message("   + Protein-level FDR filter")
    msnid <- compute_num_peptides_per_1000aa(msnid, path_to_FASTA = fasta_file)

    msnid <- filter_msgf_data(msnid, level = "accession", fdr.max = 0.01)
  }

  if(verbose) message("   + Remove decoy sequences")
  msnid <- apply_filter(msnid, "!isDecoy")

  if (verbose) message("   + Concatenating redundant RefSeq matches")
  msnid <- assess_redundant_protein_matches(msnid, collapse = ",")

  if (verbose) message("   + Assessing non-inferable proteins")
  msnid <- assess_noninferable_proteins(msnid, collapse = ",")

  if(verbose) message("   + Inference of parsimonius set")
  if (proteomics == "pr") {
    prior <- character(0)
  } else if (is.null(global_results)) {
    if(verbose) message("     > Reference global proteomics dataset NOT provided")
    prior <- character(0)
  } else {
    global_ratios <- read.table(global_results, header=TRUE, sep="\t")
    prior <- unique(global_ratios$protein_id)
  }

  msnid <- infer_parsimonious_accessions(msnid,
                                         unique_only = unique_only,
                                         prior = prior,
                                         refine_prior = refine_prior)

  if (proteomics == "pr") {
    if (verbose) message("   + Compute protein coverage")
    msnid <- compute_accession_coverage(msnid, fst)
  }

  if (proteomics != "pr") {
    if(verbose) message("   + Mapping sites to protein sequence")
    if (proteomics == "ph") {
      mod_char = "*"
    } else if (proteomics %in% c("ac", "ub")) {
      mod_char = "#"
    } else {
      stop("Proteomics variable not supported.")
    }

    msnid <- map_mod_sites(msnid,
                           fasta           = fst,
                           accession_col   = "accession",
                           peptide_mod_col = "peptide",
                           mod_char        = mod_char,
                           site_delimiter  = "lower")

    if(verbose) message("   + Map flanking sequences")
    msnid <- extract_sequence_window(msnid, fasta = fst)
  }

  if (verbose) message("- Filtering MASIC results.")
  masic_data <- filter_masic_data(masic_data, 0.5, 0)

  args <- list(msnid      = msnid,
               masic_data = masic_data,
               fractions  = study_design$fractions,
               samples    = study_design$samples,
               references = study_design$references,
               annotation = annotation,
               org_name   = species)

  if (verbose) message("- Making results tables.")
  if (proteomics == "pr") {
    suppressMessages(rii_peptide <- do.call(make_rii_peptide_gl, args))
    suppressMessages(results_ratio <- do.call(make_results_ratio_gl, args))
  } else {
    suppressMessages(rii_peptide   <- do.call(make_rii_peptide_ph,   args))
    suppressMessages(results_ratio <- do.call(make_results_ratio_ph, args))
  }

  if (verbose) message("- Saving results.")

  if(is.null(output_folder)){
    output_folder <- getwd()
  }

  if(write_results_to_file){
    if (!is.null(output_folder)) {
      if (!dir.exists(file.path(output_folder))) {
        dir.create(file.path(output_folder), recursive = TRUE)
      }

      filename = paste0(file_prefix, "-results_RII-peptide.txt")
      write.table(rii_peptide,
                  file      = file.path(output_folder, filename),
                  sep       = "\t",
                  row.names = FALSE,
                  quote     = FALSE)
      if(verbose) message("- RII file save to ",
                          file.path(output_folder, filename))

      filename = paste0(file_prefix, "-results_ratio.txt")
      write.table(results_ratio,
                  file      = file.path(output_folder, filename),
                  sep       = "\t",
                  row.names = FALSE,
                  quote     = FALSE)
      if(verbose) message("- RATIO file save to ",
                          file.path(output_folder, filename))
    }
  }
  if (save_env) {
    fileenv = paste0(file_prefix, "-env.RData")
    save.image(file = file.path(output_folder, fileenv))
    if(verbose) message("- R environment saved to ",
                        file.path(output_folder, fileenv))
  }

  if (verbose) message("Done!")

  if(return_results){
    return(list(rii_peptide = rii_peptide, results_ratio = results_ratio))
  }
}

