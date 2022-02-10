test_that("global MoTrPAC BIC functions work", {
  temp_dir <- tempdir()

  n <- 3

  # Prepare MS/MS IDs
  path_to_MSGF_results <- system.file("extdata/global/msgf_output",
                                      package = "PlexedPiperTestData")
  msgf_file <- list.files(path_to_MSGF_results, full.names = TRUE)[1:n]
  file.copy(msgf_file, temp_dir)
  msnid <- read_msgf_data(temp_dir)

  msnid <- MSnID::correct_peak_selection(msnid)
  msnid <- filter_msgf_data_peptide_level(msnid, 0.01)
  path_to_FASTA <- system.file(
    "extdata/Rattus_norvegicus_NCBI_RefSeq_2018-04-10.fasta.gz",
    package = "PlexedPiperTestData"
  )
  msnid <- compute_num_peptides_per_1000aa(msnid, path_to_FASTA)
  msnid <- filter_msgf_data_protein_level(msnid, 0.01)
  msnid <- assess_redundant_protein_matches(msnid)
  msnid <- assess_noninferable_proteins(msnid)
  fst <- Biostrings::readAAStringSet(path_to_FASTA)
  names(fst) <- gsub(" .*", "", names(fst)) # make names match accessions
  msnid <- MSnID::compute_accession_coverage(msnid, fst)

  # Prepare table with reporter ion intensities
  path_to_MASIC_results <- system.file("extdata/global/masic_output",
                                       package = "PlexedPiperTestData")
  masic_file <- list.files(path_to_MASIC_results, full.names = TRUE)[1:(2*n)]
  file.copy(masic_file, temp_dir)
  masic_data <- read_masic_data(temp_dir, interference_score = TRUE)
  masic_data <- filter_masic_data(masic_data, 0.5, 0)


  # Creating cross-tab
  library(readr)
  suppressMessages({
    fractions <- read_tsv(system.file("extdata/study_design/fractions.txt",
                                      package = "PlexedPiperTestData"))
    samples <- read_tsv(system.file("extdata/study_design/samples.txt",
                                    package = "PlexedPiperTestData"))
    references <- read_tsv(system.file("extdata/study_design/references.txt",
                                       package = "PlexedPiperTestData"))
  })
  # Modify study design to avoid
  # Error in eval(parse(text = ref_i$Reference)) : object 'ref' not found
  fractions <- fractions[1, ]
  samples <- samples[samples$PlexID == "S1", ]
  references <- references[references$PlexID == "S1", ]

  ## RII peptide table ----
  rii_peptide <- make_rii_peptide_gl(msnid, masic_data,
                                     fractions, samples, references,
                                     annotation = "RefSeq",
                                     org_name = "Rattus norvegicus")

  expect_true(all(c("protein_id", "sequence", "organism_name",
                    "gene_symbol", "entrez_id", "redundant_ids",
                    "peptide_score", "is_contaminant") %in%
                    names(rii_peptide)))
  expect_equal(dim(rii_peptide), c(17674, 18))

  ## Results ratio ----
  results_ratio <- make_results_ratio_gl(msnid, masic_data,
                                         fractions, samples, references,
                                         annotation = "RefSeq",
                                         org_name = "Rattus norvegicus")

  expect_true(all(c("protein_id", "organism_name", "gene_symbol",
                    "entrez_id", "redundant_ids", "percent_coverage",
                    "protein_score", "num_peptides", "is_contaminant") %in%
                    names(results_ratio)))
  expect_equal(dim(results_ratio), c(4703, 18))

  # Clean up
  temp_files <- list.files(temp_dir, full.names = TRUE)
  temp_files <- temp_files[!grepl("\\.fasta$", temp_files)]
  file.remove(temp_files)
  unlink(".Rcache", recursive = TRUE)
  unlink("../../.Rcache", recursive = TRUE)
})

