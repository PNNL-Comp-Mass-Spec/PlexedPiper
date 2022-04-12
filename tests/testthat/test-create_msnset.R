test_that("create_msnset works", {
  temp_dir <- tempdir()

  # Prepare MS/MS IDs
  path_to_MSGF_results <- system.file("extdata/global/msgf_output",
                                      package = "PlexedPiperTestData")
  msgf_file <- list.files(path_to_MSGF_results, full.names = TRUE)[1]
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

  # Prepare table with reporter ion intensities
  path_to_MASIC_results <- system.file("extdata/global/masic_output",
                                       package = "PlexedPiperTestData")
  masic_file <- list.files(path_to_MASIC_results, full.names = TRUE)[1:2]
  file.copy(masic_file, temp_dir)
  masic_data <- read_masic_data(temp_dir, interference_score = TRUE)
  masic_data <- filter_masic_data(masic_data, 0.5, 0)

  # Creating cross-tab
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

  crosstab <- create_crosstab(msnid, masic_data, "accession",
                              fractions, samples, references)

  m <- create_msnset(crosstab, samples)

  expect_true(validObject(m))

  # Clean up
  temp_files <- list.files(temp_dir, full.names = TRUE)
  temp_files <- temp_files[!grepl("\\.fasta$", temp_files)]
  file.remove(temp_files)
  unlink(".Rcache", recursive = TRUE)
  unlink("../../.Rcache", recursive = TRUE)
})

