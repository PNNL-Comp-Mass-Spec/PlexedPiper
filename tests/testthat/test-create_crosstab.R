test_that("create_crosstab works", {
  temp_dir <- tempdir()

  n <- 3 # number of files to use (1 or 2 produces unstable results)

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

  out <- create_crosstab(msnid, masic_data, "accession",
                         fractions, samples, references)
  dim(out)

  expect_equal(dim(out), c(4703, 9))
  expect_equal(head(rownames(out)),
               c("AP_004895.1", "NP_001001514.1", "NP_001001515.1",
                 "NP_001002016.2", "NP_001002805.1", "NP_001002807.1"))
  expect_equal(as.numeric(out[1, ]),
               c(0.45092842, -0.16592137, -0.28276820, 0.08842209,
                 -0.20495048, -1.29913791, -0.22078874, -0.53929615,
                 -0.93980622))

  # Clean up
  temp_files <- list.files(temp_dir, full.names = TRUE)
  temp_files <- temp_files[!grepl("\\.fasta$", temp_files)]
  file.remove(temp_files)
  unlink(".Rcache", recursive = TRUE)
  unlink("../../.Rcache", recursive = TRUE)
})

