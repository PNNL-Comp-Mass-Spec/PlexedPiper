test_that("assess_noninferable_proteins works", {
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

  ## assess_noninferable_proteins
  msnid <- assess_noninferable_proteins(msnid)

  expect_true("noninferableProteins" %in% names(msnid))
  expect_equal(head(unique(msnid$noninferableProteins)),
               c("NP_001102567.1|XP_006251717.1",
                 "XP_017447817.1",
                 "XP_017456475.1",
                 "NP_445967.1",
                 "NP_476443.1",
                 "XP_008771999.1"))

  # Clean up
  temp_files <- list.files(temp_dir, full.names = TRUE)
  temp_files <- temp_files[!grepl("\\.fasta$", temp_files)]
  file.remove(temp_files)
  unlink(".Rcache", recursive = TRUE)
  unlink("../../.Rcache", recursive = TRUE)
})

