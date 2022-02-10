# Check for numMissCleavages column and that filter worked properly
test_that("filter_missed_cleavages adds numMissCleavages", {
  path_to_MSGF_results <- system.file("extdata/global/msgf_output",
                                      package = "PlexedPiperTestData")
  temp_dir <- tempdir()
  file_to_copy <- list.files(path_to_MSGF_results, full.names = TRUE)[1]
  file.copy(file_to_copy, temp_dir)
  msnid <- read_msgf_data(temp_dir)

  msnid <- filter_msgf_data_missed_cleavages(msnid,
                                             missed_cleavages_threshold = 2)

  expect_true("numMissCleavages" %in% names(msnid))
  expect_equal(nrow(msnid), 94123)

  # Clean up
  temp_files <- list.files(temp_dir, full.names = TRUE)
  temp_files <- temp_files[!grepl("\\.fasta$", temp_files)]
  file.remove(temp_files)
  unlink(".Rcache", recursive = TRUE)
  unlink("../../.Rcache", recursive = TRUE)
})

