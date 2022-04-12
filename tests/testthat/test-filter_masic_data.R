test_that("filter_masic_data works", {
  temp_dir <- tempdir()
  path_to_MASIC_results <- system.file("extdata/global/masic_output",
                                       package = "PlexedPiperTestData")
  masic_file <- list.files(path_to_MASIC_results, full.names = TRUE)[1:2]
  file.copy(masic_file, temp_dir)
  masic_data <- read_masic_data(temp_dir, interference_score = TRUE)

  expect_true(all(c("Dataset", "ScanNumber", "InterferenceScore") %in%
                    names(masic_data)))

  masic_data <- filter_masic_data(masic_data, 0.5, 0)

  expect_equal(dim(masic_data), c(30893, 12))

  # Clean up
  temp_files <- list.files(temp_dir, full.names = TRUE)
  temp_files <- temp_files[!grepl("\\.fasta$", temp_files)]
  file.remove(temp_files)
  unlink(".Rcache", recursive = TRUE)
  unlink("../../.Rcache", recursive = TRUE)
})

