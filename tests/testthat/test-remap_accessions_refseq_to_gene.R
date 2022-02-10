test_that("remap_accessions_refseq_to_gene works", {
  path_to_MSGF_results <- system.file("extdata/global/msgf_output",
                                      package = "PlexedPiperTestData")
  temp_dir <- tempdir()
  file_to_copy <- list.files(path_to_MSGF_results, full.names = TRUE)[1]
  file.copy(file_to_copy, temp_dir)
  msnid <- read_msgf_data(temp_dir)

  suppressMessages(
    out <- remap_accessions_refseq_to_gene(msnid, "Rattus norvegicus")
  )
  expect_equal(head(accessions(out)),
               c("Myoz1", "Myom1", NA, "Slc25a4", "Slc25a5", "Pkm"))

  # Clean up
  temp_files <- list.files(temp_dir, full.names = TRUE)
  temp_files <- temp_files[!grepl("\\.fasta$", temp_files)]
  file.remove(temp_files)
  unlink(".Rcache", recursive = TRUE)
  unlink("../../.Rcache", recursive = TRUE)
})

