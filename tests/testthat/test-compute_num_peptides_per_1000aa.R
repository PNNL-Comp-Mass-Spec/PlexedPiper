test_that("compute peptides per 1000aa works", {
  path_to_MSGF_results <- system.file("extdata/global/msgf_output",
                                      package = "PlexedPiperTestData")
  temp_dir <- tempdir()
  file_to_copy <- list.files(path_to_MSGF_results, full.names = TRUE)[1]
  file.copy(file_to_copy, temp_dir)
  msnid <- read_msgf_data(temp_dir)
  msnid <- filter_msgf_data_peptide_level(msnid)
  path_to_FASTA <- system.file(
    "extdata/Rattus_norvegicus_NCBI_RefSeq_2018-04-10.fasta.gz",
    package = "PlexedPiperTestData"
  )
  msnid <- compute_num_peptides_per_1000aa(msnid, path_to_FASTA)
  out <- head(unique(msnid$peptides_per_1000aa))

  expect_equal(out, 1000 * c(4/296, 14/1689, 305/34232,
                             307/34695, 5/298, 3/273))

  # Clean up
  temp_files <- list.files(temp_dir, full.names = TRUE)
  temp_files <- temp_files[!grepl("\\.fasta$", temp_files)]
  file.remove(temp_files)
  unlink(".Rcache", recursive = TRUE)
  unlink("../../.Rcache", recursive = TRUE)
})

