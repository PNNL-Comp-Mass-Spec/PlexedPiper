# Check that number of rows after filtering is correct
test_that("Filter MSGF data at peptide level works", {
    path_to_MSGF_results <- system.file("extdata/global/msgf_output",
                                        package = "PlexedPiperTestData")
    temp_dir <- tempdir()
    file_to_copy <- list.files(path_to_MSGF_results, full.names = TRUE)[1]
    file.copy(file_to_copy, temp_dir)
    msnid <- read_msgf_data(temp_dir)

    msnid <- filter_msgf_data(msnid, level = "peptide")
    expect_equal(length(peptides(msnid)), 4197)

    # Clean up
    temp_files <- list.files(temp_dir, full.names = TRUE)
    temp_files <- temp_files[!grepl("\\.fasta$", temp_files)]
    file.remove(temp_files)
    unlink(".Rcache", recursive = TRUE)
    unlink("../../.Rcache", recursive = TRUE)
})
