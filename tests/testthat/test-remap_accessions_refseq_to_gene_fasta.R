test_that("remap_accessions_refseq_to_gene_fasta works", {
  library(Biostrings)

  path_to_FASTA <- system.file(
    "extdata/Rattus_norvegicus_NCBI_RefSeq_2018-04-10.fasta.gz",
    package = "PlexedPiperTestData"
  )
  temp_dir <- tempdir()
  file.copy(path_to_FASTA, temp_dir)
  path_to_FASTA <- file.path(temp_dir, basename(path_to_FASTA))
  suppressMessages(
    path_to_new_FASTA <- remap_accessions_refseq_to_gene_fasta(
      path_to_FASTA, organism_name = "Rattus norvegicus"
    )
  )
  fst <- readAAStringSet(path_to_new_FASTA) # gene symbols

  expect_equal(head(names(fst)),
               c("A1bg", "A1cf", "A2m", "A3galt2", "A4galt", "A4gnt"))
  expect_equal(head(width(fst)),
               c(513, 594, 1472, 339, 360, 351))

  # Clean up
  temp_files <- list.files(temp_dir, full.names = TRUE)
  temp_files <- temp_files[!grepl("\\.fasta$", temp_files)]
  file.remove(temp_files)
  unlink(".Rcache", recursive = TRUE)
  unlink("../../.Rcache", recursive = TRUE)
})

