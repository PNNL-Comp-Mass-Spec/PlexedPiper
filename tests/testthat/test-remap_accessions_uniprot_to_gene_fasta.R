test_that("remap_uniprot_to_gene_fasta works", {
  path_to_FASTA <- system.file("extdata", "uniprot_rat_small.fasta.gz",
                               package = "MSnID")
  temp_dir <- tempdir()
  file.copy(path_to_FASTA, temp_dir)
  path_to_FASTA <- file.path(temp_dir, basename(path_to_FASTA))

  # Create new FASTA file and output file path
  path_to_FASTA_gene <- remap_accessions_uniprot_to_gene_fasta(
    path_to_FASTA
  )

  fst_gene <- readAAStringSet(path_to_FASTA_gene)

  expect_equal(width(fst_gene[c("Aatf", "Acsl4", "Alkbh5", "Andpro",
                                "Angpt1", "Atp2b2", "B4galt1", "Calml3",
                                "Ccnf", "Cct3")]),
               c(523, 711, 395, 176, 497, 1205, 32, 149, 780, 545))
  expect_equal(length(fst_gene), 97)

  # Clean up
  try(expr = {
    temp_files <- list.files(temp_dir, full.names = TRUE)
    temp_files <- temp_files[!grepl("\\.fasta$", temp_files)]
    file.remove(temp_files)
    unlink(".Rcache", recursive = TRUE)
    unlink("../../.Rcache", recursive = TRUE)
  })
})

