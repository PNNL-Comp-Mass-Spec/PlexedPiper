# test_that("remap_uniprot_to_gene_fasta works", {
#   path_to_FASTA <- system.file("extdata", "uniprot_rat_small.fasta",
#                                package = "MSnID")
#   temp_dir <- tempdir()
#   file.copy(path_to_FASTA, temp_dir)
#   path_to_FASTA <- file.path(temp_dir, basename(path_to_FASTA))
#
#   # Create new FASTA file and output file path
#   path_to_FASTA_gene <- remap_accessions_uniprot_to_gene_fasta(
#     path_to_FASTA
#   )
#
#   fst_gene <- readAAStringSet(path_to_FASTA_gene)
#
#   expect_equal(width(fst_gene[c("A1bg", "A1cf", "A1i3", "A1m", "A2m",
#                                 "A3galt2", "A4galt", "Aacs",
#                                 "Aadac", "Aadat")]),
#                c(513, 594, 1477, 1500, 1472, 339, 360, 672, 398, 425))
#   expect_equal(length(fst_gene), 7986)
#
#   # Clean up
#   temp_files <- list.files(temp_dir, full.names = TRUE)
#   temp_files <- temp_files[!grepl("\\.fasta$", temp_files)]
#   file.remove(temp_files)
#   unlink(".Rcache", recursive = TRUE)
#   unlink("../../.Rcache", recursive = TRUE)
# })

