## !!! The Dataset column of the ascore file is incorrect. The
## resulting MSnID has 0 rows because it cannot join ascore.

# testthat("best_PTM_location_by_ascore works", {
#   temp_dir <- tempdir()
#
#   # Get AScore results
#   path_to_ascore <- system.file("extdata/phospho/ascore_output",
#                                 package = "PlexedPiperTestData")
#   ascore <- read_AScore_results(path_to_ascore)
#
#   # MS/MS IDs
#   path_to_MSGF_results <- system.file("extdata/phospho/msgf_output",
#                                       package = "PlexedPiperTestData")
#   msgf_file <- list.files(path_to_MSGF_results, full.names = TRUE)[1:2]
#   file.copy(msgf_file, temp_dir)
#   msnid <- read_msgf_data(temp_dir)
#
#   file.remove(list.files(temp_dir, full.names = TRUE)) # clean temp dir
#
#   # Improve phosphosite localization
#   out <- best_PTM_location_by_ascore(msnid, ascore)
#   dim(psms(out))
# })

