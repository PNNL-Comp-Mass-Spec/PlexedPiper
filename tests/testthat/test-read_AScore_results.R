test_that("read_AScore_results works", {
  path_to_AScore_results <- system.file("extdata/phospho/ascore_output",
                                        package = "PlexedPiperTestData")
  ascore <- read_AScore_results(path_to_AScore_results)
  expect_equal(nrow(ascore), 1177348)
})

