#' Filtering MASIC Data
#'
#' Filtering MASIC data by interference score and signal-to-noise ratio (S/N) of
#' individual channels. Note: function in its current form also drops 40 columns
#' out of the `data.frame`. E.g. for TMT10 the table contains 62 columns. After
#' filtering, the remaining columns are the 10 TMT channels, `Dataset`, and
#' `ScanNumber`.
#'
#' Default values `interference_score_threshold = 0.9` and `s2n_threshold = 4`
#' are fairly stringent. We recommend at least setting
#' `interference_score_threshold = 0.5` and `s2n_threshold = 0`. If filtering is
#' not intended, use `filter_masic_data(x,0,0)`.
#'
#' @param x (data.frame) collated MASIC output
#' @param interference_score_threshold (numeric) The interference score reflects
#'   the proportion of the ion population isolated for fragmentation that is due
#'   to the targeted ion. `1 - InferferenceScore` is due to co-isolated species.
#'   The higher the interference score, the cleaner the parent ion at the MS1
#'   level. Default is 0.9.
#' @param s2n_threshold (numeric) S/N calculated by vendor and extracted MASIC
#'   from raw files. Default is 4.
#'
#' @return (data.frame) filtered MASIC output
#'
#' @md
#'
#' @importFrom dplyr filter select inner_join mutate %>% contains starts_with
#'   all_of
#' @importFrom tidyr pivot_longer pivot_wider
#'
#' @export filter_masic_data
#'
#' @examples
#' \dontrun{
#' path_to_MASIC_results <- system.file("extdata/global/masic_output",
#'                                      package = "PlexedPiperTestData")
#' x <- read_masic_data(path_to_MASIC_results, interference_score = TRUE)
#' dim(x)
#' x1 <- filter_masic_data(x,0,0)
#' dim(x1)
#' x2 <- filter_masic_data(x)
#' dim(x2)
#' }


filter_masic_data <- function(x,
                              interference_score_threshold = 0.9,
                              s2n_threshold = 4){
   reporter_ions <- colnames(x)[grepl("SignalToNoise", colnames(x))] %>%
      sub("_SignalToNoise", "", .)

   x <- x %>%
      filter(InterferenceScore >= interference_score_threshold) %>%
      select(Dataset, ScanNumber, all_of(reporter_ions),
             contains("SignalToNoise"))

   if (s2n_threshold == 0 | is.na(s2n_threshold)) {
      x <- x %>% select(-contains("SignalToNoise"))
   } else if (s2n_threshold > 0) {
      selected <- x %>%
         select(Dataset, ScanNumber, contains("SignalToNoise")) %>%
         # gather(channel, s2n, -c(Dataset, ScanNumber)) %>%
         pivot_longer(-c(Dataset, ScanNumber),
                      names_to = "channel", values_to = "s2n") %>%
         mutate(s2n = ifelse(is.na(s2n), 0, s2n)) %>%
         filter(s2n >= s2n_threshold) %>%
         select(-s2n) %>%
         mutate(channel = sub("_SignalToNoise", "", channel))

      x <- x %>% select(Dataset, ScanNumber, all_of(reporter_ions)) %>%
         # gather(channel, intensity, -c(Dataset, ScanNumber)) %>%
         pivot_longer(-c(Dataset, ScanNumber),
                      names_to = "channel", values_to = "intensity") %>%
         inner_join(selected, by = c("Dataset","ScanNumber", "channel")) %>%
         # spread(channel, intensity)
         pivot_wider(names_from = channel, values_from = intensity)
   }
}

utils::globalVariables(c(".", "InterferenceScore", "ScanNumber", "channel",
                         "s2n", "intensity"))

