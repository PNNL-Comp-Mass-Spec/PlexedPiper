#' Filtering MASIC Data
#'
#' Filtering MASIC data by InterferenceScore and S/N of individual channels.
#' Note, function in its current form also drops 40 columns out of the data.frame.
#' E.g. for TMT10 the table contains 62 columns. After filtering, the left columns
#' are 10 TMT channels, dataset and scan.
#'
#' Default values interference_score_threshold = 0.9 and s2n_threshold = 4 are
#' fairly stringent. The least recomended setting is interference_score_threshold = 0.5
#' and s2n_threshold = 0. If not filtering intended do filter_masic_data(x,0,0) call.
#'
#' @param x (data.frame) collated MASIC output
#' @param interference_score_threshold (numeric) fetch extra metrics that MASIC extracts from dataset or not. Higher the number, the cleaner parent ion at MS1 level. Default is 0.9.
#' @param s2n_threshold (numeric) S/N calculated by vendor and extracted MASIC from raw files. Default is 4.
#' @return (data.frame) filtered MASIC output
#' @importFrom dplyr filter select inner_join mutate %>%
#' @importFrom tidyselect contains starts_with all_of
#' @importFrom tidyr gather spread
#' @export filter_masic_data
#' @examples \dontrun{
#' path_to_MASIC_results <- system.file("extdata/global/masic_output", 
#' package = "PlexedPiperTestData")
#' x <- read_masic_data(path_to_MASIC_results, extra_metrics=TRUE)
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

   x <- x %>% filter(InterferenceScore >= interference_score_threshold) %>%
      select(Dataset, ScanNumber, all_of(reporter_ions), contains("SignalToNoise"))

   if (s2n_threshold == 0 | is.na(s2n_threshold)) {
      x <- x %>% select(-contains("SignalToNoise"))
   } else if (s2n_threshold > 0) {
      selected <- x %>% select(Dataset, ScanNumber, contains("SignalToNoise")) %>%
         gather(channel, s2n, -c(Dataset, ScanNumber)) %>%
         mutate(s2n = ifelse(is.na(s2n), 0, s2n)) %>%
         filter(s2n >= s2n_threshold) %>% select(-s2n) %>%
         mutate(channel = sub("_SignalToNoise", "", channel))

      x <- x %>% select(Dataset, ScanNumber, all_of(reporter_ions)) %>%
         gather(channel, intensity, -c(Dataset, ScanNumber)) %>%
         inner_join(selected, by = c("Dataset","ScanNumber", "channel")) %>%
         spread(channel, intensity)
   }
}


#
# load("masicData_original.RData")
#
# library(dplyr)
# library(tidyr)
#
# s2n_min_threshold = 4
# interference_score_threshold = 0.90
#
# # retain only those MS2 scans that do not have much interference at MS1 levels
# masicData <- masicData %>%
#    filter(InterferenceScore > interference_score_threshold)
#
# selected <- masicData %>%
#    select(Dataset, ScanNumber, contains("SignalToNoise")) %>%
#    gather(channel, s2n, -c(Dataset, ScanNumber)) %>%
#    mutate(s2n = ifelse(is.na(s2n), 0, s2n)) %>% # impute NA with 0 values
#    # group_by(Dataset, ScanNumber) %>%
#    # summarise(s2n_min = min(s2n)) %>%
#    filter(s2n > s2n_min_threshold) %>% # key filter step
#    select(-s2n) %>%
#    mutate(channel = sub("_SignalToNoise","",channel))
# # QC viz
# selected %>% group_by(Dataset, ScanNumber) %>% summarise(n = n()) %>% .$n %>% table %>% barplot
#
# # linked filtered Dataset/ScanNumber/channel
# masicData <- masicData %>%
#    select(Dataset, ScanNumber, starts_with("Ion"), -contains("SignalToNoise")) %>%
#    gather(channel,intensity,-c(Dataset,ScanNumber)) %>%
#    inner_join(selected) %>%
#    spread(channel, intensity)
#
# save(masicData, file="masicData_filtered.RData", compress = 'xz')
#
