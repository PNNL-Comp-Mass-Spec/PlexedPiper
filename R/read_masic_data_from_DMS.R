#' Reading MASIC results from PNNL's DMS
#'
#' @param data_package_num (Numeric or Character vector) containing Data Package ID(s) located in DMS
#' @param interference_score (logical) read interference score. Default is FALSE.
#' @param id_quant_table (logical) 
#' @return (data.frame) with reporter ion intensities and other metrics
#' @importFrom dplyr select
#' @importFrom plyr llply
#' @importFrom data.table rbindlist
#' @export read_masic_data_from_DMS

read_masic_data_from_DMS <- function(data_package_num, interference_score=FALSE, 
                                     id_quant_table = FALSE){
   on.exit(gc(), add = TRUE)

   # Fetch job records for data package(s)
   if(length(data_package_num) > 1){
      jobRecords <- lapply(data_package_num, get_job_records_by_dataset_package)
      jobRecords <- Reduce(rbind, jobRecords)
   }else{
      jobRecords <- get_job_records_by_dataset_package(data_package_num)
   }

   jobRecords <- jobRecords[grepl("MASIC", jobRecords$Tool),]

   # select relevant columns from df, drop redundant dataset column (select(-2))
   masicList <- get_results_for_multiple_jobs.dt(jobRecords)
   masicData <- masicList[["_ReporterIons.txt"]] %>% 
      select(-2) %>% 
      select(Dataset, 
             ScanNumber, 
             starts_with("Ion"), 
             -contains("Resolution"))
   
   if (id_quant_table) {
      quantIDTable <- make_id_to_quant_scan_link_table(masicList[["_ScanStatsEx.txt"]][,-2]) %>%
         rename(ScanNumber = QuantScan)
      masicData <- inner_join(quantIDTable, masicData, 
                         by = c("Dataset", "ScanNumber")) %>%    
         select(-ScanNumber) %>%    
         rename(ScanNumber = IDScan)
   }

   if (interference_score) {
      masicStats <- masicList[["_SICstats.txt"]] %>% 
         select(-2) %>%
         select(Dataset, 
                ScanNumber = FragScanNumber, 
                contains("InterferenceScore"))
      masicData <- inner_join(masicData, masicStats, 
                              by = c("Dataset", "ScanNumber"))
   }
   return(masicData)
}


# fetch_masic_data_for_single_datset <- function(pathToFile, fileNamePttrn ){
#   pathToFile = list.files( path=as.character(pathToFile),
#                            pattern=fileNamePttrn,
#                            full.names=T)
#   if(length(pathToFile) == 0){
#     stop("can't find the results file")
#   }
#   if(length(pathToFile) > 1){
#     stop("ambiguous results files")
#   }
#   results = read.delim( pathToFile, header=T, stringsAsFactors = FALSE)
#   dataset = strsplit( basename(pathToFile), split=fileNamePttrn)[[1]]
#   out = data.table(Dataset=dataset, results)
#   return(out)
# }





#
#
#
#
# library("devtools")
# source_url("https://raw.githubusercontent.com/vladpetyuk/PNNL_misc/master/PNNL_DMS_utils.R")
#
#
# # 1
# jobRecords <- lapply(c(2930,2988,3087), get_job_records_by_dataset_package)
# jobRecords <- Reduce(rbind, jobRecords)
# jobRecords <- subset(jobRecords, Tool == "MASIC_Finnigan")
#
#
#
# system.time({
#     masicData <- get_results_for_multiple_jobs.dt(jobRecords)
# })
#
#
# library(plyr)
# library(data.table)
# results = llply( jobRecords[["Folder"]],
#                  get_results_for_single_job.dt,
#                  fileNamePattern="_SICstats.txt",
#                  .progress = "text")
# results.dt <- rbindlist(results)
# masicStats <- as.data.frame(results.dt)
#
# # hack to remove redundant Dataset column
# masicData  <- masicData[,-2] # use make.names and then remove Dataset.1
# masicStats <- masicStats[,-2]
# colnames(masicStats)[colnames(masicStats) == "FragScanNumber"] <- "ScanNumber" # dplyr::rename
#
# library(dplyr)
#
# x <- select(masicData, Dataset, ScanNumber, starts_with("Ion"), -contains("Resolution"))
# y <- select(masicStats, Dataset, ScanNumber, InterferenceScore)
# z <- inner_join(x, y)
# masicData <- z
# save(masicData, file="masicData_original.RData")
#
#
#
#
#
#
