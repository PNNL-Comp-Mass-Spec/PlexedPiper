#' Finding Links from PNNL's DMS
#'
#' scanStats should have only one Dataset column
#' 
#' @param scanStats Scan stats
#' @param idScanPattern Default `ms2`
#' @param quantScanPattern Default `ms3`
#' @param threshold Default `NULL`
#' @return (data.table) with reporter ion intensities and other metrics
#' @importFrom dplyr select
#' @importFrom data.table data.table
#' @importFrom plyr llply
make_id_to_quant_scan_link_table <- function(scanStats, 
                                             idScanPattern = "ms2", 
                                             quantScanPattern = "ms3", 
                                             threshold = NULL) {
  # Select only relevant columns
  scanStats <- select(scanStats, 
                      Dataset, 
                      ScanNumber, 
                      ScanFilterText = `Scan Filter Text`)
  
  # If no quant pattern is found, warn user about used patterns
  if (!any(grepl(quantScanPattern, scanStats$ScanFilterText))){
    message <- paste("The quantification scan pattern '",
                     quantScanPattern,
                     "' was not found in the input. Please check that the right pattern is supplied",
                     sep = "")
    stop(message)
  }   
  
  # Making patterns more specific and able to capture the mass
  idScanPattern <- paste("^.*Full", idScanPattern, "([0-9.]*)@.*$")
  quantScanPattern <- paste("^.*Full", quantScanPattern, "([0-9.]*)@.*$")
  
  # Split scanStats into a named list according to dataset
  scanStatsList <- split(scanStats, by = "Dataset")
  linkingList <- list()
  
  pb = txtProgressBar(min = 0, max = length(scanStatsList), initial = 0,
                      style = 3)
  for (i in 1:length(scanStatsList)) {
    dataset_i <- names(scanStatsList)[[i]]
    
    # If quantification pattern is present in the dataset
    if (any(grepl(quantScanPattern, scanStatsList[[dataset_i]]))) {
      linkingList[[dataset_i]] <- 
        make_id_to_quant_scan_link_table_for_dataset(
          scanStatsList[[dataset_i]],
          quantScanPattern = quantScanPattern, 
          idScanPattern = idScanPattern, 
          threshold = threshold)
    } else {
      # If quantification pattern is NOT present in the dataset 
      # then quantification = identification
      linkingList[[dataset_i]] <- scanStatsList[[dataset_i]] %>% 
        select(Dataset, 
               IDScan = ScanNumber, 
               QuantScan = ScanNumber)
    }
    setTxtProgressBar(pb, i)
  }
  
  # Output as list instead?
  idToQuantScansLinkTable <- rbindlist(linkingList)
  return(idToQuantScansLinkTable)
}


# Works to link id scans to quant scans for a single dataset.
# comment treshold parameter
make_id_to_quant_scan_link_table_for_dataset <- function(scanStats, 
                                                         quantScanPattern, 
                                                         idScanPattern, 
                                                         threshold = NULL) {
  dataset <- unique(scanStats$Dataset)
  if (length(dataset) > 1){
    stop("Data from more than one dataset supplied, please check input.")
  }
  
  quantIdx <- grepl(quantScanPattern, scanStats$ScanFilterText)
  identificationIdx <- grepl(idScanPattern, scanStats$ScanFilterText)
  
  # Stop with an error whenever id and quant pattern matches coincide
  if (any(quantIdx & identificationIdx)){
    message <- paste("Quantification and identification patterns matched the same scan in dataset ",
                     dataset,
                     ". Please check further.",
                     sep = "")
    stop(message)
  }
  
  # The scans matching the quantification pattern, and extracting mass.
  quantScans <- scanStats %>% 
    filter(quantIdx) %>%
    mutate(ParentMass = sub(quantScanPattern, "\\1", ScanFilterText)) %>%
    rename(QuantScan = ScanNumber)
  
  # The scans matching the identification pattern, and extracting mass.
  idScans <- scanStats %>% 
    filter(identificationIdx) %>%
    mutate(ParentMass = sub(idScanPattern, "\\1", ScanFilterText)) %>%
    rename(IDScan = ScanNumber)
  
  # Identification and quantification scans are linked by IDMass
  idToQuantScansLink <- inner_join(idScans, quantScans, by = c("ParentMass")) %>% 
    select(IDScan, QuantScan)
  
  # idToQuantScansLink <- idToQuantScansLink[ , .(QuantScan = min(link_min_vectorized(IDScan, QuantScan))), 
  #                                           by = c("IDScan")]
  
  # We wish to link IDScan to the first scan matching its mass, with QuantScan > IDScan.
  # The inner_join above does not follow this requirement. We group by 
  # IDScan, then select the smallest QuantScan still larger than IDScan.
  # If for some IDScan no larger QuantScan is found, the grouping above 
  # gives infinity, so we filter those out.

  idToQuantScansLink <- idToQuantScansLink %>%
    mutate(d = abs(QuantScan - IDScan)) %>%
    group_by(IDScan) %>% 
    slice(which.min(d)) %>%
    ungroup()
  
  # Dynamic determination of threshold and filtering
  if (is.null(threshold)){
    threshold <- quantile(idToQuantScansLink$d, probs = 0.995)[[1]]
  }
  
  idToQuantScansLink <- idToQuantScansLink %>%
    filter(d < threshold) %>%
    mutate(Dataset = dataset) %>%
    select(Dataset, everything(), -d)
  
  return(idToQuantScansLink)
}

make_id_to_quant_scan_table_for_dataset_old <- function(scanStats, quantScanPattern, 
                                                    idScanPattern) {
  dataset <- unique(scanStats$Dataset)
  if (length(dataset) > 1){
    stop("Data from more than one dataset supplied, please check input.")
  }
  
  quantIdx <- grepl(quantScanPattern, scanStats$ScanFilterText)
  identificationIdx <- grepl(idScanPattern, scanStats$ScanFilterText)
  
  # Stop with an error whenever id and quant pattern matches coincide
  if (any(quantIdx & identificationIdx)){
    message <- paste("Quantification and identification patterns matched the same scan in dataset ",
                     dataset,
                     ". Please check further.",
                     sep = "")
    stop(message)
  }
  
  # The scans matching the quantification pattern, and extracting mass.
  quantScans <- scanStats %>% 
    filter(quantIdx) %>%
    mutate(IDMass = sub(quantScanPattern, "\\1", ScanFilterText)) %>%
    rename(QuantScan = ScanNumber)
  
  # The scans matching the identification pattern, and extracting mass.
  idScans <- scanStats %>% 
    filter(identificationIdx) %>%
    mutate(IDMass = sub(idScanPattern, "\\1", ScanFilterText)) %>%
    rename(IDScan = ScanNumber)
  
  # Identification and quantification scans are linked by IDMass
  idToQuantScansLink <- inner_join(idScans, quantScans, by = c("IDMass")) %>% 
    select(IDScan, QuantScan) 
  
  # We wish to link IDScan to the first scan matching its mass, with QuantScan > IDScan.
  # The inner_join above does not follow this requirement. We group by 
  # IDScan, then select the smallest QuantScan still larger than IDScan.
  idToQuantScansLink <- idToQuantScansLink[ , .(QuantScan = min(link_min_vectorized(IDScan, QuantScan))), 
                                            by = c("IDScan")]
  
  # If for some IDScan no larger QuantScan is found, the grouping above 
  # gives infinity, so we filter those out.
  idToQuantScansLink <- idToQuantScansLink %>%
    filter(!(QuantScan== Inf)) %>%  
    mutate(Dataset = dataset) %>%
    select(Dataset, everything())
  return(idToQuantScansLink)
}









