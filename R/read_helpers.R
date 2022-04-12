# helper functions that work behind the scenes for both
# MASIC and MSGF output

#' @importFrom plyr llply
#' @importFrom dplyr bind_rows
#' @importFrom utils read.delim

read_with_name <- function(x, suffix){
   nm <- strsplit(basename(x), suffix)[[1]][1]
   x <- read.delim(x, header = T, stringsAsFactors = F)
   x <- x[,setdiff(colnames(x),"Dataset")]
   cbind.data.frame(Dataset = nm, x, stringsAsFactors=FALSE)
}

collate_files <- function(path_to_MASIC_results, suffix){
   files <- list.files(path_to_MASIC_results,
                       pattern = suffix,
                       full.names = TRUE)
   # using llply just for progress bar
   if(interactive()){

   }
   extra_metrics_list <- llply(files,
                               read_with_name,
                               suffix = suffix,
                               .progress = ifelse(interactive(),"text","none"))
   out <- bind_rows(extra_metrics_list)
}


