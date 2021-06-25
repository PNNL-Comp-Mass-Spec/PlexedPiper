## ----setup, echo=FALSE-----------------------------------------------------
knitr::opts_chunk$set(message=FALSE, warning=FALSE)
# knitr::opts_chunk$set(echo=T, message=F, warning=F, fig.align='center', out.width='10cm')

## ----libraries, eval=FALSE-------------------------------------------------
#  library(PlexedPiper)
#  library(PNNL.DMS.utils)
#  library(MSnID)
#  library(tidyverse)

## ----msms, eval=FALSE------------------------------------------------------
#  # phospho
#  msnid <- read_msms_data_from_DMS(3432)
#  ascore <- get_AScore_results(3432)
#  msnid <- best_PTM_location_by_ascore(msnid, ascore)
#  # at some point I need to filter non-phospho away
#  msnid <- filter_msgf_data_peptide_level(msnid, 0.01)
#  msnid <- remap_accessions_uniprot_to_gene(msnid, organism_name="Homo sapiens")
#  
#  path_to_FASTA <- path_to_FASTA_used_by_DMS(3432)
#  path_to_FASTA_gene <- remap_accessions_uniprot_to_gene_fasta(path_to_FASTA)
#  
#  # todo. or done? seems to be working
#  msnid <- compute_num_peptides_per_1000aa(msnid, path_to_FASTA_gene)
#  
#  
#  
#  show(msnid)
#  msnid <- filter_msgf_data_protein_level(msnid, 0.01)
#  show(msnid)
#  
#  # no need for parsimonius inference with phospho
#  # msnid <- infer_parsimonious_accessions(msnid, unique_only=TRUE)
#  # show(msnid)
#  
#  msnid <- apply_filter(msnid, "!isDecoy")
#  
#  
#  
#  # remap to site notation
#  # needed: 1) accession,
#  #         2) fasta with accession,
#  #         3) peptide sequence with asterisc location
#  # Let's borrow vp.misc remapping. !!! Ultimately it should be moved to MSnID
#  fst <- readAAStringSet(path_to_FASTA_gene, format="fasta",
#                         nrec=-1L, skip=0L, use.names=TRUE)
#  ids <- psms(msnid) %>%
#     distinct(accession, peptide)
#  ids_with_sites <- map_PTM_sites(ids, fst, "accession", "peptide", "*")
#  
#  ids_with_sites <- ids_with_sites %>%
#     mutate(idx = map(PepLoc, length) %>% unlist) %>%
#     filter(idx != 0) %>%
#     dplyr::select(-idx)
#  
#  ids_with_sites_simplified <- ids_with_sites %>%
#     dplyr::select(accession, peptide, SiteCollapsedFirst) %>%
#     mutate(site = unlist(SiteCollapsedFirst) %>%
#               gsub(",","_",.) %>%
#               paste(accession, ., sep="-")) %>%
#     dplyr::select(-SiteCollapsedFirst)
#  
#  psms(msnid) <- inner_join(psms(msnid), ids_with_sites_simplified)
#  
#  psms(msnid) <- psms(msnid) %>%
#     inner_join(ids_with_sites_simplified)
#  
#  
#  
#  
#  
#  
#  masic_data <- read_masic_data_from_DMS(3432, interference_score = T)
#  masic_data <- filter_masic_data(masic_data, 0.5, 0)
#  
#  library(dplyr)
#  fractions <- masic_data %>%
#     distinct(Dataset) %>%
#     mutate(PlexID = sub("SHSY5Y_IIS_3x3_P_(B\\d).*","\\1",Dataset))
#  head(fractions)
#  
#  library(readr)
#  samples <- read_tsv("../../3432/samples.txt")
#  head(samples,10)
#  
#  ref <- "(`10min_control1`*`10min_control2`*`10min_Insulin`*`10min_IGF1`*`10min_IGF2`*`60min_control1`*`60min_control2`*`60min_Insulin`*`60min_IGF1`*`60min_IGF2`)^(1/10)"
#  references <- samples %>%
#     distinct(PlexID, QuantBlock) %>%
#     mutate(Reference = ref)
#  
#  
#  aggregation_level <- c("site")
#  quant_cross_tab <- create_crosstab(msnid,
#                                     masic_data,
#                                     aggregation_level,
#                                     fractions, samples, references)
#  dim(quant_cross_tab)
#  head(quant_cross_tab)
#  
#  
#  

