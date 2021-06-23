## ----setup, echo=FALSE-----------------------------------------------------
knitr::opts_chunk$set(message=FALSE, warning=FALSE)
# knitr::opts_chunk$set(echo=T, message=F, warning=F, fig.align='center', out.width='10cm')

## ----test_connection-------------------------------------------------------
# if there is no connection to PNNL DMS, vignette is not compiled
if(!is_PNNL_DMS_connection_successful()){
   message("There is no connection to PNNL DMS. This code in this vignette won't be evaluated.")
   knitr::opts_chunk$set(eval = FALSE)
}

## ----libraries-------------------------------------------------------------
library(PlexedPiper)
library(PNNL.DMS.utils)
library(MSnID)

## ----msms------------------------------------------------------------------
msnid <- read_msms_data_from_DMS(3442)
show(msnid)

## ----correct_isotopes------------------------------------------------------
msnid <- correct_peak_selection(msnid)

## ----peptide_level_filter--------------------------------------------------
show(msnid)
msnid <- filter_msgf_data_peptide_level(msnid, 0.01)
show(msnid)

## ----refseq_2_gene---------------------------------------------------------
msnid <- remap_accessions_refseq_to_gene(msnid, 
                                         organism_name="Rattus norvegicus")

## ----fasta_refseq_2_gene---------------------------------------------------
path_to_FASTA <- path_to_FASTA_used_by_DMS(3442)
path_to_FASTA_gene <- remap_accessions_refseq_to_gene_fasta(
   path_to_FASTA,
   organism_name = "Rattus norvegicus")

## ----protein_level_filter--------------------------------------------------
msnid <- compute_num_peptides_per_1000aa(msnid, path_to_FASTA_gene)
show(msnid)
msnid <- filter_msgf_data_protein_level(msnid, 0.01)
show(msnid)

## ----inference-------------------------------------------------------------
msnid <- infer_parsimonious_accessions(msnid, unique_only=TRUE)
show(msnid)

## ----remove_decoy----------------------------------------------------------
msnid <- apply_filter(msnid, "!isDecoy")

## ----read_masic------------------------------------------------------------
masic_data <- read_masic_data_from_DMS(3442, interference_score=TRUE)
head(masic_data)

## ----filter_masic----------------------------------------------------------
nrow(masic_data)
masic_data <- filter_masic_data(masic_data, 0.5, 0)
nrow(masic_data)

## ----read_fractions--------------------------------------------------------
library(dplyr)
fractions <- masic_data %>%
   distinct(Dataset) %>%
   mutate(PlexID =sub("MoTrPAC_Pilot_TMT_W_(S\\d).*","\\1",Dataset))
head(fractions)

## ----read_samples----------------------------------------------------------
library(readr)
samples <- read_tsv(system.file("extdata/study_design/samples.txt", package = "PlexedPiperTestData"))
head(samples,10)

## ----read_references-------------------------------------------------------
references <- read_tsv(system.file("extdata/study_design/references.txt", package = "PlexedPiperTestData"))
head(references)

## ----crosstab_compillation-------------------------------------------------
aggregation_level <- c("accession")
quant_cross_tab <- create_crosstab(msnid, 
                                   masic_data, 
                                   aggregation_level, 
                                   fractions, samples, references)
dim(quant_cross_tab)
head(quant_cross_tab)

## ----cleanup, echo=FALSE---------------------------------------------------
unlink(".Rcache", recursive=TRUE)

