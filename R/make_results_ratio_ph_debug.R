

make_results_ratio_ph <- function(msnid, masic_data, fractions, samples, 
                                  references, org_name = "Rattus norvegicus", sep="_") {
  
  aggregation_level <- c("accession", "SiteID")
  crosstab <- create_crosstab(msnid, masic_data, aggregation_level, fractions,
                              samples, references)
  crosstab <- as.data.frame(crosstab) %>% 
    rownames_to_column('Specie')
  
  ## Create RII peptide table
  results_ratio <- crosstab %>%
    select(Specie) %>%
    mutate(protein_id = sub("(^.*)@(.*)", "\\1", Specie),
           ptm_id = sub("(^.*)@(.*)", "\\2", Specie),
           organism_name = org_name) %>%
    mutate(REFSEQ = sub("(^.*)\\.\\d+", "\\1", protein_id))
  
  ## Add Genes + EntrezID
  conv <- fetch_conversion_table(org_name, from = "REFSEQ", "SYMBOL")
  
  conv <- fetch_conversion_table(org_name, from = "REFSEQ", "ENTREZID") %>%
    inner_join(., conv)
  
  results_ratio <- left_join(results_ratio, conv) %>%
    rename(gene_symbol = SYMBOL,
           entrez_id = ENTREZID) %>%
    select(-REFSEQ)
  
  ## Additional info from MS/MS
  ids <- psms(msnid) %>%
    select(accession, peptide, SiteID,
           noninferableProteins, flankingSequence,
           MSGFDB_SpecEValue, maxAScore) %>%
    rename(protein_id = accession,
           sequence = peptide,
           ptm_id = SiteID,
           redundant_ids = noninferableProteins,
           flanking_sequence = flankingSequence) %>%
    # group at peptide level to calculate peptide score, confident score
    group_by(protein_id, sequence, ptm_id, flanking_sequence, redundant_ids) %>%
    summarize(peptide_score = min(MSGFDB_SpecEValue),
              confident_score = max(maxAScore)) %>%
    # regroup at siteID level and recalculate ptm score
    group_by(protein_id, ptm_id, flanking_sequence, redundant_ids) %>%
    summarize(ptm_score = min(peptide_score),
              confident_score = max(confident_score)) %>%
    mutate(confident_site = case_when(confident_score >= 17 ~ TRUE,
                                      confident_score < 17 ~ FALSE),
           is_contaminant = grepl("Contaminant", protein_id))
  
  results_ratio <- inner_join(results_ratio, ids) %>%
    mutate(ptm_id = gsub("-", sep, ptm_id))
  
  ## Join with crosstab
  results_ratio <- inner_join(results_ratio, crosstab) %>%
    select(-Specie)
  
  return(results_ratio)
}