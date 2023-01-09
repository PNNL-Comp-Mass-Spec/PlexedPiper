# PlexedPiper 0.4.0 (2023-01-09)

-   Updated handling of UniProt-to-gene conversion in the `make_rii_peptide_*` and `make_results_ratio_*` functions. Now uses gene symbols from FASTA header and the *AnnotationDbi* Bioconductor package to convert to Entrez gene IDs, rather than converting from protein to gene symbol and Entrez gene.
-   Added "pkgdown" GitHub Action to automatically build website.

# PlexedPiper 0.3.7 (2022-10-22)

-   Updated MoTrPAC BIC functions to extract information from GENCODE FASTA headers and include them as columns in results tables.
-   Added checks for `create_crosstab` input.
-   Updated `run_plexedpiper` to accept multiple folder paths. Now capable of integrating multiple datasets, such as MoTrPAC PASS1A/1C.

# PlexedPiper 0.3.6 (2022-07-01)

-   Added use_mzIdentML argument to read_msgf_data to read mzid files from a local folder.
-   Added prefix argument to read_study_design and run_plexedpiper.
-   Expanded the functionality of filter_msgf_data to allow for filtering at the SiteID level. This performs the same optimization procedure as when filtering at the peptide level. Also added the filter_msgf_data_SiteID_level wrapper.

# PlexedPiper 0.3.5 (2022-05-04)

-   Added `unique_only` and `refine_prior` arguments to `run_plexedpiper`. Requires an update to a version of [MSnID](https://github.com/PNNL-Comp-Mass-Spec/MSnID) built on or after 2022-04-26--when the `refine_prior` argument was added to `infer_parsimonious_accessions`.
-   Minor documentation changes.

# PlexedPiper 0.3.4 (2022-04-15)

-   Added this *NEWS* file.
-   Added website built with the [pkgdown](https://pkgdown.r-lib.org/) package.

# PlexedPiper 0.3.3 (2022-04-14)

-   Added the `run_plexedpiper` wrapper function, which performs all the steps necessary to go from MS-GF+ and MASIC data to Reporter Ion Intensity (RII) and Ratio Results tables used by the MoTrPAC Bioinformatics Center (BIC).

# PlexedPiper 0.3.2 (2022-04-13)

-   Added TMT6 and TMT18 tables to `reporter_converter`.
-   Removed functions that had previously been migrated to [PNNL.DMS.utils](https://github.com/PNNL-Comp-Mass-Spec/PNNL.DMS.utils) during the switch to [v0.3.0](https://github.com/PNNL-Comp-Mass-Spec/PlexedPiper/releases/tag/0.3.0).
-   Reduced computation time and memory usage of `filter_msgf_data`. This has the added benefit of increasing the precision of the optimization algorithms. **As a result, filtered MSnIDs may be marginally different.**
-   Added UniProt support for MoTrPAC BIC functions.
-   Fixed check issues.

# PlexedPiper 0.3.1 (2021-09-22)

-   Removed odbc package from imports. Leftover from changes in [v0.3.0](https://github.com/PNNL-Comp-Mass-Spec/PlexedPiper/releases/tag/0.3.0).

# PlexedPiper 0.3.0 (2021-09-16)

-   Migrated DMS-related tools to [PNNL.DMS.utils](https://github.com/PNNL-Comp-Mass-Spec/PNNL.DMS.utils).

# PlexedPiper 0.2.0 (2021-06-07)

-   First release since migration from [vladpetyuk/PlexedPiper](https://github.com/vladpetyuk/PlexedPiper).
-   Updated `read_study_design`.
