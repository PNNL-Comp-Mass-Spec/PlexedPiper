# PlexedPiper 0.3.5 (2022-04-28)

* Added `unique_only` and `refine_prior` arguments to `run_plexedpiper`. 
  Requires an update to a version of [MSnID](https://github.com/PNNL-Comp-Mass-Spec/MSnID) built on or after 2022-04-26--when the `refine_prior` argument was added to `infer_parsimonious_accessions`.
* Minor documentation changes.


# PlexedPiper 0.3.4 (2022-04-15)

* Added this *NEWS* file.
* Added website built with the [pkgdown](https://pkgdown.r-lib.org/) package.


# PlexedPiper 0.3.3 (2022-04-14)

* Added the `run_plexedpiper` wrapper function, which performs all the steps necessary to go from MS-GF+ and MASIC data to Reporter Ion Intensity (RII) and Ratio Results tables used by the MoTrPAC Bioinformatics Center (BIC).


# PlexedPiper 0.3.2 (2022-04-13)

* Added TMT6 and TMT18 tables to `reporter_converter`.
* Removed functions that had previously been migrated to [PNNL.DMS.utils](https://github.com/PNNL-Comp-Mass-Spec/PNNL.DMS.utils) during the switch to [v0.3.0](https://github.com/PNNL-Comp-Mass-Spec/PlexedPiper/releases/tag/0.3.0).
* Reduced computation time and memory usage of `filter_msgf_data`. This has the added benefit of increasing the precision of the optimization algorithms. **As a result, filtered MSnIDs may be marginally different.**
* Added UniProt support for MoTrPAC BIC functions.
* Fixed check issues.


# PlexedPiper 0.3.1 (2021-09-22)

* Removed odbc package from imports. Leftover from changes in [v0.3.0](https://github.com/PNNL-Comp-Mass-Spec/PlexedPiper/releases/tag/0.3.0).


# PlexedPiper 0.3.0 (2021-09-16)

* Migrated DMS-related tools to [PNNL.DMS.utils](https://github.com/PNNL-Comp-Mass-Spec/PNNL.DMS.utils).


# PlexedPiper 0.2.0 (2021-06-07)

* First release since migration from [vladpetyuk/PlexedPiper](https://github.com/vladpetyuk/PlexedPiper).
* Updated `read_study_design`.

