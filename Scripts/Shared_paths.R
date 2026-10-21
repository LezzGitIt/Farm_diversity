# Shared input paths for the pipeline ----

### Sourced by every script that reads the data-paper tables or the Scripts/00 diversity exports (00, 01a, 02a, 02b, 04a, 04b). No package dependencies. The modelling helpers are separate, in Scripts/Shared_fns.R.

## The bird data-paper tables (Bird_pcs/Bird_pcs_analysis.csv, Site_covs.csv, Event_covs.csv), read from the sibling data-wrangling repo. Column names are the data paper's English ones (Id_scr = SCR farm ID, Id_survey = point-count survey, Year_grp, Date, ...).
WRANGLING_EXCELS <- "../Ssp-bird-data-wrangling/Derived/Excels"

## The Scripts/00 diversity-estimate export every downstream script reads, pinned by its date stamp (Tax_div_<type>_<version>.csv). Pinned rather than "newest file" so that re-running 00 does not silently change the inputs to the cached brms fits; bump this deliberately after a re-run of 00.
TAX_DIV_VERSION <- "10.07.26"

## Path of one Scripts/00 export -- `type` is "incidence", "all_farms" or "coverage65"
tax_div_path <- function(type, version = TAX_DIV_VERSION) {
  file.path("Derived/Excels", sprintf("Tax_div_%s_%s.csv", type, version))
}
