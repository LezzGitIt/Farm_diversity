# Bird taxonomic-diversity (Hill number) estimates per assemblage ----

### Ported from `Scripts/qmd/02_Analysis_iNEXT.qmd`, which is kept as the reference write-up (Cmax rationale, the iNEXT 4-steps framework, the sample-completeness / rarefaction plots). This is the runnable pipeline version -- data prep, the iNEXT4steps call, and the two exports, nothing else.

### An assemblage = [data collector . farm . year-group . season]. For each assemblage the raw point-count counts are summed (across days, repeat surveys and point counts) into a species-by-assemblage matrix, then Anne Chao's iNEXT 4-steps workflow (Chao et al. 2020) gives asymptotic (q = 1, 2) and coverage-standardised non-asymptotic (q = 0, at Cmax) diversity with bootstrap uncertainty.

### Two passes, in sequence:
###   1. all assemblages            -> Derived/Excels/Tax_div_all_farms_<date>.csv   (used for q = 1 / q = 2 asymptotic)
###   2. drop the low-coverage ones -> Derived/Excels/Tax_div_coverage65_<date>.csv  (used for q = 0 non-asymptotic; raises Cmax)
### The pass-2 removal set is taken from pass 1's own SC_obs (assemblages with observed sample coverage < 0.65 at q = 0) -- self-contained, no external file.

# Setup ----
library(tidyverse)
library(iNEXT)
library(iNEXT.4steps)
library(conflicted)

conflicts_prefer(dplyr::select, dplyr::filter, .quiet = TRUE)

wrangling_excels <- "../Ssp-bird-data-wrangling/Derived/Excels"
dir.create("Derived/Excels", recursive = TRUE, showWarnings = FALSE)

## Row identifier for an assemblage
row_id <- c("Uniq_db", "Id_gcs", "Ano_grp", "Season")

## iNEXT bootstrap replicates -- 500 so the per-assemblage SEs are stable enough to use as measurement error downstream (Scripts/04_farm_mgmt_mod.R)
nboot <- 500

sc_threshold <- 0.65
today <- format(Sys.Date(), "%m.%d.%y")

# Load and link the point-count data ----

Bird_pcs <- read_csv(file.path(wrangling_excels, "Bird_pcs/Bird_pcs_analysis.csv"),
                     show_col_types = FALSE)
Site_covs <- read_csv(file.path(wrangling_excels, "Site_covs.csv"), show_col_types = FALSE)
Event_covs <- read_csv(file.path(wrangling_excels, "Event_covs.csv"), show_col_types = FALSE)

## One row per point-count survey with the assemblage identifier attached
Event_covs2 <- Site_covs %>%
  select(Id_muestreo_no_dc, Id_gcs, Habitat) %>%
  left_join(
    Event_covs %>% select(Id_muestreo_no_dc, Id_muestreo, Uniq_db, Fecha, Pc_start,
                          Ano_grp, Season, Rep_ano_grp),
    by = "Id_muestreo_no_dc"
  ) %>%
  mutate(
    Identifier = pmap_chr(across(all_of(row_id)), ~ paste(..., sep = ".")),
    Identifier = str_replace_all(Identifier, " |-", "_")
  )

Bird_pcs_linked <- Bird_pcs %>%
  left_join(Event_covs2, by = c("Id_muestreo", "Id_muestreo_no_dc", "Fecha", "Pc_start")) %>%
  mutate(Id_muestreo = str_replace_all(Id_muestreo, "-", "_"))

## Number of distinct habitat types per assemblage -- a sampling-design covariate carried in the export
num_hab <- Event_covs2 %>%
  filter(!is.na(Id_gcs)) %>%
  distinct(Identifier, Id_muestreo, Habitat) %>%
  summarize(Num.hab = as.factor(n_distinct(Habitat)), .by = Identifier)

# Species-by-assemblage matrix ----

## Sum counts within [assemblage x day x repeat], then across days, then widen to species x assemblage. Summing (not averaging) keeps the abundances iNEXT needs consistent with the accumulating species pool.
build_matrix <- function(bird_data) {
  bird_data %>%
    filter(!is.na(Id_gcs), Season != "Late") %>%
    mutate(Identifier = pmap_chr(across(all_of(row_id)), ~ paste(..., sep = ".")),
           Identifier = str_replace_all(Identifier, " |-", "_")) %>%
    summarize(Count = sum(Count),
              .by = c(Identifier, Species_ayerbe, Fecha, Rep_ano_grp)) %>%
    summarize(Count = sum(Count), .by = c(Identifier, Species_ayerbe)) %>%
    pivot_wider(names_from = Identifier, values_from = Count, values_fill = 0) %>%
    rename_with(~ str_replace_all(.x, " |-", "_")) %>%
    column_to_rownames("Species_ayerbe")
}

# One iNEXT 4-steps pass -> tidy per-assemblage summary ----

run_inext_pass <- function(bird_data) {
  bw <- build_matrix(bird_data)
  message("  iNEXT4steps on ", ncol(bw), " assemblages x ", nrow(bw), " species (nboot = ", nboot, ") ...")
  ## suppressWarnings: iNEXT4steps builds ggplot objects internally and emits cosmetic shape-palette / dropped-row warnings for > 6 assemblages
  i4 <- suppressWarnings(iNEXT4steps(data = bw, nboot = nboot, details = TRUE, datatype = "abundance"))

  ## Cmax: the standardised coverage STEP 3 compares assemblages at
  cmax <- i4$summary[["STEP 3. Non-asymptotic coverage-based rarefaction and extrapolation analysis"]] %>%
    colnames() %>% pluck(1) %>% str_split_i(" ", 3) %>% as.numeric() %>% round(2)

  ## STEP 1 / 3 / 4: observed sample completeness, non-asymptotic diversity at Cmax, evenness -- long by Order.q
  long_bits <- map2(
    i4$summary[c(1, 3, 4)], c("SC_obs", "No_Asy_TD", "Evenness"),
    function(df, value_col) {
      df %>%
        rename_with(~ "Assemblage", .cols = 1) %>%
        pivot_longer(2:4, names_to = "Order.q", values_to = value_col) %>%
        mutate(Order.q = str_remove(Order.q, "q = "))
    }
  )

  ## STEP 2: observed + asymptotic diversity with bootstrap SE / CI (used for q = 1, 2)
  summary_df <- i4$summary[[2]] %>%
    mutate(Order.q = as.character(case_when(
      qTD == "Species richness"  ~ 0,
      qTD == "Shannon diversity" ~ 1,
      qTD == "Simpson diversity" ~ 2
    ))) %>%
    left_join(long_bits[[1]], by = c("Assemblage", "Order.q")) %>%
    left_join(long_bits[[2]], by = c("Assemblage", "Order.q")) %>%
    full_join(long_bits[[3]], by = c("Assemblage", "Order.q")) %>%
    left_join(num_hab, by = c("Assemblage" = "Identifier"))

  ## Bootstrap CI of the coverage-standardised (non-asymptotic) estimate at Cmax. iNEXT4steps computes the coverage-based rarefaction/extrapolation curve with CIs (`details$iNEXT_coverage_based`) but STEP 3 surfaces only the point estimate, so interpolate the curve's CI bounds at exactly SC = Cmax. SE from the CI half-width. If iNEXT clamps qTD.LCL at 0 the SE is a slight under-estimate -- only for a few very sparse assemblages.
  cov_ci <- purrr::list_flatten(i4$details)$iNEXT_coverage_based %>%
    as_tibble() %>%
    mutate(Order.q = as.character(Order.q)) %>%
    summarize(
      No_Asy_TD_LCL = approx(SC, qTD.LCL, xout = cmax, rule = 2)$y,
      No_Asy_TD_UCL = approx(SC, qTD.UCL, xout = cmax, rule = 2)$y,
      .No_Asy_TD_check = approx(SC, qTD, xout = cmax, rule = 2)$y,
      .by = c(Assemblage, Order.q)
    ) %>%
    mutate(No_Asy_TD_se = (No_Asy_TD_UCL - No_Asy_TD_LCL) / (2 * qnorm(0.975)))

  out <- summary_df %>%
    left_join(cov_ci, by = c("Assemblage", "Order.q"))

  discrepancy <- with(out, max(abs(No_Asy_TD - .No_Asy_TD_check), na.rm = TRUE))
  if (discrepancy > 0.5) warning("coverage-based CI matched to a point estimate ", round(discrepancy, 2), " off No_Asy_TD")

  out %>%
    select(-.No_Asy_TD_check) %>%
    separate_wider_delim(Assemblage, delim = ".", names = row_id,
                         too_few = "align_start", cols_remove = FALSE) %>%
    relocate(Order.q, .after = qTD)
}

# Pass 1: all assemblages ----

message("Pass 1 (all assemblages)")
td_all <- run_inext_pass(Bird_pcs_linked)
write_csv(td_all, sprintf("Derived/Excels/Tax_div_all_farms_%s.csv", today))
cat("Wrote Tax_div_all_farms_", today, ".csv: ",
    n_distinct(td_all$Assemblage), " assemblages.\n", sep = "")

# Pass 2: drop the low-coverage assemblages, re-estimate ----

low_sc <- td_all %>%
  filter(Order.q == "0", SC_obs < sc_threshold) %>%
  distinct(Assemblage)

cat("\n", nrow(low_sc), " assemblages with observed sample coverage < ", sc_threshold,
    " at q = 0 -- removed for pass 2:\n", sep = "")
print(sort(low_sc$Assemblage))

Bird_pcs_high_sc <- Bird_pcs_linked %>%
  mutate(Assemblage = pmap_chr(across(all_of(row_id)), ~ paste(..., sep = ".")),
         Assemblage = str_replace_all(Assemblage, " |-", "_")) %>%
  anti_join(low_sc, by = "Assemblage") %>%
  select(-Assemblage)

message("\nPass 2 (low-coverage assemblages removed)")
td_c65 <- run_inext_pass(Bird_pcs_high_sc)
write_csv(td_c65, sprintf("Derived/Excels/Tax_div_coverage65_%s.csv", today))
cat("Wrote Tax_div_coverage65_", today, ".csv: ",
    n_distinct(td_c65$Assemblage), " assemblages.\n", sep = "")
