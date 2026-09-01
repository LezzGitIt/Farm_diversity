# Bird taxonomic-diversity (Hill number) estimates per assemblage ----

### An assemblage = [data collector . farm . year-group . season]. This script produces TWO diversity estimates per assemblage, from the same point-count data:

### ESTIMATE A -- point-count-standardised (iNEXT incidence). Treated as equal evidence to Estimate B, not primary/sensitivity (Aaron 2026-08-31); one goes in the main text, one in the supplement.
###   Sampling unit = a point-count LOCATION (`Id_muestreo`); a species is scored present (1) if it was detected there on any visit. Diversity (Hill q = 0/1/2) is then rarefied / extrapolated to a common NUMBER OF POINT COUNTS, m* (`estimateD(base = "size", level = m*)`), i.e. "expected diversity if this assemblage had been surveyed with m* point counts". This standardises SPATIAL survey effort directly, so the estimate carries essentially no `Num_pc` signal (effort R^2 ~ 0.02 vs ~0.34 for the coverage-based richness) and Scripts/04a does not need `Num_pc` as a covariate. m* = 6: iNEXT's own `min(2 x reference size)` rule once assemblages with < `min_pc` = 3 point counts are dropped (keeps ~102 of 112), and also the median effort.
###   CAVEATS (carry these into any write-up):
###     * INCIDENCE needs 0/1 cell values -- so with the LOCATION as the sampling unit, repeat visits are POOLED to presence/absence: a species detected once and one detected on every visit both score 1. Within-point-count and within-visit abundance is discarded, so q = 1 / q = 2 are INCIDENCE-based Hill numbers (built from among-point-count occupancy frequencies), not abundance-weighted diversity. (Pooling is a consequence of the unit choice, not a requirement of incidence -- visit-as-unit keeps visit-level info, see next caveat.)
###     * ALL POINT COUNTS ARE TREATED AS EQUIVALENT UNITS regardless of protocol. The 6 datasets differ in point-count DURATION (CIPAV / Gaica-mbd variable, up to 90 min; the rest 10 min), RADIUS (25 vs 50 m), and NUMBER OF VISITS (CIPAV 1; the rest 3-5). Standardising the NUMBER of point counts does not standardise the effort PER point count -- a longer / wider / more-revisited count detects more species, so its incidence column is fuller and the assemblage reads higher. Choosing the visit as the unit would let iNEXT's incidence model use detection frequency, but then the incidence sample size = locations x visits, which is not comparable across the single-visit vs multi-visit datasets. Between-dataset protocol differences are therefore left to `(1 | CollectorXyear)` in Scripts/04a, NOT removed at estimation. (The abundance / Cmax estimate handles per-count duration/radius somewhat better -- individuals / coverage is a finer effort currency than "number of counts" -- but leaves the number-of-counts residual. Each approach trades one confound for another; hence report both.)
###     * NO coverage standardisation. Assemblages are equalised on effort (6 point counts), not on completeness -- sample coverage at m* is ~0.78 median (0.36-0.99). A species-rich assemblage sits further from its asymptote at 6 point counts, so the estimate partly reflects detectability / density, and it reads LOWER for richly-sampled farms than the near-asymptotic coverage-based estimate.
###     * ~10 assemblages with < 3 point counts are dropped (the coverage-based export keeps them, with wide CIs). ~1/3 of the rest (3-5 point counts) are extrapolated to 6 -- within the 2x limit, but model-dependent.
###     * m* is a researcher choice. Lower m* keeps more assemblages but sits further from the asymptote.

### ESTIMATE B -- abundance, coverage-standardised (the original approach; Chao et al. 2020 iNEXT 4-steps).
###   For each assemblage the raw point-count counts are POOLED across days, repeat surveys and point counts into one species-by-assemblage abundance vector, then iNEXT4steps gives asymptotic (q = 1, 2) and coverage-standardised non-asymptotic (q = 0, at Cmax) diversity with bootstrap uncertainty. Uses all the data and equalises completeness, but is blind to spatial vs temporal replication (pooling loses it) and leaves a `Num_pc` signal in the estimate. Two passes: all assemblages, then the low-coverage ones dropped (removal set taken from pass 1's own SC_obs).
###   This section is the slow part (~1 h). Its logic is UNCHANGED from the previous version of this script (only moved / relabelled). Set `run_abundance_sensitivity` = FALSE to re-do only the fast incidence estimate.

### EXPORTS (one row per [assemblage x Order.q], keyed by `Assemblage` + the row_id parts):
###   Derived/Excels/Tax_div_incidence_<date>.csv    -- ESTIMATE A. qD + CI + delta-method SE, coverage at m*, interp/extrap flag, Num_pc, Num.hab.
###   Derived/Excels/Tax_div_all_farms_<date>.csv     -- ESTIMATE B. asymptotic q = 1 / q = 2 (`TD_asy`, `s.e.`) + observed / SC / evenness.
###   Derived/Excels/Tax_div_coverage65_<date>.csv    -- ESTIMATE B. non-asymptotic q = 0 at Cmax (`No_Asy_TD`) + coverage-based SE, low-coverage assemblages removed.
###   Figures/Incidence_vs_iNEXT.png                  -- the two estimates against each other + their residual dependence on Num_pc.

### Reference write-up (Cmax rationale, the iNEXT 4-steps framework, sample-completeness / rarefaction plots): Scripts/qmd/_archive/02_Analysis_iNEXT.qmd.

# Setup ----
library(tidyverse)
library(iNEXT)
library(iNEXT.4steps)
library(cowplot)

source("Scripts/Model_fns.R")   # latest_file()

### `conflicted` not used -- iNEXT / iNEXT.4steps do not mask the dplyr verbs. If a package added later does, qualify the call (dplyr::select, etc.) rather than reaching for conflicted.

ggplot2::theme_set(theme_cowplot())

wrangling_excels <- "../Ssp-bird-data-wrangling/Derived/Excels"
dir.create("Derived/Excels", recursive = TRUE, showWarnings = FALSE)
dir.create("Figures", showWarnings = FALSE)
set.seed(1989)

## Row identifier for an assemblage
row_id <- c("Uniq_db", "Id_gcs", "Ano_grp", "Season")

## Bootstrap replicates -- 500 so the per-assemblage SEs are stable enough to use as measurement error downstream (Scripts/04a). Shared by the incidence and the abundance passes.
nboot <- 500

## ESTIMATE A (incidence) standardisation level and the minimum effort to keep
target_level <- 6
min_pc <- 3

## ESTIMATE B (abundance) low-coverage cutoff for pass 2
sc_threshold <- 0.65

## The abundance sensitivity is the slow part (iNEXT4steps x 2 passes, ~1 h). Its logic is unchanged from the previous version of this script, so set FALSE to re-do only the fast incidence estimate; the comparison section then reads the latest existing abundance exports. Set TRUE for a full refresh (new point-count data, or new date stamp for all three).
run_abundance_sensitivity <- TRUE

today <- format(Sys.Date(), "%m.%d.%y")

# Load and link the point-count data ----

Bird_pcs <- read_csv(file.path(wrangling_excels, "Bird_pcs/Bird_pcs_analysis.csv"),
                     show_col_types = FALSE)
Site_covs <- read_csv(file.path(wrangling_excels, "Site_covs.csv"), show_col_types = FALSE)
Event_covs <- read_csv(file.path(wrangling_excels, "Event_covs.csv"), show_col_types = FALSE)

## One row per point-count survey with the assemblage identifier and habitat attached
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
  filter(!is.na(Id_gcs), Season != "Late") %>%
  mutate(Id_muestreo = str_replace_all(Id_muestreo, "-", "_"))

## Number of distinct habitat types surveyed per assemblage -- a sampling-scope covariate carried in both exports
num_hab <- Event_covs2 %>%
  filter(!is.na(Id_gcs)) %>%
  distinct(Identifier, Id_muestreo, Habitat) %>%
  summarize(Num.hab = as.factor(n_distinct(Habitat)), .by = Identifier)

# ============================================================================
# ESTIMATE A: point-count-standardised diversity (iNEXT incidence) ----
# ============================================================================

## [assemblage x point count x species] -- detections summed across repeat visits, then binarised
Pc_species <- Bird_pcs_linked %>%
  summarize(Count = sum(Count), .by = c(Identifier, Id_muestreo, Species_ayerbe))

## Per-assemblage species x point-count incidence (0/1) matrix -- the format estimateD() wants for datatype = "incidence_raw"
incidence_matrix <- function(df) {
  df %>%
    mutate(present = 1L) %>%
    select(Id_muestreo, Species_ayerbe, present) %>%
    pivot_wider(names_from = Id_muestreo, values_from = present, values_fill = 0L) %>%
    column_to_rownames("Species_ayerbe") %>%
    as.matrix()
}

inc_list <- Pc_species %>%
  group_split(Identifier) %>%
  set_names(map_chr(., ~ .x$Identifier[[1]])) %>%
  map(incidence_matrix)

Num_pc <- map_int(inc_list, ncol)
kept <- names(Num_pc)[Num_pc >= min_pc]

message("ESTIMATE A (incidence): ", length(inc_list), " assemblages, ", length(kept),
        " kept (>= ", min_pc, " point counts), standardised to ", target_level, " point counts")

inc_est <- estimateD(inc_list[kept], q = c(0, 1, 2), datatype = "incidence_raw",
                     base = "size", level = target_level, nboot = nboot)

Tax_div_incidence <- inc_est %>%
  as_tibble() %>%
  transmute(
    Assemblage,
    Order.q = as.integer(Order.q),
    Hill = c("0" = "richness", "1" = "shannon", "2" = "simpson")[as.character(Order.q)],
    qD, qD_LCL = qD.LCL, qD_UCL = qD.UCL,
    ## delta-method-style SE from the symmetric CI half-width (as for the coverage-based export)
    qD_se = (qD.UCL - qD.LCL) / (2 * qnorm(0.975)),
    SC_at_m = SC,
    method = Method,
    m_star = target_level
  ) %>%
  left_join(tibble(Assemblage = names(Num_pc), Num_pc = as.integer(Num_pc)), by = "Assemblage") %>%
  left_join(num_hab, by = c("Assemblage" = "Identifier")) %>%
  separate_wider_delim(Assemblage, delim = ".", names = row_id,
                       too_few = "align_start", cols_remove = FALSE) %>%
  relocate(Assemblage, Hill, Order.q)

write_csv(Tax_div_incidence, sprintf("Derived/Excels/Tax_div_incidence_%s.csv", today))
cat("Wrote Tax_div_incidence_", today, ".csv: ",
    n_distinct(Tax_div_incidence$Assemblage), " assemblages.\n", sep = "")

# ============================================================================
# ESTIMATE B: abundance diversity, coverage-standardised (iNEXT 4-steps) ----
# ============================================================================

## Skipped when run_abundance_sensitivity = FALSE -- the comparison section then reads the latest existing Tax_div_all_farms / _coverage65 exports.
if (run_abundance_sensitivity) {

## Sum counts within [assemblage x day x repeat], then across days, then widen to species x assemblage. Summing (not averaging) keeps the abundances iNEXT needs consistent with the accumulating species pool.
build_matrix <- function(bird_data) {
  bird_data %>%
    mutate(Identifier = pmap_chr(across(all_of(row_id)), ~ paste(..., sep = ".")),
           Identifier = str_replace_all(Identifier, " |-", "_")) %>%
    summarize(Count = sum(Count),
              .by = c(Identifier, Species_ayerbe, Fecha, Rep_ano_grp)) %>%
    summarize(Count = sum(Count), .by = c(Identifier, Species_ayerbe)) %>%
    pivot_wider(names_from = Identifier, values_from = Count, values_fill = 0) %>%
    rename_with(~ str_replace_all(.x, " |-", "_")) %>%
    column_to_rownames("Species_ayerbe")
}

## One iNEXT 4-steps pass -> tidy per-assemblage summary
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

## Pass 1: all assemblages
message("ESTIMATE B pass 1 (all assemblages)")
td_all <- run_inext_pass(Bird_pcs_linked)
write_csv(td_all, sprintf("Derived/Excels/Tax_div_all_farms_%s.csv", today))
cat("Wrote Tax_div_all_farms_", today, ".csv: ",
    n_distinct(td_all$Assemblage), " assemblages.\n", sep = "")

## Pass 2: drop the low-coverage assemblages, re-estimate
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

message("\nESTIMATE B pass 2 (low-coverage assemblages removed)")
td_c65 <- run_inext_pass(Bird_pcs_high_sc)
write_csv(td_c65, sprintf("Derived/Excels/Tax_div_coverage65_%s.csv", today))
cat("Wrote Tax_div_coverage65_", today, ".csv: ",
    n_distinct(td_c65$Assemblage), " assemblages.\n", sep = "")

}   # end if (run_abundance_sensitivity)

# ============================================================================
# Compare the primary (incidence) and sensitivity (Cmax) estimates ----
# ============================================================================

## The current-pipeline (abundance) value for each Hill, from the latest exports (just-written if the passes ran, else the most recent on disk): q = 0 non-asymptotic at Cmax, q = 1 / 2 asymptotic.
abund <- bind_rows(
  read_csv(latest_file("Derived/Excels", "^Tax_div_coverage65_.*\\.csv$"), show_col_types = FALSE) %>%
    filter(Order.q == 0) %>% transmute(Assemblage, Hill = "richness", abund_qD = No_Asy_TD),
  read_csv(latest_file("Derived/Excels", "^Tax_div_all_farms_.*\\.csv$"), show_col_types = FALSE) %>%
    filter(Order.q %in% c(1, 2)) %>%
    transmute(Assemblage, Hill = if_else(Order.q == 1, "shannon", "simpson"), abund_qD = TD_asy)
)

Compare <- Tax_div_incidence %>%
  select(Assemblage, Hill, Num_pc, incidence_qD = qD, SC_at_m) %>%
  left_join(abund, by = c("Assemblage", "Hill")) %>%
  filter(is.finite(incidence_qD), is.finite(abund_qD), incidence_qD > 0, abund_qD > 0)

Cor_summary <- Compare %>%
  summarize(n = n(),
            spearman = round(cor(incidence_qD, abund_qD, method = "spearman"), 3),
            pearson_log = round(cor(log(incidence_qD), log(abund_qD)), 3),
            .by = Hill)

## Residual dependence on Num_pc: R^2 of log(estimate) ~ log(Num_pc), incidence vs Cmax
Effort_summary <- Compare %>%
  pivot_longer(c(incidence_qD, abund_qD), names_to = "source", values_to = "value") %>%
  summarize(r2_vs_log_num_pc = round(summary(lm(log(value) ~ log(Num_pc)))$r.squared, 3),
            .by = c(Hill, source)) %>%
  pivot_wider(names_from = source, values_from = r2_vs_log_num_pc) %>%
  rename(incidence = incidence_qD, cmax = abund_qD)

hill_lab <- c(richness = "Richness (q = 0)", shannon = "Shannon (q = 1)", simpson = "Simpson (q = 2)")

## per-panel correlation label for p_scatter (replaces a standalone table in the qmd); position at each facet's bottom-right from the data range
Cor_labels <- Cor_summary %>%
  left_join(Compare %>% summarize(xpos = max(abund_qD), ypos = min(incidence_qD), .by = Hill), by = "Hill") %>%
  mutate(Hill = factor(recode(Hill, !!!hill_lab), levels = unname(hill_lab)),
         lab = sprintf("r(log) = %.2f\nrho = %.2f   (n = %d)", pearson_log, spearman, n))

p_scatter <- Compare %>%
  mutate(Hill = factor(recode(Hill, !!!hill_lab), levels = unname(hill_lab))) %>%
  ggplot(aes(abund_qD, incidence_qD, colour = Num_pc <= target_level)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
  geom_point(alpha = 0.75, size = 1.7) +
  geom_text(data = Cor_labels, aes(x = xpos, y = ypos, label = lab), inherit.aes = FALSE,
            hjust = 1, vjust = 0, size = 3, colour = "grey25", lineheight = 0.9) +
  scale_x_log10() + scale_y_log10() +
  scale_colour_manual(values = c(`FALSE` = "#08519c", `TRUE` = "#d95f02"),
                      name = sprintf("Num_pc <= %d\n(interp / extrap)", target_level)) +
  facet_wrap(~ Hill, scales = "free") +
  labs(x = "Sensitivity: abundance estimate, coverage-standardised (Cmax)",
       y = sprintf("Primary: incidence estimate, standardised to %d point counts", target_level),
       title = "Primary vs sensitivity diversity estimate")

p_effort <- Compare %>%
  pivot_longer(c(incidence_qD, abund_qD), names_to = "source", values_to = "value") %>%
  mutate(Hill = factor(recode(Hill, !!!hill_lab), levels = unname(hill_lab)),
         source = recode(source, incidence_qD = "primary (incidence @ m*)", abund_qD = "sensitivity (abundance, Cmax)")) %>%
  ggplot(aes(Num_pc, value, colour = source)) +
  geom_point(alpha = 0.5, size = 1.3) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
  scale_x_log10() + scale_y_log10() +
  scale_colour_manual(values = c(`sensitivity (abundance, Cmax)` = "#d95f02",
                                 `primary (incidence @ m*)` = "#08519c"), name = NULL) +
  facet_wrap(~ Hill, scales = "free_y") +
  labs(x = "Number of point counts (Num_pc)", y = "Diversity estimate",
       title = "Residual dependence on survey effort",
       subtitle = "The primary (incidence-at-m*) estimate is flat in Num_pc; the coverage-based one is not")

p_compare <- plot_grid(p_scatter, p_effort, ncol = 1)
ggsave("Figures/Incidence_vs_iNEXT.png", p_compare, bg = "white", width = 11, height = 9, dpi = 150)
print(p_compare)

# Report ----

cat("\n== ESTIMATE A: sample coverage of the incidence estimates at m* =", target_level, "point counts ==\n")
Tax_div_incidence %>%
  filter(Order.q == 0) %>%
  summarize(min = round(min(SC_at_m), 2), median = round(median(SC_at_m), 2),
            max = round(max(SC_at_m), 2), n_extrapolated = sum(method == "Extrapolation")) %>%
  print()

cat("\n== Correlation, primary vs sensitivity ==\n")
print(Cor_summary)

cat("\n== R^2 of log(estimate) ~ log(Num_pc): primary (incidence) vs sensitivity (Cmax) ==\n")
print(Effort_summary)
