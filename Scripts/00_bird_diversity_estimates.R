# Bird taxonomic-diversity (Hill number) estimates per assemblage ----

### An assemblage = [data collector . farm . year-group . season]. Two diversity estimates per assemblage from the same point-count data; methods detail and caveats for the write-up are in Project_notes.md ("For the paper -- methods & discussion points").
## ESTIMATE A -- point-count-standardised (iNEXT incidence): diversity at m* = 6 point counts; assemblages with < 3 point counts dropped.
## ESTIMATE B -- abundance, coverage-standardised (iNEXT 4-steps): asymptotic q = 1 / 2, and q = 0 at Cmax with assemblages of coverage < 0.65 dropped. The slow part (~1 h); set `run_abundance_sensitivity` = FALSE to skip it.

### EXPORTS (one row per [assemblage x Order.q]) -> Derived/Excels/:
## Tax_div_incidence_<date>.csv -- Estimate A
## Tax_div_all_farms_<date>.csv -- Estimate B, all assemblages (q = 1 / 2 used downstream)
## Tax_div_coverage65_<date>.csv -- Estimate B, low-coverage assemblages removed (q = 0 used downstream)
## Figures/Incidence_vs_iNEXT.png -- the two estimates against each other + their residual dependence on Num_pc

# Setup ----
library(tidyverse)
library(iNEXT)
library(iNEXT.4steps)
library(cowplot)

source("Scripts/Shared_paths.R")   # WRANGLING_EXCELS, tax_div_path()

### `conflicted` not used -- iNEXT / iNEXT.4steps do not mask the dplyr verbs. If a package added later does, qualify the call (dplyr::select, etc.) rather than reaching for conflicted.

ggplot2::theme_set(theme_cowplot())

dir.create("Derived/Excels", recursive = TRUE, showWarnings = FALSE)
dir.create("Figures", showWarnings = FALSE)
set.seed(1989)

## Row identifier for an assemblage
row_id <- c("Uniq_db", "Id_scr", "Year_grp", "Season")

## Bootstrap replicates -- 500 so the per-assemblage SEs are stable enough to use as measurement error downstream (Scripts/04a). Shared by the incidence and the abundance passes.
nboot <- 500

## ESTIMATE A (incidence) standardisation level and the minimum effort to keep
target_level <- 6
min_pc <- 3

## ESTIMATE B (abundance) low-coverage cutoff for pass 2
sc_threshold <- 0.65

## The abundance sensitivity is the slow part (iNEXT4steps x 2 passes, ~1 h). Its logic is unchanged from the previous version of this script, so set FALSE to re-do only the fast incidence estimate; the comparison section then reads the abundance exports pinned by TAX_DIV_VERSION (Shared_paths.R). Set TRUE for a full refresh (new point-count data, or new date stamp for all three), then bump TAX_DIV_VERSION to the new date.
run_abundance_sensitivity <- TRUE

today <- format(Sys.Date(), "%m.%d.%y")

# Load and link the point-count data ----

Bird_pcs <- read_csv(file.path(WRANGLING_EXCELS, "Bird_pcs/Bird_pcs_analysis.csv"),
                     show_col_types = FALSE)
Site_covs <- read_csv(file.path(WRANGLING_EXCELS, "Site_covs.csv"), show_col_types = FALSE)
Event_covs <- read_csv(file.path(WRANGLING_EXCELS, "Event_covs.csv"), show_col_types = FALSE)

## One row per point-count survey with the assemblage identifier and habitat attached
Event_covs2 <- Site_covs %>%
  select(Id_survey_no_dc, Id_scr, Habitat) %>%
  left_join(
    Event_covs %>% select(Id_survey_no_dc, Id_survey, Uniq_db, Date, Pc_start,
                          Year_grp, Season, Rep_year_grp),
    by = "Id_survey_no_dc"
  ) %>%
  mutate(
    Identifier = pmap_chr(across(all_of(row_id)), ~ paste(..., sep = ".")),
    Identifier = str_replace_all(Identifier, " |-", "_")
  )

Bird_pcs_linked <- Bird_pcs %>%
  left_join(Event_covs2, by = c("Id_survey", "Id_survey_no_dc", "Date", "Pc_start")) %>%
  filter(!is.na(Id_scr), Season != "Late") %>%
  mutate(Id_survey = str_replace_all(Id_survey, "-", "_"))

## Number of distinct habitat types surveyed per assemblage -- a sampling-scope covariate carried in both exports
num_hab <- Event_covs2 %>%
  filter(!is.na(Id_scr)) %>%
  distinct(Identifier, Id_survey, Habitat) %>%
  summarize(Num.hab = as.factor(n_distinct(Habitat)), .by = Identifier)

# ESTIMATE A: point-count-standardised diversity (iNEXT incidence) ----

## [assemblage x point count x species] -- detections summed across repeat visits, then binarised
Pc_species <- Bird_pcs_linked %>%
  summarize(Count = sum(Count), .by = c(Identifier, Id_survey, Species_ayerbe))

## Per-assemblage species x point-count incidence (0/1) matrix -- the format estimateD() wants for datatype = "incidence_raw"
incidence_matrix <- function(df) {
  df %>%
    mutate(present = 1L) %>%
    select(Id_survey, Species_ayerbe, present) %>%
    pivot_wider(names_from = Id_survey, values_from = present, values_fill = 0L) %>%
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

# ESTIMATE B: abundance diversity, coverage-standardised (iNEXT 4-steps) ----

## Skipped when run_abundance_sensitivity = FALSE -- the comparison section then reads the pinned (TAX_DIV_VERSION) Tax_div_all_farms / _coverage65 exports.
if (run_abundance_sensitivity) {

## Sum counts within [assemblage x day x repeat], then across days, then widen to species x assemblage. Summing (not averaging) keeps the abundances iNEXT needs consistent with the accumulating species pool.
build_matrix <- function(bird_data) {
  bird_data %>%
    mutate(Identifier = pmap_chr(across(all_of(row_id)), ~ paste(..., sep = ".")),
           Identifier = str_replace_all(Identifier, " |-", "_")) %>%
    summarize(Count = sum(Count),
              .by = c(Identifier, Species_ayerbe, Date, Rep_year_grp)) %>%
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

# Compare incidence and Cmax estimates ----

## The current-pipeline (abundance) value for each Hill: q = 0 non-asymptotic at Cmax, q = 1 / 2 asymptotic. From the passes just run, else from the pinned exports on disk.
if (!run_abundance_sensitivity) {
  td_c65 <- read_csv(tax_div_path("coverage65"), show_col_types = FALSE)
  td_all <- read_csv(tax_div_path("all_farms"), show_col_types = FALSE)
}
abund <- bind_rows(
  td_c65 %>%
    filter(Order.q == 0) %>% transmute(Assemblage, Hill = "richness", abund_qD = No_Asy_TD),
  td_all %>%
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
