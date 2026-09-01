# Shared helpers for the pipeline ----

### Sourced by every numbered script. `latest_file()` (below) has no dependencies -- the data-prep scripts use only that. The modelling pieces (priors, formula builder, Stan-compile-sharing grid fitter, tidiers) are used by Scripts/04a / 04b / 05 and assume `brms` + the tidyverse are already attached (those scripts load them); building `MGMT_PRIORS` needs the `brms` package installed but not attached.

# latest_file ----

## The most recently modified file in `dir` matching `pattern` -- lets consumers pick up the newest date-stamped Tax_div_* export from Scripts/00 without hard-coding the date
latest_file <- function(dir, pattern) {
  hits <- list.files(dir, pattern = pattern, full.names = TRUE)
  if (length(hits) == 0) stop("No file matching '", pattern, "' in ", dir, call. = FALSE)
  hits[which.max(file.mtime(hits))]
}

# Priors ----

## Weakly-informative priors on the log-diversity scale, shared by every farm-management model. The intercept prior sits near log(20) diversity; normal(0, 0.75) on b lightly regularises the standardised slopes; exponential(1) on the group SDs and the residual SD.
MGMT_PRIORS <- c(
  brms::prior(student_t(3, 3, 2.5), class = "Intercept"),
  brms::prior(normal(0, 0.75), class = "b"),
  brms::prior(exponential(1), class = "sd"),
  brms::prior(exponential(1), class = "sigma")
)

# Formula builder ----

## Assemble one brms formula from string parts.
##   lhs   -- the response side, including any `| resp_se(...)` term, e.g. "log_response | resp_se(se_log, sigma = TRUE)"
##   focal -- the focal predictor column name, or NULL for a no-index baseline
##   rhs   -- character vector of the remaining right-hand-side terms (region adjustment, sampling covariates, random effects), in the order they should appear
## Fits whose resulting formula, data and priors all match reuse a cached `.rds` (see `fit_model_grid()`), so keep `rhs` term order stable.
mgmt_bf <- function(lhs, focal, rhs) {
  terms <- c(if (!is.null(focal)) focal, rhs)
  brms::bf(stats::as.formula(paste(lhs, "~", paste(terms, collapse = " + "))))
}

# Grid fitter ----

## Fit a grid of brms models, each cached to `<model_dir>/<key>.rds` (`file_refit = "on_change"`).
##   grid       -- data frame with a unique `key` column (-> the `.rds` filename)
##   build_bf   -- function(row_list) -> a brmsformula (via `mgmt_bf()`)
##   build_data -- function(row_list) -> the model frame for that fit
##   prior      -- brms prior object, or function(row_list, data) -> prior
##   family     -- NULL (brms default gaussian), a brms family, or function(row_list) -> family
##   the sampler settings match Scripts/04a's Cmax-arm defaults; override per script as needed
## Each fit is an independent `brm()` call; rstan caches the compiled Stan model by code hash across
## calls in a session (`auto_write`), so distinct model structures still compile only once. (Sharing a
## fitted object via `update(recompile = FALSE)` breaks when the cached `.rds` -- which brms saves
## without the compiled model -- has to be refit in a later session.)
## Returns a list of fits in grid order, named by `key`.
fit_model_grid <- function(grid, build_bf, build_data, prior = MGMT_PRIORS, family = NULL,
                           chains = 4, iter = 4000, warmup = 1500,
                           adapt_delta = 0.995, seed = 1989,
                           model_dir = "Derived/models") {
  dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
  ## cache the compiled Stan model to disk, so identical model structures compile once even across independent brm() calls / sessions
  if (requireNamespace("rstan", quietly = TRUE)) try(rstan::rstan_options(auto_write = TRUE), silent = TRUE)
  fits <- vector("list", nrow(grid))
  names(fits) <- grid$key
  for (i in seq_len(nrow(grid))) {
    row <- as.list(grid[i, ])
    message(sprintf("[%d/%d] %s", i, nrow(grid), row$key))
    frame <- build_data(row)
    brm_args <- list(
      formula = build_bf(row), data = frame,
      prior = if (is.function(prior)) prior(row, frame) else prior,
      chains = chains, iter = iter, warmup = warmup, seed = seed,
      control = list(adapt_delta = adapt_delta), refresh = 0, silent = 2,
      file = file.path(model_dir, row$key), file_refit = "on_change"
    )
    fam <- if (is.function(family)) family(row) else family
    if (!is.null(fam)) brm_args$family <- fam
    fits[[i]] <- do.call(brms::brm, brm_args)
  }
  fits
}

# Fixed-effects tidier ----

## One row per [fit x population-level coefficient], plus per-fit fit / convergence stats.
##   fits -- named list from `fit_model_grid()` (names are the grid `key`s)
##   grid -- the grid data frame; joined back on `key` so hill / index / spec / ... travel with the coefficients
##   ci   -- central credible-interval mass (default 0.90 -> `conf_low` / `conf_high` at the 5th / 95th percentiles)
## The caller does any rounding, column selection and extra flags.
tidy_model_fits <- function(fits, grid, ci = 0.90) {
  lo <- (1 - ci) / 2
  hi <- 1 - lo
  purrr::imap(fits, function(fit, key) {
    draws <- posterior::as_draws_matrix(fit, variable = "^b_", regex = TRUE)
    coefs <- purrr::imap(asplit(draws, 2), ~ tibble::tibble(
      term = sub("^b_", "", .y),
      estimate = stats::median(.x),
      conf_low = stats::quantile(.x, lo, names = FALSE),
      conf_high = stats::quantile(.x, hi, names = FALSE),
      p_direction_pos = mean(.x > 0)
    )) |> purrr::list_rbind()
    coefs$key <- key
    coefs$n_obs <- stats::nobs(fit)
    coefs$bayes_R2 <- brms::bayes_R2(fit)[, "Estimate"]
    coefs$max_rhat <- round(max(brms::rhat(fit), na.rm = TRUE), 3)
    coefs$n_divergent <- sum(subset(brms::nuts_params(fit), Parameter == "divergent__")$Value)
    coefs
  }) |>
    purrr::list_rbind() |>
    dplyr::left_join(grid, by = "key")
}

## One row per fit: the focal coefficient (NA for a no-index baseline) plus per-fit fit / convergence stats. Use when only the focal effect matters (the incidence arm, the robustness checks), rather than every population-level coefficient.
focal_fit_summary <- function(fits, grid, focal = "focal_z", ci = 0.90) {
  lo <- (1 - ci) / 2
  hi <- 1 - lo
  bvar <- paste0("b_", focal)
  purrr::imap(fits, function(fit, key) {
    has_focal <- bvar %in% brms::variables(fit)
    d <- if (has_focal) as.numeric(posterior::as_draws_matrix(fit, variable = bvar)) else NA_real_
    tibble::tibble(
      key = key,
      n_obs = stats::nobs(fit),
      bayes_R2 = brms::bayes_R2(fit)[, "Estimate"],
      focal_est = if (has_focal) stats::median(d) else NA_real_,
      focal_lo = if (has_focal) stats::quantile(d, lo, names = FALSE) else NA_real_,
      focal_hi = if (has_focal) stats::quantile(d, hi, names = FALSE) else NA_real_,
      p_direction_pos = if (has_focal) mean(d > 0) else NA_real_,
      max_rhat = round(max(brms::rhat(fit), na.rm = TRUE), 3),
      n_divergent = sum(subset(brms::nuts_params(fit), Parameter == "divergent__")$Value)
    )
  }) |>
    purrr::list_rbind() |>
    dplyr::left_join(grid, by = "key")
}
