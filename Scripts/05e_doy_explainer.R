# What the cyclic day-of-year terms (doy_sin + doy_cos) do ----

### The models include doy_sin = sin(2*pi*doy/365) and doy_cos = cos(2*pi*doy/365), where doy is the assemblage's mean survey day (1-365). Together these two terms let the model fit ONE smooth sine wave of any amplitude and any phase, with no discontinuity at 31 Dec / 1 Jan -- a control for the seasonal wave of Nearctic migrants (present ~Oct-Mar), which inflates diversity estimates when surveys fall in that window.

### This script makes an explainer figure: (A) the two basis waves; (B) the fitted seasonal effect from a real model, on the diversity scale, with the actual survey days shown. Sourced by nothing; run standalone or for a report.

# Setup ----
library(tidyverse)
library(brms)
library(cowplot)
ggplot2::theme_set(theme_cowplot(11))
dir.create("Figures", showWarnings = FALSE)

model_path <- "Derived/models/mod_shannon__baseline__climate__primary.rds"

# (A) the two basis functions ----

days <- tibble(doy = 1:365) %>%
  mutate(`sin(2 pi doy / 365)` = sin(2 * pi * doy / 365),
         `cos(2 pi doy / 365)` = cos(2 * pi * doy / 365)) %>%
  pivot_longer(-doy, names_to = "term", values_to = "value")

month_breaks <- yday(seq(as.Date("2021-01-01"), as.Date("2021-12-01"), by = "month"))
month_labs   <- month.abb

p_basis <- ggplot(days, aes(doy, value, colour = term)) +
  annotate("rect", xmin = yday(as.Date("2021-09-01")), xmax = 365, ymin = -Inf, ymax = Inf, alpha = 0.08, fill = "#2166ac") +
  annotate("rect", xmin = 0, xmax = yday(as.Date("2021-04-30")), ymin = -Inf, ymax = Inf, alpha = 0.08, fill = "#2166ac") +
  geom_line(linewidth = 1) +
  scale_x_continuous(breaks = month_breaks, labels = month_labs, expand = expansion(0)) +
  scale_colour_manual(values = c("sin(2 pi doy / 365)" = "#1b7837", "cos(2 pi doy / 365)" = "#d95f02"), name = NULL) +
  labs(x = NULL, y = "basis value",
       title = "(A) The two terms are phase-shifted waves",
       subtitle = "b1 * sin + b2 * cos = one sine wave of any amplitude and phase; continuous across the year-end. Shaded ~ when Nearctic migrants are in Colombia (Sep-Apr)") +
  theme(legend.position = "bottom")

# (B) the fitted seasonal effect from a real model ----

fit <- readRDS(model_path)
dr  <- as_draws_df(fit)

grid <- tibble(doy = 1:365) %>%
  mutate(s = sin(2 * pi * doy / 365), c = cos(2 * pi * doy / 365))

## posterior of the seasonal contribution to log diversity, centred so the annual mean is 0
season <- map_dfr(seq_len(nrow(grid)), function(i) {
  lp <- dr$b_doy_sin * grid$s[i] + dr$b_doy_cos * grid$c[i]
  tibble(doy = grid$doy[i], est = median(lp), lo = quantile(lp, 0.05), hi = quantile(lp, 0.95))
})
yr_mean <- mean(season$est)                                 # centre on the annual mean
season <- season %>%
  mutate(est = exp(est - yr_mean), lo = exp(lo - yr_mean), hi = exp(hi - yr_mean))  # -> multiplicative effect

## the actual assemblage survey days that this model was fit on
md <- read_csv("Derived/Excels/Farm_mgmt_model_data.csv", show_col_types = FALSE) %>%
  filter(Hill == "shannon") %>% distinct(Assemblage, doy)

amp <- with(dr, sqrt(b_doy_sin^2 + b_doy_cos^2))
p_fit <- ggplot(season, aes(doy, est)) +
  annotate("rect", xmin = yday(as.Date("2021-09-01")), xmax = 365, ymin = -Inf, ymax = Inf, alpha = 0.08, fill = "#2166ac") +
  annotate("rect", xmin = 0, xmax = yday(as.Date("2021-04-30")), ymin = -Inf, ymax = Inf, alpha = 0.08, fill = "#2166ac") +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey60") +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2, fill = "#762a83") +
  geom_line(linewidth = 1, colour = "#762a83") +
  geom_rug(data = md, aes(x = doy), sides = "b", alpha = 0.5, inherit.aes = FALSE) +
  scale_x_continuous(breaks = month_breaks, labels = month_labs, expand = expansion(0)) +
  labs(x = NULL, y = "seasonal multiplier on diversity",
       title = "(B) The fitted seasonal effect (Shannon, no-index baseline, climate spec)",
       subtitle = sprintf("Amplitude ~ %.0f%% peak-to-mean (90%% CrI %.0f-%.0f%%); ticks = the %d assemblages' mean survey day",
                          100 * (exp(median(amp)) - 1),
                          100 * (exp(quantile(amp, 0.05)) - 1), 100 * (exp(quantile(amp, 0.95)) - 1),
                          nrow(md)))

p_doy <- plot_grid(p_basis, p_fit, ncol = 1, align = "v", rel_heights = c(1, 1.05))
ggsave("Figures/doy_terms_explainer.png", p_doy, width = 10, height = 8, bg = "white")
print(p_doy)
