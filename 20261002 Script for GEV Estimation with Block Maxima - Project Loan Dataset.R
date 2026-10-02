# GEV Estimation for Cohort-based Hazard Matrix with Annual Block Maxima, MLE
# (c) Andreas Jobst, July 2026

pkgs <- c("readxl", "dplyr", "tidyr", "ggplot2", "stringr")
missing_pkgs <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
if (length(missing_pkgs)) install.packages(missing_pkgs, repos = "https://cloud.r-project.org")
invisible(lapply(pkgs, library, character.only = TRUE))

LOAN_FILE  <- "20260703 Project Loan Dataset - LET AAJ22 Andreas Jobst.xlsx"
CLASS_FILE <- "Categorization.xlsx"

ASOF_YEAR   <- 2024   # data as-of date, active loans get censored here
MAX_AGE     <- 10
MIN_AT_RISK <- 30     # ignore age cells with too few loans at risk
Z90 <- qnorm(0.95)
Z95 <- qnorm(0.975)

# --- load the loan book -----------------------------------------------------

raw <- read_excel(LOAN_FILE, sheet = LOAN_SHEET)

df <- raw %>%
  rename(
    orig_year    = `Origination Year`,
    event_year   = `Maturity/ End Year`,
    loan_status  = `Loan status`,
    default_flag = Default,
    country      = Country,
    region       = Region,
    industry     = Industry,
    sub_industry = `Sub-Industry`
  ) %>%
  mutate(
    orig_year    = as.integer(orig_year),
    event_year   = as.integer(event_year),
    default_flag = tolower(trimws(as.character(default_flag))) == "yes"
  )

# loan_status: 0 active (censor at ASOF_YEAR), -1 prepaid, -2 matured, -3 defaulted
df <- df %>% mutate(censor_year = if_else(loan_status == 0, pmin(event_year, ASOF_YEAR), event_year))

n_raw <- nrow(df)
df <- df %>% mutate(age_at_event = censor_year - orig_year)
n_excluded <- sum(df$age_at_event < 0, na.rm = TRUE)
df <- df %>% filter(age_at_event >= 0)   # a handful of records have an event before origination -- junk, drop them

cat(sprintf("Loaded %d loans, dropped %d bad records (%.2f%%)\n",
            n_raw, n_excluded, 100 * n_excluded / n_raw))

# --- country grouping --------------------------------------------------------
# Country_Classification has 4 columns (EEA, OECD, EMDE-A, EMDE-B), each just
# a list of country names. A few rows are region names instead of countries
# (Nordics, Western Europe, Asia, Africa...) -- those get used as a fallback
# when the loan's own Country field is itself a region.

country_class_raw <- read_excel(CLASS_FILE, sheet = "Country_Classification",
                                 col_names = c("EEA", "OECD", "EMDE_A", "EMDE_B"),
                                 skip = 2)

region_tokens <- c("Nordics", "Western Europe", "Asia", "Africa", "Europe", "Oceania",
                    "Middle East", "Eastern Europe", "Latin America", "North America")

group_cols <- c("EEA", "OECD", "EMDE_A", "EMDE_B")
country_sets <- list()
region_sets  <- list()
for (g in group_cols) {
  vals <- na.omit(trimws(country_class_raw[[g]]))
  vals <- vals[vals != ""]
  country_sets[[g]] <- tolower(vals[!(vals %in% region_tokens)])
  region_sets[[g]]  <- tolower(vals[vals %in% region_tokens])
}

# a few spelling variants that show up in the loan file but not the classification sheet
aliases <- c(
  "netherlands the"                  = "netherlands",
  "netherlands  the"                 = "netherlands",
  "russian federation"                = "russia",
  "ivory coast"                        = "côte d'ivoire",
  "lao people's democratic republic"  = "laos",
  "timor-leste"                        = "timor-leste (east timor)",
  "myanmar (burma)"                   = "myanmar"
)

norm_country <- function(x) {
  s <- str_squish(tolower(trimws(as.character(x))))
  ifelse(s %in% names(aliases), aliases[s], s)
}

df <- df %>%
  mutate(country_norm = norm_country(country),
         region_norm  = tolower(trimws(as.character(region))))

in_group <- function(country_norm, region_norm, g) {
  (country_norm %in% country_sets[[g]]) |
    (country_norm %in% region_sets[[g]]) |
    (region_norm  %in% region_sets[[g]])
}

for (g in group_cols) df[[paste0("grp_", g)]] <- in_group(df$country_norm, df$region_norm, g)
df$grp_EEA_OECD <- df$grp_EEA | df$grp_OECD

n_country_unclassified <- sum(!(df$grp_EEA_OECD | df$grp_EMDE_A | df$grp_EMDE_B))
cat(sprintf("Loans with no country group match: %d (%.2f%%) -- mostly GCC countries not in the classification sheet\n",
            n_country_unclassified, 100 * n_country_unclassified / nrow(df)))

# --- sector grouping ---------------------------------------------------------
# Industry_Classification gives Green/Non-Green at a pretty granular
# sub-industry level. Our loan file's Sub-Industry field is coarser, so below
# is a manual mapping -- checked so each coarse label only ever maps to
# classification rows that agree on the flag. Anything not in either list
# below (public services & admin, telecom, mining, manufacturing, etc.) is
# left unclassified rather than guessed at.

green_subinds <- tolower(c(
  "Power - Renewable Electricity Generation",
  "Infrastructure - Environment (water systems; water desalination, waste treatment)",
  "Infrastructure - Environment (waste to energy)",
  "Oil & Gas - Biofuels"
))
nongreen_subinds <- tolower(c(
  "Power - Non-renewable Electricity Generation",
  "Power - Electricity Transmission & Distribution",
  "Infrastructure - Transport (bridges, tunnels, roads, railroads)",
  "Infrastructure - Ports & terminals",
  "Oil & Gas - Exploration & Production",
  "Oil & Gas - Distribution (pipelines)",
  "Oil & Gas - LNG",
  "Oil & Gas - Storage",
  "Oil & Gas - Refining",
  "Oil & Gas - Other",
  "Oil & Gas"
))

df <- df %>%
  mutate(
    sub_industry_norm = tolower(trimws(coalesce(sub_industry, ""))),
    is_green    = sub_industry_norm %in% green_subinds,
    is_nongreen = sub_industry_norm %in% nongreen_subinds,
    grp_infra_power = industry %in% c("Power", "Infrastructure")
  )

n_sector_unclassified <- sum(!(df$is_green | df$is_nongreen))
cat(sprintf("Loans with no sector flag: %d (%.2f%%)\n",
            n_sector_unclassified, 100 * n_sector_unclassified / nrow(df)))

# --- hazard matrix / block maxima (unchanged from the single-population script) ---

build_hazard_matrix <- function(sub) {
  vintages <- sort(unique(sub$orig_year))
  out <- list()
  for (v in vintages) {
    vsub <- sub %>% filter(orig_year == v)
    at_risk <- nrow(vsub)
    for (a in 1:MAX_AGE) {
      cal_year <- v + a
      if (cal_year > ASOF_YEAR) break
      defaults_a <- sum(vsub$age_at_event == a & vsub$default_flag)
      exits_a    <- sum(vsub$age_at_event == a & !vsub$default_flag)
      hz <- if (at_risk > 0) defaults_a / at_risk else NA_real_
      out[[length(out) + 1]] <- data.frame(
        vintage = v, age = a, calendar_year = cal_year, at_risk = at_risk,
        defaults = defaults_a, exits = exits_a, hazard = hz
      )
      at_risk <- at_risk - defaults_a - exits_a
    }
  }
  bind_rows(out)
}

block_maxima_from <- function(hz_df, min_at_risk = MIN_AT_RISK) {
  valid <- hz_df %>% filter(at_risk >= min_at_risk, !is.na(hazard))
  if (nrow(valid) == 0) return(data.frame(calendar_year = integer(0), hazard = numeric(0)))
  valid %>% group_by(calendar_year) %>% summarise(hazard = max(hazard), .groups = "drop") %>% arrange(calendar_year)
}

# --- GEV distribution + fitting ----------------------------------------------

pgev <- function(q, mu, sigma, xi) {
  z <- (q - mu) / sigma
  if (abs(xi) < 1e-8) return(exp(-exp(-z)))
  w <- pmax(1 + xi * z, 0)
  ifelse(w == 0 & xi > 0, 0, ifelse(w == 0 & xi < 0, 1, exp(-w^(-1 / xi))))
}
dgev <- function(x, mu, sigma, xi) {
  z <- (x - mu) / sigma
  if (abs(xi) < 1e-8) return((1 / sigma) * exp(-z - exp(-z)))
  w <- 1 + xi * z
  out <- rep(0, length(x)); ok <- w > 0
  out[ok] <- (1 / sigma) * w[ok]^(-1 / xi - 1) * exp(-w[ok]^(-1 / xi))
  out
}
qgev <- function(p, mu, sigma, xi) {
  if (abs(xi) < 1e-8) return(mu - sigma * log(-log(p)))
  mu - sigma / xi * (1 - (-log(p))^(-xi))
}

# reparametrize sigma = exp(phi) so the optimizer can't wander into sigma <= 0
negloglik_repar <- function(theta, x) {
  mu <- theta[1]; phi <- theta[2]; xi <- theta[3]
  sigma <- exp(phi)
  z <- (x - mu) / sigma
  if (abs(xi) < 1e-8) {
    t <- exp(-z)
    return(-sum(-phi - z - t))
  }
  w <- 1 + xi * z
  if (any(w <= 1e-10)) return(1e10)   # outside the support, just penalize hard
  t <- w^(-1 / xi)
  -sum(-phi - (1 + 1 / xi) * log(w) - t)
}
negloglik <- function(par, x) {
  mu <- par[1]; sigma <- par[2]; xi <- par[3]
  if (sigma <= 0) return(1e10)
  negloglik_repar(c(mu, log(sigma), xi), x)
}

# plain numerical Hessian/gradient -- no numDeriv on this box
num_hessian <- function(f, par, eps = 1e-5, ...) {
  k <- length(par); H <- matrix(0, k, k)
  for (i in 1:k) for (j in 1:k) {
    pp <- par; pp[i] <- pp[i] + eps; pp[j] <- pp[j] + eps
    pm <- par; pm[i] <- pm[i] + eps; pm[j] <- pm[j] - eps
    mp <- par; mp[i] <- mp[i] - eps; mp[j] <- mp[j] + eps
    mm <- par; mm[i] <- mm[i] - eps; mm[j] <- mm[j] - eps
    H[i, j] <- (f(pp, ...) - f(pm, ...) - f(mp, ...) + f(mm, ...)) / (4 * eps^2)
  }
  H
}
num_grad <- function(f, par, eps = 1e-6, ...) {
  k <- length(par); g <- numeric(k)
  for (i in 1:k) {
    pp <- par; pp[i] <- pp[i] + eps
    pm <- par; pm[i] <- pm[i] - eps
    g[i] <- (f(pp, ...) - f(pm, ...)) / (2 * eps)
  }
  g
}

fit_gev_robust <- function(x) {
  n <- length(x)
  sd_x <- sd(x)
  sigma0 <- max(sqrt(6) * sd_x / pi, 1e-4)
  mu0 <- mean(x) - 0.5772157 * sigma0
  sigma_floor_rel <- 0.15

  candidates <- list()
  for (xi0 in c(-0.3, -0.1, 0.0, 0.1, 0.2, 0.3, 0.5)) {
    theta0 <- c(mu0, log(sigma0), xi0)
    res <- tryCatch(
      optim(theta0, negloglik_repar, x = x, method = "Nelder-Mead",
            control = list(maxit = 40000, reltol = 1e-14)),
      error = function(e) NULL
    )
    if (!is.null(res) && is.finite(res$value) && res$value < 1e9) candidates[[length(candidates) + 1]] <- res
  }

  if (length(candidates) == 0) {
    # every start blew up -- just run one more so we have *something* to
    # report, and it'll get flagged degenerate below regardless
    theta0 <- c(mu0, log(sigma0), 0.0)
    candidates <- list(optim(theta0, negloglik_repar, x = x, method = "Nelder-Mead", control = list(maxit = 40000)))
  }

  legit <- Filter(function(c) exp(c$par[2]) >= sigma_floor_rel * sd_x, candidates)
  boundary_only <- length(legit) == 0
  pool <- if (length(legit) > 0) legit else candidates
  best <- pool[[which.min(sapply(pool, function(c) c$value))]]

  mu <- best$par[1]; sigma <- exp(best$par[2]); xi <- best$par[3]
  loglik <- -best$value

  H <- num_hessian(negloglik, c(mu, sigma, xi), x = x)
  se <- rep(NA_real_, 3); cov_mat <- matrix(NA_real_, 3, 3)
  degenerate <- boundary_only
  cov_try <- tryCatch(solve(H), error = function(e) NULL)
  if (is.null(cov_try) || any(diag(cov_try) <= 0) || any(!is.finite(diag(cov_try)))) {
    degenerate <- TRUE
  } else {
    se <- sqrt(diag(cov_try)); cov_mat <- cov_try
  }
  if (sigma < 1e-4 || abs(xi) > 1.0 || any(!is.finite(se))) degenerate <- TRUE

  gev_mean_fn <- function(p) { mu <- p[1]; sigma <- p[2]; xi <- p[3]; mu + sigma * (gamma(1 - xi) - 1) / xi }
  gev_sd_fn <- function(p) {
    mu <- p[1]; sigma <- p[2]; xi <- p[3]
    g1 <- gamma(1 - xi); g2 <- gamma(1 - 2 * xi)
    sqrt(max(sigma^2 * (g2 - g1^2) / xi^2, 0))
  }
  gev_q_fn <- function(p, pr) qgev(pr, p[1], p[2], p[3])

  par_hat <- c(mu, sigma, xi)
  mean_hat <- NA_real_; sd_hat <- NA_real_
  if (abs(xi) < 0.98) mean_hat <- tryCatch(gev_mean_fn(par_hat), error = function(e) NA_real_)
  if (abs(xi) < 0.48) sd_hat   <- tryCatch(gev_sd_fn(par_hat),   error = function(e) NA_real_)

  se_mean <- se_sd <- se_rl90 <- se_rl95 <- NA_real_
  rl90 <- tryCatch(gev_q_fn(par_hat, 0.90), error = function(e) NA_real_)
  rl95 <- tryCatch(gev_q_fn(par_hat, 0.95), error = function(e) NA_real_)
  if (!degenerate) {
    # delta method: SE of a derived quantity = sqrt(grad' * cov * grad)
    if (is.finite(mean_hat)) { g <- num_grad(gev_mean_fn, par_hat); se_mean <- sqrt(max(as.numeric(t(g) %*% cov_mat %*% g), 0)) }
    if (is.finite(sd_hat))   { g <- num_grad(gev_sd_fn, par_hat);   se_sd   <- sqrt(max(as.numeric(t(g) %*% cov_mat %*% g), 0)) }
    g90 <- num_grad(function(p) gev_q_fn(p, 0.90), par_hat)
    g95 <- num_grad(function(p) gev_q_fn(p, 0.95), par_hat)
    se_rl90 <- sqrt(max(as.numeric(t(g90) %*% cov_mat %*% g90), 0))
    se_rl95 <- sqrt(max(as.numeric(t(g95) %*% cov_mat %*% g95), 0))
  }

  ks_test <- suppressWarnings(tryCatch(ks.test(x, pgev, mu = mu, sigma = sigma, xi = xi),
                                        error = function(e) list(statistic = NA, p.value = NA)))
  xs <- sort(x)
  Fx <- pmin(pmax(pgev(xs, mu, sigma, xi), 1e-12), 1 - 1e-12)
  i <- seq_len(n)
  A2 <- -n - mean((2 * i - 1) * (log(Fx) + log(1 - rev(Fx))))   # Anderson-Darling stat
  pp <- (i - 0.5) / n
  qth <- qgev(pp, mu, sigma, xi)
  ppcc <- suppressWarnings(cor(xs, qth))   # probability-plot correlation coefficient

  k <- 3
  list(mu = mu, sigma = sigma, xi = xi, se_mu = se[1], se_sigma = se[2], se_xi = se[3],
       mean = mean_hat, sd = sd_hat, se_mean = se_mean, se_sd = se_sd,
       rl90 = rl90, se_rl90 = se_rl90, rl95 = rl95, se_rl95 = se_rl95,
       ks_d = unname(ks_test$statistic), ks_p = ks_test$p.value, ad = A2, ppcc = ppcc,
       loglik = loglik, aic = 2 * k - 2 * loglik, bic = k * log(n) - 2 * loglik,
       n_blocks = n, degenerate = degenerate,
       emp_mean = mean(x), emp_sd = sd(x), emp_p90 = quantile(x, 0.90, names = FALSE),
       emp_p95 = quantile(x, 0.95, names = FALSE), n_zero = sum(x == 0))
}

# --- run all 16 cells ---------------------------------------------------------

country_groups <- list(
  All      = rep(TRUE, nrow(df)),
  EEA_OECD = df$grp_EEA_OECD,
  EMDE_A   = df$grp_EMDE_A,
  EMDE_B   = df$grp_EMDE_B
)
sector_groups <- list(
  All        = rep(TRUE, nrow(df)),
  InfraPower = df$grp_infra_power,
  Green      = df$is_green,
  NonGreen   = df$is_nongreen
)

all_results <- list()
results_rows <- list()

for (cg in names(country_groups)) {
  for (sg in names(sector_groups)) {
    sub <- df[country_groups[[cg]] & sector_groups[[sg]], ]
    hz  <- build_hazard_matrix(sub)
    bm  <- block_maxima_from(hz)
    x   <- bm$hazard
    fit <- if (length(x) >= 8) fit_gev_robust(x) else NULL   # too few block maxima to bother fitting

    key <- paste(cg, sg, sep = "_")
    all_results[[key]] <- list(hazard_matrix = hz, block_maxima = bm, fit = fit,
                                n_loans = nrow(sub), n_defaults = sum(sub$default_flag))

    row <- data.frame(country_group = cg, sector_group = sg, n_loans = nrow(sub),
                       n_defaults = sum(sub$default_flag),
                       n_blocks = if (!is.null(fit)) fit$n_blocks else length(x))
    if (!is.null(fit)) {
      row$mu <- fit$mu; row$sigma <- fit$sigma; row$xi <- fit$xi
      row$se_mu <- fit$se_mu; row$se_sigma <- fit$se_sigma; row$se_xi <- fit$se_xi
      row$mean <- fit$mean; row$sd <- fit$sd
      row$rl90 <- fit$rl90; row$rl95 <- fit$rl95
      row$ks_p <- fit$ks_p; row$ppcc <- fit$ppcc; row$aic <- fit$aic
      row$degenerate <- fit$degenerate
      row$emp_mean <- fit$emp_mean; row$emp_p90 <- fit$emp_p90; row$emp_p95 <- fit$emp_p95
    } else {
      row$degenerate <- NA
    }
    results_rows[[key]] <- row

    flag <- if (is.null(fit)) "insufficient data" else if (fit$degenerate) "DEGENERATE" else "ok"
    cat(sprintf("%-10s x %-10s  N=%5d  blocks=%3d  %s\n", cg, sg, nrow(sub), length(x), flag))
  }
}

results_table <- bind_rows(results_rows)
write.csv(results_table, "permutation_results_R.csv", row.names = FALSE)

cat("\n--- Full results table ---\n")
print(results_table)

# --- charts -------------------------------------------------------------------

theme_report <- theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())

# hazard age profile, a handful of sample vintages, baseline population
hz_all <- all_results[["All_All"]]$hazard_matrix
sample_vintages <- c(1995, 2000, 2005, 2008, 2010, 2013)
fig1_data <- hz_all %>% filter(vintage %in% sample_vintages)
fig1 <- ggplot(fig1_data, aes(x = age, y = 100 * hazard, colour = factor(vintage), group = vintage)) +
  geom_line() + geom_point(size = 1.3) +
  labs(title = "Hazard-rate age profile by origination vintage (all loans, all countries)",
       x = "Years since inception (age)", y = "Marginal hazard rate (%)", colour = "Vintage") +
  theme_report
ggsave("chart_vintage_hazard.png", fig1, width = 7, height = 4.2, dpi = 170)

# block maxima faceted by country group (sector = All)
bm_by_country <- bind_rows(lapply(names(country_groups), function(cg) {
  all_results[[paste0(cg, "_All")]]$block_maxima %>% mutate(country_group = cg)
}))
cg_labels <- c(All = "All", EEA_OECD = "EEA+OECD", EMDE_A = "EMDE-A", EMDE_B = "EMDE-B")
bm_by_country$country_group <- factor(cg_labels[bm_by_country$country_group], levels = cg_labels[names(country_groups)])
fig2 <- ggplot(bm_by_country, aes(x = calendar_year, y = 100 * hazard)) +
  geom_col(fill = "#2b4a6f") +
  facet_wrap(~country_group, nrow = 2) +
  labs(title = "Annual block-maximum hazard rate by country grouping (all sectors)",
       x = "Calendar year", y = "Block-maximum hazard rate (%)") +
  theme_report
ggsave("chart_block_maxima_by_country.png", fig2, width = 9, height = 6, dpi = 170)

# shape parameter across all 16 cells, with degenerate ones marked separately
sg_labels <- c(All = "All sectors", InfraPower = "Infra+Power", Green = "Green", NonGreen = "Non-Green")
xi_df <- results_table %>%
  mutate(label = paste(cg_labels[country_group], sg_labels[sector_group], sep = " – "),
         lo = xi - Z95 * se_xi, hi = xi + Z95 * se_xi)
xi_df$label <- factor(xi_df$label, levels = rev(xi_df$label))

fig3 <- ggplot(xi_df, aes(x = xi, y = label)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(data = filter(xi_df, !degenerate), aes(xmin = lo, xmax = hi), height = 0, colour = "#9fbfdf", linewidth = 1) +
  geom_point(data = filter(xi_df, !degenerate), colour = "#2b4a6f", size = 2.5) +
  geom_point(data = filter(xi_df, degenerate), shape = 4, colour = "#c0392b", size = 3.5, stroke = 1.5) +
  labs(title = "GEV shape parameter by country x sector permutation",
       x = "Shape parameter xi (95% CI where identified)", y = NULL) +
  theme_report
ggsave("chart_xi_comparison.png", fig3, width = 8.2, height = 6.2, dpi = 170)

# 1-in-10-year return level, same layout, empirical fallback for degenerate cells
rl_df <- results_table %>%
  mutate(label = paste(cg_labels[country_group], sg_labels[sector_group], sep = " – "),
         rl90_show = if_else(degenerate, emp_p90, rl90))
rl_df$label <- factor(rl_df$label, levels = rev(xi_df$label))
fig4 <- ggplot(rl_df, aes(x = 100 * rl90_show, y = label, fill = degenerate)) +
  geom_col() +
  scale_fill_manual(values = c(`FALSE` = "#2b4a6f", `TRUE` = "#c0392b"),
                     labels = c(`FALSE` = "GEV 90% return level", `TRUE` = "Empirical 90th pct. (fallback)")) +
  labs(title = "1-in-10-year stressed hazard rate by permutation",
       x = "1-in-10-year block-maximum hazard rate (%)", y = NULL, fill = NULL) +
  theme_report + theme(legend.position = "bottom")
ggsave("chart_rl90_comparison.png", fig4, width = 8.2, height = 6.2, dpi = 170)

# fitted density vs histogram, baseline population
f0 <- all_results[["All_All"]]$fit
x0 <- all_results[["All_All"]]$block_maxima$hazard
xs <- seq(max(0, min(x0) - 0.02), max(x0) + 0.04, length.out = 400)
dens_df <- data.frame(x = xs, y = dgev(xs, f0$mu, f0$sigma, f0$xi))
fig5 <- ggplot() +
  geom_histogram(data = data.frame(x = x0), aes(x = x, y = after_stat(density)), bins = 10, fill = "#9fbfdf", colour = "white") +
  geom_line(data = dens_df, aes(x = x, y = y), colour = "#c0392b", linewidth = 1) +
  labs(title = "Fitted GEV distribution (all loans, all countries)",
       x = "Block-maximum hazard rate", y = "Density") +
  theme_report
ggsave("chart_gev_density.png", fig5, width = 6.5, height = 4, dpi = 170)

cat("\nDone -- wrote permutation_results_R.csv and 5 chart PNGs.\n")
