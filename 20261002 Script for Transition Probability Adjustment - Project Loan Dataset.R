# Adjustment of Infrastructure Project Loan PDs with Markov Chain of Rated Infra Debt
# (c) Andreas Jobst, July 2026

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(ggplot2)
  library(writexl)
})

xlsx_path <- "DATA_Moodys.xlsx"
ratings <- c("AAA","AA","A","BBB","BB","B","CCC to C")
years <- 1:10
state_names <- c(ratings, "D")

gamma <- 2.0  # dispersion sensitivity
k_sd <- 10.0  # scaling: how strongly SD maps into incremental PD relative to direct PD

parse_infra_transition_mats <- function(df_raw){
  year_rows <- which(str_detect(as.character(df_raw[[1]]), "^Year\\s*\\d+"))
  mats <- list()
  for(yr_idx in year_rows){
    yr_label <- as.character(df_raw[[1]][yr_idx])
    yr_num <- as.integer(str_extract(yr_label, "\\d+"))

    data_start <- yr_idx + 2
    data_end <- data_start + 6

    block <- df_raw[data_start:data_end, 1:8]
    colnames(block) <- c("From", ratings)
    block <- block %>% filter(From %in% ratings)

    mat <- as.matrix(block[, ratings])
    mat <- apply(mat, 2, as.numeric) / 100
    rownames(mat) <- block$From
    mat <- mat[ratings, ratings]
    mats[[as.character(yr_num)]] <- mat
  }
  mats
}

interpolate_mats_6_9 <- function(mats){
  mat5 <- mats[["5"]]
  mat10 <- mats[["10"]]
  for(yr in 6:9){
    w <- (yr - 5) / (10 - 5)
    mats[[as.character(yr)]] <- (1 - w) * mat5 + w * mat10
  }
  mats
}

# Option A calibration: enforce P_t(i->D) = rated marginal PD(i,t)
# Scale the 7x7 non-default transitions to sum to (1 - rated_marg)
build_calibrated_P <- function(tm7, rated_marg_year){
  P <- matrix(0, nrow = 8, ncol = 8)
  colnames(P) <- state_names
  rownames(P) <- state_names

  for(i in ratings){
    row_vals <- tm7[i, ]
    row_sum <- sum(row_vals)
    surv <- as.numeric(1 - rated_marg_year[i])
    if(row_sum <= 0){
      scaled <- rep(0, length(ratings))
    } else {
      scaled <- as.numeric(row_vals) * (surv / row_sum)
    }
    P[i, ratings] <- scaled
    P[i, "D"] <- 1 - sum(P[i, ratings])
  }

  P["D","D"] <- 1
  P
}

# Dispersion penalty: sd of off-diagonal non-default probs (excluding stay prob)
row_sd <- function(P){
  out <- rep(NA_real_, length(ratings))
  names(out) <- ratings
  for(i in ratings){
    row <- P[i, ratings]
    off <- row[names(row) != i]
    out[i] <- sd(off)
  }
  out
}

# Dynamic bracket selection based on unrated cumulative PD location in rated cumulative PDs
get_brackets <- function(rated_cum_long, unrated_cum_tbl){
  map_dfr(years, function(yr){
    rv <- rated_cum_long %>% filter(Year == yr) %>% select(Rating, cum) %>% arrange(cum)
    u <- unrated_cum_tbl %>% filter(Year == yr) %>% pull(cum)

    if(u <= rv$cum[1]){
      low <- rv$Rating[1]; high <- rv$Rating[2]
    } else if(u >= rv$cum[nrow(rv)]){
      low <- rv$Rating[nrow(rv)-1]; high <- rv$Rating[nrow(rv)]
    } else {
      hi_idx <- which(rv$cum >= u)[1]
      high <- rv$Rating[hi_idx]
      low <- rv$Rating[hi_idx-1]
    }

    low_val <- rv %>% filter(Rating == low) %>% pull(cum)
    high_val <- rv %>% filter(Rating == high) %>% pull(cum)

    if(high_val == low_val){
      w_high <- 0.5
    } else {
      w_high <- as.numeric((u - low_val) / (high_val - low_val))
    }
    w_low <- 1 - w_high

    tibble(Year = yr, Low = low, High = high, w_low = w_low, w_high = w_high)
  })
}

#-------------------------
# Read Default Rates
#-------------------------

def_raw <- read_excel(xlsx_path, sheet = "Default Rates", col_names = FALSE)

rated_hdr_idx <- which(str_detect(as.character(def_raw[[1]]), "Rated Infrastructure"))[1]
rated_block <- def_raw[(rated_hdr_idx + 1):(rated_hdr_idx + 7), 2:12]
colnames(rated_block) <- c("Rating", paste0("Year ", years))

rated_cum <- rated_block %>%
  filter(Rating %in% ratings) %>%
  mutate(across(starts_with("Year"), ~ as.numeric(.x) / 100))

rated_cum_long <- rated_cum %>%
  pivot_longer(cols = starts_with("Year"), names_to = "YearLabel", values_to = "cum") %>%
  mutate(Year = as.integer(str_extract(YearLabel, "\\d+"))) %>%
  select(Rating, Year, cum) %>%
  arrange(Rating, Year)

rated_marg_long <- rated_cum_long %>%
  group_by(Rating) %>%
  arrange(Year) %>%
  mutate(marg = if_else(Year == 1, cum, cum - lag(cum))) %>%
  ungroup()

hazard_1y <- rated_marg_long %>% filter(Year == 1) %>% select(Rating, h = marg)

unrated_idx <- which(str_detect(as.character(def_raw[[1]]), "Unrated"))[1]
unrated_vals <- def_raw[unrated_idx, 3:12] %>% unlist() %>% as.numeric() / 100
unrated_cum <- tibble(Year = years, cum = unrated_vals) %>%
  arrange(Year) %>%
  mutate(marg = if_else(Year == 1, cum, cum - lag(cum)))

#-------------------------
# Read transition matrices and interpolate
#-------------------------

trans_raw <- read_excel(xlsx_path, sheet = "Infra Trans Rates", col_names = FALSE)
trans_mats <- parse_infra_transition_mats(trans_raw)
trans_mats <- interpolate_mats_6_9(trans_mats)

#-------------------------
# Build calibrated Markov matrices P_t + penalty
#-------------------------

P_list <- list()
sd_tbl <- list()
for(yr in years){
  tm7 <- trans_mats[[as.character(yr)]]
  rated_marg_year <- rated_marg_long %>% filter(Year == yr) %>% select(Rating, marg)
  rated_marg_map <- setNames(rated_marg_year$marg, rated_marg_year$Rating)

  P <- build_calibrated_P(tm7, rated_marg_map)
  P_list[[as.character(yr)]] <- P

  sdv <- row_sd(P)
  sd_tbl[[as.character(yr)]] <- tibble(Year = yr, Rating = names(sdv), sd_offdiag = as.numeric(sdv))
}
sd_long <- bind_rows(sd_tbl)

#-------------------------
# Markov decomposition for rated: total marginal default vs direct vs indirect
#-------------------------

markov_decomp_one_start <- function(start_rating){
  dist <- rep(0, length(state_names))
  names(dist) <- state_names
  dist[start_rating] <- 1

  cum_prev <- 0
  out <- map_dfr(years, function(yr){
    h <- hazard_1y %>% filter(Rating == start_rating) %>% pull(h)
    direct <- dist[start_rating] * h

    dist2 <- as.numeric(dist %*% P_list[[as.character(yr)]])
    names(dist2) <- state_names

    cum <- dist2["D"]
    marg <- cum - cum_prev
    cum_prev <<- cum

    # Net migration pressure (downgrades minus upgrades) computed as FLOWS this year
    idx_start <- match(start_rating, ratings)
    P_nd <- P_list[[as.character(yr)]][ratings, ratings, drop = FALSE]

    # Expected flow mass from each current rating k into other ratings j
    # downgrades: j worse than start_rating; upgrades: j better than start_rating
    down_mass <- 0
    up_mass <- 0

    if (idx_start < length(ratings)) {
      worse_set <- ratings[(idx_start+1):length(ratings)]
      down_mass <- sum(dist[ratings] * rowSums(P_nd[, worse_set, drop = FALSE]), na.rm = TRUE)
    }
    if (idx_start > 1) {
      better_set <- ratings[1:(idx_start-1)]
      up_mass <- sum(dist[ratings] * rowSums(P_nd[, better_set, drop = FALSE]), na.rm = TRUE)
    }

    net_mig <- abs(down_mass - up_mass)

    # Proxy indirect default as net migration pressure times base hazard
    indirect <- net_mig * h

    # Keep accounting identity: cap so that direct + indirect cannot exceed total marginal
    indirect <- pmin(indirect, pmax(marg - direct, 0))

    dist <<- dist2

    tibble(Rating = start_rating, Year = yr, total_marg = marg, direct = direct, indirect = indirect, cum_default = cum)
  })
  out
}

rated_decomp <- map_dfr(ratings, markov_decomp_one_start) %>%
  mutate(indirect_share = if_else(total_marg == 0, NA_real_, indirect / total_marg))

#-------------------------
# Dynamic brackets + unrated adjustment (vectorized via joins)
#-------------------------

brackets <- get_brackets(rated_cum_long, unrated_cum)

unrated_adj <- unrated_cum %>%
  left_join(brackets, by = "Year") %>%
  left_join(
    rated_decomp %>% select(Rating, Year, indirect_share, direct),
    by = c("Low" = "Rating", "Year" = "Year")
  ) %>%
  rename(share_low = indirect_share, direct_low = direct) %>%
  left_join(
    rated_decomp %>% select(Rating, Year, indirect_share, direct),
    by = c("High" = "Rating", "Year" = "Year")
  ) %>%
  rename(share_high = indirect_share, direct_high = direct) %>%
  left_join(
    sd_long %>% select(Rating, Year, sd_offdiag),
    by = c("Low" = "Rating", "Year" = "Year")
  ) %>%
  rename(sd_low = sd_offdiag) %>%
  left_join(
    sd_long %>% select(Rating, Year, sd_offdiag),
    by = c("High" = "Rating", "Year" = "Year")
  ) %>%
  rename(sd_high = sd_offdiag) %>%
  mutate(
    share_w = w_low * share_low + w_high * share_high,
    sd_w = w_low * sd_low + w_high * sd_high,

    # Variant 1: no dispersion penalty
    added_indirect_nopen = marg * share_w,
    adj_marg_nopen = marg + added_indirect_nopen,

    # Variant 2: scale SD to direct PD (hazard) rather than an exogenous multiplier
    # We treat SD as increasing the indirect component in proportion to the direct PD level.
    sd_scaled_factor = 1 + sd_w,
    added_indirect_sds = marg * share_w * sd_scaled_factor,
    adj_marg_sds = marg + added_indirect_sds
  ) %>%
  arrange(Year) %>%
  mutate(adj_cum_nopen = cumsum(adj_marg_nopen), adj_cum_sds = cumsum(adj_marg_sds))

#-------------------------
# Outputs: table + chart + files
#-------------------------

out_table <- unrated_adj %>%
  transmute(
    Year,
    Unrated_Cum_PD = cum,
    Unrated_Marg_PD = marg,
    Bracket_Low = Low,
    Bracket_High = High,
    w_low,
    w_high,
    IndirectShare_w = share_w,
    DispersionSD_w = sd_w,
    AddedIndirect_NoPenalty = added_indirect_nopen,
    Adjusted_Marg_PD_NoPenalty = adj_marg_nopen,
    Adjusted_Cum_PD_NoPenalty = adj_cum_nopen,
    SDScaledFactor = sd_scaled_factor,
    AddedIndirect_SDScaled = added_indirect_sds,
    Adjusted_Marg_PD_SDScaled = adj_marg_sds,
    Adjusted_Cum_PD_SDScaled = adj_cum_sds
  )

print(out_table)

plot_df <- out_table %>%
  select(Year, Original = Unrated_Cum_PD, NoPenalty = Adjusted_Cum_PD_NoPenalty, SDScaled = Adjusted_Cum_PD_SDScaled) %>%
  pivot_longer(cols = c(Original, NoPenalty, SDScaled), names_to = "Series", values_to = "cum")

p <- ggplot(plot_df, aes(x = Year, y = cum * 100, color = Series)) +
  geom_line(linewidth = 1) +
  geom_point() +
  labs(
    title = "Unrated infrastructure cumulative PD: original vs adjusted (Markov Chain)",
    x = "Year",
    y = "Cumulative default probability (%)"
  ) +
  theme_minimal()

print(p)

ggsave("unrated_pd_original_vs_adjusted_variants.png", p, width = 7, height = 4, dpi = 200)

write.csv(out_table, "unrated_pd_adjustment_table_variants.csv", row.names = FALSE)
write.csv(rated_decomp, "rated_markov_decomposition.csv", row.names = FALSE)
write.csv(brackets, "unrated_brackets_by_year.csv", row.names = FALSE)
write.csv(sd_long, "transition_dispersion_sd.csv", row.names = FALSE)

write_xlsx(list(
  Unrated_Adjustment_Variants = out_table,
  Rated_Decomposition = rated_decomp,
  Brackets = brackets,
  Dispersion_SD = sd_long
), "unrated_pd_markov_outputs_variants.xlsx")
