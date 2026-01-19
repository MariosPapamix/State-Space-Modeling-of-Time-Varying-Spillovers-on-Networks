# ======================================================================
# OPTIONAL EXPERIMENT (MUMPS): Joint node–edge uncertainty with latent W_t
# - Scrape GOV.UK mumps confirmed cases by region (annual)
# - Poisson log model with time-varying coefficients (random walks)
# - Latent row-stochastic W_t (diag=0) via softmax(u_raw), u_raw random walk
# ======================================================================

# ---- 0) Packages ----
pkgs_cran <- c("rvest","dplyr","tidyr","stringr","readr","posterior","ggplot2")
for (p in pkgs_cran) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)

# cmdstanr install (not always on CRAN): try r-universe, then mc-stan repo
if (!requireNamespace("cmdstanr", quietly = TRUE)) {
  try(
    install.packages("cmdstanr",
                     repos = c("https://stan-dev.r-universe.dev", "https://cloud.r-project.org")
    ),
    silent = TRUE
  )
}
if (!requireNamespace("cmdstanr", quietly = TRUE)) {
  install.packages("cmdstanr",
                   repos = c("https://mc-stan.org/r-packages/", "https://cloud.r-project.org")
  )
}

invisible(lapply(c(pkgs_cran,"cmdstanr"), library, character.only = TRUE))

# Ensure CmdStan is installed (one-time; needs C++ toolchain)
ver <- tryCatch(cmdstanr::cmdstan_version(), error = function(e) NA_character_)
if (is.na(ver)) {
  message("CmdStan not found. Installing CmdStan (one-time)...")
  cmdstanr::install_cmdstan()
}

set.seed(1)

# ---- Helpers ----
clean_table <- function(tb) {
  tb <- as.data.frame(tb, stringsAsFactors = FALSE, check.names = FALSE)
  nms <- stringr::str_squish(names(tb))
  nms <- sub("^X(20\\d{2})$", "\\1", nms)                   # X2012 -> 2012
  nms <- vapply(nms, function(nm) {
    yr <- stringr::str_extract(nm, "20\\d{2}")
    if (!is.na(yr)) yr else nm
  }, character(1))
  names(tb) <- make.unique(nms, sep = "_")
  tb
}

detect_year_cols <- function(tb) {
  names(tb)[grepl("^20\\d{2}(_\\d+)?$", names(tb))]
}

# ---- 1) Scrape GOV.UK page ----
url  <- "https://www.gov.uk/government/publications/mumps-confirmed-cases-and-notifications/mumps-confirmed-cases-and-notifications"
page <- rvest::read_html(url)
tabs_raw <- rvest::html_elements(page, "table") |> rvest::html_table(fill = TRUE)
tabs <- lapply(tabs_raw, clean_table)

is_region_table <- function(tb) {
  if (!is.data.frame(tb) || ncol(tb) < 5) return(FALSE)
  firstcol <- tolower(stringr::str_squish(as.character(tb[[1]])))
  
  # NOTE: na.rm belongs to any(), not grepl()
  has_regions <- any(grepl("north\\s*east", firstcol), na.rm = TRUE) ||
    any(grepl("london", firstcol), na.rm = TRUE)
  
  year_cols <- detect_year_cols(tb)
  has_years <- length(year_cols) >= 5
  
  has_regions && has_years
}

scores <- vapply(tabs, function(tb) if (is_region_table(tb)) length(detect_year_cols(tb)) else 0L, integer(1))
if (max(scores) == 0L) stop("Could not find the mumps-by-region table on the GOV.UK page.")
tab_region <- tabs[[which.max(scores)]]

# ---- 2) Pivot to long format (fix mixed column types) ----
tab_region <- tab_region |> dplyr::rename(region = 1)
year_cols  <- detect_year_cols(tab_region)
if (length(year_cols) == 0) stop("Found region table but no year columns were detected.")

# critical: coerce year cols to character BEFORE pivot_longer()
tab_region <- tab_region |>
  dplyr::mutate(
    region = stringr::str_squish(as.character(region)),
    dplyr::across(dplyr::all_of(year_cols), as.character)
  )

mumps_long <- tab_region |>
  dplyr::filter(!grepl("Total|Unknown", region, ignore.case = TRUE)) |>
  tidyr::pivot_longer(
    cols = dplyr::all_of(year_cols),
    names_to = "year",
    values_to = "cases_raw"
  ) |>
  dplyr::mutate(
    year  = suppressWarnings(as.integer(stringr::str_extract(year, "20\\d{2}"))),
    cases = suppressWarnings(readr::parse_number(cases_raw)),
    cases = dplyr::coalesce(cases, 0),
    cases = pmax(cases, 0),
    cases = as.integer(round(cases))
  ) |>
  dplyr::filter(!is.na(year)) |>
  dplyr::select(region, year, cases) |>
  dplyr::arrange(year, region)

regions <- sort(unique(mumps_long$region))
years   <- sort(unique(mumps_long$year))

Y <- mumps_long |>
  tidyr::pivot_wider(names_from = region, values_from = cases, values_fill = 0) |>
  dplyr::arrange(year) |>
  dplyr::select(dplyr::all_of(regions)) |>
  as.matrix()

storage.mode(Y) <- "integer"
Tt <- nrow(Y)
N  <- ncol(Y)
cat(sprintf("Mumps data loaded: N=%d regions, T=%d years (%d–%d)\n", N, Tt, min(years), max(years)))
if (Tt < 2) stop("Need at least 2 time points.")

# ---- 3) Stan data ----
Y_lag  <- log1p(Y[1:(Tt - 1), , drop = FALSE])   # numeric
Y_curr <- Y[2:Tt,       , drop = FALSE]          # integer

stan_data <- list(
  N    = N,
  T    = Tt - 1,
  y    = array(Y_curr, dim = c(Tt - 1, N)),      # matches: array[T, N] int y;
  ylag = Y_lag                                   # matches: matrix[T, N] ylag;
)

# ---- 4) Stan model (FIXED: use array[T, N] vector[...] not nested arrays) ----
stan_code <- "
data {
  int<lower=2> N;
  int<lower=1> T;
  array[T, N] int<lower=0> y;
  matrix[T, N] ylag;  // log1p(y_{t-1})
}
parameters {
  vector[T] beta0;
  vector[T] beta1;
  vector[T] beta2;

  real<lower=0> s_b0;
  real<lower=0> s_b1;
  real<lower=0> s_b2;

  real<lower=0> s_u;

  // FIX: combined dims array[T, N], not array[T] array[N]
  array[T, N] vector[N-1] u_raw;
}
transformed parameters {
  array[T] matrix[N, N] W;

  for (t in 1:T) {
    for (i in 1:N) {
      vector[N-1] wrow;
      int k;

      wrow = softmax(u_raw[t, i] - mean(u_raw[t, i]));

      k = 1;
      for (j in 1:N) {
        if (j == i) {
          W[t][i, j] = 0;
        } else {
          W[t][i, j] = wrow[k];
          k += 1;
        }
      }
    }
  }
}
model {
  // hyperpriors
  s_b0 ~ normal(0, 0.5);
  s_b1 ~ normal(0, 0.5);
  s_b2 ~ normal(0, 0.5);
  s_u  ~ normal(0, 0.5);

  // initial states
  beta0[1] ~ normal(0, 1);
  beta1[1] ~ normal(0, 1);
  beta2[1] ~ normal(0, 1);

  // random walks for betas
  for (t in 2:T) {
    beta0[t] ~ normal(beta0[t-1], s_b0);
    beta1[t] ~ normal(beta1[t-1], s_b1);
    beta2[t] ~ normal(beta2[t-1], s_b2);
  }

  // random walk for latent network logits (vectorized over N-1 entries)
  for (i in 1:N) {
    u_raw[1, i] ~ normal(0, 1);
    for (t in 2:T) {
      u_raw[t, i] ~ normal(u_raw[t-1, i], s_u);
    }
  }

  // likelihood
  for (t in 1:T) {
    vector[N] ylag_vec;
    vector[N] netlag;
    vector[N] eta;

    ylag_vec = (ylag[t])';      // row_vector -> vector
    netlag   = W[t] * ylag_vec;

    eta = rep_vector(beta0[t], N)
          + beta1[t] * netlag
          + beta2[t] * ylag_vec;

    for (i in 1:N) {
      y[t, i] ~ poisson_log(eta[i]);
    }
  }
}
generated quantities {
  matrix[T, N] lambda_hat;
  for (t in 1:T) {
    vector[N] ylag_vec;
    vector[N] netlag;
    vector[N] eta;

    ylag_vec = (ylag[t])';
    netlag   = W[t] * ylag_vec;

    eta = rep_vector(beta0[t], N)
          + beta1[t] * netlag
          + beta2[t] * ylag_vec;

    for (i in 1:N) {
      lambda_hat[t, i] = exp(eta[i]);
    }
  }
}
"

# ---- 5) Compile + fit ----
mod <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(stan_code))

fit <- mod$sample(
  data = stan_data,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  seed = 1,
  adapt_delta = 0.9,
  max_treedepth = 12,
  refresh = 500
)

print(fit$summary(variables = c("s_b0","s_b1","s_b2","s_u")))

# ---- 6) Plot beta1_t ----
beta1_draws <- fit$draws("beta1")
beta1_mat   <- posterior::as_draws_matrix(beta1_draws)
beta1_cols  <- grep("^beta1\\[", colnames(beta1_mat), value = TRUE)

beta1_q <- apply(beta1_mat[, beta1_cols, drop = FALSE], 2, quantile, probs = c(0.1, 0.5, 0.9))

df_beta <- data.frame(
  year = years[-1],
  q10  = as.numeric(beta1_q[1, ]),
  q50  = as.numeric(beta1_q[2, ]),
  q90  = as.numeric(beta1_q[3, ])
)

p <- ggplot2::ggplot(df_beta, ggplot2::aes(x = year)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = q10, ymax = q90), alpha = 0.2) +
  ggplot2::geom_line(ggplot2::aes(y = q50), linewidth = 0.7) +
  ggplot2::labs(x = "Year", y = "beta1,t (latent-network spillover)") +
  ggplot2::theme_minimal()

print(p)


# ============================================================
# Produce "Table 4/5 style" numbers from YOUR cmdstanr fit
# + compare to the attached paper's Table 4 (static) and Table 5 (dynamic).
#
# Assumptions:
#   - You already ran your model and you have:
#       fit    : CmdStanMCMC object
#       years  : integer vector of observation years (length Tt)
#       regions: character vector of region names (length N)
#   - Your Stan model has parameters beta0[1:T], beta1[1:T], beta2[1:T]
#     and hyperparameters s_b0, s_b1, s_b2, s_u.
#
# If you don't have `years`/`regions` in your environment, set them from
# your data prep objects (e.g., years <- sort(unique(mumps_long$year))).
# ============================================================

pkgs <- c("dplyr","tibble","tidyr","stringr","posterior","cmdstanr","knitr")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if(length(to_install)) install.packages(to_install)
invisible(lapply(pkgs, library, character.only = TRUE))

# ---------- Helper formatting ----------
fmt_ci <- function(mean, l, u, digits = 2) {
  sprintf(paste0("%.",digits,"f (%.",digits,"f, %.",digits,"f)"), mean, l, u)
}

summ_1d <- function(x, probs = c(0.025, 0.975)) {
  x <- as.numeric(x)
  c(mean = mean(x), l = unname(quantile(x, probs[1])), u = unname(quantile(x, probs[2])))
}

summ_vars_from_fit <- function(fit, vars, probs = c(0.025, 0.975), digits = 2) {
  # vars: character vector of variable names (e.g. c("s_b0","s_b1"))
  out <- lapply(vars, function(v) {
    d <- fit$draws(v)
    m <- posterior::as_draws_matrix(d)[, 1]
    s <- summ_1d(m, probs = probs)
    tibble::tibble(
      Parameter = v,
      Mean = s["mean"],
      L = s["l"],
      U = s["u"],
      `95% credible interval` = fmt_ci(s["mean"], s["l"], s["u"], digits = digits)
    )
  })
  dplyr::bind_rows(out)
}

# ============================================================
# 1) YOUR PAPER: main posterior summaries (Table-style)
# ============================================================

# --- Hyperparameters (your paper / your model) ---
tbl_hyp <- summ_vars_from_fit(
  fit,
  vars = c("s_b0","s_b1","s_b2","s_u"),
  probs = c(0.025, 0.975),
  digits = 2
)

# --- Time-averaged betas (to make a compact "Table 5 style" summary) ---
#     Because beta0,beta1,beta2 are time-varying in your model, we summarise
#     their posterior distribution of the TIME AVERAGE across t.
beta_draws <- posterior::as_draws_matrix(fit$draws(c("beta0","beta1","beta2")))

beta0_cols <- grep("^beta0\\[", colnames(beta_draws), value = TRUE)
beta1_cols <- grep("^beta1\\[", colnames(beta_draws), value = TRUE)
beta2_cols <- grep("^beta2\\[", colnames(beta_draws), value = TRUE)

beta0_bar <- rowMeans(beta_draws[, beta0_cols, drop = FALSE])
beta1_bar <- rowMeans(beta_draws[, beta1_cols, drop = FALSE])
beta2_bar <- rowMeans(beta_draws[, beta2_cols, drop = FALSE])

s0 <- summ_1d(beta0_bar); s1 <- summ_1d(beta1_bar); s2 <- summ_1d(beta2_bar)

tbl_bars <- tibble::tibble(
  Parameter = c("beta0_bar","beta1_bar","beta2_bar"),
  Mean = c(s0["mean"], s1["mean"], s2["mean"]),
  L    = c(s0["l"],    s1["l"],    s2["l"]),
  U    = c(s0["u"],    s1["u"],    s2["u"])
) %>%
  dplyr::mutate(`95% credible interval` = fmt_ci(Mean, L, U, digits = 2)) %>%
  dplyr::select(Parameter, Mean, `95% credible interval`)

# --- Combine into ONE table (like "numbers in the paper") ---
tbl_yourpaper <- dplyr::bind_rows(
  tbl_bars,
  tbl_hyp %>% dplyr::select(Parameter, Mean, `95% credible interval`)
)

# Print as LaTeX-ready table
knitr::kable(
  tbl_yourpaper,
  format = "latex",
  booktabs = TRUE,
  digits = 2,
  caption = "Your model (latent W_t PNAR): posterior means with 95\\% credible intervals."
)

# ============================================================
# 2) ATTACHED PAPER: Table 4 (static) and Table 5 (dynamic)
#    (Hard-coded from the PDF so you can compare.)
# ============================================================

# Table 4 (static model) from the attached paper:
tbl_attached_T4 <- tibble::tribble(
  ~Parameter, ~Mean,  ~L,    ~U,
  "alpha",     2.18,  2.12,  2.23,
  "beta_EM",   0.27,  0.22,  0.33,
  "beta_EE",   0.32,  0.16,  0.46,
  "beta_LO",   0.34,  0.28,  0.41,
  "beta_NE",   0.47,  0.43,  0.51,
  "beta_NW",   0.51,  0.48,  0.56,
  "beta_SE",   0.40,  0.33,  0.46,
  "beta_WM",   0.32,  0.25,  0.39,
  "beta_YH",   0.18,  0.13,  0.24,
  "beta_SW",   0.29,  0.24,  0.33
) %>%
  dplyr::mutate(`95% credible interval` = fmt_ci(Mean, L, U, digits = 2)) %>%
  dplyr::select(Parameter, Mean, `95% credible interval`)

# Table 5 (dynamic model) from the attached paper:
tbl_attached_T5 <- tibble::tribble(
  ~Parameter, ~Mean,   ~L,     ~U,
  "alpha",     1.45,   1.02,   1.83,
  "beta_EM",  -0.31,  -0.49,  -0.13,
  "beta_EE",  -0.60,  -0.79,  -0.17,
  "beta_LO",  -0.21,  -0.34,  -0.08,
  "beta_NE",  -0.12,  -0.25,   0.03,
  "beta_NW",   0.02,  -0.09,   0.14,
  "beta_SE",  -0.03,  -0.14,   0.07,
  "beta_WM",  -0.46,  -0.66,  -0.12,
  "beta_YH",  -0.03,  -0.16,   0.10,
  "beta_SW",  -0.22,  -0.34,  -0.08
) %>%
  dplyr::mutate(`95% credible interval` = fmt_ci(Mean, L, U, digits = 2)) %>%
  dplyr::select(Parameter, Mean, `95% credible interval`)

knitr::kable(tbl_attached_T4, format="latex", booktabs=TRUE, digits=2,
             caption="Attached paper: Table 4 (static) posterior means with 95\\% credible intervals.")
knitr::kable(tbl_attached_T5, format="latex", booktabs=TRUE, digits=2,
             caption="Attached paper: Table 5 (dynamic) posterior means with 95\\% credible intervals.")

# ============================================================
# 3) COMPARISON TABLES
#    3a) Intercept comparison: your beta0_bar vs their alpha
#    3b) Own-lag comparison: your beta2_bar (global) vs their beta_i (by region)
# ============================================================

# Pull your beta0_bar and beta2_bar summary numbers as scalars:
your_beta0_bar <- tibble::tibble(
  Parameter = "beta0_bar",
  Mean = s0["mean"], L = s0["l"], U = s0["u"],
  `95% credible interval` = fmt_ci(s0["mean"], s0["l"], s0["u"], digits=2)
)

your_beta2_bar <- tibble::tibble(
  Parameter = "beta2_bar",
  Mean = s2["mean"], L = s2["l"], U = s2["u"],
  `95% credible interval` = fmt_ci(s2["mean"], s2["l"], s2["u"], digits=2)
)

# --- (3a) alpha vs beta0_bar ---
cmp_intercept <- tibble::tibble(
  Quantity = c("Intercept (your model)","Intercept (attached T4)","Intercept (attached T5)"),
  Parameter = c("beta0_bar","alpha","alpha"),
  Mean = c(your_beta0_bar$Mean, 2.18, 1.45),
  CrI = c(your_beta0_bar$`95% credible interval`, "2.18 (2.12, 2.23)", "1.45 (1.02, 1.83)")
)

knitr::kable(cmp_intercept, format="latex", booktabs=TRUE, digits=2,
             caption="Intercept comparison: your time-averaged $\\beta_{0,t}$ vs the attached paper's $\\alpha$.")

# --- (3b) region-level betas: your beta2_bar replicated vs their beta_i ---
# Keep only beta_* rows (exclude alpha)
t4_betas <- tbl_attached_T4 %>%
  dplyr::filter(grepl("^beta_", Parameter)) %>%
  dplyr::select(Parameter, Mean_T4 = Mean, CrI_T4 = `95% credible interval`)

t5_betas <- tbl_attached_T5 %>%
  dplyr::filter(grepl("^beta_", Parameter)) %>%
  dplyr::select(Parameter, Mean_T5 = Mean, CrI_T5 = `95% credible interval`)

cmp_betas <- dplyr::full_join(t4_betas, t5_betas, by = "Parameter") %>%
  dplyr::mutate(
    Your_beta2_bar_Mean = as.numeric(your_beta2_bar$Mean),
    Your_beta2_bar_CrI  = as.character(your_beta2_bar$`95% credible interval`),
    Diff_T4_minus_yours = Mean_T4 - Your_beta2_bar_Mean,
    Diff_T5_minus_yours = Mean_T5 - Your_beta2_bar_Mean
  ) %>%
  dplyr::arrange(Parameter)

knitr::kable(cmp_betas, format="latex", booktabs=TRUE, digits=2,
             caption=paste0("Own-lag comparison: attached paper's region-specific $\\beta_i$ ",
                            "vs your model's global time-averaged $\\overline{\\beta_2}$. ",
                            "Differences shown as (attached minus yours)."))

# ============================================================
# 4) OPTIONAL: Numbers behind your beta1_t plot (by year)
#    (If you want a numeric table to accompany the figure)
# ============================================================

# beta1 per time index -> per calendar year (years[-1] corresponds to t=1..T in your Stan data)
beta1_mat <- beta_draws[, beta1_cols, drop = FALSE]
beta1_mean <- apply(beta1_mat, 2, mean)
beta1_l    <- apply(beta1_mat, 2, quantile, probs = 0.025) %>% unname()
beta1_u    <- apply(beta1_mat, 2, quantile, probs = 0.975) %>% unname()

tbl_beta1_by_year <- tibble::tibble(
  Year = years[-1],
  Mean = as.numeric(beta1_mean),
  L    = as.numeric(beta1_l),
  U    = as.numeric(beta1_u)
) %>%
  dplyr::mutate(`95% credible interval` = fmt_ci(Mean, L, U, digits = 2)) %>%
  dplyr::select(Year, Mean, `95% credible interval`)

knitr::kable(tbl_beta1_by_year, format="latex", booktabs=TRUE, digits=2,
             caption="Your model: yearly posterior summaries for the spillover coefficient $\\beta_{1,t}$.")




















stan_code_nc <- "
data {
  int<lower=2> N;
  int<lower=1> T;
  array[T, N] int<lower=0> y;
  matrix[T, N] ylag;
}
parameters {
  real beta0_1;
  real beta1_1;
  real beta2_1;

  vector[T-1] z_b0;
  vector[T-1] z_b1;
  vector[T-1] z_b2;

  real<lower=0> s_b0;
  real<lower=0> s_b1;
  real<lower=0> s_b2;

  real<lower=0> s_u;

  // Non-centered RW increments for logits
  array[N] matrix[T, N-1] z_u;   // per row i: T x (N-1) std normals
}
transformed parameters {
  vector[T] beta0;
  vector[T] beta1;
  vector[T] beta2;

  array[T] matrix[N,N] W;

  beta0[1] = beta0_1;
  beta1[1] = beta1_1;
  beta2[1] = beta2_1;
  for (t in 2:T) {
    beta0[t] = beta0[t-1] + s_b0 * z_b0[t-1];
    beta1[t] = beta1[t-1] + s_b1 * z_b1[t-1];
    beta2[t] = beta2[t-1] + s_b2 * z_b2[t-1];
  }

  for (t in 1:T) {
    for (i in 1:N) {
      vector[N-1] u;
      vector[N-1] wrow;
      int k;
      // build non-centered RW logits: u_t = s_u * cumsum(z_u)
      // for simplicity: treat z_u as innovations with cumulative sum in time
      // u(t) = s_u * sum_{r<=t} z_u(r)
      for (j in 1:(N-1)) {
        real acc = 0;
        for (r in 1:t) acc += z_u[i][r,j];
        u[j] = s_u * acc;
      }
      wrow = softmax(u - mean(u));
      k = 1;
      for (j in 1:N) {
        if (j==i) W[t][i,j] = 0;
        else { W[t][i,j] = wrow[k]; k += 1; }
      }
    }
  }
}
model {
  // stronger regularization (key for short annual panels)
  s_b0 ~ normal(0, 0.3);
  s_b1 ~ normal(0, 0.3);
  s_b2 ~ normal(0, 0.3);
  s_u  ~ normal(0, 0.3);

  beta0_1 ~ normal(0, 1);
  beta1_1 ~ normal(0, 1);
  beta2_1 ~ normal(0, 1);

  z_b0 ~ std_normal();
  z_b1 ~ std_normal();
  z_b2 ~ std_normal();

  for (i in 1:N) to_vector(z_u[i]) ~ std_normal();

  for (t in 1:T) {
    vector[N] ylag_vec = (ylag[t])';
    vector[N] netlag   = W[t] * ylag_vec;
    vector[N] eta = rep_vector(beta0[t], N) + beta1[t]*netlag + beta2[t]*ylag_vec;
    for (i in 1:N) y[t,i] ~ poisson_log(eta[i]);
  }
}
generated quantities {
  // one-step predictive lambda for scoring
  matrix[T, N] lambda_hat;
  for (t in 1:T) {
    vector[N] ylag_vec = (ylag[t])';
    vector[N] netlag   = W[t] * ylag_vec;
    vector[N] eta = rep_vector(beta0[t], N) + beta1[t]*netlag + beta2[t]*ylag_vec;
    for (i in 1:N) lambda_hat[t,i] = exp(eta[i]);
  }
}
"


pkgs <- c("cmdstanr","posterior","dplyr","tibble")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if(length(to_install)) install.packages(to_install)
invisible(lapply(pkgs, library, character.only = TRUE))

# Compile latent-W model (non-centered)
mod_nc <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(stan_code_nc))

fit_nc <- mod_nc$sample(
  data = stan_data,
  chains = 4, parallel_chains = 4,
  iter_warmup = 1500, iter_sampling = 1500,
  adapt_delta = 0.95, max_treedepth = 14,
  seed = 1, refresh = 500
)

print(fit_nc$summary(variables=c("s_b0","s_b1","s_b2","s_u")))

# ---- baseline: no-network dynamic Poisson RW (beta1_t == 0, no W) ----
stan_code_baseline <- "
data {
  int<lower=2> N;
  int<lower=1> T;
  array[T, N] int<lower=0> y;
  matrix[T, N] ylag;
}
parameters {
  real beta0_1;
  real beta2_1;
  vector[T-1] z_b0;
  vector[T-1] z_b2;
  real<lower=0> s_b0;
  real<lower=0> s_b2;
}
transformed parameters {
  vector[T] beta0;
  vector[T] beta2;
  beta0[1] = beta0_1;
  beta2[1] = beta2_1;
  for (t in 2:T) {
    beta0[t] = beta0[t-1] + s_b0*z_b0[t-1];
    beta2[t] = beta2[t-1] + s_b2*z_b2[t-1];
  }
}
model {
  s_b0 ~ normal(0, 0.3);
  s_b2 ~ normal(0, 0.3);
  beta0_1 ~ normal(0, 1);
  beta2_1 ~ normal(0, 1);
  z_b0 ~ std_normal();
  z_b2 ~ std_normal();

  for (t in 1:T) {
    vector[N] ylag_vec = (ylag[t])';
    vector[N] eta = rep_vector(beta0[t], N) + beta2[t]*ylag_vec;
    for (i in 1:N) y[t,i] ~ poisson_log(eta[i]);
  }
}
generated quantities {
  matrix[T, N] lambda_hat;
  for (t in 1:T) {
    vector[N] ylag_vec = (ylag[t])';
    vector[N] eta = rep_vector(beta0[t], N) + beta2[t]*ylag_vec;
    for (i in 1:N) lambda_hat[t,i] = exp(eta[i]);
  }
}
"
mod_base <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(stan_code_baseline))

fit_base <- mod_base$sample(
  data = stan_data,
  chains = 4, parallel_chains = 4,
  iter_warmup = 1500, iter_sampling = 1500,
  adapt_delta = 0.95, max_treedepth = 14,
  seed = 2, refresh = 500
)

# ---- rolling predictive log score proxy (one-step, in-sample years) ----
# Use posterior mean lambda_hat[t,i] as plug-in for scoring; for stronger scoring,
# average dpois over draws (mixture), but this is already a clean head-to-head.

lambda_nc <- posterior::as_draws_matrix(fit_nc$draws("lambda_hat"))
lambda_ba <- posterior::as_draws_matrix(fit_base$draws("lambda_hat"))

# posterior mean lambda by t,i:
get_lambda_mean <- function(lambda_draws, T, N) {
  # columns like lambda_hat[1,1] ...
  out <- matrix(NA_real_, nrow=T, ncol=N)
  for (t in 1:T) for (i in 1:N) {
    nm <- sprintf("lambda_hat[%d,%d]", t, i)
    out[t,i] <- mean(lambda_draws[, nm])
  }
  out
}

Tstan <- stan_data$T
Nstan <- stan_data$N
lam_nc_mean <- get_lambda_mean(lambda_nc, Tstan, Nstan)
lam_ba_mean <- get_lambda_mean(lambda_ba, Tstan, Nstan)

y_obs <- array(stan_data$y, dim=c(Tstan,Nstan))
ls_nc <- sapply(1:Tstan, function(t) sum(dpois(y_obs[t,], pmax(lam_nc_mean[t,],1e-12), log=TRUE)))
ls_ba <- sapply(1:Tstan, function(t) sum(dpois(y_obs[t,], pmax(lam_ba_mean[t,],1e-12), log=TRUE)))

tbl_mumps_score <- tibble::tibble(
  year = years[-1],  # corresponds to t=1..T
  logscore_latentW = ls_nc,
  logscore_nonet   = ls_ba,
  d_logscore       = ls_nc - ls_ba
)
print(tbl_mumps_score)
cat("Mean d_logscore (latentW - nonet):", mean(tbl_mumps_score$d_logscore), "\n")


