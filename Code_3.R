# =============================================================================
# EMPIRICAL II: Chicago burglary counts (Poisson NSSM) — FULL SCRIPT 
# Data: nick3703/Chicago-Data
#   - crime.csv is ALREADY MONTHLY counts (N=552 locations x T=72 months)
#   - neighborhood.mtx is adjacency (552x552, MatrixMarket)
#
# Implements:
#   (1) Poisson NSSM network TVP (EKF/DGLM RW filter)
#   (2) Poisson TVP ablation with no network (beta1,t = 0)
#   (3) Static PNAR baseline (constant spillover) via rolling-origin GLM
#   (4) Static seasonal Poisson regression (region FE + month-of-year) via rolling-origin GLM
#
# Reports:
#   rolling-origin MAE, mean log score, and 90% PI coverage on last test_h months.
#
# =============================================================================

# ---------- 0) Packages ----------
pkgs <- c("readr","dplyr","Matrix","zoo","ggplot2","tibble")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install)) install.packages(to_install, dependencies = TRUE)
invisible(lapply(pkgs, require, character.only = TRUE))

options(timeout = 600)

# ---------- 1) Download data ----------
base_url <- "https://raw.githubusercontent.com/nick3703/Chicago-Data/master"
crime_url <- paste0(base_url, "/crime.csv")
W_url     <- paste0(base_url, "/neighborhood.mtx")

cache_dir <- "nssm_cache"
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

crime_path <- file.path(cache_dir, "crime.csv")
W_path     <- file.path(cache_dir, "neighborhood.mtx")

if (!file.exists(crime_path)) download.file(crime_url, crime_path, mode="wb", quiet=TRUE)
if (!file.exists(W_path))     download.file(W_url,     W_path,     mode="wb", quiet=TRUE)

# ---------- 2) Read monthly counts (N x T in file) ----------
crime_df <- readr::read_csv(crime_path, show_col_types = FALSE)

# Keep only monthly columns count.YYYYMM
count_cols <- grep("^count\\.[0-9]{6}$", names(crime_df), value = TRUE)
if (length(count_cols) == 0) {
  stop("Did not find any columns like count.YYYYMM in crime.csv. Check file format.")
}
count_cols <- count_cols[order(count_cols)]  # chronological

Y_loc_month <- as.matrix(crime_df[, count_cols, drop=FALSE])
storage.mode(Y_loc_month) <- "numeric"

N  <- nrow(Y_loc_month)
Tt <- ncol(Y_loc_month)

if (N  != 552) message("Note: N=", N,  " (expected 552).")
if (Tt != 72)  message("Note: T=", Tt, " (expected 72).")

# Convert to T x N (time x region)
Y <- t(Y_loc_month)

# Build monthly index from column names
yyyymm <- sub("^count\\.", "", count_cols)
dates  <- zoo::as.yearmon(
  paste0(substr(yyyymm, 1, 4), "-", substr(yyyymm, 5, 6)),
  "%Y-%m"
)

Yz <- zoo::zoo(Y, order.by = dates)

cat(sprintf("Loaded burglary panel: N=%d regions, T=%d months (%s .. %s)\n",
            ncol(Yz), nrow(Yz), format(start(Yz)), format(end(Yz))))

# ---------- 3) Load adjacency and build row-stochastic W  ----------
A0 <- Matrix::readMM(W_path)

# IMPORTANT:
#   Do NOT call Matrix::as() (it's not exported).
#   Use methods::as() / as() instead.
#
# readMM() may return a symmetric triplet sparse class (e.g. dsTMatrix).
# Safest conversion chain:
#   1) to CsparseMatrix (column-compressed; may still be symmetric)
#   2) to generalMatrix (drops symmetric wrapper if present)
#   3) to dgCMatrix (numeric general sparse)
A0 <- methods::as(A0, "CsparseMatrix")
A0 <- methods::as(A0, "generalMatrix")
A  <- methods::as(A0, "dgCMatrix")

# Sanity checks
if (nrow(A) != ncol(A)) stop("Adjacency matrix is not square.")
if (nrow(A) != N) stop(sprintf("Adjacency dimension (%d) != N in panel (%d).", nrow(A), N))

# Remove self-loops
Matrix::diag(A) <- 0
A <- Matrix::drop0(A)

# Row-stochastic normalization (each non-isolate row sums to 1)
rs <- Matrix::rowSums(A)
inv_rs <- ifelse(rs > 0, 1/rs, 0)   # isolates keep all-zero row
W <- Matrix::Diagonal(x = inv_rs) %*% A
W <- Matrix::drop0(W)

cat("Adjacency nnz:", length(A@x),
    " | Isolates:", sum(rs == 0),
    " | Max |rowSum(W)-1| (non-isolates):",
    ifelse(any(rs > 0), max(abs(Matrix::rowSums(W)[rs > 0] - 1)), NA_real_),
    "\n")

# ---------- 4) Build Poisson NSSM regressors ----------
# Observation model (per time t, region i):
# y_{i,t} ~ Pois(exp(beta0_t + beta1_t * (W log(1+y_{t-1}))_i + beta2_t * log(1+y_{i,t-1})))
Y_mat  <- zoo::coredata(Yz)              # T x N numeric
Y_lag  <- Y_mat[1:(Tt-1), , drop=FALSE]
Y_curr <- Y_mat[2:Tt,   , drop=FALSE]
dates_curr <- dates[2:Tt]

X_list <- vector("list", Tt-1)
for (t in 1:(Tt-1)) {
  ylag  <- log1p(Y_lag[t, ])
  x_net <- as.numeric(W %*% ylag)
  X_list[[t]] <- cbind(1, x_net, ylag)   # N x 3
}

# no-network ablation: beta1,t = 0 => design [1, selflag]
X_list_nonet <- lapply(1:(Tt-1), function(t) {
  ylag <- log1p(Y_lag[t, ])
  cbind(1, ylag)                         # N x 2
})

# ---------- 5) EKF / DGLM-style RW filter for Poisson(log) ----------
ekf_rw_poisson <- function(Y, X_list, Q_list, m0, C0, clip_eta = 20) {
  TT <- nrow(Y); N <- ncol(Y); K <- length(m0)
  stopifnot(length(X_list) == TT, length(Q_list) == TT)
  
  m_pred   <- matrix(NA_real_, TT, K)
  m_filt   <- matrix(NA_real_, TT, K)
  C_filt   <- array(NA_real_, dim = c(TT, K, K))
  yhat     <- matrix(NA_real_, TT, N)
  logscore <- numeric(TT)
  
  m_prev <- m0
  C_prev <- C0
  
  for (t in 1:TT) {
    # predict
    m_t_pred <- m_prev
    C_t_pred <- C_prev + Q_list[[t]]
    
    X_t <- X_list[[t]]
    y_t <- as.numeric(Y[t, ])
    
    # one-step mean
    eta0 <- as.numeric(X_t %*% m_t_pred)
    eta0 <- pmax(pmin(eta0, clip_eta), -clip_eta)
    mu0  <- exp(eta0)
    
    yhat[t, ] <- mu0
    logscore[t] <- sum(dpois(y_t, lambda = pmax(mu0, 1e-12), log = TRUE))
    
    # IRLS pseudo-observation for Poisson log link
    z <- eta0 + (y_t - mu0) / pmax(mu0, 1e-12)
    w <- pmax(mu0, 1e-12)
    
    sw <- sqrt(w)
    y_star <- sw * z
    X_star <- X_t * sw
    
    # Gaussian update in pseudo-model
    C_post <- solve(solve(C_t_pred) + crossprod(X_star))
    resid  <- y_star - as.numeric(X_star %*% m_t_pred)
    m_t_filt <- m_t_pred + as.numeric(C_post %*% crossprod(X_star, resid))
    
    m_pred[t, ] <- m_t_pred
    m_filt[t, ] <- m_t_filt
    C_filt[t, , ] <- C_post
    
    m_prev <- m_t_filt
    C_prev <- C_post
  }
  
  list(m_pred = m_pred, m_filt = m_filt, C_filt = C_filt, yhat = yhat, logscore = logscore)
}

# ---------- 6) Fit dynamic Poisson NSSM + no-network TVP ----------
Q_full  <- replicate(Tt-1, diag(c(5e-3, 5e-3, 5e-3)), simplify = FALSE) # small RW variance
Q_nonet <- replicate(Tt-1, diag(c(5e-3, 5e-3)),       simplify = FALSE)

fit_pois_net   <- ekf_rw_poisson(Y_curr, X_list,       Q_full,  m0 = c(0, 0, 0), C0 = diag(1, 3))
fit_pois_nonet <- ekf_rw_poisson(Y_curr, X_list_nonet, Q_nonet, m0 = c(0, 0),    C0 = diag(1, 2))

# ---------- 7) Rolling-origin baselines (static models refit each test month) ----------
test_h   <- 12
idx_test <- (nrow(Y_curr) - test_h + 1):nrow(Y_curr)
stopifnot(length(idx_test) == test_h)

poisson_metrics <- function(Y_true, mu_hat, idx_test) {
  mu_hat <- pmax(mu_hat, 1e-12)
  
  MAE <- mean(abs(Y_true[idx_test, , drop=FALSE] - mu_hat[idx_test, , drop=FALSE]))
  
  ls_by_t <- sapply(idx_test, function(t)
    sum(dpois(Y_true[t, ], lambda = mu_hat[t, ], log = TRUE))
  )
  mean_logscore <- mean(ls_by_t)
  
  low <- qpois(0.05, lambda = mu_hat[idx_test, , drop=FALSE])
  hi  <- qpois(0.95, lambda = mu_hat[idx_test, , drop=FALSE])
  Yt  <- Y_true[idx_test, , drop=FALSE]
  cover90 <- mean((Yt >= low) & (Yt <= hi))
  
  list(MAE = MAE, mean_logscore = mean_logscore, cover90 = cover90)
}

# --- (i) Static PNAR baseline: constant coefficients on [1, netlag, selflag]
yhat_static_pnar <- matrix(NA_real_, nrow = nrow(Y_curr), ncol = ncol(Y_curr))
logscore_static_pnar <- rep(NA_real_, nrow(Y_curr))

for (k in idx_test) {
  X_train <- do.call(rbind, X_list[1:(k-1)])
  y_train <- as.vector(t(Y_curr[1:(k-1), , drop=FALSE]))
  
  df_train <- data.frame(netlag = X_train[,2], selflag = X_train[,3], y = y_train)
  glm_fit <- glm(y ~ netlag + selflag, family = poisson(), data = df_train)
  
  Xk <- X_list[[k]]
  df_test <- data.frame(netlag = Xk[,2], selflag = Xk[,3])
  mu <- as.numeric(predict(glm_fit, newdata = df_test, type = "response"))
  
  yhat_static_pnar[k, ] <- mu
  logscore_static_pnar[k] <- sum(dpois(Y_curr[k, ], lambda = pmax(mu, 1e-12), log = TRUE))
}

# --- (ii) Static no-network Poisson autoregression baseline: constant coefficients on [1, selflag]
yhat_static_nonet <- matrix(NA_real_, nrow = nrow(Y_curr), ncol = ncol(Y_curr))
logscore_static_nonet <- rep(NA_real_, nrow(Y_curr))

for (k in idx_test) {
  X_train <- do.call(rbind, X_list_nonet[1:(k-1)])
  y_train <- as.vector(t(Y_curr[1:(k-1), , drop=FALSE]))
  
  df_train <- data.frame(selflag = X_train[,2], y = y_train)
  glm_fit <- glm(y ~ selflag, family = poisson(), data = df_train)
  
  Xk <- X_list_nonet[[k]]
  df_test <- data.frame(selflag = Xk[,2])
  mu <- as.numeric(predict(glm_fit, newdata = df_test, type = "response"))
  
  yhat_static_nonet[k, ] <- mu
  logscore_static_nonet[k] <- sum(dpois(Y_curr[k, ], lambda = pmax(mu, 1e-12), log = TRUE))
}

# --- (iii) Static Poisson regression with seasonal covariates (region FE + month-of-year)
month_fac <- factor(format(dates, "%m"), levels = sprintf("%02d", 1:12))
month_fac_curr <- month_fac[2:Tt]

yhat_static_season <- matrix(NA_real_, nrow = nrow(Y_curr), ncol = ncol(Y_curr))
logscore_static_season <- rep(NA_real_, nrow(Y_curr))

for (k in idx_test) {
  Y_train <- Y_curr[1:(k-1), , drop=FALSE]
  y_train <- as.vector(t(Y_train))
  
  df_train <- data.frame(
    y      = y_train,
    region = factor(rep(1:N, times = (k-1))),
    month  = factor(rep(month_fac_curr[1:(k-1)], each = N), levels = levels(month_fac_curr))
  )
  
  glm_fit <- glm(y ~ 0 + region + month, family = poisson(), data = df_train)
  
  df_test <- data.frame(
    region = factor(1:N, levels = levels(df_train$region)),
    month  = factor(rep(month_fac_curr[k], N), levels = levels(df_train$month))
  )
  
  mu <- as.numeric(predict(glm_fit, newdata = df_test, type = "response"))
  
  yhat_static_season[k, ] <- mu
  logscore_static_season[k] <- sum(dpois(Y_curr[k, ], lambda = pmax(mu, 1e-12), log = TRUE))
}

# ---------- 8) Summaries on the last test_h months ----------
met_net  <- poisson_metrics(Y_curr, fit_pois_net$yhat, idx_test)
met_tvp0 <- poisson_metrics(Y_curr, fit_pois_nonet$yhat, idx_test)
met_pnar <- poisson_metrics(Y_curr, yhat_static_pnar, idx_test)
met_s0   <- poisson_metrics(Y_curr, yhat_static_nonet, idx_test)
met_seas <- poisson_metrics(Y_curr, yhat_static_season, idx_test)

res_tbl <- tibble::tibble(
  model = c(
    "Poisson NSSM (network TVP)",
    "Poisson TVP (no network)",
    "Static PNAR (const spillover, rolling)",
    "Static Poisson AR (no network, rolling)",
    "Static Poisson seasonal (region FE + month, rolling)"
  ),
  test_MAE = c(met_net$MAE, met_tvp0$MAE, met_pnar$MAE, met_s0$MAE, met_seas$MAE),
  test_mean_logscore = c(
    met_net$mean_logscore,
    met_tvp0$mean_logscore,
    mean(logscore_static_pnar[idx_test],  na.rm = TRUE),
    mean(logscore_static_nonet[idx_test], na.rm = TRUE),
    mean(logscore_static_season[idx_test], na.rm = TRUE)
  ),
  test_cover90 = c(met_net$cover90, met_tvp0$cover90, met_pnar$cover90, met_s0$cover90, met_seas$cover90)
)

print(res_tbl)

# ---------- 9) Plot spillover intensity beta1,t with approx bands ----------
df_beta <- tibble::tibble(
  date = dates_curr,
  beta1 = fit_pois_net$m_filt[, 2],
  se1   = sqrt(pmax(fit_pois_net$C_filt[, 2, 2], 0))
) |>
  dplyr::mutate(lo = beta1 - 1.96 * se1, hi = beta1 + 1.96 * se1)

p_beta1 <- ggplot(df_beta, aes(x = date, y = beta1)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2) +
  geom_line(linewidth = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Month", y = expression(hat(beta)[1*t] ~ "(network spillover)")) +
  theme_minimal()

print(p_beta1)

# ---------- 10) Save outputs ----------
dir.create("figs", showWarnings = FALSE, recursive = TRUE)
ggsave("figs/crime_beta1_poisson_nssm.pdf", p_beta1, width = 7, height = 4)
write.csv(res_tbl, "figs/crime_model_comparison.csv", row.names = FALSE)

cat("Saved: figs/crime_beta1_poisson_nssm.pdf and figs/crime_model_comparison.csv\n")

















# <-- paste the PATCH block from section 0 right here -->

# ----------------------------
# Align objects
# ----------------------------
Y_full <- as.matrix(Y_mat)

# If dates length matches columns, transpose
if (exists("dates_full") && length(dates_full) == ncol(Y_full) && length(dates_full) != nrow(Y_full)) {
  Y_full <- t(Y_full)
}

Tt <- nrow(Y_full)
N  <- ncol(Y_full)

W_full <- W  # expected N x N dgCMatrix / matrix
stopifnot(all(dim(W_full) == c(N, N)))

theta_net_filt <- as.matrix(fit_pois_net$m_filt)     # (Tt-1)x3 or Tt x 3
theta_non_filt <- as.matrix(fit_pois_nonet$m_filt)   # (Tt-1)x2 or Tt x 2

theta_at <- function(t, theta_mat) {
  nr <- nrow(theta_mat)
  if (nr == Tt)     return(theta_mat[t, , drop=TRUE])
  if (nr == Tt - 1) return(theta_mat[t - 1, , drop=TRUE])  # row1 corresponds to t=2
  stop("theta_mat nrow must be Tt or Tt-1 (got ", nr, ").")
}

make_pairs <- function(idx_target, h, t_origin_min = 2L, t_target_max = Tt) {
  idx_target <- idx_target[idx_target > h & idx_target <= t_target_max]
  idx_origin <- idx_target - h
  keep <- idx_origin >= t_origin_min
  idx_target <- idx_target[keep]
  idx_origin <- idx_origin[keep]
  if (length(idx_target) == 0L) stop("No valid origin/target pairs for h=", h)
  list(origin = idx_origin, target = idx_target)
}

# ----------------------------
# A) Multi-step horizons h = 1,2,4,8
# ----------------------------
pois_fc_h_mean <- function(t_origin, h, th, W_use) {
  b0 <- th[1]; b1 <- th[2]; b2 <- th[3]
  y_prev <- as.numeric(Y_full[t_origin, ])
  
  lam <- rep(NA_real_, length(y_prev))
  for (j in seq_len(h)) {
    x_self <- log1p(pmax(y_prev, 0))
    x_net  <- as.numeric(W_use %*% x_self)
    eta <- b0 + b1 * x_net + b2 * x_self
    lam <- exp(pmax(pmin(eta, 20), -20))
    y_prev <- lam
  }
  lam
}

pois_fc_h_mean_nonet <- function(t_origin, h, th) {
  b0 <- th[1]; b2 <- th[2]
  y_prev <- as.numeric(Y_full[t_origin, ])
  
  lam <- rep(NA_real_, length(y_prev))
  for (j in seq_len(h)) {
    x_self <- log1p(pmax(y_prev, 0))
    eta <- b0 + b2 * x_self
    lam <- exp(pmax(pmin(eta, 20), -20))
    y_prev <- lam
  }
  lam
}

# test targets = last test_h months (in Y_full indexing)
idx_target_test <- (Tt - test_h + 1):Tt
h_grid <- c(1,2,4,8)

eval_pois_h <- function(h, W_use, labelW = "oracleW") {
  pr <- make_pairs(idx_target_test, h, t_origin_min = 2L, t_target_max = Tt)
  idx_origin <- pr$origin
  idx_target <- pr$target
  
  y_true <- Y_full[idx_target, , drop=FALSE]
  
  lam_net <- matrix(NA_real_, nrow=length(idx_origin), ncol=N)
  lam_non <- matrix(NA_real_, nrow=length(idx_origin), ncol=N)
  
  for (k in seq_along(idx_origin)) {
    t0 <- idx_origin[k]
    lam_net[k,] <- pois_fc_h_mean(t0, h, theta_at(t0, theta_net_filt), W_use)
    lam_non[k,] <- pois_fc_h_mean_nonet(t0, h, theta_at(t0, theta_non_filt))
  }
  
  mae_net_t <- rowMeans(abs(y_true - lam_net))
  mae_non_t <- rowMeans(abs(y_true - lam_non))
  
  ls_net_t <- vapply(seq_len(nrow(y_true)), function(k) {
    sum(dpois(y_true[k, ], lambda = pmax(lam_net[k, ], 1e-12), log = TRUE))
  }, numeric(1))
  
  ls_non_t <- vapply(seq_len(nrow(y_true)), function(k) {
    sum(dpois(y_true[k, ], lambda = pmax(lam_non[k, ], 1e-12), log = TRUE))
  }, numeric(1))
  
  dm_mae <- dm_test_nw(mae_net_t, mae_non_t, h = h)
  ci_mae <- block_boot_ci(mae_net_t - mae_non_t, B = 4000, block_len = max(4, h), seed = 1)
  
  dm_ls <- dm_test_nw(-ls_net_t, -ls_non_t, h = h)  # DM on negative log score (loss)
  ci_ls <- block_boot_ci((-ls_net_t) - (-ls_non_t), B = 4000, block_len = max(4, h), seed = 1)
  
  tibble::tibble(
    horizon = h,
    n_pairs = length(mae_net_t),
    W_scenario = labelW,
    
    MAE_net = mean(mae_net_t),
    MAE_nonet = mean(mae_non_t),
    dMAE = mean(mae_net_t - mae_non_t),
    dMAE_lo = ci_mae[1],
    dMAE_hi = ci_mae[2],
    DM_MAE_lag = dm_mae$lag,
    DM_MAE_p = dm_mae$pval,
    
    mean_logscore_net = mean(ls_net_t),
    mean_logscore_nonet = mean(ls_non_t),
    dNegLogScore = mean((-ls_net_t) - (-ls_non_t)),
    dNLS_lo = ci_ls[1],
    dNLS_hi = ci_ls[2],
    DM_NLS_lag = dm_ls$lag,
    DM_NLS_p = dm_ls$pval
  )
}

tab_crime_h <- dplyr::bind_rows(lapply(h_grid, eval_pois_h, W_use = W_full, labelW = "oracleW"))
print(tab_crime_h)

# ----------------------------
# B) Nodewise improvements
# ----------------------------
nodewise_crime <- function(h) {
  pr <- make_pairs(idx_target_test, h, t_origin_min = 2L, t_target_max = Tt)
  idx_origin <- pr$origin
  idx_target <- pr$target
  y_true <- Y_full[idx_target, , drop=FALSE]
  
  lam_net <- matrix(NA_real_, nrow=length(idx_origin), ncol=N)
  lam_non <- matrix(NA_real_, nrow=length(idx_origin), ncol=N)
  
  for (k in seq_along(idx_origin)) {
    t0 <- idx_origin[k]
    lam_net[k,] <- pois_fc_h_mean(t0, h, theta_at(t0, theta_net_filt), W_full)
    lam_non[k,] <- pois_fc_h_mean_nonet(t0, h, theta_at(t0, theta_non_filt))
  }
  
  mae_net_i <- colMeans(abs(y_true - lam_net))
  mae_non_i <- colMeans(abs(y_true - lam_non))
  
  tibble::tibble(node = seq_len(N), h = h, dMAE = mae_net_i - mae_non_i)
}

df_nodes <- dplyr::bind_rows(lapply(h_grid, nodewise_crime))
print(df_nodes %>% dplyr::group_by(h) %>%
        dplyr::summarize(frac_improve = mean(dMAE < 0),
                         median_dMAE = median(dMAE),
                         .groups="drop"))

ggplot(df_nodes, aes(x = factor(h), y = dMAE)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title="Crime: distribution of network benefit across regions",
       x="Horizon h", y="Regionwise dMAE (net - nonet)") +
  theme_minimal()

# ----------------------------
# C) Calibration: randomized PIT (1-step test window)
# ----------------------------
idx_test_curr <- (nrow(Y_curr) - test_h + 1):nrow(Y_curr)

y_true_1 <- as.vector(Y_curr[idx_test_curr, , drop=FALSE])
lam_1    <- as.vector(fit_pois_net$yhat[idx_test_curr, , drop=FALSE])

pit_vals <- pit_poisson_randomized(y_true_1, lam_1, seed = 1)

ggplot(tibble::tibble(pit = pit_vals), aes(x = pit)) +
  geom_histogram(bins = 30) +
  labs(title="Crime: randomized PIT (1-step, Poisson plug-in)",
       x="randomized PIT", y="count") +
  theme_minimal()

# ----------------------------
# D) Sensitivity curve (static W perturbation)
# ----------------------------
sens_curve_crime <- function(h = 4, alpha_grid = seq(0, 0.6, by = 0.1),
                             mode = c("mix","delete","rewire")) {
  mode <- match.arg(mode)
  
  pr <- make_pairs(idx_target_test, h, t_origin_min = 2L, t_target_max = Tt)
  idx_origin <- pr$origin
  idx_target <- pr$target
  y_true <- Y_full[idx_target, , drop=FALSE]
  
  purrr::map_dfr(alpha_grid, function(a) {
    W_hat <- if (mode == "mix") {
      as(Matrix::Matrix(W_mix(as.matrix(W_full), alpha = a, seed = 100), sparse = TRUE), "dgCMatrix")
    } else if (mode == "delete") {
      W_delete_edges(W_full, p_delete = a, seed = 100)
    } else {
      as(Matrix::Matrix(W_rewire(as.matrix(W_full), p_rewire = a, seed = 100), sparse = TRUE), "dgCMatrix")
    }
    
    dW <- op_norm2(W_hat - W_full)
    
    lam_hat <- matrix(NA_real_, nrow=length(idx_origin), ncol=N)
    for (k in seq_along(idx_origin)) {
      t0 <- idx_origin[k]
      lam_hat[k,] <- pois_fc_h_mean(t0, h, theta_at(t0, theta_net_filt), W_hat)
    }
    
    mae_t <- rowMeans(abs(y_true - lam_hat))
    
    tibble::tibble(alpha = a, mode = mode, h = h, dW_sup = dW, test_MAE = mean(mae_t))
  })
}

df_sens <- dplyr::bind_rows(
  sens_curve_crime(h = 4, mode = "mix"),
  sens_curve_crime(h = 4, mode = "delete"),
  sens_curve_crime(h = 4, mode = "rewire")
)

ggplot(df_sens, aes(x = dW_sup, y = test_MAE, shape = mode)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  labs(
    title = "Crime: sensitivity of multi-step MAE to network perturbations",
    x = "sup_t ||W_hat,t - W_t||_op",
    y = "Test MAE (h=4)"
  ) +
  theme_minimal()
