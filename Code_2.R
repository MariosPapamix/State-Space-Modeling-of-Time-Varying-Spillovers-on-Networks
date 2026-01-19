# =============================================================================
# GDP + trade NSSM experiment: adds baselines + spectral-radius stability
# Requires: nssm_cache/gdp_trade_gvar2019_1980_2016.rds from your working snippet.
# =============================================================================

pkgs <- c("zoo","RSpectra","ggplot2","tibble","dplyr")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install)) install.packages(to_install)
invisible(lapply(pkgs, library, character.only = TRUE))

# ---------- Load prepared data ----------
obj <- readRDS("nssm_cache/gdp_trade_gvar2019_1980_2016.rds")
Xz <- obj$X
W_t <- obj$W_t
countries <- obj$countries

Tt <- nrow(Xz)
N  <- ncol(Xz)
qtrs <- zoo::index(Xz)

cat("Loaded panel: T=", Tt, " quarters, N=", N, " countries\n")

# Train/test windows (per your design)
train_end <- zoo::as.yearqtr("2014 Q4")
test_start <- zoo::as.yearqtr("2015 Q1")
test_end   <- zoo::as.yearqtr("2016 Q4")

# For Y_curr we use qtrs_curr = qtrs[2:Tt]
qtrs_curr <- qtrs[2:Tt]
idx_train_curr <- which(qtrs_curr <= train_end)
idx_test <- which(qtrs_curr >= test_start & qtrs_curr <= test_end)
stopifnot(length(idx_train_curr) > 0, length(idx_test) > 0)

# ---------- Helpers ----------
op_norm2 <- function(M) {
  M <- as.matrix(M)
  if (nrow(M) <= 200) svd(M, nu = 0, nv = 0)$d[1] else RSpectra::svds(M, k = 1)$d[1]
}
stationary_dist <- function(W, tol=1e-10, maxit=5000) {
  W <- as.matrix(W)
  n <- nrow(W)
  pi <- rep(1/n, n)
  for (it in 1:maxit) {
    pi_new <- as.numeric(t(pi) %*% W)
    pi_new <- pi_new / sum(pi_new)
    if (max(abs(pi_new - pi)) < tol) break
    pi <- pi_new
  }
  pi
}
spectral_radius <- function(M) {
  ev <- eigen(M, only.values = TRUE)$values
  max(Mod(ev))
}

# Random-walk Kalman filter for multivariate dynamic regression:
# y_t (N-vector), X_t (N x K), theta_t = theta_{t-1} + u_t, u_t ~ N(0,Q_t)
# eps_t ~ N(0, sigma2 I_N)
kalman_rw_gaussian <- function(Y, X_list, sigma2, Q_list, m0, C0) {
  TT <- nrow(Y); N <- ncol(Y); K <- length(m0)
  stopifnot(length(X_list) == TT, length(Q_list) == TT)
  
  m_pred <- matrix(NA_real_, TT, K)
  m_filt <- matrix(NA_real_, TT, K)
  yhat   <- matrix(NA_real_, TT, N)
  
  m_prev <- m0
  C_prev <- C0
  
  for (t in 1:TT) {
    # predict
    m_t_pred <- m_prev
    C_t_pred <- C_prev + Q_list[[t]]
    
    X_t <- X_list[[t]]
    y_t <- as.numeric(Y[t, ])
    
    # one-step mean E[y_t | y_{1:t-1}]
    yhat[t, ] <- as.numeric(X_t %*% m_t_pred)
    
    # conjugate update
    XtX <- crossprod(X_t)
    C_post <- solve(solve(C_t_pred) + XtX / sigma2)
    
    resid <- y_t - as.numeric(X_t %*% m_t_pred)
    s <- crossprod(X_t, resid)
    m_t_filt <- m_t_pred + as.numeric(C_post %*% (s / sigma2))
    
    m_pred[t, ] <- m_t_pred
    m_filt[t, ] <- m_t_filt
    
    m_prev <- m_t_filt
    C_prev <- C_post
  }
  
  list(m_pred=m_pred, m_filt=m_filt, yhat=yhat)
}

# ---------- Build Y, lags, and design matrices ----------
Y_mat  <- coredata(Xz)                 # T x N
Y_lag  <- Y_mat[1:(Tt-1), , drop=FALSE]
Y_curr <- Y_mat[2:Tt, , drop=FALSE]

W_lag <- lapply(2:Tt, function(t) W_t[t, , ])  # list length T-1

# Full network design: [1, (W_t Y_{t-1}), Y_{t-1}]
X_list <- vector("list", Tt-1)
for (t in 1:(Tt-1)) {
  x_net <- as.numeric(W_lag[[t]] %*% Y_lag[t, ])
  X_list[[t]] <- cbind(1, x_net, Y_lag[t, ])
}

# No-network ablation design: [1, Y_{t-1}]  (beta1,t ≡ 0)
X_list_nonet <- lapply(1:(Tt-1), function(t) cbind(1, Y_lag[t, ]))

# ---------- Fit models ----------
sigma2 <- 0.02^2
Q_full  <- replicate(Tt-1, diag(c(1e-4, 1e-4, 1e-4)), simplify = FALSE)
Q_nonet <- replicate(Tt-1, diag(c(1e-4, 1e-4)), simplify = FALSE)

fit_full <- kalman_rw_gaussian(Y_curr, X_list, sigma2, Q_full,  m0=c(0,0,0), C0=diag(1,3))
fit_nonet <- kalman_rw_gaussian(Y_curr, X_list_nonet, sigma2, Q_nonet, m0=c(0,0), C0=diag(1,2))

# ---------- Static network autoregression baseline (training-only OLS) ----------
X_stack_train <- do.call(rbind, X_list[idx_train_curr])
y_stack_train <- as.vector(t(Y_curr[idx_train_curr, , drop=FALSE]))
beta_static_net <- qr.solve(X_stack_train, y_stack_train)

X_stack_all <- do.call(rbind, X_list)
yhat_static_net_stack <- as.numeric(X_stack_all %*% beta_static_net)
yhat_static_net <- matrix(yhat_static_net_stack, nrow=(Tt-1), ncol=N, byrow=TRUE)

# ---------- VAR(1) benchmark baseline (training-only OLS) ----------
# VAR(1): Y_t = c + A Y_{t-1} + eps.  Estimate (c, A) on training only.
X_var_train <- cbind(1, Y_lag[idx_train_curr, , drop=FALSE])         # (n_train x (1+N))
Y_var_train <- Y_curr[idx_train_curr, , drop=FALSE]                  # (n_train x N)
B_var <- qr.solve(X_var_train, Y_var_train)                          # ((1+N) x N)

X_var_all <- cbind(1, Y_lag)                                         # ((T-1) x (1+N))
yhat_var <- X_var_all %*% B_var                                      # ((T-1) x N)

# ---------- Full model with carry-forward network forecast (test network unknown) ----------
# carry-forward: use a 2014 network for all test periods
i2014 <- which(as.integer(floor(as.numeric(qtrs))) == 2014)
stopifnot(length(i2014) > 0)
W_cf <- W_t[i2014[1], , ]

W_lag_cf <- W_lag
for (k in seq_along(W_lag_cf)) {
  qt <- qtrs_curr[k]
  if (qt > train_end) W_lag_cf[[k]] <- W_cf
}
X_list_cf <- vector("list", Tt-1)
for (t in 1:(Tt-1)) {
  x_net <- as.numeric(W_lag_cf[[t]] %*% Y_lag[t, ])
  X_list_cf[[t]] <- cbind(1, x_net, Y_lag[t, ])
}
fit_full_cf <- kalman_rw_gaussian(Y_curr, X_list_cf, sigma2, Q_full, m0=c(0,0,0), C0=diag(1,3))

dW_sup <- max(sapply(idx_test, function(ii) op_norm2(W_lag_cf[[ii]] - W_lag[[ii]])))

# ---------- Metrics ----------
mse <- function(A, B, idx) mean((A[idx, , drop=FALSE] - B[idx, , drop=FALSE])^2)

mse_full_oracle <- mse(Y_curr, fit_full$yhat, idx_test)
mse_full_cf     <- mse(Y_curr, fit_full_cf$yhat, idx_test)
mse_nonet_tvp    <- mse(Y_curr, fit_nonet$yhat, idx_test)
mse_static_net   <- mse(Y_curr, yhat_static_net, idx_test)
mse_var          <- mse(Y_curr, yhat_var, idx_test)

# Aggregation recursion (Theorem agg-scalar)
train_quarters_idx <- which(as.integer(floor(as.numeric(qtrs))) <= 2014)
W_bar_train <- Reduce("+", lapply(train_quarters_idx, function(t) W_t[t, , ])) / length(train_quarters_idx)
pi_hat <- stationary_dist(W_bar_train)

agg_y <- as.numeric(Y_curr %*% pi_hat)
agg_y_lag <- c(sum(pi_hat * Y_mat[1, ]), agg_y[1:(length(agg_y)-1)])
beta_sum <- fit_full$m_pred[,2] + fit_full$m_pred[,3]
agg_pred <- fit_full$m_pred[,1] + beta_sum * agg_y_lag

agg_mae_all  <- mean(abs(agg_y - agg_pred))
agg_mae_test <- mean(abs(agg_y[idx_test] - agg_pred[idx_test]))

# Stability diagnostic: rho(beta1_t W_t + beta2_t I)
rho_hat <- sapply(1:(Tt-1), function(t) {
  b1 <- fit_full$m_filt[t,2]
  b2 <- fit_full$m_filt[t,3]
  spectral_radius(b1 * W_lag[[t]] + b2 * diag(N))
})

# ---------- Summary table ----------
summary_tbl <- tibble::tibble(
  model = c("NTVP-Network (oracle W)",
            "NTVP-Network (carry-forward W)",
            "NTVP No-network (beta1,t = 0)",
            "Static network autoreg (const beta)",
            "VAR(1) OLS"),
  test_MSE = c(mse_full_oracle, mse_full_cf, mse_nonet_tvp, mse_static_net, mse_var)
) |>
  dplyr::mutate(delta_vs_oracle = test_MSE - mse_full_oracle)

print(summary_tbl)

cat(sprintf("sup_t ||W_hat - W||_op (test): %.6f\n", dW_sup))
cat(sprintf("Aggregation MAE (all):  %.6f\n", agg_mae_all))
cat(sprintf("Aggregation MAE (test): %.6f\n", agg_mae_test))
cat(sprintf("Stability: max rho_t over time: %.3f\n", max(rho_hat)))
cat(sprintf("Stability: fraction rho_t < 1: %.3f\n", mean(rho_hat < 1)))

# ---------- Plots ----------
df_plot <- tibble::tibble(
  qtr = qtrs_curr,
  beta1 = fit_full$m_filt[,2],
  beta2 = fit_full$m_filt[,3],
  beta_sum = fit_full$m_filt[,2] + fit_full$m_filt[,3],
  rho = rho_hat
)

p_beta1 <- ggplot(df_plot, aes(x=qtr, y=beta1)) +
  geom_line(linewidth=0.6) +
  labs(x="Quarter", y=expression(hat(beta)[1*t]~"(network spillover)")) +
  theme_minimal()

p_stab_proxy <- ggplot(df_plot, aes(x=qtr, y=beta_sum)) +
  geom_line(linewidth=0.6) +
  geom_hline(yintercept=1, colour="red") +
  labs(x="Quarter", y=expression(hat(beta)[1*t]+hat(beta)[2*t]~"(proxy)")) +
  theme_minimal()

p_rho <- ggplot(df_plot, aes(x=qtr, y=rho)) +
  geom_line(linewidth=0.6) +
  geom_hline(yintercept=1, colour="red") +
  labs(x="Quarter", y=expression(rho(hat(beta)[1*t]*W[t]+hat(beta)[2*t]*I))) +
  theme_minimal()

df_agg <- tibble::tibble(qtr=qtrs_curr, agg_actual=agg_y, agg_pred=agg_pred)
p_agg <- ggplot(df_agg, aes(x=qtr)) +
  geom_line(aes(y=agg_actual), linewidth=0.6) +
  geom_line(aes(y=agg_pred), linewidth=0.6, linetype="dashed") +
  labs(x="Quarter", y=expression(pi^top*Y[t]~": actual vs recursion")) +
  theme_minimal()

print(p_beta1)
print(p_stab_proxy)
print(p_rho)
print(p_agg)

# ---------- OPTIONAL: save figures for LaTeX ----------
# dir.create("figs", showWarnings = FALSE)
ggsave("figs/gdp_beta1.pdf", p_beta1, width=7, height=4)
ggsave("figs/gdp_stability_proxy.pdf", p_stab_proxy, width=7, height=4)
ggsave("figs/gdp_rho_spectralradius.pdf", p_rho, width=7, height=4)
ggsave("figs/gdp_aggregation_fit.pdf", p_agg, width=7, height=4)

# ---------- OPTIONAL: export summary table ----------
write.csv(summary_tbl, "figs/gdp_model_comparison.csv", row.names = FALSE)






# =============================================================================
# GDP MULTI-STEP EVALUATION (Gaussian)
#   - GDP/crime mixing: W_t (148x33x33) must pair with Y_gdp (148x33)
#   - Builds W_full from W_t dims (not from Y_mat)
#   - Auto-detects the GDP Y matrix in your workspace by matching dims to W_t
#   - DM test + block bootstrap are made "small-sample safe" (no hard stop)
# =============================================================================

suppressPackageStartupMessages({
  library(tibble); library(dplyr); library(ggplot2)
})

# ----------------------------
# 1) Safe DM + safe block bootstrap (no hard stops on short samples)
# ----------------------------
dm_test_nw_safe <- function(lossA, lossB, h = 1, lag = NULL, min_T = 5) {
  d <- as.numeric(lossA - lossB)
  d <- d[is.finite(d)]
  Tn <- length(d)
  
  if (Tn < min_T) {
    return(list(T = Tn, mean_diff = mean(d), stat = NA_real_, pval = NA_real_, lag = NA_integer_))
  }
  
  if (is.null(lag)) lag <- max(h - 1, 0)
  lag <- min(lag, Tn - 1)
  
  dbar <- mean(d)
  x <- d - dbar
  
  gamma0 <- sum(x^2) / Tn
  vhat <- gamma0
  if (lag > 0) {
    for (k in 1:lag) {
      gamma_k <- sum(x[(k + 1):Tn] * x[1:(Tn - k)]) / Tn
      w_k <- 1 - k / (lag + 1)
      vhat <- vhat + 2 * w_k * gamma_k
    }
  }
  
  se <- sqrt(vhat / Tn)
  if (!is.finite(se) || se <= 0) {
    return(list(T = Tn, mean_diff = dbar, stat = NA_real_, pval = NA_real_, lag = lag))
  }
  
  stat <- dbar / se
  pval <- 2 * pt(-abs(stat), df = Tn - 1)
  list(T = Tn, mean_diff = dbar, stat = stat, pval = pval, lag = lag)
}

block_boot_ci_core <- function(d, B = 2000, block_len = 4,
                               probs = c(0.025, 0.975), seed = 1) {
  set.seed(seed)
  d <- as.numeric(d); d <- d[is.finite(d)]
  n <- length(d)
  
  nb <- ceiling(n / block_len)
  starts <- 1:(n - block_len + 1)
  
  boot_means <- replicate(B, {
    idx <- unlist(lapply(sample(starts, nb, replace = TRUE),
                         function(s) s:(s + block_len - 1)))
    idx <- idx[1:n]
    mean(d[idx])
  })
  
  quantile(boot_means, probs = probs, names = FALSE)
}

block_boot_ci_safe <- function(d, B = 2000, block_len = 4,
                               probs = c(0.025, 0.975), seed = 1) {
  d <- as.numeric(d); d <- d[is.finite(d)]
  n <- length(d)
  
  # too short -> return NA CI (do NOT stop)
  if (n < 6) return(c(NA_real_, NA_real_))
  
  # ensure n >= block_len + 5 by shrinking block_len if needed
  bl <- min(block_len, n - 5L)
  bl <- max(1L, bl)
  if (n < bl + 5L) return(c(NA_real_, NA_real_))
  
  block_boot_ci_core(d, B = B, block_len = bl, probs = probs, seed = seed)
}

# ----------------------------
# 2) Coerce W_t into a list of matrices W_full[[t]] (t=1..T)
# ----------------------------
coerce_W_list <- function(W_t_obj) {
  if (is.list(W_t_obj)) {
    W_list <- W_t_obj
    Tw <- length(W_list)
    Nw <- nrow(W_list[[1]])
    stopifnot(all(vapply(W_list, function(W) all(dim(W) == c(Nw, Nw)), logical(1))))
    return(list(W = W_list, Tw = Tw, Nw = Nw))
  }
  
  if (is.array(W_t_obj) && length(dim(W_t_obj)) == 3) {
    dd <- dim(W_t_obj)
    
    # case A: T x N x N
    if (dd[2] == dd[3] && dd[1] != dd[2]) {
      Tw <- dd[1]; Nw <- dd[2]
      W_list <- lapply(seq_len(Tw), function(t) W_t_obj[t, , ])
      return(list(W = W_list, Tw = Tw, Nw = Nw))
    }
    
    # case B: N x N x T
    if (dd[1] == dd[2] && dd[3] != dd[1]) {
      Nw <- dd[1]; Tw <- dd[3]
      W_list <- lapply(seq_len(Tw), function(t) W_t_obj[, , t])
      return(list(W = W_list, Tw = Tw, Nw = Nw))
    }
    
    stop("W_t has ambiguous/unexpected dims: ", paste(dd, collapse = " x "),
         ". Expected (T x N x N) with T!=N, or (N x N x T) with T!=N.")
  }
  
  stop("W_t must be a list of matrices or a 3D array.")
}

# ----------------------------
# 3) Find the GDP Y matrix by matching dims to W_t (Tw x Nw)
#    This prevents you from accidentally using the crime matrix (72x552).
# ----------------------------
find_Y_matching_W <- function(Tw, Nw, env = .GlobalEnv) {
  obj_names <- ls(envir = env)
  
  score_name <- function(nm) {
    nm2 <- tolower(nm)
    score <- 0
    if (grepl("gdp", nm2))   score <- score + 5
    if (grepl("trade", nm2)) score <- score + 2
    if (grepl("y", nm2))     score <- score + 1
    score
  }
  
  cands <- list()
  for (nm in obj_names) {
    obj <- get(nm, envir = env)
    if (is.matrix(obj) || is.data.frame(obj)) {
      dm <- dim(obj)
      if (!is.null(dm) && length(dm) == 2) {
        if (all(dm == c(Tw, Nw))) {
          cands[[nm]] <- list(Y = as.matrix(obj), score = 10 + score_name(nm))
        } else if (all(dm == c(Nw, Tw))) {
          cands[[nm]] <- list(Y = t(as.matrix(obj)), score = 10 + score_name(nm))
        }
      }
    }
  }
  
  if (length(cands) == 0) {
    # show a compact inventory to help you spot the correct object name
    mats <- lapply(obj_names, function(nm) {
      obj <- get(nm, envir = env)
      if (is.matrix(obj) || is.data.frame(obj)) {
        dm <- dim(obj)
        if (!is.null(dm) && length(dm) == 2) return(c(name = nm, nrow = dm[1], ncol = dm[2]))
      }
      NULL
    })
    mats <- do.call(rbind, mats[!vapply(mats, is.null, logical(1))])
    
    stop(
      "No GDP Y matrix found that matches W_t dims.\n",
      "W_t implies you need a matrix with dim ", Tw, " x ", Nw, " (or ", Nw, " x ", Tw, ").\n",
      "But your current Y_mat is likely the crime panel (72 x 552).\n\n",
      "load/restore your GDP panel into an object like Y_gdp (", Tw, "x", Nw, "),\n",
      "then rerun. Available matrix/data.frame objects and dims:\n",
      paste(apply(mats, 1, function(r) paste0("  - ", r[1], ": ", r[2], " x ", r[3])), collapse="\n")
    )
  }
  
  pick <- names(which.max(vapply(cands, function(z) z$score, numeric(1))))
  message("Using GDP Y matrix from object: ", pick,
          " (dim ", nrow(cands[[pick]]$Y), " x ", ncol(cands[[pick]]$Y), ")")
  cands[[pick]]$Y
}

# ----------------------------
# 4) Align W length to Y length (allow +/-1 only)
# ----------------------------
align_W_to_Y <- function(W_list, Ty, prefer = c("drop_first", "drop_last")) {
  prefer <- match.arg(prefer)
  Tw <- length(W_list)
  
  if (Tw == Ty) return(W_list)
  
  if (Tw == Ty + 1L) {
    message("Aligning W to Y: dropping one W (Tw=", Tw, " -> Ty=", Ty, ").")
    if (prefer == "drop_first") return(W_list[-1]) else return(W_list[1:Ty])
  }
  
  if (Tw == Ty - 1L) {
    message("Aligning W to Y: padding first W (Tw=", Tw, " -> Ty=", Ty, ").")
    return(c(list(W_list[[1]]), W_list))
  }
  
  stop("Can't align W to Y: length(W)=", Tw, " but nrow(Y)=", Ty,
       ". This usually means you're mixing GDP and crime objects.")
}

# ----------------------------
# 5) BUILD GDP OBJECTS (this is the part you were failing)
# ----------------------------
W_info <- coerce_W_list(W_t)
W_list_raw <- W_info$W
Tw <- W_info$Tw
Nw <- W_info$Nw

# Find the *GDP* Y matrix (148x33), NOT the crime matrix (72x552)
Y_full <- find_Y_matching_W(Tw = Tw, Nw = Nw, env = .GlobalEnv)

# Align W to Y (should be equal for GDP; if your Y is 147x33 this will drop 1 W)
W_full <- align_W_to_Y(W_list_raw, Ty = nrow(Y_full), prefer = "drop_first")

Tt <- nrow(Y_full)
N  <- ncol(Y_full)

stopifnot(length(W_full) %in% c(Tt, Tt - 1))
stopifnot(all(vapply(W_full, function(W) all(dim(W) == c(N, N)), logical(1))))

cat(sprintf("GDP objects ready: Y_full dim = %d x %d; W_full length = %d; each W is %d x %d\n",
            Tt, N, length(W_full), N, N))

# ----------------------------
# 6) Helpers: theta_at and W_at (support T or T-1 storage)
# ----------------------------
theta_full_filt  <- as.matrix(fit_full$m_filt)
theta_nonet_filt <- as.matrix(fit_nonet$m_filt)

theta_at <- function(t, theta_mat, Ty) {
  nr <- nrow(theta_mat)
  if (nr == Ty)     return(theta_mat[t, , drop = TRUE])
  if (nr == Ty - 1) return(theta_mat[t - 1, , drop = TRUE])  # row 1 corresponds to time 2
  stop("theta_mat nrow must be Ty or Ty-1. Got ", nr, " with Ty=", Ty)
}

W_at <- function(t, W_list, Ty) {
  Tw <- length(W_list)
  if (Tw == Ty) idx <- t
  else if (Tw == Ty - 1L) idx <- t - 1L
  else stop("W_list length must be Ty or Ty-1 (got ", Tw, ").")
  if (idx < 1L || idx > Tw) stop("W_at: time t=", t, " not available.")
  W_list[[idx]]
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
# 7) Define test targets for GDP
#     If you already have idx_target_test for GDP, keep it.
#     Otherwise default to the last 8 quarters (2015Q1–2016Q4 is 8 quarters in the paper). 
# ----------------------------
if (!exists("idx_target_test", inherits = TRUE) || length(idx_target_test) == 0) {
  idx_target_test <- (Tt - 7L):Tt
  message("idx_target_test not found; using last 8 times: ", paste(idx_target_test, collapse = ", "))
} else {
  # sanity check (common failure: idx_target_test built from crime dates)
  if (max(idx_target_test) > Tt) stop("idx_target_test exceeds GDP length Tt=", Tt,
                                      ". You likely used the crime date vector.")
}

# ----------------------------
# 8) Forecast recursions
# ----------------------------
gauss_fc_h_full <- function(t_origin, h, th, W_use, Ty) {
  b0 <- th[1]; b1 <- th[2]; b2 <- th[3]
  y_prev <- as.numeric(Y_full[t_origin, ])
  for (j in seq_len(h)) {
    Wtj <- W_at(t_origin + j, W_use, Ty)
    y_prev <- b0 + b1 * as.numeric(Wtj %*% y_prev) + b2 * y_prev
  }
  y_prev
}

gauss_fc_h_nonet <- function(t_origin, h, th) {
  b0 <- th[1]; b2 <- th[2]
  y_prev <- as.numeric(Y_full[t_origin, ])
  for (j in seq_len(h)) y_prev <- b0 + b2 * y_prev
  y_prev
}

gauss_fc_h_static_net <- function(t_origin, h, beta_static, W_use, Ty) {
  y_prev <- as.numeric(Y_full[t_origin, ])
  for (j in seq_len(h)) {
    Wtj <- W_at(t_origin + j, W_use, Ty)
    y_prev <- beta_static[1] +
      beta_static[2] * as.numeric(Wtj %*% y_prev) +
      beta_static[3] * y_prev
  }
  y_prev
}

gauss_fc_h_var <- function(t_origin, h, B_var) {
  cvec <- as.numeric(B_var[1, ])
  A <- t(B_var[-1, , drop = FALSE])  # N x N
  y_prev <- as.numeric(Y_full[t_origin, ])
  for (j in seq_len(h)) y_prev <- cvec + as.numeric(A %*% y_prev)
  y_prev
}

# ----------------------------
# 9) Horizon evaluation (oracle W vs carry-forward W if you have W_cf)
# ----------------------------
eval_gdp_h <- function(h, W_use_full, labelW = "oracleW") {
  pr <- make_pairs(idx_target_test, h, t_origin_min = 2L, t_target_max = Tt)
  idx_origin <- pr$origin
  idx_target <- pr$target
  
  y_true <- Y_full[idx_target, , drop = FALSE]
  
  yhat_full <- matrix(NA_real_, nrow = length(idx_origin), ncol = N)
  yhat_non  <- matrix(NA_real_, nrow = length(idx_origin), ncol = N)
  
  have_static <- exists("beta_static_net", inherits = TRUE) && length(beta_static_net) == 3
  have_var    <- exists("B_var", inherits = TRUE)
  
  yhat_stat <- if (have_static) matrix(NA_real_, nrow = length(idx_origin), ncol = N) else NULL
  yhat_varh <- if (have_var)    matrix(NA_real_, nrow = length(idx_origin), ncol = N) else NULL
  
  for (k in seq_along(idx_origin)) {
    t0 <- idx_origin[k]
    th_full <- theta_at(t0, theta_full_filt,  Tt)
    th_non  <- theta_at(t0, theta_nonet_filt, Tt)
    
    yhat_full[k, ] <- gauss_fc_h_full(t0, h, th_full, W_use_full, Tt)
    yhat_non[k, ]  <- gauss_fc_h_nonet(t0, h, th_non)
    
    if (have_static) yhat_stat[k, ] <- gauss_fc_h_static_net(t0, h, beta_static_net, W_use_full, Tt)
    if (have_var)    yhat_varh[k, ] <- gauss_fc_h_var(t0, h, B_var)
  }
  
  loss_full <- rowMeans((y_true - yhat_full)^2)
  loss_non  <- rowMeans((y_true - yhat_non)^2)
  
  dm_full_vs_non <- dm_test_nw_safe(loss_full, loss_non, h = h, min_T = 5)
  ci_full_vs_non <- block_boot_ci_safe(loss_full - loss_non, B = 4000, block_len = max(4, h), seed = 1)
  
  out <- tibble::tibble(
    horizon = h,
    W_scenario = labelW,
    n_eval = length(loss_full),
    MSE_full = mean(loss_full),
    MSE_nonet = mean(loss_non),
    diff_full_minus_nonet = mean(loss_full - loss_non),
    diff_CI_lo = ci_full_vs_non[1],
    diff_CI_hi = ci_full_vs_non[2],
    DM_stat = dm_full_vs_non$stat,
    DM_p = dm_full_vs_non$pval
  )
  
  if (have_static) out$MSE_static_net <- mean(rowMeans((y_true - yhat_stat)^2))
  else            out$MSE_static_net <- NA_real_
  
  if (have_var)    out$MSE_VAR <- mean(rowMeans((y_true - yhat_varh)^2))
  else             out$MSE_VAR <- NA_real_
  
  out
}

h_grid <- c(1, 2, 4, 8)

# --- oracle W
tab_gdp_oracle <- dplyr::bind_rows(lapply(h_grid, eval_gdp_h,
                                          W_use_full = W_full,
                                          labelW = "oracleW"))
print(tab_gdp_oracle)

# --- carry-forward W (optional)
if (exists("W_cf", inherits = TRUE)) {
  W_cf_list <- W_full
  train_end_idx <- min(idx_target_test) - 1L   # assumes test starts right after train
  for (t in seq_len(Tt)) {
    if (t > train_end_idx) W_cf_list[[t]] <- W_cf
  }
  
  tab_gdp_cf <- dplyr::bind_rows(lapply(h_grid, eval_gdp_h,
                                        W_use_full = W_cf_list,
                                        labelW = "carryforwardW"))
  print(tab_gdp_cf)
  
  tab_gdp_all <- dplyr::bind_rows(tab_gdp_oracle, tab_gdp_cf)
  ggplot(tab_gdp_all, aes(x = horizon, y = MSE_full, linetype = W_scenario)) +
    geom_line() + geom_point() +
    labs(x = "Horizon h", y = "Test MSE",
         title = "GDP: multi-step MSE for full network TVP") +
    theme_minimal()
} else {
  message("W_cf not found; skipping carry-forward scenario.")
}
