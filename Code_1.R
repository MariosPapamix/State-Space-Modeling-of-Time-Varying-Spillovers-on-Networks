# ======================================================================
# EXPERIMENT 1 (Simulation): stability, aggregation, sparse increments,
#                            and network misspecification sensitivity
# Produces:
#   - One-step MSE comparisons (full / no-network / static)
#   - Threshold rule TPR/FPR for sparse increments in beta1,t
#   - Sensitivity curve: mean||mu_oracle - mu_plugin||^2 vs sup_t||W_hat-W||_op
#   - Aggregation recursion MAE
#   - Two plots (stability diagnostics + sensitivity curve)
# ======================================================================

# ---------- Packages ----------
pkgs <- c("Matrix","RSpectra","ggplot2","dplyr","tibble","purrr")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if(length(to_install)) install.packages(to_install)
invisible(lapply(pkgs, library, character.only = TRUE))

set.seed(1)

# ---------- Helpers ----------
row_stochastic <- function(A, diag_zero = TRUE, eps = 1e-12) {
  A <- as.matrix(A)
  if (diag_zero) diag(A) <- 0
  rs <- rowSums(A)
  zero_rows <- which(rs < eps)
  if (length(zero_rows)) {
    A[zero_rows, ] <- 1
    if (diag_zero) diag(A)[zero_rows] <- 0
    rs <- rowSums(A)
  }
  A / rs
}

spectral_radius <- function(M) {
  ev <- eigen(M, only.values = TRUE)$values
  max(Mod(ev))
}

op_norm2 <- function(M) {
  M <- as.matrix(M)
  if (nrow(M) <= 200) {
    svd(M, nu = 0, nv = 0)$d[1]
  } else {
    RSpectra::svds(M, k = 1)$d[1]
  }
}

stationary_dist <- function(W, tol = 1e-12, maxit = 20000) {
  W <- as.matrix(W)
  N <- nrow(W)
  pi <- rep(1 / N, N)
  for (it in 1:maxit) {
    pi_new <- as.numeric(t(pi) %*% W)
    pi_new <- pi_new / sum(pi_new)
    if (max(abs(pi_new - pi)) < tol) break
    pi <- pi_new
  }
  pi
}

# ---------- Numerically stable Kalman filter for dynamic regression ----------
# Model (at each t):
#   y_t (N x 1) = X_t (N x K) m_t + eps_t,   eps_t ~ N(0, sigma2 I_N)
#   m_t = m_{t-1} + u_t,                    u_t ~ N(0, Q_t)
#
# Inputs:
#   Y:      T x N
#   X_list: list of length T, each N x K
#   Q_list: list of length T, each K x K
#
kalman_rw_gaussian <- function(Y, X_list, sigma2, Q_list, m0, C0, jitter = 1e-10) {
  Tt <- nrow(Y); N <- ncol(Y); K <- length(m0)
  stopifnot(length(X_list) == Tt, length(Q_list) == Tt)
  
  m_pred <- matrix(NA_real_, Tt, K)
  m_filt <- matrix(NA_real_, Tt, K)
  C_filt <- array(NA_real_, dim = c(K, K, Tt))
  yhat   <- matrix(NA_real_, Tt, N)
  loglik <- numeric(Tt)
  
  m_prev <- as.numeric(m0)
  C_prev <- as.matrix(C0)
  
  for (t in 1:Tt) {
    # prediction
    m_t_pred <- m_prev
    C_t_pred <- C_prev + Q_list[[t]]
    C_t_pred <- C_t_pred + diag(jitter, K)  # guard against singularity
    
    X_t <- as.matrix(X_list[[t]])  # N x K
    y_t <- as.numeric(Y[t, ])      # N
    
    # predictive mean
    yhat[t, ] <- as.numeric(X_t %*% m_t_pred)
    resid <- y_t - yhat[t, ]
    
    # predictive covariance for y: S = sigma2 I + X C X'
    S <- sigma2 * diag(N) + X_t %*% C_t_pred %*% t(X_t)
    
    # Cholesky for stability
    cholS <- chol(S)
    
    # Kalman gain: K = C X' S^{-1}
    K_gain <- C_t_pred %*% t(X_t)  # K x N
    # Multiply by S^{-1} without explicitly inverting S
    # K_gain <- K_gain %*% solve(S)
    K_gain <- t(backsolve(cholS, forwardsolve(t(cholS), t(K_gain))))  # K x N
    
    # update
    m_t_filt <- as.numeric(m_t_pred + K_gain %*% resid)
    C_post <- C_t_pred - K_gain %*% X_t %*% C_t_pred
    C_post <- 0.5 * (C_post + t(C_post))  # symmetrize
    
    # log predictive density (optional)
    Sinv_resid <- backsolve(cholS, forwardsolve(t(cholS), resid))
    quad <- sum(resid * Sinv_resid)
    logdetS <- 2 * sum(log(diag(cholS)))
    loglik[t] <- -0.5 * (N * log(2 * pi) + logdetS + quad)
    
    # store
    m_pred[t, ] <- m_t_pred
    m_filt[t, ] <- m_t_filt
    C_filt[,,t] <- C_post
    
    # roll
    m_prev <- m_t_filt
    C_prev <- C_post
  }
  
  list(m_pred = m_pred, m_filt = m_filt, C_filt = C_filt, yhat = yhat, loglik = loglik)
}

# ======================================================================
# 1) Simulate dynamic logistic edges -> W_t
# ======================================================================
N  <- 50
TT <- 200

# latent positions -> structured connectivity
d <- 2
z <- matrix(rnorm(N * d), N, d)
dist_mat <- as.matrix(dist(z))

# random walk controlling density
eta <- cumsum(c(0, rnorm(TT - 1, sd = 0.05))) + 1.0

W_list <- vector("list", TT)
for (t in 1:TT) {
  P <- plogis(eta[t] - dist_mat)
  diag(P) <- 0
  A <- matrix(rbinom(N * N, size = 1, prob = as.vector(P)), N, N)
  diag(A) <- 0
  W_list[[t]] <- row_stochastic(A)
}

# ======================================================================
# 2) Simulate piecewise-constant beta1,t (sparse increments) + beta0 RW
# ======================================================================
beta0 <- cumsum(rnorm(TT, sd = 0.02))
beta2 <- rep(0.35, TT)

beta1 <- rep(0.25, TT)
jump_times <- c(80, 120, 160)
beta1[jump_times[1]:jump_times[2]] <- 0.60  # near-instability window
beta1[jump_times[2]:jump_times[3]] <- 0.40
beta1[jump_times[3]:TT]            <- 0.20

rho_true <- sapply(1:TT, function(t) {
  spectral_radius(beta1[t] * W_list[[t]] + beta2[t] * diag(N))
})

# ======================================================================
# 3) Simulate Gaussian NTVP–VAR(1)
# ======================================================================
sigma <- 0.5
Y <- matrix(0, TT, N)
Y[1, ] <- rnorm(N, sd = 1)

for (t in 2:TT) {
  mu <- beta0[t] + beta1[t] * (W_list[[t]] %*% Y[t - 1, ]) + beta2[t] * Y[t - 1, ]
  Y[t, ] <- as.numeric(mu) + rnorm(N, sd = sigma)
}

# build design for t = 2..TT
Y_lag  <- Y[1:(TT - 1), , drop = FALSE]
Y_curr <- Y[2:TT,       , drop = FALSE]
W_use  <- W_list[2:TT]  # W_2,...,W_T

X_list <- vector("list", TT - 1)
for (t in 1:(TT - 1)) {
  x_net <- as.numeric(W_use[[t]] %*% Y_lag[t, ])
  X_list[[t]] <- cbind(1, x_net, Y_lag[t, ])
}

# ======================================================================
# 4) Fit: (a) full TVP network, (b) TVP no-network, (c) static OLS
# ======================================================================
sigma2 <- sigma^2

Q_full <- replicate(TT - 1, diag(c(1e-3, 5e-3, 1e-3)), simplify = FALSE)
fit_full <- kalman_rw_gaussian(
  Y = Y_curr, X_list = X_list, sigma2 = sigma2, Q_list = Q_full,
  m0 = c(0, 0, 0), C0 = diag(1, 3)
)

X_list_nonet <- lapply(X_list, function(X) X[, c(1, 3), drop = FALSE])
Q_nonet <- replicate(TT - 1, diag(c(1e-3, 1e-3)), simplify = FALSE)
fit_nonet <- kalman_rw_gaussian(
  Y = Y_curr, X_list = X_list_nonet, sigma2 = sigma2, Q_list = Q_nonet,
  m0 = c(0, 0), C0 = diag(1, 2)
)

# static OLS on stacked data (matches your construction)
X_stack <- do.call(rbind, X_list)         # (TT-1)*N x 3
y_stack <- as.vector(t(Y_curr))           # stacked by time-blocks
beta_static <- solve(crossprod(X_stack), crossprod(X_stack, y_stack))
yhat_static_stack <- as.numeric(X_stack %*% beta_static)
yhat_static <- matrix(yhat_static_stack, nrow = (TT - 1), ncol = N, byrow = TRUE)

mse_full   <- mean((Y_curr - fit_full$yhat)^2)
mse_nonet  <- mean((Y_curr - fit_nonet$yhat)^2)
mse_static <- mean((Y_curr - yhat_static)^2)

cat(sprintf("One-step MSE: full=%.4f | no-network=%.4f | static=%.4f\n",
            mse_full, mse_nonet, mse_static))

# ======================================================================
# 5) Sparse increments: simple threshold rule on |Δ beta1_hat|
# ======================================================================
beta1_hat <- fit_full$m_filt[, 2]
d_thr <- as.numeric(quantile(abs(diff(beta1_hat)), probs = 0.95))
s_hat <- c(FALSE, abs(diff(beta1_hat)) > d_thr)

beta1_true_sub <- beta1[2:TT]
s_true <- c(FALSE, diff(beta1_true_sub) != 0)

TPR <- mean(s_hat[s_true], na.rm = TRUE)
FPR <- mean(s_hat[!s_true], na.rm = TRUE)
cat(sprintf("Threshold rule: TPR=%.3f, FPR=%.3f (d=%.4f)\n", TPR, FPR, d_thr))

# ======================================================================
# 6) Network misspecification sensitivity curve (plug-in vs oracle means)
# ======================================================================
alpha_grid <- seq(0, 0.5, by = 0.05)

sens_df <- purrr::map_dfr(alpha_grid, function(a) {
  W_hat_list <- lapply(W_use, function(Wt) {
    Wnoise <- row_stochastic(matrix(runif(N * N), N, N))
    row_stochastic((1 - a) * Wt + a * Wnoise)
  })
  
  yhat_oracle <- fit_full$yhat
  yhat_plugin <- matrix(NA_real_, TT - 1, N)
  
  for (t in 1:(TT - 1)) {
    Xhat <- cbind(1, as.numeric(W_hat_list[[t]] %*% Y_lag[t, ]), Y_lag[t, ])
    yhat_plugin[t, ] <- as.numeric(Xhat %*% fit_full$m_pred[t, ])
  }
  
  dW <- max(sapply(1:(TT - 1), function(t) op_norm2(W_hat_list[[t]] - W_use[[t]])))
  dmu2 <- mean(rowSums((yhat_oracle - yhat_plugin)^2))
  
  data.frame(alpha = a, dW = dW, mean_sq_diff = dmu2)
})

print(sens_df)

# ======================================================================
# 7) Aggregation check (stationary weights from average network)
# ======================================================================
W_bar <- Reduce("+", W_use) / length(W_use)
pi_hat <- stationary_dist(W_bar)

agg_y <- as.numeric(Y_curr %*% pi_hat)
agg_y_lag <- c(sum(pi_hat * Y[1, ]), agg_y[1:(length(agg_y) - 1)])

agg_pred <- numeric(TT - 1)
for (t in 1:(TT - 1)) {
  b0 <- fit_full$m_pred[t, 1]
  b1 <- fit_full$m_pred[t, 2]
  b2 <- fit_full$m_pred[t, 3]
  agg_pred[t] <- b0 + (b1 + b2) * agg_y_lag[t]
}
cat(sprintf("Aggregation recursion MAE: %.4f\n", mean(abs(agg_y - agg_pred))))

# ======================================================================
# 8) Stability diagnostics (true vs estimated spectral radius)
# ======================================================================
rho_hat <- sapply(1:(TT - 1), function(t) {
  spectral_radius(fit_full$m_filt[t, 2] * W_use[[t]] + fit_full$m_filt[t, 3] * diag(N))
})

diag_df <- tibble::tibble(t = 2:TT, rho_true = rho_true[2:TT], rho_hat = rho_hat)

p1 <- ggplot(diag_df, aes(x = t)) +
  geom_line(aes(y = rho_true, linetype = "true")) +
  geom_line(aes(y = rho_hat,  linetype = "estimated")) +
  geom_hline(yintercept = 1, colour = "red") +
  labs(y = "spectral radius of beta1*W + beta2*I", linetype = "", x = "time")

p2 <- ggplot(sens_df, aes(x = dW, y = mean_sq_diff)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  labs(x = "sup_t ||W_hat - W||_op", y = "mean ||mu_oracle - mu_plugin||^2")

print(p1)
print(p2)

# Optional: save figures
ggsave("fig_stability.pdf", p1, width = 7, height = 4)
ggsave("fig_network_sensitivity.pdf", p2, width = 7, height = 4)


















# -------------------------------------------------------------------------
# A1) Multi-step horizons h = 1,2,4,8  (oracle W at forecast times)
#     "origin t uses filtered state at time t" with correct indexing.
#     Robust to Y being N x TT (will transpose to TT x N).
# -------------------------------------------------------------------------

Y_full <- as.matrix(Y)      # your simulation object (should be TT x N)
W_full <- W_list            # list of W_t

stopifnot(is.list(W_full), length(W_full) >= 2)

# --- Use W_list as the source of truth for TT and N ---
TT <- length(W_full)
N  <- nrow(W_full[[1]])

# Validate W_list
stopifnot(all(vapply(W_full, function(W) {
  (is.matrix(W) || inherits(W, "Matrix")) && all(dim(W) == c(N, N))
}, logical(1))))

# Ensure Y_full is TT x N (if user supplied N x TT, transpose)
if (nrow(Y_full) == TT && ncol(Y_full) == N) {
  # OK
} else if (nrow(Y_full) == N && ncol(Y_full) == TT) {
  Y_full <- t(Y_full)
} else {
  stop(sprintf(
    "Y has dim %d x %d, but expected TT x N = %d x %d (or N x TT = %d x %d).",
    nrow(Y_full), ncol(Y_full), TT, N, N, TT
  ))
}

# Filtered states (expected (TT-1) x K with rows corresponding to times 2..TT,
# but we also support TT x K if your code stored times 1..TT)
theta_full_filt  <- as.matrix(fit_full$m_filt)    # K=3 (b0,b1,b2)
theta_nonet_filt <- as.matrix(fit_nonet$m_filt)   # K=2 (b0,b2)

# Helper: grab filtered theta at time t (post-observation at time t)
theta_at <- function(t, theta_mat, TT) {
  nr <- nrow(theta_mat)
  if (nr == TT) {
    idx <- t
  } else if (nr == TT - 1) {
    idx <- t - 1   # row 1 corresponds to t=2
  } else {
    stop(sprintf("theta_mat nrow (%d) must be TT (%d) or TT-1 (%d).", nr, TT, TT-1))
  }
  if (idx < 1 || idx > nr) stop(sprintf("theta_at: time t=%d not available (idx=%d).", t, idx))
  theta_mat[idx, , drop = TRUE]
}

# Helper: grab W_t; supports list length TT (times 1..TT) or TT-1 (times 2..TT)
W_at <- function(t, W_list, TT) {
  L <- length(W_list)
  if (L == TT) {
    idx <- t
  } else if (L == TT - 1) {
    idx <- t - 1   # row 1 corresponds to t=2
  } else {
    stop(sprintf("W_list length (%d) must be TT (%d) or TT-1 (%d).", L, TT, TT-1))
  }
  if (idx < 1 || idx > L) stop(sprintf("W_at: time t=%d not available (idx=%d).", t, idx))
  W_list[[idx]]
}

# Full network multi-step recursion with oracle W_{t+1},...,W_{t+h}
gauss_fc_h <- function(t_origin, h, theta, W_use, TT) {
  stopifnot(h >= 1)
  b0 <- theta[1]; b1 <- theta[2]; b2 <- theta[3]
  y_prev <- as.numeric(Y_full[t_origin, ])
  for (j in seq_len(h)) {
    Wtj <- W_at(t_origin + j, W_use, TT)   # oracle future network at forecast time
    y_prev <- b0 + b1 * as.numeric(Wtj %*% y_prev) + b2 * y_prev
  }
  y_prev
}

# Non-network ablation (no W term)
gauss_fc_h_nonet <- function(t_origin, h, theta) {
  stopifnot(h >= 1)
  b0 <- theta[1]; b2 <- theta[2]
  y_prev <- as.numeric(Y_full[t_origin, ])
  for (j in seq_len(h)) {
    y_prev <- b0 + b2 * y_prev
  }
  y_prev
}

h_grid <- c(1,2,4,8)
h_max  <- max(h_grid)

# Safe evaluation window:
# - need t_origin >= 2 if theta_mat has TT-1 rows (times 2..TT)
# - need t_origin + h_max <= TT
t_eval_default_start <- 120
t_start <- max(2, min(t_eval_default_start, TT - h_max))
t_end   <- TT - h_max

if (t_start > t_end) {
  # fallback if TT is small
  t_start <- max(2, floor(0.6 * TT))
  t_end   <- TT - h_max
}
if (t_start > t_end) stop("Not enough time points to evaluate multi-step horizons with your TT/h_grid.")

t_eval <- t_start:t_end

eval_horizon <- function(h) {
  stopifnot(h >= 1)
  stopifnot(max(t_eval) + h <= TT)
  
  targets <- t_eval + h
  y_true  <- Y_full[targets, , drop = FALSE]
  
  yhat_full   <- matrix(NA_real_, nrow=length(t_eval), ncol=N)
  yhat_non    <- matrix(NA_real_, nrow=length(t_eval), ncol=N)
  yhat_static <- matrix(NA_real_, nrow=length(t_eval), ncol=N)
  
  if (!exists("beta_static")) {
    stop("beta_static not found. Fit your static OLS first and store coefficients in beta_static (length 3).")
  }
  
  for (k in seq_along(t_eval)) {
    t0 <- t_eval[k]
    
    th_full <- theta_at(t0, theta_full_filt, TT)
    th_non  <- theta_at(t0, theta_nonet_filt, TT)
    
    yhat_full[k,] <- gauss_fc_h(t0, h, th_full, W_full, TT)
    yhat_non[k,]  <- gauss_fc_h_nonet(t0, h, th_non)
    
    # Static OLS multi-step recursion (oracle W at forecast times)
    y_prev <- as.numeric(Y_full[t0, ])
    for (j in seq_len(h)) {
      Wtj <- W_at(t0 + j, W_full, TT)
      y_prev <- beta_static[1] +
        beta_static[2] * as.numeric(Wtj %*% y_prev) +
        beta_static[3] * y_prev
    }
    yhat_static[k,] <- y_prev
  }
  
  loss_full   <- rowMeans((y_true - yhat_full)^2)
  loss_nonet  <- rowMeans((y_true - yhat_non)^2)
  loss_static <- rowMeans((y_true - yhat_static)^2)
  
  d_full_vs_non <- loss_full - loss_nonet
  dm <- dm_test_nw(loss_full, loss_nonet, h = h)
  ci <- block_boot_ci(d_full_vs_non, B = 2000, block_len = max(4, h), seed = 1)
  
  tibble::tibble(
    h = h,
    MSE_full   = mean(loss_full),
    MSE_nonet  = mean(loss_nonet),
    MSE_static = mean(loss_static),
    diff_full_minus_nonet = mean(d_full_vs_non),
    diff_CI_lo = ci[1],
    diff_CI_hi = ci[2],
    DM_stat = dm$stat,
    DM_p    = dm$pval
  )
}

tab_sim_horizons <- dplyr::bind_rows(lapply(h_grid, eval_horizon))
print(tab_sim_horizons)

# -------------------------------------------------------------------------
# B) Distribution across nodes (referee request)
# -------------------------------------------------------------------------
node_improve <- function(h) {
  targets <- t_eval + h
  y_true  <- Y_full[targets, , drop = FALSE]
  
  yhat_full <- matrix(NA_real_, nrow=length(t_eval), ncol=N)
  yhat_non  <- matrix(NA_real_, nrow=length(t_eval), ncol=N)
  
  for (k in seq_along(t_eval)) {
    t0 <- t_eval[k]
    yhat_full[k,] <- gauss_fc_h(
      t0, h,
      theta_at(t0, theta_full_filt, TT),
      W_full, TT
    )
    yhat_non[k,] <- gauss_fc_h_nonet(
      t0, h,
      theta_at(t0, theta_nonet_filt, TT)
    )
  }
  
  mse_full_i <- colMeans((y_true - yhat_full)^2)
  mse_non_i  <- colMeans((y_true - yhat_non)^2)
  
  tibble::tibble(node = seq_len(N), h = h, dMSE = mse_full_i - mse_non_i)
}

df_nodes <- dplyr::bind_rows(lapply(h_grid, node_improve))

print(df_nodes %>%
        dplyr::group_by(h) %>%
        dplyr::summarize(
          frac_improve = mean(dMSE < 0),
          median_dMSE  = median(dMSE),
          .groups = "drop"
        ))

ggplot(df_nodes, aes(x = factor(h), y = dMSE)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    x = "Horizon h",
    y = "Nodewise dMSE (full - nonet)",
    title = "Simulation: distribution of network benefit across nodes"
  ) +
  theme_minimal()
