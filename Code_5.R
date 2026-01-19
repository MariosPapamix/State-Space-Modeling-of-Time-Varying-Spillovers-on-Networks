## ============================================================
## SSNR "reviewer-upgrade" evaluation & diagnostics (R)
## Multi-step rolling forecasts + uncertainty + PIT + stress tests
## ============================================================

## Packages
pkgs <- c("data.table","ggplot2","Matrix","RSpectra","scales")
to_install <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]
if(length(to_install)) install.packages(to_install)
lapply(pkgs, library, character.only=TRUE)

## -----------------------------
## 0) YOU PROVIDE THESE OBJECTS
## -----------------------------
## Y: T x N numeric matrix (Gaussian GDP) or integer matrix (Poisson Chicago)
## W_list: list length T, each element N x N matrix
## fit: list with either posterior means or draws. Minimal required fields:
##   fit$beta0_hat (T), fit$beta1_hat (T), fit$beta2_hat (T)
##   fit$q_sd: length-3 innovation SDs for RW on (beta0,beta1,beta2) (or 3xT)
##   fit$obs_sd: (Gaussian) scalar or length-N obs sd; optional
##
## For Poisson:
##   same betas, but obs model is Poisson(log link); no obs_sd needed.

stopifnot(is.matrix(Y), length(W_list) == nrow(Y))

T_ <- nrow(Y); N <- ncol(Y)
I_N <- Diagonal(N)

## -----------------------------
## 1) Linear algebra diagnostics
## -----------------------------
op_norm <- function(A) {
  ## Largest singular value (operator norm)
  ## RSpectra is fast for big N; fall back to svd for small N
  if(nrow(A) <= 200) {
    return(max(svd(as.matrix(A), nu=0, nv=0)$d))
  } else {
    return(RSpectra::svds(A, k=1, nu=0, nv=0)$d[1])
  }
}

spectral_radius <- function(A) {
  ## Largest modulus eigenvalue
  if(nrow(A) <= 200) {
    return(max(Mod(eigen(as.matrix(A), only.values=TRUE)$values)))
  } else {
    ev <- RSpectra::eigs(A, k=1, which="LM", opts=list(retvec=FALSE))$values
    return(max(Mod(ev)))
  }
}

## -----------------------------
## 2) Block bootstrap CI (for forecast-origin dependence)
## -----------------------------
block_boot_ci <- function(x, B=2000, block_len=6, conf=0.95, seed=1) {
  set.seed(seed)
  x <- as.numeric(x)
  n <- length(x)
  nb <- ceiling(n / block_len)
  idx_blocks <- lapply(1:(n - block_len + 1), function(s) s:(s+block_len-1))
  
  boot_mean <- replicate(B, {
    starts <- sample(1:length(idx_blocks), nb, replace=TRUE)
    idx <- unlist(idx_blocks[starts])[1:n]
    mean(x[idx], na.rm=TRUE)
  })
  
  alpha <- (1-conf)/2
  c(lo=quantile(boot_mean, alpha, na.rm=TRUE),
    hi=quantile(boot_mean, 1-alpha, na.rm=TRUE))
}

## -----------------------------
## 3) Monte Carlo predictive log score (stable)
## -----------------------------
logmeanexp <- function(v) {
  m <- max(v)
  m + log(mean(exp(v - m)))
}

## -----------------------------
## 4) Forecast simulators (Gaussian / Poisson)
##    These do rolling origins and multi-step horizons.
## -----------------------------
get_qsd_t <- function(fit, t) {
  ## allow constant q_sd (length 3) or time-varying (3 x T)
  if(is.matrix(fit$q_sd)) fit$q_sd[,t] else fit$q_sd
}

get_beta_hat_t <- function(fit, t) {
  c(b0=fit$beta0_hat[t], b1=fit$beta1_hat[t], b2=fit$beta2_hat[t])
}

W_scenario_at <- function(W_list, t, origin, scenario=c("oracle","carry_forward")) {
  scenario <- match.arg(scenario)
  if(scenario == "oracle") {
    W_list[[t]]
  } else {
    ## unknown future W: carry-forward W_{origin}
    W_list[[origin]]
  }
}

## ---- 4a) Gaussian SSNR simulator: Y_{t} = b0 + b1*W_t Y_{t-1} + b2*Y_{t-1} + eps
sim_forecast_gaussian <- function(Y, W_list, fit, origin, h,
                                  scenario=c("oracle","carry_forward"),
                                  S=1000, seed=1) {
  set.seed(seed + origin + 100*h)
  scenario <- match.arg(scenario)
  
  y_prev <- as.numeric(Y[origin, ])
  beta0 <- fit$beta0_hat[origin]
  beta1 <- fit$beta1_hat[origin]
  beta2 <- fit$beta2_hat[origin]
  qsd   <- get_qsd_t(fit, origin)  # (b0,b1,b2) RW innovation SDs
  
  ## observation noise
  obs_sd <- if(!is.null(fit$obs_sd)) fit$obs_sd else 1.0
  if(length(obs_sd) == 1) obs_sd <- rep(obs_sd, N)
  
  ## store S predictive draws for Y_{origin+h}
  y_draws <- matrix(NA_real_, nrow=S, ncol=N)
  
  for(s in 1:S) {
    b0 <- beta0; b1 <- beta1; b2 <- beta2
    y  <- y_prev
    for(k in 1:h) {
      ## propagate coefficient RW
      b0 <- b0 + rnorm(1, 0, qsd[1])
      b1 <- b1 + rnorm(1, 0, qsd[2])
      b2 <- b2 + rnorm(1, 0, qsd[3])
      
      t_next <- origin + k
      Wt <- W_scenario_at(W_list, t_next, origin, scenario=scenario)
      mu <- as.numeric(b0 + (b1 * (Wt %*% y)) + b2 * y)
      y  <- mu + rnorm(N, 0, obs_sd)
    }
    y_draws[s,] <- y
  }
  
  list(
    mean = colMeans(y_draws),
    draws = y_draws
  )
}

## ---- 4b) Poisson SSNR simulator: Y ~ Poisson(exp(eta)), eta = b0 + b1*W Y_{t-1} + b2*Y_{t-1}
sim_forecast_poisson <- function(Y, W_list, fit, origin, h,
                                 scenario=c("oracle","carry_forward"),
                                 S=1000, seed=1) {
  set.seed(seed + origin + 100*h)
  scenario <- match.arg(scenario)
  
  y_prev <- as.numeric(Y[origin, ])
  beta0 <- fit$beta0_hat[origin]
  beta1 <- fit$beta1_hat[origin]
  beta2 <- fit$beta2_hat[origin]
  qsd   <- get_qsd_t(fit, origin)
  
  y_draws <- matrix(NA_integer_, nrow=S, ncol=N)
  
  for(s in 1:S) {
    b0 <- beta0; b1 <- beta1; b2 <- beta2
    y  <- y_prev
    for(k in 1:h) {
      b0 <- b0 + rnorm(1, 0, qsd[1])
      b1 <- b1 + rnorm(1, 0, qsd[2])
      b2 <- b2 + rnorm(1, 0, qsd[3])
      
      t_next <- origin + k
      Wt <- W_scenario_at(W_list, t_next, origin, scenario=scenario)
      eta <- as.numeric(b0 + (b1 * (Wt %*% y)) + b2 * y)
      lam <- pmax(exp(eta), 1e-12)
      y   <- rpois(N, lam)
    }
    y_draws[s,] <- y
  }
  
  list(
    mean = colMeans(y_draws),
    draws = y_draws
  )
}

## -----------------------------
## 5) Rolling-origin evaluation (multi-step + uncertainty)
## -----------------------------
eval_rolling <- function(Y, W_list, fit, origins, horizons=c(1,2,4,8),
                         scenario=c("oracle","carry_forward"),
                         model=c("gaussian","poisson"),
                         S=1000, seed=1, block_len=6) {
  scenario <- match.arg(scenario)
  model <- match.arg(model)
  
  rows <- list()
  
  for(h in horizons) {
    per_origin_metric <- data.table(origin=origins,
                                    mse=NA_real_, mae=NA_real_,
                                    logscore=NA_real_, cov90=NA_real_)
    for(ii in seq_along(origins)) {
      o <- origins[ii]
      y_true <- as.numeric(Y[o+h, ])
      
      fc <- if(model=="gaussian") {
        sim_forecast_gaussian(Y, W_list, fit, origin=o, h=h, scenario=scenario, S=S, seed=seed)
      } else {
        sim_forecast_poisson(Y, W_list, fit, origin=o, h=h, scenario=scenario, S=S, seed=seed)
      }
      
      y_mean <- fc$mean
      draws  <- fc$draws
      
      per_origin_metric$mse[ii] <- mean((y_true - y_mean)^2)
      per_origin_metric$mae[ii] <- mean(abs(y_true - y_mean))
      
      ## predictive 90% interval coverage (empirical)
      q05 <- apply(draws, 2, quantile, probs=0.05)
      q95 <- apply(draws, 2, quantile, probs=0.95)
      per_origin_metric$cov90[ii] <- mean(y_true >= q05 & y_true <= q95)
      
      ## log score:
      ## Gaussian: use empirical predictive density via draws (kernel-free via Normal approx per node)
      ## Poisson: use Monte Carlo mixture exactly via draws' implied lambdas is hard if you only stored Y draws
      ## => use draws as predictive distribution and score by empirical pmf:
      if(model=="poisson") {
        ## approximate log p(y) by log(mean(1{Y_draw==y})) with smoothing
        ## better: store lambdas if you have them; this is the minimal “works now” version.
        p_hat <- pmax(colMeans(draws == matrix(y_true, nrow=nrow(draws), ncol=N, byrow=TRUE)), 1e-6)
        per_origin_metric$logscore[ii] <- sum(log(p_hat))
      } else {
        ## normal approx from draws
        mu <- y_mean
        sd <- apply(draws, 2, sd); sd <- pmax(sd, 1e-6)
        per_origin_metric$logscore[ii] <- sum(dnorm(y_true, mean=mu, sd=sd, log=TRUE))
      }
    }
    
    ## aggregate across origins
    mse_bar <- mean(per_origin_metric$mse)
    mae_bar <- mean(per_origin_metric$mae)
    ls_bar  <- mean(per_origin_metric$logscore)
    cov_bar <- mean(per_origin_metric$cov90)
    
    ## uncertainty on horizon-wise mean metrics (block bootstrap over origins)
    ci_mse <- block_boot_ci(per_origin_metric$mse, block_len=block_len, seed=seed)
    ci_mae <- block_boot_ci(per_origin_metric$mae, block_len=block_len, seed=seed)
    ci_ls  <- block_boot_ci(per_origin_metric$logscore, block_len=block_len, seed=seed)
    
    rows[[length(rows)+1]] <- data.table(
      horizon=h, scenario=scenario, model=model,
      mse=mse_bar, mse_lo=ci_mse["lo"], mse_hi=ci_mse["hi"],
      mae=mae_bar, mae_lo=ci_mae["lo"], mae_hi=ci_mae["hi"],
      logscore=ls_bar, logscore_lo=ci_ls["lo"], logscore_hi=ci_ls["hi"],
      cov90=cov_bar
    )
  }
  
  rbindlist(rows)
}

## -----------------------------
## 6) “Effect size is small” 
## -----------------------------
eval_diff <- function(Y, W_list, fit_net, fit_nonet, origins, horizons=c(1,2,4,8),
                      model=c("gaussian","poisson"), S=1000, seed=1, block_len=6) {
  model <- match.arg(model)
  out <- list()
  
  for(h in horizons) {
    d_mae <- numeric(length(origins))
    d_ls  <- numeric(length(origins))
    
    for(i in seq_along(origins)) {
      o <- origins[i]
      y_true <- as.numeric(Y[o+h, ])
      
      fc_net <- if(model=="gaussian") {
        sim_forecast_gaussian(Y,W_list,fit_net,origin=o,h=h,scenario="oracle",S=S,seed=seed)
      } else {
        sim_forecast_poisson(Y,W_list,fit_net,origin=o,h=h,scenario="oracle",S=S,seed=seed)
      }
      fc_nn <- if(model=="gaussian") {
        sim_forecast_gaussian(Y,W_list,fit_nonet,origin=o,h=h,scenario="oracle",S=S,seed=seed)
      } else {
        sim_forecast_poisson(Y,W_list,fit_nonet,origin=o,h=h,scenario="oracle",S=S,seed=seed)
      }
      
      ## MAE difference (net - no-net)
      d_mae[i] <- mean(abs(y_true - fc_net$mean)) - mean(abs(y_true - fc_nn$mean))
      
      ## Log score difference (net - no-net)
      if(model=="poisson") {
        p_net <- pmax(colMeans(fc_net$draws == matrix(y_true, nrow=nrow(fc_net$draws), ncol=N, byrow=TRUE)), 1e-6)
        p_nn  <- pmax(colMeans(fc_nn$draws  == matrix(y_true, nrow=nrow(fc_nn$draws),  ncol=N, byrow=TRUE)), 1e-6)
        d_ls[i] <- sum(log(p_net)) - sum(log(p_nn))
      } else {
        muN <- fc_net$mean; sdN <- pmax(apply(fc_net$draws,2,sd),1e-6)
        mu0 <- fc_nn$mean;  sd0 <- pmax(apply(fc_nn$draws,2,sd),1e-6)
        d_ls[i] <- sum(dnorm(y_true, muN, sdN, log=TRUE)) - sum(dnorm(y_true, mu0, sd0, log=TRUE))
      }
    }
    
    ci_mae <- block_boot_ci(d_mae, block_len=block_len, seed=seed)
    ci_ls  <- block_boot_ci(d_ls,  block_len=block_len, seed=seed)
    
    out[[length(out)+1]] <- data.table(
      horizon=h,
      dMAE=mean(d_mae), dMAE_lo=ci_mae["lo"], dMAE_hi=ci_mae["hi"],
      dLogScore=mean(d_ls), dLogScore_lo=ci_ls["lo"], dLogScore_hi=ci_ls["hi"]
    )
  }
  
  rbindlist(out)
}

## -----------------------------
## 7) Randomized PIT (Poisson calibration figure)
## -----------------------------
randomized_pit <- function(y_true, y_draws, seed=1) {
  ## y_true: length N
  ## y_draws: S x N predictive draws
  set.seed(seed)
  S <- nrow(y_draws); N <- length(y_true)
  u <- runif(N)
  
  pit <- numeric(N)
  for(i in 1:N) {
    yy <- y_draws[,i]
    pit[i] <- mean(yy < y_true[i]) + u[i]*mean(yy == y_true[i])
  }
  pit
}

## -----------------------------
## 8) Network perturbation stress test (Chicago “sensitivity curve”)
## -----------------------------
perturb_W <- function(W, eps=0.1, seed=1) {
  ## randomly drop eps fraction of edges then renormalize rows
  set.seed(seed)
  A <- as.matrix(W)
  nz <- which(A != 0, arr.ind=TRUE)
  drop_n <- floor(eps * nrow(nz))
  if(drop_n > 0) {
    drop_idx <- nz[sample(1:nrow(nz), drop_n), , drop=FALSE]
    A[drop_idx] <- 0
  }
  rs <- rowSums(A)
  rs[rs==0] <- 1
  A <- A / rs
  A
}

stress_test_network <- function(Y, W_list, fit, origins, h=4, eps_grid=seq(0,0.5,by=0.05),
                                model=c("poisson","gaussian"), S=500, seed=1) {
  model <- match.arg(model)
  baseW <- W_list[[origins[1]]]  
  out <- lapply(eps_grid, function(eps) {
    Wp <- perturb_W(baseW, eps=eps, seed=seed)
    ## operator-norm error
    dW <- op_norm(Wp - baseW)
    
    ## evaluate MAE at horizon h with perturbed W used as "oracle"
    Wp_list <- W_list
    for(t in 1:length(Wp_list)) Wp_list[[t]] <- Wp
    
    tab <- eval_rolling(Y, Wp_list, fit, origins=origins, horizons=c(h),
                        scenario="oracle", model=model, S=S, seed=seed)
    data.table(eps=eps, op_err=dW, mae=tab$mae)
  })
  rbindlist(out)
}

## -----------------------------
## 9) Stability figure: spectral radius + operator norm over time
## -----------------------------
stability_over_time <- function(W_list, fit) {
  rho <- numeric(T_)
  on  <- numeric(T_)
  for(t in 2:T_) {
    b <- get_beta_hat_t(fit, t)
    Bt <- b["b1"] * W_list[[t]] + b["b2"] * I_N
    rho[t] <- spectral_radius(Bt)
    on[t]  <- op_norm(Bt)
  }
  data.table(t=1:T_, spectral_radius=rho, op_norm=on)
}

## ============================================================
## 10) EXAMPLE CALLS (edit indices to match your application)
## ============================================================

## ----- GDP (Gaussian) -----
## Longer holdout + multi-step horizons (addresses “short window” critique)
## Example: evaluate rolling one-step and multi-step over last 40 quarters:
gdp_origins <- (T_-40):(T_-8)  # ensures origin+h <= T_ for h up to 8
gdp_tab_oracle <- eval_rolling(Y, W_list, fit, origins=gdp_origins,
                               horizons=c(1,2,4,8),
                               scenario="oracle", model="gaussian",
                               S=1000, block_len=8)
gdp_tab_carry  <- eval_rolling(Y, W_list, fit, origins=gdp_origins,
                               horizons=c(1,2,4,8),
                               scenario="carry_forward", model="gaussian",
                               S=1000, block_len=8)

gdp_tab <- rbind(gdp_tab_oracle, gdp_tab_carry)
print(gdp_tab)

## Plot “Figure 5”-style: multi-step MSE oracle vs carry-forward
ggplot(gdp_tab, aes(x=horizon, y=mse, linetype=scenario)) +
  geom_line() + geom_point() +
  scale_y_continuous(labels=label_number()) +
  labs(title="GDP: multi-step test MSE (oracle vs carry-forward W)",
       x="Horizon h", y="Test MSE") +
  theme_minimal()

## Stability diagnostic (spectral radius + operator norm)
stab <- stability_over_time(W_list, fit)
ggplot(stab[!is.na(spectral_radius)], aes(x=t)) +
  geom_line(aes(y=spectral_radius, color="spectral radius")) +
  geom_line(aes(y=op_norm, color="operator norm")) +
  geom_hline(yintercept=1, linetype=2) +
  labs(title="Stability diagnostics: spectral radius and operator norm",
       x="time", y="value", color="") +
  theme_minimal()

## ----- Chicago (Poisson) -----
## Multi-step + uncertainty on differences + calibration + stress test
chi_origins <- (T_-36):(T_-8)  # e.g., last 3 years, enough for h up to 8
chi_tab <- eval_rolling(Y, W_list, fit, origins=chi_origins,
                        horizons=c(1,2,4,8),
                        scenario="oracle", model="poisson",
                        S=1000, block_len=6)
print(chi_tab)

## If you have a no-network baseline fit (set beta1_hat=0, q_sd[2]=0)
## build fit_nonet by copying fit and forcing no network spillover:
fit_nonet <- fit
fit_nonet$beta1_hat <- rep(0, length(fit$beta1_hat))
if(is.matrix(fit_nonet$q_sd)) fit_nonet$q_sd[2,] <- 0 else fit_nonet$q_sd[2] <- 0

## Paired uncertainty on differences (this is what kills “0.001 MAE” critiques)
chi_diff <- eval_diff(Y, W_list, fit_net=fit, fit_nonet=fit_nonet,
                      origins=chi_origins, horizons=c(1,2,4,8),
                      model="poisson", S=800, block_len=6)
print(chi_diff)

## Randomized PIT histogram (“Figure 6b”-style)
## Use 1-step forecasts pooled across origins and nodes:
pit_vals <- c()
for(o in chi_origins) {
  fc1 <- sim_forecast_poisson(Y,W_list,fit,origin=o,h=1,scenario="oracle",S=800)
  pit_vals <- c(pit_vals, randomized_pit(y_true=as.numeric(Y[o+1,]), y_draws=fc1$draws))
}
ggplot(data.table(pit=pit_vals), aes(x=pit)) +
  geom_histogram(bins=20) +
  labs(title="Randomized PIT (1-step) – Poisson SSNR", x="PIT", y="count") +
  theme_minimal()

## Network perturbation sensitivity curve (“Figure 6c”-style)
stress <- stress_test_network(Y, W_list, fit, origins=chi_origins, h=4,
                              eps_grid=seq(0,0.5,by=0.05),
                              model="poisson", S=400)
ggplot(stress, aes(x=op_err, y=mae)) +
  geom_line() + geom_point() +
  labs(title="Sensitivity of multi-step MAE to network perturbations",
       x="||W_pert - W||_op", y="MAE (h=4)") +
  theme_minimal()
