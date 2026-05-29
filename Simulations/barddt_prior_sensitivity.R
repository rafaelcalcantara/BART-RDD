# =============================================================================
# BARDDT sig.prior Sensitivity Analysis
# DGPs: 1 and 3 | sig.prior: 0.01, 0.1, 1 | 50 iterations
# =============================================================================
# Key modification vs. original simulation_estimator_functions.R:
#   sigma2_leaf_init = diag(rep(sig.prior/150, 4))   [was hardcoded 0.1/150]
# =============================================================================

library(parallel)
library(MASS)

# ── Global settings ───────────────────────────────────────────────────────────
no_cores         <- 10
n                <- 4000
s                <- 100
c_val            <- 0       # cutoff (using c_val to avoid clash with base::c)
Owidth           <- 0.1
sig.prior.values <- c(0.01, 0.05, 0.1, 0.5, 1)

# ── DGP definitions (DGP 1 = row 1, DGP 3 = row 3 of original matrix) ────────
dgp_mat <- rbind(
  c(k1 = 1, k2 = 1, k3 = 0, k4 = 0.1, k5 = 0, p = 2, rho = 0.5),  # DGP 1
  c(k1 = 1, k2 = 1, k3 = 0, k4 = 0.1, k5 = 1, p = 2, rho = 0.0)   # DGP 3
)
dgp_labels <- c("DGP1", "DGP3")

# ── Directory setup ───────────────────────────────────────────────────────────
for (d in c("Data", "Results", "Results/RMSE", "Time")) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# ── Pre-generate fixed w matrices (seed 007, one per unique p) ───────────────
# Mirrors the set.seed(007) / mvrnorm block in the original master scripts.
p_vals <- unique(dgp_mat[, "p"])
w_list <- setNames(
  lapply(p_vals, function(p) {
    K <- 2 * toeplitz(seq(1, 0, length.out = p))
    set.seed(007)
    MASS::mvrnorm(n, rep(0, p), K)
  }),
  as.character(p_vals)
)

# ── Worker function: one simulation sample ────────────────────────────────────
# Everything needed is passed as arguments so parLapply can serialise it
# cleanly without relying on a shared global environment.
run_one_sample <- function(sample, dgp_params, w_fixed, c_val, Owidth, sig.prior) {

  library(MASS)
  library(stochtree)
  library(rdrobust)

  # Unpack DGP parameters
  k1  <- dgp_params[["k1"]];  k2  <- dgp_params[["k2"]]
  k3  <- dgp_params[["k3"]];  k4  <- dgp_params[["k4"]]
  k5  <- dgp_params[["k5"]];  p   <- dgp_params[["p"]]
  rho <- dgp_params[["rho"]]
  n_obs <- nrow(w_fixed)

  # ── Derived quantities (mirrors simulation_data.R) ──────────────────────────
  m        <- 1
  beta_raw <- rep(1, p) / sqrt(p)
  K        <- 2 * toeplitz(seq(1, 0, length.out = p))

  # When rho = 0 the covariate contribution vanishes
  if (rho == 0) {
    beta <- rep(0, p)
  } else {
    beta <- rho * beta_raw / sqrt(as.numeric(beta_raw %*% K %*% beta_raw))
  }

  wstar  <- function(wval) rowSums(wval) / sqrt(ncol(wval))
  mu0.x  <- function(x)    (x + 1)^3
  mu0.w  <- function(wval) (wstar(wval) + 2)^2
  mu_fn  <- function(xval, wval, k1, k3, sf)
    sf * (k1 * mu0.x(xval) + mu0.w(wval) * (sign(xval + 1) * sqrt(abs(xval + 1)))^k3)
  tau0   <- function(wv)   pnorm(2 * wv[, 1] + 3, 0, 1) / 2 + dnorm(wv[, 1], 0, 1)
  tau_fn <- function(wval, minval, k2_adj, tau.bar)
    k2_adj * (tau0(wval) - tau.bar) + minval

  # Population calibration — set seed FIRST to exactly mirror fit_general():
  #   set.seed(sample); source("simulation_data.R")  [wcond drawn inside]
  set.seed(sample)
  Kcond   <- K - K %*% beta %*% t(beta) %*% K
  wcond   <- MASS::mvrnorm(10000, -beta * m, Kcond)
  sf      <- 1 / sd(mu_fn(c_val, wcond, k1, k3, 1))
  k2.new  <- k2 / sd(tau0(wcond))
  tau.bar <- mean(tau0(wcond))
  mintau  <- -k2.new * min(tau0(wcond) - tau.bar) + k5
  sigma_y <- k4

  # ── Sample data (w is fixed; x, y, z are drawn from the seeded state) ───────
  x    <- rnorm(n_obs, m + w_fixed %*% beta, sqrt(max(1 - rho^2, 0)))
  z    <- as.numeric(x >= c_val)
  y    <- mu_fn(x, w_fixed, k1, k3, sf) +
    tau_fn(w_fixed, mintau, k2.new, tau.bar) * z +
    rnorm(n_obs, 0, sigma_y)
  cate <- tau_fn(w_fixed, mintau, k2.new, tau.bar)

  # ── Test window ───────────────────────────────────────────────────────────────
  test <- (c_val - Owidth) <= x & x <= (c_val + Owidth)

  # ── BARDDT (modified: sigma2_leaf_init uses sig.prior instead of 0.1) ────────
  global.parms <- list(
    standardize         = TRUE,
    sample_sigma_global = TRUE,
    sigma2_global_init  = 0.1
  )
  mean.parms <- list(
    num_trees          = 50,
    min_samples_leaf   = 20,
    alpha              = 0.95,
    beta               = 2,
    max_depth          = 20,
    sample_sigma2_leaf = FALSE,
    sigma2_leaf_init   = diag(rep(sig.prior / 150, 4))   # <-- KEY CHANGE
  )

  B  <- cbind(z * x,       (1 - z) * x, z,       rep(1, n_obs))
  B1 <- cbind(rep(c_val, n_obs), rep(0, n_obs),   rep(1, n_obs), rep(1, n_obs))
  B0 <- cbind(rep(0, n_obs),     rep(c_val, n_obs), rep(0, n_obs), rep(1, n_obs))

  barddt_fit <- stochtree::bart(
    X_train            = as.matrix(cbind(x, w_fixed)),
    y_train            = y,
    leaf_basis_train   = B,
    mean_forest_params = mean.parms,
    general_params     = global.parms,
    num_mcmc           = 1000,
    num_gfr            = 30
  )

  xmat_test <- as.matrix(cbind(rep(0, n_obs), w_fixed)[test, ])
  pred1     <- predict(barddt_fit, xmat_test, B1[test, ])$y_hat
  pred0     <- predict(barddt_fit, xmat_test, B0[test, ])$y_hat
  post      <- pred1 - pred0   # matrix: n_test_obs × n_mcmc_draws

  # ── RMSE ──────────────────────────────────────────────────────────────────────
  cate_test <- cate[test]
  sqrt(mean((rowMeans(post) - cate_test)^2))
}

# ── Main loop: DGPs × sig.prior values ───────────────────────────────────────
all_results <- list()

for (dgp_idx in seq_len(nrow(dgp_mat))) {

  params  <- as.list(dgp_mat[dgp_idx, ])
  lbl     <- dgp_labels[dgp_idx]
  w_fixed <- w_list[[as.character(params[["p"]])]]

  for (sp in sig.prior.values) {

    cat(sprintf("\n>>> %s | sig.prior = %g  (%d iterations) ...\n", lbl, sp, s))

    cl <- makeCluster(no_cores, type = "SOCK")

    # Export only what parLapply workers need
    clusterExport(cl,
                  varlist = c("run_one_sample", "params", "w_fixed", "c_val", "Owidth", "sp"),
                  envir   = environment()
    )

    rmse_vec <- parLapply(
      cl, 1:s,
      fun       = run_one_sample,
      dgp_params = params,
      w_fixed    = w_fixed,
      c_val      = c_val,
      Owidth     = Owidth,
      sig.prior  = sp
    )
    stopCluster(cl)

    key                <- paste0(lbl, "_sig", sp)
    all_results[[key]] <- unlist(rmse_vec)
    cat(sprintf("    Mean RMSE = %.4f  |  Median = %.4f\n",
                mean(all_results[[key]]), median(all_results[[key]])))

    # Persist individual RMSE files (mirrors original output structure)
    dgp_id <- paste(paste0(names(params), "_", unlist(params)), collapse = "_")
    rmse_dir <- file.path("Results", "RMSE",
                          paste0(dgp_id, "_sigprior_", sp))
    if (!dir.exists(rmse_dir)) dir.create(rmse_dir, recursive = TRUE)
    for (j in seq_along(all_results[[key]])) {
      write.table(all_results[[key]][j],
                  file.path(rmse_dir, paste0("barddt_sample_", j, ".csv")),
                  row.names = FALSE, col.names = FALSE)
    }
  }
}

# ── Boxplot: RMSE by sig.prior, faceted by DGP ───────────────────────────────
cols_sp <- c("0.01"  = "#E69F00",   # amber
             "0.05"  = "#D55E00",   # vermillion
             "0.1"   = "#56B4E9",   # sky blue
             "0.5"   = "#0072B2",   # deep blue
             "1"     = "#009E73")   # green

n_sp   <- length(sig.prior.values)
n_dgps <- length(dgp_labels)

# Build ordered list for boxplot()
bp_data  <- vector("list",  n_dgps * n_sp)
bp_names <- character(n_dgps * n_sp)
bp_cols  <- character(n_dgps * n_sp)

idx <- 1
for (dgp_lbl in dgp_labels) {
  for (sp in sig.prior.values) {
    key            <- paste0(dgp_lbl, "_sig", sp)
    bp_data[[idx]] <- all_results[[key]]
    bp_names[idx]  <- paste0(dgp_lbl, "\nsig=", sp)
    bp_cols[idx]   <- cols_sp[as.character(sp)]
    idx            <- idx + 1
  }
}

# Vertical dashed separator between DGPs
# Extra bottom margin to accommodate the horizontal legend
pdf("Results/Figures/boxplots_prior_sensitivity.pdf")
op <- par(mar = c(7.5, 5, 4.5, 2),bty="l")
boxplot(
  bp_data,
  col      = bp_cols,
  border   = "gray30",
  main     = "BARDDT prior sensitivity",
  ylab     = "RMSE",
  xaxt     = "n",
  las      = 1,
  cex.axis = 0.85,
  outline  = TRUE,
  whisklty = 1
)
axis(1,at=c(3,8),labels=c("DGP 1","DGP 3"),tick=FALSE)
# Vertical line separating the two DGP groups
abline(v = n_sp + 0.5, lty = 2, col = "gray50")
# Horizontal legend centred below the plot
legend(
  x      = "bottom",
  inset  = c(0, -0.25),   # shift down below the x-axis labels
  legend = sig.prior.values,
  fill   = cols_sp,
  bty    = "n",
  horiz  = TRUE,
  title  = "Prior scale",
  xpd    = TRUE            # allow drawing outside the plot region
)
par(op)
dev.off()
