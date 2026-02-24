#####################################################
### R Code from Fan Y and Bandyopadhyay D. (2026+). Individualized treatment rules for bivariate outcomes with 
### informative follow-up: Application to periodontal health
### Authors: Yiwei Fan and Dipankar Bandyopadhyay 
### Last edited: 02/23/2026
### "myfunctions.R" consists of a list of functions to be called in "analysis.R"
#######################################################

# Vector of required packages
required_packages <- c(
  "dplyr", "lme4", "survival", "splines2"
)


# Identify packages not yet installed
new_packages <- required_packages[!(
  required_packages %in% rownames(installed.packages())
)]

# Install missing packages
if (length(new_packages) > 0) {
  install.packages(new_packages, dependencies = TRUE)
}

# Load all required packages
invisible(lapply(required_packages, library, character.only = TRUE))




fit_stage1_empirical <- function(dat,
                                 y      = "COMP",
                                 id     = "STUDY_ID",
                                 covars = c("GENDER1","RACE1","RACE2","AGE1","DM_IND1",
                                            "DENT_COMMERCIAL1","TOBACCO_IND1","TOBACCO_IND2",
                                            "TX_IND1","BASE_CAL","BASE_PPD"),
                                 scale_continuous = TRUE) {
  
  keep <- c(id, y, covars)
  dat0 <- dat[, keep, drop = FALSE]
  dat0 <- dat0[!is.na(dat0[[y]]), , drop = FALSE]
  
  dat0[[id]] <- factor(dat0[[id]])
  
  # --- Standardization: apply z-scoring only to variables that appear continuous ---
  scaler <- NULL
  covars_use <- covars
  
  if (scale_continuous) {
    is_cont <- sapply(covars, function(v) {
      x <- dat0[[v]]
      if (!is.numeric(x)) return(FALSE)
      ux <- unique(x[is.finite(x)])
      length(ux) > 2   # treat as continuous if it has >2 distinct values
    })
    
    cont_vars <- covars[is_cont]
    if (length(cont_vars) > 0) {
      scaler <- lapply(cont_vars, function(v) {
        x <- dat0[[v]]
        m <- mean(x, na.rm = TRUE)
        s <- sd(x, na.rm = TRUE)
        if (!is.finite(s) || s == 0) s <- 1
        dat0[[paste0(v, "_z")]] <<- (x - m) / s
        list(var = v, mean = m, sd = s)
      })
      names(scaler) <- cont_vars
      
      # Replace original variables by standardized versions in the formula
      covars_use <- covars
      covars_use[match(cont_vars, covars_use)] <- paste0(cont_vars, "_z")
    }
  }
  
  form <- as.formula(
    paste0(y, " ~ ", paste(covars_use, collapse = " + "), " + (1 | ", id, ")")
  )
  
  fit_lmm <- lmer(form, data = dat0, REML = TRUE)
  
  re <- ranef(fit_lmm)[[id]]
  W_hat_df <- data.frame(W_hat = re[, "(Intercept)"], stringsAsFactors = FALSE)
  W_hat_df[[id]] <- rownames(re)
  W_hat_df <- W_hat_df[, c(id, "W_hat"), drop = FALSE]
  
  dat2 <- dat
  dat2[[id]] <- as.character(dat2[[id]])
  W_hat_df[[id]] <- as.character(W_hat_df[[id]])
  dat2$W_hat <- W_hat_df$W_hat[ match(dat2[[id]], W_hat_df[[id]]) ]
  
  list(dat = dat2, fit_lmm = fit_lmm, W_hat_df = W_hat_df, scaler = scaler)
}



build_intensity_data_from_visits_emp <- function(subdf,
                                                 id_col    = "STUDY_ID",
                                                 time_col  = NULL,    # prefer cumdelta; otherwise accumulate from delta
                                                 delta_col = "delta",
                                                 a_col     = "TX_IND1",
                                                 covars    = c("GENDER1","RACE1","RACE2","AGE1","DM_IND1",
                                                               "DENT_COMMERCIAL1","TOBACCO_IND1","TOBACCO_IND2","COMP",
                                                               "BASE_CAL","BASE_PPD"),
                                                 w_col     = "W_hat") {
  
  dat <- subdf
  dat$.rowid <- seq_len(nrow(dat))   # used to merge predictions back to the original data (most robust)
  
  # Choose time: prefer cumdelta; otherwise construct cumulative time from delta
  if (is.null(time_col)) {
    if ("cumdelta" %in% names(dat)) time_col <- "cumdelta"
  }
  if (is.null(time_col)) {
    # Construct cumulative time from delta
    dat <- dat[order(dat[[id_col]], dat[[delta_col]]), ]
    dat[[time_col]] <- NA_real_
    split_list <- split(dat, dat[[id_col]])
    dat_list <- lapply(split_list, function(d) {
      d[[time_col]] <- cumsum(as.numeric(d[[delta_col]]))
      d
    })
    dat <- do.call(rbind, dat_list)
    rownames(dat) <- NULL
  }
  
  dat <- dat[order(dat[[id_col]], dat[[time_col]]), ]
  
  split_list <- split(dat, dat[[id_col]])
  
  int_list <- lapply(split_list, function(d) {
    stop  <- as.numeric(d[[time_col]])
    start <- c(0, head(stop, -1))
    event <- rep(1L, length(stop))   # empirical: each row corresponds to a visit => all are events
    
    # Assemble output (including covariates, A, W_hat, and .rowid)
    out <- data.frame(
      id    = d[[id_col]],
      start = start,
      stop  = stop,
      event = event,
      rowid = d$.rowid
    )
    
    # Add covariates
    for (v in covars) out[[v]] <- d[[v]]
    out[[a_col]] <- d[[a_col]]
    out[[w_col]] <- d[[w_col]]
    
    # Drop invalid intervals (often due to baseline time=0 leading to start==stop)
    out <- out[is.finite(out$stop) & out$stop > out$start, , drop = FALSE]
    out
  })
  
  int_dat <- do.call(rbind, int_list)
  rownames(int_dat) <- NULL
  int_dat
}

library(survival)

fit_intensity_cox_emp <- function(int_dat,
                                  a_col  = "TX_IND1",
                                  w_col  = "W_hat",
                                  covars = c("GENDER1","RACE1","RACE2","AGE1","DM_IND1",
                                             "DENT_COMMERCIAL1","TOBACCO_IND1","TOBACCO_IND2","COMP",
                                             "BASE_CAL","BASE_PPD")) {
  
  rhs <- paste(
    unique(c(
      covars,
      a_col,
      w_col,
      paste0(a_col, ":", w_col),
      "cluster(id)"
    )),
    collapse = " + "
  )
  
  form <- as.formula(paste0("Surv(start, stop, event) ~ ", rhs))
  
  fit_cox <- coxph(form, data = int_dat, ties = "efron")
  list(intensity_data = int_dat, fit_cox = fit_cox)
}

add_eta_rate_to_visits_emp <- function(subdf, intensity_data, fit_cox) {
  lp <- as.numeric(predict(fit_cox, newdata = intensity_data, type = "lp"))
  intensity_data$eta_rate <- lp
  intensity_data$hr <- exp(lp)  # relative intensity exp(lp), sometimes more convenient
  
  out <- subdf
  if (!(".rowid" %in% names(out))) out$.rowid <- seq_len(nrow(out))
  out$eta_rate <- NA_real_
  out$hr <- NA_real_
  
  idx <- match(intensity_data$rowid, out$.rowid)
  out$eta_rate[idx] <- intensity_data$eta_rate
  out$hr[idx] <- intensity_data$hr
  
  out
}



fit_treatment <- function(visits) {
  fit_trt <- glm(
    TX_IND1 ~ GENDER1 + RACE1 + RACE2 + AGE1 + DM_IND1 + DENT_COMMERCIAL1 +
      TOBACCO_IND1 + TOBACCO_IND2 + COMP + BASE_CAL + BASE_PPD,
    family = binomial(),
    data = visits
  )
  visits$pi_hat <- predict(fit_trt, type = "response")
  list(dat = visits, fit_trt = fit_trt)
}


trim_by_quantile <- function(w, q = c(0.01, 0.99)) {
  qq <- quantile(w, probs = q, na.rm = TRUE)
  pmin(pmax(w, qq[1]), qq[2])
}


build_weights_simple <- function(dat,
                                 a_col  = "TX_IND1",
                                 pi_col = "pi_hat",
                                 hr_col = "hr",
                                 trim_q = c(0.01, 0.99),
                                 stabilized = FALSE) {
  out <- dat
  
  A  <- out[[a_col]]
  pi <- out[[pi_col]]
  hr <- out[[hr_col]]
  
  # treatment weight
  if (stabilized) {
    pA <- mean(A == 1, na.rm = TRUE)
    out$w_trt <- ifelse(A == 1, pA / pi, (1 - pA) / (1 - pi))
  } else {
    out$w_trt <- ifelse(A == 1, 1 / pi, 1 / (1 - pi))
  }
  
  # intensity weight (simple)
  out$w_int <- 1 / hr
  
  # combined weight
  out$w_all <- out$w_trt * out$w_int
  
  # trim
  out$w_trt_trim <- trim_by_quantile(out$w_trt, q = trim_q)
  out$w_int_trim <- trim_by_quantile(out$w_int, q = trim_q)
  out$w_all_trim <- trim_by_quantile(out$w_all, q = trim_q)
  
  out
}


##############################
library(splines2)
fit_blip_spline <- function(dat, dat_raw, weights,
                            tau1, tau2,
                            X_vars = c("X1", "X2", "X3"),
                            K1 = 5, K2 = 5,
                            n_y1 = 5, n_y2 = 5,   # kept for backward compatibility
                            degree = 3) {
  
  ## Use post-baseline visits only (here visit >= 1 is all rows)
  dat_use <- dat
  
  if (length(weights) != nrow(dat_use)) {
    stop("Length of 'weights' must match nrow(dat_use).")
  }
  dat_use$omega <- weights
  
  ## 1. Construct (y1, y2) grid using tau1_vec, tau2_vec from DGP ----
  ##    For simulate_dataset_5pts_blip, these are returned in the list.
  if (is.list(dat_raw) && all(c("tau1_vec", "tau2_vec") %in% names(dat_raw))) {
    tau1_vec <- dat_raw$tau1_vec
    tau2_vec <- dat_raw$tau2_vec
  } else {
    tau1_vec <- attr(dat_raw, "tau1_vec")
    tau2_vec <- attr(dat_raw, "tau2_vec")
  }
  
  if (is.null(tau1_vec) || is.null(tau2_vec)) {
    stop("tau1_vec or tau2_vec not found. Pass the full object from simulate_dataset_5pts_blip, or a data frame with these attributes.")
  }
  
  y1_grid <- sort(unique(as.numeric(tau1_vec)))
  y2_grid <- sort(unique(as.numeric(tau2_vec)))
  G1 <- length(y1_grid)
  G2 <- length(y2_grid)
  
  ## Helper: build 1D B-spline basis on a grid with df = K ----
  build_bs_basis <- function(y_grid, K, degree = 3) {
    rng <- range(y_grid)
    # number of interior knots
    n_int <- max(K - degree - 1, 0L)
    if (n_int > 0) {
      probs <- seq(0, 1, length.out = n_int + 2)[-c(1, n_int + 2)]
      knots <- as.numeric(quantile(y_grid, probs = probs))
    } else {
      knots <- NULL
    }
    splines2::bSpline(
      x = y_grid,
      knots = knots,
      degree = degree,
      Boundary.knots = rng,
      intercept = TRUE
    )
  }
  
  ## 2. Build 1D B-splines on y1_grid and y2_grid, then form tensor products ----
  B1 <- build_bs_basis(y1_grid, K = K1, degree = degree)  # G1 x K1_eff
  B2 <- build_bs_basis(y2_grid, K = K2, degree = degree)  # G2 x K2_eff
  
  K1_eff <- ncol(B1)
  K2_eff <- ncol(B2)
  Ktot   <- K1_eff * K2_eff
  
  # Tensor-product basis Φ(y1_g1, y2_g2) flattened to length Ktot
  G <- G1 * G2
  Phi <- matrix(NA_real_, nrow = G, ncol = Ktot)
  grid_pairs <- matrix(NA_integer_, nrow = G, ncol = 2)  # (g1, g2)
  
  idx <- 1L
  for (g1 in 1:G1) {
    for (g2 in 1:G2) {
      phi_vec <- as.numeric(kronecker(B1[g1, ], B2[g2, ]))  # length Ktot
      Phi[idx, ]       <- phi_vec
      grid_pairs[idx,] <- c(g1, g2)
      idx <- idx + 1L
    }
  }
  
  ## Locate (tau1, tau2) on the (y1_grid, y2_grid) tensor grid ----
  g1_tau <- which.min(abs(y1_grid - tau1))
  g2_tau <- which.min(abs(y2_grid - tau2))
  idx_tau <- which(grid_pairs[,1] == g1_tau & grid_pairs[,2] == g2_tau)
  if (length(idx_tau) != 1L) {
    stop("Failed to uniquely locate (tau1, tau2) on the (y1_grid, y2_grid) grid.")
  }
  # Tensor-product B-spline at (tau1, tau2), length Ktot
  phi_tau <- Phi[idx_tau, ]
  
  ## 3. Build stacked design matrix over (i, visit, y-grid) --------------------
  N_obs <- nrow(dat_use)                  # number of visit times
  X_mod <- as.matrix(dat_use[, X_vars, drop = FALSE])
  X_mod <- cbind(1, X_mod)                # add intercept
  colnames(X_mod)[1] <- "(Intercept)"
  q_x <- ncol(X_mod)                      # dimension of X_mod (with intercept)
  
  A_vec <- dat_use$A
  Y1    <- dat_use$Y1_obs
  Y2    <- dat_use$Y2_obs
  w_obs <- dat_use$omega
  
  N_total <- N_obs * G
  P       <- 2 * q_x * Ktot  # alpha and blip parts
  
  Z     <- matrix(0, nrow = N_total, ncol = P)
  U_vec <- integer(N_total)
  w_vec <- numeric(N_total)
  
  for (o in seq_len(N_obs)) {
    row_offset <- (o - 1L) * G
    x_mod_o    <- X_mod[o, ]
    A_o        <- A_vec[o]
    Y1_o       <- Y1[o]
    Y2_o       <- Y2[o]
    w_o        <- w_obs[o]
    
    for (g in 1:G) {
      r  <- row_offset + g
      g1 <- grid_pairs[g, 1]
      g2 <- grid_pairs[g, 2]
      
      y1_thr <- y1_grid[g1]
      y2_thr <- y2_grid[g2]
      
      # Joint binary outcome U_i(t; y1, y2)
      U_vec[r] <- as.integer(Y1_o <= y1_thr && Y2_o <= y2_thr)
      w_vec[r] <- w_o
      
      phi_g <- Phi[g, ]                # length Ktot
      # Baseline linear term: x_mod ⊗ phi_g (length q_x * Ktot)
      base_vec <- as.numeric(kronecker(phi_g, x_mod_o))
      
      # First half: alpha part; second half: blip part multiplied by A
      Z[r, 1:(q_x * Ktot)] <- base_vec
      Z[r, (q_x * Ktot + 1):P] <- A_o * base_vec
    }
  }
  
  ## 4. Weighted logistic regression to estimate {α_{kl}, β_{kl}} --------------
  fit <- glm.fit(x = Z, y = U_vec, family = binomial(), weights = w_vec)
  
  coef_all  <- coef(fit)
  alpha_vec <- coef_all[1:(q_x * Ktot)]
  beta_vec  <- coef_all[(q_x * Ktot + 1):P]
  
  # Reshape β into (q_x x Ktot) matrix for evaluation at (tau1, tau2)
  beta_mat <- matrix(beta_vec, nrow = q_x, ncol = Ktot, byrow = FALSE)
  
  # Blip coefficients at (tau1, tau2):
  #   β(tau1, tau2) = β_mat %*% φ(tau1, tau2)
  blip_coef_tau <- beta_mat %*% phi_tau  # length q_x, for (Intercept, X_vars)
  
  list(
    fit            = fit,
    alpha_vec      = alpha_vec,
    beta_vec       = beta_vec,
    beta_mat       = beta_mat,
    X_coef_names   = colnames(X_mod),   # (Intercept, X1, X2, X3, ...)
    y1_grid        = y1_grid,
    y2_grid        = y2_grid,
    Phi_grid       = Phi,
    phi_tau        = phi_tau,
    blip_coef_tau  = as.numeric(blip_coef_tau)  # β(tau1, tau2) w.r.t X_mod
  )
}




solve_ridge <- function(H, b,
                        lambda0 = 1e-8, mult = 10, max_try = 8, tol_qr = 1e-12) {
  p <- nrow(H)
  Hs <- 0.5 * (H + t(H))
  dbar <- mean(diag(Hs), na.rm = TRUE)
  if (!is.finite(dbar) || dbar <= 0) dbar <- 1
  lam <- lambda0 * dbar
  
  for (k in 0:max_try) {
    Hlam <- Hs + diag(lam, p)
    sol <- tryCatch(solve(Hlam, b), error = function(e) NULL)
    if (!is.null(sol) && all(is.finite(sol))) return(sol)
    lam <- lam * mult
  }
  qr.solve(Hs, b, tol = tol_qr)
}

fit_blip_spline_fast <- function(dat, dat_raw, weights,
                                 tau1, tau2,
                                 X_vars = c("X1","X2","X3"),
                                 K1 = 5, K2 = 5,
                                 n_y1 = 5, n_y2 = 5,
                                 degree = 3,
                                 maxit = 25, tol = 1e-6, verbose = FALSE) {
  
  ## Use post-baseline visits only
  dat_use <- dat
  
  ## Grid
  if (is.list(dat_raw) && all(c("tau1_vec","tau2_vec") %in% names(dat_raw))) {
    tau1_vec <- dat_raw$tau1_vec
    tau2_vec <- dat_raw$tau2_vec
  } else {
    tau1_vec <- attr(dat_raw, "tau1_vec")
    tau2_vec <- attr(dat_raw, "tau2_vec")
  }
  if (is.null(tau1_vec) || is.null(tau2_vec)) stop("tau1_vec or tau2_vec not found.")
  
  y1_grid <- sort(unique(as.numeric(tau1_vec)))
  y2_grid <- sort(unique(as.numeric(tau2_vec)))
  G1 <- length(y1_grid); G2 <- length(y2_grid)
  G  <- G1 * G2
  
  ## weights must be N x G: each column corresponds to a grid point g
  N <- nrow(dat_use)
  if (!is.matrix(weights)) {
    stop("'weights' must be a matrix with nrow = nrow(dat_use) and ncol = G (one column per grid point).")
  }
  if (nrow(weights) != N || ncol(weights) != G) {
    stop(sprintf("Dim(weights) must be %d x %d (N x G). Got %d x %d.",
                 N, G, nrow(weights), ncol(weights)))
  }
  omega <- weights
  if (any(!is.finite(omega))) stop("weights contains NA/Inf.")
  
  build_bs_basis <- function(y_grid, K, degree = 3) {
    rng <- range(y_grid)
    n_int <- max(K - degree - 1, 0L)
    if (n_int > 0) {
      probs <- seq(0, 1, length.out = n_int + 2)[-c(1, n_int + 2)]
      knots <- as.numeric(quantile(y_grid, probs = probs))
    } else {
      knots <- NULL
    }
    splines2::bSpline(
      x = y_grid, knots = knots, degree = degree,
      Boundary.knots = rng, intercept = TRUE
    )
  }
  
  B1 <- build_bs_basis(y1_grid, K = K1, degree = degree)
  B2 <- build_bs_basis(y2_grid, K = K2, degree = degree)
  K1_eff <- ncol(B1); K2_eff <- ncol(B2); Ktot <- K1_eff * K2_eff
  
  ## Phi grid (G x Ktot), order: g1 major, g2 minor
  Phi <- matrix(NA_real_, nrow = G, ncol = Ktot)
  idx <- 1L
  for (g1 in 1:G1) {
    for (g2 in 1:G2) {
      Phi[idx, ] <- as.numeric(kronecker(B1[g1, ], B2[g2, ]))
      idx <- idx + 1L
    }
  }
  
  ## Locate (tau1, tau2)
  g1_tau <- which.min(abs(y1_grid - tau1))
  g2_tau <- which.min(abs(y2_grid - tau2))
  idx_tau <- (g1_tau - 1L) * G2 + g2_tau
  phi_tau <- Phi[idx_tau, ]
  
  ## Data matrices
  X_mod <- as.matrix(dat_use[, X_vars, drop = FALSE])
  X_mod <- cbind(1, X_mod)
  colnames(X_mod)[1] <- "(Intercept)"
  q_x <- ncol(X_mod)
  
  A_vec <- as.numeric(dat_use$A)
  A_vec <- 1.0 * (A_vec != 0)
  
  Y1 <- as.numeric(dat_use$Y1_obs)
  Y2 <- as.numeric(dat_use$Y2_obs)
  
  if (any(!is.finite(Y1)) || any(!is.finite(Y2))) stop("Y1_obs/Y2_obs contains NA/Inf.")
  if (any(!is.finite(X_mod))) stop("X_vars contains NA/Inf.")
  
  qK <- q_x * Ktot
  P  <- 2 * qK
  
  ## Thresholds in Phi-row order: g1 major, g2 minor
  y1_thr <- rep(y1_grid, each = G2)
  y2_thr <- rep(y2_grid, times = G1)
  
  ## U: N x G (0/1)
  U <- (outer(Y1, y1_thr, "<=") & outer(Y2, y2_thr, "<=")) * 1.0
  
  ## Precompute Phi outer products (Ktot x Ktot) for each g
  Phi_outer_list <- lapply(1:G, function(g) tcrossprod(Phi[g, ]))
  
  ## Init theta
  theta <- rep(0, P)
  
  for (it in 1:maxit) {
    alpha_mat <- matrix(theta[1:qK], nrow = q_x, ncol = Ktot, byrow = FALSE)
    beta_mat  <- matrix(theta[(qK + 1):P], nrow = q_x, ncol = Ktot, byrow = FALSE)
    
    ## Coefs on grid: q_x x G
    AlphaCoef <- alpha_mat %*% t(Phi)
    BetaCoef  <- beta_mat  %*% t(Phi)
    
    eta <- (X_mod %*% AlphaCoef) + (X_mod %*% BetaCoef) * A_vec
    mu  <- plogis(eta)
    
    ## Score: use omega[, g] to weight column g
    R <- (U - mu) * omega
    S_alpha <- crossprod(X_mod, R) %*% Phi
    S_beta  <- crossprod(X_mod, R * A_vec) %*% Phi
    score <- c(as.vector(S_alpha), as.vector(S_beta))
    
    ## Fisher info blocks
    Haa <- matrix(0, nrow = qK, ncol = qK)
    Hab <- matrix(0, nrow = qK, ncol = qK)
    Hbb <- matrix(0, nrow = qK, ncol = qK)
    
    ## Wmat: N x G, column-g weights are omega[, g]
    Wmat <- mu * (1 - mu) * omega
    
    for (g in 1:G) {
      wg <- Wmat[, g]
      
      XWX <- crossprod(X_mod, X_mod * wg)
      
      wga   <- wg * A_vec
      XWX_a <- crossprod(X_mod, X_mod * wga)
      
      Pout <- Phi_outer_list[[g]]
      
      Haa <- Haa + kronecker(Pout, XWX)
      Hab <- Hab + kronecker(Pout, XWX_a)
      Hbb <- Hbb + kronecker(Pout, XWX_a)  # A^2 = A for 0/1
    }
    
    H <- rbind(cbind(Haa, Hab),
               cbind(t(Hab), Hbb))
    H <- 0.5 * (H + t(H))
    
    step <- solve_ridge(H, score)
    if (any(!is.finite(step))) stop("Newton step has NA/Inf even after ridge.")
    
    theta_new <- theta + step
    if (verbose) cat(sprintf("iter %d: max|step|=%.3e\n", it, max(abs(step))))
    theta <- theta_new
    if (max(abs(step)) < tol) break
  }
  
  ## Final unpack
  alpha_vec <- theta[1:qK]
  beta_vec  <- theta[(qK+1):P]
  
  alpha_mat <- matrix(alpha_vec, nrow = q_x, ncol = Ktot, byrow = FALSE)
  beta_mat  <- matrix(beta_vec,  nrow = q_x, ncol = Ktot, byrow = FALSE)
  
  ## Coefs on grid: q_x x G
  AlphaCoef <- alpha_mat %*% t(Phi)
  BetaCoef  <- beta_mat  %*% t(Phi)
  
  ## Grid-level quantities (N x G)
  kappa_grid <- X_mod %*% BetaCoef                         # x^T beta_g
  linear_predict_grid <- (X_mod %*% AlphaCoef) + kappa_grid * A_vec
  
  # Prob and kappa on the same N_obs x G layout
  prob_grid  <- 1 / (1 + exp(-linear_predict_grid))
  kappa_grid <- prob_grid * (1 - prob_grid)
  
  ## Tau-specific (vector length N)
  kappa_tau   <- kappa_grid[, idx_tau]
  prob_tau  <- prob_grid[, idx_tau]
  
  ## Coef at tau (q_x)
  blip_coef_tau <- beta_mat %*% phi_tau
  
  list(
    fit                 = list(converged = TRUE, coef = theta),
    alpha_vec           = alpha_vec,
    beta_vec            = beta_vec,
    beta_mat            = beta_mat,
    X_coef_names        = colnames(X_mod),
    y1_grid             = y1_grid,
    y2_grid             = y2_grid,
    Phi_grid            = Phi,
    idx_tau             = idx_tau,
    phi_tau             = phi_tau,
    blip_coef_tau       = as.numeric(blip_coef_tau),
    
    ## NEW outputs you asked for
    linear_predict_grid = linear_predict_grid,  # N x G
    prob_grid           = prob_grid,            # N x G
    kappa_grid          = kappa_grid,           # N x G
    kappa_tau           = as.numeric(kappa_tau),# N
    prob_tau            = as.numeric(prob_tau)  # N (optional)
  )
  
}





## Helper: extract blip coefficients at a given (tau1_new, tau2_new)
predict_blip_on_grid <- function(fit_obj, tau1_new, tau2_new) {
  y1_grid <- fit_obj$y1_grid
  y2_grid <- fit_obj$y2_grid
  G2      <- length(y2_grid)
  
  # Find the nearest grid indices g1, g2
  g1 <- which.min(abs(y1_grid - tau1_new))
  g2 <- which.min(abs(y2_grid - tau2_new))
  
  # Row index: idx = (g1 - 1)*G2 + g2
  idx <- (g1 - 1L) * G2 + g2
  
  # Corresponding tensor-product B-spline basis vector phi(tau1, tau2)
  phi_new <- fit_obj$Phi_grid[idx, ]  # length Ktot
  
  # Blip coefficients: beta(tau1, tau2) = beta_mat %*% phi(tau1, tau2)
  blip_coef <- as.numeric(fit_obj$beta_mat %*% phi_new)
  names(blip_coef) <- fit_obj$X_coef_names
  
  blip_coef
}

## Wrapper: compute blip parameters on the 5x5 DGP grid
get_blip_grid <- function(fit_obj) {
  tau1_vec <- fit_obj$y1_grid  # length 5
  tau2_vec <- fit_obj$y2_grid  # length 5
  
  grid_25 <- expand.grid(
    tau1 = tau1_vec,
    tau2 = tau2_vec
  )
  
  # Use lapply to ensure a list output (rather than a matrix)
  blip_list <- lapply(seq_len(nrow(grid_25)), function(i) {
    tau1_new <- grid_25$tau1[i]
    tau2_new <- grid_25$tau2[i]
    predict_blip_on_grid(fit_obj, tau1_new, tau2_new)
  })
  
  blip_mat <- do.call(rbind, blip_list)          # 25 x q_x
  blip_df  <- cbind(grid_25, as.data.frame(blip_mat))
  rownames(blip_df) <- NULL
  
  blip_df
}



## Predict optimal treatment using a given blip_coef (named vector)
predict_itr_blip <- function(blip_coef, newdata, X_vars,
                             rule = c("gt0","ge0")) {
  rule <- match.arg(rule)
  
  Xmat <- as.matrix(newdata[, X_vars, drop = FALSE])
  Xmat <- cbind(1, Xmat)
  colnames(Xmat)[1] <- "(Intercept)"
  
  # align columns with blip_coef names
  Xmat <- Xmat[, names(blip_coef), drop = FALSE]
  
  score <- as.numeric(Xmat %*% blip_coef)
  if (rule == "gt0") as.integer(score > 0) else as.integer(score >= 0)
}




