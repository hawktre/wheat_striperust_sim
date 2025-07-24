## ---------------------------
##
## Script name: 04a_SimFunc.R
##
## Purpose of script: Compile all necessary functions for simulation.
##
## Author: Trent VanHawkins
##
## Date Created: 2025-07-22
##
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

## ---------------------------

## view outputs in non-scientific notation

options(scipen = 6, digits = 4) 

## ---------------------------

## load up the packages we will need:  (uncomment as required)

# Function to simulate the data -------------------------------------------
disease_sim <- function(pars, mu){
  phi <- pars[["phi"]]
  
  #convert back to shape and scale 
  alpha <- mu * phi
  beta <- phi*(1-mu)
  
  #Draw from beta distribution
  map2_dbl(alpha, beta, ~rbeta(1, shape1 = .x, shape2 = .y)) 
}

# Forward Fit -------------------------------------------------------------
forward_fit <- function(blk, trt, vst, mod_dat, dist, kappa_try) {
  ## Extract Needed Model Data
  intensity <- mod_dat$intensity[, blk, trt, vst]
  intensity_prev <- mod_dat$intensity[, blk, trt, as.numeric(vst) - 1]
  wind <- mod_dat$wind[, , blk, trt, vst]
  pi <- mean(intensity == 0)
  non_zero <- which(intensity != 0)
  
  ## Set up strage for temproary results
  results_list <- list()
  
  for (kappa in kappa_try) {
    init_theta <- mod_dat$inits[, kappa, blk, trt, vst] #Current inits
    
    fit <- try(
      optim(
        par = init_theta,
        fn = neg_loglik,
        gr = neg_grad,
        method = "BFGS",
        hessian = TRUE,
        control = list(maxit = 5000, reltol = 1e-8),
        y_current = intensity[non_zero],
        y_prev = intensity_prev[non_zero],
        wind_matrix = wind[non_zero, non_zero],
        dist_matrix = dist[non_zero, non_zero]
      ),
      silent = TRUE
    )
    
    if (!inherits(fit, "try-error")) {
      se <- sqrt(diag(ginv(fit$hessian)))
      results_list[[length(results_list) + 1]] <- list(
        block = blk,
        treat = as.numeric(trt),
        visit = as.numeric(vst),
        iters = fit$counts[1],
        init_kappa = unname(init_theta["kappa"]),
        neg_loglik = fit$value,
        converged = fit$convergence == 0,
        theta = list(fit$par),
        theta_se = list(se),
        pi = pi
      )
    }
  }
  
  if (length(results_list) > 0) {
    visit_dt <- rbindlist(results_list)
    visit_dt <- visit_dt[neg_loglik > -1e5]
    if (nrow(visit_dt) > 0) {
      return(visit_dt[which.min(neg_loglik)])
    }
  }
  
  return(NULL)
}


# Backward Fit ------------------------------------------------------------
backward_fit <- function(config, blk, trt, vst, mod_dat, forward_fits) {
  # Setup
  intensity <- mod_dat$intensity[, blk, trt, vst]
  intensity_prev <- mod_dat$intensity[, blk, trt, as.numeric(vst) - 1]
  wind <- mod_dat$wind[, , blk, trt, vst]
  dist <- mod_dat$dist
  groups <- mod_dat$groups[, config]
  non_zero <- which(intensity != 0)
  n_groups <- readr::parse_number(config)
  prior <- 1 / n_groups
  
  # Starting values from forward fit
  theta <- forward_fits[
    block == blk & treat == as.numeric(trt) & visit == as.numeric(vst)
  ][["theta"]][[1]]
  
  max_em_iter <- 1000
  tol <- 1e-8
  q_track <- numeric(max_em_iter)
  
  # Initial E-step
  p_mat <- e_step(
    par = theta,
    y_current = intensity,
    y_prev = intensity_prev,
    wind_matrix = wind,
    dist_matrix = dist,
    group_id = groups,
    prior = prior
  )
  
  # EM loop
  for (em_iter in seq_len(max_em_iter)) {
    fit <- try(
      optim(
        par = theta,
        fn = wrapped_obj,
        gr = mstep_grad_em,
        method = "BFGS",
        control = list(maxit = 1000, reltol = 1e-8),
        p_mat = p_mat[non_zero, ],
        y_current = intensity[non_zero],
        y_prev = intensity_prev[non_zero],
        wind_matrix = wind[non_zero, non_zero],
        dist_matrix = dist[non_zero, non_zero],
        group_id = groups[non_zero]
      ),
      silent = TRUE
    )
    
    if (inherits(fit, "try-error")) {
      warning(paste("EM step", em_iter, "failed"))
      return(NULL)
    }
    
    theta_new <- fit$par
    
    # Update E-step
    p_mat <- e_step(
      par = theta_new,
      y_current = intensity,
      y_prev = intensity_prev,
      wind_matrix = wind,
      dist_matrix = dist,
      group_id = groups,
      prior = prior
    )
    
    # Compute Q
    mu_mat <- get_mu(
      par = theta_new,
      y_prev = intensity_prev,
      wind_matrix = wind,
      dist_matrix = dist,
      group_id = groups
    )
    
    q_val <- Q_fun(
      y = intensity,
      mu_mat = mu_mat,
      phi = theta_new[["phi"]],
      p_mat = p_mat
    )
    
    q_track[em_iter] <- q_val
    
    if (em_iter > 1 && abs(q_track[em_iter] - q_track[em_iter - 1]) < tol) {
      break
    }
    
    theta <- theta_new
  }
  
  # Return a single-row data.table
  result <- data.table(
    config = config,
    block = blk,
    treat = as.numeric(trt),
    visit = as.numeric(vst),
    em_iters = em_iter,
    converged = em_iter < max_em_iter,
    neg_loglik = q_track[em_iter],
    final_theta = list(theta_new),
    p_mat = list(p_mat)
  )
  
  return(result)
}

# Backward_Preds ----------------------------------------------------------
source_pred <- function(config, blk, trt, vst, p_mat, mod_dat) {
  
  # Ensure groups is a factor with levels = 1:n_groups
  groups <- as.factor(mod_dat$groups[,config])
  levels(groups) <- sort(unique(groups))
  
  # 1. Compute p̄ matrix (average p_mat rows by destination group)
  group_ids <- sort(unique(groups))
  n_groups <- length(group_ids)
  
  p_bar <- matrix(NA, n_groups, n_groups)
  rownames(p_bar) <- colnames(p_bar) <- as.character(group_ids)
  
  for (g in group_ids) {
    rows_in_group <- which(groups == g)
    p_bar[as.character(g), ] <- colMeans(p_mat[rows_in_group, , drop = FALSE])
  }
  
  # 2. Predict source group for each destination group (row-wise argmax)
  predicted_source <- which(p_bar == max(p_bar), arr.ind = TRUE)[[2]]
  
  # 3. Compare to ground truth
  true_source <- mod_dat$truth[blk, config]  # vector of length n_groups
  
  # 4. Compute distance-weighted accuracy
  dist <- mod_dat$grid_dist[[config]]
  error <- dist[predicted_source, true_source]
  dist_acc <- 1 - (error / max(dist[predicted_source,]))  # normalized distance accuracy
  
  result <- list(
    config = config,
    block = blk,
    treat = trt,
    visit = vst,
    predicted_source = predicted_source,
    true_source = true_source,
    dist_acc = dist_acc
  )
  # 5. Compute other metrics 
  return(as.data.table(result))
}


# Wrapper Function for the simulation -------------------------------------
single_sim <- function(sim_id, dat, forward_mod, output_dir = here("DataProcessed/results/simulation")) {
  tryCatch({
    browser()
    # Set up indices
    blocks <- dimnames(dat$intensity)[["block"]]
    treats <- dimnames(dat$intensity)[["treat"]]
    visits <- dimnames(dat$intensity)[["visit"]]
    kappa_try <- dimnames(dat$inits)[["kappa_try"]]
    configs <- dimnames(dat$groups)[["config"]]
    
    # 1. Simulate new intensity values
    intensity_sim <- dat$intensity
    for (blk in blocks) {
      for (trt in treats) {
        for (vst in visits[-1]) {
          fit <- forward_mod[block == blk & treat == as.numeric(trt) & visit == as.numeric(vst)]
          pars <- fit[["theta"]][[1]]
          fitted <- fit[["fitted"]][[1]]
          intensity_sim[, blk, trt, vst] <- disease_sim(pars, fitted)
        }
      }
    }
    dat$intensity <- intensity_sim
    
    # 2. Fit forward model
    combos_forward <- expand.grid(block = blocks, treat = treats, visit = visits[-1], stringsAsFactors = FALSE)
    forward <- pmap(combos_forward, ~forward_fit(..1, ..2, ..3, dat, dat$dist, kappa_try)) %>% rbindlist()
    
    # 3. Fit backward model
    combos_backward <- expand.grid(config = configs, bk = blocks, trt = treats, vst = visits[-1], stringsAsFactors = FALSE)
    backward <- pmap(combos_backward, ~backward_fit(..1, ..2, ..3, ..4, dat, forward)) %>% rbindlist()
    
    # 4. Source prediction for treat == 1
    backward_t1 <- backward[treat == 1]
    sources_predicted <- pmap(backward_t1, ~source_pred(config = ..1,
                                                        blk = ..2,
                                                        trt = ..3,
                                                        vst = ..4,
                                                        p_mat = ..9,
                                                        mod_dat = dat)) %>% 
      rbindlist()
    
    results_merge <- left_join(backward, forward, by = c("block", "treat", "visit"), suffix = c(".backward", ".forward")) %>% 
      dplyr::select(-c(final_theta, p_mat)) %>% 
      left_join(sources_predicted, by = c("config", "block", "treat", "visit")) %>% 
      mutate(sim = sim_id) %>% 
      dplyr::select(sim, everything())
    # 5. Return results
    return(results_merge)
    
  }, error = function(e) {
    # On failure, save what we can and keep going
    err_file <- file.path(output_dir, paste0("sim_error_", sim_id, ".txt"))
    writeLines(c("Simulation failed:", conditionMessage(e)), con = err_file)
    return(NULL)
  })
}