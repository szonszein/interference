make_exposure_prob_individual <-
  function(prob_exposure,
           i_start = NULL,
           i_end = NULL) {
    if (is.null(i_start)) {
      i_start <- 1
    }
    if (is.null(i_end)) {
      i_end <- dim(prob_exposure$I_exposure[[1]])[1]
    }
    
    R <- dim(prob_exposure$I_exposure[[1]])[2]
    
    
    #  For each of the exposure conditions in prob_exposure$I_exposure
    prob_exposure_cond <-
      do.call(rbind, lapply(
        prob_exposure$I_exposure,
        #  Use drop=FALSE to keep a matrix (not a vector) so that 
        #  rowSums works if i_start==i_end
        FUN = function(x)
          (rowSums(x[i_start:i_end, , drop=FALSE]) + 1) / (R + 1)
      ))
    
    return(prob_exposure_cond)
  }


make_kl_names <- function(exposure_names) {
  kl_names <- c()
  for (k in exposure_names) {
    for (l in exposure_names) {
      if (k!=l) {
        kl_names <-c(kl_names, paste(k, l, sep=','))
      }
    }
  }
  return(list(kl_names))
}

estimates_batch <- function(obs_exposure,
                            obs_outcome,
                            n_var_permutations = 10,
                            
                            effect_estimators = c('hajek', 'horvitz-thompson'),
                            variance_estimators = c('hajek', 'horvitz-thompson'),
                            
                            control_condition=NULL,
                            treated_conditions=NULL,
                            
                            potential_tr_vector,
                            adj_matrix,
                            exposure_map_fn,
                            exposure_map_fn_add_args = NULL,
                            
    
                            batch_size = NA) {
  if (!is.null(control_condition) & is.null(treated_conditions)) {
    treated_conditions <- setdiff(colnames(obs_exposure), control_condition)
  }
  k_to_include <- treated_conditions
  l_to_include <- control_condition
  
  
  if (('constant_effect' %in% variance_estimators) | ('max_ht_const' %in% variance_estimators)) {
    stop("constant_effect and max_ht_const variance estimators are not supported with batch estimation")
  }
  if (('hajek' %in% variance_estimators) & !('hajek' %in% effect_estimators)) {
    effect_estimators  <- c(effect_estimators, 'hajek')
  }
  
  out <- list()
  
  N <- nrow(obs_exposure)
  if (is.na(batch_size)) {
    batch_size <- N
  }
  
  obs_outcome_by_exposure <- t(obs_exposure) %*% diag(obs_outcome)
  obs_outcome_by_exposure[t(obs_exposure) == 0] <- NA
  
  # TODO: Figure out why k and l are switched relative to package

  running_yT_ht <- vector(mode = 'numeric', ncol(obs_exposure))
  running_mu_h_divisor <- vector(mode = 'numeric', ncol(obs_exposure))
  names(running_mu_h_divisor) <- colnames(obs_exposure)
  var_yT <-
    matrix(NA,
           nrow = ncol(obs_exposure),
           dimnames = list(colnames(obs_exposure)))
  #  Matrix of adjustments to the variance, initalized at zero because
  #  it is only ever summed with the variance; if the variance is undefined,
  #  the sum will be too

  var_yT_A2 <-
    matrix(0,
           nrow = ncol(obs_exposure),
           dimnames = list(colnames(obs_exposure)))
  cov_yT_A <-
    matrix(NA,
           nrow = 2 * choose(ncol(obs_exposure), 2),
           dimnames = make_kl_names(colnames(obs_exposure))
    )
  if (is.null(k_to_include)) {
    k_to_include <- colnames(obs_exposure)
  }
  if (is.null(l_to_include)) {
    l_to_include <- colnames(obs_exposure)
  }
  
  
  #  Calculate the point estimates, variances, and covariances for one batch at
  #  a time.
  for (i in seq(1, (N), by = batch_size)) {
    #  The batch goes from unit i to unit j, so loop over the starting `i`s
    #  and calculate the ending `j`s, which are always `batch_size` more
    #  except for the final one, which cannot be larger than N
    j <- i + batch_size - 1
    if (j > N) {
      j <- N
    }
  
    

    obs_prob_exposure <-
      make_exposure_prob(
        potential_tr_vector,
        adj_matrix,
        exposure_map_fn,
        exposure_map_fn_add_args = exposure_map_fn_add_args,
        i_start = i,
        i_end = j
      )
    
    # TODO: calculate obs_prob_exposure_individual_kk only once for all i, j, and subset here
    obs_prob_exposure_individual_kk <- make_exposure_prob_individual(obs_prob_exposure, i, j)
    
    sums <-
      rowSums(obs_outcome_by_exposure[, i:j] / obs_prob_exposure_individual_kk,
              na.rm = T)
    running_yT_ht <- running_yT_ht + sums
    
    
    running_mu_h_divisor <- running_mu_h_divisor + rowSums(t(obs_exposure[i:j,,drop=FALSE]) / obs_prob_exposure_individual_kk, na.rm = T)
    
    
    #  TODO: reuse the single obs_prob_exposure_individual_kk above for all_ind_kk
    all_ind_kk <- make_exposure_prob_individual(obs_prob_exposure)
    
    # k indexes the exposure conditions
    for (k in colnames(obs_exposure)) {
      ## Unadjusted variance
      pi_k <-
        obs_prob_exposure$prob_exposure_k_k[[paste(k, k, sep = ',')]]
      ind_kk <- all_ind_kk[k, ]
      
      cond_indicator <- obs_exposure[, k]
      
      ind_kk_sq <- (ind_kk %o% ind_kk)[i:j,]
      mm <-
        (cond_indicator[i:j] %o% cond_indicator) * (pi_k - ind_kk_sq) / pi_k * (obs_outcome[i:j] %o% obs_outcome) / (ind_kk_sq)
      mm[!is.finite(mm)] <- 0
      
      #  If there are no units in this condition, the variance remains NA
      if (is.na(var_yT[k, ] )& any(cond_indicator[i:j] == 1)) {var_yT[k, ] <-0}
      var_yT[k, ] <- var_yT[k, ] + sum(mm)
      
      
      # variance adjustment  A2
      m <- cond_indicator*(obs_outcome^2)/(2*ind_kk)
      
      A2_part_sum <-
        sum(outer(m[i:j], m, FUN = '+') * (pi_k == 0) * (!diag(length(m)))[i:j,] )
      
    
      var_yT_A2[k, ] <- var_yT_A2[k, ] + A2_part_sum
      
      
      # covariance
      for (l in l_to_include) {
        if ((k %in% k_to_include) & (k != l)) {
          kl <- paste(k, l, sep = ',')
          pi_l <-
            obs_prob_exposure$prob_exposure_k_k[[paste(l, l, sep = ',')]]
          ind_ll <- all_ind_kk[l, ]
          pi_k_l <- obs_prob_exposure$prob_exposure_k_l[[kl]]
          cond_indicator_l <- obs_exposure[, l]
          
          ind_kk_ll <- (ind_kk %o% ind_ll)[i:j,]
          
          mm <-
            cond_indicator[i:j] %o% cond_indicator_l * (pi_k_l - ind_kk_ll) / pi_k_l * (obs_outcome[i:j] %o% obs_outcome) / (ind_kk_ll)
          mm[!is.finite(mm)] <- 0
          
          if (is.na(cov_yT_A[kl,]) ) {
            cov_yT_A[kl,] <- 0
          }
          cov_yT_A[kl,] <- cov_yT_A[kl,] + sum(mm)

          
          # TODO: Batch this part
          second_part_cov <- 0
          # Using x and y instead of i and j because those are already being used for the batching
          for (x in 1:nrow(pi_k_l)) {
            for (y in 1:ncol(pi_k_l)) {
              if (pi_k_l[x, y] == 0) {

                second_part_cov_x_y <-
                  (
                    (cond_indicator[i+x-1] * obs_outcome[i+x-1] ^ 2 / (2 * ind_kk[i+x-1])) + 
                    (cond_indicator_l[y] * obs_outcome[y] ^ 2 / (2 * ind_ll[y]))
                  )
                
                second_part_cov <-
                  second_part_cov + second_part_cov_x_y
                
              }
            }
            
          }

          cov_yT_A[kl,] <- cov_yT_A[kl,] - second_part_cov
          
        }
      }
      
      
    }
    

    
  }
  
  cov_yT_ht <- cov_yT_A
  var_yT_ht <- var_yT + var_yT_A2
  
  yT_ht <- running_yT_ht
  yT_ht[rowSums(!is.na(obs_outcome_by_exposure)) == 0] <- NA
  
  if ('horvitz-thompson' %in% effect_estimators) {
    out[['yT_ht']] <- yT_ht
  }
  

  
  if (('hajek' %in% effect_estimators)) {
    
    mu_h <-
      yT_ht / running_mu_h_divisor
    mu_h[rowSums(!is.na(obs_outcome_by_exposure)) == 0] <- NA
    
    yT_h <- mu_h * N
    yT_h[rowSums(!is.na(obs_outcome_by_exposure)) == 0] <- NA
    
    out[['yT_h']] <- yT_h
  }
  
  
  #  The variance of the hajek estimator is calculated in the same way as the
  #  variance of the HT estimator, but substituting the residuals (`resid_h`)
  #  for the observed outcomes. This means that we first need to calculate
  #  the yT_ht and from that the residuals, and then call estimates_batch
  #  again with the residuals in place of the outcomes.
  
  if (('hajek' %in% variance_estimators)) {


    resid_h <-
      colSums((obs_outcome_by_exposure - mu_h), na.rm = T) # pass resid_h instead of obs_outcome to var, cov

    hajek_estimates <- estimates_batch(obs_exposure=obs_exposure,
                                obs_outcome=resid_h,
                                n_var_permutations = n_var_permutations,
                                
                                effect_estimators = c('horvitz-thompson'),
                                variance_estimators = c( 'horvitz-thompson'),
                                
                                control_condition=control_condition,
                                treated_conditions=treated_conditions,
                                
                                potential_tr_vector=potential_tr_vector,
                                adj_matrix=adj_matrix,
                                exposure_map_fn=exposure_map_fn,
                                exposure_map_fn_add_args = exposure_map_fn_add_args,
                                
                
                                batch_size = batch_size)
    
    var_yT_h <-  hajek_estimates[['var_yT_ht']]
    cov_yT_h <- hajek_estimates[['cov_yT_ht']]
    out[['var_yT_h']] <-var_yT_h
    out[['cov_yT_h']] <- cov_yT_h
    var_tau_h <-
      (1 / N ^ 2) * (var_yT_h[!rownames(var_yT_h) %in% control_condition,] + var_yT_h[rownames(var_yT_h) %in% control_condition,] -
                       2 * cov_yT_h[rownames(cov_yT_h) %in% paste(treated_conditions, control_condition, sep=','),])
    
    out[['var_tau_h']] <- var_tau_h
  }
  
  

  
  
  if ('horvitz-thompson' %in% variance_estimators) {
    out[['var_yT_ht']] <- var_yT_ht
    out[['cov_yT_ht']] <- cov_yT_ht
    
    var_tau_ht <-
      (1 / N ^ 2) * (var_yT_ht[!rownames(var_yT_ht) %in% control_condition,] + var_yT_ht[rownames(var_yT_ht) %in% control_condition,] -
                       2 * cov_yT_ht[rownames(cov_yT_ht) %in% paste(treated_conditions, control_condition, sep=','),])
    out[['var_tau_ht']] <- var_tau_ht
    
  }
  
  
  if (('horvitz-thompson' %in% effect_estimators)) {
    tau_ht <- (1 / N) * (yT_ht - yT_ht['no'])[names(yT_ht) != 'no']
    out[['tau_ht']] <- tau_ht
  }
  if (('hajek' %in% effect_estimators)) {
    tau_h <- (mu_h - mu_h['no'])[names(mu_h) != 'no']
    out[['tau_h']] <- tau_h
  }
  
  return(out)
  
}