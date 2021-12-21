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
    
    prob_exposure_cond <-
      do.call(rbind, lapply(
        prob_exposure$I_exposure,
        FUN = function(x)
          (rowSums(x[i_start:i_end, ]) + 1) / (R + 1)
      ))
    
    return(prob_exposure_cond)
  }


estimates_batch <- function(obs_exposure,
                            obs_outcome,
                            obs_prob_exposure,
                            potential_tr_vector,
                            adj_matrix,
                            exposure_map_fn,
                            exposure_map_fn_add_args = NULL,
                            n_var_permutations = 10,
                            batch_size = NA,
                            hop) {
  out <- list()
  
  N <- nrow(obs_exposure)
  if (is.na(batch_size)) {
    batch_size <- N
  }
  
  obs_outcome_by_exposure <- t(obs_exposure) %*% diag(obs_outcome)
  obs_outcome_by_exposure[t(obs_exposure) == 0] <- NA
  
  # TODO: Figure out why k and l are switched relative to package
  # TODO :remove requirement to specify k_to_include and l_to_include
  k_to_include <- c('dir_ind1', 'isol_dir', 'ind1')
  l_to_include <- 'no'
  
  running_yT_ht <- vector(mode = 'numeric', ncol(obs_exposure))
  var_yT <-
    matrix(NA,
           nrow = ncol(obs_exposure),
           dimnames = list(colnames(obs_exposure)))
  var_yT_A2 <-
    matrix(NA,
           nrow = ncol(obs_exposure),
           dimnames = list(colnames(obs_exposure)))
  cov_yT_A <-
    matrix(NA,
           nrow = 2 * choose(ncol(obs_exposure), 2),
           dimnames = list(names(obs_prob_exposure$prob_exposure_k_l)))
  if (is.null(k_to_include)) {
    k_to_include <- colnames(obs_exposure)
  }
  if (is.null(l_to_include)) {
    l_to_include <- colnames(obs_exposure)
  }
  for (i in seq(1, (N), by = batch_size)) {
    j <- i + batch_size - 1
    if (j > N) {
      j <- N
    }
    #TODO: deal with non-divisible batch sizes
    
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
    obs_prob_exposure_individual_kk <-
      make_exposure_prob_individual(obs_prob_exposure, i, j)
    
    sums <-
      rowSums(obs_outcome_by_exposure[, i:j] / obs_prob_exposure_individual_kk,
              na.rm = T)
    
    running_yT_ht <- running_yT_ht + sums
    
    
    
    #  TODO: reuse the single obs_prob_exposure_individual_kk above for all_ind_kk
    all_ind_kk <- make_exposure_prob_individual(obs_prob_exposure)
    
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
      m <- cond_indicator[i:j] * (obs_outcome[i:j] ^ 2) / (2 * ind_kk[i:j])
      A2_part_sum <-
        sum(outer(m, m, FUN = '+') * (pi_k[i:j, i:j] == 0) * (!diag(length(m))))
      
      #  If there are no units in this condition, the variance remains NA
      if (is.na(var_yT_A2[k, ]) & any(cond_indicator[i:j] ==1)) {var_yT_A2[k, ] <-0}
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
          
          if (is.na(cov_yT_A[kl,])) {
            cov_yT_A[kl,] <- 0
          }
          cov_yT_A[kl,] <- cov_yT_A[kl,] + sum(mm)
          
          # TODO: Batch this part
          second_part_cov <- 0
          for (x in 1:nrow(pi_k_l)) {
            for (y in 1:ncol(pi_k_l)) {
              if (pi_k_l[x, y] == 0) {
                second_part_cov_x_y <-
                  ((cond_indicator[x] * obs_outcome[x] ^ 2 / (2 * ind_kk[x])) + (
                    cond_indicator_l[y] * obs_outcome[y] ^ 2 / (2 * ind_ll[y])
                  )
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
  var_yT_ht_adjusted <- var_yT + var_yT_A2
  
  yT_ht <- running_yT_ht
  yT_ht[rowSums(!is.na(obs_outcome_by_exposure)) == 0] <- NA
  
  return(list(
    yT_ht = yT_ht,
    var_yT_ht = var_yT_ht_adjusted,
    cov_yT_ht = cov_yT_A
  ))
  
}