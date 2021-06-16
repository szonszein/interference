#' Estimate exposure-specific causal effects.
#'
#' Estimate exposure-specific causal effects and their variance.
#'
#' \code{estimates} produces values for the estimator of the average unit-level
#' causal effect of exposure k versus l and its variance estimator for the
#' exposure mappings returned by function \code{\link{make_exposure_map_AS}},
#' using a Horvitz-Thompson and a Hajek estimator. It also computes
#' Horvitz-Thompson and Hajek estimators of the total of potential outcomes,
#' which are inputs in the computation of average unit-level causal effect of
#' exposure k versus l.
#'
#' @param obs_exposure an `N`\eqn{*} K named numeric matrix of indicators for
#'   whether units `N` are in exposure condition k, where K is the total number
#'   of exposure conditions and names correspond to the exposure conditions.
#'   Such matrix is returned by function \code{\link{make_exposure_map_AS}}.
#' @param obs_outcome a vector length `N` of outcome data.
#' @param obs_prob_exposure a list of 3 lists containing exposure probabilities:
#'   \describe{ \item{`I_exposure`:}{A list of K `N` \eqn{*} `R` numeric
#'   matrices of indicators for whether units `N` are in exposure condition k
#'   over each of the possible `R` treatment assignment vectors. The number of
#'   numeric matrices K corresponds to the number of exposure conditions.}
#'   \item{`prob_exposure_k_k`:}{A list of K symmetric `N` \eqn{*} `N` numeric
#'   matrices each containing individual exposure probabilities to condition k
#'   on the diagonal, and joint exposure probabilities to condition k on the
#'   off-diagonals.} \item{`prob_exposure_k_l`:}{A list of
#'   \eqn{permutation(K,2)} nonsymmetric `N` \eqn{*} `N` numeric matrices each
#'   containing joint probabilities across exposure conditions k and l on the
#'   off-diagonal, and zeroes on the diagonal. When K = 4, the number of numeric
#'   matrices is 12; \eqn{permutation(4,2)}.} } Such list is returned by
#'   function \code{\link{make_exposure_prob}}.
#' @param n_var_permutations if `'constant_effect'` (derived in Aronow (2013)) 
#'   is one of the variance estimators specified, the number of treatment 
#'   permutations used to estimate this estimator. Default is
#'   `10`, but must be smaller or equal to `R`, the number of permutations to
#'   compute exposure probabilities. Recommended is `1000`, when `R` \eqn{>
#'   1000}.
#' @param effect_estimators string vector with names of estimators to be estimated
#' among 'hajek', 'horvitz-thompson'. Default is both.
#' @param variance_estimators string vector with names of variance estimators
#' to be estimated among 'hajek', 'horvitz-thompson', 'constant_effect',
#' 'max_ht_const'. Default includes the first two. Estimating 'constant_effect'
#' or 'max_ht_const' signficantly increases the running time.
#' @param hop number; either `1` or `2`. Must be `1` if argument `hop = 1` in
#'   function \code{\link{make_exposure_map_AS}} which assumes first-degree
#'   interference and produces four exposure conditions. Must be `2` if argument
#'   `hop = 2` in function \code{\link{make_exposure_map_AS}} which assumes
#'   second-degree interference and produces eight exposure conditions.
#' @export
#' @references Aronow, P. M. (2013). [Model assisted causal
#'   inference](https://search.proquest.com/docview/1567045106?accountid=12768).
#'   *PhD thesis, Department of Political Science, Yale University, New Haven,
#'   CT*.
#'
#'   Aronow, P.M. & Samii, C. (2017). [Estimating average causal effects under
#'   general interference, with application to a social network
#'   experiment](https://doi.org/10.1214/16-AOAS1005). *The Annals of Applied
#'   Statistics*, 11(4), 1912--1947.
#'
#'   Aronow, P.M. et al. (2020). [Spillover effects in experimental
#'   data](https://arxiv.org/abs/2001.05444). *arXiv preprint*,
#'   arXiv:2001.05444.
#' @examples
#' # Create adjacency matrix and treatment assignment vector
#' # to produce observed exposure conditions:
#'
#' adj_matrix <- make_adj_matrix(N = 9, model = 'sq_lattice')
#'
#' tr_vector <- make_tr_vec_permutation(N = 9, p = 0.2,
#'                                      R = 1, seed = 357)
#'
#' obs_exposure <- make_exposure_map_AS(adj_matrix, tr_vector,
#'                                      hop = 1)
#'
#' # Simulate a vector of outcome data:
#'
#' potential_outcome <- make_dilated_out(adj_matrix, make_corr_out,
#'                                       seed = 357, hop = 1)
#'
#' obs_outcome <- rowSums(obs_exposure*t(potential_outcome))
#'
#' # Create exposure probabilities:
#'
#' potential_tr_vector <- make_tr_vec_permutation(N = 9, p = 0.2,
#'                                                R = 36,
#'                                                seed = 357)
#'
#' obs_prob_exposure <- make_exposure_prob(potential_tr_vector,
#'                                         adj_matrix,
#'                                         make_exposure_map_AS,
#'                                         list(hop=1))
#'
#' # Estimate exposure-specific causal effects and their variance:
#'
#' estimates(obs_exposure, obs_outcome, obs_prob_exposure,
#'                                      n_var_permutations = 30,
#'                                      hop = 1)
#' @return A list of 13 lists: \describe{ \item{`yT_ht`:}{A named numeric vector
#'   which contains the values of the Horvitz-Thompson estimator of the total of
#'   potential outcomes under each exposure condition as derived in Equation 1
#'   of Aronow and Samii (2017).} \item{`yT_h`:}{A named numeric vector which
#'   contains the values of the Hajek estimator of the total of potential
#'   outcomes under each exposure condition as derived in Equation 15 of Aronow
#'   and Samii (2017).} \item{`var_yT_ht`:}{A named numeric K \eqn{*} 1 matrix
#'   which contains the values of the variance estimator of the Horvitz-Thompson
#'   estimator of the total of potential outcomes under each exposure condition
#'   as derived in Equation 7 and Proposition 5.1 of Aronow and Samii (2017).}
#'   \item{`var_yT_h`:}{A named numeric K \eqn{*} 1 matrix which contains the
#'   values of the variance estimator of the Hajek estimator of the total of
#'   potential outcomes under each exposure condition as explained in the first
#'   paragraph of page 1929 of Aronow and Samii (2017).} \item{`cov_yT_ht`:}{A
#'   named numeric \eqn{permutation(K,2)} \eqn{*} 1 matrix which contains the
#'   values of the covariance estimator of the Horvitz-Thompson estimator of the
#'   total of potential outcomes across exposures conditions k and l as derived
#'   in Equation 10 of Aronow and Samii (2017). When the number of exposure
#'   conditions K = 4, then the number of rows of this matrix is 12;
#'   \eqn{permutation(4,2)}.} \item{`cov_yT_h`:}{A named numeric
#'   \eqn{permutation(K,2)} \eqn{*} 1 matrix which contains the values of the
#'   covariance estimator of the Hajek estimator of the total of potential
#'   outcomes across exposures conditions k and l as explained in the first
#'   paragraph of page 1929 of Aronow and Samii (2017). When the number of
#'   exposure conditions K = 4, then the number of rows of this matrix is 12;
#'   \eqn{permutation(4,2)}.} \item{`tau_ht`:}{A named numeric vector which
#'   contains the  values of the Horvitz-Thompson estimator of the average
#'   unit-level causal effect of exposure k versus exposure l as derived in
#'   Equation 3 of Aronow and Samii (2017). Here exposure l is fixed to the
#'   \eqn{No Exposure} condition (i.e. no direct or indirect exposure).}
#'   \item{`tau_h`:}{A named numeric vector which contains the values of the
#'   Hajek estimator of the average unit-level causal effect of exposure k
#'   versus exposure l. Here exposure l is fixed to the \eqn{No Exposure}
#'   condition (i.e. no direct or indirect exposure).} \item{`tau_dsm`:}{A named
#'   numeric vector which contains the values of the difference in sample means
#'   estimator of the total observed outcomes across exposures k and l. Here
#'   exposure l is fixed to the \eqn{No Exposure} condition (i.e. no direct or
#'   indirect exposure).} \item{`var_tau_ht`:}{A named numeric vector which
#'   contains the values of the conservative variance estimator of the variance
#'   of the Horvitz-Thompson estimator of the average unit-level causal effect
#'   of exposure k versus exposure l as derived in Equation 11 of Aronow and
#'   Samii (2017). Here exposure l is fixed to the \eqn{No Exposure} condition
#'   (i.e. no direct or indirect exposure).} \item{`var_tau_h`:}{A named numeric
#'   vector which contains the values of the linearized variance estimator of
#'   the variance of the Hajek estimator of the average unit-level causal effect
#'   of exposure k versus exposure l as derived in Equation 11 of Aronow and
#'   Samii (2017) and further explained in the first paragraph of page 1929.
#'   Here exposure l is fixed to the \eqn{No Exposure} condition (i.e. no direct
#'   or indirect exposure).} \item{`var_tau_ht_const_eff`:}{A named numeric
#'   vector which contains the values of the constant effects variance
#'   estimator of the variance of the Horvitz-Thompson estimator of the average
#'   unit-level causal effect of exposure k versus exposure l as derived in
#'   Equation 2.15 of Aronow (2013). Here exposure l is fixed to the \eqn{No
#'   Exposure} condition (i.e. no direct or indirect exposure).}
#'   \item{`var_tau_ht_max`:}{A named numeric vector which contains the maximum
#'   between `var_tau_h` and `var_tau_ht_const_eff`.} }
estimates <-
  memoise(function(obs_exposure,
                   obs_outcome,
                   obs_prob_exposure,
                   n_var_permutations = 10,
                   effect_estimators = c('hajek', 'horvitz-thompson'),
                   variance_estimators = c('hajek', 'horvitz-thompson'),
                   hop) {

    if (('hajek' %in% variance_estimators) & ! ('hajek' %in% effect_estimators)) {
     effect_estimators  <- c(effect_estimators, 'hajek')
    }

     if ((
       ('horvitz-thompson' %in% variance_estimators)  | ('constant_effect' %in% variance_estimators) | ('max_ht_const' %in% variance_estimators)
       ) & ! ('horvitz-thompson' %in% effect_estimators)) {
     effect_estimators  <- c(effect_estimators, 'horvitz-thompson')
    }

    if (('max_ht_const' %in% variance_estimators) & !('constant_effect' %in% variance_estimators)) {
      variance_estimators <- c(variance_estimators, 'constant_effect')
    }
    if (('max_ht_const' %in% variance_estimators) & !('horvitz-thompson' %in% variance_estimators)) {
      variance_estimators <- c(variance_estimators, 'horvitz-thompson')

    }

    if (length(effect_estimators) == 0) {
      stop("Must specify effect estimators")
    }
    if (length(variance_estimators) == 0) {
      stop("Must specify variance estimators")
    }


    out <- list()

    N <- nrow(obs_exposure)
    
    obs_outcome_by_exposure <- t(obs_exposure) %*% diag(obs_outcome)
    obs_outcome_by_exposure[t(obs_exposure) == 0] <- NA
    
    obs_prob_exposure_individual_kk <-
      make_prob_exposure_cond(obs_prob_exposure)
    
    yT_ht <-
      rowSums(obs_outcome_by_exposure / obs_prob_exposure_individual_kk,
              na.rm = T)
    yT_ht[rowSums(!is.na(obs_outcome_by_exposure)) == 0] <- NA

    if (('horvitz-thompson' %in% effect_estimators)) {
      out[['yT_ht']] <- yT_ht
    }
    
    if (('hajek' %in% effect_estimators)) {
      mu_h <-
        yT_ht / rowSums(t(obs_exposure) / obs_prob_exposure_individual_kk, na.rm =
                          T)
      mu_h[rowSums(!is.na(obs_outcome_by_exposure)) == 0] <- NA
      
      yT_h <- mu_h * N
      yT_h[rowSums(!is.na(obs_outcome_by_exposure)) == 0] <- NA
    
      out[['yT_h']] <- yT_h
    }

    if (('horvitz-thompson' %in% variance_estimators)) {
      var_yT_ht <-
        var_yT_ht_adjusted(obs_exposure, obs_outcome, obs_prob_exposure)
      var_yT_ht[rowSums(!is.na(obs_outcome_by_exposure)) == 0] <- NA
      out[['var_yT_ht']] <- var_yT_ht
    }

    if (('hajek' %in% variance_estimators)) {
      resid_h <-
        colSums((obs_outcome_by_exposure - mu_h), na.rm = T) # pass resid_h instead of obs_outcome to var, cov

      var_yT_h <-
        var_yT_ht_adjusted(obs_exposure, resid_h, obs_prob_exposure)
      var_yT_h[rowSums(!is.na(obs_outcome_by_exposure)) == 0] <- NA

      out[['var_yT_h']] <- var_yT_h
    }
    
    if (('constant_effect' %in% variance_estimators)) {
      const_eff <-
        var_yT_ht_const_eff_lm(obs_exposure,
                              obs_outcome,
                              obs_prob_exposure,
                              n_var_permutations)
      var_yT_ht_const_eff <- const_eff$var_yT_ht_const_eff
      diffs <- const_eff$means - const_eff$means$no
      diffs <- subset(diffs, select = -c(no))
      var_tau_ht_const_eff <- unlist(lapply(diffs, var))
      # Count all of the cells which are not NA
      obs_outcome_by_exposure_na_conditions <-
        names(which(rowSums(!is.na(
          obs_outcome_by_exposure
        )) == 0))
      var_tau_ht_const_eff[obs_outcome_by_exposure_na_conditions] <- NA
      out[['var_tau_ht_const_eff']] <- var_tau_ht_const_eff
    }
    
    
   
    if (hop == 1) {
      remove <- 'no'
      keep <- c('dir_ind1,no', 'isol_dir,no', 'ind1,no')

      if (('horvitz-thompson' %in% variance_estimators)) {
        cov_yT_ht <-
          cov_yT_ht_adjusted(
            obs_exposure,
            obs_outcome,
            obs_prob_exposure,
            k_to_include = c('dir_ind1', 'isol_dir', 'ind1'),
            l_to_include = 'no'
          )
        out[['cov_yT_ht']] <- cov_yT_ht
        var_tau_ht <-
        (1 / N ^ 2) * (var_yT_ht[!rownames(var_yT_ht) %in% remove,] + var_yT_ht[rownames(var_yT_ht) %in% remove,] -
                         2 * cov_yT_ht[rownames(cov_yT_ht) %in% keep,])
        out[['var_tau_ht']] <- var_tau_ht
      }
      if (('hajek' %in% variance_estimators)) {
        cov_yT_h <-
          cov_yT_ht_adjusted(
            obs_exposure,
            resid_h,
            obs_prob_exposure,
            k_to_include = c('dir_ind1', 'isol_dir', 'ind1'),
            l_to_include = 'no'
          )
      out[['cov_yT_h']] <- cov_yT_h
      var_tau_h <-
        (1 / N ^ 2) * (var_yT_h[!rownames(var_yT_h) %in% remove,] + var_yT_h[rownames(var_yT_h) %in% remove,] -
                         2 * cov_yT_h[rownames(cov_yT_h) %in% keep,])
 
      out[['var_tau_h']] <- var_tau_h
      }
    }
    
    if (hop == 2) {
      remove <- 'no'
      keep <- c(
        'dir_ind1_ind2,no',
        'dir_ind1,no',
        'dir_ind2,no',
        'isol_dir,no',
        'ind1_ind2,no',
        'ind1,no',
        'ind2,no'
      )
      if (('horvitz-thompson' %in% variance_estimators)) {
        cov_yT_ht <-
          cov_yT_ht_adjusted(
            obs_exposure,
            obs_outcome,
            obs_prob_exposure,
            k_to_include = c(
              'dir_ind1_ind2',
              'dir_ind1',
              'dir_ind2',
              'isol_dir',
              'ind1_ind2',
              'ind1',
              'ind2'
            ),
            l_to_include = 'no'
          )
        var_tau_ht <-
        (1 / N ^ 2) * (var_yT_ht[!rownames(var_yT_ht) %in% remove,] + var_yT_ht[rownames(var_yT_ht) %in% remove,] -
                         2 * cov_yT_ht[rownames(cov_yT_ht) %in% keep,])

        out[['var_tau_ht']] <- var_tau_ht
      }
      if (('hajek' %in% variance_estimators)) {
        cov_yT_h <-
          cov_yT_ht_adjusted(
            obs_exposure,
            resid_h,
            obs_prob_exposure,
            k_to_include =  c(
              'dir_ind1_ind2',
              'dir_ind1',
              'dir_ind2',
              'isol_dir',
              'ind1_ind2',
              'ind1',
              'ind2'
            ),
            l_to_include = 'no'
          )
          var_tau_h <-
        (1 / N ^ 2) * (var_yT_h[!rownames(var_yT_h) %in% remove,] + var_yT_h[rownames(var_yT_h) %in% remove,] -
                         2 * cov_yT_h[rownames(cov_yT_h) %in% keep,])
          out[['var_tau_h']] <- var_tau_h
      }
    }
    
    
    if (('horvitz-thompson' %in% variance_estimators)) {
      tau_ht <- (1 / N) * (yT_ht - yT_ht['no'])[names(yT_ht) != 'no']
      out[['tau_ht']] <- tau_ht
    }
    if (('hajek' %in% effect_estimators)) {
      tau_h <- (mu_h - mu_h['no'])[names(mu_h) != 'no']
      out[['tau_h']] <- tau_h
    }
    
    if ('max_ht_const' %in% variance_estimators) {
      if (!setequal(names(var_tau_ht), names(var_tau_ht_const_eff))) {
        warning(
          "var_tau_ht and var_tau_ht_const_eff do not have the same names, var_tau_ht_max may be incorrect"
        )
      }
      var_tau_ht_max <-
        pmax(var_tau_ht[order(names(var_tau_ht))], var_tau_ht_const_eff[order(names(var_tau_ht_const_eff))])
      out[['var_tau_ht_max']] <- var_tau_ht_max
    }
   
    tau_dsm <-
      (
        rowMeans(obs_outcome_by_exposure, na.rm = T) - rowMeans(obs_outcome_by_exposure, na.rm = T)['no']
      )[names(yT_ht) != 'no']
    out[['tau_dsm']] <- tau_dsm
    
    
    return(out)
  })


#' @rdname estimates
#' @noRd
#' @export
make_prob_exposure_cond <- function(prob_exposure) {
  
  k_exposure_names <- stringi::stri_split_fixed(str = names(prob_exposure$prob_exposure_k_k), pattern=',', simplify=T)[,1]
  N <- nrow(prob_exposure$I_exposure[[1]])
  
  prob_exposure_cond <- matrix(nrow = length(k_exposure_names), ncol = N)
  for (j in 1:length(prob_exposure$prob_exposure_k_k)) {
    prob_exposure_cond[j,] <- diag(prob_exposure$prob_exposure_k_k[[j]])
  }
  rownames(prob_exposure_cond) <- k_exposure_names
  
  return(prob_exposure_cond)
}


#' @rdname estimates
#' @noRd
#' @export
var_yT_ht_unadjusted <- function(obs_exposure,obs_outcome,prob_exposure) {
  
  
  var_yT <- matrix(nrow = ncol(obs_exposure), dimnames=list(colnames(obs_exposure)))
  for (k in colnames(obs_exposure)) {
    pi_k <- prob_exposure$prob_exposure_k_k[[paste(k,k,sep=',')]]
    ind_kk <- diag(pi_k)
    cond_indicator = obs_exposure[,k]
    
    mm <- cond_indicator %o% cond_indicator * (pi_k - ind_kk %o% ind_kk)/pi_k * (obs_outcome %o% obs_outcome) / (ind_kk %o% ind_kk)
    mm[!is.finite(mm)] <- 0
    
    second_part_sum <- sum(mm)
    
    var_yT[k,] <- second_part_sum
  }
  return(var_yT)
}



#' @rdname estimates
#' @noRd
#' @export
var_yT_ht_A2_adjustment <- function(obs_exposure,obs_outcome,prob_exposure) {
  
  var_yT_A2 <- matrix(nrow = ncol(obs_exposure), dimnames=list(colnames(obs_exposure)))
  
  for (k in colnames(obs_exposure)) {
    pi_k <- prob_exposure$prob_exposure_k_k[[paste(k,k,sep=',')]]
    ind_kk <- diag(pi_k)
    cond_indicator = obs_exposure[,k]
    
    m <- cond_indicator*(obs_outcome^2)/(2*ind_kk)
    A2_part_sum <- sum(outer(m, m, FUN='+') * (pi_k == 0) * (!diag(length(m))) )
    var_yT_A2[k,] <- A2_part_sum
  }
  return(var_yT_A2)
}



#' @rdname estimates
#' @noRd
#' @export
var_yT_ht_adjusted <- function(obs_exposure,obs_outcome,prob_exposure) {
  var <- var_yT_ht_unadjusted(obs_exposure,obs_outcome,prob_exposure)
  A2 <- var_yT_ht_A2_adjustment(obs_exposure,obs_outcome,prob_exposure)
  var_adjusted <- var + A2
  return(var_adjusted)
}



#' @rdname estimates
#' @noRd
#' @export
cov_yT_ht_adjusted <- function(obs_exposure,obs_outcome,prob_exposure, k_to_include=NULL, l_to_include=NULL) {
  
  cov_yT_A <- matrix(nrow = 2*choose(ncol(obs_exposure),2), dimnames=list(names(prob_exposure$prob_exposure_k_l)))
  
  if (is.null(k_to_include)) {
    k_to_include <- colnames(obs_exposure) }
  if (is.null(l_to_include)) {
    l_to_include <- colnames(obs_exposure) }
  
  for (k in k_to_include) {
    for (l in l_to_include) {
      if (k!=l) {
        
        pi_k <- prob_exposure$prob_exposure_k_k[[paste(k,k,sep=',')]]
        ind_kk <- diag(pi_k)
        pi_l <- prob_exposure$prob_exposure_k_k[[paste(l,l,sep=',')]]
        ind_ll <- diag(pi_l)
        pi_k_l <- prob_exposure$prob_exposure_k_l[[paste(k,l,sep=',')]]
        cond_indicator_k = obs_exposure[,k]
        cond_indicator_l = obs_exposure[,l]
        
        mm <- cond_indicator_k %o% cond_indicator_l * (pi_k_l - ind_kk %o% ind_ll)/pi_k_l * (obs_outcome %o% obs_outcome) / (ind_kk %o% ind_ll)
        mm[!is.finite(mm)] <- 0
        
        first_part_cov <- sum(mm)
        
        second_part_cov <- 0
        for (i in 1:length(cond_indicator_k)) {
          for (j in 1:length(cond_indicator_l)) {
            if (pi_k_l[i,j]==0) {
              second_part_cov_i_j <- ((cond_indicator_k[i]*obs_outcome[i]^2/(2*ind_kk[i])) + (cond_indicator_l[j]*obs_outcome[j]^2/(2*ind_ll[j])))
              second_part_cov <- second_part_cov + second_part_cov_i_j
              
            }
          }
        }
        cov_yT_A[paste(k,l,sep=',') , ] <- first_part_cov - second_part_cov
      }
      
    }
    
  }
  return(cov_yT_A)
}



#' @rdname estimates
#' @noRd
#' @export
fit_avg_outcome_lm <- function(obs_exposure,obs_outcome, obs_prob_exposure){
  pi_ <- lapply(obs_prob_exposure$prob_exposure_k_k, diag)
  names(pi_) <- stringi::stri_split_fixed(str = names(pi_), pattern=',', simplify=T)[,1]
  
  pi_df <- do.call(cbind.data.frame, pi_)
  names(pi_df) <- paste0('pi_', names(pi_df))
  
  
  inv_pi_df <- 1 / (obs_exposure * pi_df)
  inv_pi_df[!sapply(inv_pi_df, is.finite)] <- 0
  wts <- rowSums(inv_pi_df)
  
  model_df <- as.data.frame(obs_exposure)
  model_df <- cbind(model_df, pi_df)
  model_df$obs_outcome <- obs_outcome
  
  lm1 <- lm(obs_outcome ~ ., data=model_df,weights = wts)
  
  obs_options <- rep(list(c(0, 1)), ncol(obs_exposure))
  names(obs_options) <-   colnames(obs_exposure)
  obs_options
  
  indicator_columns <- data.frame(do.call(
    rbind,
    c(
      list(rep(0, ncol(obs_exposure)-1)),
      unique(combinat::permn(
        c(1,rep(0, ncol(obs_exposure)-2)))
      )
    )
  ))
  colnames(indicator_columns) <- colnames(obs_exposure)[1:ncol(obs_exposure)-1]
  
  average_outcomes_df <- merge(
    indicator_columns,
    pi_df,
    by=NULL
  )
  
  average_outcomes_df[[colnames(obs_exposure)[ncol(obs_exposure)]]] <- NA_real_
  
  average_outcomes_df$y_hat <- predict(lm1, newdata=average_outcomes_df)
  average_outcomes_df <- data.table(average_outcomes_df)
  
  return(list(
    means=average_outcomes_df[, .(y=mean(y_hat,na.rm = TRUE)), by=eval(colnames(indicator_columns))],
    model=lm1
  )
  
  )
}




#' @rdname estimates
#' @noRd
#' @export
var_yT_ht_const_eff_lm <- function(obs_exposure,obs_outcome,obs_prob_exposure,n_var_permutations=1000) {
  
  #if (n_var_permutations > ncol(obs_prob_exposure$I_exposure$dir_ind1)) {
  # stop('n_var_permutations must be smaller or equal to R, the number of permutations to compute exposure probabilities')
  #}
  
  lm1 <- fit_avg_outcome_lm(obs_exposure,obs_outcome,obs_prob_exposure)$model
  
  residuals1 <- obs_outcome - predict(lm1)
  
  means_iter <- list()  
  for (i in 1:n_var_permutations) {
    obs_exposure_iter <- do.call(cbind.data.frame, lapply(obs_prob_exposure$I_exposure,function(x) x[,i]))
    this_means <- fit_avg_outcome_lm(obs_exposure_iter,residuals1,obs_prob_exposure)$means
    this_means$iter <- i
    means_iter[[i]] <- this_means
  }
  means_iter <- rbindlist(means_iter)
  
  
  ix <- apply(means_iter[,1:(ncol(means_iter)-1)],1,function(x) which(as.logical(x)))
  ix <- lapply(ix, function(x) ifelse(length(x)==0, ncol(means_iter), x))
  
  means_iter$exposure <- colnames(obs_exposure)[unlist(ix)]
  means_iter <- dcast(means_iter,iter ~ exposure, value.var='y')
  means_iter <- subset(means_iter, select=-c(iter))
  
  var <- unlist(lapply(means_iter, var))
  
  return(list(means=means_iter, var_yT_ht_const_eff=var))
  
}




#' @describeIn estimates Produces values for the estimator of the average
#'   unit-level causal effect of exposure k versus l and its variance estimator
#'   for the exposure mapping returned by function
#'   \code{\link{make_exposure_map_full_neighborhood}}, using a Horvitz-Thompson
#'   and a Hajek estimator. It also computes Horvitz-Thompson and Hajek
#'   estimators of the total of potential outcomes, which are inputs in the
#'   computation of average unit-level causal effect of exposure k versus l.
#' @examples
#' # Create adjacency matrix and treatment vector to
#' # produce observed exposure conditions according to the
#' # "full neighborhood" exposure mapping:
#'
#' adj_matrix <- make_adj_matrix(N = 81, model = 'sq_lattice')
#'
#' tr_vector <- make_tr_vec_permutation(N = 81, p = 0.5,
#'                                      R = 1, seed = 357)
#'
#' obs_exposure_full_nei <- make_exposure_map_full_neighborhood(adj_matrix,
#'                                                     tr_vector)
#' # Simulate a vector of outcome data:
#'
#' potential_outcome_full_nei <-
#'   make_dilated_out_full_neighborhood(adj_matrix, make_corr_out,
#'                                      seed = 357)
#'
#' obs_outcome_full_nei <-
#'   rowSums(obs_exposure_full_nei*t(potential_outcome_full_nei))
#'
#' # Create exposure probabilities:
#'
#' potential_tr_vector <- make_tr_vec_permutation(N = 81, p = 0.5,
#'                                                R = 36,
#'                                                seed = 357)
#' obs_prob_exposure_full_nei <- make_exposure_prob(potential_tr_vector,
#'                                         adj_matrix,
#'                                         make_exposure_map_full_neighborhood)
#'
#' # Estimate exposure-specific causal effects and their variance:
#'
#' estimators_full_neighborhood(obs_exposure_full_nei, obs_outcome_full_nei,
#'           obs_prob_exposure_full_nei,
#'           n_var_permutations = 30)
#' @export
estimators_full_neighborhood <-
  memoise(function(obs_exposure,
                   obs_outcome,
                   obs_prob_exposure,
                   n_var_permutations = 10) {
    N <- nrow(obs_exposure)
    obs_outcome_by_exposure <- t(obs_exposure) %*% diag(obs_outcome)
    obs_outcome_by_exposure[t(obs_exposure) == 0] <- NA
    
    obs_prob_exposure_individual_kk <-
      make_prob_exposure_cond(obs_prob_exposure)
    
    yT_ht <-
      rowSums(obs_outcome_by_exposure / obs_prob_exposure_individual_kk,
              na.rm = T)
    yT_ht[rowSums(!is.na(obs_outcome_by_exposure)) == 0] <- NA
    
    mu_h <-
      yT_ht / rowSums(t(obs_exposure) / obs_prob_exposure_individual_kk, na.rm =
                        T)
    mu_h[rowSums(!is.na(obs_outcome_by_exposure)) == 0] <- NA
    
    yT_h <- mu_h * N
    yT_h[rowSums(!is.na(obs_outcome_by_exposure)) == 0] <- NA
    
    resid_h <-
      colSums((obs_outcome_by_exposure - mu_h), na.rm = T) # pass resid_h instead of obs_outcome to var, cov to get Hajek estimator
    
    var_yT_ht <-
      var_yT_ht_adjusted(obs_exposure, obs_outcome, obs_prob_exposure)
    
    var_yT_h <-
      var_yT_ht_adjusted(obs_exposure, resid_h, obs_prob_exposure)
    
    cov_yT_ht <-
      cov_yT_ht_adjusted(obs_exposure, obs_outcome, obs_prob_exposure)
    cov_yT_h <-
      cov_yT_ht_adjusted(obs_exposure, resid_h, obs_prob_exposure)
    
    const_eff <-
      var_yT_ht_const_eff_lm(obs_exposure,
                             obs_outcome,
                             obs_prob_exposure,
                             n_var_permutations)
    
    var_yT_ht_const_eff <- const_eff$var_yT_ht_const_eff
    
    diffs <- const_eff$means - const_eff$means$all_control
    diffs <- subset(diffs, select = -c(all_control))
    
    var_tau_ht_const_eff <- unlist(lapply(diffs, var))
    
    remove <- 'all_control'
    keep <- c('all_treat,all_control')
    
    var_tau_ht <-
      (1 / N ^ 2) * (var_yT_ht[!rownames(var_yT_ht) %in% remove,] + var_yT_ht[rownames(var_yT_ht) %in% remove,] -
                       2 * cov_yT_ht[rownames(cov_yT_ht) %in% keep,])
    
    var_tau_h <-
      (1 / N ^ 2) * (var_yT_h[!rownames(var_yT_h) %in% remove,] + var_yT_h[rownames(var_yT_h) %in% remove,] -
                       2 * cov_yT_h[rownames(cov_yT_h) %in% keep,])
    
    var_tau_ht_max <-
      pmax(var_tau_ht[order(names(var_tau_ht))], var_tau_ht_const_eff[order(names(var_tau_ht_const_eff))])
    
    tau_ht <-
      (1 / N) * (yT_ht - yT_ht['all_control'])[names(yT_ht) != 'all_control']
    tau_h <- (mu_h - mu_h['all_control'])[names(mu_h) != 'all_control']
    
    return(
      list(
        yT_ht = yT_ht,
        yT_h = yT_h,
        tau_ht = tau_ht,
        tau_h = tau_h,
        var_tau_ht = var_tau_ht,
        var_tau_h = var_tau_h,
        var_tau_ht_const_eff = var_tau_ht_const_eff,
        var_tau_ht_max = var_tau_ht_max
      )
    )
  })
