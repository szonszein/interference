#' Create potential dilated outcomes.
#'
#' Create potential dilated outcomes for simulation.
#'
#' \code{make_dilated_out} produces potential dilated outcomes for simulation
#' according to two different exposure mappings. An exposure mapping that
#' assumes that interference happens only through direct peer connections
#' (first-degree interference) which produces four exposure conditions:
#' \eqn{Direct + Indirect Exposure}, \eqn{Isolated Direct Exposure},
#' \eqn{Indirect Exposure}, \eqn{No Exposure}. And an exposure mapping that
#' assumes second-degree interference which produces eight exposure conditions.
#' The values for the baseline \eqn{No Exposure} condition are drawn from an
#' absolute standard normal distribution which is correlated with the unit's
#' first and second order degree. The values for the other exposure conditions
#' are obtained by multiplying the vector of multipliers by the baseline \eqn{No
#' Exposure} value.
#'
#' @param make_corr_out unused.
#' @param multipliers numeric vector with dilated effects multipliers. Must be
#'   of length 3 if `hop` is 1, and length 7 if `hop` is 2. Default are
#'   `c(2,1.5,1.25)` and `c(2.25,2,1.75,1.5,1.375,1.25,1.125)`.
#' @inheritParams make_tr_vec_permutation
#' @inheritParams make_exposure_map_AS
#' @examples
#' adj_matrix <- make_adj_matrix(N = 9, model = 'sq_lattice')
#' multipliers <- c(4, 2, 3)
#' make_dilated_out(adj_matrix, make_corr_out, seed = 357,
#' multipliers = multipliers, hop = 1)
#' make_dilated_out(adj_matrix, make_corr_out, seed = 357,
#' multipliers = multipliers, hop = 2)
#' @return An K \eqn{*} `N` named numeric matrix, where K corresponds to the number of
#'   exposure conditions and `N` number of units.
#' @export
#' @references Aronow, P.M. et al. (2020). [Spillover effects in experimental
#'   data](https://arxiv.org/abs/2001.05444). *arXiv preprint*,
#'   arXiv:2001.05444.
make_dilated_out <-
  function(adj_matrix,
           make_corr_out,
           seed,
           hop,
           multipliers = NULL) {
    if (hop == 1) {
      return(
        make_dilated_out_1(
          adj_matrix,
          make_corr_out,
          multipliers = multipliers,
          seed = seed
        )
      )
    }
    if (hop == 2) {
      return(
        make_dilated_out_2(
          adj_matrix,
          make_corr_out,
          multipliers = multipliers,
          seed = seed
        )
      )
    }
  }


#' @describeIn make_dilated_out Produces potential dilated outcomes for simulation according to
#'   an exposure mapping that assumes first-degree interference and that
#'   produces two exposure conditions: \eqn{All Treat}, a condition in which a
#'   unit and 100% of its first-degree neighbors are treated, and \eqn{All Control},
#'   a unit and 100% of its first-degree neighbors are in the control
#'   condition.
#' @export
#' @examples
#' multiplier <- 4
#' make_dilated_out_full_neighborhood(adj_matrix, make_corr_out,
#'                                    multipliers=multiplier,
#'                                    seed=357)    
make_dilated_out_full_neighborhood <-
  function(adj_matrix,
           make_corr_out,
           multipliers = NULL,
           seed = NULL) {
    set.seed(seed)
    if (is.null(multipliers)) {
      multipliers = c(2)
    }
    if (length(multipliers) != 1) {
      stop('Needs 1 multiplier')
    }
    degree <- rowSums(adj_matrix)
    nei2 <- adj_matrix %*% adj_matrix
    nei2[which(nei2 > 1)] <- 1
    diag(nei2) <- rep(0, nrow(nei2))
    degree2 <- rowSums(nei2)
    
    baseline_out <- make_corr_out(degree, degree2, 'yes', seed = seed)
    potential_out <- rbind(multipliers[1] * baseline_out, baseline_out)
    rownames(potential_out) <- c('all_treat', 'all_control')
    return(potential_out)
    
  }


#' @rdname dilated_out
#' @noRd
#' @export
make_corr_out <- function(degree, degree2, correlate, seed = NULL) {
  set.seed(seed)
  switch(correlate, 'yes' = return((degree + degree2 + degree * degree2) *
                                     abs(rnorm(length(degree))) + rnorm(length(degree), 1, 0.25)
  ),
  'no' = return(abs(rnorm(length(
    degree
  )))))
}


#' @rdname dilated_out
#' @noRd
#' @export
make_dilated_out_1 <-
  function(adj_matrix,
           make_corr_out,
           multipliers = NULL,
           seed = NULL) {
    set.seed(seed)
    if (is.null(multipliers)) {
      multipliers = c(2, 1.5, 1.25)
    }
    if (length(multipliers) != 3) {
      stop('Needs 3 multipliers')
    }
    
    degree <- rowSums(adj_matrix)
    nei2 <- adj_matrix %*% adj_matrix
    nei2[which(nei2 > 1)] <- 1
    diag(nei2) <- rep(0, nrow(nei2))
    degree2 <- rowSums(nei2)
    
    baseline_out <- make_corr_out(degree, degree2, 'yes', seed = seed)
    potential_out <-
      rbind(
        multipliers[1] * baseline_out,
        multipliers[2] * baseline_out,
        multipliers[3] * baseline_out,
        baseline_out
      )
    rownames(potential_out) <- c('dir_ind1', 'isol_dir', 'ind1', 'no')
    return(potential_out)
    
  }


#' @rdname dilated_out
#' @noRd
#' @export
make_dilated_out_2 <- function(adj_matrix,
                               make_corr_out,
                               multipliers = NULL,
                               seed = NULL) {
  set.seed(seed)
  if (is.null(multipliers)) {
    multipliers = c(2.25, 2, 1.75, 1.5, 1.375, 1.25, 1.125)
  }
  if (length(multipliers) != 7) {
    stop('Needs 7 multipliers')
  }
  
  degree <- rowSums(adj_matrix)
  nei2 <- adj_matrix %*% adj_matrix
  nei2[which(nei2 > 1)] <- 1
  diag(nei2) <- rep(0, nrow(nei2))
  degree2 <- rowSums(nei2)
  
  baseline_out <- make_corr_out(degree, degree2, 'yes', seed = seed)
  potential_out <-
    rbind(
      multipliers[1] * baseline_out,
      multipliers[2] * baseline_out,
      multipliers[3] * baseline_out,
      multipliers[4] * baseline_out,
      multipliers[5] * baseline_out,
      multipliers[6] * baseline_out,
      multipliers[7] * baseline_out,
      baseline_out
    )
  rownames(potential_out) <-
    c(
      'dir_ind1_ind2',
      'dir_ind1',
      'dir_ind2',
      'isol_dir',
      'ind1_ind2',
      'ind1',
      'ind2',
      'no'
    )
  return(potential_out)
  
}



#' @rdname dilated_out
#' @noRd
#' @export
make_dilated_out_from_out <-
  function(outcomes, hop, multipliers = NULL) {
    if (hop == 1) {
      if (is.null(multipliers)) {
        multipliers = c(2, 1.5, 1.25)
      }
      if (length(multipliers) != 3) {
        stop('Needs 3 multipliers')
      }
      baseline_out <- outcomes
      potential_out <-
        rbind(
          multipliers[1] * baseline_out,
          multipliers[2] * baseline_out,
          multipliers[3] * baseline_out,
          baseline_out
        )
      rownames(potential_out) <-
        c('dir_ind1', 'isol_dir', 'ind1', 'no')
      return(potential_out)
    }
    
    if (hop == 2) {
      if (is.null(multipliers)) {
        multipliers = c(2.25, 2, 1.75, 1.5, 1.375, 1.25, 1.125)
      }
      if (length(multipliers) != 7) {
        stop('Needs 7 multipliers')
      }
      baseline_out <- outcomes
      potential_out <-
        rbind(
          multipliers[1] * baseline_out,
          multipliers[2] * baseline_out,
          multipliers[3] * baseline_out,
          multipliers[4] * baseline_out,
          multipliers[5] * baseline_out,
          multipliers[6] * baseline_out,
          multipliers[7] * baseline_out,
          baseline_out
        )
      rownames(potential_out) <-
        c(
          'dir_ind1_ind2',
          'dir_ind1',
          'dir_ind2',
          'isol_dir',
          'ind1_ind2',
          'ind1',
          'ind2',
          'no'
        )
      return(potential_out)
    }
  }


