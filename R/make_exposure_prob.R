#' Create exposure probabilities.
#'
#' Create matrices of exposure probabilities.
#'
#' \code{make_exposure_prob} produces the units' probabilities of being subject
#' to each of the possible exposure conditions which are used for estimating
#' exposure-specific causal effects, and the joint exposure probabilities used
#' for variance estimators.
#'
#' @param potential_tr_vector an `R` \eqn{*}  `N` matrix of 0, 1 entries such as that
#'   produced by \code{\link{make_tr_vec_permutation}}, or an `R` \eqn{*}  `N` matrix
#'   containing `R` permuted treatment assignments for `N` units.
#' @param exposure_map_fn function which returns the exposure mapping such as
#'   \code{\link{make_exposure_map_AS}}.
#' @param exposure_map_fn_add_args list of additional arguments which are passed
#'   to `exposure_map_fn`. `adj_matrix` and `tr_vector` should not be specified
#'   here because they are already passed with values derived from the arguments
#'   to `make_exposure_prob`. If using \code{\link{make_exposure_map_AS}}, `hop`
#'   must be specified here.
#' @inheritParams make_exposure_map_AS
#' @export
#' @references Aronow, P.M. & Samii, C. (2017). [Estimating average causal
#'   effects under general interference, with application to a social network
#'   experiment](https://doi.org/10.1214/16-AOAS1005). *The Annals of Applied
#'   Statistics*, 11(4), 1912--1947.
#'
#'   Aronow, P.M. et al. (2020). [Spillover effects in experimental
#'   data](https://arxiv.org/abs/2001.05444). *arXiv preprint*,
#'   arXiv:2001.05444.
#' @examples
#' potential_tr_vector <- make_tr_vec_permutation(N = 9, p = 0.2, R = 20, seed = 357)
#' adj_matrix <- make_adj_matrix(N = 9, model = 'sq_lattice')
#' make_exposure_prob(potential_tr_vector, adj_matrix, make_exposure_map_AS, list(hop=1))
#' @return A list of 3 lists: \describe{ \item{`I_exposure`:}{A list of K
#'   `N` \eqn{*}  `R` numeric matrices of indicators for whether units `N` are in
#'   exposure condition k over each of the possible `R` treatment assignment
#'   vectors. The number of numeric matrices K corresponds to the number of
#'   exposure conditions.} \item{`prob_exposure_k_k`:}{A list of K symmetric
#'   `N` \eqn{*}  `N` numeric matrices each containing individual exposure probabilities
#'   to condition k on the diagonal, and joint exposure probabilities to
#'   condition k on the off-diagonals.} \item{`prob_exposure_k_l`:}{A list of
#'   \eqn{permutation(K,2)} nonsymmetric `N`\eqn{*} `N` numeric matrices each
#'   containing joint probabilities across exposure conditions k and l on the
#'   off-diagonal, and zeroes on the diagonal. When K = 4, the number of numeric
#'   matrices is 12; \eqn{permutation(4,2)}.}}
make_exposure_prob <-  memoise(function(potential_tr_vector,
                                adj_matrix,
                                exposure_map_fn,
                                exposure_map_fn_add_args = NULL,
                                i_start = NULL,
                                i_end = NULL) {
  exposure_map_fn_args <-
    c(list(adj_matrix, potential_tr_vector[1, ]),
      exposure_map_fn_add_args)
  exposure_names <-
    colnames(do.call(exposure_map_fn, exposure_map_fn_args))
  n_exposure_conditions <- length(exposure_names)
  
  
  R <- nrow(potential_tr_vector)
  N <- ncol(potential_tr_vector)
  
  
  if ((is.null(i_start)) & (is.null(i_end))) {
    i_start <-1
    i_end <- N
  }
  
  # TODO: do not recreate I_exposure everytime
  I_exposure <- list()
  for (i in 1:n_exposure_conditions) {
    I_exposure[[i]] <- matrix(nrow = N, ncol = R)
  }
  names(I_exposure) <- exposure_names
  
  
  for (i in 1:R) {
    exposure_map_fn_args <-
      c(list(adj_matrix, potential_tr_vector[i, ]),
        exposure_map_fn_add_args)
    potential_exposure <-
      do.call(exposure_map_fn, exposure_map_fn_args)
    for (j in 1:ncol(potential_exposure)) {
      I_exposure[[j]][, i] <- potential_exposure[, j]
    }
  }
  
  
  
  
  prob_exposure_k_k <- list()
  prob_exposure_k_l <- list()
  for (i in 1:length(I_exposure)) {
    for (j in 1:length(I_exposure)) {
      #TODO: calculate this only if j==i
      diagonal_ones <- rbind(
        matrix(0, nrow=(i_start-1), ncol=i_end-i_start+1),
        diag(i_end-i_start+1), 
        matrix(0, nrow=N-(i_end), ncol=i_end-i_start+1)
      )
      prob_exposure_k_k[[paste(names(I_exposure)[[i]], names(I_exposure)[[i]], sep =
                                 ',')]] <-
        t(( I_exposure[[i]] %*% t(I_exposure[[i]])[,i_start:i_end] + diagonal_ones) / (R + 1))
      
      if (j != i) {
        prob_exposure_k_l[[paste(names(I_exposure)[[i]], names(I_exposure)[[j]], sep =
                                   ',')]] <-
          t((I_exposure[[i]] %*% t(I_exposure[[j]])[,i_start:i_end]) / R)
        
      }
      
    }
  }
  
  return(
    list(
      I_exposure = I_exposure,
      prob_exposure_k_k = prob_exposure_k_k,
      prob_exposure_k_l = prob_exposure_k_l
    )
  )
})