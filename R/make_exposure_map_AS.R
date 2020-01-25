#' Create treatment exposure conditions.
#'
#' Create treatment exposure conditions according to selected exposure mappings.
#'
#' \code{make_exposure_map_AS} produces treatment exposure conditions according
#' to two different exposure mappings. An exposure mapping that assumes that
#' interference happens only through direct peer connections (first-degree
#' interference) which produces four exposure conditions: \eqn{Direct + Indirect
#' Exposure}, \eqn{Isolated Direct Exposure}, \eqn{Indirect Exposure}, \eqn{No
#' Exposure}. And an exposure mapping that assumes second-degree interference
#' which produces eight exposure conditions.
#'
#' @param adj_matrix an `N` \eqn{*} `N` numeric matrix of 0, 1 entries such as those
#'   retuned by \code{\link{make_adj_matrix}}, or a `N` \eqn{*} `N` matrix, where `N`
#'   is the number of units.
#' @param tr_vector a numeric vector of length `N` of 0, 1 entries such as that
#'   retuned by \code{\link{make_tr_vec_permutation}} with `R = 1`, or a vector
#'   of length `N` with the treatment assignment.
#' @param hop number; either `1` or `2`. If `1` assumes first-degree
#'   interference and produces four exposure conditions, if `2` assumes
#'   second-degree interference and produces eight exposure conditions.
#' @references Aronow, P.M. & Samii, C. (2017). [Estimating average causal
#'   effects under general interference, with application to a social network
#'   experiment](https://doi.org/10.1214/16-AOAS1005). *The Annals of Applied
#'   Statistics*, 11(4), 1912--1947.
#'   
#'   Aronow, P.M. et al. (2020). [Spillover effects in experimental
#'   data](https://arxiv.org/abs/2001.05444). *arXiv preprint*,
#'   arXiv:2001.05444.
#' @export
#' @examples
#' # Create an adjacency matrix and a treatment vector:
#' 
#' adj_matrix <- make_adj_matrix(N = 9, model = 'sq_lattice')
#' tr_vector <- make_tr_vec_permutation(N = 9, p = 0.2, R = 1,
#'                                      seed = 357)
#' 
#' # Create treatment exposure conditions:
#' 
#' make_exposure_map_AS(adj_matrix, tr_vector, hop = 1)
#' make_exposure_map_AS(adj_matrix, tr_vector, hop = 2)
#' @return An `N` \eqn{*} K numeric matrix, where K corresponds to the number of
#'   exposure conditions.
make_exposure_map_AS <- function(adj_matrix, tr_vector, hop) {
  N <- nrow(adj_matrix)
  peer_exposure <- as.numeric(tr_vector %*% adj_matrix)
  if (hop == 1) {
    return(matrix(as.numeric(
      c(
        tr_vector > 0 & peer_exposure > 0,
        tr_vector > 0 & peer_exposure == 0,
        tr_vector == 0 & peer_exposure > 0,
        tr_vector == 0 & peer_exposure == 0
      )
    ),
    N, 4, dimnames = list(
      NULL, c('dir_ind1', 'isol_dir', 'ind1', 'no')
    )))
  }
  if (hop == 2) {
    adj_matrix_sq <-
      adj_matrix %*% adj_matrix # number of length two paths from i to j
    adj_matrix_sq[which(adj_matrix_sq > 1)] <- 1
    diag(adj_matrix_sq) <- rep(0, nrow(adj_matrix_sq))
    adj_matrix_2 <- adj_matrix_sq
    
    peer_exposure_2 <- as.numeric(tr_vector %*% adj_matrix_2)
    return(matrix(as.numeric(
      c(
        tr_vector > 0 & peer_exposure > 0 & peer_exposure_2 > 0,
        tr_vector > 0 &
          peer_exposure > 0 & peer_exposure_2 == 0,
        tr_vector > 0 &
          peer_exposure == 0 & peer_exposure_2 > 0,
        tr_vector > 0 &
          peer_exposure == 0 & peer_exposure_2 == 0,
        tr_vector == 0 &
          peer_exposure > 0 & peer_exposure_2 > 0,
        tr_vector == 0 &
          peer_exposure > 0 & peer_exposure_2 == 0,
        tr_vector == 0 &
          peer_exposure == 0 & peer_exposure_2 > 0,
        tr_vector == 0 &
          peer_exposure == 0 & peer_exposure_2 == 0
      )
    ),
    N, 8, dimnames = list(
      NULL,
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
    )))
  }
}


#' @describeIn make_exposure_map_AS Create treatment exposure conditions
#'   according to an exposure mapping that assumes first-degree interference and
#'   that produces two exposure conditions: \eqn{All Treat}; a unit and 100% of
#'   its first-degree neighbors are treated, and \eqn{All Control}; a unit and
#'   100% of its first-degree neighbors are in the control condition.
#' @examples
#' make_exposure_map_full_neighborhood(adj_matrix, tr_vector)
#' @export
make_exposure_map_full_neighborhood <- function(adj_matrix,tr_vector) {
  N <- nrow(adj_matrix)
  adj_matrix_diag <- adj_matrix
  diag(adj_matrix_diag) <- rep(1, nrow(adj_matrix_diag))
  exposure <- as.numeric(tr_vector%*%adj_matrix_diag)
  B_size <- rowSums(adj_matrix_diag)
  return(matrix(as.numeric(c(exposure==B_size,
                             exposure<B_size)),
                N, 2, dimnames = list(NULL, c('all_treat', 'all_control'))))
}
