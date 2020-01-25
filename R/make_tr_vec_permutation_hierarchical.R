#'Randomly assign treatment with uniform probability under a two-stage
#'hierarchical design.
#'
#'Create possible treatment random assignments with uniform probability under a
#'two-stage hierarchical design.
#'
#'\code{make_tr_vec_permutation_hierarchical} produces all possible (or a random
#'sample of all possible) treatment assignments under a two-stage hierarchical
#'design where groups are first randomly assigned to a high or a low level of
#'treatment saturation (psi, phi), and then units within groups are randomly
#'assigned to treatment with probability equal to their group saturation rate.
#'Sampling of treatment vectors is without replacement.
#'
#'@param group integer vector of length `N` of values specifying the `N` units'
#'  group membership (or group labels).
#'@param c proportion of groups assigned to the high saturation (psi) in the
#'  first stage.
#'@param k numeric vector of length 2 of proportions of individuals whithin
#'  groups assigned to treatment in the second stage, given the first stage
#'  treatment low and high saturations (phi and psi), respectively.
#' @param R number of repetitions (treatment permutations).  
#'@inheritParams make_tr_vec_permutation
#'@export
#'@return list of `R` data frames each of `N` observations and 3 variables:
#'  \describe{ \item{`group`:}{Integer with group label.}
#'  \item{`group_tr`:}{Numeric indicator of group assignment to high saturation
#'  (psi) in the first stage.} \item{`indiv_tr`:}{Numeric indicator of
#'  individual assignment to treatment in the second stage.} }
#' @examples
#' # Simulate a vector of group membership of 30 units
#' # equally divided into 6 groups:
#'
#' group <- rep(1:6, each = 30/6)
#'
#' # Define proportion of groups assigned to
#' # saturation psi in the first stage:
#'
#' c <- 1/2
#'
#' # Define proportion of individuals whithin groups assigned
#' # to treatment in the second stage, given the assigned
#' # first-stage low and high (phi, psi) saturation rates:
#'
#' k <- c(2/5, 3/5)
#'
#' make_tr_vec_permutation_hierarchical(group, c, k, R = 1, seed = 357)
#' make_tr_vec_permutation_hierarchical(group, c, k, R = 5, seed = 357)
#'
#'@references Hudgens, M.G. & Halloran M.E. (2008). [Toward causal inference
#'  with interference](https://doi.org/10.1198/016214508000000292). *Journal of
#'  the American Statistical Association*, 103(482), 832--842.
#'
#'  Aronow, P.M. et al. (2020). [Spillover effects in experimental
#'  data](https://arxiv.org/abs/2001.05444). *arXiv preprint*, arXiv:2001.05444.
make_tr_vec_permutation_hierarchical <-
  function(group, c, k, R, seed = NULL) {
    set.seed(seed)
    tr_vec_sampled <- vector('list', R)
    G <- length(unique(group))
    
    for (j in 1:R) {
      group_tr <-
        as.list(sample(c(rep(1, round(
          G * c
        )), rep(0, G - round(
          G * c
        )))))
      names(group_tr) <- unique(group)
      
      indiv_tr <- rep(NA, length(group))
      for (i in unique(group)) {
        n <- sum(group == i)
        p <- k[group_tr[[i]] + 1]
        indiv_tr[which(group == i)] <-
          sample(c(rep(1, round(n * p)), rep(0, n - round(n * p))))
      }
      
      group_tr <- as.data.frame(unlist(group_tr))
      colnames(group_tr) <- 'group_tr'
      group_tr$group <- unique(group)
      
      group_tr_data <-
        merge(as.data.frame(group), as.data.frame(group_tr), by = 'group')
      vec <- cbind(group_tr_data, as.data.frame(indiv_tr))
      
      while (any(unlist(lapply(tr_vec_sampled, function(x) {
        !is.null(x) && all(vec == x)
      }))))
      {
        group_tr <-
          as.list(sample(c(rep(1, round(
            G * c
          )), rep(
            0, G - round(G * c)
          ))))
        names(group_tr) <- unique(group)
        
        indiv_tr <- rep(NA, length(group))
        for (i in unique(group)) {
          n <- sum(group == i)
          p <- k[group_tr[[i]] + 1]
          indiv_tr[which(group == i)] <-
            sample(c(rep(1, round(n * p)), rep(0, n - round(n * p))))
        }
        
        group_tr <- as.data.frame(unlist(group_tr))
        colnames(group_tr) <- 'group_tr'
        group_tr$group <- unique(group)
        
        group_tr_data <-
          merge(as.data.frame(group), as.data.frame(group_tr), by = 'group')
        vec <- cbind(group_tr_data, as.data.frame(indiv_tr))
      }
      tr_vec_sampled[[j]] <- vec
    }
    return(tr_vec_sampled)
  }


#' @rdname group_labels
#' @noRd
#' @export
make_group_labels <- function(N, G) {
  if ((N %% G) != 0) {
    stop('N must be divisible by G')
  }
  
  rep(1:G, each = N / G)
}

