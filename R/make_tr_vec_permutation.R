#' Randomly assign treatment with uniform probability.
#'
#' Create a permutation matrix from possible treatment random assignments with
#' uniform probability.
#'
#' \code{make_tr_vec_permutation} produces all possible (or a random sample of
#' all possible) treatment assignments when sampling without replacement, or a
#' random sample of possible treatment assignments when sampling with
#' replacement.
#'
#' @param N number of units.
#' @param p proportion of units assigned to treatment.
#' @param R number of repetitions (treatment permutations). If `allow_repetions =
#'   FALSE` and `R` is bigger than the number of possible treatment assignments,
#'   then `R` is truncated to the number of possible treatment assignements.
#' @param seed random number for result replicability.
#' @param allow_repetitions logical; if `TRUE` sampling is with replacement
#'   rather than without replacement.
#' @export
#' @examples
#' make_tr_vec_permutation(N = 10, p = 0.2, R = 1, seed = 357)
#' make_tr_vec_permutation(N = 10, p = 0.2, R = 50, seed = 357)
#' make_tr_vec_permutation(N = 10, p = 0.2, R = 50, seed = 357,
#'                         allow_repetitions = TRUE)
#' @return An `R` \eqn{*} `N` numeric matrix. Each row cooresponds to a treatment
#'   assignment vector.
#'   
make_tr_vec_permutation <-
  function(N,
           p,
           R,
           seed = NULL,
           allow_repetitions = FALSE) {
    set.seed(seed)
    n_treated <- round(N * p)
    
    max_R <- choose(N, n_treated)
    if ((R > max_R) & (allow_repetitions == FALSE)) {
      R <- max_R
      warning(
        paste(
          "R is larger than the number of possible treatment assignements truncating to",
          max_R
        )
      )
    }
    tr_vec_sampled <- matrix(nrow = R, ncol = N)
    
    for (i in 1:R) {
      vec <- sample(c(rep(1, n_treated), rep(0, N - n_treated)))
      while (any(duplicated(rbind(vec, tr_vec_sampled[1:i - 1, ]))) &
             (allow_repetitions == FALSE))
      {
        vec <- sample(c(rep(1, n_treated), rep(0, N - n_treated)))
        
      }
      tr_vec_sampled[i, ] <- vec
    }
    return(tr_vec_sampled)
  }
