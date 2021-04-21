#'Create treatment conditions under a two-stage hierarchical design.
#'
#'Create treatment conditions under two-stage hierarchical treatment assignment,
#'assuming partial and stratified interference.
#'
#'\code{make_tr_condition} produces treatment conditions according to a
#'two-stage hierarchical design where groups are first randomly assigned to a
#'high or a low level of treatment saturation (psi, phi), and then units within
#'groups are randomly assigned to treatment with probability equal to their
#'group saturation rate. Following Hudgens and Halloran (2008), there are four
#'potential outcomes which correspond to the four potential treatment
#'conditions: \eqn{Direct + Indirect Psi}--a unit is directly treated and its
#'group is assigned treatment saturation psi, \eqn{Direct + Indirect Phi}--a
#'unit is directly treated and its group is assigned treatment saturation phi,
#'\eqn{Indirect Psi}--a unit is not directly treated and its group is assigned
#'treatment saturation psi, \eqn{Indirect Phi}--a unit is not directly treated
#'and its group is assigned treatment saturation phi. This function assumes
#'stratified interference (i.e. potential outcomes of a unit are affected by its
#'own treatment assignment and only the treated proportion of its group; the
#'precise set of treated group members does not matter).
#'@param tr_assignment data frame of `N` observations and the variables:
#'  \describe{ \item{`group_tr`:}{Indicator of group assignment to high
#'  saturation (psi) in the first stage.} \item{`indiv_tr`:}{Indicator of
#'  individual assignment to treatment in the second stage.} } This data frame
#'  is returned by function \code{\link{make_tr_vec_permutation_hierarchical}}.
#'@return named numeric `N` \eqn{*} K matrix of observed treatment conditions K
#'  for units `N`.
#' @examples
#' # Create data frame with first and second stage
#' # treatment assignment indicators:
#'
#' group <- rep(1:6, each = 30/6)
#' c <- 1/2
#' k <- c(2/5, 3/5)
#'
#' tr_assignment <-
#' make_tr_vec_permutation_hierarchical(group, c, k, R = 1,
#'                                      seed = 357)[[1]]
#'
#' # Create data frame with observed treatment conditions:
#' make_tr_condition(tr_assignment)
#'@export
#'@references Hudgens, M.G. & Halloran M.E. (2008). [Toward causal inference
#'  with interference](https://doi.org/10.1198/016214508000000292). *Journal of
#'  the American Statistical Association*, 103(482), 832--842.
#'
#'  Aronow, P.M. et al. (2020). [Spillover effects in experimental
#'  data](https://arxiv.org/abs/2001.05444). *arXiv preprint*, arXiv:2001.05444.
make_tr_condition <- function(tr_assignment) {
  N <- nrow(tr_assignment)
  group_tr <- tr_assignment[, 'group_tr']
  indiv_tr <- tr_assignment[, 'indiv_tr']
  return(matrix(as.numeric(
    c(
      group_tr > 0 & indiv_tr > 0,
      group_tr == 0 & indiv_tr > 0,
      group_tr > 0 & indiv_tr == 0,
      group_tr == 0 & indiv_tr == 0
    )
  ),
  N, 4, dimnames = list(
    NULL, c('dir_indpsi', 'dir_indphi', 'indpsi', 'indphi')
  )))
}
