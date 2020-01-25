#' Create potential dilated outcomes under a two-stage hierarchical design.
#'
#' Create potential dilated outcomes for simulation under two-stage hierarchical
#' treatment assignment, assuming partial and stratified interference.
#'
#' \code{make_dilated_out_hh} produces potential dilated outcomes for simulation
#' according to a two-stage hierarchical design where groups are first randomly
#' assigned to a high or a low level of treatment saturation (psi, phi), and
#' then units within groups are randomly assigned to treatment with probability
#' equal to their group saturation rate. Following Hudgens and Halloran (2008),
#' there are four potential outcomes which correspond to the four potential
#' treatment conditions: \eqn{Direct + Indirect Psi}--a unit is directly treated
#' and its group is assigned treatment saturation psi, \eqn{Direct + Indirect
#' Phi}--a unit is directly treated and its group is assigned treatment
#' saturation phi, \eqn{Indirect Psi}--a unit is not directly treated and its
#' group is assigned treatment saturation psi, \eqn{Indirect Phi}--a unit is not
#' directly treated and its group is assigned treatment saturation phi. This
#' function assumes stratified interference (i.e. potential outcomes of a unit
#' are affected by its own treatment assignment and only the treated proportion
#' of its group; the precise set of treated group members does not matter).
#'
#' @param multipliers numeric vector with dilated effects multipliers. Must be
#'   of length 3. Default is `c(2,1.5,1.25)`.
#' @param out function which returns a vector length `N` of outcome values for
#'   treatment condition \eqn{Indirect Phi} which is the baseline (or control)
#'   condition. Default function draws from an absolute standard normal
#'   distribution.
#' @inheritParams make_tr_vec_permutation
#' @examples
#' make_dilated_out_hh(N = 10, seed = 357)
#'
#' multipliers <- c(4, 3, 2)
#' make_dilated_out_hh(N = 10, multipliers = multipliers, seed = 357)
#' @return An K \eqn{*} `N` named numeric matrix, where K corresponds to the
#'   number of exposure conditions and `N` number of units.
#' @export
#' @references Hudgens, M.G. & Halloran M.E. (2008). [Toward causal inference
#'   with interference](https://doi.org/10.1198/016214508000000292). *Journal of
#'   the American Statistical Association*, 103(482), 832--842.
#'
#'   Aronow, P.M. et al. (2020). [Spillover effects in experimental
#'   data](https://arxiv.org/abs/2001.05444). *arXiv preprint*,
#'   arXiv:2001.05444.
make_dilated_out_hh <-
  function(N,
           out = function(x)
             abs(rnorm(x)),
           multipliers = NULL,
           seed = NULL) {
    set.seed(seed)
    if (is.null(multipliers)) {
      multipliers = c(2, 1.5, 1.25)
    }
    if (length(multipliers) != 3) {
      stop('Needs 3 multipliers')
    }
    baseline_out <- out(N)
    potential_out <-
      rbind(
        multipliers[1] * baseline_out,
        multipliers[2] * baseline_out,
        multipliers[3] * baseline_out,
        baseline_out
      )
    rownames(potential_out) <-
      c('dir_indpsi', 'dir_indphi', 'indpsi', 'indphi')
    return(potential_out)
    
  }
