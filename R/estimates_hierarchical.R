#'Estimate average causal effects under a two-stage hierarchical design.
#'
#'Estimate average causal effects and their variance under two-stage
#'hierarchical treatment assignment, assuming partial and stratified
#'interference.
#'
#'\code{estimates} produces values of the estimators proposed by Hudgens and
#'Halloran (2008) of the population average direct causal effect under high
#'treatment saturation (psi), the population average direct causal effect under
#'low treatment saturation (phi), the population average indirect causal effect,
#'the population average total causal effect, and the population average overall
#'causal effect. It also produces values of the variance estimators which assume
#'stratified interference (i.e. potential outcomes of a unit are affected by its
#'own treatment assignment and only the treated proportion of its group; the
#'precise set of treated group members does not matter).
#'@param estimator_data data frame of `N` observations and the variables:
#'  \describe{ \item{`group`:}{Integer vector specifying the group label.}
#'  \item{`group_tr`:}{Numeric indicator of group assignment to high saturation
#'  (psi) in the first stage.} \item{`indiv_tr`:}{Numeric indicator of
#'  individual assignment to treatment in the second stage.}
#'  \item{`obs_outcome`:}{Numeric vector of observed outcomes.} }
#'@return A list of 2 lists: \enumerate{ \item A list of 5 scalars corresponding
#'  to the values of the estimators of the average direct causal effect under
#'  high treatment saturation (psi), the average direct causal effect under low
#'  treatment saturation (phi), the average indirect causal effect, the average
#'  total causal effect, and the average overall causal effect. \item A list of
#'  5 scalars corresponding to the values of the variance of the estimator of
#'  the average direct causal effect under high treatment saturation (psi), the
#'  average direct causal effect under low treatment saturation (phi), the
#'  average indirect causal effect, the average total causal effect, and the
#'  average overall causal effect. }
#' @examples
#' # Simulate first and second stage treatment
#' # assignment and outcomes:
#'
#' group <- rep(1:6, each = 30/6)
#' c <- 1/2
#' k <- c(2/5, 3/5)
#'
#' tr_assignment <-
#' make_tr_vec_permutation_hierarchical(group, c, k, R = 1,
#'                                      seed = 357)[[1]]
#'
#' potential_outcomes <- make_dilated_out_hh(N = 30,
#'                                           seed = 357)
#'
#' # Create data frame with group label, first and second stage
#' # treatment assignment and outcomes:
#'
#' estimator_data <- make_estimator_data(tr_assignment,
#'                                       potential_outcomes)
#'
#' # Estimate average causal effects and their variance:
#' estimates_hierarchical(estimator_data)
#'@export
#'@references Hudgens, M.G. & Halloran M.E. (2008). [Toward causal inference
#'  with interference](https://doi.org/10.1198/016214508000000292). *Journal of
#'  the American Statistical Association*, 103(482), 832--842.
#'
#'  Aronow, P.M. et al. (2020). [Spillover effects in experimental
#'  data](https://arxiv.org/abs/2001.05444). *arXiv preprint*, arXiv:2001.05444.
estimates_hierarchical <- function(estimator_data) {
  estimator_data <- data.table(estimator_data)
  
  # causal effects
  estimator_data_group_ov <-
    estimator_data[, list(y_hat_gr_ov = mean(obs_outcome)), by = c('group', 'group_tr')]
  estimator_data_pop_ov <-
    estimator_data_group_ov[, list(y_hat_ov = mean(y_hat_gr_ov)), by = c('group_tr')]
  
  estimator_data_group <-
    estimator_data[, list(y_hat_gr = mean(obs_outcome)), by = c('group', 'group_tr', 'indiv_tr')]
  estimator_data_pop <-
    estimator_data_group[, list(y_hat = mean(y_hat_gr)), by = c('group_tr', 'indiv_tr')]
  
  direct_psi <-
    estimator_data_pop[group_tr == 1 &
                         indiv_tr == 1, y_hat] - estimator_data_pop[group_tr == 1 &
                                                                      indiv_tr == 0, y_hat]
  direct_phi <-
    estimator_data_pop[group_tr == 0 &
                         indiv_tr == 1, y_hat] - estimator_data_pop[group_tr == 0 &
                                                                      indiv_tr == 0, y_hat]
  indirect <-
    estimator_data_pop[group_tr == 1 &
                         indiv_tr == 0, y_hat] - estimator_data_pop[group_tr == 0 &
                                                                      indiv_tr == 0, y_hat]
  total <-
    estimator_data_pop[group_tr == 1 &
                         indiv_tr == 1, y_hat] - estimator_data_pop[group_tr == 0 &
                                                                      indiv_tr == 0, y_hat]
  overall <-
    estimator_data_pop_ov[group_tr == 1, y_hat_ov] - estimator_data_pop_ov[group_tr ==
                                                                             0, y_hat_ov]
  
  # variance estimators
  sample_var_group <-
    estimator_data[, list(var_hat_gr = var(obs_outcome), K_gr = .N), by = c('group', 'group_tr', 'indiv_tr')] # sigma^2_i
  sample_var_pop <-
    estimator_data_group[, list(var_hat = var(y_hat_gr)), by = c('group_tr', 'indiv_tr')] # sigma^2_g
  sample_var_pop <-
    merge(sample_var_pop, estimator_data[, list(C = (.SD[unique(group), .N])), by =
                                           group_tr], by = 'group_tr')
  
  sample_var_pop_ov <-
    estimator_data_group_ov[, list(var_hat_ov = var(y_hat_gr_ov)), by = group_tr]
  sample_var_pop_ov <-
    merge(sample_var_pop_ov, estimator_data[, list(C = (.SD[unique(group), .N])), by =
                                              group_tr], by = 'group_tr')
  
  var_direct_group <-
    sample_var_group[, list(var_hat_direct_gr = sum(var_hat_gr / K_gr, na.rm =
                                                      TRUE)), by = c('group', 'group_tr')] # Var(CE_i^D)
  # 2nd component var hat CE
  second_var_direct <-
    var_direct_group[, list(sum_var_hat_direct = sum(var_hat_direct_gr),
                            C = .N), by = group_tr] # sum Var(CE_i^D)
  
  direct_group <-
    estimator_data_group[, list(direct_hat_gr = (.SD[indiv_tr == 1, y_hat_gr] -
                                                   .SD[indiv_tr == 0, y_hat_gr])), by = c('group', 'group_tr')]
  # 1st component var hat CE
  first_var_direct <-
    direct_group[, list(sample_var_hat_direct = var(direct_hat_gr)), by = group_tr] # sigma^2_D
  first_second_var_direct <-
    merge(first_var_direct, second_var_direct, by = 'group_tr')
  
  first <- first_second_var_direct[, sample_var_hat_direct / C]
  second <- first_second_var_direct[, sum_var_hat_direct / C]
  C <- first_second_var_direct[, C]
  G <-
    c(length(estimator_data[, unique(group)]), length(estimator_data[, unique(group)]))
  
  var_direct <-  (1 - C / G) * first + (1 / G) * second
  names(var_direct) <- first_second_var_direct$group_tr
  
  var_indirect <-
    sample_var_pop[group_tr == 1 &
                     indiv_tr == 0, var_hat / C] + sample_var_pop[group_tr == 0 &
                                                                    indiv_tr == 0, var_hat / C]
  var_total <-
    sample_var_pop[group_tr == 1 &
                     indiv_tr == 1, var_hat / C] + sample_var_pop[group_tr == 0 &
                                                                    indiv_tr == 0, var_hat / C]
  var_overall <-
    sample_var_pop_ov[group_tr == 1, var_hat_ov / C] + sample_var_pop_ov[group_tr ==
                                                                           0, var_hat_ov / C]
  
  
  return(list(
    list(
      direct_psi_hat = direct_psi,
      direct_phi_hat = direct_phi,
      indirect_hat = indirect,
      total_hat = total,
      overall_hat = overall
    ),
    list(
      var_direct_psi_hat = unname(var_direct['1']),
      var_direct_phi_hat = unname(var_direct['0']),
      var_indirect_hat = var_indirect,
      var_total_hat = var_total,
      var_overall_hat = var_overall
    )
  ))
  
}


#' @export
make_estimator_data <- function(tr_assignement, potential_outcomes) {
  tr_condition <- make_tr_condition(tr_assignement)
  obs_outcome <- rowSums(tr_condition * t(potential_outcomes))
  estimator_data <-
    cbind(tr_assignement, as.data.frame(obs_outcome))
  return(estimator_data)
}
