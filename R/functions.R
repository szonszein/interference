#' @import memoise
#' @import data.table
#' @import igraph
#' @import stringi
#' @import randnet
#' @import combinat


#' @export
naive_estimate <- memoise(function(obs_outcome_by_exposure, hop) {
  if (hop==1) { 
    yT_treat <- sum((1/p)*colSums(obs_outcome_by_exposure[c('dir_ind1', 'isol_dir'),], na.rm = T))
    yT_control <- sum((1/(1-p))*colSums(obs_outcome_by_exposure[c('ind1', 'no'),], na.rm = T))
    estimate <- (1/N)*(yT_treat-yT_control)
    return(list(estimate=estimate, yT_treat=yT_treat, yT_control=yT_control))
  }
  if (hop==2) {
    yT_treat <- sum((1/p)*colSums(obs_outcome_by_exposure[c('dir_ind1_ind2', 'dir_ind1', 'dir_ind2', 'isol_dir'),], na.rm = T))
    yT_control <- sum((1/(1-p))*colSums(obs_outcome_by_exposure[c('ind1_ind2', 'ind1', 'ind2', 'no'),], na.rm = T))
    estimate <- (1/N)*(yT_treat-yT_control)
    return(list(estimate=estimate, yT_treat=yT_treat, yT_control=yT_control))
  }  
})

#' @export
make_adj_matrix_trunc <-memoise(function(adj_matrix,p,seed=NULL) {
  set.seed(seed)
  m <- adj_matrix
  for (i in 1:nrow(m)) {
    for (j in 1:ncol(m)) {
      if (j>i & m[i,j]==1) {
        m[i,j] <- ifelse(runif(1) < p,0,1)
      }
    }
  }
  m[lower.tri(m,diag=FALSE)] <- m[upper.tri(m,diag=FALSE)]
  return(m)
})

#' @export
make_adj_matrix_add <-memoise(function(adj_matrix,p,seed=NULL) {
  #  p, as passed, should be the approximate proportion of the existing edges to add
  
  #  Rescale p to be a proportion of the not-currently-existing edges
  p <- p * sum(adj_matrix) / (ncol(adj_matrix) * (ncol(adj_matrix) -1) - sum(adj_matrix))
  set.seed(seed)
  m <- adj_matrix
  for (i in 1:nrow(m)) {
    for (j in 1:ncol(m)) {
      if (j>i & m[i,j]==0) {
        m[i,j] <- ifelse(runif(1) < p,1,0)
      }
    }
  }
  m[lower.tri(m,diag=FALSE)] <- m[upper.tri(m,diag=FALSE)]
  return(m)
})

#' @export
make_adj_matrix_miss_ties <- function(adj_matrix,p,type,seed=NULL) {
  if (p==0 & type!='nothing') {
    stop('If p=0, type must be nothing')
  }
  if (p!=0 & type=='nothing') {
    stop('If type is nothing, p must be 0')
  }
  switch(type, 'trunc'=return(make_adj_matrix_trunc(adj_matrix,p,seed)),
         'nothing'=return(adj_matrix),
         'add'=return(make_adj_matrix_add(adj_matrix,p,seed)))
}