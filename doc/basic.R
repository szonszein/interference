## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- include=FALSE-----------------------------------------------------------
#  This chunk creates the plausible but fake data for the vignette

library(igraph)
library(igraphdata)
library(interference)
data(USairports)
agraph <- as.undirected(USairports)
amat <- as_adj(agraph, sparse=FALSE)
igraph::write_graph(agraph, './airports_undirected.csv', format='edgelist')

outcomes <- make_dilated_out(amat, make_corr_out, seed=0, hop=1, multipliers=NULL)
tr_vector <- make_tr_vec_permutation(nrow(amat),0.2,R=1,seed=4224)

exposure <- make_exposure_map_AS(amat, tr_vector, hop=1) 
obs_outcome <- rowSums(exposure*t(outcomes))

vars <- data.frame(id=1:755, treatment=tr_vector[1,], y=obs_outcome)
fwrite(vars, './airports_vars_hop1.csv')

outcomes <- make_dilated_out(amat, make_corr_out, seed=0, hop=2, multipliers=NULL)
exposure <- make_exposure_map_AS(amat, tr_vector, hop=2) 
obs_outcome <- rowSums(exposure*t(outcomes))

vars <- data.frame(id=1:755, treatment=tr_vector[1,], y=obs_outcome)
fwrite(vars, './airports_vars_hop2.csv')


## -----------------------------------------------------------------------------
network <- read_graph('./airports_undirected.csv', format='edgelist', directed=FALSE)
adj_matrix <- as_adj(network, sparse=FALSE)

## -----------------------------------------------------------------------------
diag(adj_matrix) <- 0

## -----------------------------------------------------------------------------
isSymmetric(adj_matrix)

## -----------------------------------------------------------------------------
vars <- read.csv('./airports_vars_hop1.csv')

## -----------------------------------------------------------------------------
non_isolates <- which(rowSums(adj_matrix)>0)

print(non_isolates)
print(length(non_isolates))
print(nrow(adj_matrix))

adj_matrix <- adj_matrix[non_isolates, non_isolates]
vars <- vars[vars$id %in% non_isolates,]


## -----------------------------------------------------------------------------
num_replicates <- 1500
prop_treated <- 0.2
N <- nrow(adj_matrix)

potential_tr_vector <- make_tr_vec_permutation(N, p=prop_treated, R=num_replicates, seed=4224)
exposure <- make_exposure_map_AS(adj_matrix, vars$treatment, hop=1) 
obs_prob_exposure <- make_exposure_prob(potential_tr_vector, adj_matrix, make_exposure_map_AS, list(hop=1))

## ---- message=FALSE, warning=FALSE, cache=FALSE-------------------------------
my_estimates <- estimates(exposure, vars$y, obs_prob_exposure, n_var_permutations = 1000, control_condition='no') 

## -----------------------------------------------------------------------------
my_estimates$tau_ht
my_estimates$var_tau_ht
my_estimates$var_tau_ht_const_eff

## -----------------------------------------------------------------------------
my_estimates$tau_h
my_estimates$var_tau_h

## -----------------------------------------------------------------------------
vars <- read.csv('./airports_vars_hop2.csv')
vars <- vars[vars$id %in% non_isolates,]

## -----------------------------------------------------------------------------
exposure <- make_exposure_map_AS(adj_matrix, vars$treatment, hop=2) 
obs_prob_exposure <- make_exposure_prob(potential_tr_vector, adj_matrix, make_exposure_map_AS, list(hop=2))

## ---- message=FALSE, warning=FALSE--------------------------------------------
my_estimates <- estimates(exposure, vars$y, obs_prob_exposure, n_var_permutations = 1000, control_condition='no') 

## -----------------------------------------------------------------------------
my_estimates$tau_ht
my_estimates$var_tau_ht
my_estimates$var_tau_ht_const_eff

## -----------------------------------------------------------------------------
my_estimates$tau_h
my_estimates$var_tau_h

