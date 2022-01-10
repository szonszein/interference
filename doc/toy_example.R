## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(interference)

N <- 10
p <- 0.2

RNGkind(sample.kind = "Rounding") # Required for compatiblity with R versions > 3.6.0.
Z <- make_tr_vec_permutation(N, p, R = 1, seed = 56)
Z

## -----------------------------------------------------------------------------
adj_matrix <- make_adj_matrix(N, model = 'small_world', seed = 492) 
adj_matrix

## -----------------------------------------------------------------------------
D <- make_exposure_map_AS(adj_matrix, Z, hop = 1)
D

## -----------------------------------------------------------------------------
omega <- make_tr_vec_permutation( N, p,
                                  R = 30, seed = 420, allow_repetitions = FALSE
)
prob_exposure <- make_exposure_prob(
  omega,
  adj_matrix, make_exposure_map_AS, list(hop = 1)
)
make_prob_exposure_cond(prob_exposure)

## -----------------------------------------------------------------------------
potential_outcomes <- make_dilated_out( adj_matrix, make_corr_out, seed = 1101,
                                        multipliers = NULL, hop = 1
)
observed_outcomes <- rowSums(D*t(potential_outcomes))

estimates(D, observed_outcomes, prob_exposure, control_condition='no')$tau_ht

