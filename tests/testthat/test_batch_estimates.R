hop <- 1
n_var_permutations <- 30
exposure_map_fn <- make_exposure_map_AS
exposure_map_fn_add_args <- list(hop = 1)

adj_matrix <- make_adj_matrix(N = 9, model = 'sq_lattice')
tr_vector <-
  make_tr_vec_permutation(N = 9,
                          p = 0.2,
                          R = 1,
                          seed = 357)
obs_exposure <-
  make_exposure_map_AS(adj_matrix, tr_vector, hop = 1)

# Simulate a vector of outcome data:
potential_outcome <-
  make_dilated_out(adj_matrix, make_corr_out, seed = 357, hop = 1)
obs_outcome <- rowSums(obs_exposure * t(potential_outcome))

# Create exposure probabilities:
potential_tr_vector <-
  make_tr_vec_permutation(N = 9,
                          p = 0.2,
                          R = 36,
                          seed = 357)
obs_prob_exposure <-
  make_exposure_prob(
    potential_tr_vector,
    adj_matrix,
    exposure_map_fn = exposure_map_fn,
    exposure_map_fn_add_args = exposure_map_fn_add_args
  )


non_batch_estimates <-
  suppressWarnings(
    estimates(
      obs_exposure,
      obs_outcome,
      obs_prob_exposure,
      n_var_permutations = n_var_permutations,
      hop = hop
    )
  )

test_that("batch estimates produce same results as non-batch, with batch_size=N", {

  batch_estimates <- estimates_batch(
    obs_exposure,
    obs_outcome,
    obs_prob_exposure,
    potential_tr_vector,
    adj_matrix,
    exposure_map_fn = exposure_map_fn,
    exposure_map_fn_add_args = exposure_map_fn_add_args,
    n_var_permutations = n_var_permutations,
  
  )
  
 
  expect_equal(batch_estimates$yT_ht, non_batch_estimates$yT_ht)
  expect_equal(batch_estimates$var_yT_ht, non_batch_estimates$var_yT_ht)
  expect_equal(batch_estimates$cov_yT_ht, non_batch_estimates$cov_yT_ht)
               
})

test_that("batch estimates produce same results as non-batch, with batch_size=1", {
  
  batch_estimates <- estimates_batch(
    obs_exposure,
    obs_outcome,
    obs_prob_exposure,
    potential_tr_vector,
    adj_matrix,
    exposure_map_fn = exposure_map_fn,
    exposure_map_fn_add_args = exposure_map_fn_add_args,
    n_var_permutations = n_var_permutations,
    batch_size=1
  )
  
  expect_equal(batch_estimates$yT_ht, non_batch_estimates$yT_ht)
  expect_equal(batch_estimates$var_yT_ht, non_batch_estimates$var_yT_ht)
  expect_equal(batch_estimates$cov_yT_ht, non_batch_estimates$cov_yT_ht)
  
})

test_that("batch estimates produce same results as non-batch, with batch_size=3 (divisible by N)", {
  
  batch_estimates <- estimates_batch(
    obs_exposure,
    obs_outcome,
    obs_prob_exposure,
    potential_tr_vector,
    adj_matrix,
    exposure_map_fn = exposure_map_fn,
    exposure_map_fn_add_args = exposure_map_fn_add_args,
    n_var_permutations = n_var_permutations,
    batch_size=3
    
  )
  
  
  expect_equal(batch_estimates$yT_ht, non_batch_estimates$yT_ht)
  expect_equal(batch_estimates$var_yT_ht, non_batch_estimates$var_yT_ht)
  expect_equal(batch_estimates$cov_yT_ht, non_batch_estimates$cov_yT_ht)
  
})

test_that("batch estimates produce same results as non-batch, with batch_size=4 (non-divisible by N)", {
  
  batch_estimates <- estimates_batch(
    obs_exposure,
    obs_outcome,
    obs_prob_exposure,
    potential_tr_vector,
    adj_matrix,
    exposure_map_fn = exposure_map_fn,
    exposure_map_fn_add_args = exposure_map_fn_add_args,
    n_var_permutations = n_var_permutations,
    batch_size=4
    
  )
  
  
  expect_equal(batch_estimates$yT_ht, non_batch_estimates$yT_ht)
  expect_equal(batch_estimates$var_yT_ht, non_batch_estimates$var_yT_ht)
  expect_equal(batch_estimates$cov_yT_ht, non_batch_estimates$cov_yT_ht)
  
})

