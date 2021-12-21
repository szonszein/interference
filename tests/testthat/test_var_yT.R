
test_that("var(yT) part of var_yT", {
    obs_exposure_test <- cbind(c(1,1,0), c(0,0,1))
    colnames(obs_exposure_test) <- c('cond1', 'cond2')
    obs_outcome_test <- c(-1,2,3)
    prob_exposure_test <- list(prob_exposure_k_k=list('cond1,cond1'=cbind(c(.3,.2,.1), c(.2,.5,.4), c(.1,.4,.8)),
                                                    'cond2,cond2'=cbind(c(.7,.4,.1), c(.4,.5,.1), c(.1,.1,.2))))

    expected <- matrix(c(9.11111,180), nrow = 2, ncol = 1, dimnames = list(c('cond1', 'cond2')))

    result <- var_yT_ht_unadjusted(obs_exposure_test, obs_outcome_test, prob_exposure_test)
    expect_true(all.equal(expected,result, tolerance=1e-5))
})


test_that("A2 part of var_yT", {
    obs_exposure_test <- cbind(c(1,1,0), c(0,0,1))
    colnames(obs_exposure_test) <- c('cond1', 'cond2')
    obs_outcome_test <- c(-1,2,3)
    prob_exposure_test <- list(prob_exposure_k_k=list('cond1,cond1'=cbind(c(.3,0,0), c(0,.5,.4), c(0,.4,.8)),
                                                    'cond2,cond2'=cbind(c(.7,.4,.1), c(.4,.5,0), c(.1,0,.2))))

    expected <- matrix(c(14.66667,45), nrow = 2, ncol = 1, dimnames = list(c('cond1', 'cond2')))

    result <- var_yT_ht_A2_adjustment(obs_exposure_test, obs_outcome_test, prob_exposure_test)
    expect_true(all.equal(expected,result, tolerance=1e-5))
})

test_that("# CovA(yT_k, yT_l) part of var_yT", {
    obs_exposure_test <- cbind(c(1,1,0), c(0,0,1))
    colnames(obs_exposure_test) <- c('cond1', 'cond2')
    obs_outcome_test <- c(-1,2,3)
    prob_exposure_test <- list(prob_exposure_k_k=list('cond1,cond1'=cbind(c(.3,.2,.1), c(.2,.5,.4), c(.1,.4,.8)),
                                                    'cond2,cond2'=cbind(c(.7,.4,.1), c(.4,.5,.1), c(.1,.1,.2))),
                            prob_exposure_k_l=list('cond1,cond2'=cbind(c(0,.2,.4), c(.1,0,.3), c(0,.15,0)),
                                                    'cond2,cond1'=cbind(c(0,.1,0), c(.2,0,.15), c(.4,.3,0))))

    expected <- matrix(c(-32.33333,-32.33333), nrow = 2, ncol = 1, dimnames = list(c('cond1,cond2', 'cond2,cond1')))
    result <- cov_yT_ht_adjusted(obs_exposure_test, obs_outcome_test, prob_exposure_test)
    expect_true(all.equal(expected,result, tolerance=1e-5))

})