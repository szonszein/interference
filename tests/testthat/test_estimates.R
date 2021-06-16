
test_that("estimates()", {
    obs_exposure_test <- matrix(c(1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0), ncol=4, nrow=3, dimnames=list(c(1, 2, 3), c('dir_ind1', 'isol_dir', 'ind1', 'no') ) )
    obs_outcome_test <- c(-1,2,3)
    prob_exposure_test <- list(prob_exposure_k_k=list(
                                    'dir_ind1,dir_ind1'=cbind(c(.3,.2,.1), c(.2,.5,.4), c(.1,.4,.8)),
                                    'isol_dir,isol_dir'=cbind(c(.3,.2,.1), c(.2,.5,.4), c(.1,.4,.8)), 
                                    'ind1,ind1'=cbind(c(.3,.2,.1), c(.2,.5,.4), c(.1,.4,.8)),
                                    'no,no'=cbind(c(.7,.4,.1), c(.4,.5,.1), c(.1,.1,.2))
                                    ),
                                prob_exposure_k_l=list(
                                    'dir_ind1,isol_dir'=cbind(c(0,.2,.4), c(.1,0,.3), c(0,.15,0)),
                                    'dir_ind1,ind1'=cbind(c(0,.2,.4), c(.1,0,.3), c(0,.15,0)),
                                    'dir_ind1,no'=cbind(c(0,.2,.4), c(.1,0,.3), c(0,.15,0)),

                                    'isol_dir,dir_ind1'=cbind(c(0,.2,.4), c(.1,0,.3), c(0,.15,0)),
                                    'isol_dir,ind1'=cbind(c(0,.2,.4), c(.1,0,.3), c(0,.15,0)),
                                    'isol_dir,no'=cbind(c(0,.2,.4), c(.1,0,.3), c(0,.15,0)),

                                    'ind1,dir_ind1'=cbind(c(0,.2,.4), c(.1,0,.3), c(0,.15,0)),
                                    'ind1,isol_dir'=cbind(c(0,.2,.4), c(.1,0,.3), c(0,.15,0)),
                                    'ind1,no'=cbind(c(0,.2,.4), c(.1,0,.3), c(0,.15,0)),

                                    'no,dir_ind1'=cbind(c(0,.1,0), c(.2,0,.15), c(.4,.3,0)),
                                    'no,isol_dir'=cbind(c(0,.1,0), c(.2,0,.15), c(.4,.3,0)),
                                    'no,ind1'=cbind(c(0,.1,0), c(.2,0,.15), c(.4,.3,0))
                                    ),
                                I_exposure=list(
                                    dir_ind1=matrix(c(1, 1, 0, 0, 0, 1), ncol=2, nrow=3 ),
                                    isol_dir=matrix(c(1, 1, 0, 0, 0, 1), ncol=2, nrow=3 ),
                                    ind1=matrix(c(1, 1, 0, 0, 0, 1), ncol=2, nrow=3 ),
                                    no=matrix(c(1, 1, 0, 0, 0, 1), ncol=2, nrow=3 )
                                )
    )
    est <- estimates(obs_exposure_test, obs_outcome_test, prob_exposure_test, n_var_permutations = 2, hop=1)
    expect_setequal(names(est), c(
      'yT_ht',
      'yT_h',
      'var_yT_ht',
      'var_yT_h',
      'cov_yT_ht',
      'var_tau_ht',
      'cov_yT_h',
      'var_tau_h',
      'var_tau_ht',
      'tau_ht',
      'tau_h',
      'tau_dsm'
      ))

      
    est <- estimates(obs_exposure_test, obs_outcome_test, prob_exposure_test, n_var_permutations = 2, hop=1, effect_estimators='hajek', variance_estimators='hajek')
    expect_setequal(names(est), c(
      'yT_h',
      'var_yT_h',
      'cov_yT_h',
      'var_tau_h',
      'tau_h',
      'tau_dsm'
      ))

    est <- estimates(obs_exposure_test, obs_outcome_test, prob_exposure_test, n_var_permutations = 2, hop=1, effect_estimators='horvitz-thompson', variance_estimators='horvitz-thompson')
    expect_setequal(names(est), c(
      'yT_ht',
      'var_yT_ht',
      'cov_yT_ht',
      'var_tau_ht',
      'var_tau_ht',
      'tau_ht',
      'tau_dsm'
      ))


})
