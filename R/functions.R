#' @import memoise
#' @import data.table
#' @import igraph
#' @import stringi

#' @export
make_tr_vec_permutation <- function(N,p,R,seed=NULL){
  set.seed(seed)
  tr_vec_sampled <- matrix(nrow=R,ncol=N)
  n_treated <- round(N*p)
  
  if (R > choose(N,n_treated)) {
    stop(paste("R must be smaller than", choose(N,n_treated),", the number of possible treatment assignements"))
  }
  
  for (i in 1:R) {
    vec <- sample(c(rep(1,n_treated),rep(0, N-n_treated)))
    while (any(duplicated(rbind(vec, tr_vec_sampled[1:i-1,]))))
    {
      vec <- sample(c(rep(1,n_treated),rep(0, N-n_treated)))
      
    }
    tr_vec_sampled[i,] <- vec
  }
  return(tr_vec_sampled)
}

#' @export
make_adj_matrix_sq_lattice <- function(N){
  if (sqrt(N) != round(sqrt(N))) {
    stop(paste('N must be a square number, not', N))
  }
  return(as.matrix(igraph::as_adj(igraph::graph.lattice(c(sqrt(N),sqrt(N)), circular = F))))
  
}

#' @export
make_adj_matrix_scale_free <- function(N, seed) {
  set.seed(seed)
  g <- igraph::barabasi.game(N, power = 1.2, m = NULL, out.dist = NULL, out.seq = NULL,
                     out.pref = FALSE, zero.appeal = 1, directed = FALSE,
                     algorithm ="psumtree", start.graph = NULL)
  while (min(igraph::degree(g, igraph::V(g)))==0) {
    g <- igraph::barabasi.game(N, power = 1.2, m = NULL, out.dist = NULL, out.seq = NULL,
                       out.pref = FALSE, zero.appeal = 1, directed = FALSE,
                       algorithm ="psumtree", start.graph = NULL)
  }
  return(as.matrix(igraph::as_adj(g)))
}

#' @export
make_adj_matrix_small_world <- function(N, seed) {
  set.seed(seed)
  g <- igraph::watts.strogatz.game(1, N, 2, 0.25, loops = FALSE, multiple = FALSE)
  while (min(igraph::degree(g, igraph::V(g)))==0) {
    g <- igraph::watts.strogatz.game(1, N, 2, 0.25, loops = FALSE, multiple = FALSE)
  }
  return(as.matrix(igraph::as_adj(g)))
}

#' @export
make_adj_matrix <- function(N, model, seed=NULL) {
  switch(model, 'sq_lattice'=return(make_adj_matrix_sq_lattice(N)),
         'scale_free'=return(make_adj_matrix_sq_lattice(N, seed)),
         'small_world'=return(make_adj_matrix_small_world(N, seed)))
}

#' @export
make_exposure_map_AS <- function(adj_matrix, tr_vector, hop) {
  N <- nrow(adj_matrix)
  peer_exposure <- as.numeric(tr_vector%*%adj_matrix)
  if (hop==1) {
    return(matrix(as.numeric(c(tr_vector>0 & peer_exposure>0,
                               tr_vector>0 & peer_exposure==0,
                               tr_vector==0 & peer_exposure>0,
                               tr_vector==0 & peer_exposure==0)),
                  N, 4, dimnames = list(NULL, c('dir_ind1', 'isol_dir', 'ind1','no'))))
  }
  if (hop==2) {
    adj_matrix_sq <- adj_matrix%*%adj_matrix # number of length two paths from i to j
    diag(adj_matrix_sq) <- rep(0, nrow(adj_matrix_sq))
    adj_matrix_2 <- adj_matrix_sq
    
    peer_exposure_2 <- as.numeric(tr_vector%*%adj_matrix_2)
    return(matrix(as.numeric(c(tr_vector>0 & peer_exposure>0 & peer_exposure_2>0,
                               tr_vector>0 & peer_exposure>0 & peer_exposure_2==0,
                               tr_vector>0 & peer_exposure==0 & peer_exposure_2>0,
                               tr_vector>0 & peer_exposure==0 & peer_exposure_2==0,
                               tr_vector==0 & peer_exposure>0 & peer_exposure_2>0,
                               tr_vector==0 & peer_exposure>0 & peer_exposure_2==0,
                               tr_vector==0 & peer_exposure==0 & peer_exposure_2>0,
                               tr_vector==0 & peer_exposure==0 & peer_exposure_2==0)),
                  N, 8, dimnames = list(NULL, c('dir_ind1_ind2', 'dir_ind1', 'dir_ind2', 'isol_dir', 'ind1_ind2', 'ind1', 'ind2', 'no'))))
  }
}

#' @export
make_corr_out <- function(degree, correlate, seed=NULL) {
  set.seed(seed) 
  switch(correlate, 'yes'= return(degree*abs(rnorm(length(degree))) + rnorm(length(degree),1,0.25)),
         'no' = return(abs(rnorm(length(degree)))))
}

#' @export
make_dilated_out_1 <- function(adj_matrix,make_corr_out,multipliers=NULL,seed=NULL) {
  set.seed(seed)  
  if (is.null(multipliers)) {
    multipliers=c(2,1.5,1.25)
  }
  if (length(multipliers)!=3) {
    stop('Needs 3 multipliers')
  }
  degree <- rowSums(adj_matrix)
  baseline_out <- make_corr_out(degree, 'yes', seed=seed)
  potential_out <- rbind(multipliers[1]*baseline_out, multipliers[2]*baseline_out,
                         multipliers[3]*baseline_out, baseline_out)
  rownames(potential_out) <- c('dir_ind1', 'isol_dir', 'ind1','no')
  return(potential_out)
  
}

#' @export
make_dilated_out_2 <- function(adj_matrix,make_corr_out,
                               multipliers=NULL,seed=NULL) {
  set.seed(seed)
  if (is.null(multipliers)) {
    multipliers=c(2.25,2,1.75,1.5,1.375,1.25,1.125)
  }
  if (length(multipliers)!=7) {
    stop('Needs 7 multipliers')
  }
  degree <- rowSums(adj_matrix)
  baseline_out <- make_corr_out(degree, 'yes', seed=seed)
  potential_out <- rbind(multipliers[1]*baseline_out, multipliers[2]*baseline_out,
                         multipliers[3]*baseline_out, multipliers[4]*baseline_out,
                         multipliers[5]*baseline_out, multipliers[6]*baseline_out,
                         multipliers[7]*baseline_out, baseline_out)
  rownames(potential_out) <- c('dir_ind1_ind2', 'dir_ind1', 'dir_ind2', 'isol_dir', 'ind1_ind2', 'ind1', 'ind2', 'no')
  return(potential_out)
  
} 

#' @export
make_dilated_out <- function(adj_matrix, make_corr_out, seed, hop, multipliers=NULL) {
  
  if (hop==1) {
    return(make_dilated_out_1(adj_matrix, make_corr_out,multipliers=multipliers, seed=seed))
  }
  if (hop==2) {
    return(make_dilated_out_2(adj_matrix, make_corr_out, multipliers=multipliers, seed=seed))
  }
}

#' @export
make_exposure_prob <- memoise(function(potential_tr_vector, adj_matrix, exposure_map_fn, exposure_map_fn_add_args=NULL) {
  exposure_map_fn_args <- c(list(adj_matrix, potential_tr_vector[1,]), exposure_map_fn_add_args)
  exposure_names <- colnames(do.call(exposure_map_fn, exposure_map_fn_args))
  n_exposure_conditions <- length(exposure_names)
  R <- nrow(potential_tr_vector)
  N <- ncol(potential_tr_vector)
  I_exposure <- list()
  for (i in 1:n_exposure_conditions) {
    I_exposure[[i]] <- matrix(nrow=N, ncol=R)
  }
  names(I_exposure) <- exposure_names
  
  
  for (i in 1:R) {
    exposure_map_fn_args <- c(list(adj_matrix, potential_tr_vector[i,]), exposure_map_fn_add_args)
    potential_exposure <- do.call(exposure_map_fn, exposure_map_fn_args)
    for (j in 1:ncol(potential_exposure)){
      I_exposure[[j]][,i] <- potential_exposure[,j]
    }
  }
  
  prob_exposure_k_k <- list()
  prob_exposure_k_l <- list()
  for (i in 1:length(I_exposure)) {
    for (j in 1:length(I_exposure)){
      prob_exposure_k_k[[paste(names(I_exposure)[[i]],names(I_exposure)[[i]],sep=',')]] <- (I_exposure[[i]]%*%t(I_exposure[[i]]) + diag(N))/(R+1)
      
      if (j!=i) {
        prob_exposure_k_l[[paste(names(I_exposure)[[i]],names(I_exposure)[[j]],sep=',')]] <- (I_exposure[[i]]%*%t(I_exposure[[j]]))/R
        
      }
      
    }
  }
  
  return(list(prob_exposure_k_k=prob_exposure_k_k, prob_exposure_k_l=prob_exposure_k_l))
  
})

#' @export
make_prob_exposure_cond <- function(prob_exposure) {
  
  k_exposure_names <- stri_split_fixed(str = names(prob_exposure$prob_exposure_k_k), pattern=',', simplify=T)[,1]
  
  prob_exposure_cond <- matrix(nrow = length(k_exposure_names), ncol = N)
  for (j in 1:length(prob_exposure$prob_exposure_k_k)) {
    prob_exposure_cond[j,] <- diag(prob_exposure$prob_exposure_k_k[[j]])
  }
  rownames(prob_exposure_cond) <- k_exposure_names
  
  return(prob_exposure_cond)
}

#' @export
var_yT_ht_unadjusted <- function(obs_exposure,obs_outcome,prob_exposure) {
  
  
  var_yT <- matrix(nrow = ncol(obs_exposure), dimnames=list(colnames(obs_exposure)))
  for (k in colnames(obs_exposure)) {
    pi_k <- prob_exposure$prob_exposure_k_k[[paste(k,k,sep=',')]]
    ind_kk <- diag(pi_k)
    cond_indicator = obs_exposure[,k]
    
    mm <- cond_indicator %o% cond_indicator * (pi_k - ind_kk %o% ind_kk)/pi_k * (obs_outcome %o% obs_outcome) / (ind_kk %o% ind_kk)
    mm[!is.finite(mm)] <- 0
    
    second_part_sum <- sum(mm)
    
    var_yT[k,] <- second_part_sum
  }
  return(var_yT)
}

#' @export
var_yT_ht_A2_adjustment <- function(obs_exposure,obs_outcome,prob_exposure) {
  
  var_yT_A2 <- matrix(nrow = ncol(obs_exposure), dimnames=list(colnames(obs_exposure)))
  
  for (k in colnames(obs_exposure)) {
    pi_k <- prob_exposure$prob_exposure_k_k[[paste(k,k,sep=',')]]
    ind_kk <- diag(pi_k)
    cond_indicator = obs_exposure[,k]
    
    m <- cond_indicator*(obs_outcome^2)/(2*ind_kk)
    A2_part_sum <- sum(outer(m, m, FUN='+') * (pi_k == 0) * (!diag(length(m))) )
    var_yT_A2[k,] <- A2_part_sum
  }
  return(var_yT_A2)
}

#' @export
var_yT_ht_adjusted <- function(obs_exposure,obs_outcome,prob_exposure) {
  var <- var_yT_ht_unadjusted(obs_exposure,obs_outcome,prob_exposure)
  A2 <- var_yT_ht_A2_adjustment(obs_exposure,obs_outcome,prob_exposure)
  var_adjusted <- var + A2
  return(var_adjusted)
}

#' @export
cov_yT_ht_adjusted <- function(obs_exposure,obs_outcome,prob_exposure, k_to_include=NULL, l_to_include=NULL) {
  
  cov_yT_A <- matrix(nrow = 2*choose(ncol(obs_exposure),2), dimnames=list(names(prob_exposure$prob_exposure_k_l)))
  
  if (is.null(k_to_include)) {
    k_to_include <- colnames(obs_exposure) }
  if (is.null(l_to_include)) {
    l_to_include <- colnames(obs_exposure) }
  
  for (k in k_to_include) {
    for (l in l_to_include) {
      if (k!=l) {
        
        pi_k <- prob_exposure$prob_exposure_k_k[[paste(k,k,sep=',')]]
        ind_kk <- diag(pi_k)
        pi_l <- prob_exposure$prob_exposure_k_k[[paste(l,l,sep=',')]]
        ind_ll <- diag(pi_l)
        pi_k_l <- prob_exposure$prob_exposure_k_l[[paste(k,l,sep=',')]]
        cond_indicator_k = obs_exposure[,k]
        cond_indicator_l = obs_exposure[,l]
        
        mm <- cond_indicator_k %o% cond_indicator_l * (pi_k_l - ind_kk %o% ind_ll)/pi_k_l * (obs_outcome %o% obs_outcome) / (ind_kk %o% ind_ll)
        mm[!is.finite(mm)] <- 0
        
        first_part_cov <- sum(mm)
        
        second_part_cov <- 0
        for (i in 1:length(cond_indicator_k)) {
          for (j in 1:length(cond_indicator_l)) {
            if (pi_k_l[i,j]==0) {
              second_part_cov_i_j <- ((cond_indicator_k[i]*obs_outcome[i]^2/(2*ind_kk[i])) + (cond_indicator_l[j]*obs_outcome[j]^2/(2*ind_ll[j])))
              second_part_cov <- second_part_cov + second_part_cov_i_j
              
            }
          }
        }
        cov_yT_A[paste(k,l,sep=',') , ] <- first_part_cov - second_part_cov
      }
      
    }
    
  }
  return(cov_yT_A)
}

#' @export
estimates <- memoise(function(obs_exposure, obs_outcome, obs_prob_exposure, hop) {
  
  obs_outcome_by_exposure <- t(obs_exposure)%*%diag(obs_outcome)
  obs_outcome_by_exposure[t(obs_exposure)==0] <- NA
  
  obs_prob_exposure_individual_kk <- make_prob_exposure_cond(obs_prob_exposure)
  
  yT_ht <- rowSums(obs_outcome_by_exposure/obs_prob_exposure_individual_kk, na.rm=T)
  yT_ht[rowSums(!is.na(obs_outcome_by_exposure)) == 0] <- NA
  
  mu_h <- yT_ht/rowSums(t(obs_exposure)/obs_prob_exposure_individual_kk, na.rm=T)
  mu_h[rowSums(!is.na(obs_outcome_by_exposure)) == 0] <- NA
  
  yT_h <- mu_h*N
  yT_h[rowSums(!is.na(obs_outcome_by_exposure)) == 0] <- NA
  
  resid_h <- colSums((obs_outcome_by_exposure-mu_h), na.rm=T) # pass resid_h instead of obs_outcome to var, cov
  
  var_yT_ht <- var_yT_ht_adjusted(obs_exposure,obs_outcome,obs_prob_exposure)
  
  var_yT_h <- var_yT_ht_adjusted(obs_exposure,resid_h,obs_prob_exposure)
  
  
  
  if (hop==1) { 
    
    cov_yT_ht <- cov_yT_ht_adjusted(obs_exposure,obs_outcome,obs_prob_exposure,
                                    k_to_include = c('dir_ind1', 'isol_dir', 'ind1'), l_to_include = 'no')
    cov_yT_h <- cov_yT_ht_adjusted(obs_exposure,resid_h,obs_prob_exposure,
                                   k_to_include = c('dir_ind1', 'isol_dir', 'ind1'), l_to_include = 'no' )
    
    remove <- 'no'
    keep <- c('dir_ind1,no', 'isol_dir,no', 'ind1,no')
    
    var_tau_ht <- (1/N^2)*(var_yT_ht[!rownames(var_yT_ht) %in% remove, ] + var_yT_ht[rownames(var_yT_ht) %in% remove, ] -
                             2*cov_yT_ht[rownames(cov_yT_ht) %in% keep, ])
    
    var_tau_h <- (1/N^2)*(var_yT_h[!rownames(var_yT_h) %in% remove, ] + var_yT_h[rownames(var_yT_h) %in% remove, ] -
                            2*cov_yT_h[rownames(cov_yT_h) %in% keep, ])
  }
  
  if (hop==2) { 
    cov_yT_ht <- cov_yT_ht_adjusted(obs_exposure,obs_outcome,obs_prob_exposure,
                                    k_to_include = c('dir_ind1_ind2', 'dir_ind1', 'dir_ind2', 'isol_dir',
                                                     'ind1_ind2', 'ind1', 'ind2'), l_to_include = 'no')
    cov_yT_h <- cov_yT_ht_adjusted(obs_exposure,resid_h,obs_prob_exposure,
                                   k_to_include =  c('dir_ind1_ind2', 'dir_ind1', 'dir_ind2', 'isol_dir',
                                                     'ind1_ind2', 'ind1', 'ind2'), l_to_include = 'no' )
    
    remove <- 'no'
    keep <- c('dir_ind1_ind2,no', 'dir_ind1,no', 'dir_ind2,no',
              'isol_dir,no', 'ind1_ind2,no', 'ind1,no', 'ind2,no')
    
    var_tau_ht <- (1/N^2)*(var_yT_ht[!rownames(var_yT_ht) %in% remove, ] + var_yT_ht[rownames(var_yT_ht) %in% remove, ] -
                             2*cov_yT_ht[rownames(cov_yT_ht) %in% keep, ])
    
    var_tau_h <- (1/N^2)*(var_yT_h[!rownames(var_yT_h) %in% remove, ] + var_yT_h[rownames(var_yT_h) %in% remove, ] -
                            2*cov_yT_h[rownames(cov_yT_h) %in% keep, ])
    
  }
  
  
  tau_ht <- (1/N)*(yT_ht-yT_ht['no'])[names(yT_ht)!='no']
  tau_h <- (mu_h-mu_h['no'])[names(mu_h)!='no']
  
  tau_dsm <- (rowMeans(obs_outcome_by_exposure, na.rm = T)-rowMeans(obs_outcome_by_exposure, na.rm = T)['no'])[names(yT_ht)!='no']
  
  
  return(list(yT_ht=yT_ht, yT_h=yT_h, tau_ht=tau_ht, tau_h=tau_h, tau_dsm=tau_dsm, var_tau_ht=var_tau_ht, var_tau_h=var_tau_h))
})



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