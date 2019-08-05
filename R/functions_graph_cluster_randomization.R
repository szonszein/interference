#' @import igraph


#' @export
make_tr_vec_bernoulli_indiv <- function(adj_matrix,p,R,seed=NULL){
  set.seed(seed)
  N <- nrow(adj_matrix)
  tr_vec_sampled <- matrix(nrow=R,ncol=N)
  
  if (R > 2^N) {
    stop(paste("R must be smaller than", 2^N,", the number of possible treatment assignements"))
  }
  
  for (i in 1:R) {
    vec <- rbinom(N,1,p)
    while (any(duplicated(rbind(vec, tr_vec_sampled[1:i-1,]))))
    {
      vec <- rbinom(N,1,p)
      
    }
    tr_vec_sampled[i,] <- vec
  }
  return(tr_vec_sampled)
}

#' @export
make_tr_vec_bernoulli_cluster <- function(adj_matrix,p,R,seed=NULL){
  set.seed(seed)
  N <- nrow(adj_matrix)
  tr_vec_sampled <- matrix(nrow=R,ncol=N)
  
  for (k in 1:R) {
  #3-net clustering
  G <- igraph::graph_from_adjacency_matrix(adjmatrix = adj_matrix, mode = 'undirected')
  
  B_2 <- igraph::ego(G, 2) # list of 2-hop neighbors
  
  vertices <- igraph::V(G)
  
  #if entry x of cluster_centers contains entry y, node x is the center of cluster 0 
  cluster_centers <- rep(NA, length(vertices))
  marked <- vector(length = length(vertices))
  j <- 0
  while(any(!marked)) {
    v_j <- sample(which(!marked), 1)
    cluster_centers[v_j] <- j
    marked[ v_j ]  <- TRUE
    marked[ unlist(B_2[v_j]) ]  <- TRUE
    j <- j +1
  }
  
  cluster_center_indices <- which(!is.na(cluster_centers ))
  distances_to_cluster_centers <- igraph::distances(G, to=cluster_center_indices)
  
  cluster_assignment <- vector(length = length(vertices))
  for (i in vertices) {
    #This breaks the tie, taking the first element
    cluster_index_i <- which(
      distances_to_cluster_centers[i,] == min(distances_to_cluster_centers[i,])
    )[1]
    cluster_assignment[i] <- cluster_centers[cluster_center_indices[cluster_index_i]]
  }
  
  # treatment assignement
  C <- length(unique(cluster_assignment))
  if (R > choose(N,C)*2^C) {
    stop(paste("R must be smaller than", choose(N,C)*2^C,", the number of possible treatment assignements"))
  }
  
  
    vec_cluster <- rbinom(C,1,p)
    cluster_assignment_index_vec_cluster <- cluster_assignment
    for (u in 1:length(unique(sort(cluster_assignment)))) {
      cluster_assignment_index_vec_cluster[which(cluster_assignment %in% unique(sort(cluster_assignment))[u])] <-u
    }
    vec <- vec_cluster[cluster_assignment_index_vec_cluster]
    
    while (any(duplicated(rbind(vec, tr_vec_sampled[1:k-1,]))))
    {
      cluster_centers <- rep(NA, length(vertices))
      marked <- vector(length = length(vertices))
      j <- 0
      while(any(!marked)) {
        v_j <- sample(which(!marked), 1)
        cluster_centers[v_j] <- j
        marked[ v_j ]  <- TRUE
        marked[ unlist(B_2[v_j]) ]  <- TRUE
        j <- j +1
      }
      
      cluster_center_indices <- which(!is.na(cluster_centers ))
      distances_to_cluster_centers <- igraph::distances(G, to=cluster_center_indices)
      
      cluster_assignment <- vector(length = length(vertices))
      for (i in vertices) {
        #This breaks the tie, taking the first element
        cluster_index_i <- which(
          distances_to_cluster_centers[i,] == min(distances_to_cluster_centers[i,])
        )[1]
        cluster_assignment[i] <- cluster_centers[cluster_center_indices[cluster_index_i]]
      }
      
      C <- length(unique(cluster_assignment))
      vec_cluster <- rbinom(C,1,p)
      cluster_assignment_index_vec_cluster <- cluster_assignment
      for (u in 1:length(unique(sort(cluster_assignment)))) {
        cluster_assignment_index_vec_cluster[which(cluster_assignment %in% unique(sort(cluster_assignment))[u])] <-u
      }
      
      vec <- vec_cluster[cluster_assignment_index_vec_cluster]
      
    }
    tr_vec_sampled[k,] <- vec
  }
  return(tr_vec_sampled)
}


#' @export
make_tr_vec_bernoulli <- function(adj_matrix,p,R,cluster,seed=NULL) {
  switch(cluster, 'no'=return(make_tr_vec_bernoulli_indiv(adj_matrix,p,R,seed)),
         'yes'=return(make_tr_vec_bernoulli_cluster(adj_matrix,p,R,seed)))
}

#' @export
make_exposure_map_full_neighborhood <- function(adj_matrix,tr_vector) {
N <- nrow(adj_matrix)
adj_matrix_diag <- adj_matrix
diag(adj_matrix_diag) <- rep(1, nrow(adj_matrix_diag))
exposure <- as.numeric(tr_vector%*%adj_matrix_diag)
B_size <- rowSums(adj_matrix_diag)
  return(matrix(as.numeric(c(exposure==B_size,
                             exposure<B_size)),
                N, 2, dimnames = list(NULL, c('all_treat', 'all_control'))))
}

#' @export
make_dilated_out_full_neighborhood <- function(adj_matrix,make_corr_out,multipliers=NULL,seed=NULL) {
  set.seed(seed)  
  if (is.null(multipliers)) {
    multipliers=c(2)
  }
  if (length(multipliers)!=1) {
    stop('Needs 1 multiplier')
  }
  degree <- rowSums(adj_matrix)
  nei2 <- adj_matrix%*%adj_matrix
  nei2[which(nei2>1)] <-1
  diag(nei2) <- rep(0, nrow(nei2))
  degree2 <- rowSums(nei2)
  
  baseline_out <- make_corr_out(degree, degree2, 'yes', seed=seed)
  potential_out <- rbind(multipliers[1]*baseline_out, baseline_out)
  rownames(potential_out) <- c('all_treat', 'all_control')
  return(potential_out)
  
}

#' @export
estimators_full_neighborhood <- memoise(function(obs_exposure, obs_outcome, obs_prob_exposure, n_var_permutations=10) {
  
  obs_outcome_by_exposure <- t(obs_exposure)%*%diag(obs_outcome)
  obs_outcome_by_exposure[t(obs_exposure)==0] <- NA
  
  obs_prob_exposure_individual_kk <- make_prob_exposure_cond(obs_prob_exposure)
  
  yT_ht <- rowSums(obs_outcome_by_exposure/obs_prob_exposure_individual_kk, na.rm=T)
  yT_ht[rowSums(!is.na(obs_outcome_by_exposure)) == 0] <- NA
  
  mu_h <- yT_ht/rowSums(t(obs_exposure)/obs_prob_exposure_individual_kk, na.rm=T)
  mu_h[rowSums(!is.na(obs_outcome_by_exposure)) == 0] <- NA
  
  yT_h <- mu_h*N
  yT_h[rowSums(!is.na(obs_outcome_by_exposure)) == 0] <- NA
  
  resid_h <- colSums((obs_outcome_by_exposure-mu_h), na.rm=T) # pass resid_h instead of obs_outcome to var, cov to get Hajek estimator
  
  var_yT_ht <- var_yT_ht_adjusted(obs_exposure,obs_outcome,obs_prob_exposure)
  
  var_yT_h <- var_yT_ht_adjusted(obs_exposure,resid_h,obs_prob_exposure)
    
  cov_yT_ht <- cov_yT_ht_adjusted(obs_exposure,obs_outcome,obs_prob_exposure)
  cov_yT_h <- cov_yT_ht_adjusted(obs_exposure,resid_h,obs_prob_exposure)
  
  const_eff <- var_yT_ht_const_eff_lm(obs_exposure,obs_outcome,obs_prob_exposure, n_var_permutations)
  
  var_yT_ht_const_eff <- const_eff$var_yT_ht_const_eff
  
  diffs <- const_eff$means - const_eff$means$all_control
  diffs <- subset(diffs, select=-c(all_control))
  
  var_tau_ht_const_eff <- unlist(lapply(diffs, var))
    
    remove <- 'all_control'
    keep <- c('all_treat,all_control')
    
    var_tau_ht <- (1/N^2)*(var_yT_ht[!rownames(var_yT_ht) %in% remove, ] + var_yT_ht[rownames(var_yT_ht) %in% remove, ] -
                             2*cov_yT_ht[rownames(cov_yT_ht) %in% keep, ])
    
    var_tau_h <- (1/N^2)*(var_yT_h[!rownames(var_yT_h) %in% remove, ] + var_yT_h[rownames(var_yT_h) %in% remove, ] -
                            2*cov_yT_h[rownames(cov_yT_h) %in% keep, ])

  var_tau_ht_max <- pmax(var_tau_ht[order(names(var_tau_ht))], var_tau_ht_const_eff[order(names(var_tau_ht_const_eff))])
    
  tau_ht <- (1/N)*(yT_ht-yT_ht['all_control'])[names(yT_ht)!='all_control']
  tau_h <- (mu_h-mu_h['all_control'])[names(mu_h)!='all_control']
  
  return(list(yT_ht=yT_ht, yT_h=yT_h, tau_ht=tau_ht, tau_h=tau_h, var_tau_ht=var_tau_ht, var_tau_h=var_tau_h,
              var_tau_ht_const_eff=var_tau_ht_const_eff, var_tau_ht_max=var_tau_ht_max))
})

