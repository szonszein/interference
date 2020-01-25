#' Assign treatment via 3-net clustering randomization.
#'
#' Create a permutation matrix from possible treatment assignments via 3-net
#' clustering randomization or via unit-level randomization under a Bernoulli
#' distribution.
#'
#' \code{make_tr_vec_bernoulli} produces all possible (or a random sample of all
#' possible) treatment assignments via 3-net clustering randomization following
#' the algorithm of Ugander et al. (2013), or via unit-level randomization under
#' a Bernoulli(p) distribution. Sampling of treatment vectors is without
#' replacement.
#'
#' @param p probability that treatment takes the value 1.
#' @param R number of repetitions (treatment assignments). `R` must be smaller
#'   or equal to the number of possible treatment assignements which is
#'   \eqn{combination(N, C) * 2 ^ C} when \code{cluster = 'yes'}, and \eqn{2 ^
#'   N} when \code{cluster = 'no'}, where N corresponds to the number of units
#'   and C to the number of clusters.
#' @param cluster string; either `'yes'` or `'no'`. If `'yes'` units are
#'   assigned to treatment via 3-net clustering randomization, if `'no'` via
#'   unit-level randomization under a Bernoulli(p) distribution.
#' @inheritParams make_tr_vec_permutation
#' @inheritParams make_exposure_map_AS
#' @references Ugander, J. et al. (2013). [Graph cluster randomization: Network
#'   exposure to multiple universes](https://doi.org/10.1145/2487575.2487695).
#'   *Proceedings of the 19th ACM SIGKDD international conference on Knowledge
#'   discovery and data mining*. 329--337.
#' @export
#' @examples
#' # Create adjacency matrix:
#'
#' adj_matrix <- make_adj_matrix(N = 50, model = 'small_world',
#'                               seed = 357)
#'
#' # Assign treatment via 3-net clustering randomization:
#' make_tr_vec_bernoulli(adj_matrix, p = 0.2, R = 1,
#'                       cluster = 'yes', seed = 357)
#'
#' # Assign treatment via unit-level randomization
#' # under Bernoulli(p) distribution:
#'
#' make_tr_vec_bernoulli(adj_matrix, p = 0.2, R = 1,
#'                             cluster = 'no', seed = 357)
#'
#' @return An `R` \eqn{*} `N` numeric matrix. Each row cooresponds to a
#'   treatment assignment vector.
make_tr_vec_bernoulli <-
  function(adj_matrix, p, R, cluster, seed = NULL) {
    switch(cluster, 'no' = return(make_tr_vec_bernoulli_indiv(adj_matrix, p, R, seed)),
           'yes' = return(make_tr_vec_bernoulli_cluster(adj_matrix, p, R, seed)))
  }


#' @rdname tr_vec_bernoulli
#' @noRd  
#' @export  
make_tr_vec_bernoulli_cluster <-
  function(adj_matrix, p , R, seed = NULL) {
    set.seed(seed)
    N <- nrow(adj_matrix)
    tr_vec_sampled <- matrix(nrow = R, ncol = N)
    
    for (k in 1:R) {
      #3-net clustering
      G <-
        igraph::graph_from_adjacency_matrix(adjmatrix = adj_matrix, mode = 'undirected')
      
      B_2 <- igraph::ego(G, 2) # list of 2-hop neighbors
      
      vertices <- igraph::V(G)
      
      #if entry x of cluster_centers contains entry y, node x is the center of cluster 0
      cluster_centers <- rep(NA, length(vertices))
      marked <- vector(length = length(vertices))
      j <- 0
      while (any(!marked)) {
        v_j <- sample(which(!marked), 1)
        cluster_centers[v_j] <- j
        marked[v_j]  <- TRUE
        marked[unlist(B_2[v_j])]  <- TRUE
        j <- j + 1
      }
      
      cluster_center_indices <- which(!is.na(cluster_centers))
      distances_to_cluster_centers <-
        igraph::distances(G, to = cluster_center_indices)
      
      cluster_assignment <- vector(length = length(vertices))
      for (i in vertices) {
        #This breaks the tie, taking the first element
        cluster_index_i <- which(distances_to_cluster_centers[i, ] == min(distances_to_cluster_centers[i, ]))[1]
        cluster_assignment[i] <-
          cluster_centers[cluster_center_indices[cluster_index_i]]
      }
      
      # treatment assignement
      C <- length(unique(cluster_assignment))
      if (R > choose(N, C) * 2 ^ C) {
        stop(
          paste(
            "R must be smaller than",
            choose(N, C) * 2 ^ C,
            ", the number of possible treatment assignements"
          )
        )
      }
      
      
      vec_cluster <- rbinom(C, 1, p)
      cluster_assignment_index_vec_cluster <- cluster_assignment
      for (u in 1:length(unique(sort(cluster_assignment)))) {
        cluster_assignment_index_vec_cluster[which(cluster_assignment %in% unique(sort(cluster_assignment))[u])] <-
          u
      }
      vec <- vec_cluster[cluster_assignment_index_vec_cluster]
      
      while (any(duplicated(rbind(vec, tr_vec_sampled[1:k - 1, ]))))
      {
        cluster_centers <- rep(NA, length(vertices))
        marked <- vector(length = length(vertices))
        j <- 0
        while (any(!marked)) {
          v_j <- sample(which(!marked), 1)
          cluster_centers[v_j] <- j
          marked[v_j]  <- TRUE
          marked[unlist(B_2[v_j])]  <- TRUE
          j <- j + 1
        }
        
        cluster_center_indices <- which(!is.na(cluster_centers))
        distances_to_cluster_centers <-
          igraph::distances(G, to = cluster_center_indices)
        
        cluster_assignment <- vector(length = length(vertices))
        for (i in vertices) {
          #This breaks the tie, taking the first element
          cluster_index_i <- which(distances_to_cluster_centers[i, ] == min(distances_to_cluster_centers[i, ]))[1]
          cluster_assignment[i] <-
            cluster_centers[cluster_center_indices[cluster_index_i]]
        }
        
        C <- length(unique(cluster_assignment))
        vec_cluster <- rbinom(C, 1, p)
        cluster_assignment_index_vec_cluster <- cluster_assignment
        for (u in 1:length(unique(sort(cluster_assignment)))) {
          cluster_assignment_index_vec_cluster[which(cluster_assignment %in% unique(sort(cluster_assignment))[u])] <-
            u
        }
        
        vec <- vec_cluster[cluster_assignment_index_vec_cluster]
        
      }
      tr_vec_sampled[k, ] <- vec
    }
    return(tr_vec_sampled)
  }



#' @rdname tr_vec_bernoulli
#' @noRd
#' @export
make_tr_vec_bernoulli_indiv <-
  function(adj_matrix, p, R, seed = NULL) {
    set.seed(seed)
    N <- nrow(adj_matrix)
    tr_vec_sampled <- matrix(nrow = R, ncol = N)
    
    if (R > 2 ^ N) {
      stop(
        paste(
          "R must be smaller than",
          2 ^ N,
          ", the number of possible treatment assignements"
        )
      )
    }
    
    for (i in 1:R) {
      vec <- rbinom(N, 1, p)
      while (any(duplicated(rbind(vec, tr_vec_sampled[1:i - 1,]))))
      {
        vec <- rbinom(N, 1, p)
        
      }
      tr_vec_sampled[i,] <- vec
    }
    return(tr_vec_sampled)
  }