#' Create adjacency matrix.
#'
#' Undirected adjacency matrix of selected networks.
#' @param model string specifying the type of network. Options are
#'   `'sq_lattice'` for a circular lattice, `'scale_free'` for a scale free
#'   network with no isolates according to the Barabasi-Albert model,
#'   `'small_world'` for a network with no isolates according to the
#'   Watts-Strogatz small-world model, `'dcbm'` for a network according to the
#'   degree corrected stochastic block model.
#' @inheritParams make_tr_vec_permutation
#' @export
#' @examples
#' make_adj_matrix(N = 9, model = 'sq_lattice')
#' \dontrun{make_adj_matrix(N = 10, model = 'sq_lattice')}
#' make_adj_matrix(N = 10, model = 'small_world', seed = 357)
#' @return An `N` \eqn{*} `N` numeric matrix, where `N` corresponds to number of units.
#'   
make_adj_matrix <- function(N, model, seed = NULL) {
  switch(
    model,
    'sq_lattice' = return(make_adj_matrix_sq_lattice(N)),
    'scale_free' = return(make_adj_matrix_scale_free(N, seed)),
    'small_world' = return(make_adj_matrix_small_world(N, seed)),
    'dcbm' = return(make_adj_matrix_dcbm(N, seed))
  )
}



#' @rdname adj_matrix
#' @noRd
#' @export
make_adj_matrix_sq_lattice <- function(N) {
  if (sqrt(N) != round(sqrt(N))) {
    stop(paste('N must be a square number, not', N))
  }
  return(as.matrix(igraph::as_adj(
    igraph::graph.lattice(c(sqrt(N), sqrt(N)), circular = TRUE)
  )))
  
}


#' @rdname adj_matrix
#' @noRd
#' @export
make_adj_matrix_scale_free <- function(N, seed) {
  set.seed(seed)
  g <-
    igraph::barabasi.game(
      N,
      power = 0.6,
      m = 5,
      out.dist = NULL,
      out.seq = NULL,
      out.pref = FALSE,
      zero.appeal = 1,
      directed = FALSE,
      algorithm = "psumtree",
      start.graph = NULL
    )
  while (min(igraph::degree(g, igraph::V(g))) == 0) {
    g <-
      igraph::barabasi.game(
        N,
        power = 0.6,
        m = 5,
        out.dist = NULL,
        out.seq = NULL,
        out.pref = FALSE,
        zero.appeal = 1,
        directed = FALSE,
        algorithm = "psumtree",
        start.graph = NULL
      )
  }
  return(as.matrix(igraph::as_adj(g)))
}


#' @rdname adj_matrix
#' @noRd
#' @export
make_adj_matrix_small_world <- function(N, seed) {
  set.seed(seed)
  g <-
    igraph::watts.strogatz.game(1, N, 2, 0.25, loops = FALSE, multiple = FALSE)
  while (min(igraph::degree(g, igraph::V(g))) == 0) {
    g <-
      igraph::watts.strogatz.game(1, N, 2, 0.25, loops = FALSE, multiple = FALSE)
  }
  return(as.matrix(igraph::as_adj(g)))
}


#' @rdname adj_matrix
#' @noRd
#' @export
make_adj_matrix_dcbm <- function(N, seed) {
  set.seed(seed)
  dt <-
    randnet::BlockModel.Gen(
      4,
      N,
      K = 4,
      beta = (3 / 7),
      rho = 0.9,
      simple = FALSE,
      power = TRUE,
      alpha = 3.5
    )
  return(dt$A)
}

