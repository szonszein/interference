#' @import data.table
#' @import memoise

#' @export
estimands_hierarchical <- function(potential_outcomes, group, k) {
  
  estimand <- cbind(as.data.frame(potential_outcomes['dir_indpsi', ]-potential_outcomes['indpsi', ]),
                    as.data.frame(potential_outcomes['dir_indphi', ]-potential_outcomes['indphi', ]),
                    as.data.frame(potential_outcomes['indpsi', ]-potential_outcomes['indphi', ]),
                    as.data.frame(potential_outcomes['dir_indpsi', ]-potential_outcomes['indphi', ]),
                    as.data.frame((potential_outcomes['dir_indpsi', ]*k[2]+potential_outcomes['indpsi', ]*k[1])-(potential_outcomes['dir_indphi', ]*k[1]+potential_outcomes['indphi', ]*k[2])),
                    group)
  colnames(estimand) <- NULL
  colnames(estimand) <- c('direct_psi', 'direct_phi', 'indirect', 'total', 'overall', 'group')
  estimand <- data.table(estimand)
  estimand <- estimand[, list(direct_psi_group=mean(direct_psi), direct_phi_group=mean(direct_phi), indirect_group=mean(indirect),
                           total_group=mean(total), overall_group=mean(overall)), by=group]
  
  estimand <- estimand[, list(direct_psi=mean(direct_psi_group), direct_phi=mean(direct_phi_group), indirect=mean(indirect_group),
                           total=mean(total_group), overall=mean(overall_group))]
  return(estimand)
}

#' @export
make_dilated_out_hh_miss <- function(N, out=function(x) abs(rnorm(x)), multipliers=NULL,seed=NULL) {
  set.seed(seed)  
  if (is.null(multipliers)) {
    multipliers=c(2.25,2,2,1.5,1.375,1.25,1.25)
  }
  if (length(multipliers)!=7) {
    stop('Needs 7 multipliers')
  }
  baseline_out <- out(N)
  potential_out <- rbind(multipliers[1]*baseline_out, multipliers[2]*baseline_out,
                         multipliers[3]*baseline_out, multipliers[4]*baseline_out,
                         multipliers[5]*baseline_out, multipliers[6]*baseline_out,
                         multipliers[7]*baseline_out, baseline_out)
  rownames(potential_out) <- c('dir_indpsipsi', 'dir_indpsiphi', 'dir_indphipsi', 'dir_indphiphi',
                               'indpsipsi','indpsiphi', 'indphipsi', 'indphiphi')
  return(potential_out)
  
}

#' @export
make_region_labels <- function(N,G,Re) {
  # Re is number of groups inside a region, not total number of regions
  if ((N%%G) !=0 | (G%%Re) != 0){
    stop('N must be divisible by G, G must be divisible by Re') 
  }
  
  rep(1:(G/Re), each=N/(G/Re))
}

#' @export
make_tr_vec_permutation_hierarchical_miss <- function(region,group,c,k,R,seed=NULL) {
  set.seed(seed)
  tr_vec_sampled <- vector('list', R)
  G <- length(unique(group))
  
  for (j in 1:R) {
    group_tr <- sample(c(rep(1,round(G*c)),rep(0, G-round(G*c))))
    names(group_tr) <- unique(group)
    
    indiv_tr <- rep(NA, length(group))
    for (i in unique(group)) {
      n <- sum(group==i)
      p <- k[group_tr[[i]] + 1]
      indiv_tr[which(group==i)] <- sample(c(rep(1,round(n*p)),rep(0, n-round(n*p))))
    }
    
    group_tr <- as.data.frame(group_tr)
    group_tr$group <- unique(group)
    
    group_tr_data <- merge(as.data.frame(cbind(region,group)), as.data.frame(group_tr), by='group')
    group_indiv_tr_data <- data.table(cbind(group_tr_data, as.data.frame(indiv_tr)))
    group_indiv_tr_data[,region_tr:=mean(group_tr),by=region]
    vec <- as.data.frame(group_indiv_tr_data)
    
    while(any(unlist(lapply(tr_vec_sampled, function(x) {!is.null(x) && all(vec==x)}))))
    {
      group_tr <- sample(c(rep(1,round(G*c)),rep(0, G-round(G*c))))
      names(group_tr) <- unique(group)
      
      indiv_tr <- rep(NA, length(group))
      for (i in unique(group)) {
        n <- sum(group==i)
        p <- k[group_tr[[i]] + 1]
        indiv_tr[which(group==i)] <- sample(c(rep(1,round(n*p)),rep(0, n-round(n*p))))
      }
      
      group_tr <- as.data.frame(group_tr)
      group_tr$group <- unique(group)
      
      group_tr_data <- merge(as.data.frame(cbind(region,group)), as.data.frame(group_tr), by='group')
      group_indiv_tr_data <- data.table(cbind(group_tr_data, as.data.frame(indiv_tr)))
      group_indiv_tr_data[,region_tr:=mean(group_tr),by=region]
      vec <- as.data.frame(group_indiv_tr_data)
      
    }
    tr_vec_sampled[[j]] <- vec
  }
  return(tr_vec_sampled)
}

#' @export
make_tr_condition_miss <- function(tr_assignment) {
  N <- nrow(tr_assignment)
  region_tr <- tr_assignment[, 'region_tr']
  group_tr <- tr_assignment[, 'group_tr']
  indiv_tr <- tr_assignment[, 'indiv_tr']
  return(matrix(as.numeric(c(region_tr==1 & group_tr>0 & indiv_tr>0,
                             region_tr==0.5 & group_tr>0 & indiv_tr>0,
                             region_tr==0.5 & group_tr==0 & indiv_tr>0,
                             region_tr==0 & group_tr==0 & indiv_tr>0,
                             region_tr==1 & group_tr>0 & indiv_tr==0,
                             region_tr==0.5 & group_tr>0 & indiv_tr==0,
                             region_tr==0.5 & group_tr==0 & indiv_tr==0,
                             region_tr==0 & group_tr==0 & indiv_tr==0)),
                N, 8, dimnames = list(NULL, c('dir_indpsipsi', 'dir_indpsiphi', 'dir_indphipsi', 'dir_indphiphi',
                                              'indpsipsi','indpsiphi', 'indphipsi', 'indphiphi'))))
}

#' @export
make_estimator_data_miss <- function(tr_assignment, potential_outcomes) {
  tr_condition <- make_tr_condition_miss(tr_assignment)
  obs_outcome <- rowSums(tr_condition*t(potential_outcomes))
  estimator_data <- cbind(tr_assignment, as.data.frame(obs_outcome)) 
  return(estimator_data)
}

#' @export
estimands_hierarchical_miss <- function(potential_outcomes, region, k) {
  # when groups and regions are of equal size, estimands from "miss" function and regular function should be same
  estimand <- cbind(as.data.frame(potential_outcomes['dir_indpsipsi', ]-potential_outcomes['indpsipsi', ]),
                    as.data.frame(potential_outcomes['dir_indphiphi', ]-potential_outcomes['indphiphi', ]),
                    as.data.frame(potential_outcomes['indpsipsi', ]-potential_outcomes['indphiphi', ]),
                    as.data.frame(potential_outcomes['dir_indpsipsi', ]-potential_outcomes['indphiphi', ]),
                    as.data.frame((potential_outcomes['dir_indpsipsi', ]*k[2]+potential_outcomes['indpsipsi', ]*k[1])-(potential_outcomes['dir_indphiphi', ]*k[1]+potential_outcomes['indphiphi', ]*k[2])),
                    region)
  colnames(estimand) <- NULL
  colnames(estimand) <- c('direct_psi', 'direct_phi', 'indirect', 'total', 'overall', 'region')
  estimand <- data.table(estimand)
  estimand <- estimand[, list(direct_psi_region=mean(direct_psi), direct_phi_region=mean(direct_phi), indirect_region=mean(indirect),
                              total_region=mean(total), overall_region=mean(overall)), by=region]
  
  estimand <- estimand[, list(direct_psi=mean(direct_psi_region), direct_phi=mean(direct_phi_region), indirect=mean(indirect_region),
                              total=mean(total_region), overall=mean(overall_region))]
  return(estimand)
}

