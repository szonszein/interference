#' @import data.table
#' @import memoise

#' @export
make_dilated_out_hh <- function(N, out=function(x) abs(rnorm(x)), multipliers=NULL,seed=NULL) {
  set.seed(seed)  
  if (is.null(multipliers)) {
    multipliers=c(2,1.5,1.25)
  }
  if (length(multipliers)!=3) {
    stop('Needs 3 multipliers')
  }
  baseline_out <- out(N)
  potential_out <- rbind(multipliers[1]*baseline_out, multipliers[2]*baseline_out,
                         multipliers[3]*baseline_out, baseline_out)
  rownames(potential_out) <- c('dir_indpsi', 'dir_indphi', 'indpsi','indphi')
  return(potential_out)
  
}


#' @export
make_group_labels <- function(N,G) {
  if ((N%%G) != 0 ){
    stop('N must be divisible by G') 
  }
  
  rep(1:G, each=N/G)
}


#' @export
make_tr_vec_permutation_hierarchical <- function(group,c,k,R,seed=NULL) {
  set.seed(seed)
  tr_vec_sampled <- vector('list',R)
  G <- length(unique(group))
  
  for (j in 1:R) {
    group_tr <- as.list(sample(c(rep(1,round(G*c)),rep(0, G-round(G*c)))))
    names(group_tr) <- unique(group)
    
    indiv_tr <- rep(NA, length(group))
    for (i in unique(group)) {
      n <- sum(group==i)
      p <- k[group_tr[[i]] + 1]
      indiv_tr[which(group==i)] <- sample(c(rep(1,round(n*p)),rep(0, n-round(n*p))))
    }
    
    group_tr <- as.data.frame(unlist(group_tr))
    colnames(group_tr) <- 'group_tr'
    group_tr$group <- unique(group)
    
    group_tr_data <- merge(as.data.frame(group), as.data.frame(group_tr), by='group')
    vec <- cbind(group_tr_data, as.data.frame(indiv_tr))
    
    while(any(unlist(lapply(tr_vec_sampled, function(x) {!is.null(x) && all(vec==x)}))))
    {
      group_tr <- as.list(sample(c(rep(1,round(G*c)),rep(0, G-round(G*c)))))
      names(group_tr) <- unique(group)
      
      indiv_tr <- rep(NA, length(group))
      for (i in unique(group)) {
        n <- sum(group==i)
        p <- k[group_tr[[i]] + 1]
        indiv_tr[which(group==i)] <- sample(c(rep(1,round(n*p)),rep(0, n-round(n*p))))
      }
      
      group_tr <- as.data.frame(unlist(group_tr))
      colnames(group_tr) <- 'group_tr'
      group_tr$group <- unique(group)
      
      group_tr_data <- merge(as.data.frame(group), as.data.frame(group_tr), by='group')
      vec <- cbind(group_tr_data, as.data.frame(indiv_tr))
    }
    tr_vec_sampled[[j]] <- vec
  }
  return(tr_vec_sampled)
}


#' @export
make_tr_condition <- function(tr_assignement) {
  N <- nrow(tr_assignement)
  group_tr <- tr_assignement[, 'group_tr']
  indiv_tr <- tr_assignement[, 'indiv_tr']
  return(matrix(as.numeric(c(group_tr>0 & indiv_tr>0,
                             group_tr==0 & indiv_tr>0,
                             group_tr>0 & indiv_tr==0,
                             group_tr==0 & indiv_tr==0)),
                  N, 4, dimnames = list(NULL, c('dir_indpsi', 'dir_indphi', 'indpsi','indphi'))))
  }

#' @export
make_estimator_data <- function(tr_assignement, potential_outcomes) {
  tr_condition <- make_tr_condition(tr_assignement)
  obs_outcome <- rowSums(tr_condition*t(potential_outcomes))
  estimator_data <- cbind(tr_assignement, as.data.frame(obs_outcome)) 
  return(estimator_data)
}

#' @export
estimates_hierarchical <- function(estimator_data) {
  estimator_data <- data.table(estimator_data)
  
  # causal effects
  estimator_data_group_ov <- estimator_data[, list(y_hat_gr_ov=mean(obs_outcome)), by=c('group', 'group_tr')]  
  estimator_data_pop_ov <- estimator_data_group_ov[, list(y_hat_ov=mean(y_hat_gr_ov)), by=c('group_tr')]

  estimator_data_group <- estimator_data[, list(y_hat_gr=mean(obs_outcome)), by=c('group', 'group_tr', 'indiv_tr')]
  estimator_data_pop <- estimator_data_group[, list(y_hat=mean(y_hat_gr)), by=c('group_tr', 'indiv_tr')]

  direct_psi <- estimator_data_pop[group_tr==1 & indiv_tr==1, y_hat]-estimator_data_pop[group_tr==1 & indiv_tr==0, y_hat]
  direct_phi <- estimator_data_pop[group_tr==0 & indiv_tr==1, y_hat]-estimator_data_pop[group_tr==0 & indiv_tr==0, y_hat]
  indirect <- estimator_data_pop[group_tr==1 & indiv_tr==0, y_hat]-estimator_data_pop[group_tr==0 & indiv_tr==0, y_hat]
  total <- estimator_data_pop[group_tr==1 & indiv_tr==1, y_hat]-estimator_data_pop[group_tr==0 & indiv_tr==0, y_hat]
  overall <- estimator_data_pop_ov[group_tr==1, y_hat_ov]-estimator_data_pop_ov[group_tr==0, y_hat_ov]
  
  # variance estimators
  sample_var_group <- estimator_data[, list(var_hat_gr=var(obs_outcome), K_gr=.N), by=c('group', 'group_tr', 'indiv_tr')] # sigma^2_i
  sample_var_pop <- estimator_data_group[, list(var_hat=var(y_hat_gr)), by=c('group_tr', 'indiv_tr')] # sigma^2_g
  sample_var_pop <- merge(sample_var_pop, estimator_data[ , list(C=(.SD[unique(group), .N])), by=group_tr], by='group_tr')
  
  sample_var_pop_ov <- estimator_data_group_ov[, list(var_hat_ov=var(y_hat_gr_ov)), by=group_tr]
  sample_var_pop_ov <- merge(sample_var_pop_ov, estimator_data[ , list(C=(.SD[unique(group), .N])), by=group_tr], by='group_tr')
  
  var_direct_group <- sample_var_group[, list(var_hat_direct_gr=sum(var_hat_gr/K_gr, na.rm=TRUE)), by=c('group', 'group_tr')] # Var(CE_i^D)
  # 2nd component var hat CE
  second_var_direct <- var_direct_group[, list(sum_var_hat_direct=sum(var_hat_direct_gr), C=.N), by=group_tr] # sum Var(CE_i^D)
  
  direct_group <- estimator_data_group[, list(direct_hat_gr=(.SD[indiv_tr==1, y_hat_gr]-.SD[indiv_tr==0, y_hat_gr])), by=c('group', 'group_tr')]
  # 1st component var hat CE
  first_var_direct <- direct_group[, list(sample_var_hat_direct=var(direct_hat_gr)), by=group_tr] # sigma^2_D
  first_second_var_direct <- merge(first_var_direct, second_var_direct, by='group_tr')
  
  first <- first_second_var_direct[, sample_var_hat_direct/C]
  second <- first_second_var_direct[, sum_var_hat_direct/C]
  C <- first_second_var_direct[, C]
  G <- c(length(estimator_data[, unique(group)]), length(estimator_data[, unique(group)]))
  
  var_direct <-  (1- C/G)*first + (1/G)*second
  names(var_direct) <- first_second_var_direct$group_tr
  
  var_indirect <- sample_var_pop[group_tr==1 & indiv_tr==0, var_hat/C] + sample_var_pop[group_tr==0 & indiv_tr==0, var_hat/C] 
  var_total <- sample_var_pop[group_tr==1 & indiv_tr==1, var_hat/C] + sample_var_pop[group_tr==0 & indiv_tr==0, var_hat/C]
  var_overall <- sample_var_pop_ov[group_tr==1, var_hat_ov/C] + sample_var_pop_ov[group_tr==0, var_hat_ov/C] 

  
  return(list(list(direct_psi_hat=direct_psi, direct_phi_hat=direct_phi, indirect_hat=indirect, total_hat=total, overall_hat=overall),
              list(var_direct_psi_hat=unname(var_direct['1']), var_direct_phi_hat=unname(var_direct['0']), var_indirect_hat=var_indirect,
                   var_total_hat=var_total, var_overall_hat=var_overall)))

  }


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
make_tr_condition_miss <- function(tr_assignement) {
  N <- nrow(tr_assignement)
  region_tr <- tr_assignement[, 'region_tr']
  group_tr <- tr_assignement[, 'group_tr']
  indiv_tr <- tr_assignement[, 'indiv_tr']
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
make_estimator_data_miss <- function(tr_assignement, potential_outcomes) {
  tr_condition <- make_tr_condition_miss(tr_assignement)
  obs_outcome <- rowSums(tr_condition*t(potential_outcomes))
  estimator_data <- cbind(tr_assignement, as.data.frame(obs_outcome)) 
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

