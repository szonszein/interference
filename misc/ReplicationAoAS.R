## Code to replicate Table 1 from Aronow & Samii (2017,AoAS) using package 'interference'

rm(list=ls())

library(interference)
library(data.table)


setwd('~/Dropbox (Samii-research)/elements-interference/aoas-replicate/')

## Get data
covs <- read.csv('datasets/covariates/covs_y.csv') # data with outcome variables
adjNames <- dir('datasets/networks') # data with adjacency matrices of all schools in AddHealth 
adjNames <- adjNames[c(-66,-67)] # get rid of problematic schools [too small: 'sch-079' (N=2) & 'sch-080' (N=19)] 
numschools <- length(adjNames) # total of 139 networks

## Set parameters for simulation:
p <- .10 # % assigned to treatment
numperm <- 1500 #number of treatment assignements to compute exposure probabilities
numiter <- 500 #replications of experiment

estimates_school <- vector(mode = "list", length = numschools)

for (j in 1:numschools) {

## Clean adjacency matrix---undirected and remove degree=0 individuals
adjmatS <- as.matrix(read.csv(paste('datasets/networks/',adjNames[j],sep="")))

# switch to an undirected graph, zero out diagonal
adjmatS[,2:ncol(adjmatS)]  <- adjmatS[,2:ncol(adjmatS)] + t(adjmatS[,2:ncol(adjmatS)])
adjmatS[adjmatS == 2]  <- 1
diag(adjmatS[,2:ncol(adjmatS)]) <- 0

degRaw <- colSums(adjmatS[,2:ncol(adjmatS)])
idRaw <- adjmatS[,1]

# match Y values
yRaw <- as.numeric(covs$Y[match(idRaw,covs$AID)]) 
yRaw[is.na(yRaw)] <- 0	#(missings for 0's--weird)

# getting zeroed because degree=0
eligible <- degRaw > 0 & !is.na(idRaw)

# adj matrix, outcomes of eligible and N eligible
adj_matrix <- unname(as.matrix(adjmatS[,2:ncol(adjmatS)][eligible,eligible]))
outcomes <- yRaw[eligible]
N <- nrow(adj_matrix)

## Generate potential outcomes
potential_outcome <- make_dilated_out_from_out(outcomes, hop=1, multipliers=c(2,1.5,1.25))

## Compute exposure probabilities
potential_tr_vector <- make_tr_vec_permutation(N,p,R=numperm,seed=4224)
obs_prob_exposure <- make_exposure_prob(potential_tr_vector, adj_matrix, make_exposure_map_AS, list(hop=1))

## Generate observed treatment vectors for simulation
obs_tr_vector_runs <- make_tr_vec_permutation(N,p,numiter, seed = 8765, allow_repetitions = TRUE) # allow_repetitions=TRUE allows to create more treatment vector permutations than the number of possible treatment permutations


run_parameters <- data.table(expand.grid(tr=1:numiter))

for (i in 1:nrow(run_parameters)) {
  if ((i %% 500) == 0) print(paste(Sys.time(), i, nrow(run_parameters)))
  
  tr <- unlist(run_parameters[i, 'tr'])
  
  obs_tr_vector <- obs_tr_vector_runs[tr, ] # user passes treatment vector
  exposure <- make_exposure_map_AS(adj_matrix, obs_tr_vector, hop=1) # returns matrix with dummies for exposure conditions
  obs_outcome <- rowSums(exposure*t(potential_outcome)) # user passes vector of observed outcomes
  
  all_estimates <- estimates(exposure, obs_outcome, obs_prob_exposure, n_var_permutations = 50, hop=1) # returns all the estimators. note: n_var_permutations defines number of iterations to compute constant effects variance estimator (should be around 1000, but we don't care about it here)
  
  outcome_ht <- as.list(all_estimates$yT_ht)
  run_parameters[i, paste0('yT_ht_', names(outcome_ht)):=outcome_ht]
  
  outcome_h <- as.list(all_estimates$yT_h)
  run_parameters[i, paste0('yT_h_', names(outcome_h)):=outcome_h]
  
  estimate_ht <- as.list(all_estimates$tau_ht)
  run_parameters[i, paste0('tau_ht_', names(estimate_ht)):=estimate_ht]
  
  estimate_h <- as.list(all_estimates$tau_h)
  run_parameters[i, paste0('tau_h_', names(estimate_h)):=estimate_h]
  
  var_yT_ht <- as.list(setNames(c(all_estimates$var_yT_ht),rownames(all_estimates$var_yT_ht)))
  run_parameters[i, paste0('var_yT_ht_', names(var_yT_ht)):=var_yT_ht]
  
  var_yT_h <- as.list(setNames(c(all_estimates$var_yT_h),rownames(all_estimates$var_yT_h)))
  run_parameters[i, paste0('var_yT_h_', names(var_yT_h)):=var_yT_h]
  
  cov_yT_ht <- as.list(setNames(c(all_estimates$cov_yT_ht),rownames(all_estimates$cov_yT_ht)))
  run_parameters[i, paste0('cov_yT_ht_', names(cov_yT_ht)):=cov_yT_ht]
  
  cov_yT_h <- as.list(setNames(c(all_estimates$cov_yT_h),rownames(all_estimates$cov_yT_h)))
  run_parameters[i, paste0('cov_yT_h_', names(cov_yT_h)):=cov_yT_h]
  
  var_estimate_ht <- as.list(all_estimates$var_tau_ht) # variance estimator for 1 experiment 
  run_parameters[i, paste0('var_ht_', names(var_estimate_ht)):=var_estimate_ht]
  
  var_estimate_h <- as.list(all_estimates$var_tau_h) # variance estimator for 1 experiment
  run_parameters[i, paste0('var_h_', names(var_estimate_h)):=var_estimate_h]
  
  var_estimate_ht_const_eff <- as.list(all_estimates$var_tau_ht_const_eff) # variance estimator for 1 experiment
  run_parameters[i, paste0('var_ht_const_eff_', names(var_estimate_ht_const_eff)):=var_estimate_ht_const_eff]
  
  var_estimate_ht_max <- as.list(all_estimates$var_tau_ht_max) # variance estimator for 1 experiment
  run_parameters[i, paste0('var_ht_max_', names(var_estimate_ht_max)):=var_estimate_ht_max]
  
  yT <- as.list(rowSums(potential_outcome))
  run_parameters[i, paste0('yT_', names(yT)):=yT]
  
  estimand_tau <- as.list((1/N)*(unlist(yT)-unlist(yT['no']))[names(unlist(yT))!='no'])
  run_parameters[i, paste0('tau_', names(estimand_tau)):=estimand_tau]
  
  run_parameters[i, N:=ncol(potential_outcome)]
  run_parameters[i, school:=adjNames[j]]
  
  }

estimates_school[[j]] <- run_parameters

}

estimates_all_schools <- rbindlist(estimates_school)

fwrite(estimates_all_schools, file = '/Volumes/GoogleDrive/My Drive/interference/replicationAoAS.csv')



## Compute aggregate estimators for all networks
# This aggregation method removes rows for which the estimator yT for a particular network and a particular treatment permutation is NA. yT=NA happens when there's 0 observations in a particular exposure condition, a particular network and a particular treatment permutation.

estimates_all_schools <- fread('/Volumes/GoogleDrive/My Drive/interference/replicationAoAS.csv')

# HT
tau11_ht <- estimates_all_schools[, sum(yT_ht_dir_ind1[!is.na(yT_ht_dir_ind1)])/sum(N[!is.na(yT_ht_dir_ind1)])-sum(yT_ht_no[!is.na(yT_ht_dir_ind1)])/sum(N[!is.na(yT_ht_dir_ind1)]), by=tr]
tau10_ht <- estimates_all_schools[, sum(yT_ht_isol_dir[!is.na(yT_ht_isol_dir)])/sum(N[!is.na(yT_ht_isol_dir)])-sum(yT_ht_no[!is.na(yT_ht_isol_dir)])/sum(N[!is.na(yT_ht_isol_dir)]), by=tr]
tau01_ht <- estimates_all_schools[, sum(yT_ht_ind1)/sum(N)-sum(yT_ht_no)/sum(N), by=tr]

# H
tau11_h <- estimates_all_schools[, sum(yT_h_dir_ind1[!is.na(yT_ht_dir_ind1)])/sum(N[!is.na(yT_ht_dir_ind1)])-sum(yT_h_no[!is.na(yT_ht_dir_ind1)])/sum(N[!is.na(yT_ht_dir_ind1)]), by=tr]
tau10_h <- estimates_all_schools[, sum(yT_h_isol_dir[!is.na(yT_ht_isol_dir)])/sum(N[!is.na(yT_ht_isol_dir)])-sum(yT_h_no[!is.na(yT_ht_isol_dir)])/sum(N[!is.na(yT_ht_isol_dir)]), by=tr]
tau01_h <- estimates_all_schools[, sum(yT_h_ind1)/sum(N)-sum(yT_h_no)/sum(N), by=tr]

# Alternatively, take weighted average of tau_ht_dir_ind1...
tau11_ht_a <- estimates_all_schools[,weighted.mean(tau_ht_dir_ind1,N,na.rm = T),by=tr]
tau10_ht_a <- estimates_all_schools[,weighted.mean(tau_ht_isol_dir,N,na.rm = T),by=tr] 
tau01_ht_a <- estimates_all_schools[,weighted.mean(tau_ht_ind1,N,na.rm = T),by=tr]

tau11_h_a <- estimates_all_schools[,weighted.mean(tau_h_dir_ind1,N,na.rm = T),by=tr] 
tau10_h_a <- estimates_all_schools[,weighted.mean(tau_h_isol_dir,N,na.rm = T),by=tr] 
tau01_h_a <- estimates_all_schools[,weighted.mean(tau_h_ind1,N,na.rm = T),by=tr] 

# var(yT)HT
varyT11_ht <- estimates_all_schools[,sum(var_yT_ht_dir_ind1[!is.na(yT_ht_dir_ind1)])/(sum(N[!is.na(yT_ht_dir_ind1)]))^2, by=tr]
varyT10_ht <- estimates_all_schools[,sum(var_yT_ht_isol_dir[!is.na(yT_ht_isol_dir)])/(sum(N[!is.na(yT_ht_isol_dir)]))^2, by=tr]
varyT01_ht <- estimates_all_schools[,sum(var_yT_ht_ind1)/(sum(N))^2, by=tr]
varyT00_ht <- estimates_all_schools[,sum(var_yT_ht_no)/(sum(N))^2, by=tr]

# var(yT)H
varyT11_h <- estimates_all_schools[,sum(var_yT_h_dir_ind1[!is.na(yT_ht_dir_ind1)])/(sum(N[!is.na(yT_ht_dir_ind1)]))^2, by=tr]
varyT10_h <- estimates_all_schools[,sum(var_yT_h_isol_dir[!is.na(yT_ht_isol_dir)])/(sum(N[!is.na(yT_ht_isol_dir)]))^2, by=tr]
varyT01_h <- estimates_all_schools[,sum(var_yT_h_ind1)/(sum(N))^2, by=tr]
varyT00_h <- estimates_all_schools[,sum(var_yT_h_no)/(sum(N))^2, by=tr]

# cov(yTk,yTl)HT
covyT11yT00_ht <- estimates_all_schools[,sum(`cov_yT_ht_dir_ind1,no`[!is.na(yT_ht_dir_ind1)])/(sum(N[!is.na(yT_ht_dir_ind1)]))^2, by=tr]
covyT10yT00_ht <- estimates_all_schools[,sum(`cov_yT_ht_isol_dir,no`[!is.na(yT_ht_isol_dir)])/(sum(N[!is.na(yT_ht_isol_dir)]))^2, by=tr]
covyT01yT00_ht <- estimates_all_schools[,sum(`cov_yT_ht_ind1,no`)/(sum(N))^2, by=tr]

# cov(yTk,yTl)H
covyT11yT00_h <- estimates_all_schools[,sum(`cov_yT_h_dir_ind1,no`[!is.na(yT_ht_dir_ind1)])/(sum(N[!is.na(yT_ht_dir_ind1)]))^2, by=tr]
covyT10yT00_h <- estimates_all_schools[,sum(`cov_yT_h_isol_dir,no`[!is.na(yT_ht_isol_dir)])/(sum(N[!is.na(yT_ht_isol_dir)]))^2, by=tr]
covyT01yT00_h <- estimates_all_schools[,sum(`cov_yT_h_ind1,no`)/(sum(N))^2, by=tr]

## Compute estimands
tau11 <- estimates_all_schools[, sum(yT_dir_ind1[!is.na(yT_ht_dir_ind1)])/sum(N[!is.na(yT_ht_dir_ind1)])-sum(yT_no[!is.na(yT_ht_dir_ind1)])/sum(N[!is.na(yT_ht_dir_ind1)]), by=tr]
tau10 <- estimates_all_schools[, sum(yT_isol_dir[!is.na(yT_ht_isol_dir)])/sum(N[!is.na(yT_ht_isol_dir)])-sum(yT_no[!is.na(yT_ht_isol_dir)])/sum(N[!is.na(yT_ht_isol_dir)]), by=tr]
tau01 <- estimates_all_schools[, sum(yT_ind1)/sum(N)-sum(yT_no)/sum(N), by=tr]

#Alternatively, weighted avg tau_dir_ind1...
tau11_a <- estimates_all_schools[,weighted.mean(tau_dir_ind1[!is.na(yT_ht_dir_ind1)],N[!is.na(yT_ht_dir_ind1)],na.rm = T),by=tr] 
tau10_a <- estimates_all_schools[,weighted.mean(tau_isol_dir[!is.na(yT_ht_isol_dir)],N[!is.na(yT_ht_isol_dir)],na.rm = T),by=tr] 
tau01_a <- estimates_all_schools[,weighted.mean(tau_ind1,N,na.rm = T),by=tr] 


### Compute estimators and statistics for Table 1
## Bias
#HT
bias_ht <- rbind(mean(tau11_ht$V1-tau11$V1),
mean(tau10_ht$V1-tau10$V1),
mean(tau01_ht$V1-tau01$V1))

#H
bias_h <- rbind(mean(tau11_h$V1-tau11$V1),
mean(tau10_h$V1-tau10$V1),
mean(tau01_h$V1-tau01$V1))

## SD
#HT
sd_ht <- rbind(sd(tau11_ht$V1),
sd(tau10_ht$V1),
sd(tau01_ht$V1))

#H
sd_h <- rbind(sd(tau11_h$V1),
sd(tau10_h$V1),
sd(tau01_h$V1))

# RMSE
rmse <- function(x,tau) mean((x-tau)^2,na.rm=TRUE)^.5
#HT
rmse_ht <- rbind(rmse(tau11_ht$V1,tau11$V1),
rmse(tau10_ht$V1,tau10$V1),
rmse(tau01_ht$V1,tau01$V1))

#H
rmse_h <- rbind(rmse(tau11_h$V1,tau11$V1),
rmse(tau10_h$V1,tau10$V1),
rmse(tau01_h$V1,tau01$V1))

## Compute variance estimators
#HT
var_tau11_ht <- varyT11_ht$V1 + varyT00_ht$V1 - 2*covyT11yT00_ht$V1
var_tau10_ht <- varyT10_ht$V1 + varyT00_ht$V1 - 2*covyT10yT00_ht$V1
var_tau01_ht <- varyT01_ht$V1 + varyT00_ht$V1 - 2*covyT01yT00_ht$V1

#H
var_tau11_h <- varyT11_h$V1 + varyT00_h$V1 - 2*covyT11yT00_h$V1
var_tau10_h <- varyT10_h$V1 + varyT00_h$V1 - 2*covyT10yT00_h$V1
var_tau01_h <- varyT01_h$V1 + varyT00_h$V1 - 2*covyT01yT00_h$V1

## Mean SE
meanse <- function(x) mean(max(0,x)^.5)
#HT
meanse_ht <- rbind(meanse(var_tau11_ht),
meanse(var_tau10_ht),
meanse(var_tau01_ht))

#H
meanseh <- rbind(meanse(var_tau11_h),
meanse(var_tau10_h),
meanse(var_tau01_h))

#95% coverage
alphanorm <- qnorm(.975) # 95%
coverage <- function(est,var,tau) mean(abs(est-tau) < alphanorm*abs(var)^.5,na.rm=TRUE)

#HT
cov95_ht <- rbind(coverage(tau11_ht$V1,var_tau11_ht,tau11$V1),
coverage(tau10_ht$V1,var_tau10_ht,tau10$V1),
coverage(tau01_ht$V1,var_tau01_ht,tau01$V1))

#H
cov95_h <- rbind(coverage(tau11_h$V1,var_tau11_h,tau11$V1),
coverage(tau10_h$V1,var_tau10_h,tau10$V1),
coverage(tau01_h$V1,var_tau01_h,tau01$V1))

#90% coverage
alphanorm <- qnorm(.95) # 90%
coverage <- function(est,var,tau) mean(abs(est-tau) < alphanorm*abs(var)^.5,na.rm=TRUE)

#HT
cov90_ht <- rbind(coverage(tau11_ht$V1,var_tau11_ht,tau11$V1),
coverage(tau10_ht$V1,var_tau10_ht,tau10$V1),
coverage(tau01_ht$V1,var_tau01_ht,tau01$V1))

#H
cov90_h <- rbind(coverage(tau11_h$V1,var_tau11_h,tau11$V1),
coverage(tau10_h$V1,var_tau10_h,tau10$V1),
coverage(tau01_h$V1,var_tau01_h,tau01$V1))

replicationTable1 <- data.table(data.frame(Estimator=c('HT','HT','HT','H','H','H'), Estimand=(rep(c('tau11','tau10','tau01'),2)),
                                Bias=rbind(bias_ht,bias_h),SD=rbind(sd_ht,sd_h),RMSE=rbind(rmse_ht,rmse_h),
                                MeanSE=rbind(meanse_ht,meanseh), CI95cov=rbind(cov95_ht,cov95_h),
                                CI90cov=rbind(cov90_ht,cov90_h)))

setorderv(replicationTable1, cols=c('Estimator', 'Estimand'), order=c(-1, 1))

write.csv(replicationTable1, file = '/Volumes/GoogleDrive/My Drive/interference/replicationAoAStable1.csv')


### Compute aggregate estimators with a different aggregation method
# Converts yT's=NA into 0

estimates_all_schools <- fread('/Volumes/GoogleDrive/My Drive/interference/replicationAoAS.csv')
estimates_all_schools[is.na(yT_ht_dir_ind1), yT_ht_dir_ind1:=0]
estimates_all_schools[is.na(yT_h_dir_ind1), yT_h_dir_ind1:=0]
estimates_all_schools[is.na(yT_ht_isol_dir), yT_ht_isol_dir:=0]
estimates_all_schools[is.na(yT_h_isol_dir), yT_h_isol_dir:=0]
estimates_all_schools[is.na(var_yT_ht_dir_ind1), var_yT_ht_dir_ind1:=0]
estimates_all_schools[is.na(var_yT_h_dir_ind1), var_yT_h_dir_ind1:=0]
estimates_all_schools[is.na(var_yT_ht_isol_dir), var_yT_ht_isol_dir:=0]
estimates_all_schools[is.na(var_yT_h_isol_dir), var_yT_h_isol_dir:=0]

# HT
tau11_ht <- estimates_all_schools[, sum(yT_ht_dir_ind1)/sum(N)-sum(yT_ht_no)/sum(N), by=tr]
tau10_ht <- estimates_all_schools[, sum(yT_ht_isol_dir)/sum(N)-sum(yT_ht_no)/sum(N), by=tr]
tau01_ht <- estimates_all_schools[, sum(yT_ht_ind1)/sum(N)-sum(yT_ht_no)/sum(N), by=tr]

# H
tau11_h <- estimates_all_schools[, sum(yT_h_dir_ind1)/sum(N[!is.na(yT_ht_dir_ind1)])-sum(yT_h_no[!is.na(yT_ht_dir_ind1)])/sum(N[!is.na(yT_ht_dir_ind1)]), by=tr]
tau10_h <- estimates_all_schools[, sum(yT_h_isol_dir[!is.na(yT_ht_isol_dir)])/sum(N[!is.na(yT_ht_isol_dir)])-sum(yT_h_no[!is.na(yT_ht_isol_dir)])/sum(N[!is.na(yT_ht_isol_dir)]), by=tr]
tau01_h <- estimates_all_schools[, sum(yT_h_ind1)/sum(N)-sum(yT_h_no)/sum(N), by=tr]


# var(yT)HT
varyT11_ht <- estimates_all_schools[,sum(var_yT_ht_dir_ind1)/(sum(N))^2, by=tr]
varyT10_ht <- estimates_all_schools[,sum(var_yT_ht_isol_dir)/(sum(N))^2, by=tr]
varyT01_ht <- estimates_all_schools[,sum(var_yT_ht_ind1)/(sum(N))^2, by=tr]
varyT00_ht <- estimates_all_schools[,sum(var_yT_ht_no)/(sum(N))^2, by=tr]

# var(yT)H
varyT11_h <- estimates_all_schools[,sum(var_yT_h_dir_ind1)/(sum(N))^2, by=tr]
varyT10_h <- estimates_all_schools[,sum(var_yT_h_isol_dir)/(sum(N))^2, by=tr]
varyT01_h <- estimates_all_schools[,sum(var_yT_h_ind1)/(sum(N))^2, by=tr]
varyT00_h <- estimates_all_schools[,sum(var_yT_h_no)/(sum(N))^2, by=tr]

# cov(yTk,yTl)HT
covyT11yT00_ht <- estimates_all_schools[,sum(`cov_yT_ht_dir_ind1,no`)/(sum(N))^2, by=tr]
covyT10yT00_ht <- estimates_all_schools[,sum(`cov_yT_ht_isol_dir,no`)/(sum(N))^2, by=tr]
covyT01yT00_ht <- estimates_all_schools[,sum(`cov_yT_ht_ind1,no`)/(sum(N))^2, by=tr]

# cov(yTk,yTl)H
covyT11yT00_h <- estimates_all_schools[,sum(`cov_yT_h_dir_ind1,no`)/(sum(N))^2, by=tr]
covyT10yT00_h <- estimates_all_schools[,sum(`cov_yT_h_isol_dir,no`)/(sum(N))^2, by=tr]
covyT01yT00_h <- estimates_all_schools[,sum(`cov_yT_h_ind1,no`)/(sum(N))^2, by=tr]

## Compute estimands
tau11 <- estimates_all_schools[, sum(yT_dir_ind1)/sum(N)-sum(yT_no)/sum(N), by=tr]
tau10 <- estimates_all_schools[, sum(yT_isol_dir)/sum(N)-sum(yT_no)/sum(N), by=tr]
tau01 <- estimates_all_schools[, sum(yT_ind1)/sum(N)-sum(yT_no)/sum(N), by=tr]

### Compute estimators and statistics for Table 1
## Bias
#HT
bias_ht <- rbind(mean(tau11_ht$V1-tau11$V1),
                 mean(tau10_ht$V1-tau10$V1),
                 mean(tau01_ht$V1-tau01$V1))

#H
bias_h <- rbind(mean(tau11_h$V1-tau11$V1),
                mean(tau10_h$V1-tau10$V1),
                mean(tau01_h$V1-tau01$V1))

## SD
#HT
sd_ht <- rbind(sd(tau11_ht$V1),
               sd(tau10_ht$V1),
               sd(tau01_ht$V1))

#H
sd_h <- rbind(sd(tau11_h$V1),
              sd(tau10_h$V1),
              sd(tau01_h$V1))

# RMSE
rmse <- function(x,tau) mean((x-tau)^2,na.rm=TRUE)^.5
#HT
rmse_ht <- rbind(rmse(tau11_ht$V1,tau11$V1),
                 rmse(tau10_ht$V1,tau10$V1),
                 rmse(tau01_ht$V1,tau01$V1))

#H
rmse_h <- rbind(rmse(tau11_h$V1,tau11$V1),
                rmse(tau10_h$V1,tau10$V1),
                rmse(tau01_h$V1,tau01$V1))

## Compute variance estimators
#HT
var_tau11_ht <- varyT11_ht$V1 + varyT00_ht$V1 - 2*covyT11yT00_ht$V1
var_tau10_ht <- varyT10_ht$V1 + varyT00_ht$V1 - 2*covyT10yT00_ht$V1
var_tau01_ht <- varyT01_ht$V1 + varyT00_ht$V1 - 2*covyT01yT00_ht$V1

#H
var_tau11_h <- varyT11_h$V1 + varyT00_h$V1 - 2*covyT11yT00_h$V1
var_tau10_h <- varyT10_h$V1 + varyT00_h$V1 - 2*covyT10yT00_h$V1
var_tau01_h <- varyT01_h$V1 + varyT00_h$V1 - 2*covyT01yT00_h$V1

## Mean SE
meanse <- function(x) mean(max(0,x)^.5)
#HT
meanse_ht <- rbind(meanse(var_tau11_ht),
                   meanse(var_tau10_ht),
                   meanse(var_tau01_ht))

#H
meanseh <- rbind(meanse(var_tau11_h),
                 meanse(var_tau10_h),
                 meanse(var_tau01_h))

#95% coverage
alphanorm <- qnorm(.975) # 95%
coverage <- function(est,var,tau) mean(abs(est-tau) < alphanorm*abs(var)^.5,na.rm=TRUE)

#HT
cov95_ht <- rbind(coverage(tau11_ht$V1,var_tau11_ht,tau11$V1),
                  coverage(tau10_ht$V1,var_tau10_ht,tau10$V1),
                  coverage(tau01_ht$V1,var_tau01_ht,tau01$V1))

#H
cov95_h <- rbind(coverage(tau11_h$V1,var_tau11_h,tau11$V1),
                 coverage(tau10_h$V1,var_tau10_h,tau10$V1),
                 coverage(tau01_h$V1,var_tau01_h,tau01$V1))

#90% coverage
alphanorm <- qnorm(.95) # 90%
coverage <- function(est,var,tau) mean(abs(est-tau) < alphanorm*abs(var)^.5,na.rm=TRUE)

#HT
cov90_ht <- rbind(coverage(tau11_ht$V1,var_tau11_ht,tau11$V1),
                  coverage(tau10_ht$V1,var_tau10_ht,tau10$V1),
                  coverage(tau01_ht$V1,var_tau01_ht,tau01$V1))

#H
cov90_h <- rbind(coverage(tau11_h$V1,var_tau11_h,tau11$V1),
                 coverage(tau10_h$V1,var_tau10_h,tau10$V1),
                 coverage(tau01_h$V1,var_tau01_h,tau01$V1))

replicationTable1AltAgg<- data.table(data.frame(Estimator=c('HT','HT','HT','H','H','H'), Estimand=(rep(c('tau11','tau10','tau01'),2)),
                                           Bias=rbind(bias_ht,bias_h),SD=rbind(sd_ht,sd_h),RMSE=rbind(rmse_ht,rmse_h),
                                           MeanSE=rbind(meanse_ht,meanseh), CI95cov=rbind(cov95_ht,cov95_h),
                                           CI90cov=rbind(cov90_ht,cov90_h)))

setorderv(replicationTable1AltAgg, cols=c('Estimator', 'Estimand'), order=c(-1, 1))

write.csv(replicationTable1AltAgg, file = '/Volumes/GoogleDrive/My Drive/interference/replicationAoAStable1AltAgg.csv')

