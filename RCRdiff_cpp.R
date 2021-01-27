setwd("/Users/Echo_Xu/Dropbox/Courses/Research/Code/RCRdiff_code/cpp_code/")
library(Rcpp)
library("matrixStats")
library("truncnorm")
library("RCRnorm")
library("statmod")

sourceCpp('Gibbs.cpp')
#get coefficient of linear regression from positive ctrl
fitWithPosCtrl = function(y, x)
{
  mod1 = stats::lm(y ~ x)
  coefs = stats::coef(mod1)
  unname(coefs)
}

#get prior range of uniform distribution
get_range = function(x, mm = 5)
{
  c(mean(x)-mm*stats::sd(x), mean(x)+mm*stats::sd(x))
}

#get residual from positive ctrl fitted with simple linear regression
get_residual = function(log_dat, RNA_conc, coefs)
  #in matrix format
{
  log_dat - sweep(sweep(RNA_conc, 2, coefs[2, ], '*'), 2, coefs[1, ], '+')
}



#'@title An Integrated a fully integrated Bayesian method for differential expression analysis using raw NanoString nCounter data
#'@description NanoString nCounter is a medium-throughput platform that measures gene or microRNA expression levels.
#'Here is a publication that introduces this platform: Malkov (2009) <doi:10.1186/1756-0500-2-80>. Here is the webpage of NanoString
#'nCounter where you can find detailed information about this platform <https://www.nanostring.com/scientific-content/technology-overview/ncounter-technology>.
#'It has great clinical application, such as diagnosis and prognosis of cancer. Based on RCRnorm developed for normalizing NanoString
#'nCounter data (Jia et al., 2019) and adaptive Bayesian LASSO for variable selection (Park and Casella, 2008; Leng et al., 2014), 
#'we propose a fully integrated Bayesian method, called RCRdiff, to detect differentially expressed genes between different
#'groups of tissue samples.
#'@param dat A list containing data for the 4 probe types: positive control, negative control, housekeeping gene and regular gene.
#'The names for the 4 elements in the list should exactly be: pos_dat, neg_dat, hk_dat and reg_dat, respectively. For an example
#'of the input data format, please refer to the FFPE_dat included in the dataset.
#'The data for each probe type should be a dataframe with rows being genes and column being patients. The number of columns (patients)
#'should be the same for data of all four probe types. The rows of positive control data should have the same order as the postive control
#'RNA amount vector supplied to the function.
#'@param pos_conc A vector of log10 RNA amount of the positive controls. The order of these controls should be the same as the rows of positive control
#'data in dat. The defaut is: log10(c(128, 32, 8, 2, 0.5, 0.125)).
#'@param iter Total number of iterations for Monte Carlo simulation. Default is 8000.
#'@param warmup Number of burnin cycles for Monte Carlo simulation. Default is 5000.
#'@param seed Seed for the MCMC sampling for reproducibility. Default is 1.
#'@param random_init Whether to estimate the starting point from data
#'@param mm Number of standard deviations for the prior uniform range.
#'@param m_ab Number of variance for the prior distribution of mu_a and mu_b.
#'@param Z a binary covariate, indicating which group the ith sample is in, 
#'for example,0 for the control group, and 1 for the treatment group.
#'@details Unlike existing methods that often require normalization performed beforehand, 
#'RCRdiff directly handles raw read counts and jointly models the behaviors of different 
#'types of internal controls along with DE and non-DE gene patterns. Doing so would avoid 
#'efficiency loss caused by ignoring estimation uncertainty from the normalization step in 
#'a sequential approach and thus can offer more reliable statistical inference.
#'@return The function returns a list of elements including: a list of MCMC samples. 
#'The number of MCMC samples quals iter-warmup. 
#'@export
#'@examples
#'data(FFPE_dat)
#'result = RCRdiff_BaLASSO(FFPE_dat)
#'@import truncnorm
#'@import statmod
#'@import matrixStats

RCRdiff_BaLASSO = function(dat, pos_conc = log10(c(128, 32, 8, 2, 0.5, 0.125)),
                           iter = 8000, warmup = 6000, random_init = F, 
                           seed = 1, mm = 3, m_ab = 9,
                           Z = c(rep(0,14),rep(1,14)))
{
  ptm <- proc.time()
  
  set.seed(seed*3723)
  
  #log10 transform the original count data
  pos_dat = log10(dat$pos_dat + 1)
  neg_dat = log10(dat$neg_dat + 1)
  hk_dat = log10(dat$hk_dat + 1)
  reg_dat = log10(dat$reg_dat + 1)
  
  #inverse gamma parameter
  u = v = .01
  #LASSO hyperprior parameter
  eta = nu = 0.78
  #number of each class of genes
  n_hk = dim(hk_dat)[1]
  n_reg = dim(reg_dat)[1]
  n_neg = dim(neg_dat)[1]
  n_pos = dim(pos_dat)[1]
  #number of patients or samples
  n_patient = dim(pos_dat)[2]
  n_group1 = length(which(Z==0))
  n_group2 = length(which(Z==1))
  
  #number of MCMC iteration used to calculate posterior after convergence
  iter_keep = iter - warmup
  #calculate the coefficient for each patient; note: with positive controls.
  all_coef = apply(pos_dat, 2, fitWithPosCtrl, pos_conc)
  
  mu_a_itm = mean(all_coef[1,])
  mu_b_itm = mean(all_coef[2,])
  
  #Jacknife to estimate mean and variance of mu_a and mu_b
  mu_a = numeric()
  mu_b = numeric()
  for (i in 1:500)
  {
    mu_a[i] = mean(sample(all_coef[1, ], n_patient - 2))
    mu_b[i] = mean(sample(all_coef[2, ], n_patient - 2))
  }
  
  mu_a_mu = mean(mu_a)
  cat(mu_a_mu, '\n')
  mu_b_mu = mean(mu_b)
  cat(mu_b_mu, '\n')
  sigma2_mu_a = m_ab * stats::var(mu_a)
  cat(sigma2_mu_a, '\n')
  sigma2_mu_b = m_ab * stats::var(mu_b)
  cat(sigma2_mu_b, '\n')
  
  #estimate genes' mean expression level 
  
  hk_RNA = sweep(sweep(hk_dat, 2, all_coef[1,], '-'), 2, all_coef[2,], '/')
  reg_RNA = sweep(sweep(reg_dat, 2, all_coef[1,], '-'), 2, all_coef[2,], '/')
  alpha_reg = sweep(sweep(reg_dat[,which(Z==0)], 2, all_coef[1,which(Z==0)], '-'), 2, all_coef[2,which(Z==0)], '/')
  beta_reg = sweep(sweep(reg_dat[,which(Z==1)], 2, all_coef[1,which(Z==1)], '-'), 2, all_coef[2,which(Z==1)], '/') - rowMeans(alpha_reg)
  
  #estimate genes' mean expression level range
  lambda_hk_range = apply(hk_RNA, 1, get_range, mm = mm)
  alpha_reg_range = apply(alpha_reg, 1, get_range, mm = mm)
  beta_reg_range = apply(beta_reg, 1, get_range, mm = mm)
  
  #estimate patient effect range by two way ANOVA with patient's regular gene expression level.
  gene = factor(rep(1:n_reg, n_patient))
  patient = factor(rep(1:n_patient, each = n_reg))
  
  mod = stats::lm(unlist(reg_RNA) ~ patient + gene, contrasts = list(patient = 'contr.sum', gene = 'contr.sum'))
  
  phi = numeric(n_patient)
  phi[1:(n_patient - 1)] = summary(mod)$coefficients[2:n_patient, 1]
  phi[n_patient] = -sum(phi)
  
  phi_L = phi - mm * summary(mod)$coefficients[2, 2]
  phi_U = phi + mm * summary(mod)$coefficients[2, 2]
  
  
  #initialize all the parameters
  aa = matrix(NA, ncol = n_patient, nrow = iter_keep)
  bb = matrix(NA, ncol = n_patient, nrow = iter_keep)
  cc = numeric(iter_keep)
  phi_return = matrix(NA, ncol = n_patient, nrow = iter_keep) #patient effect
  kappa_hk = matrix(NA, ncol = n_patient * n_hk, nrow = iter_keep)
  kappa_reg = matrix(NA, ncol = n_patient * n_reg, nrow = iter_keep)
  alpha_r = matrix(NA, ncol = n_reg, nrow = iter_keep)
  beta_r = matrix(NA, ncol = n_reg, nrow = iter_keep)
  invTau2_r = matrix(NA, ncol = n_reg, nrow = iter_keep)
  lambda_tau_r = matrix(NA, ncol = n_reg, nrow = iter_keep)
  lambda_hk = matrix(NA, ncol = n_hk, nrow = iter_keep)
  #lambda_reg = matrix(NA, ncol = n_reg, nrow = iter_keep)
  d_neg = matrix(NA, ncol = n_neg, nrow = iter_keep)
  d_pos = matrix(NA, ncol = n_pos, nrow = iter_keep)
  d_hk = matrix(NA, ncol = n_hk, nrow = iter_keep)
  d_reg = matrix(NA, ncol = n_reg, nrow = iter_keep)
  mu_a = numeric(iter_keep)
  mu_b = numeric(iter_keep)
  sigma2e_neg = numeric(iter_keep)
  sigma2e_phr = numeric(iter_keep)
  sigma2a = numeric(iter_keep)
  sigma2b = numeric(iter_keep)
  sigma2kappa_hk = numeric(iter_keep)
  sigma2kappa_reg = numeric(iter_keep)
  sigma2d_neg = numeric(iter_keep)
  sigma2d_phr = numeric(iter_keep)
  
  
  if (random_init == T)
  {
    mu_a_itm = stats::runif(1, 1, 4)
    mu_b_itm = stats::runif(1, 0, 2)
    sigma2a_itm = stats::runif(1, 0, .01)
    sigma2b_itm = stats::runif(1, 0, .01)
    a_itm = stats::rnorm(n_patient, 2.5, .1)
    b_itm = stats::rnorm(n_patient, .9, .1)
    cc_itm = stats::runif(1, -6, -1)
    
    phi_itm = stats::rnorm(n_patient, 0, 2)
    phi_itm[n_patient] = -sum(phi_itm[1:(n_patient - 1)])
    
    kappa_hk_itm = stats::rnorm(n_hk * n_patient, 0, 1)
    sigma2kappa_hk_itm = stats::runif(1, 0, 1)
    kappa_reg_itm = stats::rnorm(n_reg * n_patient, 0, 1)
    sigma2kappa_reg_itm = stats::runif(1, 0, 1)
    
    
    lambda_hk_itm = stats::rnorm(n_hk, 0, 1)
    #lambda_reg_itm = stats::rnorm(n_reg, 0, 1)
    
    alpha_r_itm = stats::rnorm(n_reg, 0, 1)
    beta_r_itm = stats::rnorm(n_reg, 0, 1)
    invTau2_itm = statmod::rinv.gaussian(n_reg, 1, 1)
    lambda_tau_itm = stats::rgamma(n_reg, 1, 1)
    
    d_neg_itm = stats::rnorm(n_neg, 0, .01)
    sigma2d_neg_itm = stats::runif(1, 0, .1)
    d_pos_itm = stats::rnorm(n_pos, 0, .01)
    sigma2d_phr_itm = stats::runif(1, 0, .1)
    d_hk_itm = stats::rnorm(n_hk, 0, .01)
    d_reg_itm = stats::rnorm(n_reg, 0, .01)
    
    
    sigma2e_neg_itm = stats::runif(1, 0, .1)
    sigma2e_phr_itm = stats::runif(1, 0, .1)
  }
  
  
  #get initial values; 
  if (random_init == F)
  {
    sigma2a_itm = stats::var(all_coef[1,])
    sigma2b_itm = stats::var(all_coef[2,])
    a_itm = all_coef[1,]
    b_itm = all_coef[2,]
    cc_itm = mean(unlist(sweep(sweep(neg_dat, 2, all_coef[1,], '-'), 2, all_coef[2,], '/')))
    
    phi_itm = phi
    
    lambda_hk_itm = apply(lambda_hk_range, 2, mean)
    alpha_r_itm = apply(alpha_reg_range, 2, mean)
    beta_r_itm = apply(beta_reg_range, 2, mean)
    
    estimate_kappa = sweep(rbind(hk_RNA, reg_RNA), 2, phi_itm, '-')
    estimate_kappa_var = sweep(sweep(rbind(hk_RNA, reg_RNA), 2, phi_itm, '-'), 1, c(lambda_hk_itm, alpha_r_itm), '-')
    
    kappa_hk_itm = as.vector(unlist(estimate_kappa[1:n_hk,]))
    sigma2kappa_hk_itm = stats::var(as.vector(unlist(estimate_kappa_var[1:n_hk,])))     #0.02497646
    kappa_reg_itm = as.vector(unlist(estimate_kappa[(1+n_hk):(n_hk+n_reg),]))
    sigma2kappa_reg_itm = stats::var(as.vector(unlist(estimate_kappa_var[(1+n_hk):(n_hk+n_reg),])) )   #0.1220459
    
    invTau2_itm <- 1 / (beta_r_itm * beta_r_itm)
    lambda_tau_itm <- sqrt(sigma2kappa_reg_itm) / abs(beta_r_itm)
    
    pos_RNA = matrix(rep(pos_conc, n_patient), ncol = n_patient)
    neg_RNA = matrix(rep(cc_itm, n_neg * n_patient), ncol = n_patient)
    
    d_neg_itm = apply(get_residual(neg_dat, neg_RNA, all_coef), 1, mean)
    sigma2d_neg_itm = stats::var(d_neg_itm)
    d_pos_itm = apply(get_residual(pos_dat, pos_RNA, all_coef), 1, mean)
    sigma2d_phr_itm = stats::var(d_pos_itm)
    d_hk_itm = rep(0, n_hk) 
    d_reg_itm = rep(0, n_reg) 
    
    
    sigma2e_neg_itm = stats::var(unlist(sweep(get_residual(neg_dat, neg_RNA, all_coef), 1, d_neg_itm, '-')))
    sigma2e_phr_itm = stats::var(unlist(sweep(get_residual(pos_dat, pos_RNA, all_coef), 1, d_pos_itm, '-')))
  }
  
  for (i in 1:iter)
  {
    pos_dat=as.matrix(pos_dat)
    neg_dat = as.matrix(neg_dat)
    hk_dat = as.matrix(hk_dat)
    reg_dat =as.matrix(reg_dat)
    pos_conc = as.vector(pos_conc)
    # Call cpp function
    res = RCRdiff_cpp(pos_dat,neg_dat,hk_dat,reg_dat,pos_conc,a_itm,b_itm,d_pos_itm,
                      cc_itm,d_neg_itm,kappa_hk_itm,d_hk_itm,phi_itm,kappa_reg_itm,
                      d_reg_itm,sigma2e_neg_itm,sigma2e_phr_itm,sigma2a_itm,
                      sigma2_mu_a,sigma2_mu_b,mu_a_itm,sigma2b_itm,mu_b_itm,
                      phi_L,phi_U,lambda_hk_itm,sigma2kappa_hk_itm,alpha_r_itm,
                      beta_r_itm,invTau2_itm,lambda_tau_itm,sigma2kappa_reg_itm,
                      lambda_hk_range[1,],lambda_hk_range[2,],alpha_reg_range[1,],
                      alpha_reg_range[2,],eta,nu,sigma2d_neg_itm,sigma2d_phr_itm,
                      mu_a_mu,mu_b_mu,u,v,n_hk,n_reg,n_neg,n_pos,n_patient,Z)
    # Update parameters
    a_itm = res$a
    b_itm = res$b
    cc_itm = res$c
    phi_itm = res$phi
    kappa_hk_itm = res$kappa_hk
    kappa_reg_itm = res$kappa_reg
    lambda_hk_itm = res$lambda_hk
    alpha_r_itm = res$alpha_reg
    beta_r_itm = res$beta_reg
    invTau2_itm = res$invTau2
    lambda_tau_itm = res$lambda_tau
    d_neg_itm = res$d_neg
    d_pos_itm = res$d_pos
    d_hk_itm = res$d_hk
    d_reg_itm = res$d_reg
    mu_a_itm = res$mu_a
    mu_b_itm = res$mu_b
      
    sigma2e_neg_itm = res$sigmas["sigma2e_neg"]
    sigma2e_phr_itm = res$sigmas["sigma2e_phr"]
    sigma2a_itm = res$sigmas["sigma2a"]
    sigma2b_itm = res$sigmas["sigma2b"]
    sigma2kappa_hk_itm = res$sigmas["sigma2kappa_hk"]
    sigma2kappa_reg_itm = res$sigmas["sigma2kappa_reg"]
    sigma2d_neg_itm = res$sigmas["sigma2d_neg"]
    sigma2d_phr_itm = res$sigmas["sigma2d_phr"]
    
    if(i > warmup){
      
      j = i - warmup
      
      aa[j,] = a_itm
      bb[j,] = b_itm
      cc[j] = cc_itm
      phi_return[j,] = phi_itm
      kappa_hk[j,] = kappa_hk_itm
      kappa_reg[j,] = kappa_reg_itm
      lambda_hk[j,] = lambda_hk_itm
      alpha_r[j,] = alpha_r_itm
      beta_r[j,] = beta_r_itm
      invTau2_r[j,] = invTau2_itm
      lambda_tau_r[j,] = lambda_tau_itm
      
      d_neg[j,] = d_neg_itm
      d_pos[j,] = d_pos_itm
      d_hk[j,] = d_hk_itm
      d_reg[j,] = d_reg_itm
      mu_a[j] = mu_a_itm
      mu_b[j] = mu_b_itm
      sigma2e_neg[j] = sigma2e_neg_itm
      sigma2e_phr[j] = sigma2e_phr_itm
      sigma2a[j] = sigma2a_itm
      sigma2b[j] = sigma2b_itm
      sigma2kappa_hk[j] = sigma2kappa_hk_itm
      sigma2kappa_reg[j] = sigma2kappa_reg_itm
      sigma2d_neg[j] = sigma2d_neg_itm
      sigma2d_phr[j] = sigma2d_phr_itm
    }
  }
  
  #get the mcmc samples of all parameters.
  mcmc.samples = list(aa = aa,
                      bb = bb,
                      d_pos = d_pos,
                      d_neg = d_neg,
                      cc = cc,
                      mu_a = mu_a,
                      mu_b = mu_b,
                      phi = phi_return,
                      alpha_r = alpha_r,
                      beta_r = beta_r,
                      invTau2_r = invTau2_r,
                      lambda_tau_r = lambda_tau_r,
                      kappa_reg = kappa_reg,
                      d_hk = d_hk,
                      d_reg = d_reg,
                      sigma2a = sigma2a,
                      sigma2b = sigma2b,
                      sigma2kappa_reg = sigma2kappa_reg,
                      sigma2kappa_hk = sigma2kappa_hk,
                      sigma2e_neg = sigma2e_neg,
                      sigma2e_phr = sigma2e_phr,
                      sigma2d_neg = sigma2d_neg,
                      sigma2d_phr = sigma2d_phr)
  
  
  
  print(proc.time() - ptm)
  
  return(mcmc.samples)
  
}

