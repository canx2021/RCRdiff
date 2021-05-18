## Simulate data
simulate_FFPE_follow_real_data_diff = function(mu_a = real_param$mu_a, mu_b = real_param$mu_b, cc = real_param$cc, lambda_hk = real_param$lambda_hk, 
                                               lambda_reg = real_param$lambda_reg,phi = real_param$phi,
                                               sigma2a = real_param$sigma2a, sigma2b = real_param$sigma2b, pos_conc = log10(c(128, 32, 8, 2, 0.5, 0.125)),
                                               n_patient = 28, n_pos = 6, n_neg = 8, n_hk = 7, n_reg = 83, sigma2kappa_hk = real_param$sigma2kappa_hk, 
                                               sigma2kappa_reg = real_param$sigma2kappa_reg, sd_extra = 0, prop_diff=0.1, Z = c(rep(0,14),rep(1,14)),
                                               sigma2d_neg = real_param$sigma2d_neg, sigma2d_phr = real_param$sigma2d_phr, sigma2e_neg = real_param$sigma2e_neg, 
                                               sigma2e_phr = real_param$sigma2e_phr, mean_diff = 1, seed = 2, dr_perc = 0.05, random_dist = "F")
{
  ## Simulate a and b in RCR model
  aa = rnorm(n_patient, mu_a, sqrt(sigma2a))
  bb = rnorm(n_patient, mu_b, sqrt(sigma2b))
  ## Generate housekeeping gene expression
  kappa_hk_vector = rnorm(rep(1,n_patient*n_hk), rep(lambda_hk,n_patient), sqrt(sigma2kappa_hk))
  n_grp1 = length(which(Z==0))
  n_grp2 = length(which(Z==1))
  n_diff = round(n_reg * prop_diff)
  dr_diff= round(n_diff * dr_perc)
  ur_diff = n_diff- dr_diff
  ## Sample alpha of regular genes from real data
  alpha_reg = sample(lambda_reg,n_reg,replace = TRUE)
  ## Generate beta of regular genes
  org_ind = c(rep(1,ur_diff),rep(-1,dr_diff),rep(0,n_reg-n_diff))
  ind_beta = sample(org_ind)
  labels = ind_beta
  labels[which(labels!=0)]=1
  beta_reg = rep(0,n_reg)
  beta_reg[ind_beta==1]=rnorm(length(which(ind_beta==1)),mean_diff,0.1)
  beta_reg[ind_beta==-1]=rnorm(length(which(ind_beta==-1)),-mean_diff,0.1)
  ## Generate regular gene expression 
  kappa_reg_g1 = rnorm(rep(1,n_reg*n_grp1),rep(alpha_reg,n_grp1),sqrt(sigma2kappa_reg))
  kappa_reg_g2 = rnorm(rep(1,n_reg*n_grp2),rep(alpha_reg+beta_reg,n_grp2),sqrt(sigma2kappa_reg))
  kappa_reg = c(kappa_reg_g1,kappa_reg_g2)
  ## Add sample effect
  kappa_hk = matrix(kappa_hk_vector, nrow = n_hk)
  kappa_reg = matrix(kappa_reg,nrow=n_reg)
  kappa = rbind(kappa_hk, kappa_reg)
  phi_s = sample(phi,n_patient-1,replace = TRUE)
  phi_a = c(phi_s,-sum(phi_s))
  ## Add white noise
  extra_kappa = matrix(rnorm((n_reg + n_hk)*n_patient, 0, 0.4), ncol = n_patient)
  kappa_patient = sweep(kappa, 2, phi_a, '+')+ extra_kappa  
  ##Generate pos and neg control data
  pos_conc_mat = matrix(rep(pos_conc, n_patient), ncol = n_patient)
  neg_conc_mat = matrix(rep(cc, n_neg*n_patient), ncol = n_patient, nrow = n_neg)
  ctrl_mat = rbind(neg_conc_mat, pos_conc_mat)
  ## Combine all parts
  full_conc_mat = rbind(ctrl_mat, kappa_patient)
  ## Add probe-specific effect
  dd = c(rnorm(n_neg, 0, sqrt(sigma2d_neg)), rnorm(n_pos+n_hk+n_reg, 0, sqrt(sigma2d_phr)))
  ## Add all effects to data
  norm_conc_mat = sweep(sweep(full_conc_mat, 2, bb, '*'), 2, aa, '+')
  norm1_conc_mat = sweep(norm_conc_mat, 1, dd, '+')
  neg_error = matrix(rnorm(n_patient*n_neg, 0, sqrt(sigma2e_neg)), nrow = n_neg)
  phr_error = matrix(rnorm(n_patient*(n_pos+n_hk+n_reg), 0, sqrt(sigma2e_phr)), ncol = n_patient)
  error_mat = rbind(neg_error, phr_error)
  dat = data.frame(codeclass = c(rep('Negative', n_neg), rep('Positive', n_pos), rep('Housekeeping', n_hk), 
                                 rep('Endogenous', n_reg)), round(10^(norm1_conc_mat + error_mat)))
  
  return(list(dat = dat, pos_dat = dat[dat$codeclass=='Positive', -1], neg_dat = dat[dat$codeclass=='Negative', -1], 
              hk_dat = dat[dat$codeclass=='Housekeeping', -1], reg_dat = dat[dat$codeclass=='Endogenous', -1],
              true_kappa = kappa, sigma2kappa_reg = sigma2kappa_reg, alpha_r = alpha_reg, beta_r = beta_reg,
              sigma2kappa_hk = sigma2kappa_hk, aa = aa, bb = bb, phi = phi_a, kappa_hk_vector = kappa_hk_vector, 
              kappa_reg_vector = as.vector(kappa_reg), labels = labels,
              dd = dd, mu_a = mu_a, mu_b = mu_b, sigma2a = sigma2a, sigma2b = sigma2b, sigma2d_neg = sigma2d_neg, 
              sigma2d_phr = sigma2d_phr, sigma2e_phr = sigma2e_phr, sigma2e_neg = sigma2e_neg))
}
# estimate parameters using FFPE_data
get_real_parameter<- function(dat=FFPE_dat,pos_conc = log10(c(128, 32, 8, 2, 0.5, 0.125))){
  #log10 transform the original count data
  pos_dat = log10(dat$pos_dat + 1)
  neg_dat = log10(dat$neg_dat + 1)
  hk_dat = log10(dat$hk_dat + 1)
  reg_dat = log10(dat$reg_dat + 1)
  #number of each class of genes.
  n_hk = dim(hk_dat)[1]
  n_reg = dim(reg_dat)[1]
  n_neg = dim(neg_dat)[1]
  n_pos = dim(pos_dat)[1]
  #number of patients or samples
  n_patient = dim(pos_dat)[2]
  all_coef = apply(pos_dat, 2, fitWithPosCtrl, pos_conc)
  
  cc = mean(unlist(sweep(sweep(neg_dat, 2, all_coef[1,], '-'), 2, all_coef[2,], '/')))
  
  hk_RNA = sweep(sweep(hk_dat, 2, all_coef[1,], '-'), 2, all_coef[2,], '/')
  reg_RNA = sweep(sweep(reg_dat, 2, all_coef[1,], '-'), 2, all_coef[2,], '/')
  
  ##estimate genes' mean expression level range
  lambda_hk_range = apply(hk_RNA, 1, get_range, mm = 3)
  #lambda_reg_range = apply(reg_RNA, 1, get_range, mm = mm)
  
  lambda_reg_range = apply(reg_RNA, 1, get_range, mm = 3)
  
  ##estimate genes' mean expression level range
  lambda_hk = apply(lambda_hk_range, 2, mean)
  lambda_reg = apply(lambda_reg_range, 2, mean)
  lambda = c(lambda_hk,lambda_reg)
  
  gene = factor(rep(1:n_reg, n_patient))
  patient = factor(rep(1:n_patient, each = n_reg))
  
  mod = stats::lm(unlist(reg_RNA) ~ patient + gene , contrasts = list(patient = 'contr.sum', gene = 'contr.sum'))
  
  phi = numeric(n_patient)
  phi[1:(n_patient - 1)] = summary(mod)$coefficients[2:n_patient, 1]
  phi[n_patient] = -sum(phi)
  
  estimate_kappa = sweep(rbind(hk_RNA, reg_RNA), 2, phi, '-')
  estimate_kappa_var = sweep(sweep(rbind(hk_RNA, reg_RNA), 2, phi, '-'), 1, c(lambda_hk, lambda_reg), '-')
  
  kappa_hk = as.vector(unlist(estimate_kappa[1:n_hk,]))
  sigma2kappa_hk = stats::var(as.vector(unlist(estimate_kappa_var[1:n_hk,])))     #0.02497646
  kappa_reg = as.vector(unlist(estimate_kappa[(1+n_hk):(n_hk+n_reg),]))
  sigma2kappa_reg = stats::var(as.vector(unlist(estimate_kappa_var[(1+n_hk):(n_hk+n_reg),])) )   #0.1220459
  
  pos_RNA = matrix(rep(pos_conc, n_patient), ncol = n_patient)
  neg_RNA = matrix(rep(cc, n_neg * n_patient), ncol = n_patient)
  
  d_neg = apply(get_residual(neg_dat, neg_RNA, all_coef), 1, mean)
  sigma2d_neg = stats::var(d_neg)
  d_pos = apply(get_residual(pos_dat, pos_RNA, all_coef), 1, mean)
  sigma2d_phr = stats::var(d_pos)
  
  sigma2e_neg = stats::var(unlist(sweep(get_residual(neg_dat, neg_RNA, all_coef), 1, d_neg, '-')))
  sigma2e_phr = stats::var(unlist(sweep(get_residual(pos_dat, pos_RNA, all_coef), 1, d_pos, '-')))
  
  real_param = list(mu_a = mean(all_coef[1,]),
                    mu_b = mean(all_coef[2,]),
                    sigma2a = var(all_coef[1,]),
                    sigma2b = var(all_coef[2,]),
                    cc = cc,
                    lambda_hk = lambda_hk,
                    lambda_reg = lambda_reg,
                    phi = phi,
                    sigma2kappa_hk = sigma2kappa_hk,
                    sigma2kappa_reg = sigma2kappa_reg,
                    sigma2d_neg = sigma2d_neg,
                    sigma2d_phr = sigma2d_phr,
                    sigma2e_neg = sigma2e_neg,
                    sigma2e_phr = sigma2e_phr
  )
  
  return(real_param)
}
real_param <- get_real_parameter()
# Simulate data with DE mean 0.5 and DE proportion 50%
Data_sim <- simulate_FFPE_follow_real_data_diff(prop_diff=0.5,sigma2e_neg =(4* real_param$sigma2e_neg),
                                                  sigma2e_phr = (4*real_param$sigma2e_phr),mean_diff = 0.5,
                                                  dr_perc = (1/3),n_patient = 28, n_reg = 83,
                                                  Z = c(rep(0,14),rep(1,14)))
# Results of RCRdiff
result_BaL <- RCRdiff_BaLASSO(Data_sim,iter = 6000, warmup = 2000,seed=1,Z = c(rep(0,14),rep(1,14))) 





