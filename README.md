# RCRdiff
A fully integrated Bayesian method for differential expression analysis using raw NanoString nCounter data
## Description
The medium-throughput mRNA abundance platform NanoString nCounter has gained great popularity in 
the past decade, due to its high sensitivity and technical reproducibility as well as remarkable applicability 
to formalin fixed paraffin embedded (FFPE) tissue samples. 
Based on RCRnorm developed for normalizing NanoString nCounter data (Jia et al., 2019) 
and adaptive Bayesian LASSO for variable selection (Park and Casella, 2008; Leng et al., 2014), 
we propose a fully integrated Bayesian method, called RCRdiff, to detect differentially 
expressed genes between different groups of tissue samples. Unlike existing methods that often 
require normalization performed beforehand, RCRdiff directly handles raw read counts and jointly 
models the behaviors of different types of internal controls along with DE and non-DE gene patterns. 

### Depends
 R (>= 2.15.0),  truncnorm, statmod, matrixStats, Rcpp

## Implementation
* Run RCRdiff_cpp.R
* An example based on simulated data is given in Example.R

## References
<a id="1">[1]</a> 
C. Xu, X. Wang, J. Lim, G. Xiao, Y. Xie, 
RCRdiff: a fully integrated Bayesian method for differential
expression analysis using raw NanoString nCounter data, 
Under review, 2020.

<a id="2">[2]</a> 
G. Jia, X. Wang, Q. Li, W. Lu, X. Tang, I. Wistuba, Y. Xie, et al., Rcrnorm: An
integrated system of random-coefficient hierarchical regression models for normalizing
nanostring ncounter data, The Annals of Applied Statistics 13 (2019) 1617–1647.

<a id="3">[3]</a> 
T. Park, G. Casella, The bayesian lasso, Journal of the American Statistical Association
103 (2008) 681–686.

<a id="4">[4]</a> 
C. Leng, M.-N. Tran, D. Nott, Bayesian adaptive lasso, Annals of the Institute of
Statistical Mathematics 66 (2014) 221–244.

