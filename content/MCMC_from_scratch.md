---
title: "MCMC: An Implementation from Scratch"
Date: 2023-12-28
tags: 
- Bayesian 
---


<img src="pictures/mcmc_chains.png" width="500" height="350">

### Lead up to writing the Metropolis Hastings MCMC Algorithem
1. Read in Data
2. Specify Likilihood Function 
3. Use Optimation to Get Starting Values
4. Specify Prior Function 
5. Specify Posterior Function

After which we can Begin to write the MCMC algorithem. 

### Read in Data
```r
dat = read.csv("/data/bodyfat.csv");n = nrow(dat)
y = dat[,1]; X = as.matrix(cbind(Int = rep(1,n),dat[,-1]))
p = ncol(X); ee = 1e-16 
```

```r
samp.o <- function(t.s) {
  round(
    c(mean = mean(t.s),
      sd = sd(t.s),
      lower = quantile(t.s, 0.025, names = F),
      upper = quantile(t.s, 0.975, names = F)), digits = 6)
}
```

### Specify Likilihood Function 
```r
llike <- function(pars) {
  beta <- pars[1:(length(pars)-1)]
  sigma_sq <- max(pars[(length(pars))],ee)
  sum(log(dnorm(y, mean = X %*% beta, sd = sqrt(sigma_sq))))}
```

### Use Optimation to Get Starting Values
```r
par0 = c(mean(y),rep(0,p-1),3)
opt <- optim(par = par0,fn = llike,method = "BFGS",
             hessian = TRUE,control = list(fnscale=-1))
np <- length(opt$par)
```

### Specify Prior Function 
```r
lprior <- function(pars) {sigma_sq = max(pars[np],ee); log(1 /(sigma_sq))}
```

### Specify Prior Function 

```r
lpost <- function(pars) llike(pars) + lprior(pars)
```

```r
library(mvtnorm)
mh.mcmc <- function(n.s, start.p, start.hessian, burnin, seed = 23, initial_scale_par = 1, n.chain = 3, thinning = 1, par_names = NULL, target_acc_rate = 0.234, learning_rate = 0.05){
  np <- length(start.p)
  draws = matrix(NA, n.s, np)
  scale_par <- initial_scale_par # Initialize scale parameter
  set.seed(seed) 
  
  chains <- list()
  scale_par_vec <- c()
  for (j in 1:n.chain){
    chain_draws = matrix(NA, n.s, np)
    chain_draws[1, ] = start.p
    acc_count = 0 # To count the number of acceptances
    C = -solve(start.hessian) * scale_par # Scale the covariance matrix
    
    for (t in 2:n.s) {
      u = runif(1) 
      prop <- rmvnorm(1, mean = chain_draws[t-1, ], sigma = C)
      if (u < exp(lpost(prop) - lpost(chain_draws[t-1,]))) { # acceptance ratio
        chain_draws[t, ] = prop # accept
        acc_count <- acc_count + 1
      } else { 
        chain_draws[t, ] = chain_draws[t -1, ] } # reject
      
      # Adaptive scaling with learning rate
      if((t > burnin)) {
        scale_par_vec <- c(scale_par_vec, scale_par)
        acc_rate = acc_count 
        scale_adj = learning_rate * (acc_rate - target_acc_rate)
        scale_par <- scale_par * exp(scale_adj)
        C = -solve(start.hessian) * scale_par # Update the covariance matrix
        acc_count = 0 # Reset acceptance count for the next window
      }
    }
    
    # Thinning and removing burn-in samples
    chains[[j]] <- chain_draws[(burnin+1):n.s,][seq(1, n.s - burnin, thinning),]
    
    if(!is.null(par_names)){
      colnames(chains[[j]]) <- par_names
    }
  }
  
  tab <- t(apply(X=do.call(rbind, chains),MARGIN = 2,FUN = samp.o))
  
  return(list(chains = chains, tab = tab, final_acceptance_rate = scale_par, scale_par_vec = scale_par_vec))
}
```

$$
f(x) = \int_{-\infty}^\infty
    f\hat(\xi),e^{2 \pi i \xi x}
    \,d\xi
$$

$$
\begin{table}[!h]
\centering
\begin{tabular}{r|r|r|r}
\hline
mean & sd & lower & upper\\
\hline
-17.984 & 20.872 & -58.703 & 22.852\\
\hline
0.057 & 0.030 & -0.003 & 0.116\\
\hline
-0.086 & 0.058 & -0.201 & 0.027\\
\hline
-0.036 & 0.167 & -0.360 & 0.293\\
\hline
-0.432 & 0.220 & -0.871 & 0.000\\
\hline
-0.017 & 0.098 & -0.209 & 0.175\\
\hline
0.889 & 0.084 & 0.723 & 1.055\\
\hline
-0.195 & 0.137 & -0.468 & 0.072\\
\hline
0.237 & 0.136 & -0.028 & 0.507\\
\hline
-0.024 & 0.233 & -0.475 & 0.427\\
\hline
0.165 & 0.207 & -0.244 & 0.566\\
\hline
0.156 & 0.162 & -0.160 & 0.475\\
\hline
0.433 & 0.188 & 0.061 & 0.804\\
\hline
-1.470 & 0.503 & -2.448 & -0.480\\
\hline
16.120 & 1.489 & 13.479 & 19.287\\
\hline
\end{tabular}
\end{table}
$$
