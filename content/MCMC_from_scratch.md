---
title: "MCMC: An Implementation from Scratch"
date: 2023-12-28
tags: 
- Bayesian 
---


<figure>
  <img src="pictures/mcmc_chains.png" width="600" height="375" alt="Description of Image">
  <figcaption>Generated with DALL-E</figcaption>
</figure>
<!--
> [!warning] Advisory to Reader
> This Post Assumes a basic knowledge of multiple regression, a limited exposure to probability and statistical theory, as well as coding in R
-->


## But first, What is Bayesian Statistics?

If not already clear, Markov Chain Monte Carlo (MCMC) is a tool for implementing Bayesian analyses. Pardon me? If that doesn't make sense, please stay with me; otherwise, please skip ahead. 

> Bayesian statistics is just an alternative paradigm to how students are usually first taught statistics, called *Frequentist Statistics*. 

While some might debate you on this, a common perspective is that neither Bayesian nor Frequentist statistics is necessarily better; what's most important is one recognizes that they are different in implementation but, more importantly, different in their philosophy, the nuances of which will not be covered here. While Frequentist vs Bayesian Statistics is an important discussion, it is also not something you can knock out in 5 minutes. For an extremely brief comparison, see the table below on the difference in Frequentist vs. Bayesian statistics.


<!--
<figure>
  <img src="pictures/bays_vs_freq2.png" width="600" height="500" alt="Description of Image">
  <figcaption>Introna, Michele, et al. "Bayesian statistics in anesthesia practice: a tutorial for anesthesiologists." Journal of anesthesia 36.2 (2022): 294-302.</figcaption>
</figure>
-->

&nbsp;

<!DOCTYPE html>
<html>
<head>
    <title>Comparison Table</title>
    <style>
        table {
            margin-left: auto;
            margin-right: auto;
            width: 80%;
            border-collapse: collapse;
            font-size: 12.5px;
        }
        th, td {
            text-align: left;
            padding: 5px;
        }
        caption {
            caption-side: top; /* Caption position: top or bottom */
            font-size: 18px; /* Adjust caption font size as needed */
            font-weight: bold; /* Optional: makes the caption text bold */
            padding: 5px; /* Optional: adds padding around the caption */
        }
    </style>
</head>
<body>

<table border="1">
  <caption>Bayesian vs. Frequentist Statistics</caption>
  <colgroup>
    <col style="width: 20%;">
    <col style="width: 40%;">
    <col style="width: 40%;">
  </colgroup>
  <tr>
    <th></th>
    <th>Frequentist</th>
    <th>Bayesian</th>
  </tr>
  <tr>
    <td>Answer given</td>
    <td>Probability of the observed data given an underlying truth</td>
    <td>Probability of the truth given the observed data</td>
  </tr>
  <tr>
    <td>Population parameter</td>
    <td>Fixed, but unknown</td>
    <td>Prob. distribution of values</td>
  </tr>
  <tr>
    <td>Outcome measure</td>
    <td>Probability of extreme results, assuming null hypothesis (P value)</td>
    <td>Posterior probability of the hypothesis</td>
  </tr>
  <tr>
    <td>Weaknesses</td>
    <td>Logical inconsistency with clinical decision-making</td>
    <td>Subjectivity in priors' choice; complexity in PKPD modeling</td>
  </tr>
  <tr>
    <td>Strengths</td>
    <td>No need for priors; well-known methods</td>
    <td>Consistency with clinical decision-making</td>
  </tr>
  <tr>
    <td>PKPD application</td>
    <td>Good estimates with large data</td>
    <td>Adaptation of data to individuals</td>
  </tr>
</table>

</body>
</html>

(Introna, Michele, et al, 2022)[^2] 

&nbsp;

## What is MCMC?

MCMC -  for the sake of keeping things simple - is a means of sampling a target distribution to arrive at empirical posterior marginal distributions for our parameters of interest, where the target distribution is the distribution that we assume the population response to be following. What's Notable about MCMC is that it does not require an envelope to sample the target distribution due to the Markov property. In addition, with MCMC, we can arrive at empirical posterior marginal distributions for the parameters of interest without needing conjugate priors. 

Practically, this means there is ample software people can use to do Bayesian analysis using MCMC without needing to know the statistical theory that allows for it to work in the first place. With regard to reducing barriers to entry, this is a big deal! However, understanding how algorithms work is crucial! It helps data scientists make informed decisions about the choice of models and algorithms, ensuring the robustness and accuracy of their analyses. Moreover, Understanding the underlying mechanics of MCMC also aids in diagnosing convergence issues and interpreting the results correctly. 

**Now, on with the post!**

&nbsp;

> [!abstract] Agenda
>
> The purpose of this post is to demonstrate how to implement the Metropolis-Hastings MCMC algorithm from scratch for a multiple linear regression. 

The First implementation will be designated implementation **(a)**; after that, I will add more complexity but also realism in implementations **(b)** through to **(c)**,  mirroring something closer to canned software such as Jags and BUGS. Let's break down some of the first steps.

&nbsp;
## MCMC Outline 

1. Sample a candidate value $\mathbf{X}^*$ from a proposal distribution $g\left(\cdot \mid \mathbf{x}^{(t)}\right)$.
2. Compute the Metropolis-Hastings ratio $R\left(\mathrm{x}^{(t)}, \mathbf{X}^*\right)$, where
$$
R(\mathbf{u}, \mathbf{v})=\frac{f(\mathbf{v}) g(\mathbf{u} \mid \mathbf{v})}{f(\mathbf{u}) g(\mathbf{v} \mid \mathbf{u})} .
$$

Note that $R\left(\mathrm{x}^{(t)}, \mathbf{X}^*\right)$ is always defined, because the proposal $\mathbf{X}^*=\mathrm{x}^*$ can only occur if $f\left(\mathrm{x}^{(t)}\right)>0$ and $g\left(\mathrm{x}^* \mid \mathrm{x}^{(t)}\right)>0$.
3. Sample a value for $\mathbf{X}^{(t+1)}$ according to the following:
$$
\mathbf{X}^{(t+1)}= \begin{cases}\mathbf{X}^* & \text { with probability } \min \left\{R\left(\mathbf{x}^{(t)}, \mathbf{X}^*\right), 1\right\}, \\ \mathbf{x}^{(t)} & \text { otherwise. }\end{cases}
$$
4. Increment $t$ and return to step 1.

**(Givens and Hoeting, 2006)**[^1]

&nbsp;
## Writing our first MCM for implementation (a)
1. Read in Data
2. Write our own Sample Function 
3. Specify Likelihood Function 
4. Use Optimisation to Get Starting Values
5. Specify Prior Function 
6. Specify Posterior Function

Once these five steps are completed, we can begin to write the MCMC algorithm. 

### 1. Reading in the Data
```r
dat = read.csv("data/bodyfat.csv");n = nrow(dat)
y = dat[,1]; X = as.matrix(cbind(Int = rep(1,n),dat[,-1]))
p = ncol(X); ee = 1e-16 
```

&nbsp;

The data is on the percent body fat for 252 adult males, where the objective is to describe 13 simple body measurements in a multiple regression model; the data was collected from **Lohman, T, 1992**[^3] .

<!--
<!DOCTYPE html>
<html>
<head>
<style>
.body {
    width: 100%;
    border-collapse: collapse;
    margin: 20px 0;
    font-family: Arial, sans-serif;
  }
  .body td {
    padding: 8px;
    text-align: left;
    border-bottom: 1px solid #ddd;
  }
  .body tr:nth-child(even) td {
    background-color: #f9f9f9;
  }
  .body tr:hover td {
    background-color: #f1f1f1;
  }
  .body caption {
    padding: 8px;
    font-size: larger;
    caption-side: top;
  }
</style>
</head>
<body>

<table class="body">
<caption>Potential Predictors of Body Fat - Data from Lohman, T, 1992</caption>
<tbody>
  <tr>
    <td>1. Age (years)</td>
    <td>8. Thigh (cm)</td>
  </tr>
  <tr>
    <td>2. Weight (pounds)</td>
    <td>9. Knee (cm)</td>
  </tr>
  <tr>
    <td>3. Height (inches)</td>
    <td>10. Ankle (cm)</td>
  </tr>
  <tr>
    <td>4. Neck (cm)</td>
    <td>11. Extended biceps (cm)</td>
  </tr>
  <tr>
    <td>5. Chest (cm)</td>
    <td>12. Forearm (cm)</td>
  </tr>
  <tr>
    <td>6. Abdomen (cm)</td>
    <td>13. Wrist (cm)</td>
  </tr>
  <tr>
    <td>7. Hip (cm)</td>
    <td></td> 
  </tr>
</tbody>
</table>

</body>
</html>
-->



<!DOCTYPE html>
<html>
<head>
<style>
/* Default styles (Light mode) */
.body {
    width: 100%;
    border-collapse: collapse;
    margin: 20px 0;
    font-family: Arial, sans-serif;
    background-color: #ffffff; /* Light background */
    color: #000000; /* Dark text */
}
.body td {
    padding: 8px;
    text-align: left;
    border-bottom: 1px solid #F6F4D6;
}
.body tr:nth-child(even) td {
    background-color: #f9f9f9;
}
.body tr:hover td {
    background-color: #f1f1f1;
}
.body caption {
    padding: 8px;
    font-size: larger;
    caption-side: top;
}

/* Dark mode specific styles */
@media (prefers-color-scheme: dark) {
  .body, .body caption {
    background-color: #4F483E; /* Dark background */
    color: #F6F4D6; /* Light text */
  }
  .body td {
    border-bottom: 1px solid #464D58; /* Darker border color */
  }
  .body tr:nth-child(even) td {
    background-color: #1a1a1a; /* Darker shade for even rows */
  }
  .body tr:hover td {
    background-color: #4F483E; /* Darker shade on hover */
  }
}
</style>
</head>
<body>

<table class="body">
<caption>Potential Predictors of Body Fat - Data from Lohman, T, 1992</caption>
<tbody>
<tr>
    <td>1. Age (years)</td>
    <td>8. Thigh (cm)</td>
  </tr>
  <tr>
    <td>2. Weight (pounds)</td>
    <td>9. Knee (cm)</td>
  </tr>
  <tr>
    <td>3. Height (inches)</td>
    <td>10. Ankle (cm)</td>
  </tr>
  <tr>
    <td>4. Neck (cm)</td>
    <td>11. Extended biceps (cm)</td>
  </tr>
  <tr>
    <td>5. Chest (cm)</td>
    <td>12. Forearm (cm)</td>
  </tr>
  <tr>
    <td>6. Abdomen (cm)</td>
    <td>13. Wrist (cm)</td>
  </tr>
  <tr>
    <td>7. Hip (cm)</td>
    <td></td> <!-- Empty cell for alignment -->
  </tr>
</tbody>
</table>

</body>
</html>


&nbsp;
### 2. Write our own Sample Function 

Later, we will need a function to summarise the results from the posterior distributions, see `samp.o`.

```r
samp.o <- function(t.s) {
  round(
    c(mean = mean(t.s),
      sd = sd(t.s),
      lower = quantile(t.s, 0.025, names = F),
      upper = quantile(t.s, 0.975, names = F)), digits = 6)
}
```

Now we have our data, we need to specify a Likelihood Function.

&nbsp;
### 3. Specify Log-Likelihood Function 

Because we have a numerical response, it is reasonable to assume that our population resposne is normally distributed:

$$
Y \sim \text{N}(\mu, \sigma^2)
$$

Note, that we are using Log Likelihood for numerical precision. 

```r
llike <- function(pars) {
  beta <- pars[1:(length(pars)-1)]
  sigma_sq <- max(pars[(length(pars))],ee)
  sum(log(dnorm(y, mean = X %*% beta, sd = sqrt(sigma_sq))))}
```

Next, we need to get estimates to use as our starting values as MCMC is sent to where they start chains for each parameter. 

&nbsp;
### 4. Use Optimisation to Get Starting Values
```r
par0 = c(mean(y),rep(0,p-1),3)
opt <- optim(par = par0,fn = llike,method = "BFGS",
             hessian = TRUE,control = list(fnscale=-1))
np <- length(opt$par)
opt$par
```


<!DOCTYPE html>
<html>
<head>
<style>
  .opt td {
    border-top: 2px solid #57B9E9;
    text-align: right;
    font-size: 15px; /* Adjust this value to change the font size of table data */
  }
  .opt th {
    text-align: right; /* Style for headers */
    font-size: 14.5px; /* Adjust this value to change the font size of headers */
  }
</style>
</head>
<body>
<table class="opt">
<thead>
  <tr>
    <th>beta_0</th>
    <th>beta_1</th>
    <th>beta_2</th>
    <th>beta_3</th>
    <th>beta_4</th>
    <th>beta_5</th>
    <th>beta_6</th>
    <th>beta_7</th>
    <th>beta_8</th>
    <th>beta_9</th>
    <th>beta_10</th>
    <th>beta_11</th>
    <th>beta_12</th>
    <th>beta_13</th>
    <th>sigma_sq</th>
  </tr>
</thead>
<tbody>
  <tr>
   <td style="text-align:right;"> -17.8 </td>
   <td style="text-align:right;"> 0.057 </td>
   <td style="text-align:right;"> -0.086 </td>
   <td style="text-align:right;"> -0.037 </td>
   <td style="text-align:right;"> -0.43 </td>
   <td style="text-align:right;"> -0.018 </td>
   <td style="text-align:right;"> 0.89 </td>
   <td style="text-align:right;"> -0.196 </td>
   <td style="text-align:right;"> 0.236 </td>
   <td style="text-align:right;"> -0.021 </td>
   <td style="text-align:right;"> 0.167 </td>
   <td style="text-align:right;"> 0.157 </td>
   <td style="text-align:right;"> 0.429 </td>
   <td style="text-align:right;"> -1.47 </td>
   <td style="text-align:right;"> 15.1 </td>
  </tr>
</tbody>
</table>
</html>

> [!note] Note
>
> These estimates are very similar to what you get with `coef(lm(y ~ ., data = dat))`

&nbsp;

### 5. Specify Prior Function 

Using the Jeffreys prior (an improper prior), we can implement a Random Walk Chain with Metropolis-Hastings Markov Chain Monte Carlo to sample the Normal target distribution where Jeffreys prior is defined as 

$$
f\left(\underline{\beta}, \sigma^2\right) \propto 1 / \sigma^2
$$

which we can write in code simply as 

```r
lprior <- function(pars) {sigma_sq = max(pars[np],ee); log(1 /(sigma_sq))}
```

&nbsp;

### 5. Specify Posterior Function

The Posterior distribution is defined as 
$$
f(\underline{\theta} \mid y)=\frac{L(\underline{\theta} \mid y) f(\underline{\theta})}{\int_0^{\infty} L(\underline{\theta} \mid y) f(\underline{\theta}) d \underline{\theta}} \propto L(\underline{\theta} \mid y) f(\underline{\theta})
$$

where the posterior is simply the prior times the likelihood, while on a log scale, it is simply addition, hence: 

```r
lpost <- function(pars) llike(pars) + lprior(pars)
```

## Implementation (a) 

<!--
	Let's break down its components:

1. **Function Definition**: `mh.mcmc` is the name of the function. It takes several arguments:
    - `n.s`: Number of samples to draw.
    - `start.p`: Starting parameter values.
    - `start.hessian`: The Hessian matrix at the starting point, used to define the proposal distribution's covariance.
    - `burnin`: Number of burn-in iterations to discard.
    - `seed`: Seed for random number generation.
    - `initial_scale_par`: Initial scale parameter for the proposal distribution.
    - `n.chain`: Number of MCMC chains to run.
    - `thinning`: Thinning interval for the MCMC samples.
    - `par_names`: Names of the parameters.
    - `target_acc_rate`: Target acceptance rate for the adaptive algorithm.
    - `learning_rate`: Learning rate for adapting the scale parameter.

2. **Initialization**:
    - The function begins by determining the number of parameters (`np`) and setting up matrices to store the MCMC draws.
    - Initializes the `scale_par` for the proposal distribution.
    - Sets the seed for reproducibility.

3. **MCMC Algorithm**:
    - The function runs `n.chain` separate MCMC chains.
    - For each chain, it initializes the chain's draws and sets up an acceptance count.
    - `C` is the negative inverse of the Hessian matrix times the scale parameter, defining the covariance of the proposal distribution.

4. **Sampling Loop**:
    - For each iteration after the first, it proposes a new sample (`prop`) using a multivariate normal distribution.
    - It calculates an acceptance ratio based on the difference in log-posterior (`lpost`) between the proposed and current sample.
    - The proposed sample is accepted or rejected based on this ratio.

5. **Adaptive Scaling**:
    - After the burn-in period, the algorithm adapts the scale parameter based on the acceptance rate, aiming to achieve the target acceptance rate.
    - It updates the scale parameter and the covariance matrix of the proposal distribution.

6. **Thinning and Burn-in Removal**:
    - After completing all iterations, it thins the samples and removes the burn-in samples.

7. **Post-Processing**:
    - If parameter names are provided, they are assigned to the columns of the chains.
    - `tab` is a table of summary statistics calculated using `samp.o` function (which is not defined in this snippet).
    - The function returns a list containing the MCMC chains, the summary table, the final acceptance rate, and the vector of scale parameters used.

8. **Dependencies**:
    - The function uses the `mvtnorm` library for drawing samples from a multivariate normal distribution.
-->

We will start off with a simpler implementation, which we will denote **(a)**
with manual scaling, Single Chain, Burn-in, No thinning,  and a Jeffrey's Prior.

First, I will present the code and then explain the details see the foldable code chunk **(a)**, below:

> [!info]- (a): manual scaling, Single Chain, Burn-in, No thinning,  and a Jeffrey's Prior
> ```r
> library(mvtnorm)
> np <- length(opt$par)
> lprior <- function(pars) {sigma_sq = max(pars[np],ee); log(1 /(sigma_sq))}
> lpost <- function(pars) llike(pars) + lprior(pars)
> n.s = 100000; draws = matrix(NA, n.s, np)
> draws[1, ] = opt$par; C = -solve(opt$hessian)
> scale_par = .45; burnin = 500; set.seed(23) 
> for (t in 2:n.s) {
>   u = runif(1) 
>   prop <- rmvnorm(1, mean = draws[t-1, ], sigma = scale_par*C)
>   if (u < exp(lpost(prop) - lpost(draws[t-1,]))) { # acceptance ratio
>     draws[t, ] = prop # accept #
>   } else { 
>     draws[t, ] = draws[t -1, ] }} # reject #
> tab.MCMC.a <- t(apply(X=draws[-c(1:burnin),],MARGIN = 2,FUN = samp.o))
> row.names(tab.MCMC.a) <- par_names; tab.MCMC.a
> (acc = (apply(draws[-c(1:burnin),],2,function(x) length(unique(x))/length(x)))[1])
> ```

To start, **(a)** uses Jeffrey's Prior, which is an improper, vague prior; it is improper because its integral does not sum to 1, which is a requirement for valid PDFs.  However, its density function is proportional to the square root of the determinant of the Fisher information matrix:

$$
p(\vec{\theta}) \propto \sqrt{\operatorname{det} \mathcal{I}(\vec{\theta})}
$$

As already stated, Jeffrey's Prior is vague prior, which is also known as a flat, or non-informative prior. Vague priors are guaranteed to play a minimal role in the posterior distribution. 

We are using the multivariate normal distribution with `rmvnorm` to generate proposal values for each of the parameters, and scaling the variance covariance matrix `C` by  `scale_par` to help adjust the acceptance ratio towards the optimal acceptance rate of 23.4 that holds for inhomogeneous target distributions $\pi\left(x^{(d)}\right)=\Pi_{i=1}^d C_i f\left(C_i x_i\right)$ (Roberts and Rosenthal, 2001)[^4].

For implementation **(a)**, we get an acceptance ratio of `0.2267`; this was achieved in part by manually choosing `scale_par` to get an acceptance ratio closer to the target of 0.234. 

Below, we can print our results

```r
tab.MCMC.a
```

<!DOCTYPE html>
<html>
<head>
<style>
  .tabA table {
    margin-left: auto; /* Centers the table horizontally by automatically adjusting the left margin */
    margin-right: auto; /* Centers the table horizontally by automatically adjusting the right margin */
    width: 80%; /* Sets the table width to 80% of the parent element's width */
  }
  .tabA th {
    text-align: left; /* For left-aligned table headers */
    padding: 5px;
    font-size: 17px; /* Set font size for table headers */
  }
  .tabA td {
    text-align: right;
    padding: 5px;
    font-size: 14.5px; /* Set font size for table cells */
  }
</style>
</head>
<body>
<table class="tabA">
<thead>
  <tr>
   <th style="text-align:right;"> Mean </th>
   <th style="text-align:right;"> SD </th>
   <th style="text-align:right;"> Lower </th>
   <th style="text-align:right;"> Upper </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:right;"> -17.890911 </td>
   <td style="text-align:right;"> 20.698141 </td>
   <td style="text-align:right;"> -57.689702 </td>
   <td style="text-align:right;"> 23.027826 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 0.056487 </td>
   <td style="text-align:right;"> 0.030255 </td>
   <td style="text-align:right;"> -0.003148 </td>
   <td style="text-align:right;"> 0.115854 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> -0.086330 </td>
   <td style="text-align:right;"> 0.058003 </td>
   <td style="text-align:right;"> -0.197201 </td>
   <td style="text-align:right;"> 0.028019 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> -0.036785 </td>
   <td style="text-align:right;"> 0.167412 </td>
   <td style="text-align:right;"> -0.362524 </td>
   <td style="text-align:right;"> 0.289261 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> -0.434688 </td>
   <td style="text-align:right;"> 0.222016 </td>
   <td style="text-align:right;"> -0.871323 </td>
   <td style="text-align:right;"> -0.004617 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> -0.016554 </td>
   <td style="text-align:right;"> 0.098487 </td>
   <td style="text-align:right;"> -0.208411 </td>
   <td style="text-align:right;"> 0.176035 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 0.890505 </td>
   <td style="text-align:right;"> 0.084205 </td>
   <td style="text-align:right;"> 0.725644 </td>
   <td style="text-align:right;"> 1.055462 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> -0.194815 </td>
   <td style="text-align:right;"> 0.137902 </td>
   <td style="text-align:right;"> -0.468699 </td>
   <td style="text-align:right;"> 0.074886 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 0.232567 </td>
   <td style="text-align:right;"> 0.138491 </td>
   <td style="text-align:right;"> -0.040740 </td>
   <td style="text-align:right;"> 0.497574 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> -0.020001 </td>
   <td style="text-align:right;"> 0.231950 </td>
   <td style="text-align:right;"> -0.471828 </td>
   <td style="text-align:right;"> 0.433188 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 0.171464 </td>
   <td style="text-align:right;"> 0.207194 </td>
   <td style="text-align:right;"> -0.230577 </td>
   <td style="text-align:right;"> 0.574587 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 0.158884 </td>
   <td style="text-align:right;"> 0.161505 </td>
   <td style="text-align:right;"> -0.155462 </td>
   <td style="text-align:right;"> 0.474175 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 0.432406 </td>
   <td style="text-align:right;"> 0.187672 </td>
   <td style="text-align:right;"> 0.053277 </td>
   <td style="text-align:right;"> 0.797462 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> -1.480241 </td>
   <td style="text-align:right;"> 0.498895 </td>
   <td style="text-align:right;"> -2.459252 </td>
   <td style="text-align:right;"> -0.498583 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 16.226183 </td>
   <td style="text-align:right;"> 1.464961 </td>
   <td style="text-align:right;"> 13.654965 </td>
   <td style="text-align:right;"> 19.337703 </td>
  </tr>
</tbody>
</table>
</html>

> [!tip] Reminder
> 
>  The posterior distribution is what we get out of a Bayesian analysis; we get one distribution for each parameter.

Thus, We can plot the posterior distribution for Weight ($\beta_3$), in which case we can see that the mode of the distribution is similar to the mean point estimate presented in the table above (-0.036). 

<iframe src="https://chart-studio.plotly.com/~dhintz1/4/#/" width="640"
height="480" frameborder="0" allowfullscreen></iframe>

Next in **(b)**, instead of using a Jeffries Prior, we can create `priors_list`, and specify a prior for each parameter. 

## Implementation (b)

> [!info]- (b):  Manual scaling, Single chain, Burn-in, No thinning, Uniform and Normal priors
> ```r
> priors_list <- list(
>   function(theta) dnorm(theta, mean = 0, sd = 1000, log = TRUE),  # beta_0 Prior
>   function(theta) dnorm(theta, mean = 0, sd = 1000, log = TRUE),  # beta_1 Prior 
>   function(theta) dnorm(theta, mean = 0, sd = 1000, log = TRUE),  # beta_2 Prior 
>   function(theta) dnorm(theta, mean = 0, sd = 1000, log = TRUE),  # beta_3 Prior 
>   function(theta) dnorm(theta, mean = 0, sd = 1000, log = TRUE),  # beta_4 Prior
>   function(theta) dnorm(theta, mean = 0, sd = 1000, log = TRUE),  # beta_5 Prior
>   function(theta) dnorm(theta, mean = 0, sd = 1000, log = TRUE),  # beta_6 Prior
>   function(theta) dnorm(theta, mean = 0, sd = 1000, log = TRUE),  # beta_7 Prior
>   function(theta) dnorm(theta, mean = 0, sd = 1000, log = TRUE),  # beta_8 Prior
>   function(theta) dnorm(theta, mean = 0, sd = 1000, log = TRUE),  # beta_9 Prior
>   function(theta) dnorm(theta, mean = 0, sd = 1000, log = TRUE),  # beta_10 Prior
>   function(theta) dnorm(theta, mean = 0, sd = 1000, log = TRUE),  # beta_11 Prior
>   function(theta) dnorm(theta, mean = 0, sd = 1000, log = TRUE),  # beta_12 Prior
>   function(theta) dnorm(theta, mean = 0, sd = 1000, log = TRUE),  # beta_13 Prior
>   function(theta) dunif(theta, min = 1e-6, max = 40,  log = TRUE) # sigma_sq Prior
> )
> # Function to calculate the total log prior
> log_prior_func <- function(theta, priors_list) {
>   log_prior_values <- mapply(function(theta_val, prior_func) prior_func(theta_val), theta, priors_list)
>   sum(log_prior_values)
> }
> # Modify lprior to use the log_prior_func
> lprior <- function(pars) log_prior_func(pars, priors_list)
> # Modify lpost to include the new lprior
> lpost <- function(pars) llike(pars) + lprior(pars)
> for (t in 2:n.s) {
>   u = runif(1) 
>   prop <- rmvnorm(1, mean = draws[t-1, ], sigma = scale_par*C)
>   if (u < exp(lpost(prop) - lpost(draws[t-1,]))) { # acceptance ratio
>     draws[t, ] = prop # accept #
>   } else { 
>     draws[t, ] = draws[t -1, ] }} # reject #
> tab.MCMC.b <- t(apply(X=draws[-c(1:burnin),],MARGIN = 2,FUN = samp.o))
> row.names(tab.MCMC.b) <- par_names; tab.MCMC.b
> (acc = (apply(draws[-c(1:burnin),],2,function(x) length(unique(x))/length(x)))[1])
> ```

Notice we print a summary table of the posterior distributions, the point estimates we get our different from **(a)**.

```r
tab.MCMC.b
```


<!DOCTYPE html>
<html>
<head>
<style>
  .tabB table {
    margin-left: auto; /* Centers the table horizontally by automatically adjusting the left margin */
    margin-right: auto; /* Centers the table horizontally by automatically adjusting the right margin */
    width: 80%; /* Sets the table width to 80% of the parent element's width */
  }
  .tabB th {
    text-align: left; /* For left-aligned table headers */
    padding: 5px;
    font-size: 17px; /* Set font size for table headers */
  }
  .tabB td {
    text-align: right;
    padding: 5px;
    font-size: 14.5px; /* Set font size for table cells */
  }
</style>
</head>
<body>
<table class="tabB">
 <thead>
  <tr>
   <th style="text-align:right;"> Mean </th>
   <th style="text-align:right;"> SD </th>
   <th style="text-align:right;"> Lower </th>
   <th style="text-align:right;"> Upper </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:right;"> -17.890911 </td>
   <td style="text-align:right;"> 20.698141 </td>
   <td style="text-align:right;"> -57.689702 </td>
   <td style="text-align:right;"> 23.027826 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 0.056487 </td>
   <td style="text-align:right;"> 0.030255 </td>
   <td style="text-align:right;"> -0.003148 </td>
   <td style="text-align:right;"> 0.115854 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> -0.086330 </td>
   <td style="text-align:right;"> 0.058003 </td>
   <td style="text-align:right;"> -0.197201 </td>
   <td style="text-align:right;"> 0.028019 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> -0.036785 </td>
   <td style="text-align:right;"> 0.167412 </td>
   <td style="text-align:right;"> -0.362524 </td>
   <td style="text-align:right;"> 0.289261 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> -0.434688 </td>
   <td style="text-align:right;"> 0.222016 </td>
   <td style="text-align:right;"> -0.871323 </td>
   <td style="text-align:right;"> -0.004617 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> -0.016554 </td>
   <td style="text-align:right;"> 0.098487 </td>
   <td style="text-align:right;"> -0.208411 </td>
   <td style="text-align:right;"> 0.176035 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 0.890505 </td>
   <td style="text-align:right;"> 0.084205 </td>
   <td style="text-align:right;"> 0.725644 </td>
   <td style="text-align:right;"> 1.055462 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> -0.194815 </td>
   <td style="text-align:right;"> 0.137902 </td>
   <td style="text-align:right;"> -0.468699 </td>
   <td style="text-align:right;"> 0.074886 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 0.232567 </td>
   <td style="text-align:right;"> 0.138491 </td>
   <td style="text-align:right;"> -0.040740 </td>
   <td style="text-align:right;"> 0.497574 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> -0.020001 </td>
   <td style="text-align:right;"> 0.231950 </td>
   <td style="text-align:right;"> -0.471828 </td>
   <td style="text-align:right;"> 0.433188 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 0.171464 </td>
   <td style="text-align:right;"> 0.207194 </td>
   <td style="text-align:right;"> -0.230577 </td>
   <td style="text-align:right;"> 0.574587 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 0.158884 </td>
   <td style="text-align:right;"> 0.161505 </td>
   <td style="text-align:right;"> -0.155462 </td>
   <td style="text-align:right;"> 0.474175 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 0.432406 </td>
   <td style="text-align:right;"> 0.187672 </td>
   <td style="text-align:right;"> 0.053277 </td>
   <td style="text-align:right;"> 0.797462 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> -1.480241 </td>
   <td style="text-align:right;"> 0.498895 </td>
   <td style="text-align:right;"> -2.459252 </td>
   <td style="text-align:right;"> -0.498583 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 16.226183 </td>
   <td style="text-align:right;"> 1.464961 </td>
   <td style="text-align:right;"> 13.654965 </td>
   <td style="text-align:right;"> 19.337703 </td>
  </tr>
</tbody>
 </table>
</html>

The next layers of complexity we will add for **(c)** are thinning, multiple chains, and Adaptive Scaling. 

## Implementation (c) 

> [!info]- (c): LR Adaptive Scaling, n-chain, Burn-in, Thinning, Jeffrey's Prior
> ```r
> mh.mcmc <- function(n.s, start.p, start.hessian, burnin, seed = 23, initial_scale_par = 1, n.chain = 3, thinning = 1, par_names = NULL, learning_rate = 0.05, target_acc_rate = 0.234){
>   np <- length(start.p)
>   set.seed(seed) 
>   
>   chains <- list()
>   scale_par_history <- vector("list", n.chain) # List to store scale_par history for each chain
>   
>   for (j in 1:n.chain){
>     chain_draws = matrix(NA, n.s, np)
>     chain_draws[1, ] = start.p
>     acc_count = 0 # To count the number of acceptances
>     scale_par <- initial_scale_par # Initialize scale parameter for each chain
>     C = -solve(start.hessian) * scale_par # Scale the covariance matrix
>     scale_par_vec <- numeric(n.s) # Store scale_par values for each iteration
>     scale_par_vec[1] <- scale_par
>     
>     for (t in 2:n.s) {
>       u = runif(1) 
>       prop <- rmvnorm(1, mean = chain_draws[t-1, ], sigma = C)
>       accept_ratio = exp(lpost(prop) - lpost(chain_draws[t-1,]))
>       
>       if (u < accept_ratio) {
>         chain_draws[t, ] = prop # accept
>         acc_count <- acc_count + 1
>       } else { 
>         chain_draws[t, ] = chain_draws[t -1, ] # reject
>       }
>       
>       # Adapt scale_par during burn-in
>       if (t <= burnin) {
>         current_acc_rate = acc_count / t
>         scale_par <- scale_par * exp(learning_rate * (current_acc_rate - target_acc_rate))
>         C = -solve(start.hessian) * scale_par
>       }
>       scale_par_vec[t] <- scale_par
>     }
>     
>     scale_par_history[[j]] <- scale_par_vec # Store scale_par history for the chain
>     
>     # Thinning and removing burn-in samples
>     chains[[j]] <- chain_draws[(burnin+1):n.s,][seq(1, n.s - burnin, thinning),]
>     
>     if(!is.null(par_names)){
>       colnames(chains[[j]]) <- par_names
>     }
>   }
>   
>   tab <- t(apply(X=do.call(rbind, chains),MARGIN = 2,FUN = samp.o))
>   
>   return(list(chains = chains, tab = tab, final_acceptance_rate = acc_count / n.s, scale_par_history = scale_par_history))
> }
> # Example usage
> set.seed(23)
> mcmc.out <- mh.mcmc(n.s = n.s, start.p = opt$par, start.hessian = opt$hessian, burnin = 500, initial_scale_par = 1, n.chain = 3, thinning = 20)
> mcmc.out$tab
> # scale_par_history
> plot(1, type = "n", ylim = c(0,1), xlim = c(1,500), ylab = "scale_par")
> lines(mcmc.out$scale_par_history[[1]][1:500], type = "l", col = "red")
> lines(mcmc.out$scale_par_history[[2]][1:500], type = "l", col = "blue")
> lines(mcmc.out$scale_par_history[[3]][1:500], type = "l", col = "green")
> tail(mcmc.out$scale_par_history[[1]])
> mcmc.out$final_acceptance_rate
> coda::autocorr.plot(mcmc.out$chains[[1]],lag.max=40)
> traceplot_mcmc <- function(mcmc_chains, param_indices = NULL, main_title = "", colors = NULL, par_names = NULL, ask = TRUE, ...){
>   n_chains <- length(mcmc_chains)
>   
>   # If specific parameters are not specified, plot all
>   if (is.null(param_indices)) {
>     param_indices <- 1:ncol(mcmc_chains[[1]])
>   }
>   
>   # Generate colors if not provided
>   if (is.null(colors) || length(colors) < n_chains) {
>     colors <- rainbow(n_chains)
>   }
>   
>   if(ask){
>     par(ask = TRUE)
>   }
>   
>   if(is.null(par_names)){
>     par_names <- paste("Parameter", param_indices)
>   }
>   
>   
>   for (i in param_indices) {
>     plot(NULL, xlim = c(1, nrow(mcmc_chains[[1]])), ylim = range(sapply(mcmc_chains, function(x) x[, i])), 
>          xlab = "Iteration", ylab = par_names[i], main = paste(main_title, par_names[i]), ...)
>     
>     for (j in 1:n_chains) {
>       lines(mcmc_chains[[j]][, i], col = colors[j], ...)
>     }
>   }
>   if(ask){
>     par(ask = FALSE)
>   }
> }
> # Example usage
> par(mfrow = c(2,4))
> traceplot_mcmc(mcmc.out$chains, param_indices = 1:8,par_names = par_names, ask = FALSE)
> traceplot_mcmc(mcmc.out$chains, param_indices = 9:15,par_names = par_names, ask = FALSE)  
> ```

### Thinning 

Thinning is a process that is done to help reduce the autocorrelation present in the chains; if there is severe autocorrelation, then the results are not usable for inference. Autocorrelation in the chains is usually diagnosed with an autocorrelation plot. Looking at the autocorrelation plot, we can see that there is still some concerning autocorrelation (we want almost no spikes for the first half dozen lags); thus, via bumping up the argument `thinning` to be something higher than 20, say 40, we might resolve the issue of autocorrelation. 

```r
par(mfrow = c(2, 2))
for (i in 1:4) {
  coda::autocorr.plot(mcmc.out$chains[[1]][, i],lag.max = 30,auto.layout = FALSE)
  title(main = paste0(par_names[i]),line = 0.8,cex.main = 0.95)
}
mtext(
  "Autocorrelation for first 4 Parameters",
  outer = TRUE,cex = 1,line = -1.5,font = 2
  )
par(mfrow = c(1, 1))
```

<figure>
  <img src="pictures/autocor.png" width="800" height="500" alt="Description of Image">
  <figcaption>Autocorrelation for first 4 Parameters, first Chain</figcaption>
</figure>


### Mixing of Chains 

With adding multiple chains, we have just added an extra loop (in this case, iterating over `j`) to repeat the same sampling procedure for each chain. 

<figure>
  <img src="pictures/mixing.png" width="800" height="500" alt="Description of Image">
  <figcaption>Mixing of 3 chains, first 4 Parameters</figcaption>
</figure>

For reference, the image below shows a good mixing of the first plot, indicating convergence, while plots 2 and 3 show that the MCMC has not converged. 

<figure>
  <img src="pictures/convergence.png" width="800" height="500" alt="Description of Image">
</figure>

(Taboga, Marco, 2021)[^5]

&nbsp;

> [!note] Note
>
> If MCMC did not converge, the results are not usable; the same goess for if there is high autocorrelation shown in the autocorrelation plot.

### Adaptive Scaling 

The complex to explain is Adaptive Scaling. 

<figure>
  <img src="pictures/adapt.png" width="800" height="500" alt="Description of Image">
</figure>

<!--FOOTNOTES-->


[^1]: [Givens and Hoeting, 2006](https://onlinelibrary.wiley.com/doi/book/10.1002/9781118555552)
[^2]: [Introna, Michele, et al, 2022](https://pubmed.ncbi.nlm.nih.gov/35147768/)
[^3]: [Lohman, T, 1992](https://www.sciepub.com/reference/33145)
[^4]: [Roberts and Rosenthal, 2001](https://www.jstor.org/stable/3182776)
[^5]: [Taboga, Marco, 2021](https://statlect.com/fundamentals-of-statistics/Markov-Chain-Monte-Carlo-diagnostics)


