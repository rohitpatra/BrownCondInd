R Code for ‘Nonparametric Notion of Residual and Test for Conditional
Independence’
================

In this article, we discuss the R implementation of the test statistic
for the test of conditional independence of \(X, Y\) given \(Z\). In the
present version the codes provided in this page can only be used when
\(X\) and \(Y\) are real valued and \(Z \in \mathbb{R}^d\) for
\(d\geq 1\).

To download the R package use the following in R:

``` r
library(devtools)
devtools::install_github(repo = "rohitpatra/BrownCondInd",ref='main')
library(BrownCondInd)
```

In the following, we generate \(n=100\) i.i.d. copies of \((X,Y,Z)\)
when \(d=2\). Moreover, we assume that
\(Z\sim N(0, \sigma_z^2 \textbf{I}_{5\times 5})\), \(X=W+Z_1+\epsilon\),
and \(X=W+Z_1+\epsilon'\), where \(Z_1\) is the first coordinate of
\(Z\) and \(\epsilon, \epsilon^\prime,\) and \(W\) are three independent
mean zero Gaussian random variables. Moreover, we assume that
\(\epsilon, \epsilon^\prime,\) and \(W\) are independent of \(Z,\) and
\(var(\epsilon)=var(\epsilon^\prime)=\sigma^2_E,\) and
\(var(W)= \sigma^2_W\). Note that this is the simulation scenario used
in Section 4 of Patra, Sen, and Székely (2015). Note that \(X\) and
\(Y\) are conditionally independent of \(Z\) only when \(\sigma_W=0\).
In the following, we fixed \(\sigma_W\) to be \(0.1\).

``` r
n <-100
d <- 2
sigma.Z <- 0.3
sigma.W <- 0.1
sigma.E <- 0.2
Z <- matrix(rnorm(n*d,0,sigma.Z),nr=n, nc=d)
W <- rnorm(n,0,sigma.W)
eps <- rnorm(n,0,sigma.E)
eps.prime <- rnorm(n,0,sigma.E)
X <- W + Z[,1] +eps
Y <- W + Z[,1] +eps.prime
Data <- cbind(X,Y,Z)
colnames(Data) <- c('X', 'Y', paste('Z.', 1:d,sep=''))
head(Data)
```

    ##               X           Y         Z.1         Z.2
    ## [1,]  0.5670427  0.47537531  0.16580002 -0.34107133
    ## [2,] -0.2737525  0.28535997  0.11948500 -0.11797877
    ## [3,] -0.3291812 -0.72813990 -0.47449526 -0.47926210
    ## [4,] -0.2736762 -0.66972302 -0.38984717  0.07950515
    ## [5,] -0.1695244  0.29729282 -0.17216284  0.83001784
    ## [6,]  0.3335767  0.01865898  0.01082668 -0.28290124

In the following we calculate the test statistic
\(\hat{\mathcal{E}}_n\), see (3.7) of Patra, Sen, and Székely (2015).

``` r
test.stat <- npresid.statistics(Data,d)
```

As the limiting behavior of \(\hat{\mathcal{E}}_n\) is unknown, in the
following we approximate the asymptotic distribution through a model
based bootstrap procedure (see Section 3.2.1 of Patra, Sen, and Székely
(2015)) and evaluate the \(p\)-value of the proposed test. In the
following \`\`boot.replic" denotes the number of bootstrap replications.
We recommend using a bootstrap replication of size \(1000.\)

``` r
out <- npresid.boot(Data,d,boot.replic=100)
```


``` 
## [1] "100 bootstrap samples obtained"
## [1] "At bootstrap iteration 50 of 100"
## [1] "At bootstrap iteration 100 of 100"
```

``` r
str(out, max.level = 1)
```

    ## List of 7
    ##  $ statistic            : num 2.38
    ##  $ p.value              : num 0.73
    ##  $ method               : chr "Cond Indep test: p-values by inverting F_hat to get bootstrap samples"
    ##  $ bandwidth.method     : chr "least-squares cross-validation, see \"np\" package"
    ##  $ data.desrip          : chr "dimension of Z is 2, sample size 100, dimension of (X,Y,Z) is 4, bootstrap replication number100"
    ##  $ bootstrap.stat.values: num [1:100] 8.43 3.89 6.38 1.71 4.66 ...
    ##  $ cond.dist.obj        :List of 6

Here \`\`cond.dist.obj’’ is the the list conataining estimators of
\(F_{X|Z}\), \(F_{Y|Z}\), and \(F_Z\) evaluated at the data points
(denoted by F.x\_z, F.y\_z, and F.z\_hat) and the bandwidth used
(denoted by Fbw.x\_z, Fbw.y\_z, and Fbw.z\_z) to evaluate the
conditional distribution functions. Note that we use the functions
available in the "np’’ package ( see Hayfield and Racine (2008)) to
compute the optimal bandwidths as well as the estimates of the
conditional distribution functions.

``` r
str(out$cond.dist.obj, max.level = 1)
```

    ## List of 6
    ##  $ F.x_z  : num [1:100] 0.893 0.111 0.657 0.651 0.578 ...
    ##  $ F.y_z  : num [1:100] 0.88 0.654 0.271 0.153 0.879 ...
    ##  $ Fbw.x_z:List of 64
    ##   ..- attr(*, "class")= chr "condbandwidth"
    ##  $ Fbw.y_z:List of 64
    ##   ..- attr(*, "class")= chr "condbandwidth"
    ##  $ Fbw.z_z:List of 2
    ##  $ F.z_hat: num [1:100, 1:2] 0.667 0.619 0.111 0.162 0.333 ...

The estimated \(p\)-value of the test procedure is given through

``` r
p.value <- out$p.value 
```

The required R file \`\`Npres\_Fucntions.R’’ can be downloaded
[here](http://stat.columbia.edu/~rohit/Code/Npres_Fucntions.R). Note the
functions require the following R-packages: boot, data.table, energy,
np, and stats.

\# References

<div id="refs" class="references">

<div id="ref-np">

Hayfield, Tristen, and Jeffrey S. Racine. 2008. “Nonparametric
Econometrics: The Np Package.” *Journal of Statistical Software* 27 (5).
<http://www.jstatsoft.org/v27/i05/>.

</div>

<div id="ref-npres">

Patra, Rohit K., Bodhisattva Sen, and Gábor Székely. 2015. “A Consistent
Bootstrap Procedure for the Maximum Score Estimator.” *Statist. Probab.
Lett. (To Appear)*. <http://arxiv.org/abs/1409.3886>.

</div>

</div>
