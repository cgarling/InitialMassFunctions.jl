InitialMassFunctions.jl
================

Stellar initial mass functions describe the distribution of initial masses that stars are born with. This package aims to implement and provide interfaces for working with initial mass functions, including but not limited to evaluating and sampling from published distributions.

Published IMFs we include are
```julia
Salpeter1955(mmin::Real=0.4,mmax::Real=Inf)
Chabrier2001BPL(mmin::Real=0.08,mmax::Real=Inf)
Chabrier2001LogNormal(mmin::Real=0.08,mmax::Real=Inf)
Chabrier2003(mmin::Real=0.08,mmax::Real=Inf)
Kroupa2001(mmin::Real=0.08,mmax::Real=Inf)
```

These all return subtypes of [`Distributions.ContinuousUnivariateDistribution`](https://juliastats.org/Distributions.jl/latest/univariate/#univariates) and have many of the typical `Distributions` methods defined for them. These include
 * pdf, logpdf
 * cdf, ccdf
 * `quantile(d, x::Real), quantile(d,x::AbstractArray), quantile!(y::AbstractArray, d, x::AbstractArray)`
 * rand
 * minimum, maximum, extrema
 * mean, median, var, skewness, kurtosis (some of the higher moments don't work for `mmax=Inf` with instances of `BrokenPowerLaw`)
 
Note that `var`, `skewness`, and `kurtosis` are not currently defined for instances of `LogNormalBPL`, such as those returned by the `Chabrier2003` function. 

Many other functions that work on [`Distributions.ContinuousUnivariateDistribution`](https://juliastats.org/Distributions.jl/latest/univariate/#univariates) will also work transparently on these `AbstractIMF` instances.