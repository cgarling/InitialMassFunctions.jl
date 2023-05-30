InitialMassFunctions.jl
================

[![Build Status](https://github.com/cgarling/InitialMassFunctions.jl/workflows/CI/badge.svg)](https://github.com/cgarling/InitialMassFunctions.jl/actions)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://cgarling.github.io/InitialMassFunctions.jl/stable/)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://cgarling.github.io/InitialMassFunctions.jl/dev/)
[![codecov](https://codecov.io/gh/cgarling/InitialMassFunctions.jl/branch/main/graph/badge.svg?token=XY3O7IQPZH)](https://codecov.io/gh/cgarling/InitialMassFunctions.jl)


Stellar initial mass functions describe the distribution of initial masses that stars are born with. This package aims to implement and provide interfaces for working with initial mass functions, including but not limited to evaluating and sampling from published distributions. See the linked documentation above for more details.

Published IMFs we include are
```julia
Salpeter1955(mmin::Real=0.4, mmax::Real=Inf)
Chabrier2001BPL(mmin::Real=0.08, mmax::Real=Inf)
Chabrier2001LogNormal(mmin::Real=0.08, mmax::Real=Inf)
Chabrier2003(mmin::Real=0.08, mmax::Real=Inf)
Chabrier2003System(mmin::Real=0.08, mmax::Real=Inf)
Kroupa2001(mmin::Real=0.08, mmax::Real=Inf)
```

These all return subtypes of [`Distributions.ContinuousUnivariateDistribution`](https://juliastats.org/Distributions.jl/latest/univariate/#univariates) and have many of the typical methods from [`Distributions.jl`](https://github.com/JuliaStats/Distributions.jl) defined for them. These include
 * pdf, logpdf
 * cdf, ccdf
 * `quantile(d, x::Real), quantile(d,x::AbstractArray), quantile!(y::AbstractArray, d, x::AbstractArray)`
 * rand
 * minimum, maximum, extrema
 * mean, median, var, skewness, kurtosis (some of the higher moments don't work for `mmax=Inf` with instances of `BrokenPowerLaw`)
 
Note that `var`, `skewness`, and `kurtosis` are not currently defined for instances of `LogNormalBPL`, such as those returned by the `Chabrier2003` function.

Many other functions that work on [`Distributions.ContinuousUnivariateDistribution`](https://juliastats.org/Distributions.jl/latest/univariate/#univariates) will also work transparently on these `AbstractIMF` instances.

Continuous broken-power-law distributions (such as those used in [Chabrier 2001](https://ui.adsabs.harvard.edu/abs/2001ApJ...554.1274C/abstract) and [Kroupa 2001](https://ui.adsabs.harvard.edu/abs/2001MNRAS.322..231K/abstract)) are provided through a new `BrokenPowerLaw` type. The lognormal distribution for masses `m<1.0` with a power law extension for higher masses as given in [Chabrier 2003](https://ui.adsabs.harvard.edu/abs/2003PASP..115..763C/abstract) is provided by a new `LogNormalBPL` type. Simpler models (e.g., `Salpeter1955` and `Chabrier2001LogNormal`) are implemented using built-in distributions from [`Distributions.jl`](https://github.com/JuliaStats/Distributions.jl) in concert with their [`truncated`](https://juliastats.org/Distributions.jl/stable/truncate/#Distributions.truncated) function.

Efficient samplers are implemented for the new types `BrokenPowerLaw` and `LogNormalBPL` such that batched calls (e.g., `rand(Chabrier2003(),1000)`) are more efficient than single calls (e.g., `d=Chabrier2003(); [rand(d) for i in 1:1000]`). These samplers can be created explicitly by calling `Distributions.sampler(d)`, with `d` being a `BrokenPowerLaw` or `LogNormalBPL` instance.

## Versioning
`InitialMassFunctions.jl` follows Julia's [recommended versioning strategy](https://pkgdocs.julialang.org/v1/compatibility/#compat-pre-1.0), where breaking changes prior to version 1.0 will result in a bump of the minor version (e.g., 0.1.x |> 0.2.0) whereas feature additions and bug patches will increment the patch version (0.1.0 |> 0.1.1). 