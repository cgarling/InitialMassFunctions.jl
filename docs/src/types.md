# Defined Types
We provide two new types that are subtypes of [`AbstractIMF`](@ref InitialMassFunctions.AbstractIMF), which itself is a subtype of [`Distributions.ContinuousUnivariateDistribution`](https://juliastats.org/Distributions.jl/latest/univariate/#univariates) and generally follow the API provided by `Distributions.jl`. These are

```@docs
InitialMassFunctions.AbstractIMF
BrokenPowerLaw
LogNormalBPL
```