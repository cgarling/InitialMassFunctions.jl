InitialMassFunctions.jl
================

Stellar initial mass functions describe the distribution of initial masses that stars are born with. This package aims to implement and provide interfaces for working with initial mass functions, including but not limited to evaluating and sampling from published distributions.

Published IMFs we include are
```julia
Salpeter1955(mmin::Real=0.4,mmax::Real=Inf)
Chabrier2001BPL(mmin::Real=0.08,mmax::Real=Inf)
Chabrier2001LogNormal(mmin::Real=0.08,mmax::Real=Inf)
Kroupa2001(mmin::Real=0.08,mmax::Real=Inf)
```

These all return subtypes of `Distributions.ContinuousUnivariateDistribution` and have many of the typical `Distributions` methods defined for them. 