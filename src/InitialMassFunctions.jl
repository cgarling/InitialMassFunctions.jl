module InitialMassFunctions

import Distributions: ContinuousUnivariateDistribution, Pareto, LogNormal, truncated, Truncated, mean, median, var, skewness, kurtosis, pdf, logpdf, cdf, ccdf, minimum, maximum, partype, quantile, cquantile, sampler, rand, Sampleable, Univariate, Continuous, eltype, params
using IrrationalConstants: sqrt2, sqrt2π, sqrthalfπ, logten 
using QuadGK: quadgk
using Random: AbstractRNG
using SpecialFunctions: erf, erfinv
using StaticArrays: SVector, MVector, SUnitRange

""" Abstract type for IMFs; a subtype of `Distributions.ContinuousUnivariateDistribution`, as all IMF models can be described as continuous, univariate PDFs. """
abstract type AbstractIMF <: ContinuousUnivariateDistribution end

"""
    partial_mean(d::Distributions.ContinuousUnivariateDistribution, x::Real)
Returns the first incomplete moment (or partial expectation) of `d` up to `x`,

```math
\\int_{\\text{min}(d)}^{x} \\, m \\, \\frac{dn(m)}{dm} \\, dm
```

which is 0 for `x <= minimum(d)` and `mean(d)` for `x >= maximum(d)`. For an IMF, `partial_mean(d, x) / mean(d)` is the fraction of the total mass in stars with masses up to `x`. Closed forms are implemented for the IMFs in this package; other distributions are integrated numerically with `QuadGK.quadgk`.
"""
function partial_mean(d::ContinuousUnivariateDistribution, x::Real)
    lo, hi = extrema(d)
    return quadgk(m -> m * pdf(d, m), lo, clamp(x, lo, hi))[1]
end

include("powerlaw.jl")
include("lognormal.jl")
include("precompile.jl")

export PowerLawIMF, Salpeter1955, Kroupa2001, Chabrier2001BPL # power law constructors
export LogNormalIMF, Chabrier2003, Chabrier2003System, Chabrier2001LogNormal   # lognormal constructors
export BrokenPowerLaw, LogNormalBPL, partial_mean, mean, median, var, skewness, kurtosis, pdf, logpdf, cdf, ccdf, partype, minimum, maximum, quantile, quantile!, cquantile
# export normalization, slope, logslope, dndm, dndlogm, pdf, logpdf, cdf, median

end # module
