module InitialMassFunctions

import Distributions: ContinuousUnivariateDistribution, Pareto, LogNormal, Normal, Gamma, Exponential, Beta, Uniform, truncated, Truncated, mean, median, var, skewness, kurtosis, pdf, logpdf, cdf, ccdf, minimum, maximum, partype, quantile, cquantile, sampler, rand, Sampleable, Univariate, Continuous, eltype, params
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

which is 0 for `x <= minimum(d)` and `mean(d)` for `x >= maximum(d)`. For an IMF, `partial_mean(d, x) / mean(d)` is the fraction of the total mass in stars with masses up to `x`. Closed forms are implemented for the IMFs in this package and for `Normal`, `LogNormal`, `Exponential`, `Gamma`, `Beta`, and `Uniform` from Distributions.jl and their truncations; other distributions are integrated numerically with `QuadGK.quadgk`.
"""
function partial_mean(d::ContinuousUnivariateDistribution, x::Real)
    lo, hi = extrema(d)
    return quadgk(m -> m * pdf(d, m), lo, clamp(x, lo, hi))[1]
end

# ∫_a^b x pdf(d, x) dx for a <= b, in closed form; most are a constant times the probability of [a, b] under a related
# distribution
const ClosedFormPartialMean = Union{Normal, LogNormal, Exponential, Gamma, Beta, Uniform}
# P(a < X <= b), from the complementary cdf in the upper tail to avoid cancellation
Δcdf(d, a, b) = (Fa = cdf(d, a); Fa > 1/2 ? ccdf(d, a) - ccdf(d, b) : cdf(d, b) - Fa)
function first_moment(d::Normal, a, b)
    μ, σ = params(d)
    return μ * Δcdf(d, a, b) - σ^2 * (pdf(d, b) - pdf(d, a))
end
function first_moment(d::LogNormal, a, b)
    μ, σ = params(d)
    return exp(μ + σ^2 / 2) * Δcdf(LogNormal(μ + σ^2, σ), a, b)
end
first_moment(d::Exponential, a, b) = d.θ * Δcdf(Gamma(2, d.θ), a, b)
function first_moment(d::Gamma, a, b)
    k, θ = params(d)
    return k * θ * Δcdf(Gamma(k + 1, θ), a, b)
end
function first_moment(d::Beta, a, b)
    α, β = params(d)
    return α / (α + β) * Δcdf(Beta(α + 1, β), a, b)
end
function first_moment(d::Uniform, a, b)
    lo, hi = params(d)
    a, b = clamp(a, lo, hi), clamp(b, lo, hi)
    return (b - a) * (b + a) / (2 * (hi - lo))
end
# minimum(d) is Float64 for some Float32 distributions
partial_mean(d::ClosedFormPartialMean, x::Real) = (lo = partype(d)(minimum(d)); first_moment(d, lo, max(x, lo)))
# Truncated LogNormal has its own method in lognormal.jl, which this would otherwise take precedence over
function partial_mean(d::Truncated{<:Union{Normal, Exponential, Gamma, Beta, Uniform}}, x::Real)
    lo, hi = partype(d).(extrema(d))
    return first_moment(d.untruncated, lo, clamp(x, lo, hi)) / Δcdf(d.untruncated, lo, hi)
end

include("powerlaw.jl")
include("lognormal.jl")
include("precompile.jl")

export PowerLawIMF, Salpeter1955, Kroupa2001, Chabrier2001BPL # power law constructors
export LogNormalIMF, Chabrier2003, Chabrier2003System, Chabrier2001LogNormal   # lognormal constructors
export BrokenPowerLaw, LogNormalBPL, partial_mean, mean, median, var, skewness, kurtosis, pdf, logpdf, cdf, ccdf, partype, minimum, maximum, quantile, quantile!, cquantile
# export normalization, slope, logslope, dndm, dndlogm, pdf, logpdf, cdf, median

end # module
