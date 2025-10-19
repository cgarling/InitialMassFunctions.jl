module InitialMassFunctions

import Distributions: ContinuousUnivariateDistribution, Pareto, LogNormal, truncated, Truncated, mean, median, var, skewness, kurtosis, pdf, logpdf, cdf, ccdf, minimum, maximum, partype, quantile, cquantile, sampler, rand, Sampleable, Univariate, Continuous, eltype, params
using IrrationalConstants: sqrt2, sqrt2π, sqrthalfπ, logten 
using Random: AbstractRNG
using SpecialFunctions: erf, erfinv
using StaticArrays: SVector, MVector

""" Abstract type for IMFs; a subtype of `Distributions.ContinuousUnivariateDistribution`, as all IMF models can be described as continuous, univariate PDFs. """
abstract type AbstractIMF <: ContinuousUnivariateDistribution end

include("powerlaw.jl")
include("lognormal.jl")
include("precompile.jl")

export PowerLawIMF, Salpeter1955, Kroupa2001, Chabrier2001BPL # power law constructors
export LogNormalIMF, Chabrier2003, Chabrier2003System, Chabrier2001LogNormal   # lognormal constructors
export BrokenPowerLaw, LogNormalBPL, mean, median, var, skewness, kurtosis, pdf, logpdf, cdf, ccdf, partype, minimum, maximum, quantile, quantile!, cquantile
# export normalization, slope, logslope, dndm, dndlogm, pdf, logpdf, cdf, median

end # module
