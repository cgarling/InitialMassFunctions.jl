module InitialMassFunctions

import Distributions: ContinuousUnivariateDistribution, Pareto, LogNormal, truncated, Truncated, mean, median, var, skewness, kurtosis, pdf, logpdf, cdf, ccdf, minimum, maximum, partype, quantile, cquantile, sampler, rand, Sampleable, Univariate, Continuous, eltype, params
import Random: AbstractRNG
import SpecialFunctions: erf,erfinv
# import StaticArrays:MVector,SVector

""" Abstract type for IMFs; a subtype of `Distributions.ContinuousUnivariateDistribution`, as all IMF models can be described as continuous, univariate PDFs. """
abstract type AbstractIMF <: ContinuousUnivariateDistribution end

include("powerlaw.jl")
include("lognormal.jl")

export PowerLaw, Salpeter1955, Kroupa2001, Chabrier2001BPL
export LogNormalIMF, Chabrier2003, Chabrier2001LogNormal # lognormal constructors
export BrokenPowerLaw, LogNormalBPL, mean, median, var, skewness, kurtosis, pdf, logpdf, cdf, ccdf, partype, minimum, maximum, quantile, quantile!, cquantile
# export normalization, slope, logslope, dndm, dndlogm, pdf, logpdf, cdf, median

end # module
