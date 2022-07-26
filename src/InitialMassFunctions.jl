module InitialMassFunctions

import Distributions: ContinuousUnivariateDistribution, Pareto, LogNormal, truncated, mean, median, var, skewness, kurtosis, pdf, logpdf, cdf, ccdf, minimum, maximum, partype, quantile, cquantile, sampler, rand, Sampleable, Univariate, Continuous
import Random:AbstractRNG
import SpecialFunctions:erf,erfinv
# import StaticArrays:MVector,SVector

include("common.jl")
include("powerlaw.jl")
include("lognormal.jl")

export PowerLawIMF, Salpeter1955, Kroupa2001, Chabrier2001BPL
export Chabrier2003, Chabrier2001LogNormal       # lognormal constructors
export BrokenPowerLaw, LogNormalBPL, mean, median, var, skewness, kurtosis, pdf, logpdf, cdf, ccdf, partype, eltype, minimum, maximum, quantile, quantile!
# export normalization, slope, logslope, dndm, dndlogm, pdf, logpdf, cdf, median

end # module
