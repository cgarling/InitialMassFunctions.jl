module InitialMassFunctions

import Distributions: ContinuousUnivariateDistribution, Pareto, LogNormal, truncated, mean, pdf, logpdf, cdf, ccdf, minimum, maximum, partype, quantile, sampler, rand, Sampleable, Univariate, Continuous
import Random:AbstractRNG
# import StaticArrays:MVector,SVector

include("common.jl")
include("powerlaw.jl")

export PowerLawIMF,Salpeter1955,Chabrier2001LogNormal,Chabrier2003,Kroupa2001,Chabrier2001BPL
export BrokenPowerLaw, mean, pdf, logpdf, cdf, ccdf, partype, eltype, minimum, maximum, quantile
# export normalization, slope, logslope, dndm, dndlogm, pdf, logpdf, cdf, median

end # module
