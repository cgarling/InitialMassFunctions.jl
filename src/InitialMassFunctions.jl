module InitialMassFunctions

import Distributions: ContinuousUnivariateDistribution, Pareto, LogNormal, truncated, pdf, logpdf, cdf, ccdf, minimum, maximum, partype, quantile, sampler, rand, Sampleable, Univariate, Continuous
import Random:AbstractRNG

include("common.jl")
include("powerlaw.jl")

export PowerLawIMF,Salpeter1955,Chabrier2001Lognormal,Chabrier2003
export BrokenPowerLaw, pdf, logpdf, cdf, ccdf, partype, eltype, minimum, maximum, quantile
# export normalization, slope, logslope, dndm, dndlogm, pdf, logpdf, cdf, median

end # module
