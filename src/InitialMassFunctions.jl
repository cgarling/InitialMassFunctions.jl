module InitialMassFunctions

import Distributions: ContinuousUnivariateDistribution, Pareto, LogNormal, truncated, pdf, logpdf

include("common.jl")
include("powerlaw.jl")

export PowerLawIMF,Salpeter1955,Chabrier2001Lognormal,Chabrier2003
# export normalization, slope, logslope, dndm, dndlogm, pdf, logpdf, cdf, median

end # module
