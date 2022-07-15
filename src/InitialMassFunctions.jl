module InitialMassFunctions

import Distributions: ContinuousUnivariateDistribution, Pareto, truncated, pdf, logpdf

include("common.jl")
include("powerlaw.jl")

export PowerLawIMF,Salpeter1955,Chabrier2003
# export normalization, slope, logslope, dndm, dndlogm, pdf, logpdf, cdf, median

end # module
