module InitialMassFunctions

import Distributions: ContinuousUnivariateDistribution, Pareto, truncated

include("common.jl")
include("powerlaw.jl")

export PowerLawIMF,Salpeter1955
export normalization, slope, logslope, dndm, dndlogm, pdf, logpdf, cdf, median

end # module
