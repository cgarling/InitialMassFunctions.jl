module InitialMassFunctions

import Distributions: ContinuousUnivariateDistribution, Pareto, truncated

include("common.jl")
include("powerlaw.jl")

export PowerLawIMF,Salpeter1955

end # module
