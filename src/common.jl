""" Abstract type for IMFs; a subtype of `Distributions.ContinuousUnivariateDistribution`, as all IMF models can be described as continuous, univariate PDFs. """
abstract type AbstractIMF <: ContinuousUnivariateDistribution end
# """
#     distribution(x::AbstractIMF)
# Get the probability distribution underpinning the IMF model. This is essentially dN/dM normalized to integrate to 1."""
# distribution(x::AbstractIMF) = x.pdf
# Base.Broadcast.broadcastable(x::AbstractIMF) = Ref(x)
# limits(x::AbstractIMF) = x.mmin, x.mmax


# module IMFConstants
# const ln10 = log(10)


# end
