abstract type AbstractIMF end
"""
    distribution(x::AbstractIMF)
Get the probability distribution underpinning the IMF model. This is essentially dN/dM normalized to integrate to 1."""
distribution(x::AbstractIMF) = x.pdf
