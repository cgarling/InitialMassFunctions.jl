###########################################################################################
# Power Law
###########################################################################################
"""
    PowerLawIMF(α::Real, mmin::Real, mmax::Real)
Descibes a single power-law IMF with probability distribution
```math
    \\frac{dn(m)}{dm} = A \\times m^{-\\alpha}
```
truncated such that the probability distribution is 0 below `mmin` and above `mmax`. `A` is a normalization constant such that the distribution integrates to 1 from the minimum valid stellar mass `mmin` to the maximum valid stellar mass `mmax`. This is simply `Distributions.truncated(Distributions.Pareto(α-1,mmin);upper=mmax)`. See the documentation for [`Pareto`](https://juliastats.org/Distributions.jl/latest/univariate/#Distributions.Pareto) and [`truncated`](https://juliastats.org/Distributions.jl/latest/truncate/#Distributions.truncated).
"""
PowerLawIMF(α::Real, mmin::Real, mmax::Real) = truncated(Pareto(α-1, mmin); upper=mmax)
function mean(d::Truncated{Pareto{T}, Continuous, T}) where T
    mmin, mmax = extrema(d)
    α, θ = params( d.untruncated )
    return (α * θ^α / (1-α) / d.ucdf) * (mmax^(1-α) - mmin^(1-α))
end
"""
    Salpeter1955(mmin::Real=0.4, mmax::Real=Inf)
The IMF model of [Salpeter 1955](https://ui.adsabs.harvard.edu/abs/1955ApJ...121..161S/abstract), a [`PowerLawIMF`](@ref) with `α=2.35`.
"""
Salpeter1955(mmin::T, mmax::T) where T <: Real= PowerLawIMF(T(2.35), mmin, mmax)
Salpeter1955(mmin::Real=0.4, mmax::Real=Inf) = Salpeter1955(promote(mmin, mmax)...)

###########################################################################################
# Broken Power Law
###########################################################################################
"""
    pl_integral(A,α,b1,b2)

Definite integral of power law ``A*x^{-α}`` from b1 (lower) to b2 (upper).

```math
\\int_{b1}^{b2} \\, A \\times x^{-\\alpha} \\, dx = \\frac{A}{1-\\alpha} \\times \\left( b2^{1-\\alpha} - b1^{1-\\alpha} \\right)
```
"""
pl_integral(A, α, b1, b2) = A/(1-α) * (b2^(1-α) - b1^(1-α)) 
"""
    BrokenPowerLaw(α::AbstractVector{T}, breakpoints::AbstractVector{S}) where {T<:Real,S<:Real}
    BrokenPowerLaw(α::Tuple, breakpoints::Tuple)
    BrokenPowerLaw{T}(A::Vector{T}, α::Vector{T}, breakpoints::Vector{T}) where {T}

An `AbstractIMF <: Distributions.ContinuousUnivariateDistribution` that describes a broken power-law IMF with probability distribution

```math
    \\frac{dn(m)}{dm} = A \\times m^{-\\alpha}
```

that is defined piecewise with different normalizations `A` and power law slopes `α` in different mass ranges. The normalization constants `A` will be calculated automatically after you provide the power law slopes and break points.

# Arguments
 - `α`; the power-law slopes of the different segments of the broken power law.
 - `breakpoints`; the masses at which the power law slopes change. If `length(α)=n`, then `length(breakpoints)=n+1`.

# Examples
`BrokenPowerLaw([1.35,2.35],[0.08,1.0,Inf])` will instantiate a broken power law defined from a minimum mass of `0.08` to a maximum mass of `Inf` with a single switch in `α` at `m=1.0`. From `0.08 ≤ m ≤ 1.0`, `α = 1.35` and from `1.0 ≤ m ≤ Inf`, `α = 2.35`.

# Notes
There is some setup necessary for `quantile` and other derived methods, so it is more efficient to call these methods directly with an array via the call signature `quantile(d::BrokenPowerLaw{T}, x::AbstractArray{S})` rather than broadcasting over `x`. This behavior is now deprecated for `quantile(d::Distributions.UnivariateDistribution, X::AbstractArray)` in Distributions.jl.

# Methods
 - `Base.convert(::Type{BrokenPowerLaw{T}}, d::BrokenPowerLaw)`
 - `minimum(d::BrokenPowerLaw)`
 - `maximum(d::BrokenPowerLaw)`
 - `partype(d::BrokenPowerLaw)`
 - `eltype(d::BrokenPowerLaw)`
 - `mean(d::BrokenPowerLaw)`
 - `median(d::BrokenPowerLaw)`
 - `var(d::BrokenPowerLaw)`, may not function correctly for large `mmax`
 - `skewness(d::BrokenPowerLaw)`, may not function correctly for large `mmax`
 - `kurtosis(d::BrokenPowerLaw)`, may not function correctly for large `mmax`
 - `pdf(d::BrokenPowerLaw,x::Real)`
 - `logpdf(d::BrokenPowerLaw,x::Real)`
 - `cdf(d::BrokenPowerLaw,x::Real)`
 - `ccdf(d::BrokenPowerLaw,x::Real)`
 - `quantile(d::BrokenPowerLaw{S},x::T) where {S,T<:Real}`
 - `quantile!(result::AbstractArray,d::BrokenPowerLaw{S},x::AbstractArray{T}) where {S,T<:Real}`
 - `quantile(d::BrokenPowerLaw{T},x::AbstractArray{S})`
 - `cquantile(d::BrokenPowerLaw{S},x::T) where {S,T<:Real}`
 - `rand(rng::AbstractRNG, d::BrokenPowerLaw,s...)` 
 - Other methods from `Distributions.jl` should also work because `BrokenPowerLaw <: AbstractIMF <: Distributions.ContinuousUnivariateDistribution`. For example, `rand!(rng::AbstractRNG, d::BrokenPowerLaw, x::AbstractArray)`.
"""
struct BrokenPowerLaw{T} <: AbstractIMF
    A::Vector{T}      # normalization parameters
    α::Vector{T}      # power law indexes
    breakpoints::Vector{T} # bounds of each break
    BrokenPowerLaw{T}(A::Vector{T}, α::Vector{T}, breakpoints::Vector{T}) where {T} =
        new{T}(A,α,breakpoints)
end
function BrokenPowerLaw(α::Vector{T}, breakpoints::Vector{T}) where T <: Real
    @assert length(breakpoints) == length(α)+1
    @assert breakpoints[1] > 0
    nbreaks = length(α)
    A = Vector{T}(undef, nbreaks)
    A[1] = one(T)
    for i in 2:nbreaks
        A[i] = A[i-1] * breakpoints[i]^-α[i-1] / breakpoints[i]^-α[i]
    end
    # Now A contains prefactors for each power law that makes them continuous
    # with the first power law having a prefactor of 1. Now we need to normalize
    # the integral from minimum(breaks) to maximum(breaks) to equal 1 by dividing
    # the entire A array by a common factor.
    total_integral = zero(T)
    for i in 1:nbreaks
        total_integral += pl_integral(A[i], α[i], breakpoints[i], breakpoints[i+1])
    end
    A ./= total_integral
    return BrokenPowerLaw{T}(A, α, breakpoints)
end
BrokenPowerLaw(α::Tuple, breakpoints::Tuple) = BrokenPowerLaw(collect(promote(α...)), collect(promote(breakpoints...)))
BrokenPowerLaw(α::AbstractVector{T}, breakpoints::AbstractVector{T}) where T <: Real =
    BrokenPowerLaw(convert(Vector{T}, α), convert(Vector{T}, breakpoints))
function BrokenPowerLaw(α::AbstractVector{T}, breakpoints::AbstractVector{S}) where {T <: Real, S <: Real}
    X = promote_type(T, S)
    return BrokenPowerLaw(convert(Vector{X},α), convert(Vector{X},breakpoints))
end

#### Conversions
Base.convert(::Type{BrokenPowerLaw{T}}, d::BrokenPowerLaw) where T =
    BrokenPowerLaw{T}(convert(Vector{T},d.A), convert(Vector{T},d.α), convert(Vector{T},d.breakpoints))
Base.convert(::Type{BrokenPowerLaw{T}}, d::BrokenPowerLaw{T}) where T <: Real = d

#### Parameters
params(d::BrokenPowerLaw) = d.A, d.α, d.breakpoints
minimum(d::BrokenPowerLaw) = minimum(d.breakpoints)
maximum(d::BrokenPowerLaw) = maximum(d.breakpoints)
partype(d::BrokenPowerLaw{T}) where T = T
eltype(d::BrokenPowerLaw{T}) where T = T

#### Statistics
function mean(d::BrokenPowerLaw)
    A, α, breakpoints = params(d)
    return sum( (A[i]*breakpoints[i+1]^(2-α[i])/(2-α[i]) -
        A[i]*breakpoints[i]^(2-α[i])/(2-α[i]) for i in 1:length(A) ) )
end
median(d::BrokenPowerLaw{T}) where T = quantile(d, T(0.5)) # this is temporary
# mode(d::BrokenPowerLaw) = d.breakpoints[argmin(d.α)] # this is not always correct

function var(d::BrokenPowerLaw)
    A, α, breakpoints = params(d)
    return sum( (A[i]*breakpoints[i+1]^(3-α[i])/(3-α[i]) -
        A[i]*breakpoints[i]^(3-α[i])/(3-α[i]) for i in 1:length(A) ) ) - mean(d)^2
end

function skewness(d::BrokenPowerLaw)
    A, α, breakpoints = params(d)
    m = mean(d)
    v = var(d)
    return ( sum( (A[i]*breakpoints[i+1]^(4-α[i])/(4-α[i]) -
                   A[i]*breakpoints[i]^(4-α[i])/(4-α[i]) for i in 1:length(A) ) )
                   - 3 * m * v - m^3) / v^(3/2)
end

function kurtosis(d::BrokenPowerLaw)
    A, α, breakpoints = params(d)
    X4 = sum( (A[i]*breakpoints[i+1]^(5-α[i])/(5-α[i]) -
        A[i]*breakpoints[i]^(5-α[i])/(5-α[i]) for i in 1:length(A) ) )
    X3 = sum( (A[i]*breakpoints[i+1]^(4-α[i])/(4-α[i]) -
        A[i]*breakpoints[i]^(4-α[i])/(4-α[i]) for i in 1:length(A) ) )
    X2 = sum( (A[i]*breakpoints[i+1]^(3-α[i])/(3-α[i]) -
        A[i]*breakpoints[i]^(3-α[i])/(3-α[i]) for i in 1:length(A) ) )
    X = sum( (A[i]*breakpoints[i+1]^(2-α[i])/(2-α[i]) -
        A[i]*breakpoints[i]^(2-α[i])/(2-α[i]) for i in 1:length(A) ) )
    μ4 = X4 - 4*X*X3 + 6*X^2*X2 - 3*X^4
    return μ4 / var(d)^2
end

#### Evaluation
function pdf(d::BrokenPowerLaw{S}, x::T) where {S, T <: Real}
    ((x < minimum(d)) || (x > maximum(d))) && (return zero(promote_type(S,T)))
    idx = findfirst(>=(x), d.breakpoints)
    idx != 1 && (idx-=1)
    return d.A[idx] * x^-d.α[idx]
end
function logpdf(d::BrokenPowerLaw{S}, x::T) where {S, T <: Real}
    if ((x >= minimum(d)) && (x <= maximum(d)))
        A, α, breakpoints = params(d)
        idx = findfirst(>=(x), breakpoints)
        idx != 1 && (idx-=1)
        return log(A[idx]) - α[idx]*log(x)
    else
        U = promote_type(S, T)
        return -U(Inf)
    end
end
function cdf(d::BrokenPowerLaw{S}, x::T) where {S, T <: Real}
    U = promote_type(S, T)
    if x <= minimum(d)
        return zero(U) 
    elseif x >= maximum(d)
        return one(U)
    end
    A,α,breakpoints = params(d)
    idx = findfirst(>=(x),breakpoints)
    idx != 1 && (idx-=1)
    return sum(pl_integral(A[i],α[i],breakpoints[i],min(x,breakpoints[i+1])) for i in 1:idx)
end
ccdf(d::BrokenPowerLaw, x::Real) = one(partype(d)) - cdf(d,x)
function quantile(d::BrokenPowerLaw{S}, x::T) where {S, T <: Real}
    U = promote_type(S, T)
    if x <= zero(T)
        return U(minimum(d))
    elseif x >= one(T)
        return U(maximum(d))
    end
    A, α, breakpoints = params(d)
    nbreaks = length(A)
    integrals = cumsum(pl_integral(A[i],α[i],breakpoints[i],breakpoints[i+1]) for i in 1:nbreaks) # calculate the cumulative integral 
    idx = findfirst(>=(x), integrals)  # find the first breakpoint where the cumulative integral   # up to each breakpoint 
    idx != 1 && (x-=integrals[idx-1]) # is greater than x. If this is not the first breakpoint, then subtract off the cumulative integral
    a = one(S) - α[idx]
    return (x*a/A[idx] + breakpoints[idx]^a)^inv(a)
end
function quantile!(result::AbstractArray, d::BrokenPowerLaw{S}, x::AbstractArray{T}) where {S, T <: Real}
    @assert axes(result) == axes(x)
    A, α, breakpoints = params(d)
    nbreaks = length(A)
    integrals = cumsum(pl_integral(A[i],α[i],breakpoints[i],breakpoints[i+1]) for i in 1:nbreaks) # calculate the cumulative integral
    @inbounds for i in eachindex(x)
        xi = x[i]
        xi<=zero(T) && (result[i]=minimum(d); continue)
        xi>=one(T) && (result[i]=maximum(d); continue)
        idx = findfirst(>=(xi), integrals)  # find the first breakpoint where the cumulative integral   # up to each breakpoint 
        idx != 1 && (xi-=integrals[idx-1]) # is greater than x. If this is not the first breakpoint, then subtract off 
        a = one(S) - α[idx]                # the cumulative integral
        result[i] = (xi*a/A[idx] + breakpoints[idx]^a)^inv(a)
    end
    return result
end
quantile(d::BrokenPowerLaw{T}, x::AbstractArray{S}) where {T, S <: Real} =
    quantile!(Array{promote_type(T,S)}(undef,size(x)),d,x)
cquantile(d::BrokenPowerLaw, x::Real) = quantile(d, one(x)-x)

##### Random sampling
##########################################################################################
# Implementing efficient sampler. We need the cumulative integral up to each breakpoint in
# order to transform a uniform random point in (0,1] to a random mass, but we dont want to
# recompute it every time. We could just use the method of
# quantile(d::BrokenPowerLaw, x::AbstractArray)
# for this purpose but the "correct" way to do it is to implement a `sampler` method as below.
# By default (e.g., without `rand(rng::AbstractRNG, d::BrokenPowerLaw)`),
# Distributions seems to call the efficient
# `quantile(d::BrokenPowerLaw, x::AbstractArray)` method anyway, but it doesn't respect
# the type of `d`, so a bit better to do it this way anyway.
struct BPLSampler{T} <: Sampleable{Univariate, Continuous}
    A::Vector{T}           # normalization parameters
    α::Vector{T}           # power law indexes
    breakpoints::Vector{T} # bounds of each break
    integrals::Vector{T}   # cumulative integral up to each breakpoint
end
function BPLSampler(d::BrokenPowerLaw)
    A, α, breakpoints = params(d)
    nbreaks = length(A)
    integrals = cumsum(pl_integral(A[i],α[i],breakpoints[i],breakpoints[i+1]) for i in 1:nbreaks) # calculate the cumulative integral
    return BPLSampler(A, α, breakpoints, integrals)
end
function rand(rng::AbstractRNG, s::BPLSampler{T}) where T
    x = rand(rng, T)
    A, α, breakpoints,integrals = s.A,s.α,s.breakpoints,s.integrals
    idx = findfirst(>=(x), integrals)  # find the first breakpoint where the cumulative integral   # up to each breakpoint 
    idx != 1 && (x-=integrals[idx-1]) # is greater than x. If this is not the first breakpoint, then subtract off the cumulative integral
    a = one(T) - α[idx]
    return (x*a/A[idx] + breakpoints[idx]^a)^inv(a)
end
sampler(d::BrokenPowerLaw) = BPLSampler(d)
rand(rng::AbstractRNG, d::BrokenPowerLaw) = rand(rng, sampler(d))

#######################################################
# Specific types of BrokenPowerLaw
#######################################################

const kroupa2001_α = [0.3, 1.3, 2.3]
const kroupa2001_breakpoints = [0.0, 0.08, 0.50, Inf]
"""
    Kroupa2001(mmin::Real=0.08, mmax::Real=Inf)

Function to instantiate a [`BrokenPowerLaw`](@ref) IMF with the parameters from Equation 2 of [Kroupa 2001](https://ui.adsabs.harvard.edu/abs/2001MNRAS.322..231K/abstract). This is equivalent to the relation given in [Kroupa 2002](https://ui.adsabs.harvard.edu/abs/2002Sci...295...82K/abstract).
"""
function Kroupa2001(mmin::T=0.08, mmax::T=Inf) where T <: Real
    @assert mmin > zero(T)
    # idx1 = max(1, findfirst(>(mmin),kroupa2001_breakpoints)-1)
    idx1 = findfirst(>(mmin), kroupa2001_breakpoints) - 1
    idx2 = findfirst(>=(mmax), kroupa2001_breakpoints)
    bp = convert(Vector{T}, kroupa2001_breakpoints[idx1:idx2])
    bp[1] = mmin
    bp[end] = mmax
    BrokenPowerLaw(convert(Vector{T}, kroupa2001_α[idx1:idx2-1]), bp)
end
Kroupa2001(mmin::Real, mmax::Real) = Kroupa2001(promote(mmin,mmax)...)

const chabrier2001bpl_α = [1.55, 2.70]
const chabrier2001bpl_breakpoints = [0.00, 1.0, Inf]
"""
    Chabrier2001BPL(mmin::T=0.08, mmax::T=Inf)

Function to instantiate a [`BrokenPowerLaw`](@ref) IMF with the parameters from the first column of Table 1 in [Chabrier 2001](https://ui.adsabs.harvard.edu/abs/2001ApJ...554.1274C/abstract).
"""
function Chabrier2001BPL(mmin::T=0.08, mmax::T=Inf) where {T<:Real}
    @assert mmin > 0
    # idx1 = max(1, findfirst(>(mmin),kroupa2001_breakpoints)-1)
    idx1 = findfirst(>(mmin), chabrier2001bpl_breakpoints) - 1
    idx2 = findfirst(>=(mmax), chabrier2001bpl_breakpoints)
    bp = convert(Vector{T}, chabrier2001bpl_breakpoints[idx1:idx2])
    bp[1] = mmin
    bp[end] = mmax
    BrokenPowerLaw(convert(Vector{T}, chabrier2001bpl_α[idx1:idx2-1]), bp)
end
Chabrier2001BPL(mmin::Real, mmax::Real) = Chabrier2001BPL(promote(mmin,mmax)...)
