PowerLawIMF(α::Real,mmin::Real,mmax::Real) = truncated(Pareto(α-1,mmin);upper=mmax)
Salpeter1955(mmin::Real=0.4,mmax::Real=10.0) = PowerLawIMF(2.35,mmin,mmax)
Chabrier2001Lognormal(mmin::Real=0.08,mmax::Real=120.0) = truncated(LogNormal(-log(10),0.627*log(10)),mmin,mmax)

struct Chabrier2003{T<:Real} <: AbstractIMF
    A1::T # lognormal parameter
    mc::T # lognormal parameter
    σ::T  # lognormal parameter
    A2::T # power law parameter for large M
    x::T  # power law parameter for large M
end
Chabrier2003(A1::Real=0.158,mc::Real=0.079,σ::Real=0.69,A2::Real=0.0443,x::Real=1.35) = Chabrier2003(promote(A1,mc,σ,A2,x)...)
# these are not pdfs; they are not normalized to integrate to one. They are normalized to equal the hipparcos value at 0.8 solar masses
function logpdf(d::Chabrier2003,m::Real)
    if m <= one(m)
        log(d.A1) - log(m) - log(log(10)) - 0.5*(log10(m)-log10(d.mc))^2/d.σ^2
    else
        log(d.A2) - log(m) - log(log(10)) -d.x * log(m)
    end
end
function pdf(d::Chabrier2003,m::Real)
    if m <= one(m)
        d.A1 / m / log(10) * exp(-0.5*(log10(m)-log10(d.mc))^2/d.σ^2)
    else
        d.A2 / m / log(10) * m^-d.x
    end
end

###########################################################################################
# Broken Power Law 
pl_integral(A,α,b1,b2) = A/(1-α) * (b2^(1-α) - b1^(1-α)) #definite integral of power law A*x^-α from b1 (lower) to b2 (upper)

struct BrokenPowerLaw{T} <: AbstractIMF
    A::Vector{T}      # normalization parameters
    α::Vector{T}      # power law indexes
    breakpoints::Vector{T} # bounds of each break
    BrokenPowerLaw{T}(A::Vector{T},α::Vector{T},breakpoints::Vector{T}) where {T} = new{T}(A,α,breakpoints)
end
function BrokenPowerLaw(α::Vector{T},breakpoints::Vector{T}) where T<:Real
    @assert length(breakpoints) == length(α)+1
    @assert breakpoints[1]>0
    nbreaks = length(α)
    A=Vector{T}(undef,nbreaks)
    A[1]=one(T)
    @inbounds for i in 2:nbreaks
        A[i] = A[i-1]*breakpoints[i]^-α[i-1] / breakpoints[i]^-α[i]
    end
    # now A contains prefactors for each power law that makes them continuous
    # with the first power law having a prefactor of 1. Now we need to normalize
    # the integral from minimum(breaks) to maximum(breaks) to equal 1 by dividing
    # the entire A array by a common factor.
    total_integral = zero(T)
    for i in 1:nbreaks
        total_integral += pl_integral(A[i],α[i],breakpoints[i],breakpoints[i+1])
    end
    A ./= total_integral
    return BrokenPowerLaw{T}(A,α,breakpoints)
end
BrokenPowerLaw(α::Tuple,breakpoints::Tuple) = BrokenPowerLaw(collect(promote(α...)),collect(promote(breakpoints...)))
BrokenPowerLaw(α::AbstractVector{T},breakpoints::AbstractVector{T}) where {T<:Real} = BrokenPowerLaw(convert(Vector{T},α),convert(Vector{T},breakpoints))
function BrokenPowerLaw(α::AbstractVector{T},breakpoints::AbstractVector{S}) where {T<:Real,S<:Real}
    X = promote_type(T,S)
    X == T ? BrokenPowerLaw(α,convert(Vector{T},breakpoints)) : BrokenPowerLaw(convert(Vector{S},α),breakpoints)
end
# struct BrokenPowerLaw{T,N,S} <: AbstractIMF
#     A::MVector{N,T}      # normalization parameters
#     α::SVector{N,T}      # power law indexes
#     breakpoints::SVector{S,T} # bounds of each break
#     BrokenPowerLaw{T,N,S}(A::MVector{N,T},α::SVector{N,T},breakpoints::SVector{S,T}) where {T,N,S} = new{T,N,S}(A,α,breakpoints)
# end
# function BrokenPowerLaw(α::SVector{N,T},breakpoints::SVector{S,T}) where {T<:Real,N,S}
#     @assert length(breakpoints) == length(α)+1
#     @assert breakpoints[1]>0
#     # nbreaks = length(α)
#     A = MVector{N,T}(undef)
#     A[1]=1
#     @inbounds for i in 2:N
#         A[i] = A[i-1]*breakpoints[i]^-α[i-1] / breakpoints[i]^-α[i]
#     end
#     # now A contains prefactors for each power law that makes them continuous
#     # with the first power law having a prefactor of 1. Now we need to normalize
#     # the integral from minimum(breaks) to maximum(breaks) to equal 1 by dividing
#     # the entire A array by a common factor.
#     total_integral = zero(T)
#     for i in 1:N
#         total_integral += pl_integral(A[i],α[i],breakpoints[i],breakpoints[i+1])
#     end
#     A ./= total_integral
#     return BrokenPowerLaw{T,N,S}(A,α,breakpoints)
# end

#### Conversions
# convert(::Type{BrokenPowerLaw{T}}, α::Real, θ::Real) where {T<:Real} = Pareto(T(α), T(θ))
Base.convert(::Type{BrokenPowerLaw{T}}, d::BrokenPowerLaw) where T = BrokenPowerLaw{T}(convert(Vector{T},d.A), convert(Vector{T},d.α), convert(Vector{T},d.breakpoints) )
Base.convert(::Type{BrokenPowerLaw{T}}, d::BrokenPowerLaw{T}) where {T<:Real} = d

params(d::BrokenPowerLaw) = d.A,d.α,d.breakpoints
minimum(d::BrokenPowerLaw) = minimum(d.breakpoints)
maximum(d::BrokenPowerLaw) = maximum(d.breakpoints)
partype(d::BrokenPowerLaw{T}) where T = T
function pdf(d::BrokenPowerLaw,x::Real)
    ((x < minimum(d)) || (x > maximum(d))) && (return zero(partype(d)))
    idx = findfirst(>=(x),d.breakpoints)
    idx != 1 && (idx-=1)
    d.A[idx] * x^-d.α[idx]
end
function logpdf(d::BrokenPowerLaw,x::Real)
    if ((x >= minimum(d)) && (x <= maximum(d)))
        A,α,breakpoints = params(d)
        idx = findfirst(>=(x),breakpoints)
        idx != 1 && (idx-=1)
        log(A[idx]) - α[idx]*log(x)
    else
        T = partype(d)
        return -T(Inf)
    end
end
function cdf(d::BrokenPowerLaw,x::Real)
    if x <= minimum(d)
        return zero(partype(d)) 
    elseif x >= maximum(d)
        return one(partype(d))
    end
    A,α,breakpoints = params(d)
    idx = findfirst(>=(x),d.breakpoints)
    idx != 1 && (idx-=1)
    result = sum(pl_integral(A[i],α[i],breakpoints[i],min(x,breakpoints[i+1])) for i in 1:idx)
end
ccdf(d::BrokenPowerLaw,x::Real) = one(partype(d)) - cdf(d,x)
function quantile(d::BrokenPowerLaw{S},x::T) where {S,T<:Real}
    U = promote_type(S,T)
    x<=zero(T) && (return U(minimum(d)))
    x>=one(T) && (return U(maximum(d)))
    A,α,breakpoints = params(d)
    nbreaks = length(A)
    integrals = cumsum(pl_integral(A[i],α[i],breakpoints[i],breakpoints[i+1]) for i in 1:nbreaks) # calculate the cumulative integral 
    idx = findfirst(>=(x),integrals)  # find the first breakpoint where the cumulative integral   # up to each breakpoint 
    idx != 1 && (x-=integrals[idx-1]) # is greater than x. If this is not the first breakpoint, then subtract off the cumulative integral
    a = one(S) - α[idx]
    return (x*a/A[idx] + breakpoints[idx]^a)^inv(a)
end
function quantile!(result::AbstractArray,d::BrokenPowerLaw{S},x::AbstractArray{T}) where {S,T<:Real}
    @assert size(result) == size(x)
    A,α,breakpoints = params(d)
    nbreaks = length(A)
    integrals = cumsum(pl_integral(A[i],α[i],breakpoints[i],breakpoints[i+1]) for i in 1:nbreaks) # calculate the cumulative integral
    @inbounds for i in eachindex(x)
        xi = x[i]
        xi<=zero(T) && (result[i]=minimum(d); continue)
        xi>=one(T) && (result[i]=maximum(d); continue)
        idx = findfirst(>=(xi),integrals)  # find the first breakpoint where the cumulative integral   # up to each breakpoint 
        idx != 1 && (xi-=integrals[idx-1]) # is greater than x. If this is not the first breakpoint, then subtract off 
        a = one(S) - α[idx]                # the cumulative integral
        result[i] = (xi*a/A[idx] + breakpoints[idx]^a)^inv(a)
    end
    return result
end
quantile(d::BrokenPowerLaw{T},x::AbstractArray{S}) where {T,S<:Real} = quantile!(Array{promote_type(T,S)}(undef,size(x)),d,x)

##### Random sampling

# Implementing efficient sampler. We need the cumulative integral up to each breakpoint in
# order to transform a uniform random point in (0,1] to a random mass, but we dont want to
# recompute it every time. We could just use the method of quantile(d::BrokenPowerLaw,x::AbstractArray)
# for this purpose but the "correct" way to do it is to implement a `sampler` method as below.
# By default (e.g., without `rand(rng::AbstractRNG, d::BrokenPowerLaw)`), Distributions seems to call
# the efficient `quantile(d::BrokenPowerLaw,x::AbstractArray)` method anyway, but it doesn't respect
# the type of `d`, so a bit better to do it this way anyway.
struct BPLSampler{T} <: Sampleable{Univariate,Continuous}
    A::Vector{T}           # normalization parameters
    α::Vector{T}           # power law indexes
    breakpoints::Vector{T} # bounds of each break
    integrals::Vector{T}   # cumulative integral up to each breakpoint
end
function BPLSampler(d::BrokenPowerLaw)
    A,α,breakpoints = params(d)
    nbreaks = length(A)
    integrals = cumsum(pl_integral(A[i],α[i],breakpoints[i],breakpoints[i+1]) for i in 1:nbreaks) # calculate the cumulative integral
    BPLSampler(A,α,breakpoints,integrals)
end
function rand(rng::AbstractRNG, s::BPLSampler{T}) where T
    x = rand(rng,T)
    A,α,breakpoints,integrals = s.A,s.α,s.breakpoints,s.integrals
    idx = findfirst(>=(x),integrals)  # find the first breakpoint where the cumulative integral   # up to each breakpoint 
    idx != 1 && (x-=integrals[idx-1]) # is greater than x. If this is not the first breakpoint, then subtract off the cumulative integral
    a = one(T) - α[idx]
    return (x*a/A[idx] + breakpoints[idx]^a)^inv(a)
end
sampler(d::BrokenPowerLaw) = BPLSampler(d)
rand(rng::AbstractRNG, d::BrokenPowerLaw) = rand(rng, sampler(d))

#######################################################
# specific types of BrokenPowerLaw
const kroupa2001_α = [0.3, 1.3, 2.3]
const kroupa2001_breakpoints = [0.01,0.08,0.50,Inf]
"""
    Kroupa2001(mmin::Real=0.01,mmax::Real=Inf)

Function to instantiate a `BrokenPowerLaw` IMF with the parameters from Equation 2 of Kroupa 2001. This IMF is defined from `M=0.01` solar masses to `M=Inf`. This is equivalent to the relation given in Kroupa 2002, which appears to be more highly cited.
"""
function Kroupa2001(mmin::T=0.01,mmax::T=Inf) where {T<:Real}
    @assert mmin>=0.01
    # idx1 = max(1, findfirst(>(mmin),kroupa2001_breakpoints)-1)
    idx1 = findfirst(>(mmin),kroupa2001_breakpoints)-1
    idx2 = findfirst(>=(mmax),kroupa2001_breakpoints)
    bp = convert(Vector{T},kroupa2001_breakpoints[idx1:idx2])
    bp[1] = mmin
    bp[end] = mmax
    BrokenPowerLaw(convert(Vector{T},kroupa2001_α[idx1:idx2-1]),bp)
end
Kroupa2001(mmin::Real,mmax::Real) = Kroup2001(promote(mmin,mmax)...)

const chabrier2001bpl_α = [1.55,2.70]
const chabrier2001bpl_breakpoints = [0.07,1.0,100.0]
"""
    Chabrier2001BPL(mmin::T=0.08,mmax::T=Inf)

Function to instantiate a `BrokenPowerLaw` IMF with the parameters from the first column of Table 1 in Chabrier 2001. This IMF is defined from `M=0.07` solar masses to `M=100`.
"""
function Chabrier2001BPL(mmin::T=0.08,mmax::T=Inf) where {T<:Real}
    @assert mmin>=0.01
    # idx1 = max(1, findfirst(>(mmin),kroupa2001_breakpoints)-1)
    idx1 = findfirst(>(mmin),chabrier2001bpl_breakpoints)-1
    idx2 = findfirst(>=(mmax),chabrier2001bpl_breakpoints)
    bp = convert(Vector{T},chabrier2001bpl_breakpoints[idx1:idx2])
    bp[1] = mmin
    bp[end] = mmax
    BrokenPowerLaw(convert(Vector{T},chabrier2001bpl_α[idx1:idx2-1]),bp)
end
Chabrier2001BPL(mmin::Real,mmax::Real) = Chabrier2001BPL(promote(mmin,mmax)...)
