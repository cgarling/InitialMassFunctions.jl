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
pl_integral(A,α,b1,b2) = A/(1-α) * (b2^(1-α) - b1^(1-α)) #definite integral of power law A*x^-α from b1 (lower) to b2 (upper)

struct BrokenPowerLaw{T,S} <: AbstractIMF
    A::T      # normalization parameters
    α::T      # power law indexes
    breakpoints::S # bounds of each break
    BrokenPowerLaw{T,S}(A::T,α::T,breakpoints::S) where {T,S} = new{T,S}(A,α,breakpoints)
end
function BrokenPowerLaw(α::T,breakpoints::S) where {T,S}
    @assert length(breakpoints) == length(α)+1
    @assert breakpoints[1]>0
    nbreaks = length(α)
    U = eltype(α)
    A=Array{U}(undef,nbreaks)
    A[1]=1
    @inbounds for i in 2:nbreaks
        A[i] = A[i-1]*breakpoints[i]^-α[i-1] / breakpoints[i]^-α[i]
    end
    # now A contains prefactors for each power law that makes them continuous
    # with the first power law having a prefactor of 1. Now we need to normalize
    # the integral from minimum(breaks) to maximum(breaks) to equal 1 by dividing
    # the entire A array by a common factor.
    total_integral = zero(U)
    for i in 1:nbreaks
        total_integral += pl_integral(A[i],α[i],breakpoints[i],breakpoints[i+1])
    end
    A ./= total_integral
    return BrokenPowerLaw{T,S}(A,α,breakpoints)
end
params(d::BrokenPowerLaw) = d.A,d.α,d.breakpoints
minimum(d::BrokenPowerLaw) = minimum(d.breakpoints)
maximum(d::BrokenPowerLaw) = maximum(d.breakpoints)
partype(d::BrokenPowerLaw{T}) where T = T
Base.eltype(d::BrokenPowerLaw{T}) where T = eltype(T)
function pdf(d::BrokenPowerLaw,x::Real)
    ((x < minimum(d)) || (x > maximum(d))) && (return zero(eltype(d)))
    idx = findfirst(>=(x),d.breakpoints)
    idx != 1 && (idx-=1)
    d.A[idx] * x^-d.α[idx]
end
function logpdf(d::BrokenPowerLaw,x::Real)
    T = eltype(d)
    if ((x >= minimum(d)) && (x <= maximum(d)))
        A,α,breakpoints = params(d)
        idx = findfirst(>=(x),breakpoints)
        idx != 1 && (idx-=1)
        log(A[idx]) - α[idx]*log(x)
    else
        return -T(Inf)
    end
end
function cdf(d::BrokenPowerLaw,x::Real)
    if x <= minimum(d)
        return zero(eltype(d))
    elseif x >= maximum(d)
        return one(eltype(d))
    end
    A,α,breakpoints = params(d)
    idx = findfirst(>=(x),d.breakpoints)
    idx != 1 && (idx-=1)
    result = sum(pl_integral(A[i],α[i],breakpoints[i],min(x,breakpoints[i+1])) for i in 1:idx)
end
ccdf(d::BrokenPowerLaw,x::Real) = 1 - cdf(d,x)
function quantile(d::BrokenPowerLaw,x::Real)
    x<=zero(x) && (return minimum(d))
    x>=one(x) && (return maximum(d))
    A,α,breakpoints = params(d)
    nbreaks = length(A)
    integrals = cumsum(pl_integral(A[i],α[i],breakpoints[i],breakpoints[i+1]) for i in 1:nbreaks) # calculate the cumulative integral 
    idx = findfirst(>=(x),integrals)  # find the first breakpoint where the cumulative integral   # up to each breakpoint 
    idx != 1 && (x-=integrals[idx-1]) # is greater than x. If this is not the first breakpoint, then subtract off the cumulative integral
    # if idx == 1                     # up to the last breakpoitn
    #     a = 1 - α[idx]
    #     return (x*a/A[idx] + breakpoints[idx]^a)^(1/a)
    # else
    #     return
    # end
    a = 1 - α[idx]
    return (x*a/A[idx] + breakpoints[idx]^a)^(1/a)
end
function quantile!(result::AbstractArray,d::BrokenPowerLaw,x::AbstractArray{T}) where {T<:Real}
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
        a = 1 - α[idx]                     # the cumulative integral
        result[i] = (xi*a/A[idx] + breakpoints[idx]^a)^(1/a)
    end
    return result
end
quantile(d::BrokenPowerLaw,x::AbstractArray) = quantile!(Array{promote_type(eltype(d),eltype(x))}(undef,size(x)),d,x)

##### Random sampling

# Implementing efficient sampler. We need the cumulative integral up to each breakpoint in
# order to transform a uniform random point in (0,1] to a random mass, but we dont want to
# recompute it every time. We could just use the method of quantile(d::BrokenPowerLaw,x::AbstractArray)
# for this purpose but the "correct" way to do it is to implement a `sampler` method as below. 
# rand(rng::AbstractRNG,d::BrokenPowerLaw) = quantile(rand())
struct BPLSampler{T,S} <: Sampleable{Univariate,Continuous}
    A::T      # normalization parameters
    α::T      # power law indexes
    breakpoints::S # bounds of each break
    integrals::T   # cumulative integral up to each breakpoint
end
function BPLSampler(d::BrokenPowerLaw)
    A,α,breakpoints = params(d)
    nbreaks = length(A)
    integrals = cumsum(pl_integral(A[i],α[i],breakpoints[i],breakpoints[i+1]) for i in 1:nbreaks) # calculate the cumulative integral
    BPLSampler(A,α,breakpoints,integrals)
end
function rand(rng::AbstractRNG, s::BPLSampler)
    x = rand(rng)
    A,α,breakpoints,integrals = s.A,s.α,s.breakpoints,s.integrals
    idx = findfirst(>=(x),integrals)  # find the first breakpoint where the cumulative integral   # up to each breakpoint 
    idx != 1 && (x-=integrals[idx-1]) # is greater than x. If this is not the first breakpoint, then subtract off the cumulative integral
    a = 1 - α[idx]
    return (x*a/A[idx] + breakpoints[idx]^a)^(1/a)
end
sampler(d::BrokenPowerLaw) = BPLSampler(d)
