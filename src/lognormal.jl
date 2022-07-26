"""
    Chabrier2001LogNormal(mmin::Real=0.08,mmax::Real=Inf)

Function to instantiate the [Chabrier 2001](https://ui.adsabs.harvard.edu/abs/2001ApJ...554.1274C/abstract) lognormal IMF. Returns an instance of `Distributions.Truncated(Distributions.LogNormal)`. See also [`Chabrier2003`](@ref) which has the same lognormal form for masses below one solar mass, but a power law extension at higher masses. 
"""
Chabrier2001LogNormal(mmin::Real=0.08,mmax::Real=Inf) = truncated(LogNormal(log(0.1),0.627*log(10));lower=mmin,upper=mmax)
# Chabrier2001LogNormal(mmin::Real=0.08,mmax::Real=Inf) = truncated(LogNormal(log(0.1)*log(10),0.627*log(10));lower=mmin,upper=mmax)

"""
    lognormal_integral(μ,σ,b1,b2)

Definite integral of the lognormal probability distribution from `b1` to `b2`.
```math
\\int_{b1}^{b2} \\, \\frac{A}{x} \\, \\exp \\left[ \\frac{ -\\left( \\log(x) - \\mu \\right)^2}{2\\sigma^2} \\right] \\, dx
```
"""
lognormal_integral(A::T,μ::T,σ::T,b1::T,b2::T) where {T<:Number} = A * sqrt(T(π)/2) * σ * (erf( (μ-log(b1))/(sqrt(T(2))*σ)) - erf( (μ-log(b2))/(sqrt(T(2))*σ)))
lognormal_integral(A::Number,μ::Number,σ::Number,b1::Number,b2::Number) = lognormal_integral(promote(A,μ,σ,b1,b2)...)
# lognormal_integral(A,μ,σ,b1,b2) = A * sqrt(π/2) * σ * (erf( (μ-log(b1))/(sqrt(2)*σ)) - erf( (μ-log(b2))/(sqrt(2)*σ)))
"""
    LogNormalBPL(μ::Real,σ::Real,α::AbstractVector{<:Real},breakpoints::AbstractVector{<:Real})
    LogNormalBPL(μ::Real,σ::Real,α::Tuple,breakpoints::Tuple)

A LogNormal distribution at low masses, with a broken power law extension at high masses. This uses the natural log base like `Distributions.LogNormal`; if you have σ and μ in base 10, then multiply them both by `log(10)`. Must have `length(α) == length(breakpoints)-2`. The probability distribution for this IMF model is

```math
    \\frac{dn(m)}{dm} = \\frac{A}{x} \\, \\exp \\left[ \\frac{ -\\left( \\log(x) - \\mu \\right)^2}{2\\sigma^2} \\right]
```

# Examples
If you want a `LogNormalBPL` with a characteristic mass of 0.5 solar masses, log10 standard deviation of 0.6, and a single power law extension with slope `α=2.35` with a break at 1 solar mass, you would do `LogNormalBPL(log(0.5),0.6*log(10),[2.35],[0.08,1.0,Inf]` where we set the minimum mass to `0.08` and maximum mass to `Inf`. If, instead, you know that `log10(m)=x`, where `m` is the characteristic mass of the `LogNormal` component, you would do `LogNormalBPL(x*log(10),0.6*log(10),[2.35],[0.08,1.0,Inf]`
"""
struct LogNormalBPL{T} <: AbstractIMF
    μ::T
    σ::T
    A::Vector{T}      # normalization parameters
    α::Vector{T}      # power law indexes
    breakpoints::Vector{T} # bounds of each break
    LogNormalBPL{T}(μ::T,σ::T,A::Vector{T},α::Vector{T},breakpoints::Vector{T}) where {T} = new{T}(μ,σ,A,α,breakpoints)
end
function LogNormalBPL(μ::T,σ::T,α::Vector{T},breakpoints::Vector{T}) where T<:Real
    @assert length(breakpoints) == length(α)+2
    @assert breakpoints[1]>0
    nbreaks = length(α) + 1
    A=Vector{T}(undef,nbreaks)
    A[1]=one(T)
    # solve for the prefactor for the first power law after the lognormal component
    A[2]=breakpoints[2]^(α[1]-1) * exp( -(log(breakpoints[2])-μ)^2/(2*σ^2))
    # if there is more than one power law component
    if nbreaks>2
        @inbounds for i in 3:nbreaks
            A[i] = A[i-1]*breakpoints[i]^-α[i-2] / breakpoints[i]^-α[i-1]
        end
    end
    # now A contains prefactors for each distribution component that makes them continuous
    # with the first lognormal component having a prefactor of 1. Now we need to normalize
    # the integral from minimum(breaks) to maximum(breaks) to equal 1 by dividing
    # the entire A array by a common factor.
    total_integral = zero(T)
    total_integral += lognormal_integral(A[1],μ,σ,breakpoints[1],breakpoints[2])
    for i in 2:nbreaks
        total_integral += pl_integral(A[i],α[i-1],breakpoints[i],breakpoints[i+1])
    end
    A ./= total_integral
    return LogNormalBPL{T}(μ,σ,A,α,breakpoints)
end
LogNormalBPL(μ::Real,σ::Real,α::Tuple,breakpoints::Tuple) = LogNormalBPL(μ,σ,collect(promote(α...)),collect(promote(breakpoints...)))
LogNormalBPL(μ::T,σ::T,α::AbstractVector{T},breakpoints::AbstractVector{T}) where {T<:Real} = LogNormalBPL(μ,σ,convert(Vector{T},α),convert(Vector{T},breakpoints))
function LogNormalBPL(μ::A,σ::B,α::AbstractVector{C},breakpoints::AbstractVector{D}) where {A<:Real,B<:Real,C<:Real,D<:Real}
    X = promote_type(A,B,C,D)
    LogNormalBPL(convert(X,μ),convert(X,σ),convert(Vector{X},α), convert(Vector{X},breakpoints))
end

#### Conversions
Base.convert(::Type{LogNormalBPL{T}}, d::LogNormalBPL) where T = LogNormalBPL{T}(convert(T,d.μ),convert(T,d.σ),convert(Vector{T},d.A), convert(Vector{T},d.α), convert(Vector{T},d.breakpoints) )
Base.convert(::Type{LogNormalBPL{T}}, d::LogNormalBPL{T}) where {T<:Real} = d

#### Parameters
params(d::LogNormalBPL) = d.μ,d.σ,d.A,d.α,d.breakpoints
minimum(d::LogNormalBPL) = minimum(d.breakpoints)
maximum(d::LogNormalBPL) = maximum(d.breakpoints)
partype(d::LogNormalBPL{T}) where T = T

#### Statistics
function mean(d::LogNormalBPL{T}) where T
    μ,σ,A,α,breakpoints = params(d)
    m = zero(T)
    # m += A[1] * exp(μ + σ^2/2) * sqrt(π/2) * σ * (erf( (μ+σ^2-log(breakpoints[1])) / (sqrt(2)*σ) ) -
    #     erf( (μ+σ^2-log(breakpoints[2])) / (sqrt(2)*σ) ) )
    m += A[1] * exp(μ + σ^2/2) * sqrt(T(π)/2) * σ * (erf( (μ+σ^2-log(breakpoints[1])) / (sqrt(T(2))*σ) ) -
        erf( (μ+σ^2-log(breakpoints[2])) / (sqrt(T(2))*σ) ) )
    m += sum( (A[i]*breakpoints[i+1]^(2-α[i-1])/(2-α[i-1]) -
        A[i]*breakpoints[i]^(2-α[i-1])/(2-α[i-1]) for i in 2:length(A) ) )
    return m
end
median(d::LogNormalBPL{T}) where T = quantile(d,T(0.5)) # this is temporary
# mode(d::BrokenPowerLaw) = d.breakpoints[argmin(d.α)] # this is not always correct

# function var(d::BrokenPowerLaw)
#     A,α,breakpoints = params(d)
#     sum( (A[i]*breakpoints[i+1]^(3-α[i])/(3-α[i]) -
#         A[i]*breakpoints[i]^(3-α[i])/(3-α[i]) for i in 1:length(A) ) ) - mean(d)^2
# end

# function skewness(d::BrokenPowerLaw)
#     A,α,breakpoints = params(d)
#     m=mean(d)
#     v=var(d)
#     ( sum( (A[i]*breakpoints[i+1]^(4-α[i])/(4-α[i]) -
#         A[i]*breakpoints[i]^(4-α[i])/(4-α[i]) for i in 1:length(A) ) ) - 3 * m * v - m^3) / v^(3/2)
# end

# function kurtosis(d::BrokenPowerLaw)
#     A,α,breakpoints = params(d)
#     X4 = sum( (A[i]*breakpoints[i+1]^(5-α[i])/(5-α[i]) -
#         A[i]*breakpoints[i]^(5-α[i])/(5-α[i]) for i in 1:length(A) ) )
#     X3 = sum( (A[i]*breakpoints[i+1]^(4-α[i])/(4-α[i]) -
#         A[i]*breakpoints[i]^(4-α[i])/(4-α[i]) for i in 1:length(A) ) )
#     X2 = sum( (A[i]*breakpoints[i+1]^(3-α[i])/(3-α[i]) -
#         A[i]*breakpoints[i]^(3-α[i])/(3-α[i]) for i in 1:length(A) ) )
#     X = sum( (A[i]*breakpoints[i+1]^(2-α[i])/(2-α[i]) -
#         A[i]*breakpoints[i]^(2-α[i])/(2-α[i]) for i in 1:length(A) ) )
#     μ4 = X4 - 4*X*X3 + 6*X^2*X2 - 3*X^4
#     μ4 / var(d)^2
# end

#### Evaluation
function pdf(d::LogNormalBPL,x::Real)
    ((x < minimum(d)) || (x > maximum(d))) && (return zero(partype(d)))
    μ,σ,A,α,breakpoints = params(d)
    idx = findfirst(>=(x),breakpoints)
    ((idx==1) || (idx==2)) ? (return A[1] / x * exp(-(log(x)-μ)^2/(2*σ^2))) : (return A[idx-1] * x^-α[idx-2])
end
function logpdf(d::LogNormalBPL,x::Real)
    if ((x >= minimum(d)) && (x <= maximum(d)))
        μ,σ,A,α,breakpoints = params(d)
        idx = findfirst(>=(x),breakpoints)
        ((idx==1) || (idx==2)) ? (return log(A[1]) - log(x) - (log(x)-μ)^2/(2*σ^2)) : (return log(A[idx-1]) - α[idx-2]*log(x))
    else
        T = partype(d)
        return -T(Inf)
    end
end
function cdf(d::LogNormalBPL,x::Real)
    if x <= minimum(d)
        return zero(partype(d)) 
    elseif x >= maximum(d)
        return one(partype(d))
    end
    μ,σ,A,α,breakpoints = params(d)
    idx = findfirst(>=(x),breakpoints)
    result = lognormal_integral(A[1],μ,σ,breakpoints[1],min(x,breakpoints[2]))
    idx > 2 && (result += sum(pl_integral(A[i],α[i-1],breakpoints[i],min(x,breakpoints[i+1])) for i in 2:idx-1))
    return result
end
ccdf(d::LogNormalBPL,x::Real) = one(partype(d)) - cdf(d,x)
function quantile(d::LogNormalBPL{S},x::T) where {S,T<:Real}
    U = promote_type(S,T)
    x<=zero(T) && (return U(minimum(d)))
    x>=one(T) && (return U(maximum(d)))
    μ,σ,A,α,breakpoints = params(d)
    nbreaks = length(A)
    # this works but the tuple interpolation is slow and allocating, so switch to a vector
    # integrals = cumsum( (lognormal_integral(A[1],μ,σ,breakpoints[1],breakpoints[2]), (pl_integral(A[i],α[i-1],breakpoints[i],breakpoints[i+1]) for i in 2:nbreaks)...) ) # calculate the cumulative integral up to each breakpoint
    integrals = Array{U}(undef,nbreaks)
    integrals[1] = lognormal_integral(A[1],μ,σ,breakpoints[1],breakpoints[2])
    @inbounds for i in 2:nbreaks
        integrals[i] = pl_integral(A[i],α[i-1],breakpoints[i],breakpoints[i+1])
    end
    cumsum!(integrals,integrals)
    idx = findfirst(>=(x),integrals)  # find the first breakpoint where the cumulative integral
    if idx == 1
        # return exp(μ - sqrt(2) * σ * erfinv( (A[1] * π * σ * erf((μ-log(breakpoints[1]))/(sqrt(2)*σ)) - sqrt(2π)*x) / (A[1]*π*σ) ))
        return exp(μ - sqrt(U(2)) * σ * erfinv( (A[1] * π * σ * erf((μ-log(breakpoints[1]))/(sqrt(U(2))*σ)) - sqrt(U(2π))*x) / (A[1]*π*σ) ))
    else
        x -= integrals[idx-1]   # If this is not the first breakpoint, then subtract off the cumulative integral and solve 
        a = one(S) - α[idx-1] # using power law CDF inversion
        return (x*a/A[idx] + breakpoints[idx]^a)^inv(a)
    end
end
function quantile!(result::AbstractArray{U},d::LogNormalBPL{S},x::AbstractArray{T}) where {S,T<:Real,U<:Real}
    @assert size(result) == size(x)
    μ,σ,A,α,breakpoints = params(d)
    nbreaks = length(A)
    # this works but the tuple interpolation is slow and allocating, so switch to a vector
    # integrals = cumsum( (lognormal_integral(A[1],μ,σ,breakpoints[1],breakpoints[2]), (pl_integral(A[i],α[i-1],breakpoints[i],breakpoints[i+1]) for i in 2:nbreaks)...) ) # calculate the cumulative integral up to each breakpoint
    integrals = Array{eltype(result)}(undef,nbreaks)
    integrals[1] = lognormal_integral(A[1],μ,σ,breakpoints[1],breakpoints[2])
    @inbounds for i in 2:nbreaks
        integrals[i] = pl_integral(A[i],α[i-1],breakpoints[i],breakpoints[i+1])
    end
    cumsum!(integrals,integrals)
    @inbounds for i in eachindex(x)
        xi = x[i]
        xi<=zero(T) && (result[i]=minimum(d); continue)
        xi>=one(T) && (result[i]=maximum(d); continue)
        idx = findfirst(>=(xi),integrals)  # find the first breakpoint where the cumulative integral   # up to each breakpoint 
        if idx == 1
            # result[i] = exp(μ - sqrt(2) * σ * erfinv( (A[1] * π * σ * erf((μ-log(breakpoints[1]))/(sqrt(2)*σ)) - sqrt(2π)*xi) / (A[1]*π*σ) ))
            result[i] = exp(μ - sqrt(U(2)) * σ * erfinv( (A[1] * π * σ * erf((μ-log(breakpoints[1]))/(sqrt(U(2))*σ)) - sqrt(U(2π))*xi) / (A[1]*π*σ) ))
        else
            xi -= integrals[idx-1]   # If this is not the first breakpoint, then subtract off the cumulative integral and solve 
            a = one(S) - α[idx-1] # using power law CDF inversion
            result[i] = (xi*a/A[idx] + breakpoints[idx]^a)^inv(a)
        end
    end
    return result
end
quantile(d::LogNormalBPL{T},x::AbstractArray{S}) where {T,S<:Real} = quantile!(Array{promote_type(T,S)}(undef,size(x)),d,x)
cquantile(d::LogNormalBPL,x::Real) = quantile(d,1-x)

#### Random sampling
struct LogNormalBPLSampler{T} <: Sampleable{Univariate,Continuous}
    μ::T                   # lognormal mean
    σ::T                   # lognormal standard deviation
    A::Vector{T}           # normalization parameters
    α::Vector{T}           # power law indexes
    breakpoints::Vector{T} # bounds of each break
    integrals::Vector{T}   # cumulative integral up to each breakpoint
end
function LogNormalBPLSampler(d::LogNormalBPL)
    μ,σ,A,α,breakpoints = params(d)
    nbreaks = length(A)
    # this works but the tuple interpolation is slow and allocating, so switch to a vector
    # integrals = cumsum( (lognormal_integral(A[1],μ,σ,breakpoints[1],breakpoints[2]), (pl_integral(A[i],α[i-1],breakpoints[i],breakpoints[i+1]) for i in 2:nbreaks)...) ) # calculate the cumulative integral up to each breakpoint
    integrals = Array{partype(d)}(undef,nbreaks)
    integrals[1] = lognormal_integral(A[1],μ,σ,breakpoints[1],breakpoints[2])
    @inbounds for i in 2:nbreaks
        integrals[i] = pl_integral(A[i],α[i-1],breakpoints[i],breakpoints[i+1])
    end
    cumsum!(integrals,integrals)
    LogNormalBPLSampler(μ,σ,A,α,breakpoints,integrals)
end
function rand(rng::AbstractRNG, s::LogNormalBPLSampler{T}) where T
    x = rand(rng,T)
    μ,σ,A,α,breakpoints,integrals = s.μ,s.σ,s.A,s.α,s.breakpoints,s.integrals
    idx = findfirst(>=(x),integrals)  # find the first breakpoint where the cumulative integral   # up to each breakpoint 
    if idx == 1
        # return exp(μ - sqrt(2) * σ * erfinv( (A[1] * π * σ * erf((μ-log(breakpoints[1]))/(sqrt(2)*σ)) - sqrt(2π)*x) / (A[1]*π*σ) ))
        return exp(μ - sqrt(T(2)) * σ * erfinv( (A[1] * π * σ * erf((μ-log(breakpoints[1]))/(sqrt(T(2))*σ)) - sqrt(T(2π))*x) / (A[1]*π*σ) ))
    else
        x -= integrals[idx-1]   # If this is not the first breakpoint, then subtract off the cumulative integral and solve 
        a = one(T) - α[idx-1] # using power law CDF inversion
        return (x*a/A[idx] + breakpoints[idx]^a)^inv(a)
    end
end
sampler(d::LogNormalBPL) = LogNormalBPLSampler(d)
rand(rng::AbstractRNG, d::LogNormalBPL) = rand(rng, sampler(d))

#######################################################
# Specific types of LogNormalBPL
#######################################################

const chabrier2003_α = [2.3]
const chabrier2003_breakpoints = [0.0,1.0,Inf]
const chabrier2003_μ = log(0.079)#*log(10)
const chabrier2003_σ = 0.69*log(10)
"""
    Chabrier2003LogNormal(mmin::Real=0.08,mmax::Real=Inf)

Function to instantiate the [Chabrier 2003](https://ui.adsabs.harvard.edu/abs/2003PASP..115..763C/abstract) lognormal IMF, with a power-law extension for masses greater than one solar mass. This will return an instance of [`LogNormalBPL`](@ref). See also [`Chabrier2001LogNormal`](@ref) which has the same lognormal form, but without a high-mass power law extension.
"""
function Chabrier2003(mmin::T=0.08,mmax::T=Inf) where {T<:Real}
    @assert mmin>0
    mmin>1.0 && return PowerLawIMF(2.3,mmin,mmax) # if mmin>1, we are ONLY using the power law extension, so return power law IMF.
    mmax<1.0 && return truncated(LogNormal(chabrier2003_μ,chabrier2003_σ);lower=mmin,upper=mmax) # if mmax<1, we are ONLY using the lognormal component, so return lognormal IMF.
    idx1 = findfirst(>(mmin),chabrier2003_breakpoints)-1
    idx2 = findfirst(>=(mmax),chabrier2003_breakpoints)
    bp = convert(Vector{T},chabrier2003_breakpoints[idx1:idx2])
    bp[1] = mmin
    bp[end] = mmax
    LogNormalBPL(convert(T,chabrier2003_μ),convert(T,chabrier2003_σ),convert(Vector{T},chabrier2003_α[idx1:(idx2-2)]),bp)
end
Chabrier2003(mmin::Real,mmax::Real) = Chabrier2003(promote(mmin,mmax)...)
