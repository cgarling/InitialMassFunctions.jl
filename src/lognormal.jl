"""
    LogNormalIMF(μ::Real, σ::Real, mmin::Real, mmax::Real)

Describes a lognormal IMF with probability distribution

```math
    \\frac{dn(m)}{dm} = \\frac{A}{x} \\, \\exp \\left[ \\frac{ -\\left( \\log(x) - \\mu \\right)^2}{2\\sigma^2} \\right]
```

truncated such that the probability distribution is 0 below `mmin` and above `mmax`. `A` is a normalization constant such that the distribution integrates to 1 from the minimum valid stellar mass `mmin` to the maximum valid stellar mass `mmax`. This is simply `Distributions.truncated(Distributions.LogNormal(μ,σ);lower=mmin,upper=mmax)`. See the documentation for [`LogNormal`](https://juliastats.org/Distributions.jl/stable/univariate/#Distributions.LogNormal) and [`truncated`](https://juliastats.org/Distributions.jl/latest/truncate/#Distributions.truncated).

# Arguments
 - `μ`; see [Distributions.LogNormal](https://juliastats.org/Distributions.jl/stable/univariate/#Distributions.LogNormal)
 - `σ`; see [Distributions.LogNormal](https://juliastats.org/Distributions.jl/stable/univariate/#Distributions.LogNormal)
"""
LogNormalIMF(μ::Real, σ::Real, mmin::Real, mmax::Real) = truncated(LogNormal(μ,σ);lower=mmin,upper=mmax)
function mean(d::Truncated{LogNormal{T}, Continuous, T}) where T
    mmin, mmax = extrema(d)
    μ, σ = params( d.untruncated )
    # return (α * θ^α / (1-α) / d.ucdf) * (mmax^(1-α) - mmin^(1-α))
    return -exp(μ + σ^2/2) / 2 / (d.ucdf - d.lcdf) *
        ( erf( (μ + σ^2 - log(mmax)) / (sqrt(2)*σ) ) - erf( (μ + σ^2 - log(mmin)) / (sqrt(2)*σ) ) )
end
"""
    Chabrier2001LogNormal(mmin::Real=0.08, mmax::Real=Inf)

Function to instantiate the [Chabrier 2001](https://ui.adsabs.harvard.edu/abs/2001ApJ...554.1274C/abstract) lognormal IMF for single stars. Returns an instance of `Distributions.Truncated(Distributions.LogNormal)`. See also [`Chabrier2003`](@ref) which has the same lognormal form for masses below one solar mass, but a power law extension at higher masses. 
"""
Chabrier2001LogNormal(mmin::Real=0.08, mmax::Real=Inf) = LogNormalIMF(log(0.1), 0.627*log(10), mmin, mmax)

"""
    lognormal_integral(μ, σ, b1, b2)

Definite integral of the lognormal probability distribution from `b1` to `b2`.
```math
\\int_{b1}^{b2} \\, \\frac{A}{x} \\, \\exp \\left[ \\frac{ -\\left( \\log(x) - \\mu \\right)^2}{2\\sigma^2} \\right] \\, dx
```
"""
lognormal_integral(A::T,μ::T,σ::T,b1::T,b2::T) where {T<:Number} =
    A * sqrt(T(π)/2) * σ * (erf( (μ-log(b1))/(sqrt(T(2))*σ)) - erf( (μ-log(b2))/(sqrt(T(2))*σ)))
lognormal_integral(A::Number,μ::Number,σ::Number,b1::Number,b2::Number) = lognormal_integral(promote(A,μ,σ,b1,b2)...)
# lognormal_integral(A,μ,σ,b1,b2) = A * sqrt(π/2) * σ * (erf( (μ-log(b1))/(sqrt(2)*σ)) - erf( (μ-log(b2))/(sqrt(2)*σ)))

"""
    LogNormalBPL(μ::Real,σ::Real,α::AbstractVector{<:Real},breakpoints::AbstractVector{<:Real})
    LogNormalBPL(μ::Real,σ::Real,α::Tuple,breakpoints::Tuple)

A LogNormal distribution at low masses, with a broken power law extension at high masses. This uses the natural log base like [Distributions.LogNormal](https://juliastats.org/Distributions.jl/stable/univariate/#Distributions.LogNormal); if you have σ and μ in base 10, then multiply them both by `log(10)`. Must have `length(α) == length(breakpoints)-2`. The probability distribution for this IMF model is

```math
    \\frac{dn(m)}{dm} = \\frac{A}{x} \\, \\exp \\left[ \\frac{ -\\left( \\log(x) - \\mu \\right)^2}{2\\sigma^2} \\right]
```
for `m < breakpoints[2]`, with a broken power law extension above this mass. See [`BrokenPowerLaw`](@ref) for interface details; the `α` and `breakpoints` are the same here as there.

# Arguments
 - `μ`; see [Distributions.LogNormal](https://juliastats.org/Distributions.jl/stable/univariate/#Distributions.LogNormal)
 - `σ`; see [Distributions.LogNormal](https://juliastats.org/Distributions.jl/stable/univariate/#Distributions.LogNormal)
 - `α`; list of power law indices with `length(α) == length(breakpoints)-2`.
 - `breakpoints`; list of masses that signal breaks in the IMF. MUST BE SORTED and bracketed with `breakpoints[1]` being the minimum valid mass and `breakpoints[end]` being the maximum valid mass.

# Examples
If you want a `LogNormalBPL` with a characteristic mass of 0.5 solar masses, log10 standard deviation of 0.6, and a single power law extension with slope `α=2.35` with a break at 1 solar mass, you would do `LogNormalBPL(log(0.5),0.6*log(10),[2.35],[0.08,1.0,Inf]` where we set the minimum mass to `0.08` and maximum mass to `Inf`. If, instead, you know that `log10(m)=x`, where `m` is the characteristic mass of the `LogNormal` component, you would do `LogNormalBPL(x*log(10),0.6*log(10),[2.35],[0.08,1.0,Inf]`.

# Notes
There is some setup necessary for `quantile` and other derived methods, so it is more efficient to call these methods directly with an array via the call signature `quantile(d::LogNormalBPL{T}, x::AbstractArray{S})` rather than broadcasting over `x`. This behavior is now deprecated for `quantile(d::Distributions.UnivariateDistribution, X::AbstractArray)` in Distributions.jl. 

# Methods
 - `Base.convert(::Type{LogNormalBPL{T}}, d::LogNormalBPL)`
 - `minimum(d::LogNormalBPL)`
 - `maximum(d::LogNormalBPL)`
 - `partype(d::LogNormalBPL)`
 - `eltype(d::LogNormalBPL)`
 - `mean(d::LogNormalBPL)`
 - `median(d::LogNormalBPL)`
 - `pdf(d::LogNormalBPL,x::Real)`
 - `logpdf(d::LogNormalBPL,x::Real)`
 - `cdf(d::LogNormalBPL,x::Real)`
 - `ccdf(d::LogNormalBPL,x::Real)`
 - `quantile(d::LogNormalBPL{S},x::T) where {S,T<:Real}`
 - `quantile!(result::AbstractArray,d::LogNormalBPL{S},x::AbstractArray{T}) where {S,T<:Real}`
 - `quantile(d::LogNormalBPL{T},x::AbstractArray{S})`
 - `cquantile(d::LogNormalBPL{S},x::T) where {S,T<:Real}`
 - `rand(rng::AbstractRNG, d::LogNormalBPL,s...)` 
 - Other methods from `Distributions.jl` should also work because `LogNormalBPL <: AbstractIMF <: Distributions.ContinuousUnivariateDistribution`. For example, `rand!(rng::AbstractRNG, d::LogNormalBPL, x::AbstractArray)`.
"""
struct LogNormalBPL{T} <: AbstractIMF
    μ::T
    σ::T
    A::Vector{T}      # normalization parameters
    α::Vector{T}      # power law indexes
    breakpoints::Vector{T} # bounds of each break
    LogNormalBPL{T}(μ::T, σ::T, A::Vector{T}, α::Vector{T}, breakpoints::Vector{T}) where {T} =
        new{T}(μ, σ, A, α, breakpoints)
end
function LogNormalBPL(μ::T, σ::T, α::Vector{T}, breakpoints::Vector{T}) where T <: Real
    @assert length(breakpoints) == length(α)+2
    @assert breakpoints[1] > 0
    nbreaks = length(α) + 1
    A = Vector{T}(undef,nbreaks)
    A[1] = one(T)
    # solve for the prefactor for the first power law after the lognormal component
    A[2] = breakpoints[2]^(α[1]-1) * exp( -(log(breakpoints[2])-μ)^2/(2*σ^2))
    # if there is more than one power law component
    if nbreaks > 2
        for i in 3:nbreaks
            A[i] = A[i-1]*breakpoints[i]^-α[i-2] / breakpoints[i]^-α[i-1]
        end
    end
    # now A contains prefactors for each distribution component that makes them continuous
    # with the first lognormal component having a prefactor of 1. Now we need to normalize
    # the integral from minimum(breaks) to maximum(breaks) to equal 1 by dividing
    # the entire A array by a common factor.
    total_integral = zero(T)
    total_integral += lognormal_integral(A[1], μ, σ, breakpoints[1], breakpoints[2])
    for i in 2:nbreaks
        total_integral += pl_integral(A[i], α[i-1], breakpoints[i], breakpoints[i+1])
    end
    A ./= total_integral
    return LogNormalBPL{T}(μ, σ, A, α, breakpoints)
end
LogNormalBPL(μ::Real, σ::Real, α::Tuple, breakpoints::Tuple) =
    LogNormalBPL(μ,σ,collect(promote(α...)),collect(promote(breakpoints...)))
LogNormalBPL(μ::T,σ::T,α::AbstractVector{T},breakpoints::AbstractVector{T}) where T<:Real =
    LogNormalBPL(μ,σ,convert(Vector{T},α),convert(Vector{T},breakpoints))
function LogNormalBPL(μ::A, σ::B, α::AbstractVector{C}, breakpoints::AbstractVector{D}) where {A<:Real, B<:Real, C<:Real, D<:Real}
    X = promote_type(A, B, C, D)
    LogNormalBPL(convert(X,μ),convert(X,σ),convert(Vector{X},α), convert(Vector{X},breakpoints))
end

#### Conversions
Base.convert(::Type{LogNormalBPL{T}}, d::LogNormalBPL) where T <: Real =
    LogNormalBPL{T}(convert(T,d.μ), convert(T,d.σ), convert(Vector{T},d.A), convert(Vector{T},d.α), convert(Vector{T},d.breakpoints))
Base.convert(::Type{LogNormalBPL{T}}, d::LogNormalBPL{T}) where T <: Real = d

#### Parameters
params(d::LogNormalBPL) = d.μ, d.σ, d.A, d.α, d.breakpoints
minimum(d::LogNormalBPL) = minimum(d.breakpoints)
maximum(d::LogNormalBPL) = maximum(d.breakpoints)
partype(d::LogNormalBPL{T}) where T = T
eltype(d::LogNormalBPL{T}) where T = T

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
median(d::LogNormalBPL{T}) where T = quantile(d, T(0.5)) # this is temporary
# mode(d::BrokenPowerLaw) = d.breakpoints[argmin(d.α)] # this is not always correct

#### Evaluation
function pdf(d::LogNormalBPL, x::Real)
    ((x < minimum(d)) || (x > maximum(d))) && (return zero(partype(d)))
    μ, σ, A, α, breakpoints = params(d)
    idx = findfirst(>=(x), breakpoints)
    ((idx==1) || (idx==2)) ? (return A[1] / x * exp(-(log(x)-μ)^2/(2*σ^2))) : (return A[idx-1] * x^-α[idx-2])
end
function logpdf(d::LogNormalBPL, x::Real)
    if ((x >= minimum(d)) && (x <= maximum(d)))
        μ, σ, A, α, breakpoints = params(d)
        idx = findfirst(>=(x), breakpoints)
        ((idx==1) || (idx==2)) ? (return log(A[1]) - log(x) - (log(x)-μ)^2/(2*σ^2)) : (return log(A[idx-1]) - α[idx-2]*log(x))
    else
        T = partype(d)
        return -T(Inf)
    end
end
function cdf(d::LogNormalBPL, x::Real)
    if x <= minimum(d)
        return zero(partype(d)) 
    elseif x >= maximum(d)
        return one(partype(d))
    end
    μ, σ, A, α, breakpoints = params(d)
    idx = findfirst(>=(x), breakpoints)
    result = lognormal_integral(A[1], μ, σ, breakpoints[1], min(x,breakpoints[2]))
    idx > 2 && (result += sum(pl_integral(A[i],α[i-1],breakpoints[i],min(x,breakpoints[i+1])) for i in 2:idx-1))
    return result
end
ccdf(d::LogNormalBPL, x::Real) = one(partype(d)) - cdf(d,x)
function quantile(d::LogNormalBPL{S}, x::T) where {S, T <: Real}
    U = promote_type(S, T)
    x <= zero(T) && (return U(minimum(d)))
    x >= one(T) && (return U(maximum(d)))
    μ, σ, A, α, breakpoints = params(d)
    nbreaks = length(A)
    # this works but the tuple interpolation is slow and allocating, so switch to a vector
    # integrals = cumsum( (lognormal_integral(A[1],μ,σ,breakpoints[1],breakpoints[2]), (pl_integral(A[i],α[i-1],breakpoints[i],breakpoints[i+1]) for i in 2:nbreaks)...) ) # calculate the cumulative integral up to each breakpoint
    integrals = Array{U}(undef,nbreaks)
    integrals[1] = lognormal_integral(A[1],μ,σ,breakpoints[1],breakpoints[2])
    for i in 2:nbreaks
        integrals[i] = pl_integral(A[i],α[i-1],breakpoints[i],breakpoints[i+1])
    end
    cumsum!(integrals, integrals)
    idx = findfirst(>=(x), integrals)  # find the first breakpoint where the cumulative integral
    if idx == 1
        # return exp(μ - sqrt(2) * σ * erfinv( (A[1] * π * σ * erf((μ-log(breakpoints[1]))/(sqrt(2)*σ)) - sqrt(2π)*x) / (A[1]*π*σ) ))
        return exp(μ - sqrt(U(2)) * σ * erfinv( (A[1] * π * σ * erf((μ-log(breakpoints[1]))/(sqrt(U(2))*σ)) - sqrt(U(2π))*x) / (A[1]*π*σ) ))
    else
        x -= integrals[idx-1]   # If this is not the first breakpoint, then subtract off the cumulative integral and solve 
        a = one(S) - α[idx-1] # using power law CDF inversion
        return (x*a/A[idx] + breakpoints[idx]^a)^inv(a)
    end
end
function quantile!(result::AbstractArray{U}, d::LogNormalBPL{S}, x::AbstractArray{T}) where {S, T<:Real, U<:Real}
    @assert axes(result) == axes(x)
    μ, σ, A, α, breakpoints = params(d)
    nbreaks = length(A)
    # this works but the tuple interpolation is slow and allocating, so switch to a vector
    # integrals = cumsum( (lognormal_integral(A[1],μ,σ,breakpoints[1],breakpoints[2]), (pl_integral(A[i],α[i-1],breakpoints[i],breakpoints[i+1]) for i in 2:nbreaks)...) ) # calculate the cumulative integral up to each breakpoint
    integrals = Array{eltype(result)}(undef,nbreaks)
    integrals[1] = lognormal_integral(A[1],μ,σ,breakpoints[1],breakpoints[2])
    @inbounds for i in 2:nbreaks
        integrals[i] = pl_integral(A[i],α[i-1],breakpoints[i],breakpoints[i+1])
    end
    cumsum!(integrals, integrals)
    @inbounds for i in eachindex(x)
        xi = x[i]
        xi <= zero(T) && (result[i]=minimum(d); continue)
        xi >= one(T) && (result[i]=maximum(d); continue)
        idx = findfirst(>=(xi), integrals)  # find the first breakpoint where the cumulative integral   # up to each breakpoint 
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
quantile(d::LogNormalBPL{T}, x::AbstractArray{S}) where {T, S <: Real} =
    quantile!(Array{promote_type(T,S)}(undef,size(x)), d, x)
cquantile(d::LogNormalBPL,x::Real) = quantile(d, 1-x)

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
    μ, σ, A, α, breakpoints = params(d)
    nbreaks = length(A)
    # this works but the tuple interpolation is slow and allocating, so switch to a vector
    # integrals = cumsum( (lognormal_integral(A[1],μ,σ,breakpoints[1],breakpoints[2]), (pl_integral(A[i],α[i-1],breakpoints[i],breakpoints[i+1]) for i in 2:nbreaks)...) ) # calculate the cumulative integral up to each breakpoint
    integrals = Array{partype(d)}(undef,nbreaks)
    integrals[1] = lognormal_integral(A[1],μ,σ,breakpoints[1],breakpoints[2])
    @inbounds for i in 2:nbreaks
        integrals[i] = pl_integral(A[i],α[i-1],breakpoints[i],breakpoints[i+1])
    end
    cumsum!(integrals, integrals)
    LogNormalBPLSampler(μ, σ, A, α, breakpoints, integrals)
end
function rand(rng::AbstractRNG, s::LogNormalBPLSampler{T}) where T
    x = rand(rng, T)
    μ, σ, A, α, breakpoints, integrals = s.μ, s.σ, s.A, s.α, s.breakpoints, s.integrals
    idx = findfirst(>=(x), integrals)  # find the first breakpoint where the cumulative integral   # up to each breakpoint 
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
    Chabrier2003(mmin::Real=0.08, mmax::Real=Inf)

Function to instantiate the [Chabrier 2003](https://ui.adsabs.harvard.edu/abs/2003PASP..115..763C/abstract) IMF for single stars. This is a lognormal IMF with a power-law extension for masses greater than one solar mass. This IMF is valid for single stars and takes parameters from the "Disk and Young Clusters" column of Table 2 in the above paper. This will return an instance of [`LogNormalBPL`](@ref). See also [`Chabrier2003System`](@ref) which implements the IMF for general stellar systems with multiplicity from this same paper, and [`Chabrier2001LogNormal`](@ref) which has the same lognormal form as this model but without a high-mass power law extension.
"""
function Chabrier2003(mmin::T=0.08, mmax::T=Inf) where T <: Real
    @assert mmin > 0
    mmin > one(T) && return PowerLaw(2.3,mmin,mmax) # if mmin>1, we are ONLY using the power law extension, so return power law IMF.
    mmax < one(T) && return truncated(LogNormal(chabrier2003_μ,chabrier2003_σ);lower=mmin,upper=mmax) # if mmax<1, we are ONLY using the lognormal component, so return lognormal IMF.
    idx1 = findfirst(>(mmin), chabrier2003_breakpoints)-1
    idx2 = findfirst(>=(mmax), chabrier2003_breakpoints)
    bp = convert(Vector{T}, chabrier2003_breakpoints[idx1:idx2])
    bp[1] = mmin
    bp[end] = mmax
    LogNormalBPL(convert(T,chabrier2003_μ), convert(T,chabrier2003_σ), convert(Vector{T},chabrier2003_α[idx1:(idx2-2)]), bp)
end
Chabrier2003(mmin::Real, mmax::Real) = Chabrier2003(promote(mmin,mmax)...)
#######################################################
const chabrier2003_system_μ = log(0.22)
const chabrier2003_system_σ = 0.57*log(10)
"""
    Chabrier2003System(mmin::Real=0.08, mmax::Real=Inf)

Function to instantiate the [Chabrier 2003](https://ui.adsabs.harvard.edu/abs/2003PASP..115..763C/abstract) system IMF. This is a lognormal IMF with a power-law extension for masses greater than one solar mass. This IMF is valid for general star systems with stellar multiplicity (e.g., binaries) and differs from the typical single-star models. Parameters for this distribution are taken from Equation 18 in the above paper. This will return an instance of [`LogNormalBPL`](@ref). See also [`Chabrier2003`](@ref) for the single star IMF.
"""
function Chabrier2003System(mmin::T=0.08, mmax::T=Inf) where T <: Real
    @assert mmin > 0
    mmin > one(T) && return PowerLaw(2.3,mmin,mmax) # if mmin>1, we are ONLY using the power law extension, so return power law IMF.
    mmax < one(T) && return truncated(LogNormal(chabrier2003_system_μ,chabrier2003_system_σ);lower=mmin,upper=mmax) # if mmax<1, we are ONLY using the lognormal component, so return lognormal IMF.
    idx1 = findfirst(>(mmin), chabrier2003_breakpoints)-1
    idx2 = findfirst(>=(mmax), chabrier2003_breakpoints)
    bp = convert(Vector{T}, chabrier2003_breakpoints[idx1:idx2])
    bp[1] = mmin
    bp[end] = mmax
    LogNormalBPL(convert(T,chabrier2003_system_μ), convert(T,chabrier2003_system_σ), convert(Vector{T},chabrier2003_α[idx1:(idx2-2)]), bp)
end
Chabrier2003System(mmin::Real, mmax::Real) = Chabrier2003System(promote(mmin,mmax)...)
