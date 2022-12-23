"""
    PowerLaw(α::Real, mmin::Real, mmax::Real)

Descibes a single power-law IMF with probability distribution

```math
    \\frac{dn(m)}{dm} = A \\times m^{-\\alpha}
```

truncated such that the probability distribution is 0 below `mmin` and above `mmax`. `A` is a normalization constant such that the distribution integrates to 1 from the minimum valid stellar mass `mmin` to the maximum valid stellar mass `mmax`.

# Notes
The `quantile` and `rand` methods are currently fairly unoptimized. I did add a `randcoeff` field into the struct to hold the constant `mmin^(1-α)` which is necessary for `quantile` and `rand`, but it is not used by any other methods and is not returned by `params`; this should be considered a hidden, non-stable field. If the random sampling speed becomes problematic, the code for the Pareto distribution in Distributions.jl might be of some help. 
"""
struct PowerLaw{T} <: AbstractIMF
    A::T      # normalization parameters
    α::T      # power law index
    mmin::T   # minimum mass
    mmax::T   # minimum mass
    randcoeff::T # mmin^(1-α) for random sampling routine. 
    PowerLaw{T}(A::T, α::T, mmin::T, mmax::T, randcoeff::T) where {T} = new{T}(A,α,mmin,mmax,randcoeff)
end
function PowerLaw(α::Real, mmin::Real, mmax::Real)
    α, mmin, mmax = promote(α, mmin, mmax)
    T = typeof(α)
    if isfinite(mmax)
        mminα = mmin^α
        mmaxα = mmax^α
        A = (mmaxα * mminα * (α - 1)) / (mmaxα * mmin - mmax * mminα)
    else
        A = mmin^(α-1) * (α - 1)
    end
    randcoeff = mmin^(1-α)
    return PowerLaw{T}(A, α, mmin, mmax, randcoeff)
end

#### Conversions
Base.convert(::Type{PowerLaw{T}}, d::PowerLaw) where T = PowerLaw{T}(convert(T,d.A), convert(T,d.α), convert(T,d.mmin), convert(T,d.mmax), convert(T, d.randcoeff) )
Base.convert(::Type{PowerLaw{T}}, d::PowerLaw{T}) where {T<:Real} = d

#### Parameters
params(d::PowerLaw) = d.A,d.α,d.mmin,d.mmax
minimum(d::PowerLaw) = d.mmin
maximum(d::PowerLaw) = d.mmax
partype(d::PowerLaw{T}) where T = T

#### Statistics
function mean(d::PowerLaw)
    A,α,mmin,mmax = params(d)
    return A / (2-α) * (mmax^(2-α) - mmin^(2-α))
end
function median(d::PowerLaw)
    A,α,mmin,mmax = params(d)
    return (( 2A * mmin^(1-α) - α + 1) / 2A)^inv(1-α)
end
function var(d::PowerLaw)
    A,α,mmin,mmax = params(d)
    return A / (3-α) * (mmax^(3-α) - mmin^(3-α)) - mean(d)^2
end
function skewness(d::PowerLaw)
    A,α,mmin,mmax = params(d)
    σ² = var(d)
    μ = mean(d)
    return ( A / (4-α) * (mmax^(4-α) - mmin^(4-α)) - 3μ * σ² - μ^3 ) / sqrt( σ²^3 )
end
function kurtosis(d::PowerLaw)
    A,α,mmin,mmax = params(d)
    σ² = var(d)
    σ⁴ = σ²^2
    μ = mean(d)
    X2 = A / (3-α) * (mmax^(3-α) - mmin^(3-α))
    X3 = A / (4-α) * (mmax^(4-α) - mmin^(4-α))
    X4 = A / (5-α) * (mmax^(5-α) - mmin^(5-α))
    return ( X4 - 4*X3*μ + 6*X2*μ^2 - 3*μ^4 ) / σ⁴
end

#### Evaluation
function pdf(d::PowerLaw{S}, x::T) where {S, T<:Real}
    U = promote_type(S,T)
    if (x < minimum(d)) | (x > maximum(d))
        return zero(U) 
    else
        return d.A * x^-d.α
    end
end
function logpdf(d::PowerLaw{S}, x::T) where {S, T<:Real}
    U = promote_type(S,T)
    if (x < minimum(d)) | (x > maximum(d))
        return -U(Inf)
    else
        log(d.A) - d.α * log(x)
    end
end
function cdf(d::PowerLaw{S}, x::T) where {S, T<:Real}
    U = promote_type(S,T)
    if x <= minimum(d)
        return zero(U) 
    elseif x >= maximum(d)
        return one(U) 
    else
        A,α,mmin,mmax = params(d)
        return pl_integral(A, α, mmin, x)
    end
end
ccdf(d::PowerLaw, x::Real) = one(partype(d)) - cdf(d,x)
@inline function quantile(d::PowerLaw{S}, x::T) where {S, T<:Real}
    U = promote_type(S, T)
    if x <= zero(T)
        return U(minimum(d))
    elseif x >= one(T)
        return U(maximum(d))
    else
        A,α,mmin,mmax = params(d)
        # return (mmin^(1-α) + ( (1-α) * x / A))^inv(1-α)
        # return exp( inv(1-α) * log( (mmin^(1-α) + ( (1-α) * x / A)) ) ) # Log transform
        return exp( inv(1-α) * log( (d.randcoeff + ( (1-α) * x / A)) ) ) # Log transform with saved property
    end
end
function quantile!(result::AbstractArray, d::PowerLaw{S}, x::AbstractArray{T}) where {S, T<:Real}
    @assert size(result) == size(x)
    A,α,mmin,mmax = params(d)
    rcf = d.randcoeff
    @inbounds @simd for i in eachindex(x)
        xi = x[i]
        if xi <= zero(T)
            result[i] = U(minimum(d))
        elseif xi >= one(T)
            result[i] = U(maximum(d))
        else
            # result[i] = ( (A * mmin^(1-α) / (1-α) + xi) * (1-α) / A )^inv(1-α)
            result[i] = exp( inv(1-α) * log( (rcf + ( (1-α) * xi / A)) ) )
        end
    end
    return result
end
quantile(d::PowerLaw{T}, x::AbstractArray{S}) where {T, S<:Real} = quantile!(Array{promote_type(T,S)}(undef, size(x)), d, x)
cquantile(d::PowerLaw, x::Real) = quantile(d, 1-x)
# Distributions.jl defines fallbacks for other types of rand calls.
eltype(d::PowerLaw{T}) where T = T # Defines the type of random variates drawn. See docstring for eltype(::Type{Distributions.Sampleable}).
rand(rng::AbstractRNG, d::PowerLaw{T}) where T = quantile(d, rand(rng, T))


"""
    Salpeter1955(mmin::Real=0.4, mmax::Real=Inf)

The IMF model of [Salpeter 1955](https://ui.adsabs.harvard.edu/abs/1955ApJ...121..161S/abstract), a [`PowerLaw`](@ref) with `α=2.35`.
"""
Salpeter1955(mmin::T, mmax::T) where T<:Real = PowerLaw(T(2.35), mmin, mmax)
Salpeter1955(mmin::Real=0.4, mmax::Real=Inf) = Salpeter1955(promote(mmin,mmax)...)
