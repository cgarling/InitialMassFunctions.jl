# this nonsense of trying to use the Distributions objects seems ill-fated.
# Switch back to manual implementations
struct PowerLawIMF{T<:Real} <: AbstractIMF
    α::T
    mmin::T
    mmax::T
    ξ0::T
    PowerLawIMF{T}(α::T, mmin::T, mmax::T, ξ0::T) where {T} = new{T}(α, mmin, mmax, ξ0)
end
function PowerLawIMF(α::T,mmin::T,mmax::T,ξ0::T) where {T<:Real}
    @assert mmin<mmax
    @assert mmin>0
    @assert α>0
    @assert ξ0>0
    PowerLawIMF{T}(α,mmin,mmax,ξ0)
end
PowerLawIMF(α::Real,mmin::Real,mmax::Real,ξ0::Real) = PowerLawIMF(promote(α,mmin,mmax,ξ0)...)
function PowerLawIMF(α::Real,mmin::Real,mmax::Real)
    α1 = α + 1
    ξ0 = α * mmax^α1 * mmin^α1 / (mmax^α1 * mmin - mmax * mmin^α1)
    return PowerLawIMF(promote(α,mmin,mmax,ξ0)...)
end
normalization(x::PowerLawIMF) = x.ξ0
prefactor(x::PowerLawIMF) = x.A
slope(x::PowerLawIMF) = x.α + 1
logslope(x::PowerLawIMF) = x.α
function dndm(law::PowerLawIMF,m::Real)
    m <= zero(m) && return zero(m)
    α = slope(law) # linear slope for dn/dm, not dn/dlogm
    return  normalization(law) / m^α
end
dndlogm(law::PowerLawIMF,m::T) where {T<:Real} = dndm(law,m) * m * T(log(10))
# m > zero(m) ? (return normalization(law) / m^(slope(law)-1) * T(log(10))) : (return zero(m))

function pdf(law::PowerLawIMF,m::Real)
    m <= zero(m) && return zero(m)
    α = slope(law) # linear slope for dn/dm, not dn/dlogm
    mmin,mmax=limits(law) # need to recalculate the normalization
    ξ0 = (α-1) * mmax^α * mmin^α / (mmax^α * mmin - mmax * mmin^α)
    return ξ0 / m^α
end
function logpdf(law::PowerLawIMF,m::Real)
    @assert m>=zero(m)
    α = slope(law) # linear slope for dn/dm, not dn/dlogm
    mmin,mmax=limits(law) # need to recalculate the normalization
    ξ0 = (α-1) * mmax^α * mmin^α / (mmax^α * mmin - mmax * mmin^α)
    log(ξ0) - α * log(m)
end
function cdf(law::PowerLawIMF,m::Real)
    m <= zero(m) && return zero(m)
    α = slope(law) # linear slope for dn/dm, not dn/dlogm
    mmin,mmax=limits(law) # need to recalculate the normalization
    ξ0 = (α-1) * mmax^α * mmin^α / (mmax^α * mmin - mmax * mmin^α)
    α2 = 1-α
    return (ξ0 * m^α2 - ξ0 * mmin^α2) / α2
end
function median(law::PowerLawIMF)
    α = slope(law) # linear slope for dn/dm, not dn/dlogm
    mmin,mmax=limits(law) # need to recalculate the normalization
    ξ0 = (α-1) * mmax^α * mmin^α / (mmax^α * mmin - mmax * mmin^α)
    return ((0.5/ξ0) * (1 - α) + mmin^(1-α))^(1/(1-α))
end
function mode(law::PowerLawIMF)
    mmin,mmax = limits(law)
    mmin
end

# Particular types of PowerLawIMFs
Salpeter1955(mmin::Real=0.4,mmax::Real=10.0) = PowerLawIMF(1.35,mmin,mmax,0.03)
