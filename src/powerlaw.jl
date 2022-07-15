# this nonsense of trying to use the Distributions objects seems ill-fated.
# Switch back to manual implementations
struct PowerLawIMF{T<:Real,S<:ContinuousUnivariateDistribution} <: AbstractIMF
    ξ0::T
    mmin::T
    mmax::T
    pdf::S
    PowerLawIMF{T,S}(ξ0::T, mmin::T, mmax::T, pdf::S) where {T,S} = new{T,S}(ξ0, mmin, mmax, pdf)
end
function PowerLawIMF(ξ0::T,mmin::T,mmax::T,pdf::S) where {T<:Real,S<:ContinuousUnivariateDistribution}
    @assert mmin<mmax
    @assert mmin>0
    PowerLawIMF{T,S}(ξ0,mmin,mmax,pdf)
end
PowerLawIMF(ξ0::Real,mmin::Real,mmax::Real,pdf::S) where {S<:ContinuousUnivariateDistribution} = PowerLawIMF(promote(ξ0,mmin,mmax)...,pdf)

dndm(law::PowerLawIMF,x::Real) = pdf(distribution(law),x)
dndlogm(law::PowerLawIMF,x::Real) = pdf(distribution(law),x)

Salpeter1955(mmin::Real=0.4,mmax::Real=10.0) = PowerLawIMF(0.03,mmin,mmax,truncated(Pareto(1.35,mmin);upper=mmax))
