PowerLawIMF(α::Real,mmin::Real,mmax::Real) = truncated(Pareto(α-1,mmin);upper=mmax)
Salpeter1955(mmin::Real=0.4,mmax::Real=10.0) = PowerLawIMF(2.35,mmin,mmax)

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
