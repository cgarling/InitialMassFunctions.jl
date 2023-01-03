var documenterSearchIndex = {"docs":
[{"location":"constructors/#Convenience-Constructors-for-Published-IMFs","page":"Convenience Constructors for Published IMFs","title":"Convenience Constructors for Published IMFs","text":"","category":"section"},{"location":"constructors/","page":"Convenience Constructors for Published IMFs","title":"Convenience Constructors for Published IMFs","text":"We provide convenience constructors for published IMFs that can be called without arguments, or with positional arguments to set different minimum and maximum stellar masses. The provided constructors are","category":"page"},{"location":"constructors/","page":"Convenience Constructors for Published IMFs","title":"Convenience Constructors for Published IMFs","text":"Salpeter1955\nChabrier2001BPL\nKroupa2001\nChabrier2001LogNormal\nChabrier2003","category":"page"},{"location":"constructors/#InitialMassFunctions.Salpeter1955","page":"Convenience Constructors for Published IMFs","title":"InitialMassFunctions.Salpeter1955","text":"Salpeter1955(mmin::Real=0.4, mmax::Real=Inf)\n\nThe IMF model of Salpeter 1955, a PowerLawIMF with α=2.35.\n\n\n\n\n\n","category":"function"},{"location":"constructors/#InitialMassFunctions.Chabrier2001BPL","page":"Convenience Constructors for Published IMFs","title":"InitialMassFunctions.Chabrier2001BPL","text":"Chabrier2001BPL(mmin::T=0.08, mmax::T=Inf)\n\nFunction to instantiate a BrokenPowerLaw IMF with the parameters from the first column of Table 1 in Chabrier 2001.\n\n\n\n\n\n","category":"function"},{"location":"constructors/#InitialMassFunctions.Kroupa2001","page":"Convenience Constructors for Published IMFs","title":"InitialMassFunctions.Kroupa2001","text":"Kroupa2001(mmin::Real=0.08, mmax::Real=Inf)\n\nFunction to instantiate a BrokenPowerLaw IMF with the parameters from Equation 2 of Kroupa 2001. This is equivalent to the relation given in Kroupa 2002.\n\n\n\n\n\n","category":"function"},{"location":"constructors/#InitialMassFunctions.Chabrier2001LogNormal","page":"Convenience Constructors for Published IMFs","title":"InitialMassFunctions.Chabrier2001LogNormal","text":"Chabrier2001LogNormal(mmin::Real=0.08, mmax::Real=Inf)\n\nFunction to instantiate the Chabrier 2001 lognormal IMF. Returns an instance of Distributions.Truncated(Distributions.LogNormal). See also Chabrier2003 which has the same lognormal form for masses below one solar mass, but a power law extension at higher masses. \n\n\n\n\n\n","category":"function"},{"location":"constructors/#InitialMassFunctions.Chabrier2003","page":"Convenience Constructors for Published IMFs","title":"InitialMassFunctions.Chabrier2003","text":"Chabrier2003LogNormal(mmin::Real=0.08, mmax::Real=Inf)\n\nFunction to instantiate the Chabrier 2003 lognormal IMF, with a power-law extension for masses greater than one solar mass. This will return an instance of LogNormalBPL. See also Chabrier2001LogNormal which has the same lognormal form, but without a high-mass power law extension.\n\n\n\n\n\n","category":"function"},{"location":"constructors/","page":"Convenience Constructors for Published IMFs","title":"Convenience Constructors for Published IMFs","text":"We also provide a constructor for a single-power-law IMF,","category":"page"},{"location":"constructors/","page":"Convenience Constructors for Published IMFs","title":"Convenience Constructors for Published IMFs","text":"PowerLawIMF","category":"page"},{"location":"constructors/#InitialMassFunctions.PowerLawIMF","page":"Convenience Constructors for Published IMFs","title":"InitialMassFunctions.PowerLawIMF","text":"PowerLawIMF(α::Real, mmin::Real, mmax::Real)\n\nDescibes a single power-law IMF with probability distribution\n\n    fracdn(m)dm = A times m^-alpha\n\ntruncated such that the probability distribution is 0 below mmin and above mmax. A is a normalization constant such that the distribution integrates to 1 from the minimum valid stellar mass mmin to the maximum valid stellar mass mmax. This is simply Distributions.truncated(Distributions.Pareto(α-1,mmin);upper=mmax). See the documentation for Pareto and truncated.\n\n\n\n\n\n","category":"function"},{"location":"constructors/","page":"Convenience Constructors for Published IMFs","title":"Convenience Constructors for Published IMFs","text":"which internally creates a truncated Pareto distribution. Similarly, we provide a constructor for a single-component LogNormal IMF,","category":"page"},{"location":"constructors/","page":"Convenience Constructors for Published IMFs","title":"Convenience Constructors for Published IMFs","text":"LogNormalIMF","category":"page"},{"location":"constructors/#InitialMassFunctions.LogNormalIMF","page":"Convenience Constructors for Published IMFs","title":"InitialMassFunctions.LogNormalIMF","text":"LogNormalIMF(μ::Real, σ::Real, mmin::Real, mmax::Real)\n\nDescribes a lognormal IMF with probability distribution\n\n    fracdn(m)dm = fracAx  exp left frac -left( log(x) - mu right)^22sigma^2 right\n\ntruncated such that the probability distribution is 0 below mmin and above mmax. A is a normalization constant such that the distribution integrates to 1 from the minimum valid stellar mass mmin to the maximum valid stellar mass mmax. This is simply Distributions.truncated(Distributions.LogNormal(μ,σ);lower=mmin,upper=mmax). See the documentation for LogNormal and truncated.\n\nArguments\n\nμ; see Distributions.LogNormal\nσ; see Distributions.LogNormal\n\n\n\n\n\n","category":"function"},{"location":"utilities/#Utilities","page":"Utilities","title":"Utilities","text":"","category":"section"},{"location":"utilities/","page":"Utilities","title":"Utilities","text":"InitialMassFunctions.pl_integral\nInitialMassFunctions.lognormal_integral","category":"page"},{"location":"utilities/#InitialMassFunctions.pl_integral","page":"Utilities","title":"InitialMassFunctions.pl_integral","text":"pl_integral(A,α,b1,b2)\n\nDefinite integral of power law A*x^-α from b1 (lower) to b2 (upper).\n\nint_b1^b2  A times x^-alpha  dx = fracA1-alpha times left( b2^1-alpha - b1^1-alpha right)\n\n\n\n\n\n","category":"function"},{"location":"utilities/#InitialMassFunctions.lognormal_integral","page":"Utilities","title":"InitialMassFunctions.lognormal_integral","text":"lognormal_integral(μ, σ, b1, b2)\n\nDefinite integral of the lognormal probability distribution from b1 to b2.\n\nint_b1^b2  fracAx  exp left frac -left( log(x) - mu right)^22sigma^2 right  dx\n\n\n\n\n\n","category":"function"},{"location":"#Index","page":"Index","title":"Index","text":"","category":"section"},{"location":"","page":"Index","title":"Index","text":"","category":"page"},{"location":"types/#Defined-Types","page":"Defined Types","title":"Defined Types","text":"","category":"section"},{"location":"types/","page":"Defined Types","title":"Defined Types","text":"We provide two new types that are subtypes of AbstractIMF, which itself is a subtype of Distributions.ContinuousUnivariateDistribution and generally follow the API provided by Distributions.jl. These are","category":"page"},{"location":"types/","page":"Defined Types","title":"Defined Types","text":"InitialMassFunctions.AbstractIMF\nBrokenPowerLaw\nLogNormalBPL","category":"page"},{"location":"types/#InitialMassFunctions.AbstractIMF","page":"Defined Types","title":"InitialMassFunctions.AbstractIMF","text":"Abstract type for IMFs; a subtype of Distributions.ContinuousUnivariateDistribution, as all IMF models can be described as continuous, univariate PDFs. \n\n\n\n\n\n","category":"type"},{"location":"types/#InitialMassFunctions.BrokenPowerLaw","page":"Defined Types","title":"InitialMassFunctions.BrokenPowerLaw","text":"BrokenPowerLaw(α::AbstractVector{T}, breakpoints::AbstractVector{S}) where {T<:Real,S<:Real}\nBrokenPowerLaw(α::Tuple, breakpoints::Tuple)\nBrokenPowerLaw{T}(A::Vector{T}, α::Vector{T}, breakpoints::Vector{T}) where {T}\n\nAn AbstractIMF <: Distributions.ContinuousUnivariateDistribution that describes a broken power-law IMF with probability distribution\n\n    fracdn(m)dm = A times m^-alpha\n\nthat is defined piecewise with different normalizations A and power law slopes α in different mass ranges. The normalization constants A will be calculated automatically after you provide the power law slopes and break points.\n\nArguments\n\nα; the power-law slopes of the different segments of the broken power law.\nbreakpoints; the masses at which the power law slopes change. If length(α)=n, then length(breakpoints)=n+1.\n\nExamples\n\nBrokenPowerLaw([1.35,2.35],[0.08,1.0,Inf]) will instantiate a broken power law defined from a minimum mass of 0.08 to a maximum mass of Inf with a single switch in α at m=1.0. From 0.08 ≤ m ≤ 1.0, α = 1.35 and from 1.0 ≤ m ≤ Inf, α = 2.35.\n\nNotes\n\nThere is some setup necessary for quantile and other derived methods, so it is more efficient to call these methods directly with an array via the call signature quantile(d::BrokenPowerLaw{T}, x::AbstractArray{S}) rather than broadcasting over x. This behavior is now deprecated for quantile(d::Distributions.UnivariateDistribution, X::AbstractArray) in Distributions.jl.\n\nMethods\n\nBase.convert(::Type{BrokenPowerLaw{T}}, d::BrokenPowerLaw)\nminimum(d::BrokenPowerLaw)\nmaximum(d::BrokenPowerLaw)\npartype(d::BrokenPowerLaw)\neltype(d::BrokenPowerLaw)\nmean(d::BrokenPowerLaw)\nmedian(d::BrokenPowerLaw)\nvar(d::BrokenPowerLaw), may not function correctly for large mmax\nskewness(d::BrokenPowerLaw), may not function correctly for large mmax\nkurtosis(d::BrokenPowerLaw), may not function correctly for large mmax\npdf(d::BrokenPowerLaw,x::Real)\nlogpdf(d::BrokenPowerLaw,x::Real)\ncdf(d::BrokenPowerLaw,x::Real)\nccdf(d::BrokenPowerLaw,x::Real)\nquantile(d::BrokenPowerLaw{S},x::T) where {S,T<:Real}\nquantile!(result::AbstractArray,d::BrokenPowerLaw{S},x::AbstractArray{T}) where {S,T<:Real}\nquantile(d::BrokenPowerLaw{T},x::AbstractArray{S})\ncquantile(d::BrokenPowerLaw{S},x::T) where {S,T<:Real}\nrand(rng::AbstractRNG, d::BrokenPowerLaw,s...) \nOther methods from Distributions.jl should also work because BrokenPowerLaw <: AbstractIMF <: Distributions.ContinuousUnivariateDistribution. For example, rand!(rng::AbstractRNG, d::BrokenPowerLaw, x::AbstractArray).\n\n\n\n\n\n","category":"type"},{"location":"types/#InitialMassFunctions.LogNormalBPL","page":"Defined Types","title":"InitialMassFunctions.LogNormalBPL","text":"LogNormalBPL(μ::Real,σ::Real,α::AbstractVector{<:Real},breakpoints::AbstractVector{<:Real})\nLogNormalBPL(μ::Real,σ::Real,α::Tuple,breakpoints::Tuple)\n\nA LogNormal distribution at low masses, with a broken power law extension at high masses. This uses the natural log base like Distributions.LogNormal; if you have σ and μ in base 10, then multiply them both by log(10). Must have length(α) == length(breakpoints)-2. The probability distribution for this IMF model is\n\n    fracdn(m)dm = fracAx  exp left frac -left( log(x) - mu right)^22sigma^2 right\n\nfor m < breakpoints[2], with a broken power law extension above this mass. See BrokenPowerLaw for interface details; the α and breakpoints are the same here as there.\n\nArguments\n\nμ; see Distributions.LogNormal\nσ; see Distributions.LogNormal\nα; list of power law indices with length(α) == length(breakpoints)-2.\nbreakpoints; list of masses that signal breaks in the IMF. MUST BE SORTED and bracketed with breakpoints[1] being the minimum valid mass and breakpoints[end] being the maximum valid mass.\n\nExamples\n\nIf you want a LogNormalBPL with a characteristic mass of 0.5 solar masses, log10 standard deviation of 0.6, and a single power law extension with slope α=2.35 with a break at 1 solar mass, you would do LogNormalBPL(log(0.5),0.6*log(10),[2.35],[0.08,1.0,Inf] where we set the minimum mass to 0.08 and maximum mass to Inf. If, instead, you know that log10(m)=x, where m is the characteristic mass of the LogNormal component, you would do LogNormalBPL(x*log(10),0.6*log(10),[2.35],[0.08,1.0,Inf].\n\nNotes\n\nThere is some setup necessary for quantile and other derived methods, so it is more efficient to call these methods directly with an array via the call signature quantile(d::LogNormalBPL{T}, x::AbstractArray{S}) rather than broadcasting over x. This behavior is now deprecated for quantile(d::Distributions.UnivariateDistribution, X::AbstractArray) in Distributions.jl. \n\nMethods\n\nBase.convert(::Type{LogNormalBPL{T}}, d::LogNormalBPL)\nminimum(d::LogNormalBPL)\nmaximum(d::LogNormalBPL)\npartype(d::LogNormalBPL)\neltype(d::LogNormalBPL)\nmean(d::LogNormalBPL)\nmedian(d::LogNormalBPL)\npdf(d::LogNormalBPL,x::Real)\nlogpdf(d::LogNormalBPL,x::Real)\ncdf(d::LogNormalBPL,x::Real)\nccdf(d::LogNormalBPL,x::Real)\nquantile(d::LogNormalBPL{S},x::T) where {S,T<:Real}\nquantile!(result::AbstractArray,d::LogNormalBPL{S},x::AbstractArray{T}) where {S,T<:Real}\nquantile(d::LogNormalBPL{T},x::AbstractArray{S})\ncquantile(d::LogNormalBPL{S},x::T) where {S,T<:Real}\nrand(rng::AbstractRNG, d::LogNormalBPL,s...) \nOther methods from Distributions.jl should also work because LogNormalBPL <: AbstractIMF <: Distributions.ContinuousUnivariateDistribution. For example, rand!(rng::AbstractRNG, d::LogNormalBPL, x::AbstractArray).\n\n\n\n\n\n","category":"type"}]
}
