# Convenience Constructors for Published IMFs
We provide convenience constructors for published IMFs that can be called without arguments, or with positional arguments to set different minimum and maximum stellar masses. The provided constructors are

```@docs
Salpeter1955
Chabrier2001BPL
Kroupa2001
Chabrier2001LogNormal
Chabrier2003
Chabrier2003System
```

```@setup pdfcompare
# @setup ensures input and output are hidden from final output
using Distributions, InitialMassFunctions
import GR

# Display figures as SVGs
GR.inline("svg")

# Write function to make comparison plot of provided literature IMF PDFs
function pdfcompare(Mmin, Mmax, npoints; kws...)
    GR.figure(;
        title="Literature IMF PDFs",
	xlabel="M\$_\\odot\$",
	ylabel="PDF (dN/dM)",
	xlog=true,
	ylog=true,
	grid=false,
        backgroundcolor=0, # white instead of transparent background for dark Documenter scheme
	# font="Helvetica_Regular", # work around https://github.com/JuliaPlots/Plots.jl/issues/2596
	linewidth=2.0, # thicker lines
	size=(800,600),
	xlim=(Mmin, 10.0),
	ylim=(1e-4,20.0),
	kws...)

    imfs = [Salpeter1955(Mmin, Mmax),
            Chabrier2001BPL(Mmin, Mmax),
    	    Kroupa2001(Mmin, Mmax),
            Chabrier2001LogNormal(Mmin, Mmax),
	    Chabrier2003(Mmin, Mmax),
	    Chabrier2003System(Mmin, Mmax)]
    imf_labels = ["Salpeter1955",
                  "Chabrier2001BPL",
                  "Kroupa2001",
		  "Chabrier2001LogNormal",
		  "Chabrier2003",
		  "Chabrier2003System"]
    masses = exp10.(range(log10(Mmin), log10(Mmax); length=npoints))
    pdfs = [pdf.(imf, masses) for imf in imfs]

    # # Normalization test
    # using Test
    # import QuadGK: quadgk
    # let Mmin = 0.08, Mmax = 100.0,
    #     imfs = [Chabrier2001BPL(Mmin, Mmax),
    #             Kroupa2001(Mmin, Mmax),
    #             Chabrier2001LogNormal(Mmin, Mmax),
    #             Chabrier2003(Mmin, Mmax),
    #             Chabrier2003System(Mmin, Mmax)]
    #     @test all( isapprox.([quadgk(Base.Fix1(pdf, imf), Mmin, Mmax)[1] for imf in imfs],1) )
    # end
	    
    # # Plot first IMF
    # GR.plot(masses, Base.Fix1(pdf,Salpeter1955(Mmin, Mmax)))
    # # Hold open
    # GR.hold(true)
    # for imf in imfs
    #     GR.plot(masses, Base.Fix1(pdf, imf))
    # end
    # GR.hold(false)
    # # GR.legend(imf_labels...)

    # Switch to plotting all lines at same time, easier for legend
    # Call is GR.plot(x1, y1, x2, y2, ...and so on)
    # Legend labels are then keyword `labels` and legend location is
    # `location` with 1 = upper right, 2 = upper left, 
    # 3 = lower left, 4 = lower right, 11 = top right (outside of axes),
    # 12 = half right (outside of axes), 13 = bottom right (outside of axes),
    GR.plot( (collect((masses, i) for i in pdfs)...)...,
             labels=imf_labels, location=3)
end
```

Below is a comparison plot of the probability density functions of the literature IMFs we provide, all normalized to integrate to 1 over the initial mass range (0.08, 100.0) solar masses.

```@example pdfcompare
pdfcompare(0.08, 100.0, 1000) # hide
```
\

We also provide a constructor for a single-power-law IMF,

```@docs
PowerLawIMF
```

which internally creates a truncated [Pareto](https://juliastats.org/Distributions.jl/stable/univariate/#Distributions.Pareto) distribution. Similarly, we provide a constructor for a single-component `LogNormal` IMF,

```@docs
LogNormalIMF
```