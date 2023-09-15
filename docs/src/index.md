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

We also provide a constructor for a single-power-law IMF,

```@docs
PowerLawIMF
```

which internally creates a truncated [Pareto](https://juliastats.org/Distributions.jl/stable/univariate/#Distributions.Pareto) distribution. Similarly, we provide a constructor for a single-component `LogNormal` IMF,

```@docs
LogNormalIMF
```