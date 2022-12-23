using InitialMassFunctions
using QuadGK
using Test

function test_pl(d)
    mmin,mmax = extrema(d)
    integral = quadgk(x->pdf(d,x),mmin,mmax) # test that the pdf is properly normalized
    # integral = quadgk(x->exp(x)*pdf(d,exp(x)),log(mmin),log(mmax))
    @test integral[1] ≈ oneunit(integral[1]) rtol=1e-12 atol=integral[2] # test normalization
    # meanmass_gk = quadgk(x->x*pdf(d,x),mmin,mmax)
    # # meanmass_gk = quadgk(x->exp(x)^2*pdf(d,exp(x)),log(mmin),log(mmax))
    # @test meanmass_gk[1] ≈ mean(d) rtol=1e-12 atol=meanmass_gk[2] # test mean
    # var_gk = quadgk(x->x^2*pdf(d,x),mmin,mmax)
    # @test var_gk[1]-mean(d)^2 ≈ var(d) rtol=1e-12 atol=var_gk[2] # test variance
    # skew_gk = quadgk(x->x^3*pdf(d,x),mmin,mmax)
    # @test (skew_gk[1]-3*mean(d)*var(d)-mean(d)^3)/var(d)^(3/2) ≈ skewness(d) rtol=1e-5 atol=var_gk[2] # test skewness; higher error
    # dmean = mean(d)
    # var_d = var(d)
    # kurt_gk = quadgk(x->(x-dmean)^4 / var_d^2 * pdf(d,x),mmin,mmax)
    # @test kurt_gk[1] ≈ kurtosis(d) rtol=1e-10 atol=kurt_gk[2] # test kurtosis
    # test correctness of CDF, quantile, etc.
    test_points = [0.1,0.2,0.3,1.0,1.5,10.0,100.0]
    test_points = test_points[mmin .< test_points .< mmax]
    for i in test_points
        x = cdf(d,i)
        integral = quadgk(x->pdf(d,x),mmin,i)
        @test integral[1] ≈ x rtol=1e-12 atol=integral[2]
        @test quantile(d,x) ≈ i rtol=1e-12
        @test cquantile(d,x) ≈ quantile(d,1-x) rtol=1e-12
        @test logpdf(d,i) ≈ log(pdf(d,i)) rtol=1e-12
        @test x ≈ 1 - ccdf(d,i) rtol=1e-12
    end
end

function test_bpl(d::BrokenPowerLaw)
    mmin,mmax = extrema(d)
    integral = quadgk(x->pdf(d,x),mmin,mmax) # test that the pdf is properly normalized
    # integral = quadgk(x->exp(x)*pdf(d,exp(x)),log(mmin),log(mmax))
    @test integral[1] ≈ oneunit(integral[1]) rtol=1e-12 atol=integral[2] # test normalization
    meanmass_gk = quadgk(x->x*pdf(d,x),mmin,mmax)
    # meanmass_gk = quadgk(x->exp(x)^2*pdf(d,exp(x)),log(mmin),log(mmax))
    @test meanmass_gk[1] ≈ mean(d) rtol=1e-12 atol=meanmass_gk[2] # test mean
    var_gk = quadgk(x->x^2*pdf(d,x),mmin,mmax)
    @test var_gk[1]-mean(d)^2 ≈ var(d) rtol=1e-12 atol=var_gk[2] # test variance
    skew_gk = quadgk(x->x^3*pdf(d,x),mmin,mmax)
    @test (skew_gk[1]-3*mean(d)*var(d)-mean(d)^3)/var(d)^(3/2) ≈ skewness(d) rtol=1e-5 atol=var_gk[2] # test skewness; higher error
    dmean = mean(d)
    var_d = var(d)
    kurt_gk = quadgk(x->(x-dmean)^4 / var_d^2 * pdf(d,x),mmin,mmax)
    @test kurt_gk[1] ≈ kurtosis(d) rtol=1e-10 atol=kurt_gk[2] # test kurtosis
    # test correctness of CDF, quantile, etc.
    test_points = [0.1,0.2,0.3,1.0,1.5,10.0,100.0]
    test_points = test_points[mmin .< test_points .< mmax]
    for i in test_points
        x = cdf(d,i)
        integral = quadgk(x->pdf(d,x),mmin,i)
        @test integral[1] ≈ x rtol=1e-12 atol=integral[2]
        @test quantile(d,x) ≈ i rtol=1e-12
        @test cquantile(d,x) ≈ quantile(d,1-x) rtol=1e-12
        @test logpdf(d,i) ≈ log(pdf(d,i)) rtol=1e-12
        @test x ≈ 1 - ccdf(d,i) rtol=1e-12
    end
end

function test_lognormalbpl(d::LogNormalBPL)
    mmin,mmax = extrema(d)
    integral = quadgk(x->pdf(d,x),mmin,mmax)
    # integral = quadgk(x->exp(x)*pdf(d,exp(x)),log(mmin),log(mmax))
    @test integral[1] ≈ oneunit(integral[1]) rtol=1e-12 atol=integral[2] # test normalization
    meanmass_gk = quadgk(x->x*pdf(d,x),mmin,mmax)
    # meanmass_gk = quadgk(x->exp(x)^2*pdf(d,exp(x)),log(mmin),log(mmax))
    @test meanmass_gk[1] ≈ mean(d) rtol=1e-12 atol=meanmass_gk[2] # test mean
    # these are not defined yet for LogNormalBPL
    # var_gk = quadgk(x->x^2*pdf(d,x),mmin,mmax)
    # @test var_gk[1]-mean(d)^2 ≈ var(d) rtol=1e-12 atol=var_gk[2] # test variance
    # skew_gk = quadgk(x->x^3*pdf(d,x),mmin,mmax)
    # @test (skew_gk[1]-3*mean(d)*var(d)-mean(d)^3)/var(d)^(3/2) ≈ skewness(d) rtol=1e-5 atol=var_gk[2] # test skewness; higher error
    # dmean = mean(d)
    # var_d = var(d)
    # kurt_gk = quadgk(x->(x-dmean)^4 / var_d^2 * pdf(d,x),mmin,mmax)
    # @test kurt_gk[1] ≈ kurtosis(d) rtol=1e-10 atol=kurt_gk[2] # test kurtosis
    # test correctness of CDF, quantile, etc.
    test_points = [0.1,0.2,0.3,1.0,1.5,10.0,100.0]
    test_points = test_points[mmin .< test_points .< mmax]
    for i in test_points
        x = cdf(d,i)
        integral = quadgk(x->pdf(d,x),mmin,i)
        @test integral[1] ≈ x rtol=1e-12 atol=integral[2]
        @test quantile(d,x) ≈ i rtol=1e-12
        @test cquantile(d,x) ≈ quantile(d,1-x) rtol=1e-12
        @test logpdf(d,i) ≈ log(pdf(d,i)) rtol=1e-12
        @test x ≈ 1 - ccdf(d,i) rtol=1e-12
    end
end

# function test_type(T::Type,d::BrokenPowerLaw)
#     @test d isa BrokenPowerLaw{T}
#     @test partype(d) == T
#     @test convert(BrokenPowerLaw{Float32},d) isa BrokenPowerLaw{Float32}
#     @test convert(BrokenPowerLaw{T},d) === d

# end

@testset "PowerLawIMF" begin
    @testset "Float64" begin
        d = PowerLawIMF(2.35,0.08,100.0)
        # @test d isa PowerLawIMF{Float64}
        @test partype(d) == Float64
        # @test convert(PowerLawIMF{Float32},d) isa PowerLawIMF{Float32}
        # @test convert(PowerLawIMF{Float64},d) === d

        @test pdf(d,1.0) isa Float64
        @test pdf(d,1.0f0) isa Float64
        @test logpdf(d,1.0) isa Float64
        @test logpdf(d,1.0f0) isa Float64
        @test cdf(d,1.0) isa Float64
        @test cdf(d,1.0f0) isa Float64
        @test ccdf(d,1.0) isa Float64
        @test ccdf(d,1.0f0) isa Float64
        @test quantile(d,0.5) isa Float64
        @test quantile(d,0.5f0) isa Float64
        @test quantile.(d,[0.5f0,0.75f0]) isa Vector{Float64}
        @test quantile.(d,[0.5,0.75]) isa Vector{Float64}
        @test rand(d) isa Float64
        # @test mean(d) isa Float64
        @test median(d) isa Float64
    end
    @testset "Float32" begin
        d = PowerLawIMF(2.35f0,0.08f0,100.0f0)
        # @test d isa PowerLawIMF{Float32}
        @test partype(d) == Float32
        # @test convert(PowerLawIMF{Float64},d) isa PowerLawIMF{Float64}
        # @test convert(PowerLawIMF{Float32},d) === d

        @test pdf(d,1.0) isa Float64
        @test pdf(d,1.0f0) isa Float32
        @test logpdf(d,1.0) isa Float64
        @test logpdf(d,1.0f0) isa Float32
        @test cdf(d,1.0) isa Float64
        @test cdf(d,1.0f0) isa Float32
        @test ccdf(d,1.0) isa Float64
        @test ccdf(d,1.0f0) isa Float32
        @test quantile(d,0.5) isa Float64
        @test quantile(d,0.5f0) isa Float32
        @test quantile.(d,[0.5f0,0.75f0]) isa Vector{Float32}
        @test quantile.(d,[0.5,0.75]) isa Vector{Float64}
        # @test rand(d) isa Float32
        # @test mean(d) isa Float64
        @test median(d) isa Float32
    end
    @testset "Params" begin
        d = PowerLawIMF(2.35,0.08,100.0)
        @test minimum(d) == 0.08
        @test maximum(d) == 100.0
    end

    # Test named types of PLs
    @testset "Salpeter1955" begin
        d = Salpeter1955(0.08, 100.0)
        test_pl(d)
    end
end

@testset "BrokenPowerLaw" begin
    @testset "Float64" begin
        d = BrokenPowerLaw([1.3,2.35],[0.08,1.0,100.0])
        @test d isa BrokenPowerLaw{Float64}
        @test partype(d) == Float64
        @test convert(BrokenPowerLaw{Float32},d) isa BrokenPowerLaw{Float32}
        @test convert(BrokenPowerLaw{Float64},d) === d

        @test pdf(d,1.0) isa Float64
        @test pdf(d,1.0f0) isa Float64
        @test logpdf(d,1.0) isa Float64
        @test logpdf(d,1.0f0) isa Float64
        @test cdf(d,1.0) isa Float64
        @test cdf(d,1.0f0) isa Float64
        @test ccdf(d,1.0) isa Float64
        @test ccdf(d,1.0f0) isa Float64
        @test quantile(d,0.5) isa Float64
        @test quantile(d,0.5f0) isa Float64
        @test quantile(d,[0.5f0,0.75f0]) isa Vector{Float64}
        @test quantile(d,[0.5,0.75]) isa Vector{Float64}
        @test rand(d) isa Float64
        @test mean(d) isa Float64
        @test median(d) isa Float64
    end
    @testset "Float32" begin
        d = BrokenPowerLaw([1.3f0,2.35f0],[0.08f0,1.0f0,100.0f0])
        @test d isa BrokenPowerLaw{Float32}
        @test partype(d) == Float32
        @test convert(BrokenPowerLaw{Float32},d) === d
        @test convert(BrokenPowerLaw{Float64},d) isa BrokenPowerLaw{Float64}

        @test pdf(d,1.0) isa Float64
        @test pdf(d,1.0f0) isa Float32
        @test logpdf(d,1.0) isa Float64
        @test logpdf(d,1.0f0) isa Float32
        @test cdf(d,1.0) isa Float64
        @test cdf(d,1.0f0) isa Float32
        @test ccdf(d,1.0) isa Float64
        @test ccdf(d,1.0f0) isa Float32
        @test quantile(d,0.5) isa Float64
        @test quantile(d,0.5f0) isa Float32
        @test cquantile(d,0.5) isa Float64
        @test cquantile(d,0.5f0) isa Float32
        @test quantile(d,[0.5f0,0.75f0]) isa Vector{Float32}
        @test quantile(d,[0.5,0.75]) isa Vector{Float64}
        @test rand(d) isa Float32
        @test mean(d) isa Float32
        @test median(d) isa Float32
    end
    @testset "Params" begin
        d = BrokenPowerLaw([1.3,2.35],[0.08,1.0,100.0])
        @test minimum(d) == 0.08
        @test maximum(d) == 100.0
    end

    # test tuple constructor
    d = BrokenPowerLaw((1.3,2.35),(0.08,1.0,100.0))
    @test BrokenPowerLaw((1.3f0,2.35f0),(0.08,1.0,100.0)) isa BrokenPowerLaw{Float64}
    
    @test BrokenPowerLaw((1.3f0,2.35f0),(0.08f0,1.0f0,100.0f0)) isa BrokenPowerLaw{Float32}

    # test named types of BPLs
    @testset "Chabrier2001BPL" begin
        d = Chabrier2001BPL(0.08,100.0)
        test_bpl(d)
    end
    @testset "Kroupa2001" begin
        d = Kroupa2001(0.08,100.0)
        test_bpl(d)
    end
end

##############################################################################

@testset "LogNormalBPL" begin
    @testset "Float64" begin
        d = LogNormalBPL(-5.0,1.5,[2.35],[0.08,1.0,100.0])
        @test d isa LogNormalBPL{Float64}
        @test partype(d) == Float64
        @test convert(LogNormalBPL{Float32},d) isa LogNormalBPL{Float32}
        @test convert(LogNormalBPL{Float64},d) === d

        @test pdf(d,1.0) isa Float64
        @test pdf(d,1.0f0) isa Float64
        @test logpdf(d,1.0) isa Float64
        @test logpdf(d,1.0f0) isa Float64
        @test cdf(d,1.0) isa Float64
        @test cdf(d,1.0f0) isa Float64
        @test ccdf(d,1.0) isa Float64
        @test ccdf(d,1.0f0) isa Float64
        @test quantile(d,0.5) isa Float64
        @test quantile(d,0.5f0) isa Float64
        @test cquantile(d,0.5) isa Float64
        @test cquantile(d,0.5f0) isa Float64
        @test quantile(d,[0.5f0,0.75f0]) isa Vector{Float64}
        @test quantile(d,[0.5,0.75]) isa Vector{Float64}
        @test rand(d) isa Float64
        @test mean(d) isa Float64
    end
    @testset "Float32" begin
        d = LogNormalBPL(-5.0f0,1.5f0,[2.35f0],[0.08f0,1.0f0,100.0f0])
        @test d isa LogNormalBPL{Float32}
        @test partype(d) == Float32
        @test convert(LogNormalBPL{Float32},d) === d
        @test convert(LogNormalBPL{Float64},d) isa LogNormalBPL{Float64}

        @test pdf(d,1.0) isa Float64
        @test pdf(d,1.0f0) isa Float32
        @test logpdf(d,1.0) isa Float64
        @test logpdf(d,1.0f0) isa Float32
        @test cdf(d,1.0) isa Float64
        @test cdf(d,1.0f0) isa Float32
        @test ccdf(d,1.0) isa Float64
        @test ccdf(d,1.0f0) isa Float32
        @test quantile(d,0.5) isa Float64
        @test quantile(d,0.5f0) isa Float32
        @test quantile(d,[0.5f0,0.75f0]) isa Vector{Float32}
        @test quantile(d,[0.5,0.75]) isa Vector{Float64}
        @test rand(d) isa Float32
        @test mean(d) isa Float32
    end
    @testset "Params" begin
        d = LogNormalBPL(-5.0,1.5,[2.35],[0.08,1.0,100.0])
        @test minimum(d) == 0.08
        @test maximum(d) == 100.0
    end

    # test tuple constructor
    d = LogNormalBPL(-5.0,1.5,(2.35,),(0.08,1.0,100.0))
    @test LogNormalBPL(-5.0f0,1.5f0,(2.35f0,),(0.08,1.0,100.0)) isa LogNormalBPL{Float64}
    @test LogNormalBPL(-5.0,1.5f0,(2.35f0,),(0.08f0,1.0f0,100.0f0)) isa LogNormalBPL{Float64}

    # test named types of BPLs
    @testset "Chabrier2003" begin
        d = Chabrier2003(0.08,100.0)
        test_lognormalbpl(d)
    end
end
