using InitialMassFunctions
using QuadGK
using Test

function test_bpl(d::BrokenPowerLaw)
    mmin,mmax = extrema(d)
    integral = quadgk(x->pdf(d,x),mmin,mmax)
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
end

# function test_type(T::Type,d::BrokenPowerLaw)
#     @test d isa BrokenPowerLaw{T}
#     @test partype(d) == T
#     @test convert(BrokenPowerLaw{Float32},d) isa BrokenPowerLaw{Float32}
#     @test convert(BrokenPowerLaw{T},d) === d

# end

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
        @test quantile(d,[0.5f0,0.75f0]) isa Vector{Float32}
        @test quantile(d,[0.5,0.75]) isa Vector{Float64}
        @test rand(d) isa Float32
        @test mean(d) isa Float32
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
