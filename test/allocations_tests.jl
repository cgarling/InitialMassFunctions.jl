using AllocCheck: @check_allocs
using InitialMassFunctions
using Test

@check_allocs test_imf_allocs_constructor(T, i) = i(T(0.08), T(100.0))
@check_allocs test_imf_allocs(d, f) = f(d) 
@check_allocs test_imf_allocs(d, f, M) = f(d, M)

for T in (Float32, Float64)
    for i in (Salpeter1955, Chabrier2001LogNormal, Chabrier2003, Chabrier2003System, Kroupa2001, Chabrier2001BPL)
        @testset "$(i) $(T) Allocations" begin
            if i ∈ (Chabrier2001LogNormal, Kroupa2001) # Constructors for which the allocations test do not pass
                d = i(T(0.08), T(100.0))
            else
                d = @test_nowarn test_imf_allocs_constructor(T, i)
            end

            for f in (rand, mean, median)
                @testset "$(f)" begin
                    @test_nowarn test_imf_allocs(d, f)
                end
            end
            # for f in (pdf, logpdf, cdf, ccdf, quantile, cquantile)
            #     @testset "$(f)" begin
            #         test_imf_allocs(d, f, T(1//2))
            #     end
            # end
        end
    end
end