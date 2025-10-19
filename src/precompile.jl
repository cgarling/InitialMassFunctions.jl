using PrecompileTools: @setup_workload, @compile_workload

@setup_workload begin
    @compile_workload begin
        for i in (Chabrier2003,Chabrier2001BPL,Chabrier2001LogNormal,Kroupa2001,Salpeter1955)
            for T in (Float32, Float64)
                I = i(T(0.08), T(100.0))
                for f in (rand, mean, median)
                    f(I)
                end
                for f in (pdf, logpdf, cdf, ccdf, quantile, cquantile)
                    f(I, T(0.5))
                end
            end
        end
    end
end
