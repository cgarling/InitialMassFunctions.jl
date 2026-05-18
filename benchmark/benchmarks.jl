using BenchmarkTools
using PrettyTables: pretty_table
using InitialMassFunctions
using Random: rand!
using StaticArrays: SVector

function show_benchmarks(results)
    # Collect results
    sorted  = sort(collect(results), by=first)
    names   = [k for (k,_) in sorted]
    trials  = [v for (_,v) in sorted]

    # Pack into matrix
    data = hcat(
        names,
        [BenchmarkTools.prettytime(median(t).time) for t in trials],
        [BenchmarkTools.prettymemory(median(t).memory) for t in trials],
        [median(t).allocs for t in trials]
    )

    # Make pretty table
    pretty_table(data;
        column_labels = ["Benchmark", "Median Time", "Memory", "Allocs"],
        alignment     = [:l, :r, :r, :r]
    )
end

const SUITE = BenchmarkGroup()

# BrokenPowerLaw benchmarks
function bench_broken_power_law(T; samples::Int=1000)
    group = BenchmarkGroup()
    # Test with both regular arrays and static vectors to check for any performance differences
    d = BrokenPowerLaw(T[1.35, 2.35], T[0.08, 1.0, 100.0])
    d_s = BrokenPowerLaw(SVector{2,T}(1.35, 2.35), SVector{3,T}(0.08, 1.0, 100.0))
    for (name, dist) in (("", d), (" svec", d_s))
        group["rand"*name] = @benchmarkable rand($dist) samples=samples
        group["pdf"*name] = @benchmarkable pdf($dist, $(T(0.5))) samples=samples
        group["logpdf"*name] = @benchmarkable logpdf($dist, $(T(0.5))) samples=samples
        group["cdf"*name] = @benchmarkable cdf($dist, $(T(0.5))) samples=samples
        group["quantile"*name] = @benchmarkable quantile($dist, $(T(0.5))) samples=samples

        xs = T[0.5, 0.25, 0.75, 0.1, 0.9]
        ys = similar(xs)
        group["quantile!"*name] = @benchmarkable quantile!($ys, $dist, $xs) samples=samples
        group["rand! array"*name] = @benchmarkable rand!($dist, $(T[1.0,2.0])) samples=samples
    end
    tune!(group) # Some calls are very fast (~10 ns), so we need to tune to get accurate measurements
    return group
end

# LogNormalBPL benchmarks
function bench_log_normal_bpl(T; samples::Int=1000)
    group = BenchmarkGroup()
    d = LogNormalBPL(T(-5.0), T(1.5), T[2.35], T[0.08, 1.0, 100.0])
    d_s = LogNormalBPL(T(-5.0), T(1.5), SVector{1,T}(2.35), SVector{3,T}(0.08, 1.0, 100.0))
    for (name, dist) in (("", d), (" svec", d_s))
        group["rand"*name] = @benchmarkable rand($dist) samples=samples
        group["pdf"*name] = @benchmarkable pdf($dist, $(T(0.5))) samples=samples
        group["logpdf"*name] = @benchmarkable logpdf($dist, $(T(0.5))) samples=samples
        group["cdf"*name] = @benchmarkable cdf($dist, $(T(0.5))) samples=samples
        group["quantile"*name] = @benchmarkable quantile($dist, $(T(0.5))) samples=samples

        xs = T[0.5, 0.25, 0.75, 0.1, 0.9]
        ys = similar(xs)
        group["quantile!"*name] = @benchmarkable quantile!($ys, $dist, $xs) samples=samples

        group["rand! array"*name] = @benchmarkable rand!($dist, $(T[1.0,2.0])) samples=samples
    end
    tune!(group) # Some calls are very fast (~10 ns), so we need to tune to get accurate measurements
    return group
end

for T in (Float64, Float32)
    SUITE["BrokenPowerLaw $T"] = bench_broken_power_law(T)
    SUITE["LogNormalBPL $T"] = bench_log_normal_bpl(T)
end

# If not on CI, we'll show a nice table
if get(ENV, "CI", "false") == "false"
     results = run(SUITE, verbose=true)

    for (name, group) in results
        println("\n=== $name ===")
        display(show_benchmarks(group))
    end
end
