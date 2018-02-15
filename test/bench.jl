using GraphChallenge
using CSV
using DataFrames

function bench(T, nv, seed)
    println("Running experiment for $T with $nv nodes")
    srand(seed)
    _, metrics = static_partition_experiment(T, nv);
    metrics
end

function run_benchmarks(types, num_nodes, seed)
    # Warm the JIT by running the smallest size experiment on the type.
    for T in types
        _ = bench(T, 50, seed)
    end
    # Actually run and collect the metrics
    metrics = [bench(T, n, seed) for n in num_nodes for T in types]
    metrics
end

function bench(
    types = [Array{Int64, 2}, InterblockEdgeCountStinger, InterblockEdgeCountDictDict, InterblockEdgeCountVectorDict, InterblockEdgeCountPostgres],
    num_nodes = [50, 100, 500, 1000];
    seed = 50, output_prefix = "bench_results"
    )
    metrics = run_benchmarks(types, num_nodes, seed)
    df = convert(DataFrame, metrics)
    CSV.write("$(output_prefix)_$(Dates.now()).csv", df)
    df
end
