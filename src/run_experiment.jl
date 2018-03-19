#using Compose
using TimerOutputs
import GR
using StatPlots
using CSV

const BACKEND_NAMES = Dict(
    Array{Int64, 2} => "Dense matrix",
    InterblockEdgeCountStinger => "Stinger",
    InterblockEdgeCountSQLite => "SQLite",
    InterblockEdgeCountDictDict => "Dict Dict",
    InterblockEdgeCountVectorDict => "Vector Dict",
    InterblockEdgeCountPostgres => "Postgres",
)

function static_partition_experiment(T::Type, num_nodes::Int64; plot_file="blocks.svg")
    g = load_graph(num_nodes)
    timer = TimerOutput()
    result = @timeit timer "timer" partition(T, g, num_nodes, timer)
    p, count_log, df = result
    @show df
    truth = load_truth_partition(num_nodes)
    #draw(SVG(plot_file, 32cm, 32cm), plotbeliefs(g, p))
    correctness_metrics = evaluate_partition(truth, p.b)
    agg_time = TimerOutputs.time(timer["timer"]["agglomerative_updates"])
    agg_bytes = TimerOutputs.allocated(timer["timer"]["agglomerative_updates"])
    iters = TimerOutputs.ncalls(timer["timer"]["agglomerative_updates"])
    nodal_time = TimerOutputs.time(timer["timer"]["nodal_updates"])
    nodal_bytes = TimerOutputs.allocated(timer["timer"]["nodal_updates"])
    CSV.write("$(num_nodes)_sparsity.csv", df)
    gr()
    plt = @df df plot(:nodal_itr, :nnz, group=:M_size)
    savefig(plt, "$(num_nodes)_nnz")
    plt = @df df plot(:nodal_itr, :sparse_rows_changed, group=:M_size)
    savefig(plt, "$(num_nodes)_rows_changed")
    plt = @df df plot(:nodal_itr, :sparse_entries_flipped, group=:M_size)
    savefig(plt, "$(num_nodes)_flipped")
    m = BenchmarkMetrics(
        BACKEND_NAMES[T],
        nv(g),
        ne(g),
        p.B,
        PerformanceMetrics(
            TimerOutputs.tottime(timer),
            TimerOutputs.allocated(timer["timer"]),
            agg_time,
            agg_bytes,
            nodal_time,
            nodal_bytes,
            iters
        ),
        correctness_metrics,
        count_log,
        T
    )
    p, m, timer
end
