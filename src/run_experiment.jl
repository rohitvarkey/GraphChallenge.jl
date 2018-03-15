#using Compose
using TimerOutputs

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
    p, count_log = result
    truth = load_truth_partition(num_nodes)
    #draw(SVG(plot_file, 32cm, 32cm), plotbeliefs(g, p))
    correctness_metrics = evaluate_partition(truth, p.b)
    agg_time = TimerOutputs.time(timer["timer"]["agglomerative_updates"])
    agg_bytes = TimerOutputs.allocated(timer["timer"]["agglomerative_updates"])
    iters = TimerOutputs.ncalls(timer["timer"]["agglomerative_updates"])
    nodal_time = TimerOutputs.time(timer["timer"]["nodal_updates"])
    nodal_bytes = TimerOutputs.allocated(timer["timer"]["nodal_updates"])
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
