#using Compose
using TimerOutputs

function static_partition_experiment(T::Type, num_nodes::Int64; plot_file="blocks.svg")
    g = load_graph(num_nodes)
    timer = TimerOutput()
    result = @timeit timer "T $(num_nodes)" partition(T, g, num_nodes, timer)
    p, count_log = result
    truth = load_truth_partition(num_nodes)
    #draw(SVG(plot_file, 32cm, 32cm), plotbeliefs(g, p))
    correctness_metrics = evaluate_partition(truth, p.b)
    m = BenchmarkMetrics(
        T,
        nv(g),
        ne(g),
        p.B,
        PerformanceMetrics(
            TimerOutputs.tottime(timer),
            TimerOutputs.allocated(timer),
            timer
        ),
        correctness_metrics,
        count_log
    )
    p, m
end
