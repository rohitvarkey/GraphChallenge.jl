#using Compose

function static_partition_experiment(T::Type, num_nodes::Int64; plot_file="blocks.svg")
    g = load_graph(num_nodes)
    result, t, bytes, gctime, memallocs = @timed partition(T, g, num_nodes)
    p, count_log = result
    truth = load_truth_partition(num_nodes)
    #draw(SVG(plot_file, 32cm, 32cm), plotbeliefs(g, p))
    correctness_metrics = evaluate_partition(truth, p.b)
    m = BenchmarkMetrics(
        T,
        PerformanceMetrics(t, bytes, gctime, memallocs.allocd),
        correctness_metrics,
        count_log
    )
    p, m
end
