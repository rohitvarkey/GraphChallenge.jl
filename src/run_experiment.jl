using Compose

function static_partition_experiment(T::Type, num_nodes::Int64; plot_file="blocks.svg")
    g = load_graph(num_nodes)
    p, t, bytes, gctime, memallocs = @timed partition(T, g, num_nodes)
    truth = load_truth_partition(num_nodes)
    #draw(SVG(plot_file, 32cm, 32cm), plotbeliefs(g, p))
    correctness_metrics = evaluate_partition(truth, p)
    p, t, bytes, correctness_metrics
end
