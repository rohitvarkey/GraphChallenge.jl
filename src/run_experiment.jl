using Compose

function static_partition_experiment(T::Type, num_nodes::Int64; plot_file="blocks.svg")
    g = load_graph(num_nodes)
    p = partition(T, g, num_nodes)
    t = load_truth_partition(num_nodes)
    draw(SVG(plot_file, 32cm, 32cm), plotbeliefs(g, p))
    evaluate_partition(t, p)
end
