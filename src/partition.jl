

function agglomerative_updates(M, g::SimpleWeightedDiGraph, )

end


function partition(T::Type, sampling_type::String, num_nodes::Int64)
    g = SimpleWeightedDiGraph(num_nodes)
    load_graph!(g, sampling_type, num_nodes, 1)
    num_blocks = nv(g)
    partition = collect(1:num_blocks)

    β = 3
    num_agg_proposals_per_block = 10 # number of proposals per block
    num_block_reduction_rate = 0.5 # fraction of blocks to reduce until the golden ratio bracket is established

    # nodal partition updates parameters
    # maximum number of iterations
    max_num_nodal_itr = 100
    # stop iterating when the change in entropy falls below this fraction of the overall entropy
    # lowering this threshold results in more nodal update iterations and likely better performance, but longer runtime
    delta_entropy_threshold1 = 5e-4
    # threshold after the golden ratio bracket is established (typically lower to fine-tune to partition)
    delta_entropy_threshold2 = 1e-4
    # width of the moving average window for the delta entropy convergence criterion
    delta_entropy_moving_avg_window = 3

    M, d_out, d_in, d = initialize_edge_counts(
        T, g, num_blocks, partition
    )
end
