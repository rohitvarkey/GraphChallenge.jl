

function agglomerative_updates(
    M, num_blocks::Int64, b::Vector{Int64}, d::Vector{Int64},
    d_in::Vector{Int64}, d_out::Vector{Int64};
    num_agg_proposals_per_block::Int64 = 10
    )
    best_merge_for_each_block = fill(-1, num_blocks)
    delta_entropy_for_each_block = fill(Inf, num_blocks)
    for current_block = 1:num_blocks
        out_neighbors = findn(M[:, current_block])
        in_neighbors = findn(M[current_block, :])
        neighbors = collect(Set(out_neighbors) ∪ Set(in_neighbors))
        k_out = sum(M[:, out_neighbors])
        k_in = sum(M[in_neighbors, :])
        k = k_out + k_in
        for proposal_idx = 1:num_agg_proposals_per_block
            proposal = propose_new_partition(
                M, current_block, b, num_blocks,
                d, neighbors, true
            )
            ∇ = evaluate_proposal_agg(
                M, current_block, proposal, num_blocks, d, d_in, d_out,
                k, k_in, k_out
            )
            if ∇ < delta_entropy_for_each_block[current_block]
                best_merge_for_each_block[current_block] = proposal
                delta_entropy_for_each_block[current_block] = ∇
            end
        end
    end
    @show best_merge_for_each_block
    @show delta_entropy_for_each_block
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

    optimal_num_blocks_found = false
    agglomerative_updates(M, num_blocks, partition, d, d_in, d_out,
    num_agg_proposals_per_block=num_agg_proposals_per_block)

end
