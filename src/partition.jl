function carry_out_best_merges(
    delta_entropy_for_each_block::Vector{Float64},
    best_merge_for_each_block::Vector{Int64},
    b::Vector{Int64}, num_blocks::Int64, num_blocks_to_merge::Int64
    )
    @show num_blocks_to_merge
    best_blocks = sortperm(delta_entropy_for_each_block)
    blocks_merged = 0
    counter = 1
    block_map = collect(1:num_blocks)
    while blocks_merged < num_blocks_to_merge
        mergeFrom = best_blocks[counter]
        mergeTo = block_map[best_merge_for_each_block[mergeFrom]]
        counter += 1
        if mergeTo != mergeFrom
            block_map[findin(block_map, mergeFrom)] = mergeTo
            b[findin(b, mergeFrom)] = mergeTo
            blocks_merged += 1
        end
    end
    remaining_blocks = unique(b)
    @show length(remaining_blocks)
    mapping = -ones(Int64, num_blocks)
    mapping[remaining_blocks] = 1:length(remaining_blocks)
    b = mapping[b]
    @show b
    return b
end

function agglomerative_updates(
    M, num_blocks::Int64, b::Vector{Int64}, d::Vector{Int64},
    d_in::Vector{Int64}, d_out::Vector{Int64};
    num_agg_proposals_per_block::Int64 = 10,
    num_block_reduction_rate::Float64 = 0.5
    )
    best_merge_for_each_block = fill(-1, num_blocks)
    delta_entropy_for_each_block = fill(Inf, num_blocks)
    for current_block = 1:num_blocks
        neighbors, k_out, k_in, k = compute_block_neighbors_and_degrees(
            M, current_block
        )
        for proposal_idx = 1:num_agg_proposals_per_block
            proposal = propose_new_partition_agg(
                M, current_block, b, num_blocks,
                d, neighbors
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
    # Get the new block assignments
    new_num_blocks = floor(Int64, num_blocks * num_block_reduction_rate)
    b = carry_out_best_merges(
        delta_entropy_for_each_block, best_merge_for_each_block, b,
        num_blocks, new_num_blocks
    )

    b, new_num_blocks
end

function partition(T::Type, num_nodes::Int64)
    g = SimpleWeightedDiGraph(num_nodes)
    load_graph!(g, num_nodes)
    partition(T, g, num_nodes)
end

function partition(T::Type, sampling_type::String, num_nodes::Int64)
    g = SimpleWeightedDiGraph(num_nodes)
    load_graph!(g, sampling_type, num_nodes, 1)
    partition(T, g, num_nodes)
end

function partition(T::Type, g::SimpleWeightedDiGraph, num_nodes::Int64)
    num_blocks = nv(g)
    partition = collect(1:num_blocks)

    β = 3
    num_agg_proposals_per_block = 10 # number of proposals per block
    num_block_reduction_rate = 0.5 # fraction of blocks to reduce until the golden ratio bracket is established

    # nodal partition updates parameters
    # maximum number of iterations
    max_num_nodal_itr = 10
    # stop iterating when the change in entropy falls below this fraction of the overall entropy
    # lowering this threshold results in more nodal update iterations and likely better performance, but longer runtime
    delta_entropy_threshold1 = 5e-4
    # threshold after the golden ratio bracket is established (typically lower to fine-tune to partition)
    delta_entropy_threshold2 = 1e-4
    # width of the moving average window for the delta entropy convergence criterion
    delta_entropy_moving_avg_window = 3

    old_overall_entropy = [Inf, Inf, Inf]
    best_partitions = Vector{Partition{Array{Int64, 2}}}(3)
    for i=1:3
        #Create dummy partitions
        best_partitions[i] = Partition(
            zeros(Int64, num_blocks, num_blocks),
            Inf,
            zeros(Int64, nv(g)),
            zeros(Int64, nv(g)),
            zeros(Int64, nv(g)),
            zeros(Int64, nv(g)),
            typemax(Int64)
        )
    end
    M = initialize_edge_counts(T, g, num_blocks, partition)
    d_out, d_in, d = compute_block_degrees(M, num_blocks)
    optimal_num_blocks_found = false

    while optimal_num_blocks_found == false
        println("Merging down from $num_blocks to $(floor(Int64, num_blocks * num_block_reduction_rate))")
        partition, num_blocks = agglomerative_updates(
            M, num_blocks, partition, d, d_in, d_out,
            num_agg_proposals_per_block = num_agg_proposals_per_block,
            num_block_reduction_rate = num_block_reduction_rate
        )
        @show partition, maximum(partition), num_blocks
        M = initialize_edge_counts(T, g, num_blocks, partition)
        d_out, d_in, d = compute_block_degrees(M, num_blocks)

        # compute the global entropy for MCMC convergence criterion
        @show overall_entropy = compute_overall_entropy(
            M, d_out, d_in, num_blocks, nv(g), ne(g)
        )

        total_num_nodal_moves = 0
        nodal_itr_delta_entropy = zeros(max_num_nodal_itr)
        for nodal_itr in 1:max_num_nodal_itr
            num_nodal_moves = 0
            nodal_itr_delta_entropy[nodal_itr] = 0

            for current_node in vertices(g)
                current_block  = partition[current_node]
                out_n = out_neighbors(g, current_node)
                in_n = in_neighbors(g, current_node)
                out_wts = [
                    floor(Int64, get_weight(g, current_node, n)) for n in out_n
                ]
                in_wts = [
                    floor(Int64, get_weight(g, n, current_node)) for n in in_n
                ]
                proposal = propose_new_partition_nodal(
                    M, current_block, partition, num_blocks,
                    d, vcat(out_n, in_n), vcat(out_wts, in_wts)
                )
                if (proposal != current_block)
                    #info("Performing nodal move on $current_node")
                    blocks_out_count_map = countmap(
                        partition[out_n], Distributions.weights(out_wts)
                    )
                    blocks_in_count_map = countmap(
                        partition[in_n], Distributions.weights(in_wts)
                    )

                    self_edge_weight = floor(
                        Int64, get_weight(g, current_node, current_node)
                    )

                    k_in = indegree(g, current_node)
                    k_out = outdegree(g, current_node)

                    M_r_row, M_r_col, M_s_row, M_s_col, Δ, p_accept =
                    evaluate_nodal_proposal(
                        M, current_block, proposal, num_blocks, β,
                        d, d_in, d_out, k_in + k_out, k_in, k_out,
                        self_edge_weight,
                        blocks_out_count_map, blocks_in_count_map
                    )

                    if rand() <= p_accept
                        total_num_nodal_moves += 1
                        num_nodal_moves += 1
                        nodal_itr_delta_entropy[nodal_itr] += Δ
                        M = update_partition(
                            M, current_block, proposal,
                            M_r_col, M_s_col, M_r_row, M_s_row
                        )
                        d_out, d_in, d = compute_block_degrees(M, num_blocks)
                        partition[current_node] = proposal
                    end
                end
            end
            # exit MCMC if the recent change in entropy falls below a small fraction of the overall entropy
           if (nodal_itr >= delta_entropy_moving_avg_window)
               window = (nodal_itr-delta_entropy_moving_avg_window+1):nodal_itr
               if all(isfinite.(old_overall_entropy)) == false # golden ratio bracket not yet established
                   if -mean(nodal_itr_delta_entropy[window]) <
                       delta_entropy_threshold1 * overall_entropy
                       println("Stopping at $nodal_itr iterations with $(nodal_itr_delta_entropy[nodal_itr]/overall_entropy)")
                       break
                   end
               else # golden ratio bracket is established. Fine-tuning partition.
                   if mean(itr_delta_entropy[window]) <
                       delta_entropy_threshold2 * overall_entropy
                       println("Stopping at $nodal_itr iterations with $(nodal_itr_delta_entropy[nodal_itr]/overall_entropy)")
                       break
                   end
               end
           end
        end

        @show partition, maximum(partition)
        # compute the global entropy for MCMC convergence criterion
        overall_entropy = compute_overall_entropy(
            M, d_out, d_in, num_blocks, nv(g), ne(g)
        )
        println("$total_num_nodal_moves nodal moves performed with entropy of $overall_entropy")
        new_partition, best_partitions, optimal_num_blocks_found, num_blocks_to_merge =
            prepare_for_partition_on_next_num_blocks(
                Partition(M, overall_entropy, partition, d, d_out, d_in, num_blocks),
                best_partitions, num_block_reduction_rate
            )
        @show new_partition
        @show best_partitions
        M = new_partition.M
        overall_entropy = new_partition.S
        partition = new_partition.b
        d = new_partition.d
        d_out = new_partition.d_out
        d_in = new_partition.d_in
        old_overall_entropy = [x.S for x in best_partitions]
        @show maximum(partition)
        #optimal_num_blocks_found = true #FIXME: Remove once all done
    end
end
