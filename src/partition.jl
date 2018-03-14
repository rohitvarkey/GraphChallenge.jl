using Base.Threads
using TimerOutputs

# Fall through
function initial_setup(x)
    nothing
end

function initialize_edge_counts(T, g, B, b, config, count_log)
    initialize_edge_counts(T, g, B, b, count_log)
end
function zeros_interblock_edge_matrix(M, size, config::Void)
    zeros_interblock_edge_matrix(M, size)
end

function carry_out_best_merges(
    delta_entropy_for_each_block::Vector{Float64},
    best_merge_for_each_block::Vector{Int64},
    b::Vector{Int64}, num_blocks::Int64, num_blocks_to_merge::Int64
    )
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
    remaining_blocks::Vector{Int64} = sort(unique(b))
    mapping = -ones(Int64, num_blocks)
    mapping[remaining_blocks] = 1:length(remaining_blocks)
    return mapping[b]
end

function agglomerative_updates_kernel{T}(
    p::Partition{T},
    current_block::Int64,
    best_merge_for_each_block::Vector{Int64},
    delta_entropy_for_each_block::Vector{Float64},
    num_agg_proposals_per_block::Int64,
    count_log::CountLog
    )
    neighbors, k_out, k_in, k = compute_block_neighbors_and_degrees(
        p, current_block, count_log
    )
    for proposal_idx = 1:num_agg_proposals_per_block
        proposal = propose_new_partition_agg(p, current_block, neighbors, count_log)
        Δ = evaluate_proposal_agg(
            p, current_block, proposal, k, k_in, k_out, count_log
        )
        if Δ < delta_entropy_for_each_block[current_block]
            best_merge_for_each_block[current_block] = proposal
            delta_entropy_for_each_block[current_block] = Δ
        end
    end
end

function agglomerative_updates{T}(
    p::Partition{T}, num_blocks_to_merge::Int64, config, count_log::CountLog,
    num_agg_proposals_per_block::Int64 = 10,
    num_block_reduction_rate::Float64 = 0.5
    )
    num_blocks = p.B
    best_merge_for_each_block = fill(-1, num_blocks)
    delta_entropy_for_each_block = fill(Inf, num_blocks)
    @threads for current_block = 1:num_blocks
        agglomerative_updates_kernel(
            p, current_block, best_merge_for_each_block,
            delta_entropy_for_each_block, num_agg_proposals_per_block,
            count_log
        )
    end
    # Get the new block assignments
    new_num_blocks = num_blocks - num_blocks_to_merge
    b = carry_out_best_merges(
        delta_entropy_for_each_block, best_merge_for_each_block, p.b,
        num_blocks, num_blocks_to_merge
    )

    Partition(T, p.g, b, num_blocks - num_blocks_to_merge, config, count_log = count_log)
end

function prepare_for_partition_on_next_num_blocks{T}(
    current_partition::Partition{T},
    best_partitions::Vector{Partition{T}},
    B_rate::Float64
    )
    optimal_B_found = false
    if current_partition.S <= best_partitions[2].S
        if best_partitions[2].B > current_partition.B
            best_partitions[1] = best_partitions[2]
        else
            best_partitions[3] = best_partitions[2]
        end
        best_partitions[2] = current_partition
    elseif best_partitions[2].B > current_partition.B
        best_partitions[3] = current_partition
    else
        best_partitions[1] = current_partition
    end

    if (best_partitions[3].S == Inf)
        B_to_merge = floor(Int64, current_partition.B * B_rate)
        if B_to_merge == 0
            optimal_B_found = true
        end
        partition = deepcopy(best_partitions[2])
    else
        if best_partitions[1].B - best_partitions[3].B == 2
            optimal_B_found = true
            partition = deepcopy(best_partitions[2])
            B_to_merge = 0
        else
            if (best_partitions[1].B - best_partitions[2].B) >=
                (best_partitions[2].B - best_partitions[3].B)  # the higher segment in the bracket is bigger
                index = 1
            else  # the lower segment in the bracket is bigger
                index = 2
            end
            next_B_to_try = best_partitions[index + 1].B +
                round(Int64,
                (best_partitions[index].B - best_partitions[index + 1].B) * 0.618)
            B_to_merge = best_partitions[index].B - next_B_to_try
            partition = deepcopy(best_partitions[index])
        end
    end
    return partition, best_partitions, optimal_B_found, B_to_merge
end

function propose_new_partition_agg{T}(
    p::Partition{T}, r::Int64, neighbors::Vector{Int64}, count_log::CountLog
    )
    # Pick a neighbor randomly
    if length(neighbors) == 0
        candidates = Set(1:p.B)
        # Force to be different than r.
        pop!(candidates, r)
        s = sample(collect(candidates))
        return s
    end
    rand_neighbor = sample(neighbors, Distributions.weights(p.d[neighbors]./sum(p.d)))
    u = p.b[rand_neighbor]
    if rand() < p.B / (p.d[u] + p.B)
        candidates = Set(1:p.B)
        pop!(candidates, r)
        s = sample(collect(candidates))
    else
        multinomial_prob = compute_multinomial_probs(p, u, count_log)
        multinomial_prob[r] = 0
        if sum(multinomial_prob) == 0
            candidates = Set(1:p.B)
            pop!(candidates, r)
            s = sample(collect(candidates))
            return s
        else
            # Normalize back
            multinomial_prob /= sum(multinomial_prob)
        end
        s = findn(rand(Multinomial(1, multinomial_prob)))[1]
    end
    return s
end

function propose_new_partition_nodal{T}(
    p::Partition{T}, r::Int64, neighbors::Vector{Int64}, wts::Vector{Int64},
    count_log::CountLog
    )

    # Pick a neighbor randomly
    if length(neighbors) == 0
        candidates = Set(1:p.B)
        s = sample(collect(candidates))
        return s
    end
    rand_neighbor = sample(neighbors, Distributions.weights(wts./sum(wts)))
    u = p.b[rand_neighbor]
    if rand() < p.B / (p.d[u] + p.B)
        candidates = Set(1:p.B)
        pop!(candidates, r)
        s = sample(collect(candidates))
    else
        multinomial_prob = compute_multinomial_probs(p, u, count_log)
        #multinomial_prob[r] = 0
        if sum(multinomial_prob) == 0
            candidates = Set(1:B)
            pop!(candidates, r)
            s = sample(collect(candidates))
            return s
        else
            # Normalize back
            multinomial_prob /= sum(multinomial_prob)
        end
        s = findn(rand(Multinomial(1, multinomial_prob)))[1]
    end
    return s
end

function evaluate_proposal_agg{T}(
    p::Partition{T}, r::Int64, s::Int64, k::Int64, k_in::Int64, k_out::Int64,
    count_log::CountLog
    )
    M_r_row, M_r_col, M_s_row, M_s_col =
        compute_new_matrix_agglomerative(p, r, s, count_log)
    new_degrees = [copy(degrees) for degrees in [p.d_out, p.d_in, p.d]]
    for (new_d, degree) in zip(new_degrees, [k_out, k_in, k])
        new_d[r] -= degree
        new_d[s] += degree
    end
    compute_delta_entropy(
        p, r, s, M_r_col, M_s_col, M_r_row, M_s_row, new_degrees[1], new_degrees[2],
        count_log
    )
end

function evaluate_nodal_proposal{T}(
    p::Partition{T}, r::Int64, s::Int64,β::Int64,
    k::Int64, k_in::Int64, k_out::Int64, self_edge_weight::Int64,
    blocks_out_count_map, blocks_in_count_map, count_log
    )

    M_r_row, M_r_col, M_s_row, M_s_col = compute_new_matrix(
        p, r, s, blocks_out_count_map, blocks_in_count_map, self_edge_weight,
        count_log
    )

    new_degrees = [copy(degrees) for degrees in [p.d_out, p.d_in, p.d]]
    for (new_d, degree) in zip(new_degrees, [k_out, k_in, k])
        new_d[r] -= degree
        new_d[s] += degree
    end

    hastings_correction = compute_hastings_correction(
        s, p, M_r_row, M_r_col, new_degrees[3],
        blocks_out_count_map, blocks_in_count_map,
        count_log
    )

    Δ = compute_delta_entropy(
            p, r, s, M_r_col, M_s_col, M_r_row, M_s_row,
            new_degrees[1], new_degrees[2], count_log
        )

    p_accept = min(exp(-β * Δ) * hastings_correction, 1)

    #println("p_accept: $(p_accept), Δ: $Δ, β: $β, H: $(hastings_correction), exp: $(exp(-β*Δ)*hastings_correction)")

    M_r_row, M_r_col, M_s_row, M_s_col, Δ, p_accept
end

function nodal_update_kernel{T}(
    current_node::Int64,
    current_partition::Partition{T},
    g::SimpleWeightedDiGraph,
    β::Int64,
    b_new::Vector{Int64},
    out_n::Vector{Int64},
    in_n::Vector{Int64},
    num_nodal_moves::Vector{Int64},
    nodal_itr_delta_entropy::Vector{Float64},
    nodal_itr::Int64,
    count_log::CountLog
    )
    current_block = current_partition.b[current_node]
    out_wts = [
        floor(Int64, get_weight(g, current_node, n))::Int64 for n in out_n
    ]
    in_wts = [
        floor(Int64, get_weight(g, n, current_node))::Int64 for n in in_n
    ]
    proposal = propose_new_partition_nodal(
        current_partition, current_block, vcat(out_n, in_n),
        vcat(out_wts, in_wts), count_log
    )
    #println("proposal: $proposal, current:$current_block")
    if (proposal != current_block)
        #println("Performing nodal move on $current_node")
        blocks_out_count_map = countmap(
            current_partition.b[out_n], Distributions.weights(out_wts)
        )
        blocks_in_count_map = countmap(
            current_partition.b[in_n], Distributions.weights(in_wts)
        )

        self_edge_weight = floor(
            Int64, get_weight(g, current_node, current_node)
        )

        k_in::Int64 = length(in_n)
        k_out::Int64 = LightGraphs.outdegree(g, current_node)

        M_r_row, M_r_col, M_s_row, M_s_col, Δ::Float64, p_accept::Float64 =
        evaluate_nodal_proposal(
            current_partition, current_block, proposal, β,
            k_in + k_out, k_in, k_out, self_edge_weight,
            blocks_out_count_map, blocks_in_count_map,
            count_log
        )
        random_value = rand(Uniform())
        if random_value <= p_accept
            num_nodal_moves[Threads.threadid()] += 1
            nodal_itr_delta_entropy[nodal_itr] += Δ
            #update proposal
            b_new[current_node] = proposal
        end
    end
end

function update_partition{T}(
    p::Partition{T}, current_node::Int64, r::Int64, s::Int64, M_r_col,
    M_s_col, M_r_row, M_s_row, count_log
    )
    p.M = update_partition(p.M, r, s, M_r_col, M_s_col, M_r_row, M_s_row, count_log)
    p.b[current_node] = s
    p.d_out, p.d_in, p.d = compute_block_degrees(p.M, p.B, count_log)
    p.S = compute_overall_entropy(
        p.M, p.d_out, p.d_in, p.B, nv(p.g), ne(p.g), count_log
    )
    p
end

function partition(T::Type, num_nodes::Int64)
    g = load_graph(num_nodes)
    partition(T, g, num_nodes)
end

function partition(T::Type, sampling_type::String, num_nodes::Int64)
    g = SimpleWeightedDiGraph(num_nodes)
    load_graph!(g, sampling_type, num_nodes, 1)
    partition(T, g, num_nodes)
end

function partition(T::Type, g::SimpleWeightedDiGraph, num_nodes::Int64)
    timer = TimerOutput()
    @timeit timer "$T $num_nodes" p = partition(T, g, num_nodes, timer)
    p, timer
end

function partition(T::Type, g::SimpleWeightedDiGraph, num_nodes::Int64, timer::TimerOutput)

    num_blocks = nv(g)
    b = collect(1:num_blocks)

    vertex_in_neighbors = inneighbors(g)

    β = 3
    num_agg_proposals_per_block = 100 # number of proposals per block
    num_block_reduction_rate = 0.5 # fraction of blocks to reduce until the golden ratio bracket is established

    # nodal partition updates parameters
    # maximum number of iterations
    max_num_nodal_itr = 1000
    # stop iterating when the change in entropy falls below this fraction of the overall entropy
    # lowering this threshold results in more nodal update iterations and likely better performance, but longer runtime
    delta_entropy_threshold1 = 5e-8
    # threshold after the golden ratio bracket is established (typically lower to fine-tune to partition)
    delta_entropy_threshold2 = 1e-10
    # width of the moving average window for the delta entropy convergence criterion
    delta_entropy_moving_avg_window = 3

    old_overall_entropy = [Inf, Inf, Inf]
    best_partitions = Vector{Partition{T}}(3)

    count_log = CountLog()

    config = initial_setup(T)

    for i=1:3
        #Create dummy partitions
        best_partitions[i] = Partition(
            zeros_interblock_edge_matrix(T, nv(g), config),
            g,
            Inf,
            zeros(Int64, 0),
            zeros(Int64, 0),
            zeros(Int64, 0),
            zeros(Int64, 0),
            typemax(Int64),
            count_log
        )
    end


    current_partition::Partition{T} = Partition(T, g, b, length(b), config)

    optimal_num_blocks_found = false
    num_blocks_to_merge = ceil(Int64, num_blocks * num_block_reduction_rate)

    total_num_nodal_moves::Int64 = 0
    while optimal_num_blocks_found == false
        println("Merging down from $(current_partition.B) to $(current_partition.B - num_blocks_to_merge)")
        @timeit timer "agglomerative_updates" current_partition = agglomerative_updates(
            current_partition, num_blocks_to_merge, config, count_log,
            num_agg_proposals_per_block, num_block_reduction_rate
        )

        nodal_itr_delta_entropy = zeros(max_num_nodal_itr)
        @timeit timer "nodal_updates" for nodal_itr in 1:max_num_nodal_itr
            num_nodal_moves = zeros(Int64, Threads.nthreads())
            nodal_itr_delta_entropy[nodal_itr] = 0
            b_new::Vector{Int64} = copy(current_partition.b)

            @threads for current_node::Int64 in vertices(g)
                nodal_update_kernel(
                    current_node, current_partition, g, β, b_new,
                    outneighbors(g, current_node), vertex_in_neighbors[current_node],
                    num_nodal_moves, nodal_itr_delta_entropy, nodal_itr, count_log
                )
            end
            # Sequential Aggregation of updates
            for (vertex::Int64, new_block::Int64) in enumerate(b_new)
                current_block = current_partition.b[vertex]
                if new_block != current_block
                    # perform update
                    # TODO: More efficient way of doing this?
                    out_n::Vector{Int64} = outneighbors(g, vertex)
                    in_n::Vector{Int64} = vertex_in_neighbors[vertex]
                    out_wts = [
                        floor(Int64, get_weight(g, vertex, n))::Int64 for n in out_n
                    ]
                    in_wts = [
                        floor(Int64, get_weight(g, n, vertex))::Int64 for n in in_n
                    ]
                    blocks_out_count_map = countmap(
                        current_partition.b[out_n], Distributions.weights(out_wts)
                    )
                    blocks_in_count_map = countmap(
                        current_partition.b[in_n], Distributions.weights(in_wts)
                    )
                    self_edge_weight = floor(
                        Int64, get_weight(g, vertex, vertex)
                    )
                    M_r_row, M_r_col, M_s_row, M_s_col = compute_new_matrix(
                        current_partition, current_block, new_block,
                        blocks_out_count_map, blocks_in_count_map,
                        self_edge_weight, count_log
                    )
                    current_partition.M = update_partition(
                        current_partition.M, current_block, new_block,
                        M_r_col, M_s_col, M_r_row, M_s_row, count_log
                    )
                    current_partition.b[vertex] = new_block
                    #println("Updated $vertex from $(current_block) to $(new_block)")
                    #println("Partition is $(current_partition)")
                end
            end
            compute_block_degrees(current_partition, count_log)
            compute_overall_entropy(current_partition, count_log)
            #println("Partition after recomputations is $(current_partition)")
            # exit MCMC if the recent change in entropy falls below a small fraction of the overall entropy
            println("Itr: $nodal_itr, nodal moves: $(sum(num_nodal_moves)), Δ: $(nodal_itr_delta_entropy[nodal_itr]), fraction: $(-nodal_itr_delta_entropy[nodal_itr]/current_partition.S)")
            if (nodal_itr >= delta_entropy_moving_avg_window)
               window::UnitRange{Int64} = (nodal_itr-delta_entropy_moving_avg_window+1):nodal_itr
               println("$(-mean(nodal_itr_delta_entropy[window])), $(delta_entropy_threshold1 * current_partition.S), $(current_partition.S)")
               if all(isfinite.(old_overall_entropy)) == false # golden ratio bracket not yet established
                   if -mean(nodal_itr_delta_entropy[window]) <
                       delta_entropy_threshold1 * current_partition.S
                       break
                   end
               else # golden ratio bracket is established. Fine-tuning partition.
                   if -mean(nodal_itr_delta_entropy[window]) <
                       delta_entropy_threshold2 * current_partition.S
                       break
                   end
               end
           end
        end

        println("$total_num_nodal_moves nodal moves performed with entropy of $(current_partition.S)")
        println("")
        new_partition, best_partitions, optimal_num_blocks_found, num_blocks_to_merge =
            prepare_for_partition_on_next_num_blocks(
                current_partition, best_partitions, num_block_reduction_rate
            )

        @show count_log

        old_overall_entropy = [x.S for x in best_partitions]
        if optimal_num_blocks_found == true
            println("Total number of nodal moves: ", total_num_nodal_moves)
            println("Final partition: ", new_partition)
            println("Best partition :", new_partition.b, "Num blocks : ", new_partition.B)
            return new_partition, count_log
        end
        current_partition = new_partition
    end
    current_partition, count_log
end
