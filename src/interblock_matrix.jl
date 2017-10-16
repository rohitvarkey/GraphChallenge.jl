
"""
Initializes the edge count matrix M between the blocks.
Calculates the new out, in and total degrees for the updated edge count matrix.
Returns a tuple of M, d_out, d_in, d
"""
function initialize_edge_counts!(
    M::Array{Int64, 2}, g::SimpleWeightedDiGraph, B::Int64, b::Vector{Int64}
    )
    @assert size(M) == (B, B)
    for edge in edges(g)
            M[b[dst(edge)], b[src(edge)]] += weight(edge)
    end
end

function compute_block_neighbors_and_degrees(M::Array{Int64, 2}, block::Int64)
    out_neighbors = findn(M[:, block])
    in_neighbors = findn(M[block, :])
    neighbors = collect(Set(out_neighbors) ∪ Set(in_neighbors))
    k_in = sum(M[block, in_neighbors])
    k_out = sum(M[out_neighbors, block])
    k = k_out + k_in
    neighbors, k_out, k_in, k
end
function compute_block_degrees(M::Array{Int64, 2}, B::Int64)
    # Sum across rows to get the outdegrees for each block
    d_out = reshape(sum(M, 1), B)
    # Sum across cols to get the indegrees for each block
    d_in = reshape(sum(M, 2), B)
    d = d_out + d_in
    return d_out, d_in, d
end

function initialize_edge_counts(
    _::Type{Array{Int64, 2}}, g::SimpleWeightedDiGraph, B::Int64,
    b::Vector{Int64}
    )
    M = zeros(Int64, B, B)
    initialize_edge_counts!(M, g, B, b)
    M
end

function compute_new_matrix(
    M::Array{Int64, 2}, r::Int64, s::Int64, num_blocks::Int64,
    out_block_count_map, in_block_count_map, self_edge_weight::Int64
    )

    M_r_row = copy(M[r, :])
    M_r_col = copy(M[:, r])
    M_s_row = copy(M[s, :])
    M_s_col = copy(M[:, s])

    for (block, out_count) in out_block_count_map
        M_r_col[block] -= out_count
        M_s_col[block] += out_count
        if block == r
            M_r_row[r] -= out_count
            M_r_row[s] += out_count
        elseif block == s
            M_s_row[r] -= out_count
            M_s_row[s] += out_count
        end
    end

    for (block, in_count) in in_block_count_map
        M_r_row[block] -= in_count
        M_s_row[block] += in_count
        if block == r
            M_r_col[r] -= in_count
            M_r_col[s] += in_count
        elseif block == s
            M_s_col[r] -= in_count
            M_s_col[s] += in_count
        end
    end

    M_s_row[r] -= self_edge_weight
    M_s_row[s] += self_edge_weight
    M_s_col[r] -= self_edge_weight
    M_s_col[s] += self_edge_weight

    return M_r_row, M_r_col, M_s_row, M_s_col
end

"""Computes the new rows and cols in `M`, when all nodes from `r` are shifted to
block `s`."""
function compute_new_matrix_agglomerative(
    M::Array{Int64, 2}, r::Int64, s::Int64, num_blocks::Int64
    )

    M_r_row = zeros(Int64, num_blocks)
    M_r_col = zeros(Int64, num_blocks)

    M_s_row = copy(M[s, :])
    M_s_col = copy(M[:, s])

    #TODO: Optimize this
    M_s_row += M[r, :]
    M_s_row[r] = 0
    M_s_row[s] += M[r, r] + M[s, r]

    M_s_col += M[:, r]
    M_s_col[r] = 0
    M_s_col[s] += M[r, r] + M[r, s]

    #TODO : Check self edge case
    return M_r_row, M_r_col, M_s_row, M_s_col
end

function propose_new_partition_agg(
    M::Array{Int64, 2}, r::Int64, b::Vector{Int64}, B::Int64,
    d::Vector{Int64}, neighbors::Vector{Int64}
    )

    # Pick a neighbor randomly
    if length(neighbors) == 0
        candidates = Set(1:B)
        # Force to be different than r.
        pop!(candidates, r)
        s = sample(collect(candidates))
        return s
    end
    rand_neighbor = sample(neighbors, Distributions.weights(d[neighbors]./sum(d)))
    u = b[rand_neighbor]
    if rand() < B / (d[u] + B)
        candidates = Set(1:B)
        pop!(candidates, r)
        s = sample(collect(candidates))
    else
        multinomial_prob = (M[:, u] .+ M[u, :]) ./ d[u]
        multinomial_prob[r] = 0
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

function propose_new_partition_nodal(
    M::Array{Int64, 2}, r::Int64, b::Vector{Int64},
    B::Int64, d::Vector{Int64}, neighbors::Vector{Int64}, wts::Vector{Int64},
    )

    # Pick a neighbor randomly
    if length(neighbors) == 0
        candidates = Set(1:B)
        s = sample(collect(candidates))
        return s
    end
    rand_neighbor = sample(neighbors, Distributions.weights(wts./sum(wts)))
    u = b[rand_neighbor]
    if rand() < B / (d[u] + B)
        candidates = Set(1:B)
        pop!(candidates, r)
        s = sample(collect(candidates))
    else
        multinomial_prob = (M[:, u] .+ M[u, :]) ./ d[u]
        multinomial_prob[r] = 0
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

function compute_delta_entropy(
    M::Array{Int64, 2}, r::Int64, s::Int64,
    M_r_col::Vector{Int64}, M_s_col::Vector{Int64}, M_r_row::Vector{Int64},
    M_s_row::Vector{Int64}, d_out::Vector{Int64}, d_in::Vector{Int64},
    d_out_new::Vector{Int64}, d_in_new::Vector{Int64}
    )
    delta = 0.0
    # Sum over col of r in new M
    for t1 in findn(M_r_col)
        # Skip if t1 is r or s to prevent double counting
        if t1 ∈ (r, s)
            continue
        end
        delta -= M_r_col[t1] * log(M_r_col[t1] / d_in_new[t1] / d_out_new[r])
    end
    for t1 in findn(M_s_col)
        if t1 ∈ (r, s)
            continue
        end
        delta -= M_s_col[t1] * log(M_s_col[t1] / d_in_new[t1] / d_out_new[s])
    end
    # Sum over row of r in new M
    for t2 in findn(M_r_row)
        delta -= M_r_row[t2] * log(M_r_row[t2] / d_in_new[r] / d_out_new[t2])
    end
    # Sum over row of s in new M
    for t2 in findn(M_s_row)
        delta -= M_s_row[t2] * log(M_s_row[t2] / d_in_new[s] / d_out_new[t2])
    end
    # Sum over columns in old M
    for t2 in (r, s)
        for t1 in findn(M[:, t2])
            # Skip if t1 is r or s to prevent double counting
            if t1 ∈ (r, s)
                continue
            end
            delta += M[t1, t2] * log(M[t1, t2] / d_in[t1] / d_out[t2])
        end
    end
    # Sum over rows in old M
    for t1 in (r, s)
        for t2 in findn(M[t1, :])
            delta += M[t1, t2] * log(M[t1, t2] / d_in[t1] / d_out[t2])
        end
    end
    #println("Delta: $delta")
    delta
end

function compute_overall_entropy(
        M::Array{Int64, 2}, d_out::Vector{Int64}, d_in::Vector{Int64},
        B::Int64, N::Int64, E::Int64
    )
    rows, cols = findn(M)  # all non-zero entries
    summation_term = 0.0
    for (col, row) in zip(cols, rows)
        summation_term -= M[row, col] * log(M[row, col]/ d_in[row] / d_out[col])
    end
    model_S_term = B^2 / E
    model_S = E * (1 + model_S_term) * log(1 + model_S_term) -
        model_S_term * log(model_S_term) + N*log(B)
    S = model_S + summation_term
    return S
end

function compute_hastings_correction(
        s::Int64, M::Array{Int64, 2}, M_r_row::Vector{Int64},
        M_r_col::Vector{Int64}, B::Int64, d::Vector{Int64}, d_new::Vector{Int64},
        blocks_out_count_map, blocks_in_count_map
    )
    blocks = Set(keys(blocks_out_count_map)) ∪ Set(keys(blocks_in_count_map))
    p_forward = 0.0
    p_backward = 0.0
    for t in blocks
        degree = get(blocks_out_count_map, t, 0) + get(blocks_in_count_map, t, 0)
        p_forward += degree * (M[t, s] + M[s, t] + 1) / (d[t] + B)
        p_backward += degree * (M_r_row[t] + M_r_col[t] + 1) / (d_new[t] + B)
    end
    return p_backward / p_forward
end

function evaluate_proposal_agg(
    M::Array{Int64, 2}, r::Int64, s::Int64, num_blocks::Int64,
    d::Vector{Int64}, d_in::Vector{Int64}, d_out::Vector{Int64},
    k::Int64, k_in::Int64, k_out::Int64
    )
    M_r_row, M_r_col, M_s_row, M_s_col =
        compute_new_matrix_agglomerative(M, r, s, num_blocks)
    new_degrees = [copy(degrees) for degrees in [d_out, d_in, d]]
    for (new_d, degree) in zip(new_degrees, [k_out, k_in, k])
        new_d[r] -= degree
        new_d[s] += degree
    end
    #@show d_in, new_degrees[2], k_in
    #@show d_out, new_degrees[1], k_out
    #@show d, new_degrees[3], k
    compute_delta_entropy(
        M, r, s, M_r_col, M_s_col, M_r_row, M_s_row, d_out, d_in,
        new_degrees[1], new_degrees[2]
    )
end

function evaluate_nodal_proposal(
    M::Array{Int64, 2}, r::Int64, s::Int64, num_blocks::Int64, β::Int64,
    d::Vector{Int64}, d_in::Vector{Int64}, d_out::Vector{Int64},
    k::Int64, k_in::Int64, k_out::Int64, self_edge_weight::Int64,
    blocks_out_count_map, blocks_in_count_map
    )

    M_r_row, M_r_col, M_s_row, M_s_col = compute_new_matrix(
        M, r, s, num_blocks,
        blocks_out_count_map, blocks_in_count_map,
        self_edge_weight
    )

    new_degrees = [copy(degrees) for degrees in [d_out, d_in, d]]
    for (new_d, degree) in zip(new_degrees, [k_out, k_in, k])
        new_d[r] -= degree
        new_d[s] += degree
    end

    hastings_correction = compute_hastings_correction(
        s, M, M_r_row, M_r_col, num_blocks, d, new_degrees[3],
        blocks_out_count_map, blocks_in_count_map
    )

    Δ = compute_delta_entropy(
            M, r, s, M_r_col, M_s_col, M_r_row, M_s_row, d_out, d_in,
            new_degrees[1], new_degrees[2]
        )

    p_accept = min(exp(-β * Δ) * hastings_correction, 1)

    #println("p_accept: $(p_accept), Δ: $Δ, β: $β, H: $(hastings_correction), exp: $(exp(-β*Δ)*hastings_correction)")

    M_r_row, M_r_col, M_s_row, M_s_col, Δ, p_accept
end

function update_partition(
    M::Array{Int64, 2}, r::Int64, s::Int64,
    M_r_col::Vector{Int64}, M_s_col::Vector{Int64},
    M_r_row::Vector{Int64}, M_s_row::Vector{Int64}
    )
    M[:, r] = M_r_col
    M[r, :] = M_r_row
    M[:, s] = M_s_col
    M[s, :] = M_s_row
    #info("Updated partition")
    M
end

function prepare_for_partition_on_next_num_blocks(
    current_partition::Partition{Array{Int64, 2}},
    best_partitions::Vector{Partition{Array{Int64, 2}}},
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
        partition = copy(best_partitions[2])
    else
        if best_partitions[1].B - best_partitions[3].B == 2
            optimal_B_found = true
            partition = copy(best_partitions[2])
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
            partition = copy(best_partitions[index])
        end
    end
    return partition, best_partitions, optimal_B_found, B_to_merge
end
