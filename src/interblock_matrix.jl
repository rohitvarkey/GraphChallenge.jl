
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
    k_out = sum(M[:, out_neighbors])
    k_in = sum(M[in_neighbors, :])
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
    M_s_row += M_r_row
    M_s_row[r] = 0
    M_s_row[s] += M[r, r] + M[r, s]

    M_s_col += M_r_col
    M_s_col[r] = 0
    M_s_col[s] += M[r, r] + M[s, r]

    return M_r_row, M_r_col, M_s_row, M_s_col
end

function propose_new_partition(
    M::Array{Int64, 2}, r::Int64, b::Vector{Int64}, B::Int64,
    d::Vector{Int64}, neighbors::Vector{Int64}, agglomerative_move::Bool
    )

    # Pick a neighbor randomly
    if length(neighbors) == 0
        return r
    end
    u = sample(neighbors, Distributions.weights(d[neighbors]./sum(d)))
    if rand() < B / (d[u] + B)
        candidates = Set(1:B)
        if agglomerative_move
            # Force to be different than r.
            pop!(candidates, r)
        end
        s = sample(collect(candidates))
    else
        multinomial_prob = (M[:, u] .+ M[u, :]) ./ d[u]
        if agglomerative_move
            # Force to be different than r.
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
    delta
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
    compute_delta_entropy(
        M, r, s, M_r_col, M_s_col, M_r_row, M_s_col, d_out, d_in,
        new_degrees[1], new_degrees[2]
    )
end
