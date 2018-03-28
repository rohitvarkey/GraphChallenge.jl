type SparseUpdateIBEM
    M::SparseMatrixCSC{Int64, Int64}
    C::Vector{Int64}
    A::SparseMatrixCSC{Int64, Int64}
    D::Array{Int64, 2}
    length::Int64
end

function initialize_edge_counts(
    _::Type{SparseUpdateIBEM}, g::SimpleWeightedDiGraph, B::Int64,
    b::Vector{Int64}, count_log::CountLog
    )
    C = sparse(1:size(b,1), b, 1, size(b, 1), B)
    A = adjacency_matrix(g)
    M = (C' * A * C)'
    D = zeros(Int64, 3, nv(g))
    SparseUpdateIBEM(M, copy(b), A, D, 0)
end

function compute_block_neighbors_and_degrees(
    p::Partition{SparseUpdateIBEM},
    block::Int64,
    count_log::CountLog
    )
    out_neighbors = findn(p.M.M[:, block])
    in_neighbors = findn(p.M.M[block, :])
    neighbors = collect(Set(out_neighbors) ∪ Set(in_neighbors))
    k_in = sum(p.M.M[block, in_neighbors])
    k_out = sum(p.M.M[out_neighbors, block])
    k = k_out + k_in
    # FIXME: Add correct counts
    #count_log.edges_traversed += (size(p.M, 1) + size(p.M, 2))
    neighbors, k_out, k_in, k
end

function compute_block_degrees(
    M::SparseUpdateIBEM, B::Int64, count_log::CountLog
    )
    #TODO: Figure out how to count edges traversed here.
    # Sum across rows to get the outdegrees for each block
    d_out = reshape(sum(M.M, 1), B)
    # Sum across cols to get the indegrees for each block
    d_in = reshape(sum(M.M, 2), B)
    d = d_out + d_in
    # FIXME: Add accurate counts
    #count_log.edges_traversed += (size(M, 1) + size(M, 2))
    return d_out, d_in, d
end

function compute_new_matrix(
    p::Partition{SparseUpdateIBEM}, r::Int64, s::Int64,
    out_block_count_map, in_block_count_map, self_edge_weight::Int64,
    count_log::CountLog
    )
    #TODO: Figure out how to count edges traversed here.
    M_r_row = copy(p.M.M[r, :])
    M_r_col = copy(p.M.M[:, r])
    M_s_row = copy(p.M.M[s, :])
    M_s_col = copy(p.M.M[:, s])

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
    p::Partition{SparseUpdateIBEM}, r::Int64, s::Int64, count_log::CountLog
    )

    #TODO: Figure out how to count edges traversed here.
    M_r_row = sparse(zeros(Int64, p.B))
    M_r_col = sparse(zeros(Int64, p.B))

    M_s_row = copy(p.M.M[s, :])
    M_s_col = copy(p.M.M[:, s])

    #TODO: Optimize this
    M_s_row += p.M.M[r, :]
    M_s_row[r] = 0
    M_s_row[s] += p.M.M[r, r] + p.M.M[s, r]

    M_s_col += p.M.M[:, r]
    M_s_col[r] = 0
    M_s_col[s] += p.M.M[r, r] + p.M.M[r, s]

    #TODO : Check self edge case
    return M_r_row, M_r_col, M_s_row, M_s_col
end

function compute_multinomial_probs(
    p::Partition{SparseUpdateIBEM}, vertex::Int64, count_log::CountLog
    )
    #FIXME: Correct edges traversed counts.
    #count_log.edges_traversed += (size(p.M, 1) + size(p.M, 2))
    return full((p.M.M[:, vertex] .+ p.M.M[vertex, :]) ./ p.d[vertex])
end

function compute_delta_entropy(
    p::Partition{SparseUpdateIBEM}, r::Int64, s::Int64,
    M_r_col::SparseVector{Int64, Int64}, M_s_col::SparseVector{Int64, Int64},
    M_r_row::SparseVector{Int64, Int64}, M_s_row::SparseVector{Int64, Int64},
    d_out_new::Vector{Int64}, d_in_new::Vector{Int64},
    count_log::CountLog
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
        for t1 in findn(p.M.M[:, t2])
            # Skip if t1 is r or s to prevent double counting
            count_log.edges_traversed += 1
            if t1 ∈ (r, s)
                continue
            end
            delta += p.M.M[t1, t2] * log(p.M.M[t1, t2] / p.d_in[t1] / p.d_out[t2])
        end
    end
    # Sum over rows in old M
    for t1 in (r, s)
        for t2 in findn(p.M.M[t1, :])
            count_log.edges_traversed += 1
            delta += p.M.M[t1, t2] * log(p.M.M[t1, t2] / p.d_in[t1] / p.d_out[t2])
        end
    end
    #println("Delta: $delta")
    delta
end

function compute_overall_entropy(
        M::SparseUpdateIBEM, d_out::Vector{Int64}, d_in::Vector{Int64},
        B::Int64, N::Int64, E::Int64, count_log::CountLog
    )
    rows, cols = findn(M.M)  # all non-zero entries
    summation_term = 0.0
    for (col, row) in zip(cols, rows)
        summation_term -= M.M[row, col] * log(M.M[row, col]/ d_in[row] / d_out[col])
    end
    model_S_term = B^2 / E
    model_S = E * (1 + model_S_term) * log(1 + model_S_term) -
        model_S_term * log(model_S_term) + N*log(B)
    S = model_S + summation_term
    count_log.edges_traversed += size(M.M, 1) * size(M.M, 2)
    return S
end

function compute_hastings_correction(
    s::Int64, p::Partition{SparseUpdateIBEM},
    M_r_row::SparseVector{Int64, Int64}, M_r_col::SparseVector{Int64, Int64},
    d_new::Vector{Int64}, blocks_out_count_map, blocks_in_count_map,
    count_log::CountLog
    )
    blocks = Set(keys(blocks_out_count_map)) ∪ Set(keys(blocks_in_count_map))
    p_forward = 0.0
    p_backward = 0.0
    for t in blocks
        degree = get(blocks_out_count_map, t, 0) + get(blocks_in_count_map, t, 0)
        p_forward += degree * (p.M.M[t, s] + p.M.M[s, t] + 1) / (p.d[t] + p.B)
        p_backward += degree * (M_r_row[t] + M_r_col[t] + 1) / (d_new[t] + p.B)
    end
    count_log.edges_traversed += length(blocks)
    return p_backward / p_forward
end

function update_partition(
    p::Partition{SparseUpdateIBEM}, b_new::Vector{Int64},
    g::SimpleWeightedDiGraph, vertex_in_neighbors, count_log::CountLog
    )
    # Create Δ from current_partition.b and b_new.
    p.M.length = 0
    for (vertex, block) in enumerate(b_new)
        if block != p.b[vertex]
            p.M.length += 1
            p.M.D[1, p.M.length] = vertex
            p.M.D[2, p.M.length] = p.b[vertex]
            p.M.D[3, p.M.length] = block
        end
    end
    Δ = sparse(
        cat(1, p.M.D[1, 1:p.M.length], p.M.D[1, 1:p.M.length]),
        cat(1, p.M.D[2, 1:p.M.length], p.M.D[3, 1:p.M.length]),
        cat(1, fill(-1, p.M.length), ones(Int64, p.M.length)),
        size(b_new, 1),
        p.B
    )
    C = sparse(1:size(p.b, 1), p.b, 1, size(p.b, 1), p.B)
    @show size(Δ), size(C), size(p.M.A)
    p.M.M = p.M.M + (Δ' * p.M.A * C + C' * p.M.A * Δ + Δ' * p.M.A * Δ)'
    #p.M.C = copy(b_new) #Is this a necessary copy?
    dropzeros!(p.M.M)
    p.b = b_new
    compute_block_degrees(p, count_log)
    compute_overall_entropy(p, count_log)
end

function zeros_interblock_edge_matrix(::Type{SparseUpdateIBEM}, size::Int64)
    SparseUpdateIBEM(
        sparse(zeros(Int64, size, size)),
        zeros(Int64, size),
        sparse(zeros(Int64, size, size)),
        zeros(Int64, 3, size),
        0
    )
end
