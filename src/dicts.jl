immutable InterblockEdgeCountDict
    block_in_edges::Dict{Int64, Dict{Int64, Int64}}
    block_out_edges::Dict{Int64, Dict{Int64, Int64}}
end

"""
Initializes the edge count matrix M between the blocks.
Calculates the new interblock_matrix edge count matrix.
"""
function initialize_edge_counts!(
    M::InterblockEdgeCountDict, g::SimpleWeightedDiGraph, B::Int64,
    b::Vector{Int64}
    )
    for edge in edges(g)
        s, d = b[src(edge)], b[dst(edge)]
        if s in keys(M.block_out_edges)
            M.block_out_edges[s][d] = get(M.block_out_edges[s], d, 0) + 1
        else
            M.block_out_edges[s] = Dict(d=>1)
        end
        if d in keys(M.block_in_edges)
            M.block_in_edges[d][s] = get(M.block_in_edges[d], s, 0) + 1
        else
            M.block_in_edges[d] = Dict(s=>1)
        end
    end
    M
end

function initialize_edge_counts(
    _::Type{InterblockEdgeCountDict}, g::SimpleWeightedDiGraph, B::Int64,
    b::Vector{Int64}
    )
    M = InterblockEdgeCountDict(
        Dict{Int64, Dict{Int64, Int64}}(),
        Dict{Int64, Dict{Int64, Int64}}()
    )
    initialize_edge_counts!(M, g, B, b)
    M
end

function compute_block_neighbors_and_degrees(M::InterblockEdgeCountDict, block::Int64)
    out_neighbors = collect(keys(M.block_out_edges[block]))
    in_neighbors = collect(keys(M.block_in_edges[block]))
    neighbors = collect(Set(out_neighbors) ∪ Set(in_neighbors))
    k_in = sum(values(M.block_in_edges[block]))
    k_out = sum(values(M.block_out_edges[block]))
    k = k_out + k_in
    neighbors, k_out, k_in, k
end


function compute_block_degrees(M::InterblockEdgeCountDict, B::Int64)
    # Sum across rows to get the outdegrees for each block
    d_out = zeros(Int64, B)
    d_in = zeros(Int64, B)
    for (k, v) in M.block_out_edges
        d_out[k] = sum(values(v))
    end
    for (k, v) in M.block_in_edges
        d_in[k] = sum(values(v))
    end
    d = d_out + d_in
    return d_out, d_in, d
end

function compute_new_matrix(
    M::InterblockEdgeCountDict, r::Int64, s::Int64, num_blocks::Int64,
    out_block_count_map, in_block_count_map, self_edge_weight::Int64
    )

    if r in keys(M.block_out_edges)
        M_r_col = copy(M.block_out_edges[r])
    else
        M_r_col = Dict{Int64, Int64}()
    end
    if r in keys(M.block_in_edges)
        M_r_row = copy(M.block_in_edges[r])
    else
        M_r_row = Dict{Int64, Int64}()
    end
    if s in keys(M.block_out_edges)
        M_s_col = copy(M.block_out_edges[s])
    else
        M_s_col = Dict{Int64, Int64}()
    end
    if s in keys(M.block_in_edges)
        M_s_row = copy(M.block_in_edges[s])
    else
        M_s_row = Dict{Int64, Int64}()
    end

    @show r, s, out_block_count_map, in_block_count_map
    @show M_r_row
    @show M_r_col
    @show M_s_row
    @show M_s_col

    for (block, out_count) in out_block_count_map
        # Move outgoing edges from r to s.
        M_r_col[block] -= out_count
        M_s_col[block] = get(M_s_col, block, 0) + out_count
        if (block => 0) in M_r_col
            pop!(M_r_col, block)
        end
        if block == r
            # Edges in the same block.
            M_r_row[r] -= out_count # Remove from in count of r.
            # Add to the in count of edges to r from s
            M_r_row[s] = get(M_r_row, s, 0) + out_count
            if (r => 0) in M_r_row
                pop!(M_r_row, r)
            end
        elseif block == s
            # Edges from r to s.
            M_s_row[r] -= out_count
            # Add as self edges
            M_s_row[s] = get(M_s_row, s, 0) + out_count
            if (r => 0) in M_s_row
                pop!(M_s_row, r)
            end
        end
    end

    for (block, in_count) in in_block_count_map
        M_r_row[block] -= in_count
        M_s_row[block] = get(M_s_row, block, 0) + in_count
        if (block => 0) in M_r_row
            pop!(M_r_row, block)
        end
        if block == r
            M_r_col[r] -= in_count
            M_r_col[s] = get(M_r_col, s, 0) + in_count
            if (r => 0) in M_r_col
                pop!(M_r_col, r)
            end
        elseif block == s
            M_s_col[r] -= in_count
            M_s_col[s] = get(M_s_col, s, 0) + in_count
            if (s => 0) in M_s_col
                pop!(M_s_col, s)
            end
        end
    end

    if self_edge_weight > 0
        # Correct counts based on self_edge_weight
        M_s_row[r] -= self_edge_weight
        M_s_row[s] = get(M_s_row, s, 0) + self_edge_weight
        M_s_col[r] -= self_edge_weight
        M_s_col[s] = get(M_s_col, s, 0) + self_edge_weight
    end

    if (r => 0) in M_s_row
        pop!(M_s_row, r)
    end
    if (r => 0) in M_s_col
        pop!(M_s_col, r)
    end

    @show M_r_row
    @show M_r_col
    @show M_s_row
    @show M_s_col

    return M_r_row, M_r_col, M_s_row, M_s_col
end

"""Computes the new rows and cols in `M`, when all nodes from `r` are shifted to
block `s`."""
function compute_new_matrix_agglomerative(
    M::InterblockEdgeCountDict, r::Int64, s::Int64, num_blocks::Int64
    )

    M_r_row = Dict{Int64, Int64}()
    M_r_col = Dict{Int64, Int64}()

    if s in keys(M.block_in_edges)
        M_s_row = copy(M.block_in_edges[s])
    else
        M_s_row = Dict{Int64, Int64}()
    end
    if s in keys(M.block_out_edges)
        M_s_col = copy(M.block_out_edges[s])
    else
        M_s_col = Dict{Int64, Int64}()
    end

    if r in keys(M.block_out_edges)
        # Add all outgoing edges in r to s.
        for (out_neighbor, edgecount) in M.block_out_edges[r]
            M_s_col[out_neighbor] = get(M_s_col, out_neighbor, 0) + edgecount
        end
        # Add self edges within r to s
        self_edge_counts = get(M_s_col, s, 0) + get(M.block_out_edges[r], r, 0)
        if self_edge_counts > 0
            M_s_col[s] = self_edge_counts
        end
    end

    if r in keys(M_s_col)
        pop!(M_s_col, r) #Set to 0 by popping.
    end

    # Add edges that went from s to r
    if s in keys(M.block_out_edges)
        self_edge_counts  = get(M_s_col, s, 0) + get(M.block_out_edges[s], r, 0)
        if self_edge_counts > 0
            M_s_col[s] = self_edge_counts
        end
    end

    if r in keys(M.block_in_edges)
        # Add all incoming edges in r to s
        for (in_neighbor, edgecount) in M.block_in_edges[r]
            M_s_row[in_neighbor] = get(M_s_row, in_neighbor, 0) + edgecount
        end
        # Add self edges within r to s
        self_edge_counts = get(M_s_row, s, 0) + get(M.block_in_edges[r], r, 0)
        if self_edge_counts > 0
            M_s_row[s] = self_edge_counts
        end

    end
    if r in keys(M_s_row)
        pop!(M_s_row, r) #Set to 0 by popping.
    end

    # Add all edges that went from r to s.
    if s in keys(M.block_in_edges)
        self_edge_counts = get(M_s_row, s, 0) + get(M.block_in_edges[s], r, 0)
        if self_edge_counts > 0
            M_s_row[s] = self_edge_counts
        end
    end

    return M_r_row, M_r_col, M_s_row, M_s_col
end


function compute_multinomial_probs(
    M::InterblockEdgeCountDict,  degrees::Vector{Int64}, block::Int64
    )
    probabilities = zeros(length(degrees))
    for (out_neighbor, edgecount) in M.block_out_edges[block]
        probabilities[out_neighbor] += edgecount/degrees[block]
    end
    for (in_neighbor, edgecount) in M.block_in_edges[block]
        probabilities[in_neighbor] += edgecount/degrees[block]
    end
    return probabilities
end

function compute_delta_entropy(
    M::InterblockEdgeCountDict, r::Int64, s::Int64,
    M_r_col::Dict{Int64, Int64}, M_s_col::Dict{Int64, Int64},
    M_r_row::Dict{Int64, Int64}, M_s_row::Dict{Int64, Int64},
    d_out::Vector{Int64}, d_in::Vector{Int64},
    d_out_new::Vector{Int64}, d_in_new::Vector{Int64}
    )
    delta = 0.0
    # Sum over col of r in new M
    for (t1, edgecount) in M_r_col
        # Skip if t1 is r or s to prevent double counting
        if t1 ∈ (r, s)
            continue
        end
        delta -= edgecount * log(edgecount / d_in_new[t1] / d_out_new[r])
    end
    for (t1, edgecount) in M_s_col
        if t1 ∈ (r, s)
            continue
        end
        delta -= edgecount * log(edgecount / d_in_new[t1] / d_out_new[s])
    end
    # Sum over row of r in new M
    for (t2, edgecount) in M_r_row
        delta -= edgecount * log(edgecount / d_in_new[r] / d_out_new[t2])
    end
    # Sum over row of s in new M
    for (t2, edgecount) in M_s_row
        delta -= edgecount * log(edgecount / d_in_new[s] / d_out_new[t2])
    end
    # Sum over columns in old M
    for t2 in (r, s)
        for (t1, edgecount) in M.block_out_edges[t2]
            # Skip if t1 is r or s to prevent double counting
            if t1 ∈ (r, s)
                continue
            end
            delta += edgecount * log(edgecount / d_in[t1] / d_out[t2])
        end
    end
    # Sum over rows in old M
    for t1 in (r, s)
        for (t2, edgecount) in M.block_in_edges[t1]
            delta += edgecount * log(edgecount / d_in[t1] / d_out[t2])
        end
    end
    #println("Delta: $delta")
    delta
end

function compute_overall_entropy(
        M::InterblockEdgeCountDict, d_out::Vector{Int64}, d_in::Vector{Int64},
        B::Int64, N::Int64, E::Int64
    )
    summation_term = 0.0
    for (out_block, edges) in M.block_out_edges
        for (in_block, edgecount) in edges
            summation_term -= edgecount * log(edgecount/ d_in[in_block] / d_out[out_block])
        end
    end
    model_S_term = B^2 / E
    model_S = E * (1 + model_S_term) * log(1 + model_S_term) -
        model_S_term * log(model_S_term) + N*log(B)
    S = model_S + summation_term
    return S
end

function compute_hastings_correction(
        s::Int64, M::InterblockEdgeCountDict, M_r_row::Dict{Int64, Int64},
        M_r_col::Dict{Int64, Int64}, B::Int64, d::Vector{Int64}, d_new::Vector{Int64},
        blocks_out_count_map, blocks_in_count_map
    )
    blocks = Set(keys(blocks_out_count_map)) ∪ Set(keys(blocks_in_count_map))
    p_forward = 0.0
    p_backward = 0.0
    for t in blocks
        degree = get(blocks_out_count_map, t, 0) + get(blocks_in_count_map, t, 0)
        if t in keys(M.block_out_edges)
            m = get(M.block_out_edges[t], s, 0)
        end
        if s in keys(M.block_out_edges)
            m += get(M.block_out_edges[s], t, 0)
        end
        p_forward += degree * (m + 1) / (d[t] + B)
        p_backward += degree * (get(M_r_row, t, 0) + get(M_r_col, t, 0) + 1) / (d_new[t] + B)
    end
    return p_backward / p_forward
end

function evaluate_proposal_agg(
    M::InterblockEdgeCountDict, r::Int64, s::Int64, num_blocks::Int64,
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
        M, r, s, M_r_col, M_s_col, M_r_row, M_s_row, d_out, d_in,
        new_degrees[1], new_degrees[2]
    )
end

function evaluate_nodal_proposal(
    M::InterblockEdgeCountDict, r::Int64, s::Int64, num_blocks::Int64, β::Int64,
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
    M::InterblockEdgeCountDict, r::Int64, s::Int64,
    M_r_col::Dict{Int64, Int64}, M_s_col::Dict{Int64, Int64},
    M_r_row::Dict{Int64, Int64}, M_s_row::Dict{Int64, Int64}
    )
    M.block_out_edges[r] = M_r_col
    M.block_in_edges[r] = M_r_row
    M.block_out_edges[s] = M_s_col
    M.block_in_edges[s] = M_s_row
    for (out_block, edgecount) in M_r_col
        M.block_in_edges[out_block][r] = edgecount
    end
    for (out_block, edgecount) in M_s_col
        M.block_in_edges[out_block][s] = edgecount
    end
    for (in_block, edgecount) in M_r_row
        M.block_out_edges[in_block][r] = edgecount
    end
    for (in_block, edgecount) in M_s_row
        M.block_out_edges[in_block][s] = edgecount
    end
    println("Updated partition")
    @show r, s
    @show M.block_out_edges[r], M.block_in_edges[r]
    @show M.block_out_edges[s], M.block_in_edges[s]
    M
end
