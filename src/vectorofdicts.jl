immutable InterblockEdgeCountVectorDict
    block_in_edges::Vector{Dict{Int64, Int64}}
    block_out_edges::Vector{Dict{Int64, Int64}}
end

"""
Initializes the edge count matrix M between the blocks.
Calculates the new interblock_matrix edge count matrix.
"""
function initialize_edge_counts!(
    M::InterblockEdgeCountVectorDict, g::SimpleWeightedDiGraph, B::Int64,
    b::Vector{Int64}
    )
    for edge in edges(g)
        s, d = b[src(edge)], b[dst(edge)]
        M.block_out_edges[s][d] = get(M.block_out_edges[s], d, 0) + 1
        M.block_in_edges[d][s] = get(M.block_in_edges[d], s, 0) + 1
    end
    M
end

function initialize_edge_counts(
    _::Type{InterblockEdgeCountVectorDict}, g::SimpleWeightedDiGraph, B::Int64,
    b::Vector{Int64}
    )
    M = InterblockEdgeCountVectorDict(
        [Dict{Int64, Int64}() for i=1:B],
        [Dict{Int64, Int64}() for i=1:B]
    )
    initialize_edge_counts!(M, g, B, b)
    M
end

function compute_block_neighbors_and_degrees(M::InterblockEdgeCountVectorDict, block::Int64)
    neighbors = Set{Int64}()
    k_in = zero(Int64)
    k_out = zero(Int64)
    for (neighbor, out_count) in M.block_out_edges[block]
        k_out += out_count
        push!(neighbors, neighbor)
    end
    for (neighbor, in_count) in M.block_in_edges[block]
        k_in += in_count
        push!(neighbors, neighbor)
    end
    k = k_out + k_in
    collect(neighbors), k_out, k_in, k
end


function compute_block_degrees(M::InterblockEdgeCountVectorDict, B::Int64)
    # Sum across rows to get the outdegrees for each block
    d_out = [sum(values(out_blocks)) for out_blocks in M.block_out_edges]
    d_in = [sum(values(in_blocks)) for in_blocks in M.block_in_edges]
    d = d_out + d_in
    return d_out, d_in, d
end

function compute_new_matrix(
    M::InterblockEdgeCountVectorDict, r::Int64, s::Int64, num_blocks::Int64,
    out_block_count_map, in_block_count_map, self_edge_weight::Int64
    )

    M_r_col = copy(M.block_out_edges[r])
    M_r_row = copy(M.block_in_edges[r])
    M_s_col = copy(M.block_out_edges[s])
    M_s_row = copy(M.block_in_edges[s])

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

    return M_r_row, M_r_col, M_s_row, M_s_col
end

"""Computes the new rows and cols in `M`, when all nodes from `r` are shifted to
block `s`."""
function compute_new_matrix_agglomerative(
    M::InterblockEdgeCountVectorDict, r::Int64, s::Int64, num_blocks::Int64
    )

    M_r_row = Dict{Int64, Int64}()
    M_r_col = Dict{Int64, Int64}()
    M_s_row = copy(M.block_in_edges[s])
    M_s_col = copy(M.block_out_edges[s])

    for (out_neighbor, edgecount) in M.block_out_edges[r]
        M_s_col[out_neighbor] = get(M_s_col, out_neighbor, 0) + edgecount
    end
    # Add self edges within r to s
    self_edge_counts = get(M_s_col, s, 0) + get(M.block_out_edges[r], r, 0)
    if self_edge_counts > 0
        M_s_col[s] = self_edge_counts
    end

    if r in keys(M_s_col)
        pop!(M_s_col, r) #Set to 0 by popping.
    end

    # Add edges that went from s to r
    self_edge_counts  = get(M_s_col, s, 0) + get(M.block_out_edges[s], r, 0)
    if self_edge_counts > 0
        M_s_col[s] = self_edge_counts
    end

    # Add all incoming edges in r to s
    for (in_neighbor, edgecount) in M.block_in_edges[r]
        M_s_row[in_neighbor] = get(M_s_row, in_neighbor, 0) + edgecount
    end
    # Add self edges within r to s
    self_edge_counts = get(M_s_row, s, 0) + get(M.block_in_edges[r], r, 0)
    if self_edge_counts > 0
        M_s_row[s] = self_edge_counts
    end

    if r in keys(M_s_row)
        pop!(M_s_row, r) #Set to 0 by popping.
    end

    # Add all edges that went from r to s.
    self_edge_counts = get(M_s_row, s, 0) + get(M.block_in_edges[s], r, 0)
    if self_edge_counts > 0
        M_s_row[s] = self_edge_counts
    end

    return M_r_row, M_r_col, M_s_row, M_s_col
end


function compute_multinomial_probs(
    M::InterblockEdgeCountVectorDict,  degrees::Vector{Int64}, block::Int64
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
    M::InterblockEdgeCountVectorDict, r::Int64, s::Int64,
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
        for (t1, edgecount) in get(M.block_out_edges, t2, Dict{Int64, Int64}())
            # Skip if t1 is r or s to prevent double counting
            if t1 ∈ (r, s)
                continue
            end
            delta += edgecount * log(edgecount / d_in[t1] / d_out[t2])
        end
    end
    # Sum over rows in old M
    for t1 in (r, s)
        for (t2, edgecount) in get(M.block_in_edges, t1, Dict{Int64, Int64}())
            delta += edgecount * log(edgecount / d_in[t1] / d_out[t2])
        end
    end
    #println("Delta: $delta")
    delta
end

function compute_overall_entropy(
        M::InterblockEdgeCountVectorDict, d_out::Vector{Int64}, d_in::Vector{Int64},
        B::Int64, N::Int64, E::Int64
    )
    summation_term = 0.0
    for (out_block, edges) in enumerate(M.block_out_edges)
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
        s::Int64, M::InterblockEdgeCountVectorDict, M_r_row::Dict{Int64, Int64},
        M_r_col::Dict{Int64, Int64}, B::Int64, d::Vector{Int64}, d_new::Vector{Int64},
        blocks_out_count_map, blocks_in_count_map
    )
    blocks = Set(keys(blocks_out_count_map)) ∪ Set(keys(blocks_in_count_map))
    p_forward = 0.0
    p_backward = 0.0
    for t in blocks
        degree = get(blocks_out_count_map, t, 0) + get(blocks_in_count_map, t, 0)
        m = get(M.block_out_edges[t], s, 0)
        m += get(M.block_out_edges[s], t, 0)
        p_forward += degree * (m + 1) / (d[t] + B)
        p_backward += degree * (get(M_r_row, t, 0) + get(M_r_col, t, 0) + 1) / (d_new[t] + B)
    end
    return p_backward / p_forward
end

function update_partition(
    M::InterblockEdgeCountVectorDict, r::Int64, s::Int64,
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
    M
end

function zeros_interblock_edge_matrix(::Type{InterblockEdgeCountVectorDict}, size::Int64)
    return InterblockEdgeCountVectorDict(
    [Dict{Int64, Int64}() for i=1:size],
    [Dict{Int64, Int64}() for i=1:size]
    )
end
