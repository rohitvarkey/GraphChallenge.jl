import StingerGraphs: Stinger, insert_edge!, indegree, outdegree, foralledges, edgeparse, stingerconfig, edgeweight

"""
Initializes the edge count matrix M between the blocks.
Calculates the new interblock_matrix edge count matrix.
"""
function initialize_edge_counts!(
    M::Stinger, g::SimpleWeightedDiGraph, B::Int64,
    b::Vector{Int64}
    )
    #TODO: Optimize this by using batch insertions.
    for edge in edges(g)
        s, d = b[src(edge)], b[dst(edge)]
        insert_edge!(M, 0, s, d, 1, 1)
    end
    M
end

function initialize_edge_counts(
    _::Type{Stinger}, g::SimpleWeightedDiGraph, B::Int64,
    b::Vector{Int64}
    )
    M = Stinger(stingerconfig(B + 1)) # To prevent off by one errors
    initialize_edge_counts!(M, g, B, b)
    M
end

function compute_block_neighbors_and_degrees(M::Stinger, block::Int64)
    neighbors = Set{Int64}()
    k_in = zero(Int64)
    k_out = zero(Int64)
    foralledges(M, block) do edge, src, etype
        direction, neighbor = edgeparse(edge)
        push!(neighbors, neighbor)
        if direction != 1
            k_out += edge.weight
        end
        if direction != 2
            k_in += edgeweight(M, neighbor, block, 0)
        end
    end
    collect(neighbors), k_out, k_in, k_in + k_out
end


function compute_block_degrees(M::Stinger, B::Int64)
    # Sum across rows to get the outdegrees for each block
    d_out = [outdegree(M, block) for block=1:B]
    d_in = [indegree(M, block) for block=1:B]
    d = d_out + d_in
    return d_out, d_in, d
end

function compute_new_matrix(
    M::Stinger, r::Int64, s::Int64, num_blocks::Int64,
    out_block_count_map, in_block_count_map, self_edge_weight::Int64
    )

    M_r_col = zeros(Int64, num_blocks)
    M_r_row = zeros(Int64, num_blocks)
    M_s_col = zeros(Int64, num_blocks)
    M_s_row = zeros(Int64, num_blocks)

    foralledges(M, r) do edge, src, etype
        direction, block = edgeparse(edge)
        # Move outgoing edges from r to s.
        if direction != 1
            # out edge
            out_count = get(out_block_count_map, block, 0)
            M_r_col[block] = edge.weight - out_count
            M_s_col[block] = out_count
            if block == r
                # Edges in the same block.
                M_r_row[r] -= out_count # Remove from in count of r.
                # Add to the in count of edges to r from s
                M_r_row[s] += out_count
            elseif block == s
                # Edges from r to s.
                M_s_row[r] -= out_count
                # Add as self edges
                M_s_row[s] += out_count
            end
        end
        if direction != 2
            # in edge
            in_count = get(in_block_count_map, block, 0)
            M_r_row[block] = edge.weight - in_count
            M_s_row[block] = in_count
            if block == r
                M_r_col[r] -= in_count
                M_r_col[s] += in_count
            elseif block == s
                M_s_col[r] -= in_count
                M_s_col[s] += in_count
            end
        end
    end

    foralledges(M, s) do edge, src, etype
        direction, block = edgeparse(edge)
        if direction != 1
            M_s_col[block] += edge.weight
        end
        if direction != 2
            M_s_row[block] += edgeweight(M, block, s, 0)
        end
    end


    if self_edge_weight > 0
        # Correct counts based on self_edge_weight
        M_s_row[r] -= self_edge_weight
        M_s_row[s] += self_edge_weight
        M_s_col[r] -= self_edge_weight
        M_s_col[s] += self_edge_weight
    end

    return M_r_row, M_r_col, M_s_row, M_s_col
end

"""Computes the new rows and cols in `M`, when all nodes from `r` are shifted to
block `s`."""
function compute_new_matrix_agglomerative(
    M::Stinger, r::Int64, s::Int64, num_blocks::Int64
    )

    M_r_col = zeros(Int64, num_blocks)
    M_r_row = zeros(Int64, num_blocks)
    M_s_col = zeros(Int64, num_blocks)
    M_s_row = zeros(Int64, num_blocks)

    foralledges(M, s) do edge, src, etype
        direction, block = edgeparse(edge)
        if direction != 1
            # out edges
            M_s_col[block] += edge.weight
            if block == r
                # s->r edge
                M_s_row[s] += edge.weight
                M_s_col[s] += edge.weight
            end
        end
        if direction != 2
            M_s_row[block] += edgeweight(M, block, s, 0)
        end
    end

    foralledges(M, r) do edge, src, etype
        direction, block = edgeparse(edge)
        # Move outgoing edges from r to s.
        if direction != 1
            # out edge
            M_s_col[block] += edge.weight
            if block == r || block == s
                # r->r self edge and r->s edge
                M_s_row[s] += edge.weight
                if block == r
                    M_s_col[s] += edge.weight
                end
            end
        end
        # Avoid double counting s->r.
        if direction != 2 && block != s
            # in edge
            M_s_row[block] += edgeweight(M, block, r, 0)
        end
    end

    M_s_row[r] = 0
    M_s_col[r] = 0

    return M_r_row, M_r_col, M_s_row, M_s_col
end


function compute_multinomial_probs(
    M::Stinger, degrees::Vector{Int64}, block::Int64
    )
    probabilities = zeros(length(degrees))
    foralledges(M, block) do edge, src, etype
        direction, block = edgeparse(edge)
        if direction != 1
            # out edge
            probabilities[neighbor] += edge.weight
        end
        if direction != 2
            # in edge
            probabilities[neighbor] += edgeweight(M, neighbor, block, 0)
        end
    end
    return probabilities
end

function compute_delta_entropy(
    M::Stinger, r::Int64, s::Int64,
    M_r_col::Array{Int64, 1}, M_s_col::Array{Int64, 1},
    M_r_row::Array{Int64, 1}, M_s_row::Array{Int64, 1},
    d_out::Vector{Int64}, d_in::Vector{Int64},
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
    # Sum over edges in old M
    for block in (r, s)
        foralledges(M, block) do edge, src, etype
            direction, neighbor = edgeparse(edge)
            if direction != 1
                # edge is block -> neighbor
                delta += edge.weight * log(edge.weight / d_in[neighbor] / d_out[block])
            end
            # Prevent double counting the r->s, s->r edges.
            if direction != 2 && !((block, neighbor) == (r, s) || (block, neighbor) == (s, r))
                # edge is neighbor -> block
                edgecount = edgeweight(M, neighbor, block, 0)
                delta += edgecount * log(edgecount / d_in[block] / d_out[neighbor])
            end
        end
    end
    #println("Delta: $delta")
    delta
end

function compute_overall_entropy(
        M::Stinger, d_out::Vector{Int64}, d_in::Vector{Int64},
        B::Int64, N::Int64, E::Int64
    )
    summation_term = 0.0
    for block=1:B
        foralledges(M, B) do edge, src, etype
            direction, neighbor = edgeparse(edge)
            if direction != 1
                edgecount = edge.weight
                summation_term -= edgecount * log(edgecount/ d_in[neighbor] / d_out[block])
            end
        end
    end
    model_S_term = B^2 / E
    model_S = E * (1 + model_S_term) * log(1 + model_S_term) -
        model_S_term * log(model_S_term) + N*log(B)
    S = model_S + summation_term
    return S
end


function compute_hastings_correction(
        s::Int64, M::Stinger, M_r_row::Dict{Int64, Int64},
        M_r_col::Dict{Int64, Int64}, B::Int64, d::Vector{Int64},
        d_new::Vector{Int64},
        blocks_out_count_map, blocks_in_count_map
    )
    blocks = Set(keys(blocks_out_count_map)) ∪ Set(keys(blocks_in_count_map))
    p_forward = 0.0
    p_backward = 0.0
    foralledges(M, s) do edge, src, etype
        direction, t = edgeparse(edge)
        if t in blocks
            degree = get(blocks_out_count_map, t, 0) +
                get(blocks_in_count_map, t, 0)
            m = 0.0
            if direction != 1
                m += edge.weight
            end
            if direction != 2
                m += edgeweight(M, t, s, 0)
            end
            p_forward += degree * (m + 1) / (d[t] + B)
            p_backward += degree * (get(M_r_row, t, 0) + get(M_r_col, t, 0) + 1) / (d_new[t] + B)
        end
    end
    return p_backward / p_forward
end

#=
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
=#
function zeros_interblock_edge_matrix(::Type{Stinger}, size::Int64)
    return Stinger(stingerconfig(0))
end
