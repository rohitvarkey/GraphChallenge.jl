import StingerGraphs: Stinger, remove_edge!, insert_edge!, indegree, outdegree, foralledges, edgeparse, stingerconfig, edgeweight

immutable InterblockEdgeCountStinger
    s::Stinger
    self_edge_counts::Vector{Int64} #to store self edge counts
    outdegrees::Vector{Int64}
    indegrees::Vector{Int64}
end

"""
Initializes the edge count matrix M between the blocks.
Calculates the new interblock_matrix edge count matrix.
"""
function initialize_edge_counts!(
    M::InterblockEdgeCountStinger, g::SimpleWeightedDiGraph, B::Int64,
    b::Vector{Int64}
    )
    #TODO: Optimize this by using batch insertions.
    edge_counts = Dict{Pair{Int64, Int64}, Int64}()
    for edge in edges(g)
        s, d = b[src(edge)], b[dst(edge)]
        edge_counts[(s=>d)] = get(edge_counts, s=>d, 0) + 1
        M.outdegrees[s] += 1
        M.indegrees[d] += 1
    end
    for ((s, d), edgecount) in edge_counts
        if s == d
            M.self_edge_counts[s] = edgecount
        else
            # Do not add self edges to stinger
            insert_edge!(M.s, 0, s, d, edgecount, 1)
        end
    end
    M
end

function initialize_edge_counts(
    _::Type{InterblockEdgeCountStinger}, g::SimpleWeightedDiGraph, B::Int64,
    b::Vector{Int64}
    )
    M = InterblockEdgeCountStinger(
        Stinger(stingerconfig(B + 1)), # To prevent off by one errors
        zeros(Int64, B), zeros(Int64, B), zeros(Int64, B)
    )
    initialize_edge_counts!(M, g, B, b)
    M
end

function compute_block_neighbors_and_degrees(M::InterblockEdgeCountStinger, block::Int64)
    neighbors = Set{Int64}()
    k_in = M.self_edge_counts[block]
    k_out = M.self_edge_counts[block]
    foralledges(M.s, block) do edge, src, etype
        direction, neighbor = edgeparse(edge)
        push!(neighbors, neighbor)
        if direction != 1
            k_out += edge.weight
        end
        if direction != 2
            k_in += edgeweight(M.s, neighbor, block, 0)
        end
    end
    collect(neighbors), k_out, k_in, k_in + k_out
end


function compute_block_degrees(M::InterblockEdgeCountStinger, B::Int64)
    # Sum across rows to get the outdegrees for each block
    d_out = copy(M.outdegrees)
    d_in = copy(M.indegrees)
    d = d_out + d_in
    return d_out, d_in, d
end

function compute_new_matrix(
    M::InterblockEdgeCountStinger, r::Int64, s::Int64, num_blocks::Int64,
    out_block_count_map, in_block_count_map, self_edge_weight::Int64
    )

    M_r_col = zeros(Int64, num_blocks)
    M_r_row = zeros(Int64, num_blocks)
    M_s_col = zeros(Int64, num_blocks)
    M_s_row = zeros(Int64, num_blocks)

    M_r_col[r] = M.self_edge_counts[r]
    M_r_row[r] = M.self_edge_counts[r]
    M_s_col[s] = M.self_edge_counts[s]
    M_s_col[s] = M.self_edge_counts[s]

    #@show r, s
    #@show out_block_count_map
    #@show in_block_count_map
    foralledges(M.s, r) do edge, src, etype
        direction, block = edgeparse(edge)
        # Move outgoing edges from r to s.
        if direction != -1
            if direction != 1
                M_r_col[block] += edge.weight
                # out edge
                if block in keys(out_block_count_map)
                    #@show r, block, edge.weight, out_block_count_map[block]
                    out_count = out_block_count_map[block]
                    M_r_col[block] -= out_count
                    M_s_col[block] += out_count
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
            end
            if direction != 2
                # in edge
                M_r_row[block] += edgeweight(M.s, block, r, 0)
                if block in keys(in_block_count_map)
                    in_count = get(in_block_count_map, block, 0)
                    #@show r, block, edgeweight(M.s, block, r, 0), in_count
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
            end
        end
    end

    #@show M_r_row, M_r_col, M_s_row, M_s_col
    foralledges(M.s, s) do edge, src, etype
        direction, block = edgeparse(edge)
        if direction != -1 && direction != 1
            M_s_col[block] += edge.weight
        end
        if direction != -1 && direction != 2
            M_s_row[block] += edgeweight(M.s, block, s, 0)
        end
    end
    #@show M_r_row, M_r_col, M_s_row, M_s_col


    if self_edge_weight > 0
        # Correct counts based on self_edge_weight
        M_s_row[r] -= self_edge_weight
        M_s_row[s] += self_edge_weight
        M_s_col[r] -= self_edge_weight
        M_s_col[s] += self_edge_weight
    end

    #@show M_r_row, M_r_col, M_s_row, M_s_col
    if any(any.(map.(x->x < 0, [M_r_row, M_r_col, M_s_row, M_s_col])))
        @show any.(map.(x->x < 0, [M_r_row, M_r_col, M_s_row, M_s_col]))
    end
    return M_r_row, M_r_col, M_s_row, M_s_col
end

"""Computes the new rows and cols in `M`, when all nodes from `r` are shifted to
block `s`."""
function compute_new_matrix_agglomerative(
    M::InterblockEdgeCountStinger, r::Int64, s::Int64, num_blocks::Int64
    )

    M_r_col = zeros(Int64, num_blocks)
    M_r_row = zeros(Int64, num_blocks)
    M_s_col = zeros(Int64, num_blocks)
    M_s_row = zeros(Int64, num_blocks)

    foralledges(M.s, s) do edge, src, etype
        direction, block = edgeparse(edge)
        if direction != -1 && direction != 1
            # out edges
            M_s_col[block] += edge.weight
            if block == r
                # s->r edge
                M_s_row[s] += edge.weight
                M_s_col[s] += edge.weight
            end
        end
        if direction != -1 && direction != 2
            M_s_row[block] += edgeweight(M.s, block, s, 0)
        end
    end

    foralledges(M.s, r) do edge, src, etype
        direction, block = edgeparse(edge)
        # Move outgoing edges from r to s.
        if direction != -1 && direction != 1
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
        if direction != -1 && direction != 2 && block != s
            # in edge
            M_s_row[block] += edgeweight(M.s, block, r, 0)
        end
    end

    M_s_row[r] = 0
    M_s_col[r] = 0

    return M_r_row, M_r_col, M_s_row, M_s_col
end


function compute_multinomial_probs(
    M::InterblockEdgeCountStinger, degrees::Vector{Int64}, block::Int64
    )
    probabilities = zeros(length(degrees))
    foralledges(M.s, block) do edge, src, etype
        direction, neighbor = edgeparse(edge)
        if direction != -1 && direction != 1
            # out edge
            probabilities[neighbor] += edge.weight
        end
        if direction != -1 && direction != 2
            # in edge
            probabilities[neighbor] += edgeweight(M.s, neighbor, block, 0)
        end
    end
    return probabilities
end

function compute_delta_entropy(
    M::InterblockEdgeCountStinger, r::Int64, s::Int64,
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
        #@show t1, M_r_col[t1], d_in_new[t1], d_out_new[r]
        delta -= M_r_col[t1] * log(M_r_col[t1] / d_in_new[t1] / d_out_new[r])
    end
    for t1 in findn(M_s_col)
        if t1 ∈ (r, s)
            continue
        end
        #@show t1, M_s_col[t1], d_in_new[t1], d_out_new[s]
        delta -= M_s_col[t1] * log(M_s_col[t1] / d_in_new[t1] / d_out_new[s])
    end
    # Sum over row of r in new M
    for t2 in findn(M_r_row)
        #@show t2, M_r_row[t2], d_in_new[r], d_out_new[t2]
        delta -= M_r_row[t2] * log(M_r_row[t2] / d_in_new[r] / d_out_new[t2])
    end
    # Sum over row of s in new M
    for t2 in findn(M_s_row)
        #@show t2, M_s_row[t2], d_in_new[s], d_out_new[t2]
        delta -= M_s_row[t2] * log(M_s_row[t2] / d_in_new[s] / d_out_new[t2])
    end
    # Sum over edges in old M
    for block in (r, s)
        foralledges(M.s, block) do edge, src, etype
            direction, neighbor = edgeparse(edge)
            if direction != -1 && direction != 1
                # edge is block -> neighbor
                #@show block, neighbor, edge.weight , d_in[neighbor] , d_out[block]
                delta += edge.weight * log(edge.weight / d_in[neighbor] / d_out[block])
            end
            # Prevent double counting the r->s, s->r edges.
            if direction != -1 && direction != 2 && !((block, neighbor) == (r, s) || (block, neighbor) == (s, r))
                # edge is neighbor -> block
                edgecount = edgeweight(M.s, neighbor, block, 0)
                #@show block, neighbor, edgecount , d_in[block] , d_out[neighbor]
                delta += edgecount * log(edgecount / d_in[block] / d_out[neighbor])
            end
        end
    end
    #println("Delta: $delta")
    delta
end

function compute_overall_entropy(
        M::InterblockEdgeCountStinger, d_out::Vector{Int64}, d_in::Vector{Int64},
        B::Int64, N::Int64, E::Int64
    )
    summation_term = 0.0
    for block=1:B
        foralledges(M.s, B) do edge, src, etype
            direction, neighbor = edgeparse(edge)
            if direction != -1 && direction != 1
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
        s::Int64, M::InterblockEdgeCountStinger, M_r_row::Vector{Int64},
        M_r_col::Vector{Int64}, B::Int64, d::Vector{Int64},
        d_new::Vector{Int64},
        blocks_out_count_map, blocks_in_count_map
    )
    blocks = Set(keys(blocks_out_count_map)) ∪ Set(keys(blocks_in_count_map))
    p_forward = 0.0
    p_backward = 0.0
    foralledges(M.s, s) do edge, src, etype
        direction, t = edgeparse(edge)
        if t in blocks
            degree = get(blocks_out_count_map, t, 0) +
                get(blocks_in_count_map, t, 0)
            m = 0.0
            if direction != -1 && direction != 1
                m += edge.weight
            end
            if direction != -1 && direction != 2
                m += edgeweight(M.s, t, s, 0)
            end
            p_forward += degree * (m + 1) / (d[t] + B)
            p_backward += degree * (get(M_r_row, t, 0) + get(M_r_col, t, 0) + 1) / (d_new[t] + B)
        end
    end
    return p_backward / p_forward
end

function update_partition(
    M::InterblockEdgeCountStinger, r::Int64, s::Int64,
    M_r_col::Vector{Int64}, M_s_col::Vector{Int64},
    M_r_row::Vector{Int64}, M_s_row::Vector{Int64}
    )
    println("Updating partition for $r, $s")

    M.self_edge_counts[r] = M_r_col[r]
    M.self_edge_counts[s] = M_s_col[s]

    # Setting degree counts for r and s
    M.outdegrees[r] = sum(M_r_col)
    M.indegrees[r] = sum(M_r_row)
    M.outdegrees[s] = sum(M_s_col)
    M.indegrees[s] = sum(M_s_row)

    # Setting to 0 as we don't care about these anymore.
    M_r_col[r] = 0
    M_r_row[r] = 0
    M_s_col[s] = 0
    M_s_row[s] = 0

    # Updating edges that already exist
    foralledges(M.s, r) do edge, src, etype
        direction, neighbor = edgeparse(edge)
        #@show direction, neighbor, edge.weight, edgeweight(M.s, 0, neighbor, r)
        if direction != -1 && direction != 1 && neighbor != r
            # Update indegree count of neighbor and remove the edge
            M.indegrees[neighbor] -= (edge.weight - M_r_col[neighbor])
            if M_r_col[neighbor] == 0
                remove_edge!(M.s, 0, r, neighbor)
            else
                insert_edge!(M.s, 0, r, neighbor, M_r_col[neighbor], 1)
                M_r_col[neighbor] = 0
            end
        end
        if direction != -1 && direction != 2 && neighbor != r
            # Update outdegree count of neighbor and remove the edge
            M.outdegrees[neighbor] -= (edgeweight(M.s, neighbor, r, 0) - M_r_row[neighbor])
            if M_r_row[neighbor] == 0
                remove_edge!(M.s, 0, neighbor, r)
            else
                insert_edge!(M.s, 0, neighbor, r, M_r_row[neighbor], 1)
                M_r_row[neighbor] = 0
            end
        end
    end

    foralledges(M.s, s) do edge, src, etype
        direction, neighbor = edgeparse(edge)
        if direction != -1 && direction != 1 && neighbor != s
            M.indegrees[neighbor] -= (edge.weight - M_s_col[neighbor])
            if M_s_col[neighbor] == 0
                remove_edge!(M.s, 0, s, neighbor)
            else
                insert_edge!(M.s, 0, s, neighbor, M_s_col[neighbor], 1)
                M_s_col[neighbor] = 0
            end
        end
        if direction != -1 && direction != 2 && neighbor != s
            M.outdegrees[neighbor] -= (edgeweight(M.s, neighbor, s, 0) - M_s_row[neighbor])
            if M_s_row[neighbor] == 0
                remove_edge!(M.s, 0, neighbor, s)
            else
                insert_edge!(M.s, 0, neighbor, s, M_s_row[neighbor], 1)
                M_s_row[neighbor] = 0
            end
        end
    end

    # Adding new edges
    for idx in findn(M_r_col)
        insert_edge!(M.s, 0, r, idx, M_r_col[idx], 1)
        M.indegrees[idx] += M_r_col[idx]
        #@show idx, r, M_r_col[idx], edgeweight(M.s, r, idx, 0)
    end
    for idx in findn(M_r_row)
        insert_edge!(M.s, 0, idx, r, M_r_row[idx], 1)
        M.outdegrees[idx] += M_r_row[idx]
        #@show idx, r, M_r_row[idx], edgeweight(M.s, idx, r, 0)
    end
    for idx in findn(M_s_col)
        insert_edge!(M.s, 0, s, idx, M_s_col[idx], 1)
        M.indegrees[idx] += M_s_col[idx]
        #@show s, idx, M_s_col[idx], edgeweight(M.s, s, idx, 0)
    end
    for idx in findn(M_s_row)
        insert_edge!(M.s, 0, idx, s, M_s_row[idx], 1)
        M.outdegrees[idx] += M_r_row[idx]
        #@show s, idx, M_s_row[idx], edgeweight(M.s, idx, s, 0)
    end
    println("Updated partition")
    M
end

function zeros_interblock_edge_matrix(::Type{InterblockEdgeCountStinger}, size::Int64)
    return InterblockEdgeCountStinger(
        Stinger(stingerconfig(0)), zeros(Int64, 0), zeros(Int64, 0), zeros(Int64, 0))
end
