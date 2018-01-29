using NamedTuples, DataStreams, LibPQ, Missings

export InterblockEdgeCountPostgres

const edgelist_tuple = @NT(block_num::Array{Int32}, src_block::Array{Int32}, dst_block::Array{Int32}, edgecount::Array{Int32})

immutable InterblockEdgeCountPostgres
    conn::Connection

    function InterblockEdgeCountPostgres(;dbname::String="rohitvarkey")
        conn = Connection("dbname=$dbname")
        edgelist_create_stmt = """CREATE TABLE edgelist (
            block_num   integer NOT NULL,
            src_block   integer NOT NULL,
            dst_block   integer NOT NULL,
            edgecount   integer NOT NULL
        );
        """
        blockmap_create_stmt = """CREATE TABLE blockmap (
            block_num   integer NOT NULL,
            vertexid   integer NOT NULL,
            block   integer NOT NULL
        );
        """
        execute(conn, "DROP TABLE IF EXISTS edgelist, blockmap;")
        execute(conn, edgelist_create_stmt)
        execute(conn, blockmap_create_stmt)

        new(conn)
    end
end



function initialize_edge_counts!(
    M::InterblockEdgeCountPostgres, g::SimpleWeightedDiGraph, B::Int64,
    b::Vector{Int64}, count_log::CountLog
    )

    block_edge_counts = Dict{Int64, Dict{Int64, Int64}}()

    total_block_edges = 0
    # Calculate edge block counts.
    for edge in edges(g)
        s, d = b[src(edge)], b[dst(edge)]
        if s in keys(block_edge_counts)
            if !(d in keys(block_edge_counts[s]))
                total_block_edges += 1
            end
            block_edge_counts[s][d] = get(block_edge_counts[s], d, 0) + 1

        else
            block_edge_counts[s] = Dict(d=>1)
            total_block_edges += 1
        end
    end

    block_num = fill(Int32(B), total_block_edges)
    src_block = ones(Int32, total_block_edges)
    dst_block = ones(Int32, total_block_edges)
    edgecounts = ones(Int32, total_block_edges)
    counter = 0
    # Create NamedTuple with the above data.
    for (src, edges) in block_edge_counts
        for (dst, edgecount) in block_edge_counts[src]
            counter += 1
            src_block[counter] = src
            dst_block[counter] = dst
            edgecounts[counter] = edgecount
        end
    end

    data = edgelist_tuple(block_num, src_block, dst_block, edgecounts)

    stmt = Data.stream!(
        data,
        Statement,
        M.conn,
        "INSERT INTO edgelist VALUES (\$1, \$2, \$3, \$4);",
    )
end

function initialize_edge_counts(
    _::Type{InterblockEdgeCountPostgres}, g::SimpleWeightedDiGraph, B::Int64,
    b::Vector{Int64}, count_log::CountLog
    )
    M = InterblockEdgeCountPostgres()
    initialize_edge_counts!(M, g, B, b, count_log)
    M
end

function compute_block_neighbors_and_degrees(
    p::Partition{InterblockEdgeCountPostgres}, block::Int64, count_log::CountLog
    )
    out_neighbors = fetch!(
        NamedTuple,
        execute(p.M.conn, "SELECT dst_block, edgecount FROM edgelist WHERE src_block = $block and block_num=$(p.B);")
    )
    in_neighbors = fetch!(
        NamedTuple,
        execute(p.M.conn, "SELECT src_block, edgecount FROM edgelist WHERE dst_block = $block and block_num=$(p.B);")
    )
    neighbors = collect(Int64,
        Set(skipmissing(out_neighbors[:dst_block])) ∪
        Set(skipmissing(in_neighbors[:src_block]))
    )
    k_out = sum(out_neighbors[:edgecount])
    k_in = sum(in_neighbors[:edgecount])
    k = k_out + k_in
    neighbors, k_out, k_in, k
end

function compute_block_degrees(M::InterblockEdgeCountPostgres, B::Int64, count_log::CountLog)

    d_out_query = fetch!(
        NamedTuple,
        execute(M.conn, "SELECT src_block, SUM(edgecount) FROM edgelist WHERE block_num=$B GROUP BY src_block;")
    )
    d_in_query = fetch!(
        NamedTuple,
        execute(M.conn, "SELECT dst_block, SUM(edgecount) FROM edgelist WHERE block_num=$B GROUP BY dst_block;")
    )

    d_out = zeros(Int64, B)
    d_in = zeros(Int64, B)

    for (src_block, degree) in zip(d_out_query[:src_block], d_out_query[:sum])
        d_out[src_block] = degree
    end

    for (dst_block, degree) in zip(d_in_query[:dst_block], d_in_query[:sum])
        d_in[dst_block] = degree
    end

    d = d_out .+ d_in

    return d_out, d_in, d
end

function compute_new_matrix(
    p::Partition{InterblockEdgeCountPostgres}, r::Int64, s::Int64,
    out_block_count_map, in_block_count_map, self_edge_weight::Int64,
    count_log::CountLog
    )

    M_r_col = zeros(Int64, p.B)
    M_r_row = zeros(Int64, p.B)
    M_s_col = zeros(Int64, p.B)
    M_s_row = zeros(Int64, p.B)

    # Setup with initial counts
    for block in (r, s)
        if block == r
            col = M_r_col
            row = M_r_row
        else
            col = M_s_col
            row = M_s_row
        end
        out_neighbors = fetch!(
            NamedTuple,
            execute(p.M.conn, "SELECT dst_block, edgecount FROM edgelist WHERE src_block = $block and block_num=$(p.B);")
        )
        in_neighbors = fetch!(
            NamedTuple,
            execute(p.M.conn, "SELECT src_block, edgecount FROM edgelist WHERE dst_block = $block and block_num=$(p.B);")
        )
        for (dst_block, edgecount) in zip(out_neighbors[:dst_block], out_neighbors[:edgecount])
            col[dst_block] = edgecount
        end
        for (src_block, edgecount) in zip(in_neighbors[:src_block], in_neighbors[:edgecount])
            row[src_block] = edgecount
        end
    end

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
    p::Partition{InterblockEdgeCountPostgres}, r::Int64, s::Int64, count_log::CountLog
    )

    M_r_row = zeros(Int64, p.B)
    M_r_col = zeros(Int64, p.B)

    M_s_row = zeros(Int64, p.B)
    M_s_col = zeros(Int64, p.B)

    # Setup with initial counts
    for block in (r, s)
        out_neighbors = fetch!(
            NamedTuple,
            execute(p.M.conn, "SELECT dst_block, edgecount FROM edgelist WHERE src_block = $block and block_num=$(p.B);")
        )
        in_neighbors = fetch!(
            NamedTuple,
            execute(p.M.conn, "SELECT src_block, edgecount FROM edgelist WHERE dst_block = $block and block_num=$(p.B);")
        )
        for (dst_block, edgecount) in zip(out_neighbors[:dst_block], out_neighbors[:edgecount])
            M_s_col[dst_block] += edgecount
            # s->r , r->r
            if dst_block == r && block == r
                M_s_row[s] += edgecount
                M_s_col[s] += edgecount
            elseif dst_block == r && block == s
                M_s_col[s] += edgecount
            elseif dst_block == s && block == r
                M_s_row[s] += edgecount
            end
        end
        for (src_block, edgecount) in zip(in_neighbors[:src_block], in_neighbors[:edgecount])
            # r->s , r->r
            M_s_row[src_block] += edgecount
        end
    end

    M_s_row[r] = 0
    M_s_col[r] = 0

    return M_r_row, M_r_col, M_s_row, M_s_col
end

function compute_multinomial_probs(
    p::Partition{InterblockEdgeCountPostgres}, vertex::Int64, count_log::CountLog
    )
    out_counts = fetch!(
        NamedTuple,
        execute(p.M.conn, """
        SELECT dst_block as neighbor, edgecount FROM edgelist WHERE
        src_block=$vertex and block_num=$(p.B);
        """
        )
    )
    in_counts = fetch!(
        NamedTuple,
        execute(p.M.conn, """
        SELECT src_block as neighbor, edgecount FROM edgelist WHERE
        dst_block=$vertex and block_num=$(p.B);
        """
        )
    )
    probs = zeros(Int64, p.B)
    for counts in (out_counts, in_counts)
        for (neighbor, edgecount) in zip(out_counts[:neighbor], counts[:edgecount])
            probs[neighbor] += edgecount
        end
    end
    return probs/p.d[vertex]
end

function compute_delta_entropy(
    p::Partition{InterblockEdgeCountPostgres}, r::Int64, s::Int64,
    M_r_col::Array{Int64, 1}, M_s_col::Array{Int64, 1},
    M_r_row::Array{Int64, 1}, M_s_row::Array{Int64, 1},
    d_out_new::Vector{Int64}, d_in_new::Vector{Int64}, count_log::CountLog
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
    edges_query = fetch!(
        NamedTuple,
        execute(p.M.conn, """
        SELECT src_block, dst_block, edgecount FROM edgelist WHERE
        block_num=$(p.B) and src_block IN ($r, $s) or dst_block IN ($r, $s);
        """
        )
    )
    for (src_block, dst_block, edgecount) in zip(
        edges_query[:src_block],
        edges_query[:dst_block],
        edges_query[:edgecount]
    )
        delta += edgecount * log(edgecount / p.d_in[dst_block] / p.d_out[src_block])
    end
    delta
end

function compute_overall_entropy(
        M::InterblockEdgeCountPostgres, d_out::Vector{Int64}, d_in::Vector{Int64},
        B::Int64, N::Int64, E::Int64, count_log::CountLog
    )
    summation_term = 0.0
    edges_query =  fetch!(
        NamedTuple,
        execute(M.conn, """
        SELECT src_block, dst_block, edgecount FROM edgelist WHERE
        block_num=$B;
        """
        )
    )
    for (src_block, dst_block, edgecount) in zip(
        edges_query[:src_block],
        edges_query[:dst_block],
        edges_query[:edgecount]
    )
        summation_term -= edgecount * log(edgecount/ d_in[dst_block] / d_out[src_block])
    end

    model_S_term = B^2 / E
    model_S = E * (1 + model_S_term) * log(1 + model_S_term) -
        model_S_term * log(model_S_term) + N*log(B)
    S = model_S + summation_term
end

function compute_hastings_correction(
        s::Int64, p::Partition{InterblockEdgeCountPostgres}, M_r_row::Vector{Int64},
        M_r_col::Vector{Int64}, d_new::Vector{Int64},
        blocks_out_count_map, blocks_in_count_map, count_log::CountLog
    )
    blocks_out = collect(keys(blocks_out_count_map))
    blocks_in = collect(keys(blocks_in_count_map))
    blocks = Set(blocks_out) ∪ Set(blocks_in)

    out_counts_query = fetch!(
        NamedTuple,
        execute(p.M.conn, """
        SELECT src_block as neighbor, edgecount FROM edgelist WHERE
        src_block in ($(join(blocks_in, ','))) and block_num=$(p.B) and dst_block=$s;
        """
        )
    )
    in_counts_query = fetch!(
        NamedTuple,
        execute(p.M.conn, """
        SELECT dst_block as neighbor, edgecount FROM edgelist WHERE
        dst_block in ($(join(blocks_out, ','))) and block_num=$(p.B) and src_block=$s;
        """
        )
    )
    out_counts = Dict{Int64, Int64}()
    in_counts = Dict{Int64, Int64}()

    for (neighbor, edgecount) in zip(out_counts_query[:neighbor], out_counts_query[:edgecount])
        out_counts[neighbor] = edgecount
    end
    for (neighbor, edgecount) in zip(in_counts_query[:neighbor], in_counts_query[:edgecount])
        in_counts[neighbor] = edgecount
    end
    p_forward = 0.0
    p_backward = 0.0
    for t in blocks
        degree = get(blocks_out_count_map, t, 0) + get(blocks_in_count_map, t, 0)
        edgecounts = get(out_counts, t, 0) + get(in_counts, t, 0)
        p_forward += degree * (edgecounts + 1) / (p.d[t] + p.B)
        p_backward += degree * (M_r_row[t] + M_r_col[t] + 1) / (d_new[t] + p.B)
    end
    return p_backward / p_forward
end

function update_partition(
    M::InterblockEdgeCountPostgres, r::Int64, s::Int64,
    M_r_col::Vector{Int64}, M_s_col::Vector{Int64},
    M_r_row::Vector{Int64}, M_s_row::Vector{Int64},
    count_log::CountLog
    )
    B = length(M_r_col)
    execute(
        M.conn,
        """
        DELETE FROM edgelist WHERE src_block in ($r, $s) or dst_block in ($r, $s) and block_num=$B;
        """
    )

    edges_to_insert = Set{NTuple{4,Int32}}()
    for dst_block in findn(M_r_col)
        push!(edges_to_insert, (B, r, dst_block, M_r_col[dst_block]))
    end
    for dst_block in findn(M_s_col)
        push!(edges_to_insert, (B, s, dst_block, M_s_col[dst_block]))
    end
    for src_block in findn(M_r_row)
        push!(edges_to_insert, (B, src_block, r, M_r_row[src_block]))
    end
    for src_block in findn(M_s_row)
        push!(edges_to_insert, (B, src_block, s, M_s_row[src_block]))
    end

    edges = collect(edges_to_insert)
    total_block_edges = 0
    # Calculate edge block counts.
    block_nums = fill(B, length(edges))
    src_blocks = [edge[2] for edge in edges]
    dst_blocks = [edge[3] for edge in edges]
    edgecounts = [edge[4] for edge in edges]
    data = edgelist_tuple(block_nums, src_blocks, dst_blocks, edgecounts)
    #info("Updated partition")

    stmt = Data.stream!(
        data,
        Statement,
        M.conn,
        "INSERT INTO edgelist VALUES (\$1, \$2, \$3, \$4);",
    )
    M
end

function zeros_interblock_edge_matrix(::Type{InterblockEdgeCountPostgres},
    size::Int64)
    return InterblockEdgeCountPostgres()
end
