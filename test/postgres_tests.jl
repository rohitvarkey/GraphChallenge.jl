using LightGraphs
using SimpleWeightedGraphs
using GraphChallenge
using LibPQ, DataStreams, NamedTuples

function test_initialize_counts(M::InterblockEdgeCountPostgres, g::SimpleWeightedDiGraph)
    edges_query = fetch!(
        NamedTuple,
        execute(M.conn, """
        SELECT src_block, dst_block, edgecount FROM edgelist WHERE
        block_num=$(nv(g));
        """
        )
    )
    g_edges = collect(LightGraphs.edges(g))
    for (src_block, dst_block, edgecount) in zip(
        edges_query[:src_block],
        edges_query[:dst_block],
        edges_query[:edgecount]
    )
        @test SimpleWeightedEdge(src_block, dst_block, edgecount) in g_edges
    end
end

function test_compute_new_matrix_agglomerative(::Type{InterblockEdgeCountPostgres})
    M = InterblockEdgeCountPostgres()
    A = [8 3 5; 2 9 6; 4 12 10]

    g = load_graph(50)

    total_block_edges = 9
    block_num = fill(Int32(3), total_block_edges)
    src_block = ones(Int32, total_block_edges)
    dst_block = ones(Int32, total_block_edges)
    edgecounts = ones(Int32, total_block_edges)

    counter = 0
    for col in 1:3
        for row in 1:3
            counter += 1
            src_block[counter] = col
            dst_block[counter] = row
            edgecounts[counter] = A[row, col]
        end
    end

    data = GraphChallenge.edgelist_tuple(block_num, src_block, dst_block, edgecounts)

    stmt = Data.stream!(
        data,
        Statement,
        M.conn,
        "INSERT INTO edgelist VALUES (\$1, \$2, \$3, \$4);",
    )

    d_out, d_in, d = compute_block_degrees(M, 3, CountLog())
    p = Partition(M, g, Inf, collect(1:3), d, d_out, d_in, 3, CountLog())

    r = 1
    s = 2
    M_r_row, M_r_col, M_s_row, M_s_col =
        compute_new_matrix_agglomerative(p, r, s, CountLog())
    @test M_r_row == [0, 0, 0]
    @test M_r_col == [0, 0, 0]
    @test M_s_row == [0, 22, 11]
    @test M_s_col == [0, 22, 16]
end

function test_compute_new_matrix(::Type{InterblockEdgeCountPostgres})
    M = InterblockEdgeCountPostgres()
    A = [8 3 5; 2 9 6; 4 12 10]

    g = load_graph(50)

    r = 1
    s = 2
    block_out_count_map = Dict(
        1=>1, 2=>2, 3=>3
    )
    block_in_count_map = Dict(
        1=>2, 2=>1, 3=>2
    )

    total_block_edges = 9
    block_num = fill(Int32(3), total_block_edges)
    src_block = ones(Int32, total_block_edges)
    dst_block = ones(Int32, total_block_edges)
    edgecounts = ones(Int32, total_block_edges)

    counter = 0
    for col in 1:3
        for row in 1:3
            counter += 1
            src_block[counter] = col
            dst_block[counter] = row
            edgecounts[counter] = A[row, col]
        end
    end

    data = GraphChallenge.edgelist_tuple(block_num, src_block, dst_block, edgecounts)

    stmt = Data.stream!(
        data,
        Statement,
        M.conn,
        "INSERT INTO edgelist VALUES (\$1, \$2, \$3, \$4);",
    )

    @show d_out, d_in, d = compute_block_degrees(M, 3, CountLog())
    @show overall_entropy = compute_overall_entropy(
        M, d_out, d_in, 3, 3, sum(d), CountLog()
    )
    p = Partition(M, g, Inf, collect(1:3), d, d_out, d_in, 3, CountLog())

    M_r_row, M_r_col, M_s_row, M_s_col =
        compute_new_matrix(p, r, s, block_out_count_map, block_in_count_map, 0, CountLog())
    @test M_r_row == [5, 3, 3]
    @test M_r_col == [5, 2, 1]
    @test M_s_row == [2, 12, 8]
    @test M_s_col == [3, 12, 15]

    M = update_partition(
        M, r, s,
        M_r_col, M_s_col, M_r_row, M_s_row, CountLog()
    )

    #@test M.self_edge_counts == [5, 12, 10]
    #@test M.outdegrees == [5+2+1, 3+12+15, 5+6+10]
    #@test M.indegrees == [5+3+3, 2+12+8, 4+12+10]

    @show d_out, d_in, d = compute_block_degrees(M, 3, CountLog())
    @show overall_entropy = compute_overall_entropy(
        M, d_out, d_in, 3, 3, sum(d), CountLog()
    )
end
