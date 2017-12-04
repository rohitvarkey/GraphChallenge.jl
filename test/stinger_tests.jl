using LightGraphs
using SimpleWeightedGraphs
using StingerGraphs
using GraphChallenge

function test_initialize_counts(M::InterblockEdgeCountStinger, g::SimpleWeightedDiGraph)
    for edge in edges(g)
        @test edgeweight(M.s, src(edge), dst(edge), 0) == 1
    end
end

function test_compute_new_matrix_agglomerative(::Type{InterblockEdgeCountStinger})
    block_out_edges = Dict(
        1 => Dict(1=>8, 2=>2, 3=>4),
        2 => Dict(1=>3, 2=>9, 3=>12),
        3 => Dict(1=>5, 2=>6, 3=>10)
    )
    block_in_edges = Dict(
        1 => Dict(1=>8, 2=>3, 3=>5),
        2 => Dict(1=>2, 2=>9, 3=>6),
        3 => Dict(1=>4, 2=>12, 3=>10)
    )
    M = InterblockEdgeCountStinger(Stinger(stingerconfig(4)), zeros(3), zeros(3), zeros(3))
    for (src_block, edges) in block_out_edges
        for (dst_block, edgecount) in edges
            if src_block != dst_block
                insert_edge!(M.s, 0, src_block, dst_block, edgecount, 1)
            end
        end
    end

    for (block, out_edges) in  block_out_edges
        M.outdegrees[block] = sum(values(out_edges))
    end
    for (block, in_edges) in  block_in_edges
        M.indegrees[block] = sum(values(in_edges))
    end

    for i=1:3
        M.self_edge_counts[i] = block_out_edges[i][i]
    end
    @show M
    r = 1
    s = 2
    @show M_r_row, M_r_col, M_s_row, M_s_col =
        compute_new_matrix_agglomerative(M, r, s, 3)
    @test M_r_row == [0, 0, 0]
    @test M_r_col == [0, 0, 0]
    @test M_s_row == [0, 22, 11]
    @test M_s_col == [0, 22, 16]
end

function test_compute_new_matrix(::Type{InterblockEdgeCountStinger})
    block_out_edges = Dict(
        1 => Dict(1=>8, 2=>2, 3=>4),
        2 => Dict(1=>3, 2=>9, 3=>12),
        3 => Dict(1=>5, 2=>6, 3=>10)
    )
    block_in_edges = Dict(
        1 => Dict(1=>8, 2=>3, 3=>5),
        2 => Dict(1=>2, 2=>9, 3=>6),
        3 => Dict(1=>4, 2=>12, 3=>10)
    )
    M = InterblockEdgeCountStinger(Stinger(stingerconfig(4)), zeros(3), zeros(3), zeros(3))
    for (src_block, edges) in block_out_edges
        for (dst_block, edgecount) in edges
            if src_block != dst_block
                insert_edge!(M.s, 0, src_block, dst_block, edgecount, 1)
            end
        end
    end

    for (block, out_edges) in  block_out_edges
        M.outdegrees[block] = sum(values(out_edges))
    end
    for (block, in_edges) in  block_in_edges
        M.indegrees[block] = sum(values(in_edges))
    end

    for i=1:3
        M.self_edge_counts[i] = block_out_edges[i][i]
    end
    @show M

    r = 1
    s = 2
    block_out_count_map = Dict(
        1=>1, 2=>2, 3=>3
    )
    block_in_count_map = Dict(
        1=>2, 2=>1, 3=>2
    )
    M_r_row, M_r_col, M_s_row, M_s_col =
        compute_new_matrix(M, r, s, 3, block_out_count_map, block_in_count_map, 0)
    @test M_r_row == [5, 3, 3]
    @test M_r_col == [5, 2, 1]
    @test M_s_row == [2, 12, 8]
    @test M_s_col == [3, 12, 15]

    @show d_out, d_in, d = compute_block_degrees(M, 3)
    @show overall_entropy = compute_overall_entropy(
        M, d_out, d_in, 3, 3, sum(d)
    )

    M = update_partition(
        M, r, s,
        M_r_col, M_s_col, M_r_row, M_s_row
    )

    @test M.self_edge_counts == [5, 12, 10]
    @test M.outdegrees == [5+2+1, 3+12+15, 5+6+10]
    @test M.indegrees == [5+3+3, 2+12+8, 4+12+10]

    @show d_out, d_in, d = compute_block_degrees(M, 3)
    @show overall_entropy = compute_overall_entropy(
        M, d_out, d_in, 3, 3, sum(d)
    )
end
