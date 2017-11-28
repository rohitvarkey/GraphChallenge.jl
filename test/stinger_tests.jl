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
        3 => Dict(1=>5, 2=>12, 3=>10)
    )
    M = InterblockEdgeCountStinger(Stinger(stingerconfig(4)), zeros(3), zeros(3), zeros(3))
    for (src_block, edges) in block_out_edges
        for (dst_block, edgecount) in edges
            insert_edge!(M.s, 0, src_block, dst_block, edgecount, 1)
        end
    end

    r = 1
    s = 2
    M_r_row, M_r_col, M_s_row, M_s_col =
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
        3 => Dict(1=>5, 2=>12, 3=>10)
    )
    M = InterblockEdgeCountStinger(Stinger(stingerconfig(4)), zeros(3), zeros(3), zeros(3))
    for (src_block, edges) in block_out_edges
        for (dst_block, edgecount) in edges
            insert_edge!(M.s, 0, src_block, dst_block, edgecount, 1)
        end
    end

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
end
