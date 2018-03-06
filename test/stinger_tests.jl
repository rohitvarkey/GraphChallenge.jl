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

    p = test_partition(InterblockEdgeCountStinger)
    r = 1
    s = 2
    @show M_r_row, M_r_col, M_s_row, M_s_col =
        compute_new_matrix_agglomerative(p, r, s, CountLog())
    @test M_r_row == [0, 0, 0]
    @test M_r_col == [0, 0, 0]
    @test M_s_row == [0, 22, 11]
    @test M_s_col == [0, 22, 16]
end

function test_compute_new_matrix(::Type{InterblockEdgeCountStinger})

    p = test_partition(InterblockEdgeCountStinger)

    r = 1
    s = 2
    block_out_count_map = Dict(
        1=>1, 2=>2, 3=>3
    )
    block_in_count_map = Dict(
        1=>2, 2=>1, 3=>2
    )
    M_r_row, M_r_col, M_s_row, M_s_col =
        compute_new_matrix(p, r, s, block_out_count_map, block_in_count_map, 0, CountLog())
    @test M_r_row == [5, 3, 3]
    @test M_r_col == [5, 2, 1]
    @test M_s_row == [2, 12, 8]
    @test M_s_col == [3, 12, 15]

    M = update_partition(
        p.M, r, s,
        M_r_col, M_s_col, M_r_row, M_s_row, CountLog()
    )

    @test M.self_edge_counts == [5, 12, 10]
    @test M.outdegrees == [5+2+1, 3+12+15, 5+6+10]
    @test M.indegrees == [5+3+3, 2+12+8, 4+12+10]

    @show d_out, d_in, d = compute_block_degrees(M, 3, CountLog())
    @show overall_entropy = compute_overall_entropy(
        M, d_out, d_in, 3, 3, sum(d), CountLog()
    )
end
