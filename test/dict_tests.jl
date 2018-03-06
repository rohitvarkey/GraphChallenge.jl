using LightGraphs
using SimpleWeightedGraphs
using GraphChallenge

function test_initialize_counts(M::InterblockEdgeCountDictDict, g::SimpleWeightedDiGraph)
    for edge in edges(g)
        @test M.block_in_edges[dst(edge)][src(edge)] == 1
        @test M.block_out_edges[src(edge)][dst(edge)] == 1
    end
end

function test_compute_new_matrix_agglomerative(::Type{InterblockEdgeCountDictDict})
    p = test_partition(InterblockEdgeCountDictDict)
    r = 1
    s = 2
    M_r_row, M_r_col, M_s_row, M_s_col =
        compute_new_matrix_agglomerative(p, r, s, CountLog())
    @test M_r_row == Dict{Int64, Int64}()
    @test M_r_col == Dict{Int64, Int64}()
    @test M_s_row == Dict(2=>22, 3=>11)
    @test M_s_col == Dict(2=>22, 3=>16)
end

function test_compute_new_matrix(::Type{InterblockEdgeCountDictDict})
    p = test_partition(InterblockEdgeCountDictDict)
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
    @test M_r_row == Dict{Int64, Int64}(1=>5, 2=>3, 3=>3)
    @test M_r_col == Dict{Int64, Int64}(1=>5, 2=>2, 3=>1)
    @test M_s_row == Dict(1=>2, 2=>12, 3=>8)
    @test M_s_col == Dict(1=>3, 2=>12, 3=>15)

    @show d_out, d_in, d = compute_block_degrees(p.M, 3, CountLog())
    @show overall_entropy = compute_overall_entropy(
        p.M, d_out, d_in, 3, 3, sum(d), CountLog()
    )
end

function test_initialize_counts(M::InterblockEdgeCountVectorDict, g::SimpleWeightedDiGraph)
    for edge in edges(g)
        @test M.block_in_edges[dst(edge)][src(edge)] == 1
        @test M.block_out_edges[src(edge)][dst(edge)] == 1
    end
end

function test_compute_new_matrix_agglomerative(::Type{InterblockEdgeCountVectorDict})
    p = test_partition(InterblockEdgeCountVectorDict)
    r = 1
    s = 2
    M_r_row, M_r_col, M_s_row, M_s_col =
        compute_new_matrix_agglomerative(p, r, s, CountLog())
    @test M_r_row == Dict{Int64, Int64}()
    @test M_r_col == Dict{Int64, Int64}()
    @test M_s_row == Dict(2=>22, 3=>11)
    @test M_s_col == Dict(2=>22, 3=>16)
end

function test_compute_new_matrix(::Type{InterblockEdgeCountVectorDict})
    p = test_partition(InterblockEdgeCountVectorDict)
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
    @test M_r_row == Dict{Int64, Int64}(1=>5, 2=>3, 3=>3)
    @test M_r_col == Dict{Int64, Int64}(1=>5, 2=>2, 3=>1)
    @test M_s_row == Dict(1=>2, 2=>12, 3=>8)
    @test M_s_col == Dict(1=>3, 2=>12, 3=>15)

    @show d_out, d_in, d = compute_block_degrees(p.M, 3, CountLog())
    @show overall_entropy = compute_overall_entropy(
        p.M, d_out, d_in, 3, 3, sum(d), CountLog()
    )
end
