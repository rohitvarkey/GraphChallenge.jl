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
    M = InterblockEdgeCountDictDict(block_in_edges, block_out_edges)
    r = 1
    s = 2
    M_r_row, M_r_col, M_s_row, M_s_col =
        compute_new_matrix_agglomerative(M, r, s, 3)
    @test M_r_row == Dict{Int64, Int64}()
    @test M_r_col == Dict{Int64, Int64}()
    @test M_s_row == Dict(2=>22, 3=>11)
    @test M_s_col == Dict(2=>22, 3=>16)
end

function test_compute_new_matrix(::Type{InterblockEdgeCountDictDict})
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
    M = InterblockEdgeCountDictDict(block_in_edges, block_out_edges)
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
    @test M_r_row == Dict{Int64, Int64}(1=>5, 2=>3, 3=>3)
    @test M_r_col == Dict{Int64, Int64}(1=>5, 2=>2, 3=>1)
    @test M_s_row == Dict(1=>2, 2=>12, 3=>8)
    @test M_s_col == Dict(1=>3, 2=>12, 3=>15)

    @show d_out, d_in, d = compute_block_degrees(M, 3)
    @show overall_entropy = compute_overall_entropy(
        M, d_out, d_in, 3, 3, sum(d)
    )
end

function test_initialize_counts(M::InterblockEdgeCountVectorDict, g::SimpleWeightedDiGraph)
    for edge in edges(g)
        @test M.block_in_edges[dst(edge)][src(edge)] == 1
        @test M.block_out_edges[src(edge)][dst(edge)] == 1
    end
end

function test_compute_new_matrix_agglomerative(::Type{InterblockEdgeCountVectorDict})
    block_out_edges = [
        Dict(1=>8, 2=>2, 3=>4),
        Dict(1=>3, 2=>9, 3=>12),
        Dict(1=>5, 2=>6, 3=>10)
    ]
    block_in_edges = [
        Dict(1=>8, 2=>3, 3=>5),
        Dict(1=>2, 2=>9, 3=>6),
        Dict(1=>4, 2=>12, 3=>10)
    ]
    M = InterblockEdgeCountVectorDict(block_in_edges, block_out_edges)
    r = 1
    s = 2
    M_r_row, M_r_col, M_s_row, M_s_col =
        compute_new_matrix_agglomerative(M, r, s, 3)
    @test M_r_row == Dict{Int64, Int64}()
    @test M_r_col == Dict{Int64, Int64}()
    @test M_s_row == Dict(2=>22, 3=>11)
    @test M_s_col == Dict(2=>22, 3=>16)
end

function test_compute_new_matrix(::Type{InterblockEdgeCountVectorDict})
    block_out_edges = [
        Dict(1=>8, 2=>2, 3=>4),
        Dict(1=>3, 2=>9, 3=>12),
        Dict(1=>5, 2=>6, 3=>10)
    ]
    block_in_edges = [
        Dict(1=>8, 2=>3, 3=>5),
        Dict(1=>2, 2=>9, 3=>6),
        Dict(1=>4, 2=>12, 3=>10)
    ]
    M = InterblockEdgeCountVectorDict(block_in_edges, block_out_edges)
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
    @test M_r_row == Dict{Int64, Int64}(1=>5, 2=>3, 3=>3)
    @test M_r_col == Dict{Int64, Int64}(1=>5, 2=>2, 3=>1)
    @test M_s_row == Dict(1=>2, 2=>12, 3=>8)
    @test M_s_col == Dict(1=>3, 2=>12, 3=>15)

    @show d_out, d_in, d = compute_block_degrees(M, 3)
    @show overall_entropy = compute_overall_entropy(
        M, d_out, d_in, 3, 3, sum(d)
    )
end
