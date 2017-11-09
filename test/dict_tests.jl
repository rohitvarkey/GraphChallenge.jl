using LightGraphs
using SimpleWeightedGraphs
using GraphChallenge

function test_initialize_counts(M::InterblockEdgeCountDict, g::SimpleWeightedDiGraph)
    for edge in edges(g)
        @test M.block_in_edges[dst(edge)][src(edge)] == 1
        @test M.block_out_edges[src(edge)][dst(edge)] == 1
    end
end

function test_compute_new_matrix_agglomerative(::Type{InterblockEdgeCountDict})
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
    M = InterblockEdgeCountDict(block_in_edges, block_out_edges)
    r = 1
    s = 2
    @show M_r_row, M_r_col, M_s_row, M_s_col =
        compute_new_matrix_agglomerative(M, r, s, 3)
    @test M_r_row == Dict{Int64, Int64}()
    @test M_r_col == Dict{Int64, Int64}()
    @test M_s_row == Dict(2=>22, 3=>11)
    @test M_s_col == Dict(2=>22, 3=>16)
end
