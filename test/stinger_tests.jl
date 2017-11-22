using LightGraphs
using SimpleWeightedGraphs
using StingerGraphs
using GraphChallenge

function test_initialize_counts(M::Stinger, g::SimpleWeightedDiGraph)
    for edge in edges(g)
        @test edgeweight(M, src(edge), dst(edge), 0) == 1
    end
end

function test_compute_new_matrix_agglomerative(::Type{Stinger})
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
    M = Stinger(stingerconfig(4))
    for (src_block, edges) in block_out_edges
        for (dst_block, edgecount) in edges
            insert_edge!(M, 0, src_block, dst_block, edgecount, 1)
        end
    end

    for src=1:3
        println("edges for $src")
        foralledges(M, src) do edge, src, etype
            direction, neighbor = edgeparse(edge)
            @show src, neighbor, direction, edge.weight
        end
    end
    r = 1
    s = 2
    @show M_r_row, M_r_col, M_s_row, M_s_col =
        compute_new_matrix_agglomerative(M, r, s, 3)
    @test M_r_row == [0, 0, 0]
    @test M_r_col == [0, 0, 0]
    @test M_s_row == [0, 22, 11]
    @test M_s_col == [0, 22, 16]
end
