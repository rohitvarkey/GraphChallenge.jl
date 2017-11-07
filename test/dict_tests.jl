using LightGraphs
using SimpleWeightedGraphs
using GraphChallenge

function test_initialize_counts(M::InterblockEdgeCountDict, g::SimpleWeightedDiGraph)
    for edge in edges(g)
        @test M.block_in_edges[dst(edge)][src(edge)] == 1
        @test M.block_out_edges[src(edge)][dst(edge)] == 1
    end
end
