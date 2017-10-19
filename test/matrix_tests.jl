using LightGraphs
using SimpleWeightedGraphs
using GraphChallenge

function test_initialize_counts(M::Array{Int64, 2}, g::SimpleWeightedDiGraph)
    for edge in edges(g)
        @test M[dst(edge), src(edge)] == 1
    end
end
