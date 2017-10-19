using LightGraphs
using SimpleWeightedGraphs
using GraphChallenge

function test_initialize_counts(M::Array{Int64, 2}, g::SimpleWeightedDiGraph)
    for edge in edges(g)
        @test M[dst(edge), src(edge)] == 1
    end
end

#test_initialize_counts(10)

#partition(Array{Int64, 2}, "emergingEdges", 1000)

#partition(Array{Int64, 2}, 1000)
