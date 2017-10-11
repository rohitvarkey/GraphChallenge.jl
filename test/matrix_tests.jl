using LightGraphs
using SimpleWeightedGraphs
using GraphChallenge

function test_initialize_counts(num_nodes::Int64)
    g = SimpleWeightedDiGraph(LightGraphs.CompleteDiGraph(num_nodes))
    M = initialize_edge_counts(
        Array{Int64, 2}, g, num_nodes, collect(1:num_nodes)
    )
    d_out, d_in, d = compute_block_degrees(M, num_nodes)
    for edge in edges(g)
        @test M[dst(edge), src(edge)] == 1
    end
    for v in vertices(g)
        @test outdegree(g, v) == d_out[v]
        @test indegree(g, v) == d_in[v]
    end
end

test_initialize_counts(10)

partition(Array{Int64, 2}, "emergingEdges", 1000)

partition(Array{Int64, 2}, 1000)
