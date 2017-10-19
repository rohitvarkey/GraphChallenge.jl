
function test_initialize_counts(
        T::Type, g::SimpleWeightedDiGraph, num_nodes::Int64
    )
    M = initialize_edge_counts(
        T, g, num_nodes, collect(1:num_nodes)
    )
    test_initialize_counts(M, g)
    M
end

function test_compute_block_degrees(M, g::SimpleWeightedDiGraph, num_nodes::Int64)
    d_out, d_in, d = compute_block_degrees(M, num_nodes)
    for v in vertices(g)
        @test outdegree(g, v) == d_out[v]
        @test indegree(g, v) == d_in[v]
        @test degree(g, v) == d[v]
    end
end


@testset "Initialization" begin
    for T in (Array{Int64, 2}, )
        num_nodes = 10
        g = SimpleWeightedDiGraph(LightGraphs.CompleteDiGraph(num_nodes))
        M = test_initialize_counts(T, g, num_nodes)
        test_compute_block_degrees(M, g, num_nodes)
    end
end
