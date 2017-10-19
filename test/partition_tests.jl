import GraphChallenge: compute_block_neighbors_and_degrees

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
        num_nodes = 50
        g = load_graph(50)
        M = initialize_edge_counts(
            T, g, num_nodes, collect(1:num_nodes)
        )
        test_initialize_counts(M, g)
        test_compute_block_degrees(M, g, num_nodes)
    end
end

@testset "Agglomerative step" begin
    for T in (Array{Int64, 2}, )
        num_nodes = 50
        g = load_graph(50)
        M = initialize_edge_counts(
            T, g, num_nodes, collect(1:num_nodes)
        )
        d_out, d_in, d = compute_block_degrees(M, num_nodes)
        current_block = 1
        block_neighbors, k_out, k_in, k = compute_block_neighbors_and_degrees(
            M, current_block
        )
        @test Set(all_neighbors(g, current_block)) == Set(block_neighbors)
        @test k_out == outdegree(g, current_block)
        @test k_in == indegree(g, current_block)
        @test k == degree(g, current_block)
    end
end
