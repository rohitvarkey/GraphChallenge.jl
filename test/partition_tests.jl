import GraphChallenge: compute_block_neighbors_and_degrees,
                       propose_new_partition_agg,
                       evaluate_proposal_agg,
                       compute_new_matrix_agglomerative

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
        block_partition = collect(1:num_nodes)
        M = initialize_edge_counts(
            T, g, num_nodes, block_partition
        )
        d_out, d_in, d = compute_block_degrees(M, num_nodes)
        current_block = 1
        num_blocks = floor(Int64, 0.5 * num_nodes)
        block_neighbors, k_out, k_in, k = compute_block_neighbors_and_degrees(
            M, current_block
        )

        @test Set(all_neighbors(g, current_block)) == Set(block_neighbors)
        @test k_out == outdegree(g, current_block)
        @test k_in == indegree(g, current_block)
        @test k == degree(g, current_block)

        srand(42) #Seed the RNG

        proposal = propose_new_partition_agg(
            M, current_block, block_partition, 25,
            d, block_neighbors
        )
        @test proposal == 2

        Δ = evaluate_proposal_agg(
            M, current_block, proposal, num_blocks, d, d_in, d_out,
            k, k_in, k_out
        )
        @test round(Δ, 8) == round(12.922923932215921, 8)
    end
end

function test_matrix_update_agglomerative()
    M_r_row, M_r_col, M_s_row, M_s_col =
        compute_new_matrix_agglomerative(M, current_block, proposal, num_blocks)

    @test all(x->x == 0, M_r_row)
    @test all(x->x == 0, M_r_col)
    @test Set(findn(M_s_col)) == Set([6, 14, 15, 17, 18, 27, 35, 36, 48])
    @test Set(findn(M_s_row)) == Set([2, 7, 9, 12, 15, 16, 20, 21, 22, 33, 37, 40, 42, 45])
end
