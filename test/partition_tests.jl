import GraphChallenge: compute_block_neighbors_and_degrees,
                       propose_new_partition_agg,
                       evaluate_proposal_agg,
                       compute_new_matrix_agglomerative,
                       compute_new_matrix,
                       compute_delta_entropy,
                       carry_out_best_merges

import StingerGraphs: Stinger

function test_compute_block_degrees(M, g::SimpleWeightedDiGraph, num_nodes::Int64)
    d_out, d_in, d = compute_block_degrees(M, num_nodes)
    for v in vertices(g)
        @test LightGraphs.outdegree(g, v) == d_out[v]
        @test LightGraphs.indegree(g, v) == d_in[v]
        @test LightGraphs.degree(g, v) == d[v]
    end
end


@testset "Initialization" begin
    for T in (Array{Int64, 2}, InterblockEdgeCountDictDict, InterblockEdgeCountVectorDict, Stinger)
        num_nodes = 50
        g = load_graph(50)
        M = initialize_edge_counts(
            T, g, num_nodes, collect(1:num_nodes)
        )
        test_initialize_counts(M, g)
        test_compute_block_degrees(M, g, num_nodes)
    end
end


function test_carry_out_best_merges()
    delta_entropy_for_each_block = [5.2, 10.1, 2.0, 22.2, 19.1, 1.6]
    b = [6, 1, 5, 2, 1, 6, 3, 4, 5, 2]
    best_merge_for_each_block = [5, 1, 6, 5, 3, 2]
    num_blocks = 6
    num_blocks_to_merge = 3
    # 6 is merged with 2, 3 is merged with 6 (2), 1 is merged with 5
    # 2, 4, 5 -> 1, 2, 3
    # 1, 2, 3, 4, 5, 6 -> 3, 1, 1, 2, 3, 1
    new_b = carry_out_best_merges(
        delta_entropy_for_each_block, best_merge_for_each_block,
        b, num_blocks, num_blocks_to_merge
    )
    @test new_b == [1, 3, 3, 1, 3, 1, 1, 2, 3, 1]
end

@testset "Agglomerative step" begin
    for T in (Array{Int64, 2}, InterblockEdgeCountDictDict, InterblockEdgeCountVectorDict, Stinger)
        test_compute_new_matrix_agglomerative(T)
        num_nodes = 50
        g = load_graph(50)
        block_partition = collect(1:num_nodes)
        M = initialize_edge_counts(
            T, g, num_nodes, block_partition
        )
        d_out, d_in, d = compute_block_degrees(M, num_nodes)
        current_block = 1
        num_blocks = num_nodes
        new_num_blocks = floor(Int64, 0.5 * num_nodes)
        block_neighbors, k_out, k_in, k = compute_block_neighbors_and_degrees(
            M, current_block
        )

        @test Set(all_neighbors(g, current_block)) == Set(block_neighbors)
        @test k_out == LightGraphs.outdegree(g, current_block)
        @test k_in == LightGraphs.indegree(g, current_block)
        @test k == LightGraphs.degree(g, current_block)

        srand(42) #Seed the RNG

        proposal = propose_new_partition_agg(
            M, current_block, block_partition, num_blocks,
            d, block_neighbors
        )
        @test proposal == 22

        # Results computed with 2 and cross checked with Python implementation.
        proposal = 2

        Δ = evaluate_proposal_agg(
            M, current_block, proposal, num_blocks, d, d_in, d_out,
            k, k_in, k_out
        )
        @test round(Δ, 8) == round(12.922923932215921, 8)

        dels = [
            11.526476498348591,
            18.399940230473767,
            10.721552778631198,
            12.676102415605499,
            14.51057342884269,
            17.370500395409749,
            11.75383351662429,
            11.325121362841692,
            9.3352584175113105,
            10.140182137228694,
            12.375461451499252
        ]

        for (del, n) in zip(dels, sort(collect(all_neighbors(g, current_block))))
            Δ = evaluate_proposal_agg(
                M, current_block, n, num_blocks, d, d_in, d_out,
                k, k_in, k_out
            )
            @test round(Δ, 8) == round(del, 8)
        end
    end
    test_carry_out_best_merges()
end

@testset "Nodal step" begin
    for T in (Array{Int64, 2}, InterblockEdgeCountDictDict, InterblockEdgeCountVectorDict, Stinger)
        test_compute_new_matrix(T)
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
