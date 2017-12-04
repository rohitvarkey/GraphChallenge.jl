using LightGraphs
using SimpleWeightedGraphs
using GraphChallenge

function test_initialize_counts(M::Array{Int64, 2}, g::SimpleWeightedDiGraph)
    for edge in edges(g)
        @test M[dst(edge), src(edge)] == 1
    end
end

function test_compute_new_matrix_agglomerative(::Type{Array{Int64, 2}})
    M = zeros(Int64, 3, 3)
    M = [8 3 5; 2 9 6; 4 12 10]
    r = 1
    s = 2
    M_r_row, M_r_col, M_s_row, M_s_col =
        compute_new_matrix_agglomerative(M, r, s, 3)
    @test M_r_row == [0, 0, 0]
    @test M_r_col == [0, 0, 0]
    @test M_s_row == [0, 22, 11]
    @test M_s_col == [0, 22, 16]
end

function test_compute_new_matrix(::Type{Array{Int64, 2}})
    M = zeros(Int64, 3, 3)
    M = [8 3 5; 2 9 6; 4 12 10]
    r = 1
    s = 2
    block_out_count_map = Dict(
        1=>1, 2=>2, 3=>3
    )
    block_in_count_map = Dict(
        1=>2, 2=>1, 3=>2
    )
    M_r_row, M_r_col, M_s_row, M_s_col =
        compute_new_matrix(M, r, s, 3, block_out_count_map, block_in_count_map, 0)
    @test M_r_row == [5, 3, 3]
    @test M_r_col == [5, 2, 1]
    @test M_s_row == [2, 12, 8]
    @test M_s_col == [3, 12, 15]

    @show d_out, d_in, d = compute_block_degrees(M, 3)
    @show overall_entropy = compute_overall_entropy(
        M, d_out, d_in, 3, 3, sum(d)
    )
end
