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
