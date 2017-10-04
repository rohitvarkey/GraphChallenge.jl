using GraphChallenge
using LightGraphs
using SimpleWeightedGraphs
using Base.Test

function test_load_graph(num_nodes::Int64)
    g = SimpleWeightedDiGraph(num_nodes)
    edges = 0
    for chunk_num = 1:9
        load_graph!(g, "snowballSampling", num_nodes, chunk_num)
        @test nv(g) == num_nodes
        edges += countlines(joinpath(
                GraphChallenge.INPUT_PATH, "snowballSampling", "$(num_nodes)_nodes",
                "simulated_blockmodel_graph_$(num_nodes)_nodes_snowball_$(chunk_num).tsv"))
        @test ne(g) == edges
    end
end

test_load_graph(1000)

include("matrix_tests.jl")
