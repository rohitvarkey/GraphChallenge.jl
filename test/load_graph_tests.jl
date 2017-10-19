using LightGraphs
using SimpleWeightedGraphs
using GraphChallenge

function test_load_streaming_graph(num_nodes::Int64)
    g = SimpleWeightedDiGraph(num_nodes)
    edges = 0
    for chunk_num = 1:9
        load_graph!(g, "emergingEdges", num_nodes, chunk_num)
        @test nv(g) == num_nodes
        edges += countlines(joinpath(
                GraphChallenge.INPUT_PATH, "streaming", "emergingEdges", "$(num_nodes)_nodes",
                "simulated_blockmodel_graph_$(num_nodes)_nodes_edgeSample_$(chunk_num).tsv"))
        @test ne(g) == edges
    end
end

function test_load_graph(num_nodes::Int64)
    g = load_graph(num_nodes)
    @test nv(g) == num_nodes
    edges = countlines(joinpath(
                GraphChallenge.INPUT_PATH, "static",
                "simulated_blockmodel_graph_$(num_nodes)_nodes.tsv"))
    @test ne(g) == edges
end


@testset "Graph Loading Tests" begin
    test_load_graph(1000)
    test_load_streaming_graph(1000)
end
