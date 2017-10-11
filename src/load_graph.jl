using LightGraphs
using SimpleWeightedGraphs

"""
Load graph given the sampling_type, the number of vertices and the number of the
streaming peice to be added onto a given graph.
"""
function load_graph!(
    g::SimpleWeightedDiGraph, sampling_type::String, num_vertices::Int64,
    streaming_num::Int64=1
)
    if sampling_type == "snowballSampling"
        sample_name = "snowball"
    else
        sample_name = "edgeSample"
    end
    filename = joinpath(
            INPUT_PATH, "streaming", sampling_type, "$(num_vertices)_nodes",
            "simulated_blockmodel_graph_$(num_vertices)_nodes_$(sample_name)_$(streaming_num).tsv")
    edgePieces = readdlm(filename)
    for i = 1:size(edgePieces, 1)
        success = add_edge!(g, edgePieces[i, 1], edgePieces[i, 2], edgePieces[i, 3])
        if success == false
            throw("Error adding edges.")
        end
    end
end


function load_graph!(g::SimpleWeightedDiGraph, num_vertices::Int64)
    filename = joinpath(
            INPUT_PATH, "static",
            "simulated_blockmodel_graph_$(num_vertices)_nodes.tsv")
    edgePieces = readdlm(filename)
    for i = 1:size(edgePieces, 1)
        success = add_edge!(g, edgePieces[i, 1], edgePieces[i, 2], edgePieces[i, 3])
        if success == false
            throw("Error adding edges.")
        end
    end
end
