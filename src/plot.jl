using SimpleWeightedGraphs
using LightGraphs
using GraphPlot
using Colors
using Plots

function plotbeliefs(wg::SimpleWeightedDiGraph, b::AbstractVector)
    g = DiGraph(adjacency_matrix(wg))
    blocks = unique(b)
    block_colors = distinguishable_colors(length(blocks))
    block_color_map = Dict()
    for (block,col) in zip(blocks, block_colors)
        block_color_map[block] = col
    end
    nodefillc = [block_color_map[block] for block in b]
    return gplot(g,
        layout=(g)-> spring_layout(g; C=1), nodelabel=1:nv(g), nodelabelsize=1, nodefillc=nodefillc
    )
end
