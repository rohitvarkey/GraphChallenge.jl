"""
Initializes the edge count matrix M between the blocks.
Calculates the new out, in and total degrees for the updated edge count matrix.
Returns a tuple of M, d_out, d_in, d
"""
function initialize_edge_counts!(
    M::Array{Int64, 2}, g::SimpleWeightedDiGraph, B::Int64, b::Vector{Int64}
    )
    @assert size(M) == (B, B)
    M = zeros(M, Int64) # create a zero matrix of B x B
    for edge in edges(g)
            M[b[dst(edge)], b[src(edge)]] += weight(edge)
    end
    # Sum across rows to get the outdegrees for each block
    d_out = reshape(sum(M, 1), B)
    # Sum across cols to get the indegrees for each block
    d_in = reshape(sum(M, 2), B)
    d = d_out + d_in
    return M, d_out, d_in, d
end

function initialize_edge_counts(
        _::Type{Array{Int64, 2}}, g::SimpleWeightedDiGraph, B::Int64,
        b::Vector{Int64}
    )
    M = zeros(Int64, B, B)
    initialize_edge_counts!(M, g, B, b)
end
