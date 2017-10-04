"""
Initializes the edge count matrix M between the blocks.
Calculates the new out, in and total degrees for the updated edge count matrix.
Returns a tuple of M, d_out, d_in, d
"""
function initalize_edge_counts(g::DiGraph, B::Int64, b::Vector{Int64})
    M = zeros(Int64, B, B) # create a zero matrix of B x B
    for v in 1:nv(g)
        for n in out_neighbors(g, v)
            # Increment count by 1
            # NOTE: Column major instead of row major
            info("Incrementing M at ($(b[n]), $(b[v]) )")
            M[b[n], b[v]] += 1
        end
    end
    # Sum across rows to get the outdegrees for each block
    d_out = reshape(sum(M, 1), B)
    # Sum across cols to get the indegrees for each block
    d_in = reshape(sum(M, 2), B)
    d = d_out .+ d_in
    return M, d_out, d_in, d
end

"""
## Propose a new block assignment for the current node or block

### Parameters
    r : Int64
            current block assignment for the node under consideration
    neighbors_out : Array{Int64, 2}, has 2 columns.
            out neighbors for the block
    neighbors_in : Array{Int64, 2}, has 2 columns.
            in neighbors for the block
    b : Vector{Int64}
        array of block assignment for each node
    M : Array{Int64, 2}, size is (B, B)
            edge count matrix between all the blocks.
    d : Vector{Int}
            total number of edges to and from each block
    B : Int64
            total number of blocks
    agg_move : Bool
            whether the proposal is a block move

### Returns
    s : int
            proposed block assignment for the node under consideration
    k_out : int
            the out degree of the node
    k_in : int
            the in degree of the node
    k : int
            the total degree of the node

### Notes
- $d_u$: degree of block u

Randomly select a neighbor of the current node, and obtain its block assignment $u$. With probability $\frac{B}{d_u + B}$, randomly propose
a block. Otherwise, randomly selects a neighbor to block $u$ and propose its block assignment. For block (agglomerative) moves,
avoid proposing the current block.
"""
function propose_new_partition(
        r::Int64, neighbors_out::Array{Int64, 2}, neighbors_in::Array{Int64, 2}, b::Vector{Int64}, M::Array{Int64, 2},
        d::Vector{Int64}, B::Int64, agg_move::Bool
    )
    neighbors = vcat(neighbors_out, neighbors_in)
    k_out = sum(neighbors_out[:, 1])
    k_in = sum(neighbors_in[:, 1])
    k = k_out + k_in
    rand_neighbor = sample(neighbors[:,0], pweights(neighbors[:,1]/k))
    u = b[rand_neighbor]
    # propose a new block randomly
    if rand() <= B/(d[u] + B)  # chance inversely prop. to block_degree
        if agg_move  # force proposal to be different from current block
            candidates = set(1:B)
            pop!(candidates, r)
            s = sample(collect(candidates))
        else
            s = rand(1:B)
        end
    else  # propose by random draw from neighbors of block partition[rand_neighbor]
        multinomial_prob = M[:, u] + M[u, :] / d[u]
        if agg_move  # force proposal to be different from current block
            multinomial_prob[r] = 0
            if sum(multinomial_prob) == 0  # the current block has no neighbors. randomly propose a different block
                candidates = set(1:B)
                pop!(candidates, r)
                s = sample(collect(candidates))
                return s, k_out, k_in, k
            else
                multinomial_prob = multinomial_prob / sum(multinomial_prob)
            end
        candidates_vec =  findn(multinomial_prob)
        s = candidates_vec[findn(rand(Multinomial(1, multinomial_prob[candidates_vec])))[1]]
        end
    end
    return s, k_out, k_in, k
end
