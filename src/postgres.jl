using NamedTuples, DataStreams, LibPQ

export InterblockEdgeCountPostgres

const edgelist_tuple = @NT(block_num::Array{Int32}, src_block::Array{Int32}, dst_block::Array{Int32}, edgecount::Array{Int32})

immutable InterblockEdgeCountPostgres
    conn::Connection

    function InterblockEdgeCountPostgres(;dbname::String="postgres")
        conn = Connection("dbname=$dbname")
        edgelist_create_stmt = """CREATE TABLE edgelist (
            block_num   integer NOT NULL,
            src_block   integer NOT NULL,
            dst_block   integer NOT NULL,
            edgecount   integer NOT NULL
        );
        """
        blockmap_create_stmt = """CREATE TABLE blockmap (
            block_num   integer NOT NULL,
            vertexid   integer NOT NULL,
            block   integer NOT NULL
        );
        """
        execute(conn, "DROP TABLE IF EXISTS edgelist, blockmap;")
        execute(conn, edgelist_create_stmt)
        execute(conn, blockmap_create_stmt)

        new(conn)
    end
end



function initialize_edge_counts(
    M::InterblockEdgeCountPostgres, g::SimpleWeightedDiGraph, B::Int64,
    b::Vector{Int64}, count_log::CountLog
    )

    block_edge_counts = Dict{Int64, Dict{Int64, Int64}}()

    total_block_edges = 0
    # Calculate edge block counts.
    for edge in edges(g)
        s, d = b[src(edge)], b[dst(edge)]
        if s in keys(block_edge_counts)
            if !(d in keys(block_edge_counts[s]))
                total_block_edges += 1
            end
            block_edge_counts[s][d] = get(block_edge_counts[s], d, 0) + 1

        else
            block_edge_counts[s] = Dict(d=>1)
            total_block_edges += 1
        end
    end

    block_num = fill(Int32(B), total_block_edges)
    src_block = ones(Int32, total_block_edges)
    dst_block = ones(Int32, total_block_edges)
    edgecounts = ones(Int32, total_block_edges)
    counter = 0
    # Create NamedTuple with the above data.
    for (src, edges) in block_edge_counts
        for (dst, edgecount) in block_edge_counts[src]
            counter += 1
            src_block[counter] = src
            dst_block[counter] = dst
            edgecounts[counter] = edgecount
        end
    end

    data = edgelist_tuple(block_num, src_block, dst_block, edgecounts)

    stmt = Data.stream!(
        data,
        Statement,
        M.conn,
        "INSERT INTO edgelist VALUES (\$1, \$2, \$3, \$4);",
    )
end
