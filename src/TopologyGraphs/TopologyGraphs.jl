module TopologyGraphs

using DataStructures
using Distributions
using Logging
using Random
using SparseArrays: SparseMatrixCSC, sparse, spzeros, nnz, findnz, spdiagm, nzrange
using Test

using Graphs
using SimpleWeightedGraphs: SimpleWeightedGraphs
using SimpleWeightedGraphs: src, dst
using SimpleWeightedGraphs: edgetype, is_directed, nv, ne, vertices, edges
using SimpleWeightedGraphs: add_vertex!, add_vertices!, add_edge!, rem_vertex!, rem_edge!
using SimpleWeightedGraphs: has_vertex, has_edge, inneighbors, outneighbors
using SimpleWeightedGraphs: indegree, outdegree, degree, has_self_loops, num_self_loops

using SimpleWeightedGraphs: adjacency_matrix, laplacian_matrix, weights
using SimpleWeightedGraphs: connected_components, cartesian_product, induced_subgraph, pagerank

using SimpleWeightedGraphs: AbstractSimpleWeightedGraph, AbstractSimpleWeightedEdge
using SimpleWeightedGraphs: SimpleWeightedGraph
using SimpleWeightedGraphs: SimpleWeightedEdge, SimpleWeightedGraphEdge
using SimpleWeightedGraphs: weight, weighttype, get_weight, degree_matrix

export AbstractTopologyGraph, AbstractTopologyEdge
export TopologyGraph
export TopologyEdge, TopologyGraphEdge
export weight, weighttype, get_weight, degree_matrix

export bfs_shortestpath, BFShortestState, getpath, getdist

include("topologyedge.jl")
include("topologygraph.jl")
include("test.jl")

end # module
