module TopologyGeneration

using Random

export candidate_edges, adjacency_matrix, topology_id, is_admissible,
       enumerate_topologies, random_node_positions

candidate_edges(n::Integer) = [(j, i) for i in 1:n-1 for j in i+1:n]

function adjacency_matrix(mask, n::Integer)
    edges = candidate_edges(n)
    length(mask) == length(edges) ||
        throw(DimensionMismatch("mask must contain one value per possible edge"))
    adjacency = falses(n, n)
    for ((j, i), active) in zip(edges, mask)
        adjacency[i, j] = adjacency[j, i] = !iszero(active)
    end
    adjacency
end

topology_id(mask) = join(Int(x != 0) for x in mask)

function component(adjacency, start)
    visited = falses(size(adjacency, 1))
    visited[start] = true
    stack = [start]
    while !isempty(stack)
        node = pop!(stack)
        for neighbour in findall(view(adjacency, :, node))
            visited[neighbour] && continue
            visited[neighbour] = true
            push!(stack, neighbour)
        end
    end
    visited
end

"""
Check graph-level admissibility.

All clamp nodes must share one component. A branch may be unused (degree zero);
an active branch must belong to the clamp component and have at least
`minimum_branch_degree` incident beams.
"""
function is_admissible(mask; n=5, clamp_nodes=(1, 2, 5),
                       branch_nodes=(3, 4), minimum_branch_degree=2)
    sort!(collect((clamp_nodes..., branch_nodes...))) == collect(1:n) ||
        throw(ArgumentError("clamp_nodes and branch_nodes must partition 1:n"))
    adjacency = adjacency_matrix(mask, n)
    degrees = vec(sum(adjacency; dims=1))
    visited = component(adjacency, first(clamp_nodes))
    all(visited[collect(clamp_nodes)]) || return false
    for node in branch_nodes
        degrees[node] == 0 && continue
        degrees[node] >= minimum_branch_degree || return false
        visited[node] || return false
    end
    true
end

"""
Enumerate all graph-level admissible binary topologies in deterministic order.
For n=5 this examines all 2^10 = 1024 masks.
"""
function enumerate_topologies(; n=5, clamp_nodes=(1, 2, 5),
                              branch_nodes=(3, 4),
                              minimum_branch_degree=2)
    edges = candidate_edges(n)
    topologies = NamedTuple[]
    for bits in 0:(2^length(edges)-1)
        mask = BitVector([((bits >> (i-1)) & 1) == 1 for i in eachindex(edges)])
        is_admissible(mask; n, clamp_nodes, branch_nodes,
                      minimum_branch_degree) || continue
        adjacency = adjacency_matrix(mask, n)
        push!(topologies, (id=topology_id(mask), mask, adjacency,
                           elements=count(mask),
                           degrees=vec(sum(adjacency; dims=1))))
    end
    topologies
end

"""
    random_node_positions(rng; n=5, gridsize=100)

Draw `n` distinct integer node positions from the square grid
`0:gridsize × 0:gridsize`. Passing an explicit random-number generator makes
the geometry reproducible and lets every topology use the same geometry for a
given seed.
"""
function random_node_positions(rng::AbstractRNG; n=5, gridsize=100)
    n > 0 || throw(ArgumentError("n must be positive"))
    gridsize >= 0 || throw(ArgumentError("gridsize must be non-negative"))
    available = (gridsize + 1)^2
    n <= available || throw(ArgumentError("grid contains fewer than n points"))
    linear = randperm(rng, available)[1:n] .- 1
    positions = Matrix{Float64}(undef, 2, n)
    positions[1, :] .= linear .% (gridsize + 1)
    positions[2, :] .= linear .÷ (gridsize + 1)
    positions
end

end
