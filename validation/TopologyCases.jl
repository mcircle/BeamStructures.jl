module TopologyCases

import BeamStructures as BS

export topology_cases, topology_quality, topology_adjacency,
       apply_boundary_condition, prepare_load_case, admissible

"""
Reproducible ground-structure benchmarks for method 2.

Geometry uses mm, forces N and moments Nmm. An edge `(j, i)` is oriented
from node `i` to node `j`, with `j > i`, as in `BS.getindices`.
"""
function candidate_adjacency(::Type{T}, n, edges) where {T}
    adjacency = zeros(T, n, n)
    for (j, i) in edges
        1 <= i < j <= n || throw(ArgumentError("candidate edges must satisfy 1 <= i < j <= n"))
        adjacency[i, j] = adjacency[j, i] = one(T)
    end
    adjacency
end

function beams_for(nodes, edges; height, width, youngs_modulus)
    T = typeof(nodes[1].x)
    map(edges) do (j, i)
        dx = nodes[j].x - nodes[i].x
        dy = nodes[j].y - nodes[i].y
        BS.Beam{T}(hypot(dx, dy), height, width, zero(T);
                   E=youngs_modulus, θs=atan(dy, dx))
    end
end

function benchmark(name, nodes, reference_edges, load_cases, output_node,
                   minimum_degrees, description)
    T = typeof(nodes[1].x)
    edges = Tuple.(BS.getindices(length(nodes)))
    edge_to_index = Dict(edge => i for (i, edge) in enumerate(edges))
    reference = falses(length(edges))
    for edge in reference_edges
        reference[edge_to_index[edge]] = true
    end
    beams = beams_for(nodes, edges; height=T(1), width=T(5),
                      youngs_modulus=T(2.1e5))
    (; name, nodes, beams, edges,
       adjacency=candidate_adjacency(T, length(nodes), edges),
       reference, load_cases, output_node, minimum_degrees, description,
       units=(length="mm", force="N", moment="Nmm"))
end

zero6(::Type{T}) where {T} = ntuple(_ -> zero(T), 6)

function three_beam_characteristic(::Type{T}) where {T<:AbstractFloat}
    z = zero6(T)
    nodes = BS.Boundary{T}[
        BS.Clamp{T}(z...),
        BS.Clamp{T}(T(30), T(70), z[3:end]...),
        BS.Branch{T}(T(50), T(50), z[3:end]...),
        BS.Clamp{T}(T(50), T(75), z[3:end]...)]
    edges = [(3,1), (3,2), (4,3), (4,1), (4,2), (2,1)]
    loads = [
        (name="left", kind=:prescribed_displacement, node=4,
         value=(x=T(-5), y=zero(T), phi=zero(T), fx=zero(T), fy=zero(T), mz=zero(T)),
         target=(quantity=:reaction, node=4, fx=T(2.2), fy=zero(T), mz=zero(T))),
        (name="right", kind=:prescribed_displacement, node=4,
         value=(x=T(5), y=zero(T), phi=zero(T), fx=zero(T), fy=zero(T), mz=zero(T)),
         target=(quantity=:reaction, node=4, fx=T(0.2), fy=zero(T), mz=zero(T)))]
    benchmark("three_beam_characteristic", nodes, edges[1:3], loads, 4,
              Dict(3=>3),
              "Three-beam characteristic from StructureSynthesis with three bypass edges.")
end

function cross_axis_pivot(::Type{T}) where {T<:AbstractFloat}
    z = zero6(T); h = T(30); span = T(30)
    nodes = BS.Boundary{T}[
        BS.Clamp{T}(z...),
        BS.Clamp{T}(span, zero(T), z[3:end]...),
        BS.Branch{T}(span, h, z[3:end]...),
        BS.Branch{T}(zero(T), h, z[3:end]...),
        BS.Clamp{T}(span/2, h, z[3:end]...)]
    edges = [(3,1), (4,2), (5,3), (5,4), (3,2), (4,1)]
    angles = T[-T(0.15), T(0.15)]
    loads = map(enumerate(angles)) do (i, phi)
        (name="rotation$(i)", kind=:prescribed_displacement, node=5,
         value=(x=-(h/2)*sin(phi), y=(h/2)*(cos(phi)-one(T)), phi=phi,
                fx=zero(T), fy=zero(T), mz=zero(T)),
         target=(quantity=:reaction, node=5, fx=zero(T), fy=zero(T), mz=zero(T)))
    end
    benchmark("cross_axis_pivot", nodes, edges[1:4], loads, 5,
              Dict(3=>2, 4=>2),
              "Cross-axis pivot based on two inclined flexures and two output links.")
end

function snap_through(::Type{T}) where {T<:AbstractFloat}
    z = zero6(T)
    nodes = BS.Boundary{T}[
        BS.Clamp{T}(z...),
        BS.Clamp{T}(T(0.5), zero(T), z[3:end]...),
        BS.Clamp{T}(zero(T), T(2), z[3:end]...),
        BS.Clamp{T}(T(0.5), T(2), z[3:end]...),
        BS.Branch{T}(zero(T), one(T), z[3:end]...),
        BS.Clamp{T}(T(0.5), one(T), z[3:end]...)]
    edges = [(5,1), (6,2), (5,3), (6,4), (6,5), (5,2)]
    loads = [
        (name="compression", kind=:prescribed_displacement, node=6,
         value=(x=T(-0.08), y=zero(T), phi=zero(T), fx=zero(T), fy=zero(T), mz=zero(T)),
         target=(quantity=:reaction, node=6, fx=zero(T), fy=zero(T), mz=zero(T))),
        (name="extension", kind=:prescribed_displacement, node=6,
         value=(x=T(0.08), y=zero(T), phi=zero(T), fx=zero(T), fy=zero(T), mz=zero(T)),
         target=(quantity=:reaction, node=6, fx=zero(T), fy=zero(T), mz=zero(T)))]
    benchmark("snap_through", nodes, edges[1:5], loads, 6,
              Dict(5=>3),
              "Snap-through example with one diagonal shortcut as a competing edge.")
end

function force_path(::Type{T}) where {T<:AbstractFloat}
    z = zero6(T)
    nodes = BS.Boundary{T}[
        BS.Clamp{T}(z...),
        BS.Clamp{T}(zero(T), T(40), z[3:end]...),
        BS.Branch{T}(T(45), T(20), z[3:end]...),
        BS.Branch{T}(T(90), T(20), z[3:end]...)]
    edges = [(3,1), (3,2), (4,1), (4,2), (4,3), (2,1)]
    loads = [
        (name="positive_force", kind=:applied_load, node=3,
         value=(x=zero(T), y=zero(T), phi=zero(T), fx=T(10), fy=zero(T), mz=zero(T)),
         target=(quantity=:displacement, node=4, x=T(-2), y=zero(T), phi=zero(T))),
        (name="negative_force", kind=:applied_load, node=3,
         value=(x=zero(T), y=zero(T), phi=zero(T), fx=T(-10), fy=zero(T), mz=zero(T)),
         target=(quantity=:displacement, node=4, x=T(2), y=zero(T), phi=zero(T)))]
    benchmark("force_path", nodes, edges[1:5], loads, 4,
              Dict(3=>2, 4=>2),
              "Force-controlled transmission based on StructureSynthesis force boundaries.")
end

topology_cases(::Type{T}=Float64) where {T<:AbstractFloat} =
    [three_beam_characteristic(T), cross_axis_pivot(T),
     snap_through(T), force_path(T)]

function topology_adjacency(case, mask=case.reference)
    length(mask) == length(case.edges) ||
        throw(DimensionMismatch("one mask value per candidate edge required"))
    T = eltype(case.adjacency)
    adjacency = zeros(T, size(case.adjacency))
    for ((j, i), active) in zip(case.edges, mask)
        adjacency[i, j] = adjacency[j, i] = T(active)
    end
    adjacency
end

function admissible(case, mask)
    length(mask) == length(case.edges) || return false
    active = Bool.(mask .!= 0)
    adjacency = topology_adjacency(case, active)
    degrees = vec(sum(adjacency .!= 0; dims=1))
    degrees[case.output_node] > 0 || return false
    all(degrees[node] >= degree for (node, degree) in case.minimum_degrees) ||
        return false
    clamps = findall(node -> node isa BS.Clamp, case.nodes)
    visited = falses(length(case.nodes))
    stack = copy(clamps)
    visited[clamps] .= true
    while !isempty(stack)
        node = pop!(stack)
        for neighbour in findall(x -> !iszero(x), view(adjacency, :, node))
            visited[neighbour] && continue
            visited[neighbour] = true
            push!(stack, neighbour)
        end
    end
    visited[case.output_node] &&
        all(i -> !(case.nodes[i] isa BS.Branch) || visited[i], eachindex(case.nodes))
end

function topology_quality(case, mask)
    length(mask) == length(case.reference) ||
        throw(DimensionMismatch("one decision per candidate edge required"))
    predicted = Bool.(mask .!= 0)
    reference = case.reference
    tp = count(predicted .& reference)
    fp = count(predicted .& .!reference)
    fn = count(.!predicted .& reference)
    tn = count(.!predicted .& .!reference)
    precision = tp + fp == 0 ? missing : tp/(tp+fp)
    recall = tp + fn == 0 ? missing : tp/(tp+fn)
    f1 = ismissing(precision) || ismissing(recall) || precision + recall == 0 ?
         missing : 2precision*recall/(precision+recall)
    (; tp, fp, fn, tn, precision, recall, f1,
       hamming=(fp+fn)/length(reference),
       admissible=admissible(case, predicted))
end

function apply_boundary_condition(case, load_case)
    nodes = copy(case.nodes)
    node = nodes[load_case.node]
    value = load_case.value
    if load_case.kind === :prescribed_displacement
        node isa BS.Clamp ||
            throw(ArgumentError("prescribed displacement requires a Clamp"))
        nodes[load_case.node] = BS.Clamp(node.x + value.x, node.y + value.y,
            node.ϕ + value.phi, node.fx + value.fx, node.fy + value.fy,
            node.mz + value.mz)
    elseif load_case.kind === :applied_load
        node isa BS.Branch ||
            throw(ArgumentError("applied load requires a Branch"))
        nodes[load_case.node] = BS.Branch(node.x + value.x, node.y + value.y,
            node.ϕ + value.phi, node.fx + value.fx, node.fy + value.fy,
            node.mz + value.mz)
    else
        throw(ArgumentError("unknown boundary-condition kind $(load_case.kind)"))
    end
    nodes
end

function prepare_load_case(case, load_case)
    nodes = apply_boundary_condition(case, load_case)
    BS.prepare(nodes..., case.beams...)
end

end
