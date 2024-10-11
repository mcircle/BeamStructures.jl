
struct Structure{A<:AbstractMatrix,So,Se,KW}
    AdjMat::A
    Solver::So
    SensAlg::Se
    kwargs::KW
    function Structure(adj::A,solver::So,sensealg::Se,kwargs::KW) where{A,So,Se,KW}
        new{A,So,Se,KW}(adj,solver,sensealg,kwargs)
    end 
end

function Structure(adj::AbstractMatrix{T};solver = Tsit5(),sensealg = ForwardDiffSensitivity(),kwargs...) where{T}
    pos = vec(any(x->!isapprox(x,zero(T)),adj,dims = 1))
    Structure(adj[pos,pos],solver,sensealg,kwargs)
end 

function Structure(adj::Connections;solver = Tsit5(),sensalg = ForwardDiffSensitivity(),kwargs...) 
    Structure(adj,solver,sensalg,kwargs)
end 

function Base.show(io::IO,str::Structure)
    return println(io, "Structure with $(size(str.AdjMat,1)) Nodes and $(Int(sum(str.AdjMat)/2)) Beam(s).")
end 

function initialize_beam(node::Boundary,beam::Beam,parameters::AbstractVector{T}) where {T}
    x,y,θ0, = node.x,node.y,node.ϕ
    l,θs,κ = beam.l,beam.θs,beam.κ0
    m,fx,fy = parameters
    m *= normfactor_m(beam) #am Balkenelement
    fx *= normfactor_f(beam) #am Balkenelement
    fy *= normfactor_f(beam) #am Balkenelement
    return [m,θ0 + θs,x./l,y./l,fx,fy,κ*l]
end 

function initialize_beam(node::CompliantClamp,beam::Beam,parameters::AbstractVector{T}) where {T}
    x,y,θ0, = node.x,node.y,node.ϕ
    l,κ = beam.l,beam.κ0
    m = node.c * (beam.θs .- parameters[1]) .* normfactor_m(beam)
    fx = parameters[2] * normfactor_f(beam) #am Balkenelement
    fy = parameters[2] * normfactor_f(beam) #am Balkenelement
    return [m,θ0 + beam.θs + parameters[1],x./l,y./l,fx,fy,κ*l]
end

function initialize_beam(node::Clamp,beam::Beam,parameters::AbstractVector{T}) where {T}
    x,y,θ0 = node.x,node.y,node.ϕ
    l,θs,κ = beam.l,beam.θs,beam.κ0     
    m,fx,fy = parameters .* [normfactor_m(beam),normfactor_f(beam),normfactor_f(beam)] #am Balkenelement
    return [m,θ0 + θs,x/l,y/l,fx,fy,κ*l]
end 

function initialize_beam(node::Free,beam::Beam,parameters::AbstractVector{T}) where {T}
    x,y,θ0 = node.x,node.y,node.ϕ
    l,θs,κ = beam.l,beam.θs,beam.κ0     
    Δθ,Δx,Δy = parameters  #am Balkenelement
    return [zero(T),θ0 + θs + Δθ,(x + Δx)./l,(y + Δy)./l,zero(T),zero(T),κ*l]
end 

function initialize_beam(node::ExtForces,beam::Beam,parameters::AbstractVector{T}) where {T}
    x,y,θ0 = node.x,node.y,node.ϕ
    l,θs,κ = beam.l,beam.θs,beam.κ0
    Δθ,Δx,Δy = parameters
    m = node.mz .* normfactor_m(beam) #am Balkenelement
    fx = node.fx .* normfactor_f(beam) #am Balkenelement
    fy = node.fy .* normfactor_f(beam) #am Balkenelement
    return [m,θ0 + θs + Δθ,(x + Δx)./l,(y + Δy)./l,fx,fy,κ*l]
end 

function initialize_beam(beams::NamedTuple,nodes::NamedTuple,x,nodepos,i::Int)
    initialize_beam(nodes[nodepos[i]],beams[i],x[:,i])
end 

scalepos(beam::Beam,y::AbstractArray,::Val{2})  = [y[1] - beam.θe,y[2].*beam.l,y[3].*beam.l]
scalepos(beam::Beam,y::AbstractArray,::Val{1})  = [y[1] - beam.θs,y[2].*beam.l,y[3].*beam.l]

function reduceposat(node::Boundary,beams::NamedTuple,y::AbstractArray{T,3},beamnbrs) where{T}
    res = map(x-> [node.ϕ,node.x,node.y] .- scalepos(beams[x],y[2:4,2,x],Val(2)),beamnbrs[1])
    reduce(hcat,res)
end

scaleforce(beam,y) =  y ./normvector(beam)

function reduceforceat(node::Boundary,Beams::NamedTuple,y::AbstractArray{T,3},beamsidxs) where{T}
    solp = reduce((init,beampos)->init .+ scaleforce(Beams[beampos],y[[1,5,6],2,beampos]),beamsidxs[1];init = zeros(T,3))
    solm = reduce((init,beampos)->init .+ scaleforce(Beams[beampos],y[[1,5,6],1,beampos]),beamsidxs[2];init = zeros(T,3)) 
    return solp .- solm .+ node[[6,4,5]]
end 

reduceforceat(node::Clamp,Beams::NamedTuple,y::AbstractArray{T,3},beamsidxs) where {T} = Vector{T}()

function residuals!(residuals,str::Structure,y::AbstractArray{T,3},bn) where{T}
    idcs = LinearIndices(CartesianIndices((1:3,1:fld(length(residuals),3))))
    adj = str.AdjMat
    nodes = findall(x->!isapprox(x,0),LowerTriangular(adj))
    start = 1
    for n in unique(first.(Tuple.(nodes)))
        node = bn.Nodes[n]
        beams = findbeamsatnode(node,n,nodes)
        res = reduceforceat(node,bn.Beams,y,beams)
        if !isempty(res)
            residuals[idcs[:,start]] .= res
            start += 1
        end 
        res = reduceposat(node,bn.Beams,y,beams)
        if !isempty(res)
            idxs = start:start + size(res,2) - 1
            residuals[idcs[:,idxs]] .= res
            start = idxs[end] + 1
        end  
    end
    residuals 
end 

addposition(node::BT,pos::AbstractVector{T}) where{T,BT<:Boundary{T}} =node + BT(pos...,zeros(T,3)...)
addposition(node::BT,pos::AbstractVector{T}) where{T,Tb,BT<:Boundary{Tb}} =gettype(node)(T.(node[1:6])...)+ gettype(node)(pos...,zeros(T,3)...)

addposition(node::Boundary,pos::Nothing)  = node

function addpositions(nodes::NamedTuple,x::AbstractMatrix)
    branches = 1:length(nodes) |> filter((x)->!isa(values(nodes[x]),Clamp))
    newnodes = map((node,pos)->keys(nodes)[node] => addposition(nodes[node],pos),branches,eachcol(x[:,1:length(branches)]))
    merge(nodes,newnodes)
end 

function changestartnodes(nodes,x)
    anz = length(nodes) - count(x->isa(x,Clamp),values(nodes))
    nodes_ = addpositions(nodes,x)
    return anz, nodes_
end 
function (str::Structure)(x::AbstractMatrix{T},bn::NamedTuple,::Val{false}) where{T}
    nodepos = getstartnodes(str.AdjMat)

    anz,nodes_ = changestartnodes(bn.Nodes,x)
    function prob_func(prob,i,repeat) 
        u0 = initialize_beam(bn.Beams,nodes_,x[:,anz+1:end],nodepos,i)
        remake(prob;u0 = u0)
    end 

    ensprob  =  EnsembleProblem(prob;
                prob_func = prob_func,
                )

    sol = solve(ensprob,str.Solver,
                EnsembleThreads(),
                reltol = 1e-12,abstol = 1e-12,
                save_start = true,save_on = false,save_end = true,
                sensealg=str.SensAlg,
                trajectories = length(bn.Beams)
                )
    Array(sol),(;Beams = bn.Beams,Nodes = nodes_)
end

function (str::Structure)(x::AbstractMatrix{T},bn::NamedTuple, 
    ::Val{true}
    ) where{T}
    nodepos = getstartnodes(str.AdjMat)
    beams,nodes = bn
    anz,nodes_ = changestartnodes(bn.Nodes,x)
    nodes_ = addpositions(nodes,x) 
    function prob_func(prob,i,repeat) 
        u0 = initialize_beam(beams,nodes_,x[:,anz+1:end],nodepos,i)
        remake(prob;u0 = u0)
    end 

    ensprob  =  EnsembleProblem(prob;prob_func = prob_func)

    sol = solve(ensprob,str.Solver,
                EnsembleThreads(),
                reltol = 1e-12,abstol = 1e-12,
                sensealg=str.SensAlg,
                trajectories = length(bn.Beams);str.kwargs...
                )
end
(str::Structure)(x::AbstractMatrix,bn,plt::Bool = false) = str(x,bn,Val(plt))
(str::Structure)(x::AbstractVector,bn,plt::Bool = false) = str(reshape(x,3,:),bn,Val(plt))


function (str::Structure)(residuals::T,values::T,bn::NamedTuple) where{T} #new loss
    sols,bn_ = str(values,bn)
    residuals!(residuals,str,sols,bn_)
end 



function getinitials(str::Structure,nb::NamedTuple)
    nodes = nb.Nodes
    beams = nb.Beams
    adjpos = findall(x->!isapprox(x,0),LowerTriangular(str.AdjMat))
    initsize = 0
    for ap in adjpos
        initsize += isa(nodes[ap[2]],Branch) ? 2 : 1
    end 
    return  (3,initsize)
end 

function Base.zeros(::Type{T},str::Structure,nb::NamedTuple) where{T}
    zeros(T,getinitials(str,nb))
end 

function Random.rand(::Type{T},str::Structure,nb::NamedTuple) where{T}
    rand(T,getinitials(str,nb))
end 
# function Base.zeros(::Type{T},strs::Vector{Structure}) where{T}
#     inits = Vector{Vector{T}}(undef,length(strs))
#     for ind in eachindex(inits)
#         inits[ind] = zeros(T,getinitials(strs[ind]))
#     end 
#     return inits
# end 

# function Random.rand(::Type{T},strs::Vector{Structure}) where{T}
#     inits = Vector{Vector{T}}(undef,length(strs))
#     for ind in eachindex(inits)
#         inits[ind] = rand(T,getinitials(strs[ind]))
#     end 
#     return inits
# end 