
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
    m = node.cz * (beam.θs .- parameters[1]) .* normfactor_m(beam)
    fx = parameters[2] * normfactor_f(beam) #am Balkenelement
    fy = parameters[2] * normfactor_f(beam) #am Balkenelement
    return [m,θ0 + beam.θs + parameters[1],x./l,y./l,fx,fy,κ*l]
end

function initialize_beam(node::Clamp,beam::Beam,parameters::AbstractVector{T}) where {T}
    x,y,θ0 = node.x,node.y,node.ϕ
    l,θs,κ = beam.l,beam.θs,beam.κ0     
    m,fx,fy = parameters 
    m *=  normfactor_m(beam)
    fx *= normfactor_f(beam)
    fy *= normfactor_f(beam) #am Balkenelement
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

function initialize_beam(node::Movable,beam::Beam,parameters::AbstractVector{T}) where {T}
    x,y,θ0 = node.x,node.y,node.ϕ
    l,θs,κ = beam.l,beam.θs,beam.κ0
    mz,fx,fy = node.trans .* parameters
    m = mz .* normfactor_m(beam) #am Balkenelement
    fx = fx .* normfactor_f(beam) #am Balkenelement
    fy = fy .* normfactor_f(beam) #am Balkenelement
    return [m,θ0 + θs,x./l,y./l,fx,fy,κ*l]
end

function initialize_beam(beams::NamedTuple,nodes::NamedTuple,x,nodepos,i::Int)

    initialize_beam(nodes[nodepos[i]],beams[i],x[:,i])
end 

scalepos(beam::Beam,y::AbstractArray,::Val{2})  = [y[1] - beam.θe,y[2].*beam.l,y[3].*beam.l]
scalepos(beam::Beam,y::AbstractArray,::Val{1})  = [y[1] - beam.θs,y[2].*beam.l,y[3].*beam.l]

function reduceposat(node::Boundary,beams::NamedTuple,y::AbstractArray{T,3},beamnbrs) where{T}
    res = map(x-> [node.ϕ,node.x,node.y] .- scalepos(beams[x],y[2:4,2,x],Val(2)),beamnbrs[1])
    if isempty(res)
        return Vector{T}()
    end
    reduce(hcat,res)
end
dir_vector(node::LinearSlider) = [cos(node.dir),sin(node.dir)]
function reduceposat(node::LinearSlider,beams::NamedTuple,y::AbstractArray{T,3},beamnbrs) where{T}
    ϕ = node.ϕ
    x,y = node.p0 .+ s .* dir_vector(node)
    res = map(x-> [ϕ,x,y] .- scalepos(beams[x],y[2:4,2,x],Val(2)),beamnbrs[1])
    reduce(hcat,res)
end

function reduceposat(node::Joint,beams::NamedTuple,y::AbstractArray{T,3},beamnbrs) where{T}
    ϕ = node.ϕ
    x,y = node.p0 .+ node.s .* dir_vector(node)
    res = map(x-> [ϕ,x,y] .- scalepos(beams[x],y[2:4,2,x],Val(2)),beamnbrs[1])
    reduce(hcat,res)
end

scaleforce(beam,y) =  y ./normvector(beam)
function reduceforceat(node::Boundary,Beams::NamedTuple,y::AbstractArray{T,3},beamsidxs) where{T}
    solp = reduce((init,beampos)->init .+ scaleforce(Beams[beampos],y[[1,5,6],2,beampos]),beamsidxs[1];init = zeros(T,3))
    solm = reduce((init,beampos)->init .+ scaleforce(Beams[beampos],y[[1,5,6],1,beampos]),beamsidxs[2];init = zeros(T,3)) 
    return solp .- solm 
end 
function dir_matrix(node::LinearSlider{T}) where{T} 
    mat = zeros(T,3,3)
    mat[1] = one(T)
    mat[2:3,2] .= [cos(node.dir),sin(node.dir)]
    return mat
end 
function reduceforceat(node::LinearSlider{TN},Beams::NamedTuple,y::AbstractArray{T,3},beamsidxs) where{T,TN}
    solp = reduce((init,beampos)->init .+ scaleforce(Beams[beampos],y[[1,5,6],2,beampos]),beamsidxs[1];init = zeros(T,3))
    solm = reduce((init,beampos)->init .+ scaleforce(Beams[beampos],y[[1,5,6],1,beampos]),beamsidxs[2];init = zeros(T,3)) 
    return dir_matrix(node) * (solp .- solm) 
end

function reduceforceat(node::Joint{TN},Beams::NamedTuple,y::AbstractArray{T,3},beamsidxs) where{T,TN}
    solp = reduce((init,beampos)->init .+ scaleforce(Beams[beampos],y[[1,5,6],2,beampos]),beamsidxs[1];init = zeros(T,3))
    solm = reduce((init,beampos)->init .+ scaleforce(Beams[beampos],y[[1,5,6],1,beampos]),beamsidxs[2];init = zeros(T,3)) 
    return [1,0,0] .* (solp .- solm) 
end
reduceforceat(node::Clamp,Beams::NamedTuple,y::AbstractArray{T,3},beamsidxs) where {T} = Vector{T}()

function getnodeswithbeams(adj::AbstractMatrix,nodes::NamedTuple)
    nodeswithbeams = Vector{Int}()
    for ap in axes(adj,1)
        if isa(nodes[ap],Branch) 
            push!(nodeswithbeams,ap)
        elseif isa(nodes[ap],Clamp) && ap > 1 && sum(adj[1:ap,ap]) > 0 
            push!(nodeswithbeams,ap)
        end
    end 
    return nodeswithbeams
end    

function residuals!(residuals,str::Structure,y::AbstractArray{T,3},bn) where{T}
    # idcs = LinearIndices(residuals)#CartesianIndices((1:3,1:fld(length(residuals),3))))
    adj = str.AdjMat
    nodes = findall(x->!isapprox(x,0),LowerTriangular(adj))
    start = 1
    for n in getnodeswithbeams(adj,bn.Nodes)
        node = bn.Nodes[n]
        beams = findbeamsatnode(node,n,nodes)
        res = reduceforceat(node,bn.Beams,y,beams)
        if !isempty(res)
            residuals[:,start] .= res
            start += 1
        end 
        res = reduceposat(node,bn.Beams,y,beams)
        if !isempty(res)
            idxs = start:start + size(res,2) - 1
            residuals[:,idxs] .= res
            start = idxs[end] + 1
        end  
    end
    residuals 
end 

addposition(node::BT,pos::AbstractVector{T}) where{T,BT<:Boundary{T}} =node + (;x = pos[1],y = pos[2],ϕ = pos[3])

function addposition(node::BT,pos::AbstractVector{T}) where{T,Tb,BT<:Boundary{Tb}} 
    Te = promote_type(T,Tb)
    type(BT)(Te.(node[1:end])...) + (;x = pos[1],y = pos[2],ϕ = pos[3])
end 
addposition(node::Clamp,pos) = node
addposition(node::LinearSlider,pos::AbstractVector{T}) where{T}  = LinearSlider{T}(node[1:3]...,pos[1]) 
addposition(node::Joint,pos::AbstractVector{T}) where{T}  = Joint{T}(node[1:2]...,pos[1]) 
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

function toArray(x::AbstractVectorOfArray,pos = [0,1])
    vecs  = vec.(x(pos))
    Array(reshape(reduce(hcat,vecs),size(x.u[1],1),length(pos),length(x)))
end

function (str::Structure)(x::AbstractMatrix{T},bn::NamedTuple,::Val{false}) where{T}
    # @show T
    nodepos = getstartnodes(str)
    x_ = reshape(x,3,:)
    anz,nodes_ = changestartnodes(bn.Nodes,x_)
    # @show typeof(nodes_[1])
    function prob_func(prob,i,repeat) 
        
        u0 = initialize_beam(bn.Beams,nodes_,x_[:,anz+1:end],nodepos,i)
        remake(prob;u0 = u0,)
    end 

    ensprob  =  EnsembleProblem(prob;
                prob_func = prob_func,
                )

    sol = solve(ensprob,str.Solver,
                EnsembleThreads(),
                reltol = 1e-6,abstol = 1e-6,
                save_start = true,save_on = false,save_end = true,
                sensealg=str.SensAlg,
                trajectories = length(bn.Beams)
                )
              
    Array(sol),(;Beams = bn.Beams,Nodes = nodes_)
end

function (str::Structure)(x::AbstractMatrix{T},bn::NamedTuple, 
    ::Val{true}
    ) where{T}
    nodepos = getstartnodes(str)
    beams,nodes = bn
    anz,nodes_ = changestartnodes(bn.Nodes,x)
    # nodes_ = addpositions(nodes,x) 
    function prob_func(prob,i,repeat) 
        u0 = initialize_beam(beams,nodes_,x[:,anz+1:end],nodepos,i)
        remake(prob;u0 = u0)
    end 

    ensprob  =  EnsembleProblem(prob;prob_func = prob_func)

    sol = solve(ensprob,str.Solver,
                EnsembleThreads(),
                reltol = 1e-6,abstol = 1e-6,
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
    initsize = length(beams) + sum(x->isa(x,Branch),values(nodes))
    # for ap in adjpos
    #     initsize += isa(nodes[ap[2]],Branch) ? 1 : 0
    # end 

    return  (3,initsize)
end 

function Base.zeros(::Type{T},str::Structure,nb::NamedTuple) where{T}
    zeros(T,getinitials(str,nb))
end 

function Random.rand(::Type{T},str::Structure,nb::NamedTuple) where{T}
    rand(T,getinitials(str,nb))
end 