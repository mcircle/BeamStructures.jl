
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
    return [m,θ0 + θs,x./l,y./l,fx,fy,κ*l],SciMLBase.NullParameters()
end 

function initialize_beam(node::Boundary,beam::CurvedBeam,parameters::AbstractVector{T}) where {T}
    x,y,θ0, = node.x,node.y,node.ϕ
    l,θs = beam.l,beam.θs
    m,fx,fy = parameters
    m *= normfactor_m(beam) #am Balkenelement
    fx *= normfactor_f(beam) #am Balkenelement
    fy *= normfactor_f(beam) #am Balkenelement
    return [m,θ0 + θs,x./l,y./l,fx,fy,zero(T)],beam.κ0 .*l
end 

function initialize_beam(node::CompliantClamp,beam::Beam,parameters::AbstractVector{T}) where {T}
    x,y,θ0, = node.x,node.y,node.ϕ
    l,κ = beam.l,beam.κ0
    m = node.cz * (beam.θs .- parameters[1]) .* normfactor_m(beam)
    fx = parameters[2] * normfactor_f(beam) #am Balkenelement
    fy = parameters[2] * normfactor_f(beam) #am Balkenelement
    return [m,θ0 + beam.θs + parameters[1],x./l,y./l,fx,fy,κ*l],SciMLBase.NullParameters()
end

function initialize_beam(node::Clamp,beam::Beam,parameters::AbstractVector{T}) where {T}
    x,y,θ0 = node.x,node.y,node.ϕ
    l,θs,κ = beam.l,beam.θs,beam.κ0     
    m,fx,fy = parameters 
    m *=  normfactor_m(beam)
    fx *= normfactor_f(beam)
    fy *= normfactor_f(beam) #am Balkenelement
    return [m,θ0 + θs,x/l,y/l,fx,fy,κ*l],SciMLBase.NullParameters()
end 

function initialize_beam(node::Clamp,beam::CurvedBeam,parameters::AbstractVector{T}) where {T}
    x,y,θ0 = node.x,node.y,node.ϕ
    l,θs = beam.l,beam.θs
    m,fx,fy = parameters
    m *= normfactor_m(beam) #am Balkenelement
    fx *= normfactor_f(beam) #am Balkenelement
    fy *= normfactor_f(beam) #am Balkenelement
    return [m,θ0 + θs,x./l,y./l,fx,fy,zero(T)],beam.κ0 .*l
end 

function initialize_beam(node::Free,beam::Beam,parameters::AbstractVector{T}) where {T}
    x,y,θ0 = node.x,node.y,node.ϕ
    l,θs,κ = beam.l,beam.θs,beam.κ0     
    Δθ,Δx,Δy = parameters  #am Balkenelement
    return [zero(T),θ0 + θs + Δθ,(x + Δx)./l,(y + Δy)./l,zero(T),zero(T),κ*l],SciMLBase.NullParameters()
end 

function initialize_beam(node::ExtForces,beam::Beam,parameters::AbstractVector{T}) where {T}
    x,y,θ0 = node.x,node.y,node.ϕ
    l,θs,κ = beam.l,beam.θs,beam.κ0
    Δθ,Δx,Δy = parameters
    m = node.mz .* normfactor_m(beam) #am Balkenelement
    fx = node.fx .* normfactor_f(beam) #am Balkenelement
    fy = node.fy .* normfactor_f(beam) #am Balkenelement
    return [m,θ0 + θs + Δθ,(x + Δx)./l,(y + Δy)./l,fx,fy,κ*l],SciMLBase.NullParameters()
end 

function initialize_beam(node::Movable,beam::Beam,parameters::AbstractVector{T}) where {T}
    x,y,θ0 = node.x,node.y,node.ϕ
    l,θs,κ = beam.l,beam.θs,beam.κ0
    mz,fx,fy = node.trans .* parameters
    m = mz .* normfactor_m(beam) #am Balkenelement
    fx = fx .* normfactor_f(beam) #am Balkenelement
    fy = fy .* normfactor_f(beam) #am Balkenelement
    return [m,θ0 + θs,x./l,y./l,fx,fy,κ*l],SciMLBase.NullParameters()
end

function initialize_beam(beams::NamedTuple,nodes::NamedTuple,x,nodepos,i::Int)

    initialize_beam(nodes[nodepos[i]],beams[i],x[:,i])
end 

scalepos(beam::BeamElement,y::AbstractArray,::Val{2})  = @SVector [y[1] - beam.θe,y[2].*beam.l,y[3].*beam.l]
scalepos(beam::BeamElement,y::AbstractArray,::Val{1})  = @SVector [y[1] - beam.θs,y[2].*beam.l,y[3].*beam.l]

function reduceposat(node::Boundary,beams::NamedTuple,y::AbstractArray{T,3},beamnbrs) where{T}
    res = map(x-> [node.ϕ,node.x,node.y] .- scalepos(beams[x],y[2:4,2,x],Val(2)),beamnbrs)
    reduce(hcat,res)
end

function reduceposat!(res,node::Boundary,beams::NamedTuple,y::AbstractArray{T,3},beamnbrs) where{T}
    for (i,b) in enumerate(beamnbrs)
        res[:,i] .= [node.ϕ,node.x,node.y] .- scalepos(beams[b],y[2:4,2,b],Val(2))
    end
    return nothing
end

dir_vector(node::LinearSlider) = [cos(node.dir),sin(node.dir)]
function reduceposat(node::LinearSlider,beams::NamedTuple,y::AbstractArray{T,3},beamnbrs) where{T}
    ϕ = node.ϕ
    x,y = node.p0 .+ s .* dir_vector(node)
    res = map(x-> [ϕ,x,y] .- scalepos(beams[x],y[2:4,2,x],Val(2)),beamnbrs)
    reduce(hcat,res)
end

function reduceposat(node::Joint,beams::NamedTuple,y::AbstractArray{T,3},beamnbrs) where{T}
    ϕ = node.ϕ
    x,y = node.p0 .+ node.s .* dir_vector(node)
    
    res = map(x-> [ϕ,x,y] .- scalepos(beams[x],y[2:4,2,x],Val(2)),beamnbrs)
    reduce(hcat,res)
end

scaleforce(beam,y) =  y ./normvector(beam)

function reduceforceat(node::Boundary,Beams::NamedTuple,y::AbstractArray{T,3},beamsidxs) where{T}
    # scale the forces according to the beamdimensions and add them up
    solp = reduce((init,beampos)->init .+ scaleforce(Beams[beampos],y[[1,5,6],2,beampos]),beamsidxs[1];init = .-node[[6,4,5]])
    solm = reduce((init,beampos)->init .- scaleforce(Beams[beampos],y[[1,5,6],1,beampos]),beamsidxs[2];init = solp ) 
    return solm 
end 
function reduceforceat!(res,node::Boundary,Beams::NamedTuple,y::AbstractArray{T,3},beamsidxs) where{T}
    # scale the forces according to the beamdimensions and add them up
    res .-= node[[6,4,5]]
    for b in beamsidxs[1]
        res .+= scaleforce(Beams[b],y[[1,5,6],2,b])
    end 
    for b in beamsidxs[2]
        res .-= scaleforce(Beams[b],y[[1,5,6],1,b])
    end
    return nothing 
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
 

function residuals!(residuals,str::Structure,y::AbstractArray{T,3},beamstpl,nodestpl) where{T}
    # idcs = LinearIndices(residuals)#CartesianIndices((1:3,1:fld(length(residuals),3))))
    adj = str.AdjMat
    branches = count(x->!isa(x,Clamp),nodestpl)
    residuals_forces = @view residuals[:,1:branches]
    residuals_positions = @view residuals[:,branches+1:end]
    forces = 1
    positions = 1
    for (node,beams) in beamsatnode(adj,nodestpl,beamstpl)

        if forcesatnode(nodestpl[node])
            reduceforceat!(view(residuals_forces,:,forces),nodestpl[node],beamstpl,y,beams)
            # residuals_forces[:,forces] .= res
            forces += 1
        end 
        if !isempty(beams[2]) 
            # idxs = positions:positions + length(beams[2]) -1
            reduceposat!(view(residuals_positions,:,beams[2]),nodestpl[node],beamstpl,y,beams[2])
            # positions = length(beams[2]) + 1
        end  
    end
    residuals 
end 

function residuals!(residuals::Matrix,str::Structure,y::EnsembleSolution,bn) 
    y_ = toArray(y)
    residuals!(residuals,str,y_,bn.Beams,bn.Nodes)
end
function residuals!(residuals::Matrix,str::Structure,y::EnsembleSolution,beams,nodes) 
    y_ = toArray(y)
    residuals!(residuals,str,y_,beams,nodes)
end

addposition(node::BT,pos::AbstractVector{T}) where{T,BT<:Boundary{T}} =node + (;x = pos[1],y = pos[2],ϕ = pos[3])

function addposition(node::BT,pos::AbstractVector{T}) where{T,Tb,BT<:Boundary{Tb}} 
    Te = promote_type(T,Tb)
    type(BT)(Te.(node[1:end])...) + (;x = pos[1],y = pos[2],ϕ = pos[3])
end 

function addposition(node::Movable{Tb},pos::AbstractVector{T}) where{T,Tb}  
    Te = promote_type(T,Tb)
    pos_ = node.trans .* pos
    Movable{Te}(Te.(node[1:end-1])...,node.trans) + (;x = pos_[1],y = pos_[2],ϕ = pos_[3])
end

addposition(node::Clamp{T},pos::AbstractVector{T}) where{T} = node
addposition(node::LinearSlider,pos::AbstractVector{T}) where{T}  = LinearSlider{T}(node[1:3]...,pos[1]) 
addposition(node::Joint,pos::AbstractVector{T}) where{T}  = Joint{T}(node[1:2]...,pos[1]) 
addposition(node::Boundary,pos::Nothing)  = node


function addpositions(nodes::NamedTuple{names,T},xpos::AbstractMatrix) where{names,T}

    newnodes = Vector{Pair{Symbol,Boundary}}(undef, size(xpos,2))
    n = 0
    pullbacknodes = Vector{Symbol}(undef,size(xpos,2))
    for node in names
        if canchangeposition(nodes[node])
            n += 1
            newnode = addposition(nodes[node],xpos[:,n])
            @inbounds pullbacknodes[n] = node
            @inbounds newnodes[n] = node => newnode
        end 

    end
    merge(nodes,newnodes)
end 


function make_prob_func(beams::BT, nodes_::NT, xforces, nodepos) where {BT,NT}
    (prob, i, repeat) -> begin
         u0, p = initialize_beam(beams, nodes_, xforces, nodepos, i)
        remake(prob; u0=u0, p=p)
    end
end

function output_function_(sol,i)
    sol.u,false
end

function reduction_funcF!(u,data,I)
    for (d_out,i) in zip(data,I)
    
        @inbounds u[:,1,i] .= d_out[1]
        @inbounds u[:,2,i] .= d_out[2]
    end 
    u,false
end

initialize_u(::Type{T},beams,v) where{T} = v ?  [] : Array{T,3}(undef,7,2,length(beams))

reduction_funcF!(v::Bool) = v ? (u, data, I) -> (append!(u, data), false) : reduction_funcF!

output_func(v::Bool) = v ? (sol,i) -> (sol, false) : output_function_

output(x::SciMLBase.EnsembleSolution{T,N,uType},beams,nodes,::Val{false}) where{T,N,uType} = x.u,beams,nodes
output(x::SciMLBase.EnsembleSolution{T,N,uType},beams,nodes,::Val{true}) where{T,N,uType} = x

function (str::Structure)(x::AbstractArray{T,2},beams::NamedTuple,nodes::NamedTuple,v::Bool = false) where{T}
   
    nodepos = getstartnodes(str)
    movables = count(x->canchangeposition(x),values(nodes))
    xpos = @view x[:,1:movables]
    xforces = @view x[:,movables+1:end]
    nodes_ = addpositions(nodes,xpos)
    
    
    prob_func = make_prob_func(beams, nodes_, xforces, nodepos)

    ensprob  =  EnsembleProblem(prob;
                prob_func = prob_func,
                output_func = output_func(v),
                reduction = reduction_funcF!(v),
                u_init = initialize_u(T,beams,v)
                )

    sol = solve(ensprob,str.Solver,
                EnsembleThreads(),
                reltol = 1e-6,abstol = 1e-6,#saveat = 0.01,
                save_start = true,save_on = v,save_end = true,
                sensealg=str.SensAlg,
                trajectories = length(beams)
                )

    output(sol,beams,nodes_,Val(v))
end

(str::Structure)(x::AbstractArray{T,2},bn,plt::Bool = false) where{T} = str(x,bn.Beams,bn.Nodes,plt)
# (str::Structure)(x::AbstractMatrix,beams,nodes,plt::Bool = false) = str(x,beams,nodes,Val(plt))


function (str::Structure)(residuals::T,values::T,beams::NamedTuple,nodes::NamedTuple) where{T} #new loss
    sols,beams,nodes_ = str(values,beams,nodes,false)
    residuals!(residuals,str,sols,beams,nodes_)
    nothing
end 
(str::Structure)(residuals::T,values::T,bn) where{T} = str(residuals,values,bn[1],bn[2])

function getinitials(str::Structure,beams::NamedTuple,nodes::NamedTuple)
    initsize = length(beams) + count(x->!isa(x,Clamp),values(nodes))
    return  (3,initsize)
end 

function Base.zeros(::Type{T},str::Structure,beams,nodes) where{T}
    zeros(T,getinitials(str,beams,nodes))
end 


function Random.rand(::Type{T},str::Structure,beams,nodes) where{T}
    rand(T,getinitials(str,beams,nodes))
end 