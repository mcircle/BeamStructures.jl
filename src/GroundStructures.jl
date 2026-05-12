struct GroundStructure{So,Se}
    Solver::So
    SensAlg::Se
    kwargs
    function GroundStructure(solver::So,sensalg::Se,kwargs) where{So,Se}
        new{So,Se}(solver,sensalg,kwargs)
    end 
end 

function GroundStructure(;solver = Tsit5(),sensalg = ForwardDiffSensitivity(),kwargs...)
    GroundStructure(solver,sensalg,kwargs)
end 



function (str::GroundStructure)(x::AbstractMatrix{T},beamtpl::NamedTuple,nodetpl::NamedTuple,adj::AbstractMatrix{TA},plt::Bool) where{T,TA}
    # x_ = reshape(x,3,:)
    
    nodepos = getstartnodes(adj)
    cbeams = size(adj,1) * (size(adj,1)-1) ÷ 2
    movables = count(x->canchangeposition(x),values(nodetpl))
    xpos = @view x[:,1:movables]
    xforces = @view x[:,movables+1:end]
    nodes_ = addpositions(nodetpl,xpos)
    
    prob_func = make_prob_func(beamtpl, nodes_, xforces, nodepos)

    ensprob  =  EnsembleProblem(prob;
                prob_func = prob_func,
                output_func = output_func(plt),
                reduction = reduction_funcF!(plt),
                u_init = initialize_u(T,beamtpl,plt)
                )

    sol = solve(ensprob,str.Solver,
                EnsembleThreads(),
                reltol = 1e-6,abstol = 1e-6,
                save_start = true,save_on = plt,save_end = true,
                sensealg=str.SensAlg,
                trajectories = cbeams
                )

    output(sol,beamtpl,nodes_,Val(plt))
end



function reduction_func_admittance!(u,data,I,solfw,adj,beams)
    for (sol,id,d0,beam) in zip(eachslice(solfw,dims = 3),I,data,beams)
        isapprox(adj[id],0) && continue
        l = beam.l
        d = d0 .* adj[id] 

        x,y = sol[3:4,2] .- sol[3:4,1]
        # 1 = M, 2 = Fx , 3 = Fy
        i = id[1]
        j = id[2]
        i_ = 3 * (i-1)+1:3*i 
        j_ = 3 * (j-1)+1:3*j
        u[i_,i_] .+= d
        u[j_,j_] .+= d
        d[1,:] .= -(d[1,:] .+ y .* d[2,:] .+ x .* d[3,:])
        
        u[i_,j_] .-= d 
        u[j_,i_] .-= d' 
    end 
    for (idx,val) in enumerate(diag(u))
        iszero(val) || continue
        u[idx,idx] =  eps(eltype(val)) # avoid singularity
    end
    u,false
end 

function output_func_admittance(bsol::ODESolution{T,N,Q},i,fsol) where{T,N,Q}
    (Array(bsol)[[1,5,6],:,2] ,false) #.* [4/9,x,y]
end

function getidxs(sz)
    len = reduce(+,1:sz-1)
    idxs = Vector{NTuple{2,Int}}(undef,len)
    i = 1
    j = 1
    for id in eachindex(idxs)
        j += i == sz ? 1 : 0
        i += i == sz ? -(sz-1-j) : 1   
        idxs[id] = (i,j)      
    end
    idxs
end 

sigmoid(x) = 1 ./ (1 .+ exp.(-x))


function make_prob_func_admittance(sol::EnsembleSolution{T},adj,beams) where{T}
    (prob,i,repeat) -> begin 
        u0 = zeros(T,14,3)
        u0[[2],1] .= one(T) #du/dθ
        u0[[3],2] .= one(T) #du/dx
        u0[[4],3] .= one(T) #du/dy
        u0[end-6:end,:] .= sol
        remake(prob;u0 = u0,p = SciMLBase.NullParameters(),tspan = (one(T),zero(T)))
    end
end

function admittance_matrix(sol::AbstractArray{T,N},adj,str,beams) where{T,N}
    # adj = sigmoid(adj)
    idxs = getindices(size(adj,1))
    lensol = length(beams)
    
    prob_func = make_vjp_func(sol,beams)
    #dadj Gradient der Steifigkeitsmatrix berechnen  
    ensprob2  =  EnsembleProblem(vjpprob;prob_func = prob_func,
                                output_func = (bsol,i) ->  output_func_admittance(bsol,i,sol),
                                reduction = (u,data,I) -> reduction_func_admittance!(u,data,idxs[I],sol,adj,beams),
                                u_init = zeros(eltype(adj),3 .* size(adj)...),
                                )
    sol = solve(ensprob2,str.Solver,save_on = false,save_start = true,save_end = true,#save_idxs = ad,#[2,3,4,9,10,11,16,17,18,22,23,24],#
                EnsembleThreads(),
                reltol = 1e-6,abstol = 1e-6,
                trajectories = lensol
                )
    
    return sol.u
end 

function reduction_func_admittance_!((u,data_out),data,I,solfw,adj,beams)
    reduction_func_admittance!(u,data,I,solfw,adj,beams)
    [u,append!(data_out,data)],false
end 
# AV being the indices of the non-moveable nodes, 
# [1,2,3] => [ϕ,x,y] of node 1 for example 
# AB being the indices of the moveable nodes
function effective_stiffness(k::AbstractMatrix{T},u::Pair{AV,AB}) where{T,AV,AB}
    nomoveidcs = [u[1]...]
    moveables = [u[2]...]
    ids = setdiff(axes(k,1),nomoveidcs)
    ke = @view k[ids,ids]
    f = zeros(T,size(k,1))
    u = zero(f)
    f[moveables] .= [1]
    u[ids] .= ke\f[ids]
    u' * k * u
end

function effective_movement(k::AbstractMatrix{T},u::Pair{AV,AB}) where{T,AV,AB}
    nomoveidcs = [u[1]...]
    moveables = [u[2]...]
    ids = setdiff(axes(k,1),nomoveidcs)    
    ke = @view k[ids,ids]
    f = zeros(T,size(k,1))
    u = zero(f)
    u[moveables] .= [1]
    f[ids] .= ke * u[ids]
    f
end

function check_structure(beams,nodes,adj)
    cbeams = reduce(+,1:size(adj,1)-1)
    @assert cbeams == length(beams) "Missing Beams! Need $(cbeams) and got only $(length(bn.Beams))"
    @assert size(adj,1) == length(nodes) "Missing Nodes! Need $(size(adj,1)) and got only $(length(nodetpl))"
end 

function (str::GroundStructure)(x::AbstractMatrix,beams,nodes,adj,plt::Bool)
    check_structure(beams,nodes,adj)
    return str(x,beams,nodes,adj,plt)
end 

(str::GroundStructure)(x::AbstractMatrix,bn::NamedTuple,adj,plt::Bool = false) = str(x,bn.Beams,bn.Nodes,adj,plt)
(str::GroundStructure)(x::AbstractMatrix,beams::NamedTuple,nodes::NamedTuple,adj,plt::Bool = false) = str(x,beams,nodes,adj,plt)

function (str::GroundStructure)(residuals::T,values::T,bn::NamedTuple,adj) where{T} #new loss
    sols,bn_ = str(values,bn,adj)
    # sols_ = toArray(sols)
    residuals!(residuals,adj,sols,bn_)
end 

function reduceposat(node::Boundary,beams::NamedTuple,y::AbstractArray{T,N},factors::AbstractVector{TF},beamnbrs) where{T,TF,N}
    # @show beamnbrs
    res = map((ind)-> factors[ind] .* ([node.ϕ,node.x,node.y] .- scalepos(beams[ind],y[2:4,2,ind],Val(2))),beamnbrs[1])
    if isempty(res)
        return Vector{T}()
    end
    reduce(hcat,res)
end
#at Clamp, force is equal to surounding
reduceforceat(node::Clamp,Beams::NamedTuple,y::AbstractArray{T,3},facs::AbstractVector{TF},beamsidxs) where {T,TF} = Vector{T}()

function reduceforceat(node::Boundary,Beams::NamedTuple,y::AbstractArray{T,3},factors::AbstractVector{TF},beamsidxs) where{T,TF}
    # @show beamsidxs
    solp = reduce((init,beampos)->init .+ factors[beampos] * scaleforce(Beams[beampos],y[[1,5,6],2,beampos]),beamsidxs[1];init = zeros(T,3))
    solm = reduce((init,beampos)->init .+ factors[beampos] * scaleforce(Beams[beampos],y[[1,5,6],1,beampos]),beamsidxs[2];init = zeros(T,3)) 
    return solp .- solm #.+ node[[6,4,5]]
end 

function residuals!(residuals::Matrix,adj::AbstractMatrix{TA},y::AbstractArray{T,N},beamstpl,nodestpl) where{T,TA,N}
    
    ids = getindices(size(adj,1))
    adj_ = ifelse.(adj .> 1,1,adj)
    adj_ = ifelse.(adj_ .< 0,0,adj)

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

function residuals!(residuals::Matrix,adj::AbstractMatrix{T},y::EnsembleSolution,bn) where{T}
    y = toArray(y)
    residuals!(residuals,adj,y,bn.Beams,bn.Nodes)
end

function residuals!(residuals::Matrix,adj::AbstractMatrix{T},y::EnsembleSolution,beams,nodes) where{T}
    y = toArray(y)
    residuals!(residuals,adj,y,beams,nodes)
end