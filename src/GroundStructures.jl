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

function (str::GroundStructure)(x::AbstractMatrix{T},bn::NamedTuple,adj::AbstractMatrix{TA},saveat::Union{AbstractFloat,AbstractVector}) where{T,TA}
    x_ = reshape(x,3,:)
    
    nodepos = getstartnodes(adj)
    cbeams = reduce(+,1:size(adj,1)-1)
    anz,nodes_ = changestartnodes(bn.Nodes,x_)
    function prob_func(prob,i,repeat) 
        
        u0 = initialize_beam(bn.Beams,nodes_,x_[:,anz+1:end],nodepos,i)
        remake(prob;u0 = u0)
    end 

    ensprob  =  EnsembleProblem(prob;
                prob_func = prob_func,
                )

    sol = solve(ensprob,str.Solver,
                EnsembleThreads(),
                reltol = 1e-6,abstol = 1e-6,saveat = saveat,
                save_start = true,save_on = true,save_end = true,
                sensealg=str.SensAlg,
                trajectories = cbeams
                )
    sol,(;Beams = bn.Beams,Nodes = nodes_)
end

function EIz(beam::Beam{T}) where{T}
    Iz = beam.w * beam.h^3/12
    EIz =beam.E * Iz
end

function reduction_func_admittance!(u,data,I,solfw,adj,beams)
    for (sol,id,d0,beam) in zip(solfw,I,data,beams)
        isapprox(adj[id],0) && continue
        l = beam.l
        d = d0 .* adj[id] 
        # f = [4    6*l  6*l;
        #        6*l 12*l^2 12*l^2;
        #        6*l 12*l^2 12*l^2]
        # d .*=       [18 36 36; 
        #         18 12 12;
        #         18 12 12]
        for i in 1:3
            iszero(d[i,i]) || continue 
            d[i,i] += 12/beam.h^2
        end
        x,y = sol[end][3:4] .- sol[1][3:4]
        # 1 = M, 2 = Fx , 3 = Fy
        i = id[1]
        j = id[2]
        i_ = 3 * (i-1)+1:3*i 
        j_ = 3 * (j-1)+1:3*j
        u[i_,i_] .+= d #.*f
        u[j_,j_] .+= d #.*f
        d[1,:] .= (-d[1,:] .+ y .* d[2,:] .+ x .* d[3,:])
        # d[2:3,2:3] .*= -1
        u[i_,j_] .+= d #.* f 
        u[j_,i_] .+= d' #.*f'
    end 
    for (idx,val) in enumerate(diag(u))
        iszero(val) || continue
        u[idx,idx] .= one(eltype(val)) # avoid singularity
    end
    u,false
end 

function output_func_admittance(bsol::ODESolution{T,N,Q},i,fsol) where{T,N,Q}

    (bsol[end][[1,5,6],:] ,false) #.* [4/9,x,y]
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

function admittance_matrix(sol::EnsembleSolution{T,N,S},adj,str,beams) where{T,N,S}
    # adj = sigmoid(adj)
    idxs =findall(x->!isapprox(x,0),LowerTriangular(adj))
    lensol = length(sol)
    function prob_func(prob,i,repeat) 
        u0 = zeros(T,7,3)
        u0[[2],1] .= one(T) #du/dx
        u0[[3],2] .= one(T) #du/dy
        u0[[4],3] .= one(T) #du/dθ
        remake(prob;u0 = u0,p = sol[i],tspan = (one(T),zero(T)))
    end
    #dadj Gradient der Steifigkeitsmatrix berechnen  
    ensprob2  =  EnsembleProblem(vjpprob;prob_func = prob_func,
                                output_func = (bsol,i) ->  output_func_admittance(bsol,i,sol),
                                reduction = (u,data,I) -> reduction_func_admittance!(u,data,idxs[I],sol,adj,beams),
                                u_init = zeros(eltype(adj),3 .* size(adj)...),
                                )
    sol = solve(ensprob2,str.Solver,save_on = false,save_start = false,save_end = true,#save_idxs = ad,#[2,3,4,9,10,11,16,17,18,22,23,24],#
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

function check_structure(bn,adj)
    cbeams = reduce(+,1:size(adj,1)-1)
    @assert cbeams == length(bn.Beams) "Missing Beams! Need $(cbeams) and got only $(length(bn.Beams))"
    @assert size(adj,1) == length(bn.Nodes) "Missing Nodes! Need $(size(adj,1)) and got only $(length(bn.Nodes))"
end 

 function (str::GroundStructure)(x::AbstractMatrix,bn,adj,saveat,::Val{false})
    check_structure(bn,adj)
    x = str(x,bn,adj,saveat)
    sol,bn_ = x

    sol,bn_#,st
end 

function (str::GroundStructure)(x::AbstractMatrix,bn,adj,saveat,::Val{true})
    check_structure(bn,adj)
    sol,bn_ = str(x,bn,adj,saveat)
    sol
end 

(str::GroundStructure)(x::AbstractMatrix,bn,adj,plt::Bool = false,saveat::Union{Real,AbstractVector} = []) = str(x,bn,adj,saveat,Val(plt))
(str::GroundStructure)(x::AbstractVector,bn,adj,plt::Bool = false,saveat =[]) = str(reshape(x,3,:),bn,adj,saveat,Val(plt))

function (str::GroundStructure)(residuals::T,values::T,bn::NamedTuple,adj) where{T} #new loss
    sols,bn_ = str(values,bn,adj)
    # sols_ = toArray(sols)
    residuals!(residuals,adj,sols,bn_)
end 

function reduceposat(node::Boundary,beams::NamedTuple,y::AbstractArray{T,3},factors::AbstractVector{TF},beamnbrs) where{T,TF}
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

function residuals!(residuals::Matrix,adj::AbstractMatrix{TA},y::AbstractArray{T,N},bn) where{T,TA,N}
    
    ids = getindices(size(adj,1))
    adj_ = ifelse.(adj .> 1,1,adj)
    adj_ = ifelse.(adj_ .< 0,0,adj)
    branches = count(x->isa(x,Branch),bn.Nodes)
    residuals_forces = @view residuals[:,1:branches]
    residuals_positions = @view residuals[:,branches+1:end]
    # idcs = LinearIndices(residuals)
    # nodes = findall(x->!isapprox(x,0),LowerTriangular(adj))
    forces = 1
    positions = 1
    for n in axes(adj,1)
        node = bn.Nodes[n]
        beams = findbeamsatnode(node,n,ids)
        res = reduceforceat(node,bn.Beams,y,adj[ids],beams)
       
        if !isempty(res)
            residuals_forces[:,forces] .= res
            forces += 1
        end 
        res = reduceposat(node,bn.Beams,y,adj[ids],beams)
        
        if !isempty(res)
            idxs = positions:positions + size(res,2) - 1
            residuals_positions[:,idxs] .= res
            positions = idxs[end] + 1
        end  
    end
    residuals 
end 

function residuals!(residuals::Matrix,adj::AbstractMatrix{T},y::EnsembleSolution,bn) where{T}
    y = toArray(y)
    residuals!(residuals,adj,y,bn)
end