function CRC.rrule(::typeof(reducevcat),res)
    sol = reducevcat(res)
    function vcat_back(ȳ)
        if isempty(res)
            grads = NoTangent()
        else
            grads = @thunk(collect(eachcol(reshape(ȳ,3,:))))
        end 
        return NoTangent(),grads
    end 

    return sol,vcat_back
end 

function CRC.rrule(::Type{Beam},l,h,w,κ0,E,θs,θe)
    function back_beam(ȳ)
        return NoTangent(),ȳ.l,ȳ.h,ȳ.w,ȳ.κ0,ȳ.E,ȳ.θs,ȳ.θe
    end 
    return Beam(l,h,w,κ0,E,θs,θe),back_beam
end 

function CRC.rrule(::typeof(scalepos),beam::Beam{T},y,side::Val{N}) where{T,N}
    sol = scalepos(beam,y,side)
    function scalepos_back(ȳ)
        dy =  ȳ .* [1,beam.l,beam.l]
        if N == 2
            dbeam = Tangent{Beam{T}}(;l = sum(ȳ[2:3] .* y[2:3]),θe = -ȳ[1])
        else
            dbeam = Tangent{Beam{T}}(;l = sum(ȳ[2:3] .* y[2:3]),θs = -ȳ[1])
        end 
        return NoTangent(),dbeam, dy,NoTangent() 
    end 
    return sol,scalepos_back
end 

function CRC.rrule(::typeof(reduceposat),node::Bo,beams,y,beamnbrs) where{Bo<:Boundary{T}} where{T} 
    sol = reduceposat(node,beams,y,beamnbrs)
    function reduceposat_back(ȳ)
        rȳ = reshape(ȳ,3,:)
        
        ∂beams= Tangent{typeof(beams)}(;ntuple(x->keys(beams)[x]=>ZeroTangent(),length(beams))...)
        ∂y = zero(y)
        for (n,b) in enumerate(beamnbrs[1])
            _, spback = CRC.rrule(scalepos,beams[b],y[2:4,2,b],Val(2))
            _,dbeam,dy,_ = spback(rȳ[:,n])    
            ∂beams -= Tangent{typeof(beams)}(;Symbol(:Beam_,b) => dbeam)
            ∂y[2:4,2,b] .-= dy
        end
        dnode = Tangent{Bo}(;x = sum(rȳ[2,:]),y = sum(rȳ[3,:]),ϕ = sum(rȳ[1,:]))
        return NoTangent(),dnode,∂beams,∂y,NoTangent()
    end 
    return sol, reduceposat_back
end 

function CRC.rrule(::typeof(reduceposat),node::Bo,beams,y::AbstractArray{T},facs,beamnbrs) where{T,Bo<:Boundary{T}} 
    sol = reduceposat(node,beams,y,facs,beamnbrs)
    function reduceposat_back(ȳ)
        rȳ = reshape(ȳ,3,:)
        ∂beams= Tangent{typeof(beams)}(;ntuple(x->keys(beams)[x]=>ZeroTangent(),length(beams))...)
        ∂y = zero(y)
        
        ∂f = zero(facs)
        dnode = Tangent{Bo}()
        for (n,b) in enumerate(beamnbrs[1])
            sp, spback = CRC.rrule(scalepos,beams[b],y[2:4,2,b],Val(2))

            # tmp = facs[b] .* rȳ[:,n] #.*([node.ϕ,node.x,node.y] .- sp)

            _,dbeam,dy,_ = spback(rȳ[:,n])    
            
            ∂beams -= Tangent{typeof(beams)}(;Symbol(:Beam_,b) => dbeam)
            ∂y[2:4,2,b] .-= dy
            #faktor hat auch einfluss auf Anfangsknoten -> diff zwischen knoten und beam enden 
            # @show ([node.ϕ,node.x,node.y] .- sp)' * ȳ[:,n]
            ∂f[b] =   ([node.ϕ,node.x,node.y] .- sp)' * ȳ[:,n]
            # dnode = dnode + Tangent{Bo}(;x = tmp[2], y = tmp[3],ϕ = tmp[1])
            dnode = dnode + Tangent{Bo}(;x = rȳ[2,n], y = rȳ[3,n],ϕ = rȳ[1,n])
        end
        return NoTangent(),dnode,∂beams,∂y,∂f, NoTangent()
    end 
    return sol, reduceposat_back
end 

function pullback_normfactor_m(ȳ,b::Beam{T}) where{T}
    nom = 1/(b.E * b.w * b.h^3)
    ∂b = Tangent{Beam{T}}(;
        l = 12 * ȳ  * nom,
        h = -36 * b.l * ȳ * nom / b.h,
        w = -12 * b.l * ȳ * nom / b.w ,
        E = -12 * b.l * ȳ * nom / b.E
    )
    return (NoTangent(), ∂b)
end

function CRC.rrule(::typeof(normfactor_m), b::Beam{T}) where{T}
    y = normfactor_m(b)
    pullback_norm_m(ȳ) = pullback_normfactor_m(ȳ,b)  
    return y, pullback_norm_m
end

function pullback_normfactor_f(ȳ,b::Beam{T}) where{T}
    ∂m = b.l * ȳ
    _, pb_m = CRC.rrule(normfactor_m, b)
    _, ∂b_m = pb_m(∂m)
    ∂b = Tangent{Beam{T}}(;
        l = ȳ * normfactor_m(b) + ∂b_m.l,
        h = ∂b_m.h,
        w = ∂b_m.w,
        E = ∂b_m.E
    )
    return (NoTangent(), ∂b)
end

function CRC.rrule(::typeof(normfactor_f), b::Beam{T}) where{T}
    y = normfactor_f(b)
    pullback_norm_f(ȳ) = pullback_normfactor_f(ȳ,b)
    return y, pullback_norm_f
end

function pullback_normfactor(ȳ,b::Beam{T}) where{T}
    _,∂b_m = pullback_normfactor_m(ȳ[1],b)
    _,∂b_f1 = pullback_normfactor_f(ȳ[2],b)
    _,∂b_f2 = pullback_normfactor_f(ȳ[3],b)
    ∂beam = ∂b_m + ∂b_f1 + ∂b_f2
    return (NoTangent(), ∂beam)
end

function CRC.rrule(::typeof(normvector), b::Beam{T}) where{T}
    y = normvector(b)
    pullback_norm(ȳ) = pullback_normfactor(ȳ,b)
    return y, pullback_norm
end

function CRC.rrule(::typeof(scaleforce),b::Beam,y)
    nv , pb_nv = CRC.rrule(normvector, b)
    result =  y ./ nv
    
    function pullback(ȳ)
        ∂y = ȳ ./ nv
        ∂nv = -ȳ .* y ./ nv.^2
        
        _, ∂beam = pb_nv(∂nv)
        return (NoTangent(),∂beam,∂y)
    end
    
    return result, pullback
end

function pullback_init_beam(ȳ,node::No,beam::Beam{BT},pars) where{BT,No<:Boundary{T} where{T}}
    l,κ = beam.l,beam.κ0
    x,y = node.x,node.y

    _,∂beam = pullback_normfactor(pars[1:3] .* ȳ[[1,5,6]],beam)
    ∂l = -x/l^2 * ȳ[3] - y/l^2 * ȳ[4] + κ * ȳ[7]

    ∂beam += Tangent{Beam{BT}}(l = ∂l,κ0 = l*ȳ[7],θs = ȳ[2])

    ∂node = Tangent{No}(;x = ȳ[3]/l,y=ȳ[4]/l,ϕ = ȳ[2])

    ∂pars = normvector(beam) .* ȳ[[1,5,6]]
    return NoTangent(),∂node,∂beam,∂pars
end 

function CRC.rrule(::typeof(initialize_beam),node::No,beam::Beam{BT},pars::AbstractVector) where{BT,No<:Boundary}
    x,y,θ0, = node.x,node.y,node.ϕ
    l,θs,κ = beam.l,beam.θs,beam.κ0
    m,fx,fy = pars .* normvector(beam) #am Balkenelement
    result = [m,θ0 + θs,x ./l,y ./l,fx,fy,κ*l]
    init_beam_back(ȳ) = pullback_init_beam(ȳ,node,beam,pars)
    return result, init_beam_back
end 

function init_u0_pullback(ȳ,beams::NamedTuple,nodes::NamedTuple,pars,idxs,i)
    x_idxs = idxs[i]
       
    ∂nodes = Tangent{typeof(nodes)}(;ntuple(x-> keys(nodes)[x]=>ZeroTangent(),length(nodes))...) 
    ∂beams = Tangent{typeof(beams)}(;ntuple(x-> keys(beams)[x]=>ZeroTangent(),length(beams))...) 
    
    # ∂pars = zero(pars)
    _,dnode,dbeam,dpars = pullback_init_beam(ȳ,nodes[x_idxs],beams[i],pars[:,x_idxs])
    # ∂pars[:,i] += dpars 
    ∂beams += Tangent{typeof(beams)}(;Symbol(:Beam_,i) => dbeam)
    ∂nodes += Tangent{typeof(nodes)}(;Symbol(:Node_,x_idxs) => dnode)
    return NoTangent(), ∂beams,∂nodes,dpars,NoTangent(),NoTangent()
end 

function CRC.rrule(::typeof(initialize_beam),beams::NamedTuple,nodes::NamedTuple,pars,x_idxs::Number,i::Int)
    result = initialize_beam(nodes[x_idxs],beams[i],pars[:,x_idxs])
    init_u0_back(ȳ) = init_u0_pullback(ȳ,beams,nodes,pars,x_idxs,i)
    result,init_u0_back
end 

function CRC.rrule(::typeof(initialize_beam),beams::NamedTuple,nodes::NamedTuple,pars,idxs::AbstractVector,i::Int)
    x_idxs = idxs[i]
    result = initialize_beam(nodes[x_idxs],beams[i],pars[:,x_idxs])
    init_u0_back(ȳ) = init_u0_pullback(ȳ,beams,nodes,pars,idxs,i)
    result,init_u0_back
end 

function CRC.rrule(::typeof(reduceforceat),node::No,beams,y,idxs) where{No}
    
    sol = reduceforceat(node,beams,y,idxs)
    function force_bound_back(ȳ)

        y_idxs = getforceindices(idxs)
        ∂y = zero(y)
        y_ = @view y[y_idxs]
        ∂y_ = @view ∂y[y_idxs] 
        Tb = eltype(beams)
        ∂beams = Tangent{typeof(beams)}(;ntuple(x->keys(beams)[x]=>ZeroTangent(),length(beams))...)
        for (b,yb,dyb) in zip(vcat(idxs...),eachcol(y_),eachcol(∂y_))
            _, rr_norm = CRC.rrule(scaleforce,beams[b],yb)
            if b in idxs[2]
                _,dbeam,dyrr_ = rr_norm(-ȳ) 
            else
                _,dbeam,dyrr_ = rr_norm(ȳ)
            end 
            ∂beams += Tangent{typeof(beams)}(;Symbol(:Beam_,b) => dbeam)
            dyb .= dyrr_
        end 

        ∂node = Tangent{No}(;fx = ȳ[2],fy = ȳ[3],mz = ȳ[1])
        return NoTangent(),∂node,∂beams,∂y,NoTangent()
    end 
    return sol,force_bound_back
end 

function CRC.rrule(::typeof(reduceforceat),node::No,beams,y::AbstractArray{T,3},factors,idxs) where{T,No}
    
    sol = reduceforceat(node,beams,y,factors,idxs)
    function force_bound_back(ȳ)
        y_idxs = getforceindices(idxs)
        ∂y = zero(y)
        ∂fac = zero(factors)
        y_ = @view y[y_idxs]
        ∂y_ = @view ∂y[y_idxs] 
        Tb = eltype(beams)
        ∂beams = Tangent{typeof(beams)}(;ntuple(x->keys(beams)[x]=>ZeroTangent(),length(beams))...)
        for (b,yb,dyb) in zip(vcat(idxs...),eachcol(y_),eachcol(∂y_))
            
            f, rr_norm = CRC.rrule(scaleforce,beams[b],yb)
            if b in idxs[2]
                _,dbeam,dyrr_ = rr_norm(-ȳ) 
                ∂fac[b] = sum(f .* ȳ) 
            else
                _,dbeam,dyrr_ = rr_norm(ȳ)
                ∂fac[b] = -sum(f .* ȳ) 
            end 
            ∂beams += Tangent{typeof(beams)}(;Symbol(:Beam_,b) => dbeam)
            dyb .= dyrr_
        end 

        ∂node = Tangent{No}()#;fx = ȳ[2],fy = ȳ[3],mz = ȳ[1])
        return NoTangent(),∂node,∂beams,∂y,∂fac,NoTangent()
    end 
    return sol,force_bound_back
end 

@non_differentiable reduceforceat(n::Clamp,beams,y,idxs) 

function CRC.rrule(::typeof(residuals!),residuals,str::Structure,y,bn)
    ind = 1
    adj = str.AdjMat
    idcs = LinearIndices(residuals)
    nodes = findall(x->!isapprox(x,0),LowerTriangular(adj))
    start = 1
    resposdict = Dict{Int,AbstractRange{Int}}()
    resforcedict = Dict{Int,Int}()
    rrposdict = Dict{Int,Function}()
    rrforcedict = Dict{Int,Function}()
    for n in unique(first.(Tuple.(nodes)))
        node = bn.Nodes[n]
        beams = findbeamsatnode(node,n,nodes)
        res, pullback_reduceforceat = rrule(reduceforceat,node,bn.Beams,y,beams)
        if !isempty(res)
            residuals[idcs[:,start]] .= res
            resforcedict[n] = start
            rrforcedict[n] = pullback_reduceforceat
            start += 1
        end 
        res, pullback_reduceposat = rrule(reduceposat,node,bn.Beams,y,beams)
        if !isempty(res)
            idxs = start:start + size(res,2) - 1
            resposdict[n] = idxs
            rrposdict[n] = pullback_reduceposat
            residuals[idcs[:,idxs]] .= res
            start = idxs[end] + 1
        end  
    end
    function residuals!_back(ȳ)
        ∂res = InplaceableThunk(dself -> dself .+= ȳ,@thunk(copy(ȳ)))
        ∂y = zero(y)
        ∂beams = Tangent{typeof(bn.Beams)}(;ntuple(x->keys(bn.Beams)[x] =>ZeroTangent(),length(bn.Beams))...)        
        ∂nodes = Tangent{typeof(bn.Nodes)}(;ntuple(x->keys(bn.Nodes)[x] =>ZeroTangent(),length(bn.Nodes))...)

        for (ind,pos) in resforcedict
            _,dnode,dbeams,dy,_ =rrforcedict[ind](ȳ[idcs[:,pos]])
            ∂nodes += Tangent{typeof(bn.Nodes)}(;Symbol(:Node_,ind) => dnode)
            ∂beams += dbeams
            ∂y .+= dy
        end
        
        for (ind,pos) in resposdict
            _,dnode,dbeams,dy,_ = rrposdict[ind](ȳ[idcs[:,pos]])
            ∂nodes += Tangent{typeof(bn.Nodes)}(;Symbol(:Node_,ind) => dnode)
            ∂beams +=  dbeams
            ∂y .+= dy
        end 
        ∂bn = Tangent{typeof(bn)}(;Beams = ∂beams,Nodes = ∂nodes)
        return NoTangent(),∂res,NoTangent(),∂y,∂bn
    end 

    return residuals,residuals!_back
end 

function CRC.rrule(::typeof(residuals!),residuals::Matrix,adj::AbstractMatrix,y,bn)

    ind = 1
    idcs = LinearIndices(residuals)
    nodes = findall(x->!isapprox(x,0),LowerTriangular(adj))
    start = 1
    resposdict = Dict{Int,AbstractRange{Int}}()
    resforcedict = Dict{Int,Int}()
    rrposdict = Dict{Int,Function}()
    rrforcedict = Dict{Int,Function}()
    for n in unique(first.(Tuple.(nodes)))
        node = bn.Nodes[n]
        beams = findbeamsatnode(node,n,nodes)
        res, pullback_reduceforceat = rrule(reduceforceat,node,bn.Beams,y,adj[nodes],beams)
        if !isempty(res)
            residuals[idcs[:,start]] .= res
            resforcedict[n] = start
            rrforcedict[n] = pullback_reduceforceat
            start += 1
        end 
        res, pullback_reduceposat = rrule(reduceposat,node,bn.Beams,y,adj[nodes],beams)
        if !isempty(res)
            idxs = start:start + size(res,2) - 1
            resposdict[n] = idxs
            rrposdict[n] = pullback_reduceposat
            residuals[idcs[:,idxs]] .= res
            start = idxs[end] + 1
        end  
    end    
    function residuals!_back(ȳ)
        ∂res = InplaceableThunk(dself -> dself .+= ȳ,@thunk(copy(ȳ)))
        ∂y = zero(y)
        ∂fac = zero(adj)
        ∂beams = Tangent{typeof(bn.Beams)}(;ntuple(x->keys(bn.Beams)[x] =>ZeroTangent(),length(bn.Beams))...)        
        ∂nodes = Tangent{typeof(bn.Nodes)}(;ntuple(x->keys(bn.Nodes)[x] =>ZeroTangent(),length(bn.Nodes))...)

        for (ind,pos) in resforcedict
            _,dnode,dbeams,dy,df =rrforcedict[ind](ȳ[idcs[:,pos]])
            ∂nodes += Tangent{typeof(bn.Nodes)}(;Symbol(:Node_,ind) => dnode)
            ∂beams += dbeams
            ∂y .+= dy
            ∂fac[nodes] .+= df
        end
        
        for (ind,pos) in resposdict
            _,dnode,dbeams,dy,df = rrposdict[ind](ȳ[idcs[:,pos]])
            ∂nodes += Tangent{typeof(bn.Nodes)}(;Symbol(:Node_,ind) => dnode)
            ∂beams +=  dbeams
            ∂y .+= dy
            ∂fac[nodes] .+= df
        end 
        ∂bn = Tangent{typeof(bn)}(;Beams = ∂beams,Nodes = ∂nodes)
        return NoTangent(),∂res,∂fac .+ transpose(∂fac),∂y,∂bn
    end    


    return residuals,residuals!_back
end 

function CRC.rrule(::Type{BT}, x,y,θ,fx,fy,mz) where{BT<:Boundary{T}} where{T}
    Bar_pullback(ȳ) = (NoTangent(),ȳ.x,ȳ.y,ȳ.θ,ȳ.fx,ȳ.fy,ȳ.mz)
    return BT(x,y,θ,fx,fy,mz), Bar_pullback
end
function CRC.rrule(::Type{Clamp{T}}, x,y,θ,fx,fy,mz) where{T}
    Bar_pullback(ȳ) = (NoTangent(),ȳ.x,ȳ.y,ȳ.θ,ȳ.fx,ȳ.fy,ȳ.mz)
    return Clamp{T}(x,y,θ,fx,fy,mz), Bar_pullback
end

function CRC.rrule(::typeof(+),a::BT,b::BT) where{BT<:Boundary{T}} where{T}
    backaddb(ȳ) = (NoTangent(),ȳ,ȳ)
    return a + b,backaddb
end 

addpos_back(ȳ) = (NoTangent(),ȳ,[ȳ.x,ȳ.y,ȳ.ϕ])
function CRC.rrule(::typeof(addposition),node::BT,x) where{T,BT<:Boundary{T}}
    res = addposition(node,x)
    return res,addpos_back
end     

function pullback_addposition(ȳ,∂x,branches)
    for (colx,node) in zip(eachcol(∂x),branches)
        colx[1] = ȳ[node].x
        colx[2] = ȳ[node].y
        colx[3] = ȳ[node].ϕ
    end 
    return NoTangent(),ȳ,∂x
end 

function CRC.rrule(::typeof(addpositions),nodes,x)
    branches = 1:length(nodes) |> filter((x)->!isa(values(nodes[x]),Clamp))
    newnodes = map((node,pos)->keys(nodes)[node] => addposition(nodes[node],pos),branches,eachcol(x[:,1:length(branches)]))
    addpos_back(ȳ) = pullback_addposition(ȳ,zero(x),branches)
    return merge(nodes,newnodes),addpos_back
end

function CRC.rrule(::typeof(changestartnodes),nodes,x)
    anz = length(nodes) - count(x->isa(x,Clamp),values(nodes))
    nodes_,back = CRC.rrule(addpositions,nodes,x)
    function changeback(ȳ)
        _,∂nodes,∂x = back(ȳ[2])
        return NoTangent(),∂nodes,∂x
    end 
    return (anz,nodes_),changeback
end 


function output_func(u,∂sol,beams,nodes,nodepos,x,i)
    x_idxs = nodepos[i]
    # _,∂beam,∂node,∂pars,_,_ = init_u0_pullback(u(0) .+ ∂sol[:,1,i],beams,nodes,x,nodepos,i)
    _,dnode,dbeam,dpars = pullback_init_beam(u(0) .+ ∂sol[:,1,i],nodes[x_idxs],beams[i],x[:,x_idxs])
    ∂beams = Tangent{typeof(beams)}(;Symbol(:Beam_,i) => dbeam)
    ∂nodes = Tangent{typeof(nodes)}(;Symbol(:Node_,x_idxs) => dnode)
    (∂beams,∂nodes,dpars),false
end 

function reduction_func!(u,data,I)
    
    ∂pars,∂Beams,∂Nodes = u
    for ((∂beam,∂node,∂pars),i) in zip(data,I)
        u[3] += ∂node
        u[2] += ∂beam
        u[1][:,i] .+= ∂pars
    end 
    u,false
end 

function u_init(x, bn)
    [zero(x),
    Tangent{typeof(bn.Beams)}(;ntuple(x->keys(bn.Beams)[x] =>ZeroTangent(),length(bn.Beams))...),
    Tangent{typeof(bn.Nodes)}(;ntuple(x->keys(bn.Nodes)[x] =>ZeroTangent(),length(bn.Nodes))...)]
end 

function CRC.rrule(str::Structure,x,bn::NamedTuple)
    nodepos = getstartnodes(str)

    (anz,nodes_),change_pullback = CRC.rrule(changestartnodes,bn.Nodes,x)
    
    function prob_func(prob,i,repeat) 
        u0 = initialize_beam(bn.Beams,nodes_,x[:,anz+1:end],nodepos,i) #x(:,anz+1:end) Parameters for the beams
        remake(prob;u0 = u0)
    end 

    ensprob  =  EnsembleProblem(prob;prob_func = prob_func)

    sol = solve(ensprob,str.Solver,
                EnsembleThreads(),
                reltol = 1e-6,abstol = 1e-6,
                trajectories = length(bn.Beams)
                )

    function back_ode(ȳ)
        @inbounds ∂sol,∂bn = ȳ

        function prob_func(prob,i,repeat) 
            u0 = collect(∂sol[:,end,i])
            remake(prob;u0 = u0,p = sol[i])
        end
        ensprob  =  EnsembleProblem(vjpprob;prob_func = prob_func,
                                    output_func = (sol,i) -> output_func(sol,∂sol,bn.Beams,bn.Nodes,nodepos,x,i),
                                    reduction = (u,data,I) -> reduction_func!(u,data,anz .+ I),
                                    u_init = u_init(x,bn))
        sol = solve(ensprob,str.Solver,
                    EnsembleThreads(),
                    reltol = 1e-6,abstol = 1e-6,
                    trajectories = length(bn.Beams)
                    )
        ∂x,∂Beams,∂Nodes = sol

        _,∂nodes_change, ∂xchange = change_pullback((nothing,∂bn.Nodes))
        ∂Nodes += ∂nodes_change
        ∂x .+= ∂xchange
        ∂bn = Tangent{typeof(bn)}(;Beams = ∂Beams,Nodes = ∂Nodes)
        return NoTangent(),∂x,∂bn
    end
    (toArray(sol),(;Beams = bn.Beams,Nodes = nodes_)),back_ode
end 

function CRC.rrule(str::GroundStructure,x::AbstractMatrix,bn::NamedTuple,adj)
    
    nodepos = getstartnodes(adj)
    s = size(adj,1)
    cbeams = reduce(+,1:size(adj,1)-1)
    
    (anz,nodes_),change_pullback = CRC.rrule(changestartnodes,bn.Nodes,x)
    
    function prob_func(prob,i,repeat) 
        
        u0 = initialize_beam(bn.Beams,nodes_,x[:,anz+1:end],nodepos,i) #x(:,anz+1:end) Parameters for the beams
        remake(prob;u0 = u0)
    end 

    ensprob  =  EnsembleProblem(prob;prob_func = prob_func)

    sol = solve(ensprob,str.Solver,
                EnsembleThreads(),
                reltol = 1e-6,abstol = 1e-6,
                trajectories = cbeams
                )
    solst = admittance_matrix(sol,adj,str)

    function back_groundstr(ȳ)
        @inbounds ∂sol,∂bn,∂st = ȳ
        
        ∂Beams = ∂bn.Beams   
        ∂Nodes = ∂bn.Nodes
        function prob_func(prob,i,repeat) 
            u0 = collect(∂sol[:,end,i])
            remake(prob;u0 = u0,p = sol[i])
        end

        ensprob  =  EnsembleProblem(vjpprob;prob_func = prob_func,
                                            output_func = (sol,i) -> output_func(sol,∂sol,bn.Beams,bn.Nodes,nodepos,x,i),
                                            reduction = (u,data,I) -> reduction_func!(u,data,anz.+I),
                                            u_init = u_init(x,bn))
        solp = solve(ensprob,str.Solver,
                    EnsembleThreads(),
                    reltol = 1e-6,abstol = 1e-6,
                    trajectories = cbeams
                    )
        ∂x,∂Beams,∂Nodes = solp
        
        _,∂nodes_change, ∂xchange = change_pullback((nothing,∂Nodes))
        ∂Nodes += ∂nodes_change
        ∂x .+= ∂xchange
        ∂bn = Tangent{typeof(bn)}(;Beams = ∂Beams,Nodes = ∂Nodes)
        ∂adj = sum(solst.u .* ∂st ,dims = 3) 
        
        # @show 1/(1/sol[2](0)[5] + 1/sol[4](0)[5]) + sol[6](0)[5]
        return NoTangent(),∂x,∂bn,∂adj[:,:,1]
    end
    (toArray(sol),(;Beams = bn.Beams,Nodes = nodes_),solst.u),back_groundstr
end 