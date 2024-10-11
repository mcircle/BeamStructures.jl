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
        # dnode = Bo(reduce((init,x) -> init .+ rȳ[[2,3,1],x],1:length(beamnbrs[1]),init = zeros(T,3))...,zeros(T,3)...)
        return NoTangent(),dnode,∂beams,∂y,NoTangent()
    end 
    return sol, reduceposat_back
end 

function CRC.rrule(::typeof(reduceposat),node::Bo,beams,y::AbstractArray{T},facs,beamnbrs) where{T,Bo<:Boundary{T}} 
    sol = reduceposat(node,beams,y,facs,beamnbrs)
    function reduceposat_back(ȳ)
        
        ∂beams= Tangent{typeof(beams)}(;ntuple(x->keys(beams)[x]=>ZeroTangent(),length(beams))...)
        ∂y = zero(y)
        ∂f = zero(facs)
        dnode = Tangent{Bo}()
    
        for (n,b) in enumerate(beamnbrs[1])
            sp, spback = CRC.rrule(scalepos,beams[b],y[2:4,2,b],Val(2))

            tmp = 2 .* facs[b] .*([node.ϕ,node.x,node.y] .- sp) .* ȳ

            _,dbeam,dy,_ = spback(-tmp)    
            
            ∂beams += Tangent{typeof(beams)}(;Symbol(:Beam_,b) => dbeam)
            ∂y[2:4,2,b] .+= dy
            ∂f[b] = abs2.([node.ϕ,node.x,node.y] .- scalepos(beams[b],y[2:4,2,b],Val(2)))' * ȳ
            dnode = dnode + Tangent{Bo}(;x = tmp[2], y = tmp[3],ϕ = tmp[1])
        end
        return NoTangent(),dnode,∂beams,∂y,∂f, NoTangent()
    end 
    return sol, reduceposat_back
end 

function CRC.rrule(::typeof(normfactor_m), b::Beam{T}) where{T}
    y = normfactor_m(b)
    function pullback(ȳ)
        nom = 1/(b.E * b.w * b.h^3)
        ∂b = Tangent{Beam{T}}(;
            l = 12 * ȳ  * nom,
            h = -36 * b.l * ȳ * nom / b.h,
            w = -12 * b.l * ȳ * nom / b.w ,
            E = -12 * b.l * ȳ * nom / b.E
        )
        return (NoTangent(), ∂b)
    end
    return y, pullback
end

function CRC.rrule(::typeof(normfactor_f), b::Beam{T}) where{T}
    y = normfactor_f(b)
    function pullback(ȳ)
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
    return y, pullback
end

function CRC.rrule(::typeof(normvector), b::Beam{T}) where{T}
    y = normvector(b)
    function pullback(ȳ)
        _, pb_m = CRC.rrule(normfactor_m, b)
        _, pb_f = CRC.rrule(normfactor_f, b)
        _, ∂b_m = pb_m(ȳ[1])
        _, ∂b_f1 = pb_f(ȳ[2])
        _, ∂b_f2 = pb_f(ȳ[3])
        ∂beam = ∂b_m + ∂b_f1 + ∂b_f2

        return (NoTangent(), ∂beam)
    end
    return y, pullback
end

function CRC.rrule(::typeof(scaleforce),b::Beam,y)
    nv = normvector(b)
    result =  y ./ nv
    
    function pullback(ȳ)
        ∂y = ȳ ./ nv
        ∂nv = -ȳ .* y ./ nv.^2
        
        _, pb_nv = CRC.rrule(normvector, b)
        _, ∂beam = pb_nv(∂nv)
        return (NoTangent(),∂beam,∂y)
    end
    
    return result, pullback
end

function CRC.rrule(::typeof(initialize_beam),node::No,beam::Beam{BT},pars::AbstractVector) where{BT,No<:Boundary}
    x,y,θ0, = node.x,node.y,node.ϕ
    l,θs,κ = beam.l,beam.θs,beam.κ0
    m,fx,fy = pars
    norm,normback = rrule(normvector,beam)
    m *= norm[1] #am Balkenelement
    fx *= norm[2] #am Balkenelement
    fy *= norm[3] #am Balkenelement
    result = [m,θ0 + θs,x ./l,y ./l,fx,fy,κ*l]
    function init_beam_back(ȳ)
        _,∂beam = normback(pars[1:3] .* ȳ[[1,5,6]])
        ∂l = -result[3]/l * ȳ[3] - result[4]/l * ȳ[4] + κ * ȳ[7]
        ∂beam += Tangent{Beam{BT}}(l = ∂l,κ0 = l*ȳ[7],θs = ȳ[2])
        ∂node = Tangent{No}(;x = ȳ[3]/l,y=ȳ[4]/l,ϕ = ȳ[2])
        ∂pars = norm .* ȳ[[1,5,6]]
        return NoTangent(),∂node,∂beam,∂pars
    end 
    return result, init_beam_back
end 

function CRC.rrule(::typeof(initialize_beam),node::Clamp{CT},beam::Beam{BT},pars::AbstractVector) where{CT,BT}
    x,y,θ0 = node.x,node.y,node.ϕ
    l,θs,κ = beam.l,beam.θs,beam.κ0 
    norm,normback = rrule(normvector,beam)
    m,fx,fy = pars .* norm #am Balkenelement
    result = [m,θ0 + θs,x/l,y/l,fx,fy,κ*l]
    function init_back_clamp(ȳ)

        _,∂beam = normback(pars .* ȳ[[1,5,6]])
        ∂l = -result[3]/beam.l * ȳ[3] - result[4]/beam.l * ȳ[4] + beam.κ0 * ȳ[7]
        ∂beam += Tangent{Beam{BT}}(;l = ∂l,κ0 = beam.l*ȳ[7],θs = ȳ[2])

        # dn = [1/l,1/l,1] .* ȳ[[3,4,2]]
        ∂node = Tangent{Clamp{CT}}(; x = ȳ[3]/beam.l,y =ȳ[4]/beam.l ,ϕ = ȳ[2])

        ∂pars = norm .* ȳ[[1,5,6]]
        return NoTangent(),∂node,∂beam,∂pars
    end 
    result,init_back_clamp
end 

function CRC.rrule(::typeof(initialize_beam),beams::NamedTuple,nodes::NamedTuple,pars,nodepos,i::Int)
    x_idxs = nodepos[i]
    result,init_back = CRC.rrule(initialize_beam,nodes[nodepos[i]],beams[i],pars[:,x_idxs])
    function init_u0_back(ȳ)
        x_idxs = nodepos[i]
        ∂nodes = Tangent{typeof(nodes)}(;ntuple(x-> keys(nodes)[x]=>ZeroTangent(),length(nodes))...) 
        ∂beams = Tangent{typeof(beams)}(;ntuple(x-> keys(beams)[x]=>ZeroTangent(),length(beams))...) 
        
        ∂pars = zero(pars)
        _,dnode,dbeam,dpars = init_back(ȳ)
        ∂pars[:,x_idxs] += dpars 
        ∂beams += Tangent{typeof(beams)}(;Symbol(:Beam_,i) => dbeam)
        ∂nodes += Tangent{typeof(nodes)}(;Symbol(:Node_,x_idxs) => dnode)
        return NoTangent(), ∂beams,∂nodes,∂pars,NoTangent(),NoTangent()
    end 
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
                ∂fac[b] = -sum(f .* ȳ) 
            else
                _,dbeam,dyrr_ = rr_norm(ȳ)
                ∂fac[b] = sum(f .* ȳ) 
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
    idcs = LinearIndices(CartesianIndices((1:3,1:fld(length(residuals),3))))
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
    
    cis = getindices(size(adj,1))
    idx = 1
    rrposdict = Dict{Int,Function}()
    rrforcedict = Dict{Int,Function}()
    
    for n in axes(adj,1)
        node = bn.Nodes[n]
        beams = findbeamsatnode(node,n,cis)
        if any(!isempty,beams)
            res, pullback_reduceforceat = rrule(reduceforceat,node,bn.Beams,y,adj[cis],beams)
            residuals[:,idx] .= res
            rrforcedict[idx] = pullback_reduceforceat
            res, pullback_reduceposat = rrule(reduceposat,node,bn.Beams,y,adj[cis],beams)
            residuals[:,idx +1] .= res
            rrposdict[idx+1] = pullback_reduceposat
            idx += 2 
        end 
    end

    function residuals!_back(ȳ)
        ∂res = InplaceableThunk(dself -> dself .+= ȳ,@thunk(copy(ȳ)))
        ∂y = zero(y)
        ∂adj = zero(adj)
        ∂beams = Tangent{typeof(bn.Beams)}(;ntuple(x->keys(bn.Beams)[x] =>ZeroTangent(),length(bn.Beams))...)        
        ∂nodes = Tangent{typeof(bn.Nodes)}(;ntuple(x->keys(bn.Nodes)[x] =>ZeroTangent(),length(bn.Nodes))...)
        idx = 1
        for ind in axes(adj,1)
            node = bn.Nodes[ind]
            beams = findbeamsatnode(node,ind,cis)
            if any(!isempty,beams)
                pb = rrforcedict[idx]
                _,dnode,dbeams,dy,dadj,_ =pb(ȳ[:,idx])
                ∂nodes += Tangent{typeof(bn.Nodes)}(;Symbol(:Node_,ind) => dnode)
                ∂beams += dbeams
                ∂y .+= dy
                ∂adj[cis] .+= dadj
                _,dnode,dbeams,dy,dadj,_ = rrposdict[idx+1](ȳ[:,idx+1])
                ∂nodes += Tangent{typeof(bn.Nodes)}(;Symbol(:Node_,ind) => dnode)
                ∂beams +=  dbeams
                ∂y .+= dy
                ∂adj[cis] .+= dadj
                idx += 2
            end 
        end 
        ∂bn = Tangent{typeof(bn)}(;Beams = ∂beams,Nodes = ∂nodes)
        return NoTangent(),∂res,∂adj,∂y,∂bn
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

function CRC.rrule(::typeof(addposition),node::BT,x) where{T,BT<:Boundary{T}}
    res = addposition(node,x)
    addpos_back(ȳ) = (NoTangent(),ȳ,[ȳ.x,ȳ.y,ȳ.ϕ])
    return res,addpos_back
end     

function CRC.rrule(::typeof(addpositions),nodes,x)
    branches = 1:length(nodes) |> filter((x)->!isa(values(nodes[x]),Clamp))
    newnodes = map((node,pos)->keys(nodes)[node] => addposition(nodes[node],pos),branches,eachcol(x[:,1:length(branches)]))
    function addposs_back(ȳ)
        ∂x = zero(x)
        for (colx,node) in zip(eachcol(∂x),branches)
            colx[1] = ȳ[node].x
            colx[2] = ȳ[node].y
            colx[3] = ȳ[node].ϕ
        end 
        return NoTangent(),ȳ,∂x
    end 
    
    return merge(nodes,newnodes),addposs_back
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