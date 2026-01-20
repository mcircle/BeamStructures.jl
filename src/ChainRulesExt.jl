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

function CRC.rrule(::Type{B},x,y,ϕ,fx,fy,mz) where{T<:Real,B<:Boundary{T}}
    function back_boundary(ȳ)
        # @show "now in rrule back_boundary for $(B)"
        return NoTangent(),ȳ.x,ȳ.y,ȳ.ϕ,ȳ.fx,ȳ.fy,ȳ.mz
    end 
    return B(x,y,ϕ,fx,fy,mz),back_boundary
end

function CRC.rrule(::typeof(+),b::B,nt::NamedTuple) where{T<:Real,B<:Boundary{T}}
    function back_add(ȳ)
        # @show "now in rrule add to boundary for $(B)"
        ∂b = ȳ 
        ∂nt = NoTangent()#NamedTuple{keys(nt)}(ntuple(i -> getproperty(ȳ, keys(nt)[i]), length(nt)))
        return NoTangent(), ∂b, ∂nt
    end 
    return b + nt,back_add
end

function scalepos_back(ȳ,y,beam::Beam{T}) where{T}
    dy = ȳ .* [1,beam.l,beam.l]
    dbeam = Tangent{Beam{T}}(;l = sum(y[2:3] .* ȳ[2:3] ) ,θe = -ȳ[1])
    return NoTangent(),dbeam, dy,NoTangent() 
end 

function CRC.rrule(::typeof(scalepos),beam::Beam{T},y,side::Val{2}) where{T}
    sol = scalepos(beam,y,side)
    sc_back(ȳ) = scalepos_back(ȳ,y,beam)
    return sol,sc_back
end 

function CRC.rrule(::typeof(reduceposat),node::LinearSlider{T},beams,y,beamnbrs) where{T}
    sol = reduceposat(node,beams,y,beamnbrs)
    function reduceposat_back(ȳ)
        rȳ = reshape(ȳ,3,:)
        
        ∂beams= Tangent{typeof(beams)}(;ntuple(x->keys(beams)[x]=>ZeroTangent(),length(beams))...)
        ∂y = zero(y)
        for (n,b) in enumerate(beamnbrs[1])
            _,dbeam,dy,_ = scalepos_back(rȳ[:,n],y[:,n],beams[b])    
            ∂beams -= Tangent{typeof(beams)}(;Symbol(:Beam_,b) => dbeam)
            ∂y[2:4,2,b] .-= dy
        end
        dϕ =  sum(rȳ[1,:])
        dx = sum(rȳ[2,:]) 
        dy = sum(rȳ[3,:]) 
        ds = dx * cos(dϕ) + dy * sin(dϕ)
        dir = node.s*(dx * sum(x->-sin(x),rȳ[1,:]) + dy * sum(x->cos(x),rȳ[1,:]))
        dnode = Tangent{LinearSlider{T}}(;x = dx,y = dy,dir = dir,ϕ = dϕ,s = ds)
        return NoTangent(),dnode,∂beams,∂y,NoTangent()
    end 
    return sol, reduceposat_back
end 

function CRC.rrule(::typeof(reduceposat),node::Bo,beams,y,beamnbrs) where{T,Bo<:Boundary{T}} 
    sol = reduceposat(node,beams,y,beamnbrs)
    function reduceposat_back(ȳ)
        rȳ = reshape(ȳ,3,:)
        
        ∂beams= Tangent{typeof(beams)}(;ntuple(x->keys(beams)[x]=>ZeroTangent(),length(beams))...)
        ∂y = zero(y)
        dnode = Tangent{Bo}(;x = 0,y = 0,ϕ = 0)
        for (n,b) in enumerate(beamnbrs[1])
            # _, spback = CRC.rrule(scalepos,beams[b],y[2:4,2,b],Val(2))
            
            _,dbeam,dy,_ = scalepos_back(rȳ[:,n],y[2:4,2,b],beams[b])     
            
            ∂beams -= Tangent{typeof(beams)}(;Symbol(:Beam_,b) => dbeam)
            ∂y[2:4,2,b] .-= dy
            # dnode += Tangent{Bo}(;x = beams[b].l * rȳ[2,n],y = beams[b].l *  rȳ[3,n],ϕ = rȳ[1,n])
            dnode += Tangent{Bo}(;x = rȳ[2,n],y =  rȳ[3,n],ϕ = rȳ[1,n])
        end
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
            #normalerweise  mit facs[b] multipliziert, wenn adj-Matrix berücktsichtigt wird
            _,dbeam,dy,_ =  spback(facs[b] .* rȳ[:,n])  
            ∂beams -= Tangent{typeof(beams)}(;Symbol(:Beam_,b) => dbeam)
            ∂y[2:4,2,b] .-=  dy
            
            ∂f[b] = ([node.ϕ,node.x,node.y] .- sp) ⋅ rȳ[:,n] 
            # normalerweise  mit facs[b] multipliziert
            dnode += Tangent{Bo}(;x = facs[b] * rȳ[2,n], y = facs[b] * rȳ[3,n] ,ϕ =  facs[b] * rȳ[1,n]) 
        end
        # dnode = Tangent{Bo}(;x = sum(rȳ[2,:]),y = sum(rȳ[3,:]),ϕ = sum(rȳ[1,:]))
        return NoTangent(),dnode,∂beams,∂y,∂f, NoTangent()
    end 
    return sol, reduceposat_back
end 

function pullback_normfactor_m(ȳ,b::Beam{T}) where{T}
    
    nom = 12/(b.E * (b.w) * (b.h)^3)
    ∂b = Tangent{Beam{T}}(;
        l = ȳ  * nom ,
        h = (-3 * b.l * ȳ * nom / b.h) ,
        w = -b.l * ȳ * nom / b.w  ,
        E =  -b.l * ȳ * nom / b.E
    )
    return (NoTangent(), ∂b)
end

function CRC.rrule(::typeof(normfactor_m), b::Beam{T}) where{T}
    y = normfactor_m(b)
    pullback_norm_m(ȳ) = pullback_normfactor_m(ȳ,b)  
    return y, pullback_norm_m
end

function pullback_normfactor_f(ȳ,b::Beam{T}) where{T}
    ∂m = b.l * ȳ #/ (1+exp(-b.l))
    fm, pb_m = CRC.rrule(normfactor_m, b)
    _, ∂b_m = pb_m(∂m)
    
    ∂b = Tangent{Beam{T}}(;
        l = fm * ȳ + ∂b_m.l,
        h =  ∂b_m.h ,
        w =  ∂b_m.w ,
        E =  -∂b_m.E 
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

function pullback_scaleforce(ȳ,y,nv,b,::Val{false})
    ∂y = ȳ ./ nv
    return (NoTangent(),ZeroTangent(),-∂y)
end

function pullback_scaleforce(ȳ,y,nv,b,::Val{true})
    ∂y = ȳ ./ nv
    ∂nv = - ∂y .* y ./ nv 
    _, ∂beam = pullback_normfactor(∂nv,b)
    return (NoTangent(),∂beam,∂y)
end 

function CRC.rrule(::typeof(scaleforce),b::Beam,y)
    nv , pb_nv = CRC.rrule(normvector, b)
    # nv = normvector(b)
    result =  y ./ nv
    # pullback_scaleforce(ȳ) = pullback_scaleforce(ȳ,y,nv,b,Val(true))
    function pullback(ȳ)
        ∂y = ȳ ./ nv
        ∂nv = - ∂y .* y ./ nv 
        _, ∂beam = pb_nv(∂nv)
        return (NoTangent(),∂beam,∂y)
    end
    
    return result, pullback
end



function pullback_init_beam(ȳ,node::No,beam::Beam{BT},pars) where{BT,No<:Boundary{T} where{T}}
    @inbounds mb,θb,xb,yb,fxb,fyb,κb = ȳ
    @inbounds m,fx,fy = pars 
    l,κ = beam.l,beam.κ0
    x,y = node.x,node.y
    _,∂beam = pullback_normfactor([m*mb, fx*fxb,fy*fyb] ,beam)
    ∂l = -x/l^2 * xb - y/l^2 * yb + κ * κb
    ∂beam += Tangent{Beam{BT}}(l = ∂l, κ0 = l*κb,θs = θb)

    dx = xb/l
    dy = yb/l
    ∂node = Tangent{No}(;x = dx,y=dy,ϕ = θb)
    ∂pars =  [mb,fxb,fyb] .* normvector(beam)
    # if isa(node,Branch)
    #     @show ∂pars
    # end 
    return NoTangent(),∂node,∂beam,∂pars
end 

function pullback_init_beam(ȳ,node::LinearSlider{T},beam::Beam{BT},pars) where{T,BT}
    l,κ = beam.l,beam.κ0
    x,y,trans = node.x,node.y,node.trans

    _,∂beam = pullback_normfactor(pars[1:3] .* ȳ[[1,5,6]],beam)
    ∂l = -x/l^2 * ȳ[3] - y/l^2 * ȳ[4] + κ * ȳ[7]
    ∂beam += Tangent{Beam{BT}}(l = ∂l,κ0 = l*ȳ[7],θs = ȳ[2])
    ∂node = Tangent{LinearSlider{T}}(;x = ȳ[3]/l,y=ȳ[4]/l,ϕ = ȳ[2])
    ∂pars = trans .* normvector(beam) .* ȳ[[1,5,6]]
    return NoTangent(),∂node,∂beam,∂pars
end

function pullback_init_beam(ȳ,node::Joint{T},beam::Beam{BT},pars) where{T,BT}
    l,κ = beam.l,beam.κ0
    x,y,trans = node.x,node.y,node.trans

    _,∂beam = pullback_normfactor(pars[1:3] .* ȳ[[1,5,6]],beam)
    ∂l = -x/l^2 * ȳ[3] - y/l^2 * ȳ[4] + κ * ȳ[7]
    ∂beam += Tangent{Beam{BT}}(l = ∂l,κ0 = l*ȳ[7],θs = ȳ[2])
    ∂node = Tangent{Joint{T}}(;x = ȳ[3]/l,y=ȳ[4]/l,ϕ = ȳ[2])
    ∂pars = trans .* normvector(beam) .* ȳ[[1,5,6]]
    return NoTangent(),∂node,∂beam,∂pars
end

function CRC.rrule(::typeof(initialize_beam),node::No,beam::Beam{BT},pars::AbstractVector) where{BT,No<:Boundary}
    x,y,θ0, = node.x,node.y,node.ϕ
    l,θs,κ = beam.l,beam.θs,beam.κ0
    m,fx,fy = pars .* normvector(beam) #am Balkenelement
    result = [m,θ0 + θs,x./l,y./l,fx,fy,κ*l]
    init_beam_back(ȳ) = pullback_init_beam(ȳ,node,beam,pars)
    return result, init_beam_back
end 

function init_u0_pullback(ȳ,beams::NamedTuple,nodes::NamedTuple,pars,nodepos,i)
       
    ∂nodes = Tangent{typeof(nodes)}(;ntuple(x-> keys(nodes)[x]=>CRC.zero_tangent(nodes[x]),length(nodes))...) 
    ∂beams = Tangent{typeof(beams)}(;ntuple(x-> keys(beams)[x]=>CRC.zero_tangent(beams[x]),length(beams))...) 
    
    ∂pars = zero(pars)
    _,dnode,dbeam,dpars = pullback_init_beam(ȳ,nodes[nodepos[i]],beams[i],pars[:,i])
    ∂pars[:,i] += dpars 
    ∂beams += Tangent{typeof(beams)}(;Symbol(:Beam_,i) => dbeam)
    ∂nodes += Tangent{typeof(nodes)}(;Symbol(:Node_,nodepos[1]) => dnode)
    return NoTangent(), ∂beams,∂nodes,∂pars,NoTangent(),NoTangent()
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
        # @show No,idxs
        y_idxs = getforceindices(idxs)
        ∂y = zero(y)
        y_ = @view y[y_idxs]
        ∂y_ = @view ∂y[y_idxs] 
        ∂beams = Tangent{typeof(beams)}(;ntuple(x->keys(beams)[x]=>ZeroTangent(),length(beams))...)
        for (b,yb,dyb) in zip(vcat(idxs...),eachcol(y_),eachcol(∂y_))
            
            nv = normvector(beams[b])
            _,dbeam,dy = pullback_scaleforce(ȳ,yb,nv,beams[b],Val(b in idxs[1]))

            ∂beams += Tangent{typeof(beams)}(;Symbol(:Beam_,b) => dbeam)
            dyb .+= dy 
            
        end 

        ∂node = CRC.zero_tangent(No) # Tangent{No}(;fx = ȳ[2],fy = ȳ[3],mz = ȳ[1]) #CRC.zero_tangent(No) #
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
        ∂beams = Tangent{typeof(beams)}(;ntuple(x->keys(beams)[x]=>ZeroTangent(),length(beams))...)
        for (b,yb,dyb) in zip(vcat(idxs...),eachcol(y_),eachcol(∂y_))
            # scale = normvector(beams[b])
            nv = normvector(beams[b])
            _,dbeam,dy = pullback_scaleforce(factors[b].* ȳ,yb,nv,beams[b],Val(b in idxs[1]))

            ∂beams += Tangent{typeof(beams)}(;Symbol(:Beam_,b) => dbeam)
            dyb .+= dy 

            # f, rr_norm = CRC.rrule(scaleforce,beams[b],yb)
            # if b in idxs[2]
            #     _,dbeam,dyrr_ = rr_norm(-factors[b].*ȳ) 
            #     ∂fac[b] = sum(-f .* ȳ) 
                
            #     # if isa(node,Branch)
            #     #     @show b, dyrr_
            #     # end
            # else
            #     _,dbeam,dyrr_ = rr_norm(factors[b] .* ȳ)
            #     ∂fac[b] = sum(f .* ȳ) 
                
            # end 
            # ∂beams += Tangent{typeof(beams)}(;Symbol(:Beam_,b) => dbeam)
            # dyb .+= dyrr_
        end 
        
        ∂node = CRC.zero_tangent(No)#(;fx = ȳ[2],fy = ȳ[3],mz = ȳ[1])
        return NoTangent(),∂node,∂beams,∂y,∂fac,NoTangent()
    end 
    return sol,force_bound_back
end 

function CRC.rrule(::typeof(reduceforceat),node::LinearSlider{No},beams,y::AbstractArray{T,3},factors,idxs) where{T,No}
    
    sol = reduceforceat(node,beams,y,factors,idxs)
    function force_bound_back(ȳ)
        y_idxs = getforceindices(idxs)
        ∂y = zero(y)
        ∂fac = zero(factors)
        y_ = @view y[y_idxs]
        ∂y_ = @view ∂y[y_idxs] 
        
        ∂beams = Tangent{typeof(beams)}(;ntuple(x->keys(beams)[x]=>ZeroTangent(),length(beams))...)
        for (b,yb,dyb) in zip(vcat(idxs...),eachcol(y_),eachcol(∂y_))
            
            f, rr_norm = CRC.rrule(scaleforce,beams[b],yb)
            if b in idxs[2]
                _,dbeam,dyrr_ = rr_norm(-ȳ) 
                ∂fac[b] = sum(trans .* f .* ȳ) 
            else
                _,dbeam,dyrr_ = rr_norm(ȳ)
                ∂fac[b] = -sum(trans .* f .* ȳ) 
            end 
            ∂beams += Tangent{typeof(beams)}(;Symbol(:Beam_,b) => dbeam)
            dyb .= dyrr_ #./ length(idxs)
        end 
        ∂node = CRC.zero_tangent(LinearSlider{No})#;fx = ȳ[2],fy = ȳ[3],mz = ȳ[1])
        return NoTangent(),∂node,∂beams,∂y,∂fac,NoTangent()
    end 
    return sol,force_bound_back
end 

@non_differentiable reduceforceat(n::Clamp,beams,y,idxs) 

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

function CRC.rrule(::typeof(residuals!),residuals,str::Structure,y,bn)
    ind = 1
    adj = str.AdjMat
    # idcs = LinearIndices(residuals)
    nodes = findall(x->!isapprox(x,0),LowerTriangular(adj))
    branches = count(x->isa(x,Branch),bn.Nodes)
    residuals_forces = @view residuals[:,1:branches]
    residuals_positions = @view residuals[:,branches+1:end]

    forces = 1
    positions = 1
    resposdict = Dict{Int,AbstractRange{Int}}()
    resforcedict = Dict{Int,Int}()
    rrposdict = Dict{Int,Function}()
    rrforcedict = Dict{Int,Function}()
    for n in getnodeswithbeams(adj,bn.Nodes)
        node = bn.Nodes[n]
        beams = findbeamsatnode(node,n,nodes)
        res, pullback_reduceforceat = rrule(reduceforceat,node,bn.Beams,y,beams)
        if !isempty(res)
            residuals_forces[:,forces] .= res
            resforcedict[n] = forces
            rrforcedict[n] = pullback_reduceforceat
            forces += 1
        end 
        res, pullback_reduceposat = rrule(reduceposat,node,bn.Beams,y,beams)
        if !isempty(res)
            idxs = positions:positions + size(res,2) - 1
            resposdict[n] = idxs
            rrposdict[n] = pullback_reduceposat
            residuals_positions[:,idxs] .= res
            positions = idxs[end] + 1
        end  
    end
    function residuals!_back(ȳt)
        ȳ = CRC.unthunk(ȳt)
        ∂res = ZeroTangent()
        ∂y = zero(y) 
        ȳ_forces = @view ȳ[:,1:branches]
        ȳ_positions = @view ȳ[:,branches+1:end]

        ∂beams = Tangent{typeof(bn.Beams)}(;ntuple(x->keys(bn.Beams)[x] =>ZeroTangent(),length(bn.Beams))...)        
        ∂nodes = Tangent{typeof(bn.Nodes)}(;ntuple(x->keys(bn.Nodes)[x] =>ZeroTangent(),length(bn.Nodes))...)
        for (ind,pos) in resforcedict
            _,dnode,dbeams,dy,_ =rrforcedict[ind](ȳ_forces[:,pos])
            ∂nodes += Tangent{typeof(bn.Nodes)}(;Symbol(:Node_,ind) => dnode)
            ∂beams += dbeams
            # @show "F",pos,ind, dbeams
            ∂y .+= dy 
        end

        for (ind,pos) in resposdict
            _,dnode,dbeams,dy,_ = rrposdict[ind](ȳ_positions[:,pos])
            ∂nodes += Tangent{typeof(bn.Nodes)}(;Symbol(:Node_,ind) => dnode)
            ∂beams += dbeams
            # @show "P", pos, ind, dbeams
            ∂y .+= dy
        end 
        ∂bn = Tangent{typeof(bn)}(;Beams = ∂beams,Nodes = ∂nodes)
        # @show ∂beams[:Beam_6]
        # ∂y
        return NoTangent(),∂res,NoTangent(),∂y,∂bn
    end 

    return residuals,residuals!_back
end 

function CRC.rrule(::typeof(residuals!),residuals::Matrix,str::Structure,y::EnsembleSolution,bn)
    y = toArray(y)
    res,back = rrule(residuals!,residuals,str,y,bn)
    return res,back
end

function CRC.rrule(::typeof(residuals!),residuals::AbstractMatrix,adj_::AbstractMatrix{T},y::AbstractArray,bn) where{T}

    ind = 1
    idcs = getindices(size(adj_,1))

    adj = ifelse.(adj_ .> 1,one(T),adj_ )
    adj = ifelse.(adj .< 0, zero(T),adj_ )
    # start = 1
    branches = count(x->isa(x,Branch),bn.Nodes)
    residuals_forces = @view residuals[:,1:branches]
    residuals_positions = @view residuals[:,branches+1:end]

    forces = 1
    positions = 1
    resposdict = Dict{Int,AbstractRange{Int}}()
    resforcedict = Dict{Int,Int}()
    rrposdict = Dict{Int,Function}()
    rrforcedict = Dict{Int,Function}()
    for n in getnodeswithbeams(adj,bn.Nodes)#axes(adj,1)
        node = bn.Nodes[n]
        beams = findbeamsatnode(node,n,idcs)
        
        res, pullback_reduceforceat = rrule(reduceforceat,node,bn.Beams,y,adj[idcs],beams)
        if !isempty(res)
            residuals_forces[:,forces] .= res
            resforcedict[n] = forces
            rrforcedict[n] = pullback_reduceforceat
            forces += 1
        end 
        res, pullback_reduceposat = rrule(reduceposat,node,bn.Beams,y,adj[idcs],beams)
        if !isempty(res)
            idxs = positions:positions + size(res,2) - 1
            resposdict[n] = idxs
            rrposdict[n] = pullback_reduceposat
            residuals_positions[:,idxs] .= res
            positions = idxs[end] + 1
        end  
    end    
    function residuals!_back(ȳ)
        ȳ_ = CRC.unthunk(ȳ)
        ∂res = @thunk(CRC.zero_tangent(residuals))
        ∂y = zero(y)
        ∂fac = zero(adj)
        ∂beams = Tangent{typeof(bn.Beams)}(;ntuple(x->keys(bn.Beams)[x] =>CRC.zero_tangent(bn.Beams[x]),length(bn.Beams))...)        
        ∂nodes = Tangent{typeof(bn.Nodes)}(;ntuple(x->keys(bn.Nodes)[x] =>CRC.zero_tangent(bn.Nodes[x]),length(bn.Nodes))...)

        ȳ_forces = @view ȳ_[:,1:branches]
        ȳ_positions = @view ȳ_[:,branches+1:end]


        for (ind,pos) in resforcedict
            _,dnode,dbeams,dy,df =rrforcedict[ind](ȳ_forces[:,pos])
            ∂nodes += Tangent{typeof(bn.Nodes)}(;Symbol(:Node_,ind) => dnode)
            ∂beams += dbeams
            # @show ind, pos, dbeams
            ∂y .+= dy 
            # ∂fac[idcs] .+= df
        end
        
        for (ind,pos) in resposdict
            _,dnode,dbeams,dy,df = rrposdict[ind](ȳ_positions[:,pos])
            ∂nodes += Tangent{typeof(bn.Nodes)}(;Symbol(:Node_,ind) => dnode)
            ∂beams +=  dbeams
            ∂y .+= dy
            # @show ind, pos, dbeams
            # ∂fac[idcs] .+= df
        end 
        # @show ∂fac
        ∂bn = Tangent{typeof(bn)}(;Beams = ∂beams,Nodes = ∂nodes)
        return NoTangent(),∂res,∂fac,∂y,∂bn
    end    
    return residuals,residuals!_back
end 

function CRC.rrule(::typeof(residuals!),residuals::Matrix,adj::AbstractMatrix,y::EnsembleSolution,bn)
    y = toArray(y)
    res,back = rrule(residuals!,residuals,adj,y,bn)
    return res,back
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
    backaddb(ȳ) = (NoTangent(), ȳ,ȳ)
    return a + b,backaddb
end 

function CRC.rrule(::typeof(addposition),mv::LinearSlider{BT},disp) where{BT}
    res = addposition(mv,disp)
    function addpos_back_x(ȳ)
        return NoTangent(),ȳ,ȳ.s
    end
    return res,addpos_back_x
end 

function CRC.rrule(::typeof(addposition),mv::Joint{BT},disp) where{BT}
    res = addposition(mv,disp)
    function addpos_back_x(ȳ)
        return NoTangent(),ȳ,ȳ.ϕ
    end
    return res,addpos_back_x
end 

function pullback_addposition(ȳ,∂x,branches)
    for (colx,node) in zip(eachcol(∂x),branches)
        colx[1] += ȳ[node].x
        colx[2] += ȳ[node].y 
        colx[3] += ȳ[node].ϕ 
    end 

    return NoTangent(),ȳ,∂x
end 

function CRC.rrule(::typeof(addpositions),nodes,xpos)
    branches = 1:length(nodes) |> filter((x)->isa(values(nodes[x]),Branch))
    newnodes = map((node,pos)->keys(nodes)[node] => addposition(nodes[node],pos),branches,eachcol(xpos))
    addpos_back(ȳ) = pullback_addposition(ȳ,zero(xpos),branches)
    return merge(nodes,newnodes),addpos_back
end

function CRC.rrule(::typeof(changestartnodes),nodes,x)
    anz = length(nodes) - count(x->isa(x,Clamp),values(nodes))
    xpos = @view x[:,1:anz]
    xforces = @view x[:,anz+1:end]
    nodes_,back = CRC.rrule(addpositions,nodes,xpos)
    function changeback(ȳ)
        # @show ȳ
        ∂x = zero(x)
        ∂x[:,1+anz:end] .+= ȳ[1]
        _,∂nodes,∂xpos = back(ȳ[2])
        ∂x[:,1:anz] .+= ∂xpos
        return NoTangent(),∂nodes,∂x
    end 
    return (xforces,nodes_),changeback
end 

function output_func(u,∂sol,beams,nodes,nodepos,x,i)
    x_idxs = nodepos[i]
    _,dnode,dbeam,∂xforces = pullback_init_beam(u(0),nodes[x_idxs],beams[i],x[:,i])
    # @show i,dbeam
    if isa(nodes[x_idxs],Branch)
        ∂xforces .+= ∂sol[[1,5,6],1,i] .* normvector(beams[i])
    end
    ∂beams = Tangent{typeof(beams)}(;Symbol(:Beam_,i) => dbeam)
    ∂nodes = Tangent{typeof(nodes)}(;Symbol(:Node_,x_idxs) => dnode)
    (∂beams,∂nodes,∂xforces),false
end 

function reduction_func!(u,data,I)
    
    # ∂pars,∂Beams,∂Nodes = u
    for ((∂beam,∂node,∂xforces),i) in zip(data,I)
        # @show i, ∂beam
        u[1][:,i] .+= ∂xforces
        u[3] += ∂node
        u[2] += ∂beam
    end 
    u,false
end 

function totangent(::T,x) where{BT,T<:CRC.StructuralTangent{BT}}
    T((getproperty(x, fname) for fname in fieldnames(x))...)
end

function u_init(x, bn,dbeams,dnodes)
    [zero(x),
    Tangent{typeof(bn.Beams)}(;ntuple(x->keys(bn.Beams)[x] =>CRC.StructuralTangent{typeof(bn.Beams[x])}(CRC.backing(dbeams[x])),length(bn.Beams))...),
    Tangent{typeof(bn.Nodes)}(;ntuple(x->keys(bn.Nodes)[x] =>CRC.StructuralTangent{typeof(bn.Nodes[x])}(CRC.backing(dnodes[x])),length(bn.Nodes))...)]
end 

function CRC.rrule(str::Structure,x::AbstractArray{T,N},bn::NamedTuple) where{T,N}
    nodepos = getstartnodes(str.AdjMat)

    (xforces,nodes_),change_pullback = CRC.rrule(changestartnodes,bn.Nodes,x)
    
    function prob_func(prob,i,repeat) 
        u0 = initialize_beam(bn.Beams,nodes_,xforces,nodepos,i) 
        remake(prob;u0 = u0)
    end 

    ensprob  =  EnsembleProblem(prob;prob_func = prob_func)

    solforward = solve(ensprob,str.Solver,
                EnsembleThreads(),
                reltol = 1e-6,abstol = 1e-6,
                trajectories = length(bn.Beams)
                )
    bn_ = (;Beams = bn.Beams,Nodes = nodes_)
    function back_ode(ȳ)
        @inbounds ∂sol,∂bn = ȳ
        ∂Beams = ∂bn.Beams   
        ∂Nodes = ∂bn.Nodes

        function prob_func(prob,i,repeat) 
            u0 = @view ∂sol[:,2,i]
            remake(prob;u0 = u0,p = solforward[i],tspan = (one(T),zero(T)))
        end
        ensprob  =  EnsembleProblem(vjpprob;prob_func = prob_func,
                                    output_func = (sol,i) -> output_func(sol,∂sol,bn.Beams,nodes_,nodepos,xforces,i),
                                    reduction = (u,data,I) -> reduction_func!(u,data,I),
                                    u_init = u_init(xforces,bn,∂Beams,∂Nodes))
        solp = solve(ensprob,str.Solver,
                    EnsembleThreads(),
                    reltol = 1e-6,abstol = 1e-6,
                    trajectories = length(bn.Beams)
                    )
        ∂xforces,dBeams,dNodes = solp
        # @show ∂xforces
        # @show dBeams
        _,∂nodes_change, ∂x = change_pullback((∂xforces,dNodes))
        # @show ∂x
        ∂bn = Tangent{typeof(bn)}(;Beams = dBeams,Nodes = dNodes)
        return NoTangent(),∂x,∂bn,ZeroTangent(),NoTangent()
    end
    (toArray(solforward),bn_),back_ode
end 

function Base.inv(x::Tangent{T,NT}) where{T<:Beam,NT}
    fields = fieldnames(T)
    CRC.StructuralTangent{T}((;map((f) -> f => Float32(inv(getproperty(x,f))), fields)...))
end
function Base.sqrt(x::Tangent{T,NT}) where{T<:Beam,NT}
    fields = fieldnames(T)
    CRC.StructuralTangent{T}((;map((f) -> f => Float32(sqrt(getproperty(x,f))), fields)...))
end
function Base.:*(x::Tangent{T},y::Tangent{T}) where{T}
    fields = fieldnames(T)
    CRC.StructuralTangent{T}((;map((f) -> f => Float32(getproperty(x,f) * getproperty(y,f)), fields)...))
end
function Statistics.realXcY(x::Tangent{T,NT},y::Tangent{T,NT}) where{T<:Beam,NT}
    x * y 
end

function CRC.rrule(str::GroundStructure,x::AbstractMatrix{T},bn::NamedTuple,adj) where{T}
    
    nodepos = getstartnodes(adj)
    cbeams = reduce(+,1:size(adj,1)-1)
    
    (xforces,nodes_),change_pullback = CRC.rrule(changestartnodes,bn.Nodes,x)
    function prob_func(prob,i,repeat) 
        u0 = initialize_beam(bn.Beams,nodes_,xforces,nodepos,i)
        remake(prob;u0 = u0)
    end 

    ensprob  =  EnsembleProblem(prob;prob_func = prob_func)

    sol = solve(ensprob,str.Solver,
                EnsembleThreads(),
                reltol = 1e-6,abstol = 1e-6,
                trajectories = cbeams
                )   
    bn_ = (;Beams = bn.Beams,Nodes = nodes_)
    function back_groundstr(ȳ)
        @inbounds ∂sol,∂bn = ȳ
        ∂Beams = ∂bn.Beams   
        ∂Nodes = ∂bn.Nodes
        
        function prob_func(prob,i,repeat) 
            u0 = @view ∂sol[:,2,i]
            remake(prob;u0 = u0,p = sol[i],tspan = (one(T),zero(T)))
        end
        #integriere Rückwärtsproblem
        ensprob  =  EnsembleProblem(vjpprob;prob_func = prob_func,
                                            output_func = (sol,i) -> output_func(sol,∂sol,bn.Beams,nodes_,nodepos,xforces,i),
                                            reduction = (u,data,I) -> reduction_func!(u,data,I),
                                            u_init = u_init(xforces,bn,∂Beams,∂Nodes)
                                            )
        solp = solve(ensprob,str.Solver,
                    EnsembleThreads(),
                    reltol = 1e-6,abstol = 1e-6,
                    trajectories = cbeams
                    )
        ∂xforces,dBeams,dNodes = solp
        # @show dNodes
        ∂adj = ZeroTangent()#one(adj)
        _,∂nodes_change, ∂x = change_pullback((∂xforces,dNodes))
        # @show ∂x
        ∂bn = Tangent{typeof(bn)}(;Beams = dBeams,Nodes = dNodes)
        return NoTangent(),∂x,∂bn,∂adj,ZeroTangent()
    end
    (sol,bn_),back_groundstr
end 

function output_func_admittance_back(sol,i)
    sol[end],false
end

function reduction_func_admittance_back!(u,data,I,ad)
    for ((i,j),d) in zip(I,data)
        u[i,j,:] .+=  d 
    end
    u,false
end

function CRC.rrule(::typeof(admittance_matrix),solfw::EnsembleSolution{T,N,S},adj,str,beams) where{T,N,S}
    
    idxs =getindices(size(adj,1))
    lensol = length(solfw)
    function prob_func(prob,i,repeat) 
        u0 = zeros(T,7,3)
        u0[[2],1] .= one(T) #du/dx
        u0[[3],2] .= one(T) #du/dy
        u0[[4],3] .= one(T) #du/dθ
        remake(prob;u0 = u0,p = solfw[i],tspan = (one(T),zero(T)))
    end
    #dadj Gradient der Steifigkeitsmatrix berechnen  
    ensprob2  =  EnsembleProblem(vjpprob;prob_func = prob_func,
                                output_func = (bsol,i) ->  output_func_admittance(bsol,i,solfw),
                                reduction = (u,data,I) -> reduction_func_admittance_!(u,data,idxs[I],solfw,adj,beams),
                                u_init = [zeros(eltype(adj),3 .* size(adj)...),[]]
                                )
    sol = solve(ensprob2,str.Solver,save_on = false,save_start = false,save_end = true,
                EnsembleThreads(),
                reltol = 1e-9,abstol = 1e-9,
                trajectories = lensol
                )

    ad,d_out = sol.u
    function admittance_back(ȳt)
        ȳ = CRC.unthunk(ȳt)
        
        ∂adj = zeros(T,size(adj)...)
        Δd = Matrix{eltype(adj)}(undef,3,3)
        # dtmp = Matrix{eltype(adj)}(undef,3,3)
        for (id,beam,d0,sol) in zip(idxs,beams,d_out,solfw)
            i_ = 3 * (id[1]-1)+1:3*id[1] 
            j_ = 3 * (id[2]-1)+1:3*id[2]
            x,y = sol[end][3:4] .- sol[1][3:4]

            Δd .=  -(view(ȳ,j_,i_)' + view(ȳ,i_,j_))
            # Pullback der Transformationen von d

            Δd[2, :] .+= y .* Δd[1, :]
            Δd[3,:]  .-= x .* Δd[1, :] 
            Δd[1, :] .*= -1
            Δd .+= (view(ȳ,i_,i_) + view(ȳ,j_,j_))

            # Pullback der elementweisen Multiplikation mit d0
            ∂adj[id] += sum(d0 .* Δd ) 
        end
        return NoTangent(),NoTangent(),∂adj,NoTangent(),NoTangent()
    end 
    ad,admittance_back
end 

function CRC.rrule(::typeof(effective_stiffness),k::AbstractMatrix{T},u) where{T}
    nomoveidcs = [u[1]...]
    moveables = [u[2]...]
    ids = setdiff(axes(k,1),nomoveidcs)
    ke = @view k[ids,ids]
    f = zeros(T,size(k,1))
    f[moveables] .= [one(T)]
    mov = ke \ f[ids]
    
    e_out = mov' * ke * mov
    function compute_energy_pullback(Δenergy)
        Δ_dot = CRC.unthunk(Δenergy)
        # Zurück in volle Matrix
        Δ_k = zero(k)
        Δ_k[ids, ids] .= - Δ_dot * mov * mov' 
        
        return CRC.NoTangent(),Δ_k, NoTangent()
    end
    e_out,compute_energy_pullback
end 

function CRC.rrule(::typeof(changenode),bn,nodes,nt)
    bn_out = changenode(bn,nodes,nt)
    # @show bn
    # @show bn_out
    function changenodeback(ȳ)
        # @show ȳ
        return CRC.NoTangent(), ȳ, CRC.NoTangent(), CRC.NoTangent()
    end 
    return bn_out,changenodeback
end
