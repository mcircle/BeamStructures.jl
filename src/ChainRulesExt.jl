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
        return NoTangent(),ȳ.x,ȳ.y,ȳ.ϕ,ȳ.fx,ȳ.fy,ȳ.mz
    end 
    return B(x,y,ϕ,fx,fy,mz),back_boundary
end

function CRC.rrule(::typeof(+),b::B,nt::NamedTuple) where{T<:Real,B<:Boundary{T}}
    function back_add(ȳ)
        return NoTangent(), ȳ, ZeroTangent()
    end 
    return b + nt,back_add
end

function scalepos_back(ȳ,y,beam::B) where{T,B<:BeamElement{T}}
    dy = ȳ .* [1,beam.l,beam.l]
    dbeam = Tangent{B}(;l = sum(y[2:3] .* ȳ[2:3] ) ,θe = -ȳ[1])
    return NoTangent(),dbeam, dy,NoTangent() 
end 

function CRC.rrule(::typeof(scalepos),beam::B,y,side::Val{2}) where{T,B<:BeamElement{T}}
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
    
    function reduceposat_back(ȳ,∂beams = zero_tangent(beams))
        rȳ = reshape(ȳ,3,:)
        
        # ∂beams= Tangent{typeof(beams)}(;ntuple(x->keys(beams)[x]=>ZeroTangent(),length(beams))...)
        ∂y = zero(y)
        dnode = Tangent{Bo}(;x = 0,y = 0,ϕ = 0)
        for (n,b) in enumerate(beamnbrs[1])
            _,dbeam,dy,_ = scalepos_back(rȳ[:,n],y[2:4,2,b],beams[b])     
            ∂beams -= Tangent{typeof(beams)}(;Symbol(:Beam_,b) => dbeam)
            ∂y[2:4,2,b] .-= dy
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

function pullback_normfactor_m(ȳ,b::B,nom) where{T,B<:BeamElement{T}}
    
    ∂b = Tangent{B}(;
        l = ȳ  * nom ,
        h = (-3 * b.l * ȳ * nom / b.h) ,
        w = -b.l * ȳ * nom / b.w ,
        E = -b.l * ȳ * nom / b.E
    )
    return (NoTangent(), ∂b)
end

function CRC.rrule(::typeof(normfactor_m), b::B) where{T,B<:BeamElement{T}}
    nf = normfactor(b)
    pullback_norm_m(ȳ) = pullback_normfactor_m(ȳ,b,nf)  
    return y, pullback_norm_m
end

function pullback_normfactor_f(ȳ,b::B) where{T,B<:BeamElement{T}}
    ∂m = b.l * ȳ 
    fm, pb_m = CRC.rrule(normfactor_m, b)
    _, ∂b_m = pb_m(∂m)
    
    ∂b = Tangent{B}(;
        l = fm * ȳ + ∂b_m.l,
        h =  ∂b_m.h ,
        w =  ∂b_m.w ,
        E =  -∂b_m.E 
    )
    return (NoTangent(), ∂b)
end

function CRC.rrule(::typeof(normfactor_f), b::BeamElement{T}) where{T}
    y = normfactor(b)
    pullback_norm_f(ȳ) = pullback_normfactor_f(ȳ,b)
    return y, pullback_norm_f
end

function pullback_normfactor(ȳ,b::BeamElement{T},nom) where{T}
    all(ȳ .== 0) && return (NoTangent(),ZeroTangent())

    # Ableitung von m nach b-Parametern
    ∂m_l = nom
    ∂m_h = -3 * b.l * nom / b.h
    ∂m_w = -b.l * nom / b.w
    ∂m_E = -b.l * nom / b.E
    
    # Ableitung von f = b.l * m nach b-Parametern
    ∂f_l = nom + b.l * ∂m_l
    ∂f_h = b.l * ∂m_h
    ∂f_w = b.l * ∂m_w
    ∂f_E = b.l * ∂m_E
    
    # Zusammenführen (ȳ = [ȳ_m, ȳ_f1, ȳ_f2])
    ∂beam = Tangent{typeof(b)}(;
        l = ȳ[1] * ∂m_l + ȳ[2] * ∂f_l + ȳ[3] * ∂f_l,
        h = ȳ[1] * ∂m_h + ȳ[2] * ∂f_h + ȳ[3] * ∂f_h,
        w = ȳ[1] * ∂m_w + ȳ[2] * ∂f_w + ȳ[3] * ∂f_w,
        E = ȳ[1] * ∂m_E + ȳ[2] * ∂f_E + ȳ[3] * ∂f_E,
    )
    return (NoTangent(), ∂beam)
end

function pullback_normfactor(ȳ,b::BeamElement{T}) where{T}
    nom = normfactor(b)
    return pullback_normfactor(ȳ,b,nom)
end

function CRC.rrule(::typeof(normvector), b::B) where{T,B<:BeamElement{T}}
    y = normfactor(b)
    pullback_norm(ȳ) = pullback_normfactor(ȳ,b,y)
    return y, pullback_norm
end

function pullback_scaleforce(ȳ,y,b,::Val{false})
    m = normfactor_m(b)
    f = b.l * m
    ∂y = @SVector [ȳ[1] / m, ȳ[2] / f,ȳ[3] / f]
    return (NoTangent(),ZeroTangent(),∂y)
end

function pullback_scaleforce(ȳ,y,b,::Val{true})

    m = normfactor_m(b)
    f = b.l * m
    ∂y = @SVector [-ȳ[1] / m,-ȳ[2] / f,-ȳ[3] / f]
    ∂nv = @SVector [-ȳ[1] * y[1] / m^2, -ȳ[2] * y[2] / f^2, -ȳ[3] * y[3] / f^2]
    ∂beam = pullback_normfactor(∂nv,b,m)[2]
    return (NoTangent(),∂beam,∂y)
end 

function CRC.rrule(::typeof(scaleforce),b::B,y) where{T,B<:BeamElement{T}}
    nv , pb_nv = CRC.rrule(normvector, b)

    result =  y ./ nv
    
    function pullback(ȳ)
        ∂y = ȳ ./ nv
        ∂nv = - ∂y .* y ./ nv 
        _, ∂beam = pb_nv(∂nv)
        return (NoTangent(),∂beam,∂y)
    end
    
    return result, pullback
end

vecOrnum(x::AbstractVector,::Val{N}) where{N} = x
vecOrnum(x::AbstractVector,::Val{1}) = only(x)
vecOrnum(x::AbstractVector) = vecOrnum(x,Val(length(x)))

function pullback_init_beam(ȳ,node::No,beam::B,pars) where{BT,B<:BeamElement{BT},No<:Boundary{T} where{T}}
    @inbounds mb,θb,xb,yb,fxb,fyb,κb... = ȳ
    @inbounds m,fx,fy = pars 
    normf = normfactor(beam)
    _,∂beam = pullback_normfactor([m*mb, fx*fxb,fy*fyb],beam,normf)

    ∂beam += Tangent{B}(l = -node.x/beam.l^2 * xb - node.y/beam.l^2 * yb + only(beam.κ0' * κb),
                       κ0 = vecOrnum(beam.l * κb),
                       θs = θb)

    ∂node = Tangent{No}(;x = xb/beam.l,y=yb/beam.l,ϕ = θb)
    ∂pars =  @SVector [mb * normf * beam.l,fxb * normf * beam.l^2,fyb * normf * beam.l^2]
    
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

function CRC.rrule(::typeof(initialize_beam),node::No,beam::BeamElement{BT},pars::AbstractVector) where{BT,No<:Boundary}
    x,y,θ0, = node.x,node.y,node.ϕ
    l,θs,κ = beam.l,beam.θs,beam.κ0
    m,fx,fy = pars .* normvector(beam) #am Balkenelement
    result = [m,θ0 + θs,x./l,y./l,fx,fy,κ*l]
    init_beam_back(ȳ) = pullback_init_beam(ȳ,node,beam,pars)
    return result, init_beam_back
end 

function init_u0_pullback(ȳ,beams::BT,nodes::NT,pars,nodepos,i,symbolcache) where{BT<:NamedTuple,NT<:NamedTuple}
           
    ∂pars = zero(pars)
    _,dnode,dbeam,dpars = pullback_init_beam(ȳ,nodes[nodepos[i]],beams[i],pars[:,i])
    ∂pars[:,i] += dpars 
    ∂beams = Tangent{typeof(beams)}(;symbolcache[1][i] => dbeam)
    ∂nodes = Tangent{typeof(nodes)}(;symbolcache[2][nodepos[1]] => dnode)
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
    symbolcache = (keys(beams),keys(nodes))  
    init_u0_back(ȳ) = init_u0_pullback(ȳ,beams,nodes,pars,idxs,i,symbolcache)
    result,init_u0_back
end 

function CRC.rrule(::typeof(reduceforceat),node::No,beams,y,idxs) where{No}
    
    sol = reduceforceat(node,beams,y,idxs)
    function force_bound_back(ȳ,∂beams::Tangent{BT} = zero_tangent(beams)) where{BT}

        y_idxs = getforceindices(idxs)
        ∂y = CRC.zero_tangent(y)
        y_ = @view y[y_idxs]
        ∂y_ = @view ∂y[y_idxs] 
        n = 1 
        for b in idxs[1] 
            _,_,dy = pullback_scaleforce(ȳ,y_[:,n],beams[b],Val(false))
            ∂y_[:,n] .+= dy
            n += 1
        end
        for b in idxs[2]
            _,dbeam,dy = pullback_scaleforce(ȳ,y_[:,n],beams[b],Val(true))
            ∂beams += Tangent{BT}(;Symbol(:Beam_,b) => dbeam)
            ∂y_[:,n] .+= dy 
            n += 1 
        end 

        ∂node = ZeroTangent() # CRC.zero_tangent(No)
        return NoTangent(),∂node,∂beams,∂y,NoTangent()
    end 
    return sol,force_bound_back
end   

function CRC.rrule(::typeof(reduceforceat),node::No,beams,y::AbstractArray{T,3},factors,idxs) where{T,No}

    sol = reduceforceat(node,beams,y,factors,idxs)
    function force_bound_back(ȳ,∂beams::Tangent{BT} = zero_tangent(beams)) where{BT}
        
        y_idxs = getforceindices(idxs)
        ∂y = CRC.zero_tangent(y)
        ∂fac = CRC.zero_tangent(factors)
        y_ = @view y[y_idxs]
        ∂y_ = @view ∂y[y_idxs] 
        n = 1
        for b in idxs[1] 
            _,_,dy = pullback_scaleforce(factors[b].*ȳ,y_[:,n],beams[b],Val(false))
            ∂y_[:,n] .+= dy
            n += 1
        end
        for b in idxs[2]
            _,dbeam,dy = pullback_scaleforce(factors[b] .* ȳ,y_[:,n],beams[b],Val(true))
            ∂beams += Tangent{BT}(;Symbol(:Beam_,b) => dbeam)
            ∂y_[:,n] .+= dy 
            n += 1 
        end 
        
        ∂node = ZeroTangent() 
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

function CRC.rrule(::typeof(residuals!),residuals,str::Structure,y,bn)
    ind = 1
    adj = str.AdjMat

    nodes = findall(x->!isapprox(x,0),LowerTriangular(adj))
    branches = filter(x->!isa(bn.Nodes[x],Clamp),eachindex(bn.Nodes))
    positional_nodes = filter(x->any(adj[1:x,x].> 0),axes(adj,1))
    residuals_forces = @view residuals[:,1:length(branches)]
    residuals_positions = @view residuals[:,length(branches)+1:end]
    forces = 1
    positions = 1

    rrforce = Vector{Function}(undef,length(branches))
    rrpos = Vector{Pair{AbstractRange{Int},Function}}()

    for n in getnodeswithbeams(adj,bn.Nodes)
        node = bn.Nodes[n]
        beams = findbeamsatnode(node,n,nodes)
        res, pullback_reduceforceat = rrule(reduceforceat,node,bn.Beams,y,beams)
        if !isempty(res)
            residuals_forces[:,forces] .= res
            rrforce[forces] = pullback_reduceforceat
            forces += 1
        end 
        res, pullback_reduceposat = rrule(reduceposat,node,bn.Beams,y,beams)
        if !isempty(res)
            
            idxs = positions:positions + size(res,2) - 1
            residuals_positions[:,idxs] .= res
            push!(rrpos,idxs =>pullback_reduceposat)
            positions = idxs[end] + 1
        end  
    end
    function residuals!_back(ȳt)
        ȳ = CRC.unthunk(ȳt)
        ∂res = ZeroTangent()
        ∂y = CRC.zero_tangent(y) 
        ȳ_forces = @view ȳ[:,1:length(branches)]
        ȳ_positions = @view ȳ[:,length(branches)+1:end]

        ∂beams = CRC.zero_tangent(bn.Beams)       
        ∂nodes = CRC.zero_tangent(bn.Nodes)
        for (pos,rrforce_rule) in enumerate(rrforce)
            _,dnode,_,dy,_ = rrforce_rule(ȳ_forces[:,pos],∂beams)
            ∂nodes += Tangent{typeof(bn.Nodes)}(;Symbol(:Node_,ind) => dnode)
            ∂y .+= dy 
            
        end
        for (ind,(pos,rrpos_rule)) in enumerate(rrpos)    
            _,dnode,_,dy,_ = rrpos_rule(ȳ_positions[:,pos],∂beams)
            ∂nodes += Tangent{typeof(bn.Nodes)}(;Symbol(:Node_,positional_nodes[ind]) => dnode)
            ∂y .+= dy
        end 
        ∂bn = Tangent{typeof(bn)}(;Beams = ∂beams,Nodes = ∂nodes)

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

    branches = count(x->isa(x,Branch),bn.Nodes)
    positional_nodes = filter(x->any(adj[1:x,x].> 0),axes(adj,1))
    residuals_forces = @view residuals[:,1:branches]
    residuals_positions = @view residuals[:,branches+1:end]

    forces = 1
    positions = 1

    rrforce = Vector{Function}(undef,length(branches))
    rrpos = Vector{Pair{AbstractRange{Int},Function}}()

    for n in getnodeswithbeams(adj,bn.Nodes)#axes(adj,1)
        node = bn.Nodes[n]
        beams = findbeamsatnode(node,n,idcs)
        
        res, pullback_reduceforceat = rrule(reduceforceat,node,bn.Beams,y,adj[idcs],beams)
        if !isempty(res)
            residuals_forces[:,forces] .= res
            rrforce[forces] = pullback_reduceforceat
            forces += 1
        end 
        res, pullback_reduceposat = rrule(reduceposat,node,bn.Beams,y,adj[idcs],beams)
    
        if !isempty(res)
            idxs = positions:positions + size(res,2) - 1
            residuals_positions[:,idxs] .= res
            push!(rrpos,idxs =>pullback_reduceposat)
            positions = idxs[end] + 1
        end  
    end    
    function residuals!_back(ȳ)
        ȳ_ = CRC.unthunk(ȳ)
        ∂res = @thunk(CRC.zero_tangent(residuals))
        ∂y = zero(y)
        ∂fac = zero(adj)
        ∂beams = CRC.zero_tangent(bn.Beams)         
        ∂nodes = CRC.zero_tangent(bn:Nodes) 

        ȳ_forces = @view ȳ_[:,1:branches]
        ȳ_positions = @view ȳ_[:,branches+1:end]


        # for (ind,pos) in resforcedict
        for (pos,rrforce_rule) in enumerate(rrforce)
            _,dnode,_,dy,_ =rrforce_rule[ind](ȳ_forces[:,pos],∂beams)
            ∂nodes += Tangent{typeof(bn.Nodes)}(;Symbol(:Node_,ind) => dnode)
            ∂y .+= dy 
            # ∂fac[idcs] .+= df
        end
        
        # for (ind,pos) in resposdict
        for (ind,(pos,rrpos_rule)) in enumerate(rrpos)    
            _,_,dbeams,dy,_ = rrposdict[ind](ȳ_positions[:,pos],∂beams)
            ∂nodes += Tangent{typeof(bn.Nodes)}(;Symbol(:Node_,ind) => dnode)
            ∂y .+= dy

            # ∂fac[idcs] .+= df
        end 
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
    anz = count(x->!isa(x,Clamp),values(nodes))
    xpos = @view x[:,1:anz]
    xforces = @view x[:,anz+1:end]
    nodes_,back = CRC.rrule(addpositions,nodes,xpos)
    function changeback(ȳ)
        
        ∂x = zero(x)
        ∂x[:,1+anz:end] .+= ȳ[1]
        _,∂nodes,∂xpos = back(ȳ[2])
        ∂x[:,1:anz] .+= ∂xpos
        return NoTangent(),∂nodes,∂x
    end 
    return (xforces,nodes_),changeback
end 

function output_func(u,beams,nodes,nodepos,x,i)
    x_idxs = nodepos[i]
    _,dnode,dbeam,∂xforces = pullback_init_beam(u(0),nodes[x_idxs],beams[i],x[:,i])
    ∂beams = Tangent{typeof(beams)}(;Symbol(:Beam_,i) => dbeam)
    ∂nodes = Tangent{typeof(nodes)}(;Symbol(:Node_,x_idxs) => dnode)
    (∂beams,∂nodes,∂xforces),false
end 

function reduction_func!(u,data,I)
    
    for ((∂beam,∂node,∂xforces),i) in zip(data,I)
        @inbounds u[1][:,i] .+= ∂xforces
        u[2] += ∂beam
        u[3] += ∂node
    end 


    u,false
end 

function u_init(x, bn,dbeams,dnodes)
    dx = similar(x)
    for i in axes(x,2)
        dx[:,i] .= x[:,i] .* normvector(bn.Beams[i])
    end    
    [dx,
    ZeroTangent(),
    ZeroTangent(),
    ]
end 

function initialize_u0_vjp(x::AbstractVector{T},::Beam) where{T}
    x
end 

function initialize_u0_vjp(x::AbstractVector{T},beam::CurvedBeam{T}) where{T}
    u0 = zeros(T,length(x) + length(beam.κ0))
    u0[1:length(x)] .= x
    u0[length(x)+1:end] .= beam.κ0
    u0
end 

function CRC.rrule(str::Structure,x::AbstractArray{T,N},bn::NamedTuple) where{T,N}
    nodepos = getstartnodes(str)

    (xforces,nodes_),change_pullback = CRC.rrule(changestartnodes,bn.Nodes,x)
    function prob_func(prob,i,repeat) 
        u0,p = initialize_beam(bn.Beams,nodes_,xforces,nodepos,i) 
        remake(prob;u0 = u0,p = p)
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
        # tet = CRC.unthunk(∂sol)
        # @show tet[2:4,1,1] ./[1,bn.Beams[1].l,bn.Beams[1].l]
        # @show tet[2:4,1,2] ./[1,bn.Beams[2].l,bn.Beams[2].l]
        # @show tet[2:4,1,3] ./[1,bn.Beams[3].l,bn.Beams[3].l]
        # @show tet[2:4,1,4]./[1,bn.Beams[4].l,bn.Beams[4].l]
        SymbolCache = (keys(bn.Beams),keys(bn.Nodes))

        function prob_func(prob,i,repeat) 
            u0 = initialize_u0_vjp(∂sol[:,2,i],bn.Beams[i])
            remake(prob;u0 = u0,p = solforward[i],tspan = (one(T),zero(T)))
        end
        ensprob  =  EnsembleProblem(vjpprob;prob_func = prob_func,
                                    output_func = (sol,i) -> output_func(sol,bn.Beams,nodes_,nodepos,xforces,i,SymbolCache),
                                    reduction = (u,data,I) -> reduction_func!(u,data,I),
                                    u_init = u_init(∂sol[[1,5,6],1,:],bn,∂Beams,∂Nodes))
        solp = solve(ensprob,str.Solver,
                    EnsembleThreads(),
                    reltol = 1e-6,abstol = 1e-6,
                    trajectories = length(bn.Beams)
                    )
        ∂xforces,dBeams,dNodes = solp

        _,∂nodes_change, ∂x = change_pullback((∂xforces,dNodes))
        # add_startgrads!(∂x,∂sol[[3,4,2],1,:],bn.Beams,bn.Nodes,)
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
        u0,p = initialize_beam(bn.Beams,nodes_,xforces,nodepos,i)
        remake(prob;u0 = u0,p = p)
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
