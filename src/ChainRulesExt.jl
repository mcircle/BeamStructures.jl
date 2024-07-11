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
            dbeam = Beam{T}(sum(ȳ[2:3] .* y[2:3]),zeros(T,5)...,-ȳ[1])
        else
            dbeam = Beam{T}(sum(ȳ[2:3] .* y[2:3]),zeros(T,4)...,-ȳ[1],zero(T))
        end 
        return NoTangent(),dbeam, dy,NoTangent() 
    end 
    return sol,scalepos_back
end 

function CRC.rrule(::typeof(reduceposat),node::Clamp{T},beams,y,beamnbrs) where{T} 
    sol = reduceposat(node,beams,y,beamnbrs)
    function reduceposat_back(ȳ)
        
        # indbeam = Symbol.(:Beam_,beamnbrs[1])
        # beam_l = getfield.(values(beams[indbeam]),:l)
        rȳ = reshape(ȳ,3,:)
        
        ∂beams= Vector{Beam}(undef,length(beams))
        ∂y = zero(y)
        Tb = eltype(beams)
        fill!(∂beams,Tb(zeros(Int,7)...))
        for (n,b) in enumerate(beamnbrs[1])
            _, spback = CRC.rrule(scalepos,beams[b],y[2:4,2,b],Val(2))
            _,dbeam,dy,_ = spback(rȳ[:,n])    
            ∂beams[b] -= dbeam
            ∂y[2:4,2,b] .-= dy
        end
        dnode = Clamp(reduce((init,x) -> init .+ rȳ[[2,3,1],x],1:length(beamnbrs[1]),init = zeros(T,3))...,zeros(T,3)...)
        return NoTangent(),dnode,∂beams,∂y,NoTangent()
    end 
    return sol, reduceposat_back
end 

function CRC.rrule(::typeof(reduceposat),node::Boundary,beams,y,beamnbrs) 
    bnbrs = vcat(beamnbrs...)
    idxs = getpositionindices(beamnbrs)
    pos = y[idxs]
    scaledpos = map((x,y)->scalepos(beams[x],y,Val(getside(x,beamnbrs))),bnbrs,eachcol(pos))
    res = map(y -> scaledpos[1].- y,scaledpos[2:end])
    
    function reduceposBound_back(ȳ)
        indbeam = Symbol.(:Beam_,bnbrs)
        beam_l = getfield.(values(beams[indbeam]),:l)
        rȳ = reshape(ȳ,3,:)
        
        dbeams = Vector{Beam}(undef,length(beams))
        Tb = eltype(beams)        
        fill!(dbeams,Tb(zeros(Int,7)...))
        
        ∂y = zero(y)
        b = bnbrs[1]
        side = getside(b,beamnbrs)
        _,spback = CRC.rrule(scalepos,beams[b],y[2:4,side,b],Val(side))
        _,dbeam,dy,_ = spback(sum(rȳ,dims = 2))
        dbeams[b] += dbeam
        ∂y[2:4,side,b] .+= dy 
        for (n,b) in enumerate(bnbrs[2:end])
            side = getside(b,beamnbrs)
            _,spback = CRC.rrule(scalepos,beams[b],y[2:4,side,b],Val(side))
            _,dbeam,dy,_ = spback(-rȳ[:,n])
            dbeams[b] += dbeam
            ∂y[2:4,side,b] .+= dy
        end 
        return NoTangent(),NoTangent(),dbeams,∂y,NoTangent()
    end 
    return reducevcat(res),reduceposBound_back
end 

function CRC.rrule(::typeof(normfactor_m), b::Beam{T}) where{T}
    y = normfactor_m(b)
    function pullback(ȳ)
        ∂b = Beam{T}(
            12 * ȳ / (b.E * b.w * b.h^3),
            -36 * b.l * ȳ / (b.E * b.w * b.h^4),
            -12 * b.l * ȳ / (b.E * b.w^2 * b.h^3),
            zero(T),
            -12 * b.l * ȳ / (b.E^2 * b.w * b.h^3),
            zero(T),zero(T)
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
        ∂b = Beam(
            ȳ * normfactor_m(b) + ∂b_m.l,
            ∂b_m.h,
            ∂b_m.w, zero(T),
            ∂b_m.E, zero(T),zero(T)
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

function CRC.rrule(::typeof(reduceforceat),inits,b::Beam,y)
    nv = normvector(b)
    result = inits .+ y ./ nv
    
    function pullback(ȳ)
        ∂y = ȳ ./ nv
        ∂nv = -ȳ .* y ./ nv.^2
        
        _, pb_nv = CRC.rrule(normvector, b)
        _, ∂beam = pb_nv(∂nv)
        ∂inits = ȳ
        return (NoTangent(),∂inits,∂beam,∂y)
    end
    
    return result, pullback
end

function CRC.rrule(::typeof(initialize_beam),node::No,beam::Beam{BT},pars::AbstractVector) where{BT,No<:Boundary}
    x,y,θ0, = node.x,node.y,node.ϕ
    l,θs,κ = beam.l,beam.θs,beam.κ0
    m,fx,fy,Δθ,Δx,Δy = pars
    norm,normback = rrule(normvector,beam)
    m *= norm[1] #am Balkenelement
    fx *= norm[2] #am Balkenelement
    fy *= norm[3] #am Balkenelement
    result = [m,θ0 + θs + Δθ,(x + Δx)./l,(y + Δy)./l,fx,fy,κ*l]
    function init_beam_back(ȳ)
        _,∂beam = normback(pars[1:3] .* ȳ[[1,5,6]])
        ∂l = -result[3]/l * ȳ[3] - result[4]/l * ȳ[4] + κ * ȳ[7]
        ∂beam += Beam{BT}(∂l,zeros(BT,2)...,l*ȳ[7],zero(BT),ȳ[2],zero(BT))
        dn = [1/l,1/l,1] .* ȳ[[3,4,2]]
        ∂node = No(dn...,zeros(eltype(ȳ),3)...)
        ∂pars = [norm .* ȳ[[1,5,6]]...,ȳ[2],ȳ[3]/l,ȳ[4]/l]
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
        ∂l = -result[3]/l * ȳ[3] - result[4]/l * ȳ[4] + κ * ȳ[7]
        ∂beam += Beam{BT}(∂l,zeros(BT,2)...,l*ȳ[7],zero(BT),ȳ[2],zero(BT))

        dn = [1/l,1/l,1] .* ȳ[[3,4,2]]
        ∂node = Clamp{CT}(dn...,zeros(CT,3)...)

        ∂pars = norm .* ȳ[[1,5,6]]
        return NoTangent(),∂node,∂beam,∂pars
    end 
    result,init_back_clamp
end 

function CRC.rrule(::typeof(initialize_beam),bn::NamedTuple,pars,nodes,i::Int)
    x_idxs = get_index_pars(bn.Nodes,nodes[1:i])
    result,init_back = CRC.rrule(initialize_beam,bn.Nodes[nodes[i]],bn.Beams[i],pars[x_idxs])
    function init_u0_back(ȳ)
        dbeams = Vector{Beam}(undef,length(bn.Beams))
        fill!(dbeams,eltype(bn.Beams)(zeros(Int,7)...))

        f(x) = typeof(bn.Nodes[x])(zeros(Int,6)...)
        dnodes = [f.(keys(bn.Nodes))...]
        
        ∂pars = zero(pars)
        _,dnode,dbeam,dpars = init_back(ȳ)
        ∂pars[x_idxs] += dpars 
        dbeams[i] += dbeam
        dnodes[nodes[i]] += dnode
        ∂bn = prepare(dnodes...,dbeams...)
        return NoTangent(), ∂bn,∂pars,NoTangent(),NoTangent()
    end 
    result,init_u0_back
end 

function CRC.rrule(::typeof(reduceforceat),node::No,beams,y,idxs) where{No}
    sol = reduceforceat(node,beams,y,idxs)
    function force_bound_back(ȳ)

        ∂beams = Vector{Beam}(undef,length(beams))
        y_idxs = getforceindices(idxs)
        ∂y = zero(y)
        y_ = @view y[y_idxs]
        ∂y_ = @view ∂y[y_idxs] 
        Tb = eltype(beams)
        fill!(∂beams,Tb(zeros(Int,7)...))
        for (b,yb,dyb) in zip(vcat(idxs...),eachcol(y_),eachcol(∂y_))
            _, rr_norm = CRC.rrule(reduceforceat,zeros(Float32,3),beams[b],yb)
            if b in idxs[2]
                _,_,dbeam,dyrr_ = rr_norm(-ȳ) 
            else
                _,_,dbeam,dyrr_ = rr_norm(ȳ)
            end 
            ∂beams[b] = dbeam
            dyb .= dyrr_
        end 

        ∂node = No(zeros(Int,3)...,ȳ[[3,1,2]]...)
        return NoTangent(),∂node,∂beams,∂y,NoTangent()
    end 
    return sol,force_bound_back
end 

@non_differentiable reduceforceat(n::Clamp,beams,y,idxs) 

function CRC.rrule(::typeof(residuals!),residuals,str,y,bn)
    ind = 1
    adj = str.AdjMat
    nodes = findall(x->!isapprox(x,0),LowerTriangular(adj))
    start = 1
    resposdict = Dict{Int,AbstractRange{Int}}()
    resforcedict = Dict{Int,AbstractRange{Int}}()
    rrposdict = Dict{Int,Function}()
    rrforcedict = Dict{Int,Function}()
    for n in unique(first.(Tuple.(nodes)))
        node = bn.Nodes[n]
        beams = findbeamsatnode(node,n,nodes)
        res, pullback_reduceforceat = rrule(reduceforceat,node,bn.Beams,y,beams)
        if !isempty(res)
            residuals[start:start+2] .= res
            resforcedict[n] = start:start + 2
            rrforcedict[n] = pullback_reduceforceat
            start += 3
        end 
        res, pullback_reduceposat = rrule(reduceposat,node,bn.Beams,y,beams)
        if !isempty(res)
            idxs = range(start,length = length(res))
            resposdict[n] = idxs
            rrposdict[n] = pullback_reduceposat
            residuals[idxs] .= res
            start = idxs[end] + 1
        end  
    end
    function residuals!_back(ȳ)
        ∂res = InplaceableThunk(dself -> dself .+= ȳ,@thunk(copy(ȳ)))
        ∂y = zero(y)
        start = 1
        ∂beams = Vector{Beam}(undef,length(bn.Beams))
        fill!(∂beams,Beam(zeros(Int,7)...))
        
        f(x) = typeof(bn.Nodes[x])(zeros(Int,6)...)
        ∂nodes = [f.(keys(bn.Nodes))...]
        for (ind,pos) in resforcedict
            _,dnode,dbeams,dy,_ =rrforcedict[ind](ȳ[pos])
            ∂nodes[ind] += dnode
            ∂beams .+= dbeams
            ∂y .+= dy
        end
        for (ind,pos) in resposdict
            _,dnode,dbeams,dy,_ = rrposdict[ind](ȳ[pos])
            ∂nodes[ind] += dnode
            ∂beams .+=  dbeams
            ∂y .+= dy
        end 
        ∂bn = prepare(∂beams...,∂nodes...)
        return NoTangent(),∂res,NoTangent(),∂y,∂bn
    end 

    return residuals,residuals!_back
end 
                
function CRC.rrule(::Type{BT}, x,y,θ,fx,fy,mz) where{BT<:Boundary{T}} where{T}
    Bar_pullback(ȳ) = (NoTangent(),BT(ȳ[1:6]...))
    return BT(x,y,θ,fx,fy,mz), Bar_pullback
end
function CRC.rrule(::Type{Clamp{T}}, x,y,θ,fx,fy,mz) where{T}
    Bar_pullback(ȳ) = (NoTangent(),ȳ[1:6]...)
    return Clamp{T}(x,y,θ,fx,fy,mz), Bar_pullback
end

function CRC.rrule(::typeof(+),a::BT,b::BT) where{BT<:Boundary{T}} where{T}
    backaddb(ȳ) = (NoTangent(),ȳ,ȳ)
    return a + b,backaddb
end 