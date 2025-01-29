struct Beam{T<:Real} 
    l::Union{T,Nothing}
    h::Union{T,Nothing}
    w::Union{T,Nothing}
    κ0::Union{T,Nothing}
    E::Union{T,Nothing}
    θs::Union{T,Nothing}
    θe::Union{T,Nothing}
end 

function Beam(l,h,w,κ0;E = 2.1f5,θs = 0,θe = θs + l*κ0)
    p = promote(l,h,w,κ0,E,θs,θe)
    Beam(p...)
end 
function Beam{T}(l,h,w,κ0;E = 2.1f5,θs = 0f0,θe = θs + l*κ0) where{T}
    Beam{T}(l,h,w,κ0,E,θs,θe)
end 

function Beam{T}(nt::NamedTuple) where{T}
    Beam{T}(nt.l,nt.h,nt.w,nt.κ0,nt.E,nt.θs,nt.θe)
end 
    
function (b::Beam{T})(nt) where{T}
    Beam{T}(nt)
end 

function Base.show(io::IO,beam::Beam)
    return println(io, "Beam with Length: $(beam.l),width: $(beam.w), height: $(beam.h), curvature: $(beam.κ0) and E: $(beam.E)")
end 

function change_beam(beam::Beam;kwargs...)
    for (field,value) in kwargs
        beam = Setfield.@set beam.$field = value 
    end 
    beam
end  

Base.length(b::Beam) = 7
Base.getindex(b::Beam,idx::AbstractVector) = map(x->getfield(b,x),fieldnames(Beam)[idx])
Base.getindex(b::Beam,idx::Int) = getfield(b,fieldnames(Beam)[idx])
Base.iterate(b::Beam,i::Int = 1) = i > 7 ? nothing : (getfield(b,i),i+1)
Base.IteratorSize(b::Beam) = Base.HasLength()

Base.:*(a::Real,b::Beam) = Beam(a .* b...)
Base.:*(b::Beam,a::Real) = Beam(a .* b...)
Base.:*(a::Beam,b::Beam) = Beam(a .* b...)
Base.:-(a::Beam,b::Beam) = Beam(a .- b...)
Base.:-(a::Beam,b::Real) = Beam(a .- b...)
Base.:-(a::Real,b::Beam) = Beam(a .- b...)

Base.:abs2(b::Beam) = b*b
Base.:/(a::Real,b::Beam) = Beam(a ./ b...)
Base.:/(b::Beam,a::Real) = Beam(b ./ a...)
Base.:+(a::Beam,b::Beam) = Beam(a .+ b...)
# Base.:+(a::Beam,b::NamedTuple) = Beam(map(x->isnothing(getfield(b,x)) ? getfield(a,x) : getfield(a,x) + getfield(b,x),keys(b))...)
# Base.:+(b::NamedTuple,a::Beam) = a + b

Optimisers.functor(b::Beam{T}) where{T} = (NamedTuple{fieldnames(Beam)}(b[1:7]),Beam{T})
Optimisers.init(o::Adam, x::Beam{T}) where{T} = (Beam{T}(zeros(T,7)...),Beam{T}(zeros(T,7)...), T.(o.beta))
Optimisers.init(o::WeightDecay, x::Beam) = nothing
Optimisers.isnumeric(::Beam) = true
Optimisers.subtract!(a::Beam{T},b::Beam) where{T} = Beam{T}(a .- Beam{T}(merge(Optimisers.mapvalue(_->zero(T),Optimisers.functor(b)[1]),Optimisers.trainable(b)))...)
Optimisers.trainable(b::Beam) = (;l = b.l,h = b.h,w = b.w,κ0 = b.κ0,E = b.E,θs = b.θs,θe = b.θe)

Optimisers._trainable(b::Beam{T},fr) where{T} = Beam{T}(merge(Optimisers.mapvalue(_ -> nothing, Optimisers.functor(b)[1]), Optimisers.trainable(b)))

Base.zero(::Beam{T}) where{T} = Beam(zeros(T,7)...)

@generated function combine(β,mt::Beam{T},dx) where{T}
    fields = fieldnames(mt)
    # exprs = [:(β * getproperty(mt, $(QuoteNode(f))) + (1-β) * getproperty(dx, $(QuoteNode(f)))) for f in fields]
    exprs = [:(isnothing(getproperty(dx, $(QuoteNode(f)))) ? zero(T) : β * getproperty(mt, $(QuoteNode(f))) + (1-β) * getproperty(dx, $(QuoteNode(f)))) for f in fields]
    return quote
        $(Expr(:call,Beam, exprs...))
    end
end 

@generated function combineabs2(β,mt::Beam{T},dx) where{T}
    fields = fieldnames(mt)
    exprs = [:(isnothing(getproperty(dx, $(QuoteNode(f)))) ? zero(T) : β * getproperty(mt, $(QuoteNode(f))) + (1-β) * abs2(getproperty(dx, $(QuoteNode(f))))) for f in fields]
    return quote
        $(Expr(:call,Beam, exprs...))
    end
end 

@generated function combine(η,βt,mt::B,vt::B,ϵ) where{B<:Beam}
    fields = fieldnames(mt)
    exprs = [:(getproperty(mt, $(QuoteNode(f))) / (1-βt[1]) / (sqrt(getproperty(vt, $(QuoteNode(f))) / (1 -βt[2]))+ ϵ)* η) for f in fields]
    return quote
        $(Expr(:call,B, exprs...))
    end
end 

function Optimisers.apply!(o::Adam,state,b::Beam{T},dx) where{T}
    η, β, ϵ = T(o.eta), T.(o.beta), T(o.epsilon)
    mt, vt, βt = state
    
    mt = combine(β[1],mt,dx)
    vt = combineabs2(β[2],vt,dx)
    dx′ = combine(η,βt,mt,vt,ϵ) #  mt / (1 - βt[1]) / (sqrt(vt / (1 - βt[2])) + ϵ) * η
    return (mt, vt, βt .* β), dx′
end 

function Optimisers.apply!(o::WeightDecay, state, x::Beam{T}, dx) where{T}
    λ = T(o.lambda)
    dx′ = dx + λ * x
  
    return state, Beam{T}(dx′...)
end

Optimisers.init(o::OptimiserChain, x::Beam) = map(opt -> Optimisers.init(opt, x), o.opts)

Optimisers.init(o::Optimisers.ClipNorm, x::Beam) = nothing

function Optimisers.apply!(o::Optimisers.ClipNorm, state, x::Beam{T}, dx) where T
  nrm = norm(dx, o.p)
  if o.throw && !isfinite(nrm)
    throw(DomainError("gradient has $(o.p)-norm $nrm, for array $(summary(x))"))
  end
  λ = T(min(o.omega / nrm, 1))
#   @show Beam{T}(dx * λ...)
  return state, λ * Beam{T}(dx) 
end

function ode!(dT,t::AbstractVector{T},p,s) where{T}
    @inbounds m,θ,x,y,fx,fy,κ = t
    dT[1] = fx*sin(θ)-fy*cos(θ)   #dM
    dT[2] = m + κ               #dΘ
    dT[3] = cos(θ) #* (1 + T[5]*h̃^2/(12)) #dx
    dT[4] = sin(θ) #* (1 + T[6]*h̃^2/(12)) #dy
    dT[5] = zero(T)
    dT[6] = zero(T)
    dT[7] = zero(T)
    return dT
end 

function ode(t::AbstractVector{T},p,s) where{T}
    @inbounds m,θ,x,y,fx,fy,κ = t
    [fx*sin(θ)-fy*cos(θ),   #dM
    m + κ,               #dΘ
    cos(θ), #* (1 + T[5]*h̃^2/(12)) #dx
    sin(θ), #* (1 + T[6]*h̃^2/(12)) #dy
    zero(T),
    zero(T),
    zero(T)]
end 

function jac(t::AbstractArray{T,N},p,s) where{T,N} #jacobi 
    @inbounds m,θ,x,y,fx,fy,κ = t
    dT = zeros(T,7,7)
    dT[1,2] = fx*cos(θ) + fy*sin(θ)
    dT[1,5] = sin(θ)
    dT[1,6] = -cos(θ)
    dT[2,1] = one(T)
    dT[2,7] = one(T)
    dT[3,2] = -sin(θ)
    dT[4,2] = cos(θ)
    return dT
end 

function jac!(dt,t::AbstractArray{T,N},p,s) where{T,N} #jacobi 
    @inbounds m,θ,x,y,fx,fy,κ = t
    dt[1,2] = fx*cos(θ) + fy*sin(θ)
    dt[1,5] = sin(θ)
    dt[1,6] = -cos(θ)
    dt[2,1] = one(T)#
    dt[2,7] = one(T)#
    dt[3,2] = -sin(θ)
    dt[4,2] = cos(θ)
    return dt
end
# vjp!(Jv,v,p,t) = BeamStructures.vjp!(Jv,v,sol.u[end],p,t)
function vjp!(Jv,λ::AbstractArray{T,N},u,t) where{T,N} #vjp
    @inbounds δm,δθ,δx,δy,δfx,δfy,δκ = λ
    @inbounds m,θ,x,y,fx,fy,κ = u
    Jv[1] = -δθ 
    Jv[2] = -δy * cos(θ) - δx * -sin(θ) - δm * (fx*cos(θ) + fy*sin(θ)) 
    # Jv[3] = zero(T) 
    # Jv[4] = zero(T)
    Jv[5] = -δm  *  sin(θ)
    Jv[6] = -δm  * -cos(θ) 
    Jv[7] = -δθ 
    # Jv .+= 
    return nothing 
end
function vjp!(Jv,λ::AbstractArray{T,N},u::ODESolution,t) where{T,N} #vjp
    @inbounds δm,δθ,δx,δy,δfx,δfy,δκ = λ
    @inbounds m,θ,x,y,fx,fy,κ = u(t)
    Jv[1] = -δθ 
    Jv[2] = -δy * cos(θ) - δx * -sin(θ) - δm * (fx*cos(θ) + fy*sin(θ)) 
    # Jv[3] = zero(T) 
    # Jv[4] = zero(T)
    Jv[5] = -δm  *  sin(θ)
    Jv[6] = -δm  * -cos(θ) 
    Jv[7] = -δθ 

    return nothing 
end

func = ODEFunction{true}(ode!, jac = jac!)

prob = ODEProblem(func,zeros(Float32,7),(0f0,1f0))

vjpfunc = ODEFunction(vjp!)
vjpprob = ODEProblem(vjpfunc,ones(Float32,7),(1f0,0f0))

