
abstract type Boundary{T} end #extern


struct Clamp{A<:Real} <:Boundary{A}
    x::A
    y::A
    ϕ::A
    fx::A
    fy::A
    mz::A
    function Clamp(x::T,y::T,ϕ::T,fx::T,fy::T,mz::T) where{T<:Real}
        new{T}(x,y,ϕ,fx,fy,mz)
    end 
    Clamp{T}(x,y,ϕ,fx,fy,mz) where{T<:Real} = new{T}(x,y,ϕ,fx,fy,mz)
end 
Clamp(x,y,ϕ,fx,fy,mz) = Clamp(promote(x,y,ϕ,fx,fy,mz)...)
type(::Clamp{T}) where{T} = Clamp
type(::Type{Clamp{T}}) where{T} = Clamp

struct Branch{A<:Real} <:Boundary{A}
    x::A
    y::A
    ϕ::A
    fx::A
    fy::A
    mz::A
    function Branch(x::T,y::T,ϕ::T,fx::T,fy::T,mz::T) where{T<:Real}
        new{T}(x,y,ϕ,fx,fy,mz)
    end 
    Branch{T}(x,y,ϕ,fx,fy,mz) where{T<:Real} = new{T}(x,y,ϕ,fx,fy,mz)
end
(cl::Branch{T})(displacement) where{T} = cl
Branch(x,y,ϕ,fx,fy,mz) = Branch(promote(x,y,ϕ,fx,fy,mz)...)
type(::Branch{T}) where{T} = Branch
type(::Type{Branch{T}}) where{T} = Branch

struct Free{A<:Real} <:Boundary{A}
    x::A
    y::A
    ϕ::A
    fx::A
    fy::A
    mz::A
    function Free(x::T = 0,y::T = zero(T),ϕ::T = zero(T),fx::T = zero(T),fy::T = zero(T),mz::T = zero(T)) where{T<:Real}
        new{T}(x,y,ϕ,fx,fy,mz)
    end 
    Free{T}(x,y = zero(T),ϕ = zero(T),fx = zero(T),fy = zero(T),mz = zero(T)) where{T<:Real} = new{T}(x,y,ϕ,fx,fy,mz)
end
(cl::Free{T})(displacement) where{T} = cl

struct ExtForces{A<:Real} <: Boundary{A}
    x::A
    y::A
    ϕ::A
    fx::A
    fy::A
    mz::A
    function ExtForces(x::T,y::T,ϕ::T,fx::T,fy::T,mz::T) where{T}
        new{T}(x,y,ϕ,fx,fy,mz)
    end 
end
(cl::ExtForces{T})(displacement) where{T} = cl[1:6]

struct LinearSlider{A<:Real} <:Boundary{A}
    x::A # x0
    y::A # y0
    ϕ::A 
    fx::A
    fy::A
    mz::A
    function LinearSlider(x::T,y::T,ϕ::T,fx::T,fy::T,mz::T) where{T}
        new{T}(x,y,ϕ,fx,fy,mz)
    end 
end

function LinearSlider{T}(x,y,ϕ,fx,fy,mz) where{T<:Real} 
    p = T.([x,y,ϕ,fx,fy,mz])
    LinearSlider(p...)
end 
struct Joint{A<:Real} <:Boundary{A}
    x::A
    y::A
    ϕ::A
    fx::A
    fy::A
    function Joint(x::T,y::T,ϕ::T,fx::T,fy::T,mz::T) where{T}
        new{T}(x,y,ϕ,fx,fy,mz)
    end 
end

# x,y,phi geben SollPosition an 
mutable struct Movable{A<:Real} <:Boundary{A}
    x::A
    y::A
    ϕ::A
    fx::A
    fy::A
    mz::A
    trans::AbstractVector{Int}
    function Movable(x::T,y::T,ϕ::T,fx::T,fy::T,mz::T,dir) where{T}
        new{T}(x,y,ϕ,fx,fy,mz,dir)
    end 
    function Movable{T}(x,y,ϕ,fx,fy,mz,dir = ones(Int,3)) where{T<:Real} 
        new{T}(T(x),T(y),T(ϕ),T(fx),T(fy),T(mz),dir)
    end 
end

function Movable(x,y,ϕ,fx,fy,mz)
    p = promote(x,y,ϕ,fx,fy,mz)
    Branch{eltype(p)}(p...)
end 
# function Movable(x,y,ϕ,dir,fx,fy,mz,dir) 
#     p = promote(x,y,ϕ,fx,fy,mz)
#     Movable{eltype(p)}(p...,,dir)
# end 

struct CompliantClamp{A} <:Boundary{A}
    x::A
    y::A
    ϕ::A
    cx::A
    cy::A
    cz::A
    function CompliantClamp(x::T,y::T,ϕ::T,cx::T,cy::T,cz::T) where{T}
        new{T}(x,y,ϕ,cx,cy,cz)
    end 
end
(cl::CompliantClamp{T})(displacement) where{T} = cl[1:6]

gettype(::B) where{T,B<:Boundary{T}} = B

(::Type{B})(x::AbstractVector) where{T,B<:Boundary{T}} = B(T.(x)...)
(::Type{B})(x::AbstractVector) where{B<:Boundary} = B(x...)


function Base.getindex(b::B, i::Int) where B<:Boundary
    fields = fieldnames(B)
    N = fieldcount(B)
    if i > N
        throw(BoundsError(b, i))
    end
    getfield(b, fields[i])
end

function Base.iterate(b::B, i::Int=1) where B<:Boundary
    N = fieldcount(B)
    i> N && return nothing
    fields = fieldnames(B)
    return getfield(b, fields[i]), i+1
end

Base.length(::B) where{B<:Boundary} = fieldcount(B)
Base.lastindex(::B) where{B<:Boundary} = fieldcount(B)

# Base.iterate(b::B,i::Int = 1) where{B<:Boundary} = i > length(b) ? nothing : (getfield(b,i),i+1)
Base.getindex(b::B,i::AbstractVector) where{B<:Boundary} = map(j->getindex(b,j),i)
Base.IteratorSize(b::T) where{T<:Boundary} = Base.HasLength()

function (b::Type{B})(x::NamedTuple) where{B<:Boundary}
    B(map(k->getfield(x,k),fieldnames(b))...)
end 

function (b::Type{B})(x::NamedTuple{kwargs}) where{T,kwargs,B<:Boundary{T}}
    
    B(map(k->getfield(x,k),fieldnames(b))...)
end 

# function change_node(node::Boundary;kwargs...)
#     for (field,value) in kwargs
#         node = Setfield.@set node.$field = value 
#     end 
#     node
# end  

forcesatnode(::Boundary) = true
forcesatnode(::Clamp) = false

canchangeposition(::Boundary) = true
canchangeposition(::Clamp) = false

degreesofreedom(::Boundary) = [1,1,1]
degreesofreedom(::Clamp) = [0,0,0]
degreesofreedom(node::LinearSlider) = [0,cos(node.ϕ % π/2),sin(node.ϕ % π/2)]

Base.zero(::B) where{T,B<:Boundary{T}} = B(zeros(T,6)...)

Base.:*(a::Real,b::T) where{T<:Boundary} = T(a .* b...)
Base.:*(b::T,a::Real) where{T<:Boundary} = T(a .* b...)
Base.:*(a::T,b::T) where{T<:Boundary} = T(a .* b...)
Base.:+(a::T,b::T) where{T<:Boundary} = T(a .+ b...)
Base.:-(a::T,b::T) where{T<:Boundary}= T(a .- b...)
Base.:-(a::T,b::Real) where{T<:Boundary} = T(a .- b...)
Base.:-(a::Real,b::T) where{T<:Boundary} = T(a .- b...)

Base.:abs2(b::T) where{T<:Boundary} = b*b
Base.:/(a::Real,b::T) where{T<:Boundary} = T(a ./ b...)
Base.:/(b::T,a::Real)where{T<:Boundary} = T(b ./ a...)


@generated function Base.:+(a::B,b::NamedTuple{V,Tu}) where{T,N,T2,V,Tu<:NTuple{N,T2},B<:Boundary{T}} 
    # println(Tu)
    fn = fieldnames(B)
    t = promote_type(T, T2)
    a_ = type(a)
    exprs = [:( $(QuoteNode(f)) in $V ? getproperty(a, $(QuoteNode(f))) + $t(getproperty(b, $(QuoteNode(f)))) : getproperty(a, $(QuoteNode(f))) ) for f in fn]
    return quote
        $(Expr(:call, a_, exprs...))
    end
    # B( map(f -> hasproperty(b,f) ? getproperty(a,f) + getproperty(b,f) : getproperty(a,f), fn)... )
end

Base.zero(::Type{T}) where{T<:Boundary} = T(zeros(eltype(T),6)...)
Base.promote_rule(::Type{T}, ::Type{S}) where{T<:Boundary{TT},S<:Boundary{SS}} where{TT,SS} = T, S, promote_type(TT,SS)

Optimisers.functor(x::T) where{T<:Boundary} = (NamedTuple{fieldnames(T)}(x[1:end]),T)
Optimisers.init(o::Adam, x::B) where{B<:Boundary{T}} where{T}  = (B(zeros(T,6)...), B(zeros(T,6)...), T.(o.beta))
Optimisers.init(o::WeightDecay, x::Boundary) = nothing
Optimisers.isnumeric(x::T) where{T<:Boundary} = true
Optimisers.trainable(c::Clamp) = (;x = c.x,y = c.y,ϕ = c.ϕ)
Optimisers.trainable(c::Branch) = (;x = c.x,y = c.y,ϕ = c.ϕ)
Optimisers.subtract!(a::T,b::T) where{T<:Boundary{TT}} where{TT} = a - T(merge(Optimisers.mapvalue(_->zero(TT),Optimisers.functor(b)[1]),Optimisers.trainable(b)))
Optimisers.init(o::OptimiserChain, x::Boundary) = map(opt -> Optimisers.init(opt, x), o.opts)
Optimisers._trainable(b::T,fr) where{T<:Boundary} =T(merge(Optimisers.mapvalue(_ -> nothing, Optimisers.functor(b)[1]), Optimisers.trainable(b)))

BOUNDARYSCALE = (
    x  = 1,      # 1 m bleibt 1
    y  = 1,
    ϕ  = 0,      # Winkel ϕ 10x „langsamer“ updaten
    fx = 1,      # ggf. eigene Wahl
    fy = 1,
    mz = 1,
)


@generated function scale_boundary(dx,::B, scale) where{B<:Boundary} 
    fields = fieldnames(dx)
    exprs = [:( getproperty(dx, $(QuoteNode(f))) * getproperty(scale,$(QuoteNode(f)) )) for f in fields]
    return quote
        $(Expr(:call, :B, exprs...))
    end
end

@generated function invscale_boundary(dx,::B, scale) where{B<:Boundary}
    fields = fieldnames(dx)
    exprs = [:( getproperty(dx, $(QuoteNode(f))) / getproperty(scale, $(QuoteNode(f)))) for f in fields]
    return quote
        $(Expr(:call, :B, exprs...))
    end
end

@generated function combine(β,mt::B,dx) where{T,B<:Boundary{T}}
    fields = fieldnames(dx)
    exprs = [:(isnothing(getproperty(dx, $(QuoteNode(f)))) ? zero(T) : β * getproperty(mt, $(QuoteNode(f))) + (1-β) * getproperty(dx, $(QuoteNode(f)))) for f in fields]
    return quote
        $(Expr(:call, :B, exprs...))
    end
end 

@generated function combineabs2(β,mt::B,dx) where{T,B<:Boundary{T}}
    fields = fieldnames(dx)
    exprs = [:(isnothing(getproperty(dx, $(QuoteNode(f)))) ? zero(T) : β * getproperty(mt, $(QuoteNode(f))) + (1-β) * abs2(getproperty(dx, $(QuoteNode(f))))) for f in fields]
    return quote
        $(Expr(:call, :B, exprs...))
    end
end 

@generated function combine(η,βt,mt::B,vt::B,ϵ) where{B<:Boundary}
    fields = fieldnames(mt)
    exprs = [:(getproperty(mt, $(QuoteNode(f))) / (1-βt[1]) / (sqrt(getproperty(vt, $(QuoteNode(f))) / (1 -βt[2]))+ ϵ)* η) for f in fields]
    return quote
        $(Expr(:call,B, exprs...))
    end
end



function Optimisers.apply!(o::Adam,state,b::BT,dx) where{BT<:Boundary{T}} where{T}
    η, β, ϵ = T(o.eta), T.(o.beta), T(o.epsilon)
    mt, vt, βt = state
    
    # dx_scaled = scale_boundary(dx,mt,BOUNDARYSCALE)
    mt =combine(β[1],mt,dx) # β[1] * mt + (1 - β[1]) * dx
    vt = combineabs2(β[2],vt,dx) # β[2] * vt + (1 - β[2]) * abs2(dx)
    dx_scaled = combine( η,βt,mt,vt,ϵ) # mt / (1 - βt[1]) / (sqrt(vt / (1 - βt[2])) + ϵ) * η
     # 3) Schritt zurück in physikalische Einheiten transformieren
     dx′ = scale_boundary(dx_scaled,mt,BOUNDARYSCALE)
    return (mt, vt, βt .* β), dx′
end 

Optimisers.init(o::AdamW, x::Bo) where {T,Bo<:Boundary{T}} = (Bo(zeros(T,6)...), Bo(zeros(T,6)...), T.(o.beta))

function Optimisers.apply!(o::AdamW, state, x::Bo, dx) where {T,Bo<:Boundary{T}}
    η, β, ϵ, λ = T(o.eta), T.(o.beta), T(o.epsilon), T(o.lambda)
    mt, vt, βt = state
  
    # standard Adam update with learning rate eta=1
    mt = combine(β[1],mt,dx)
    vt = combineabs2(β[2],vt,dx)
    dx′ = combine(η,βt,mt,vt,ϵ)
  
    # apply learning rate and weight decay
    if o.couple
      dx′′ =  η * (dx′ + λ * x)
    else
      dx′′ =  η * dx′ + λ * x
    end
  
    return (mt, vt, βt .* β), dx′′
end

Optimisers.init(o::Optimisers.ClipNorm, x::Boundary) = nothing

function Optimisers.apply!(o::Optimisers.ClipNorm, state, x::BT, dx) where{BT<:Boundary{T}} where T
  nrm = norm(dx, o.p)
  if o.throw && !isfinite(nrm)
    throw(DomainError("gradient has $(o.p)-norm $nrm, for array $(summary(x))"))
  end
  λ = T(min(o.omega / nrm, 1))

  return state,  BT(dx * λ...)
end