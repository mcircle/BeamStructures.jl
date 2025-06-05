
# abstract type AbstractBoundary{T} end
abstract type Boundary{T} end #extern
# abstract type Movable{A<:Real} <:AbstractBoundary{A} end 

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
(cl::Clamp{T})(displacement) where{T} = cl


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
    x::A
    y::A
    dir::A
    ϕ::A
    s::A
    fx::A
    fy::A
    mz::A
    function LinearSlider(x::T,y::T,dir::T,ϕ::T,fx::T,fy::T,mz::T) where{T}
        new{T}(x,y,dir,ϕ,zero(T))
    end 
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

# Eine Funktion f(s) -> x,y,ϕ,trans notwendig
mutable struct Movable{A<:Real} <:Boundary{A}
    x::A
    y::A
    ϕ::A
    trans::AbstractMatrix{A}
    fx::A
    fy::A
    mz::A
    function Movable(x::T,y::T,ϕ::T,fx::T,fy::T,mz::T) where{T}
        new{T}(x,y,ϕ,Diagonal(ones(T,3)),fx,fy,mz)
    end 
end


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


Base.length(::B) where{B<:Boundary} = fieldcount(B)
Base.lastindex(::B) where{B<:Boundary} = fieldcount(B)
Base.getindex(b::Boundary,idx::AbstractVector) = map(x->getfield(b,x),fieldnames(Clamp)[idx])
Base.getindex(b::Boundary,idx::Int) = getfield(b,fieldnames(Clamp)[idx])
Base.iterate(b::B,i::Int = 1) where{B<:Boundary} = i > length(b) ? nothing : (getfield(b,i),i+1)
Base.IteratorSize(b::T) where{T<:Boundary} = Base.HasLength()

function (b::Type{B})(x::NamedTuple) where{B<:Boundary}
    B(map(k->getfield(x,k),fieldnames(b))...)
end 

function (b::Type{B})(x::NamedTuple) where{T,B<:Boundary{T}}
    B(map(k->getfield(x,k),fieldnames(b))...)
end 

function change_node(node::Boundary;kwargs...)
    for (field,value) in kwargs
        node = Setfield.@set node.$field = value 
    end 
    node
end  

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

function Base.:+(a::T,b::NamedTuple) where{T<:Boundary} 
    T(map(x->hasproperty(b,x) ? getfield(a,x) + getfield(b,x) : getfield(a,x),fieldnames(T))...)
end
# function CRC.rrule(::typeof(Base.:+), a::T, b::NamedTuple) where{T<:Boundary}
#     ab = T(map(x->hasproperty(b,x) ? getfield(a,x) + getfield(b,x) : getfield(a,x),fieldnames(T))...), 
#     function back_add_nt(dy)
#         da = Tangent{T}(map(x->hasproperty(b,x) ? x -> getfield(dy,x) : nothing,fieldnames(T))...)
#         db = Tangent{typeof(b)}(map(x->hasproperty(b,x) ? x -> getfield(dy,x) : nothing,fieldnames(T))...)
#         return CRC.NoTangent(),da,db
#     end
#     return ab,back_add_nt
# end

Optimisers.functor(x::T) where{T<:Boundary} = (NamedTuple{fieldnames(T)}(x[1:end]),T)
Optimisers.init(o::Adam, x::B) where{B<:Boundary{T}} where{T}  = (B(zeros(T,6)...), B(zeros(T,6)...), T.(o.beta))
Optimisers.init(o::WeightDecay, x::Boundary) = nothing
Optimisers.isnumeric(x::T) where{T<:Boundary} = true
Optimisers.trainable(c::Clamp) = (;x = c.x,y = c.y,ϕ = c.ϕ)
Optimisers.trainable(c::Branch) = (;x = c.x,y = c.y,ϕ = c.ϕ)
Optimisers.subtract!(a::T,b::T) where{T<:Boundary{TT}} where{TT} = a - T(merge(Optimisers.mapvalue(_->zero(TT),Optimisers.functor(b)[1]),Optimisers.trainable(b)))
Optimisers.init(o::OptimiserChain, x::Boundary) = map(opt -> Optimisers.init(opt, x), o.opts)
Optimisers._trainable(b::T,fr) where{T<:Boundary} =T(merge(Optimisers.mapvalue(_ -> nothing, Optimisers.functor(b)[1]), Optimisers.trainable(b)))



@generated function combine(β,mt::B,dx) where{T,B<:Boundary{T}}
    fields = fieldnames(mt)
    exprs = [:(isnothing(getproperty(dx, $(QuoteNode(f)))) ? zero(T) : β * getproperty(mt, $(QuoteNode(f))) + (1-β) * getproperty(dx, $(QuoteNode(f)))) for f in fields]
    return quote
        $(Expr(:call, :B, exprs...))
    end
end 

@generated function combineabs2(β,mt::B,dx) where{T,B<:Boundary{T}}
    fields = fieldnames(mt)
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
    
    mt =combine(β[1],mt,dx) # β[1] * mt + (1 - β[1]) * dx
    vt = combineabs2(β[2],vt,dx) # β[2] * vt + (1 - β[2]) * abs2(dx)
    dx′ = combine(η,βt,mt,vt,ϵ) # mt / (1 - βt[1]) / (sqrt(vt / (1 - βt[2])) + ϵ) * η
  
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