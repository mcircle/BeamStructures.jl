
abstract type Boundary end

struct Clamp{T} <:Boundary
    x::T
    y::T
    ϕ::T
    fx::T
    fy::T
    mz::T
end 


function Clamp(x::X,y::Y,ϕ::P,fx::FX,fy::FY,mz::MZ) where {X,Y,P,FX,FY,MZ}
    p = promote(x,y,ϕ,fx,fy,mz)
    Clamp(p...)
end  

struct Branch{T} <:Boundary
    x::T
    y::T
    ϕ::T
    fx::T
    fy::T
    mz::T
end

function Branch(x::X,y::Y,ϕ::P,fx::FX,fy::FY,mz::MZ) where {X,Y,P,FX,FY,MZ}
    p = promote(x,y,ϕ,fx,fy,mz)
    Branch(p...)
end  

struct Free{T} <:Boundary
    x::T
    y::T
    ϕ::T
    fx::T
    fy::T
    mz::T
end

function Free(x::X,y::Y,ϕ::P,fx::FX,fy::FY,mz::MZ) where {X,Y,P,FX,FY,MZ}
    p = promote(x,y,ϕ,fx,fy,mz)
    Free(p...)
end  

struct ExtForces{T}<:Boundary
    x::T
    y::T
    ϕ::T
    fx::T
    fy::T
    mz::T
end 

function ExtForces(x::X,y::Y,ϕ::P,fx::FX,fy::FY,mz::MZ) where {X,Y,P,FX,FY,MZ}
    p = promote(x,y,ϕ,fx,fy,mz)
    ExtForces(p...)
end  

struct Slider{T} <:Boundary
    x::T
    y::T
    ϕ::T
    fx::T
    fy::T
    mz::T
end

function Slider(x::X,y::Y,ϕ::P,fx::FX,fy::FY,mz::MZ) where {X,Y,P,FX,FY,MZ}
    p = promote(x,y,ϕ,fx,fy,mz)
    Slider(p...)
end  

Base.length(b::Boundary) = 6
Base.getindex(b::Boundary,idx::AbstractVector) = map(x->getfield(b,x),fieldnames(Clamp)[idx])
Base.getindex(b::Boundary,idx::Int) = getfield(b,fieldnames(Clamp)[idx])

struct CompliantClamp{T} <:Boundary
    x::T
    y::T
    ϕ::T
    fx::T
    fy::T
    c::T
end

function CompliantClamp(x::X,y::Y,ϕ::P,fx::FX,fy::FY,c::C) where {X,Y,P,FX,FY,C}
    p = promote(x,y,ϕ,fx,fy,c)
    CompliantClamp(p...)
end  

struct Beam{T} 
    l::T
    h::T
    w::T
    κ0::T
    E::T
    θs::T
    θe::T
end 

function Beam(l,h,w,κ0;E = 2.1e5,θs = 0.,θe = θs + l*κ0)
    p = promote(l,h,w,κ0,E,θs,θe)
    Beam(p...)
end 

function Base.show(io::IO,beam::Beam)
    return println(io, "Beam with Length: $(beam.l),width: $(beam.w), height: $(beam.h), curvature: $(beam.κ0) and E: $(beam.E)")
end 

Base.length(b::Beam) = 7
Base.getindex(b::Beam,idx::AbstractVector) = map(x->getfield(b,x),fieldnames(Beam)[idx])
Base.getindex(b::Beam,idx::Int) = getfield(b,fieldnames(Beam)[idx])
