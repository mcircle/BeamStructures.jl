
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


