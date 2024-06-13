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

function ode!(dT,T,p,s)
    m,θ,x,y,fx,fy,κ = T
    dT[1] = fx*sin(θ)-fy*cos(θ)   #dM
    dT[2] = m + κ               #dΘ
    dT[3] = cos(θ) #* (1 + T[5]*h̃^2/(12)) #dx
    dT[4] = sin(θ) #* (1 + T[6]*h̃^2/(12)) #dy
    dT[5] = 0.
    dT[6] = 0.
    dT[7] = 0.
end 

function jac!(dT,T,p,s) #jacobi 
    m,θ,x,y,fx,fy,κ = T
    dT[1,2] = fx*cos(θ) + fy*sin(θ)
    dT[1,5] = sin(θ)
    dT[1,6] = -cos(θ)
    dT[2,1] = 1 
    dT[2,7] = 1
    dT[3,2] = -sin(θ)
    dT[4,2] = cos(θ)
end 

func = ODEFunction(ode!, jac = jac!)

prob = ODEProblem(func,zeros(Float64,7),(0.,1.))

