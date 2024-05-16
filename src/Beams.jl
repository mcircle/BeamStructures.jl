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