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

func = ODEFunction(ode!)

prob = ODEProblem(func,zeros(Float64,7),(0.,1.))