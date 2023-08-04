#New_Str1
using Plots
plotlyjs(size = (800,600))
#parameters
DR = 30e-3 #konst
R1 = 20e-3
alpha = 1.6
d = 5e-3
#pre-calculus
L1  = R1 - DR/2
R2  = R1 + d
R3  = R2 + L1
R4  = R3 + d
#changing 
E = 2e11
h = 2e-3
w = 10e-3
Iz = w*h^3/12
ΔFy = 100/(E*Iz)  
ΔFx = 100/(E*Iz)
ΔM = 5/(E*Iz)
#Node's positions
Nx1  =  DR/2
Ny1  =  0.
Nx2  =  R1 
Ny2  =  0.
Nx3  =  R1*cos(alpha/2)# - DR/2
Ny3  =  R1*sin(alpha/2)
Nx4  =  R1*cos(alpha/2)# - DR/2
Ny4  = -R1*sin(alpha/2)
Nx5  =  R2*cos(alpha/2)# - DR/2
Ny5  =  R2*sin(alpha/2)
Nx6  =  R2*cos(alpha/2)# - DR/2
Ny6  = -R2*sin(alpha/2)
Nx7  =  R2# - DR/2
Ny7  =  0
Nx8  =  R3 #- DR/2
Ny8  =  0
Nx9  =  R3*cos(alpha/2) #- DR/2
Ny9  =  R3*sin(alpha/2)
Nx10 =  R3*cos(alpha/2)# - DR/2
Ny10 = -R3*sin(alpha/2)
Nx11 =  R4*cos(alpha/2) #- DR/2
Ny11 =  R4*sin(alpha/2)
Nx12 =  R4*cos(alpha/2)# - DR/2
Ny12 = -R4*sin(alpha/2)
Nx13 =  R4 #- DR/2
Ny13 =  0
ϕ = 0.
Nx14 =  (R4 + L1)*cos(ϕ) #- DR/2
Ny14 =  (R4 + L1)*sin(ϕ)


#               beam       1       2              3       4             5              6           7            8    9             10           11           12              13            14           15       
beamfeatures = collect([  L1   R1*alpha/2   R1*alpha/2    π*d/2         π*d/2         R2*alpha/2  R2*alpha/2    L1   R3*alpha/2    R3*alpha/2   π*d/2        π*d/2           R4*alpha/2    R4*alpha/2   L1 ; #length of the beam 
                           0.  π/2          3π/2          π/2+alpha/2   3π/2-alpha/2  alpha/2-π/2 5π/2-alpha/2  0    π/2           3π/2         π/2+alpha/2  3π/2-alpha/2    alpha/2-π/2   5π/2-alpha/2  0. ; #θs angle as going out slope 
                           0.  1/R1         -1/R1         -2/d          2/d           -1/R2       1/R2          0.   1/R3          -1/R3        -2/d         2/d             -1/R4         1/R4         0. ; #κ as curvatature of the beam 
                           0.  π/2+alpha/2  3π/2-alpha/2  alpha/2-π/2   5π/2-alpha/2  -π/2        5π/2          0    π/2+alpha/2   3π/2-alpha/2 alpha/2-π/2  5π/2-alpha/2    -π/2          5π/2          0. ; #θe angle as going in slope
]')

#nodes as intersections/Connections of the beams

    #            node     1     2      3      4      5      6     7     8    9    10    11    12    13    14
nodefeatures = collect([ 0.     0.     0.     0.     0.     0.    0.    0.   0.   0.    0.    0.    0.    ϕ    ; # Drehung des Nodes um Δθ
                    	 Nx1    Nx2    Nx3    Nx4    Nx5    Nx6   Nx7   Nx8  Nx9  Nx10  Nx11  Nx12  Nx13  Nx14-0.001 ; # position x 
                         Ny1    Ny2    Ny3    Ny4    Ny5    Ny6   Ny7   Ny8  Ny9  Ny10  Ny11  Ny12  Ny13  Ny14 ; # position y 
                         0.     0.     0.     0.     0.     0.    0.    0.   0.   0.    0.    0.    0.    0.   ; # ΔM
                         0.     0.     0.     0.     0.     0.    0.    0.   0.   0.    0.    0.    0.    0.   ; # ΔFx
                         0.     0.     0.     0.     0.     0.    0.    0.   0.   0.    0.    0.    0.    0.   ; # ΔFy
]')
#adjacence matrix defining the connected nodes
        #  1  2  3  4  5  6  7  8  9  10 11 12 13 14
adj =[ 0  1  0  0  0  0  0  0  0  0  0  0  0  0;  #1
       1  0  1  1  0  0  0  0  0  0  0  0  0  0;  #2
       0  1  0  0  1  0  0  0  0  0  0  0  0  0;  #3
       0  1  0  0  0  1  0  0  0  0  0  0  0  0;  #4
       0  0  1  0  0  0  1  0  0  0  0  0  0  0;  #5
       0  0  0  1  0  0  1  0  0  0  0  0  0  0;  #6
       0  0  0  0  1  1  0  1  0  0  0  0  0  0;  #7
       0  0  0  0  0  0  1  0  1  1  0  0  0  0;  #8
       0  0  0  0  0  0  0  1  0  0  1  0  0  0;  #9
       0  0  0  0  0  0  0  1  0  0  0  1  0  0;  #10
       0  0  0  0  0  0  0  0  1  0  0  0  1  0;  #11
       0  0  0  0  0  0  0  0  0  1  0  0  1  0;  #12
       0  0  0  0  0  0  0  0  0  0  1  1  0  1;  #13
       0  0  0  0  0  0  0  0  0  0  0  0  1  0;] #14
adjm = adj[1:10,1:10]
#defining which nodes can move/are clamped
tps = Dict{Int,Boundary}(1=>Clamp(),2=>Branch(),3=>Branch(),4=>Branch(),5=>Branch(),6=>Branch(),7=>Branch(),8=>Branch(),9=>Branch(),10=>Branch(),11=>Branch(),12=>Branch(),13=>Branch(),14=>Clamp())
tpsm = Dict{Int,Boundary}(1=>Clamp(),2=>Branch(),3=>Branch(),4=>Branch(),5=>Branch(),6=>Branch(),7=>Branch(),8=>Branch(),9=>Branch(),10=>Branch())
# 
con = Connections(adj,tps)
adj_beam = edge_adjacence(con)
str = Structure(con)

#solving the system with Trust-region mehtod schooting the moments/forces for each beam 
nsol = nlsolve((r,x) -> loss!(r,str,x,collect(nodefeatures[:,1:6]'),collect(beamfeatures')),0.1rand(Float64,90))
nsol.zero
#solved state
inits = initialize(str,nsol.zero,collect(nodefeatures[:,1:6]'),collect(beamfeatures'))
# you can use the next line to plot the structure without any forces
inits = initialize(str,zeros(Float64,90),collect(nodefeatures[:,1:6]'),collect(beamfeatures')) 
# inits ./100e-3^2

x = zeros(Float64,90)
r = similar(x)
loss!(r,str,x,collect(nodefeatures[:,1:6]'),collect(beamfeatures'))
po = findall(x->x>0.0005,r)

odesol = str(inits,true)
#plot solution 
plt = plot();
for ind in eachindex(odesol)
    plot!(odesol[ind],idxs = ((x,y)->beamfeatures[ind,1].*(x,y),3,4),label = ind)
end 
plot!(aspect_ratio = :equal, lims =:auto)

#Node 2
odesol[1](1)[[3,4]]*(L1)*1000
odesol[1](1)[[3,4]]*(L1)*1000 - [Nx2;Ny2]*1000


#Node 3
odesol[2](1)[[3,4]]*(R1*alpha/2)
odesol[2](1)[[3,4]]*(R1*alpha/2)*1000 - [Nx3;Ny3]*1000
odesol[2](1)[2] -(π/2 + alpha/2)
odesol[4](0)[2] -(π/2 + alpha/2)
#Node 4
odesol[3](1)[[3,4]]*(R1*alpha/2)*1000
odesol[3](1)[[3,4]]*(R1*alpha/2)*1000 - [Nx4;Ny4]*1000

#Node 5
odesol[4](1)[[3,4]]*(π*R1)*1000
odesol[4](1)[[3,4]]*(π*R1)*1000 - [N5x;N5y]*1000
odesol[4](1)[2]
odesol[6](0)[2]
odesol[6](1)[2]
#Node 6
odesol[5](1)[[3,4]]*π*R1*1000
odesol[5](1)[[3,4]]*π*R1*1000 - [N6x;N6y]*1000

#moments 
odesol[1](0)[[1]]*(E*Iz)*-1000/(π*R1/2)

π/2 - alpha/2 

savefig("example")