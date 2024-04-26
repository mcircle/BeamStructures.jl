
include("BeamStructures.jl")
using Plots
plotlyjs(size = (800,600))
#               beam    1   2   3   4
beamfeatures = collect([100e-3 100e-3  100e-3 100e-3 100e-3 40e-3; #length of the beam 
                        0.     0.     π/2    π/2    0.    0π/4   ; #θs angle as going out slope 
                        0π     0.     0.     0.     0.    2π ; #κ as curvatature of the beam 
                        0.     0.     π/2    π/2    0.    0 ; #θe angle as going in slope
]')
#nodes as intersections/Connections of the beams
#changing 
E = 2e11
h = 2.e-3
w = 10.e-3
Iz = w*h^3/12

beam1 = Beam(100e-3,h,w,0.)
beam2 = Beam(100e-3,h,w,0.)
beam3 = Beam(100e-3,h,w,0.,θs = π/2)
beam4 = Beam(100e-3,h,w,0.,θs = π/2)
beam5 = Beam(100e-3,h,w,0.)
beam6 = Beam(40e-3,h,w,0., θs = 0.,)

ΔFy = 100. 
ΔFx = 100. 
    #            node   1    2       3       4       5
nodefeatures = collect([0.   0.      0.      0.     0.     π/4; # Drehung des Nodes um Δθ
                        0.   100e-3  200e-3  100e-3  200e-3 240e-3; #position x 
                        0.   0.      0.      100e-3  100e-3 120e-3;# position y 
                        0.   0.      0.      0.      0.     0.; # ΔM
                        0.   0ΔFx     0ΔFx    -0ΔFx   -0ΔFx    0.; # ΔFx
                        0.   0ΔFy     -0ΔFy   0ΔFy    -0ΔFy    0.; # ΔFy
]')
#adjacence matrix defining the connected nodes
# node 1 2 3 4 5  node
no1 = CompliantClamp(zeros(Float64,5)...,100.)
@time fieldnames(typeof(no1))
no2 = Branch(100e-3,0.,0.,0.,0.,0.)
no3 = Branch(200e-3,0.,0.,100.,-100.,0.)
no4 = Branch(100e-3,100e-3,0.,-100.,100.,0.)
no5 = Branch(200e-3,100e-3,0.,0.,0.,10.)
no6 = Free(240e-3,100e-3,0.,0.,0.,0.)

test = loaddatabase() 

#closed chain 
adj =[ 0 1 0 0 0 0; #1
       1 0 1 1 0 0; #2
       0 1 0 0 1 0; #3
       0 1 0 0 1 0;
       0 0 1 1 0 1;
       0 0 0 0 1 0]

con = Connections(adj,no1,no2,no3,no4,no5,no6)
con.Nodes[1]

str = Structure(con,beam1,beam2,beam3,beam4,beam5,beam6)


str.AdjMat.Nodes[1].x
@time str = change_node(str,2;x=1,y = 2,ϕ = 3)
str.AdjMat.Nodes[1].ϕ
op = :x
str2 = @set str.AdjMat.Nodes[3].$op = 1.1

r = rand(Float64,36)
s = zeros(Float64,36)
inits = initialize(str,s)
str(r,s)
r[36]
findall(x->x>=0.00001,r)
#solving the system with Trust-region mehtod schooting the moments/forces for each beam 
nsol = nlsolve((r,x) -> str(r,x),zeros(Float64,36))
nsol.zero
#solved state

odesol = str(nsol.zero,true)
size(odesol)
#plot solution 

plt = plot();
for ind in eachindex(odesol)
    # plot!(odesol[ind],idxs = ((x,y)->beamfeatures[ind,1].*(x,y),3,4),label = ind)
    plot!(odesol[ind],idxs = ((x,y)->str.Beams[ind].l.*(x,y),3,4),label = ind)
end 
plot!(aspect_ratio = :equal, lims =:auto)

x = zeros(Float64,30)
r = similar(x)
loss!(r,str,x,collect(nodefeatures[:,1:6]'),collect(beamfeatures'))
#open chain
adj =[ 0 1 0 0 ; #1
       1 0 1 1 ; #2
       0 1 0 0 ; #3
       0 1 0 0 ;]      
tps = Dict{Int,Boundary}(1 => Clamp(), 2 => Branch(),3=>Branch(),4=>Branch())
con = Connections(adj,tps)
adj_beam = edge_adjacence(con)
str = Structure(con)

#solving the system with Trust-region mehtod schooting the moments/forces for each beam 
nsol = nlsolve((r,x) -> loss!(r,str,x,collect(nodefeatures[:,1:6]'),collect(beamfeatures')),zeros(Float64,30))
nsol.zero
#solved state
inits = initialize(str,nsol.zero,collect(nodefeatures[:,1:6]'),collect(beamfeatures'))
odesol = str(inits,true)
#plot solution 

plt = plot();
for ind in eachindex(odesol)
    plot!(odesol[ind],idxs = ((x,y)->beamfeatures[ind,1].*(x,y),3,4),label = ind)
end 
plot!(lims =:auto)

#for debugging and control
x = zeros(Float64,30)
r = collect(1.:30.)
loss!(r,str,nsol.zero,collect(nodefeatures[:,1:6]'),collect(beamfeatures'))
r
odesol[1](1)[[5,6]]   
odesol[2](1)[[5,6]]./ 80e-3^2    
odesol[3](1)[[5,6]]./ 100e-3^2    
odesol[4](1)[[5,6]]./ 80e-3^2 

rad2deg(odesol[2](1)[2]) 
odesol[5](1)[[3,4]] 
x1,y1 = odesol[4](0)[[3,4]] .*100
x2,y2 = odesol[5](0)[[3,4]] .*80 



