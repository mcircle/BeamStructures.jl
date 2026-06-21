using GeometryBasics,Dates
# using MeshIO,FileIO

# function maxheight(beams)
#     maximum(beam -> beam.h, beams)
# end 

# function createarea(x,beam)
#     θ,x,y = x
#     nx,ny = sincos(θ)
#     x1 = x * beam.l - nx * beam.h / 2
#     y1 = y * beam.l + ny * beam.h / 2
#     x2 = x * beam.l + nx * beam.h / 2
#     y2 = y * beam.l - ny * beam.h / 2
#     p1 = Point3f(x1,y1,-beam.w ÷ 2)
#     p2 = Point3f(x2,y2,-beam.w ÷ 2)
#     p3 = Point3f(x1,y1,beam.w ÷ 2)
#     p4 = Point3f(x2,y2,beam.w ÷ 2)
#     [p1,p2,p3,p4]
# end 


# function createpoints(sol::ODESolution{T,N,Q},s::AbstractVector{<:Real},height,width,len) where{T,N,Q}
#     xarr = Matrix{T}(undef,length(s),5)
#     yarr = Matrix{T}(undef,length(s),5)
#     zarr = Matrix{T}(undef,length(s),5) 

#     for idx in eachindex(s)
#         θ,x,y = sol(s[idx],idxs =2:4)
#         nx,ny = sincos(θ + π/2f0)
        
#         x *= len
#         y *= len
#         ny *= height/2
#         nx *= height/2
#         xarr[idx,[1,4,5]] .= x + ny 
#         xarr[idx,[2,3]]   .= x - ny
#         yarr[idx,[1,4,5]] .= y + nx
#         yarr[idx,[2,3]]   .= y - nx
#         zarr[idx,[1,2,5]] .= -width/2
#         zarr[idx,[3,4]]   .= width/2
#     end
#     (xarr,yarr,zarr)
# end 
# createpoints(sol,s,beam) = createpoints(sol,s,beam.h,beam.w,beam.l)
# function createpoints(sol::ODESolution{T,N,Q},s::Real,beam) where{T,N,Q}
#     createpoints(sol,[s],beam.h,beam.w,beam.l)
# end

# function createarearing(sol::ODESolution{T,N,Q},s::AbstractVector{<:Real},beam) where{T,N,Q}
#     # xarr = Matrix{T}(undef,length(s),5)
#     # yarr = Matrix{T}(undef,length(s),5)
#     # zarr = Matrix{T}(undef,length(s),5) 
#     heigth = beam.h
#     width = beam.w
#     len = beam.l
#     pts = Vector{Meshes.Point}(undef,2 * length(s))
#     for idx in eachindex(s)
#         θ,x,y = sol(s[idx],idxs =2:4)
#         nx,ny = sincos(θ + π/2f0)
        
#         x *= len
#         y *= len
#         ny *= heigth/2
#         nx *= heigth/2
#         pts[idx] = Meshes.Point(x + ny,y + nx)
#         pts[end-idx+1] = Meshes.Point(x - ny,y - nx)

#     end
#     pts
# end

# function createmesh(sol::ODESolution{T,N,Q},beam) where{T,N,Q}

#     wxy = sol.(0:0.001:1,idxs = 2:4)
#     poly = mapreduce(x->createarea(x,beam),vcat,wxy) 
#     sz = length(poly)
#     faces = Vector{TriangleFace{Int}}()
#     sizehint!(faces, 2 * length(wxy)) # 2 faces per 4 points
#     #stirnseite 1 
#     push!(faces,TriangleFace(1,3,2))
#     push!(faces,TriangleFace(2,3,4))
#     #stirnseite 2
#     push!(faces,TriangleFace(sz-3,sz-2,sz-1))
#     push!(faces,TriangleFace(sz-2,sz,sz-1))
#     #mantelfläche
#     for i in 1:4:(sz-4)

#         push!(faces,TriangleFace(i,i+1,i+5))
#         push!(faces,TriangleFace(i,i+5,i+4))
        
#         push!(faces,TriangleFace(i+1,i+3,i+7))
#         push!(faces,TriangleFace(i+1,i+7,i+5))
        
#         push!(faces,TriangleFace(i+2,i+6,i+7))
#         push!(faces,TriangleFace(i+2,i+7,i+3))
        
#         push!(faces,TriangleFace(i,i+6,i+2))
#         push!(faces,TriangleFace(i,i+4,i+6))
#     end
#     me = GeometryBasics.Mesh(poly,faces)
# end 

# function createmesh(str,beams,nodes,inits::AbstractMatrix; name = str_$(now()) )
#     sol = str(inits,beams,nodes,true) # get EnsembleSolution for all beams 
#     isfile("./meshobjects") && mkdir("./meshobjects")
#     meobjs = Vector{GeometryBasics.Mesh}(undef, length(beams))
#     for idx in eachindex(sol)
#         meobjs[idx] = createmesh(sol[idx],beams[idx])    
#     end
#     merged = GeometryBasics.merge(meobjs)
#     save("./meshobjects/$name.stl",merged)
#     merged,meobjs
# end 



# function createpoints(sol::ODESolution{T,N,Q},s::AbstractVector{<:Real},beam) where{T,N,Q} 
#     ptsfor = Vector{Point2f}(undef,length(s))
#     ptsback = Vector{Point2f}(undef,length(s))
#     len = beam.l
#     heigth = beam.h
#     # width = beam.w

#     for idx in eachindex(s)
#         θ,x,y = sol(s[idx],idxs =2:4)
#         nx,ny = sincos(θ + π/2f0)
#         x *= len
#         y *= len
#         ny *= heigth/2
#         nx *= heigth/2
#         ptsfor[idx] = Point2f(x + ny,y + nx)
#         ptsback[end-idx+1] = Point2f(x - ny,y - nx)
#     end
#     (ptsfor,ptsback)
# end
    

# function signed_dist_to_plane(n,p::Point3f, plane_point::Point3f)
#     return dot(n, p - plane_point)
# end

# function point_in_triangle_3d(p::Point3f, tri::Vector{Point3f}; eps=eps(Float32))
#     # check angles of point. schould be between angles  
#     a, b, c = tri
#     for idx in 1:3
#         v1 = tri[idx]
#         v2 = tri[mod1(idx+1,3)]
#         v3 = tri[mod1(idx+2,3)]
#         rtv1 = v2 - v1 # richtungsvektor
#         rtv2 = v3 - v1
#         rtv3 = p - v1
#         #angle between rtv1 and rtv2 should be bigger than angle between rtv1 and rtv3 then point in triangle
#         nrmvrt1 = norm(rtv1)
#         acos(dot(rtv1,rtv2) / (norm(rtv2) * nrmvrt1)) ≥ acos(dot(rtv1,rtv3) / (norm(rtv3) * nrmvrt1)) || return false
#     end
#     true 
# end

# function faces_intersecting(face1::Vector{Point3f},face2::Vector{Point3f},n1,n2)
#     # now for Trianglefaces 
#     @assert length(face1) == 3 && length(face2) == 3 "Only triangle faces are supported"

#     a1,b1,c1 = face1
#     a2,b2,c2 = face2
    
#     isapprox(cross(n1,n2), Point3f(0,0,0)) && return false #planes are parallel or coincident
    
#     # check if linesegemnents of face2 intersect with plane of face1 or vice versa by one point must be on/other sider of plane 
#     distpoints_f2 = map(x->signed_dist_to_plane(n1,x,a1), face2)
#     if all(distpoints_f2 .> 0) || all(distpoints_f2 .< 0)
#         return false # all points of face2 are on the same side of face1
#     end
#     distpoints_f1 = map(x->signed_dist_to_plane(n2,x,a2), face1)
#     if all(distpoints_f1 .> 0) || all(distpoints_f1 .< 0)
#         return false # all points of face1 are on the same side of face2
#     end

#     return true

# end


# function faces_intersecting!(ori::Vector{T},dir::Vector{T}, face1::Vector{Point{N,T}},face2::Vector{Point{N,T}}) where{T<:Real,N}
#     @assert length(ori) == N && length(dir) == N "Origin and direction vectors must have the same dimension as the points"
#     a1,b1,c1 = face1
#     a2,b2,c2 = face2
#     n1 = cross(a1 - b1,a1 - c1)
#     n2 = cross(a2 - b2,a2 - c2)
#     faces_intersecting(face1,face2,n1,n2) || return false
#     p1 = dot(n1,a1)
#     p2 = dot(n2,a2)

#     det = dot(n1,n1) * dot(n2,n2) - dot(n1, n2)^2
#     d12 = dot(n1, n2)

#     c1 = p1 * dot(n2, n2) + p2 * d12
#     c1 /= det
#     c2 = p2 * dot(n1, n1) + p1 * d12
#     c2 /= det
    
#     ori .= (c1 * n1 + c2 * n2) 
#     dir .= cross(n2,n1)
    
#     return true
# end 

# function edge_intersecting(face,ori,dir,idx)
#     v1 = face[idx]
#     v2 = face[mod1(idx+1,3)]
#     v3 = face[mod1(idx+2,3)]
#     nrm = cross(v2 - v1, v3 - v1)
#     edge_dir = v2 - v1
#     #create plane trough edge and normal of face 
#     edge_nrm = cross(nrm, edge_dir)
#     # edge_nrm dot (P2 -P1) P1 and P2 are on line
#     dem = dot(edge_nrm,dir)
#     dem ≈ 0 && return nothing # line is parallel to plane

#     # edge_nrm dot (P3 -P1) P3 is on plane 
#     nom = dot(edge_nrm,v2 - ori)
#     #  edge_nrm dot (P3 -P1) / edge_nrm dot (P2 -P1) gives u for line equation P = P1 + u * (P2 -P1)
#     u = nom / dem
#     intersect_pt = Point3f(ori + u * dir)
    
#     return intersect_pt
# end

# function points_intersecting_faces(face1::Vector{Point3f},face2::Vector{Point3f})
#     #ori and dir of line of intersection between planes of face1 and face2
#     ori = Vector{Float32}(undef, 3)
#     dir = Vector{Float32}(undef, 3)
#     a1,b1,c1 = face1
#     a2,b2,c2 = face2


#     intersecs = faces_intersecting!(ori,dir,face1,face2) 
#     #find points of intersection between line and facesedges
#     # create new plane with normal of edge and normal of face and find intersection point with ori and dir
#     intersecting_pts = Vector{Point3f}()
#     for f in (face1,face2)
#         for idx in 1:3
#             intersect_pt = edge_intersecting(f,ori,dir,idx)
#             isnothing(intersect_pt) && continue     
#             push!(intersecting_pts, intersect_pt)
#         end
#     end 
#     return intersecting_pts
# end



using BSplineKit, StaticArrays, LinearAlgebra, Gmsh
using DifferentialEquations

# ─────────────────────────────────────────────────────────
# 1. Glatte neutrale Faser aus ODESolution
# ─────────────────────────────────────────────────────────

function fit_centerline(sol::ODESolution,beam)
    # idxs 3:4 → x, y Koordinaten der ODESolution
    # Auswertung Kurve und Ableitung
    function centerline(s) 
        x,y = sol(s, idxs=3:4) .* beam.l
        SVector(x, y, 0.0)
    end 

    function dcenterline(s)
        dx,dy = sol(s,Val{1}; idxs=3:4) .* beam.l
        SVector(dx, dy, 0.0)
    end 
    return centerline, dcenterline
end

# ─────────────────────────────────────────────────────────
# 2. Lokales Frenet-Frame entlang der Faser
#    Für Balken in xy-Ebene: Binormale = z-Achse (stabil)
# ─────────────────────────────────────────────────────────

function beam_frame(dcenterline, s)
    tang  = normalize(dcenterline(s))
    binom = SVector(0.0, 0.0, 1.0)        # stabil für planare Balken
    norm_ = normalize(cross(binom, tang))  # liegt in xy-Ebene
    return tang, norm_, binom
end

# ─────────────────────────────────────────────────────────
# 3. Querschnittsecken im Weltkoordinatensystem
# ─────────────────────────────────────────────────────────

function crosssection_corners(centerline, dcenterline, s, hw, hh)
    c       = centerline(s)
    _, n, b = beam_frame(dcenterline, s)
    # Uhrzeigersinn von vorne betrachtet → wichtig für Wire-Orientierung
    return [
        c + n * hh + b * hw,   # oben rechts
        c - n * hh + b * hw,   # unten rechts
        c - n * hh - b * hw,   # unten links
        c + n * hh - b * hw,   # oben links
    ]
end

# ─────────────────────────────────────────────────────────
# 4. Einen Wire (Querschnitt) in Gmsh anlegen
#    Gibt den Wire-Tag zurück
# ─────────────────────────────────────────────────────────

function add_wire!(corners, meshsize = 0.0)
    npts = length(corners)

    # Punkte anlegen
    pt_tags = map(corners) do c
        gmsh.model.occ.addPoint(c[1], c[2], c[3], meshsize)
    end

    # Linien zwischen den Punkten (geschlossen)
    ln_tags = map(1:npts) do i
        p1 = pt_tags[i]
        p2 = pt_tags[mod1(i + 1, npts)]
        gmsh.model.occ.addLine(p1, p2)
    end

    # CurveLoop und Wire
    # gmsh.model.occ.addCurveLoop(ln_tags,loop_tag[])
    wt = gmsh.model.occ.addWire(ln_tags)

    # wt = wire_tag[]
    # loop_tag[]  += 1
    # wire_tag[]  += 1
    return wt
end

# ─────────────────────────────────────────────────────────
# 5. Einen einzelnen Balken als ThruSection-Solid anlegen
#    Gibt den Solid-Tag zurück
# ─────────────────────────────────────────────────────────

function add_beam_solid!(centerline, dcenterline, beam,
                         ;
                         nstations=20,meshsize=0.0)

    svals    = range(0.0, 1.0, length=nstations)
    hw, hh   = beam.w / 2, beam.h / 2
    wire_tags = Int[]

    for s in svals
        corners = crosssection_corners(centerline, dcenterline, s, hw, hh)
        wt = add_wire!(corners, meshsize)
        push!(wire_tags, wt)
    end

    # Solid durch alle Querschnitte legen
    st = gmsh.model.occ.addThruSections(wire_tags, -1, true, false)
    return only(st)
end

# ─────────────────────────────────────────────────────────
# 6. Alle Balken → fuse → STEP
# ─────────────────────────────────────────────────────────

function export_structure_step(str, beams, nodes;
                                filename="structure",
                                nstations=20)
    inits = zeros(Float32,beams,nodes)
    sol = str(inits, beams, nodes, true)  # EnsembleSolution

    gmsh.isInitialized() == 1 || gmsh.initialize()
    gmsh.model.add(filename)

    solid_tags = Tuple{Int32,Int32}[]

    for idx in eachindex(sol)
        cl, dcl = fit_centerline(sol[idx], beams[idx])
        st = add_beam_solid!(cl, dcl, beams[idx]; nstations=nstations,meshsize = 0.0)
        push!(solid_tags, st)
    end

    # if length(solid_tags) > 1
    #     gmsh.model.occ.fuse(solid_tags[1:1], solid_tags[2:end])
    # end

    gmsh.model.occ.synchronize()
    gmsh.write("$(filename).step")
    gmsh.finalize()
    # @info "Exportiert: $(filename).step"
end

