import trimesh
import numpy as np
import networkx as nx
import time
import math
#Discrete Differential-Geometry Operators for Triangulated 2-Manifolds
#https://arxiv.org/pdf/1706.02413.pdf
#https://computergraphics.stackexchange.com/questions/1718/what-is-the-simplest-way-to-compute-principal-curvature-for-a-mesh-triangle
import csv
def write_vertex_attributes(path,c):
	with open(path, mode='w') as file_to_write:
	    writer = csv.writer(file_to_write, delimiter=',')
	    for v,curv in enumerate(c):
	    	writer.writerow([v,curv])


def average_over_one_ring(mesh, c):
    g = nx.from_edgelist(mesh.edges_unique)
    one_rings = [list(g[i].keys()) for i in range(len(m.vertices))]

    avg = [0]*len(c)
    for vertex_id,one_ring in enumerate(one_rings):
        avg[vertex_id]=np.mean(c[one_ring])

    return np.array(avg)

def average_over_n_ring(mesh,c,n):
	g = nx.from_edgelist(mesh.edges_unique)
	one_rings = [list(g[i].keys()) for i in range(len(m.vertices))]

    avg = [0]*len(one_rings)
	for vertex_id in range(len(one_rings)):

        n_ring_vertex_ids = []
        current_ring = one_rings[vertex_id]
        for ring in range(n-1):
            #can speed up so not rechecking
            current_ring = np.unique( [one_ring[i] for i in current_ring] )
            n_ring_vertex_ids.append(current_ring)

        avg[vertex_id] = np.mean(c[n_ring_vertex_ids])

    return np.array(avg)




def my_discrete_mean_curvature_measure(mesh):
    #Area (should be easy)

    g = nx.from_edgelist(mesh.edges_unique)
    one_rings = [list(g[i].keys()) for i in range(len(m.vertices))]
    print(len(one_rings))
    print(len(mesh.vertices))
    #tangents = np.tan(m.face_angles_sparse.tocsr()) #vertices by faces
    start = time.time()
    face_angles = m.face_angles_sparse
    cotangents = { f"{vertex},{face}":1/np.tan(angle) for vertex,face,angle in zip(face_angles.row, face_angles.col, face_angles.data)}
    end = time.time()
    # print(end-start)
    # start = time.time()
    # fa = m.face_angles_sparse.todok()
    # cotangents = fa.update((key,1/np.tan(angle)) for key,value in fa.items())
    # end = time.time()
    # print(end-start)

    #adjacency_unshared_sorted = 
    print("b")
    #[print((mesh.vertices[0]- mesh.vertices[1])*cotangents[key]) for key in cotangents]
    fa = mesh.face_adjacency
    fae = mesh.face_adjacency_edges
    edge_measure = {f"{fae[i][0]},{fae[i][1]}":(mesh.vertices[fae[i][1]] - mesh.vertices[fae[i][0]])*(cotangents[f"{v[0]},{fa[i][0]}"]+cotangents[f"{v[1]},{fa[i][1]}"]) for i,v in enumerate(mesh.face_adjacency_unshared) }
    
    # for i,v in enumerate(mesh.face_adjacency_edges):
    #     if v[0]==1518 or v[1]==1518 or v[0]==21390 or v[1]==21390:
    #             print(f"fae {v}")


    # for i,v in enumerate(mesh.edges_unique):
    #     if v[0]==1518 or v[1]==1518 or v[0]==21390 or v[1]==21390:
    #         print(f"eu {v}")

    print("c")
    #get all faces in vertices 
    #face_adjacency_unshared Return the vertex index of the two vertices not in the shared edge between two adjacent faces
    mean_curv = [0]*len(mesh.vertices)
    percentage = int(len(mesh.vertices)/10)

    for vertex_id,face_ids in enumerate(mesh.vertex_faces):
        face_ids = face_ids[face_ids!=-1] #faces associated with vertex_id
        one_ring = one_rings[vertex_id]
        delta_s = 0;
        for one_ring_vertex_id in one_ring:
            if f"{vertex_id},{one_ring_vertex_id}" in edge_measure:
                delta_s += edge_measure[f"{vertex_id},{one_ring_vertex_id}"]
            elif f"{one_ring_vertex_id},{vertex_id}"  in edge_measure:
                delta_s -= edge_measure[f"{one_ring_vertex_id},{vertex_id}"]
            #else :
                #print("hey")

        delta_s *= 1/(2*sum(mesh.area_faces[face_ids])/3) #use 1/3 of the areas
        mean_curv[vertex_id] = 0.5*np.linalg.norm(delta_s)
        if vertex_id % percentage == 0:
        	print(vertex_id/len(mesh.vertices))
        #print(mean_curv[vertex_id])

        #delta_s=sum([edge_measure[f"{vertex_id},{one_ring_vertex_id}"] if f"{vertex_id},{one_ring_vertex_id}" in edge_measure else if f"{one_ring_vertex_id},{vertex_id}" in edge_measure -1*edge_measure[f"{one_ring_vertex_id},{vertex_id}"] for one_ring_vertex_id in one_ring ])

    return np.array(mean_curv)


def my_discrete_gaussian_curvature_measure(mesh):
    """
    Return the discrete gaussian curvature measure of a sphere centered
    at a point as detailed in 'Restricted Delaunay triangulations and normal
    cycle', Cohen-Steiner and Morvan.
    Parameters
    ----------
    points : (n,3) float, list of points in space
    radius : float, the sphere radius
    Returns
    --------
    gaussian_curvature: (n,) float, discrete gaussian curvature measure.
    """

    #points = np.asanyarray(points, dtype=np.float64)
    #if not util.is_shape(points, (-1, 3)):
    #    raise ValueError('points must be (n,3)!')

    g = nx.from_edgelist(mesh.edges_unique)
    #nearest = mesh.kdtree.query_ball_point(points, radius)
    one_ring = [list(g[i].keys()) for i in range(len(m.vertices))]
    #print(one_ring)

    #gauss_curv = [mesh.vertex_defects[vertices].sum() for vertices in one_ring]
    print(type(mesh.vertex_defects))
    #[print(mesh.vertex_faces[vertex]) for vertex in range(len(mesh.vertices))]
    
    #FACTOR OF 3 since using 1/3 area?
    gauss_curv = [3*mesh.vertex_defects[vertex]/mesh.area_faces[ mesh.vertex_faces[vertex][mesh.vertex_faces[vertex]!=-1] ].sum() for vertex in range(len(mesh.vertices))]

    #one_ring = [list(g[i].keys()) for i in range(len(m.vertices))]
    #gauss_curv = [gauss_curv[vertices] for vertices in one_ring]
    return np.asarray(gauss_curv)





# meshes = [trimesh.creation.uv_sphere() for i in range(10)]
# for i, m in enumerate(meshes):
# 	m.vertices*=(np.random.random(3) +1)*2
# 	m.apply_translation([0,0, i*6])
# 	#radii = np.linalg.norm(m.vertices - m.center_mass, axis=1)
# 	#m.visual.vertex_colors = trimesh.visual.interpolate(radii, color_map='viridis')
# 	c = trimesh.curvature.discrete_gaussian_curvature_measure(m, m.vertices, 2.0)/(4*np.pi)
# 	m.visual.vertex_colors = trimesh.visual.interpolate(c, color_map='viridis')
radius = 5
m= trimesh.creation.icosphere(radius=radius)

print(my_discrete_gaussian_curvature_measure(m)/(1/(radius*radius)))
print(my_discrete_mean_curvature_measure(m)/(0.5*(1/radius + 1/radius)))
c=my_discrete_mean_curvature_measure(m)

m = trimesh.load('/groups/cosem/cosem/ackermand/paperResultsWithFullPaths/collected/jrc_ctl-id8_a01.n5/mesh/1.obj')#'multiresolutionMeshes/test/mito_obj_meshes_s2/345809856042.obj')'multiresolutionMeshes/test/mito_obj_meshes_s2/345809856042.obj')#'/groups/cosem/cosem/ackermand/paperResultsWithFullPaths/collected/jrc_ctl-id8_a01.n5/mesh/1.obj')#'multiresolutionMeshes/test/mito_obj_meshes_s2/345809856042.obj')
m.remove_duplicate_faces()
m.remove_unreferenced_vertices()
#trimesh.smoothing.filter_humphrey(m)
print(m.is_watertight)
#c = np.log(np.abs(trimesh.curvature.discrete_gaussian_curvature_measure(m, m.vertices, 500.0)))
c = my_discrete_gaussian_curvature_measure(m)
write_vertex_attributes('/groups/cosem/cosem/ackermand/meshesForAubrey/gaussianCurvature.csv',c)

c=np.clip(c,np.percentile(c,0), np.percentile(c,90))
avg = average_over_one_ring(m,c)
avg = average_over_one_ring(m,c)


#c = np.log(np.abs(c))
m.visual.vertex_colors = trimesh.visual.interpolate(c, color_map='viridis')
m.export("/groups/cosem/cosem/ackermand/meshesForAubrey/1_gc_low_res.obj")

#avg = np.log(np.abs(avg))
m.visual.vertex_colors = trimesh.visual.interpolate(avg, color_map='viridis')
m.export("/groups/cosem/cosem/ackermand/meshesForAubrey/1_gc_low_res_avg.obj")

# c = np.log(np.abs(my_discrete_gaussian_curvature_measure(m)))
# m.visual.vertex_colors = trimesh.visual.interpolate(c, color_map='viridis')
# m.export("my_gauss.obj")

c = my_discrete_mean_curvature_measure(m)
write_vertex_attributes('/groups/cosem/cosem/ackermand/meshesForAubrey/meanCurvature.csv',c)


avg = average_over_one_ring(m,c)
avg = average_over_one_ring(m,avg)
avg = average_over_one_ring(m,avg)
#avg = average_over_one_ring(m,avg)
#avg = average_over_one_ring(m,avg)



c=np.clip(c,np.percentile(c,0), np.percentile(c,90))
c+=0.000001
#c = np.log(c)
m.visual.vertex_colors = trimesh.visual.interpolate(c, color_map='viridis')
m.export("/groups/cosem/cosem/ackermand/meshesForAubrey/1_mc_low_res.obj")


avg+=0.000001
#avg=10**((avg-np.min(avg))/np.max(avg)-np.min(avg))
avg=np.clip(avg,np.percentile(avg,30), np.percentile(avg,90))
#avg = np.log(avg)
m.visual.vertex_colors = trimesh.visual.interpolate(avg, color_map='viridis')
m.export("/groups/cosem/cosem/ackermand/meshesForAubrey/1_mc_low_res_avg3.obj")



#c=trimesh.curvature.discrete_mean_curvature_measure(m, m.vertices, 500.0)/(4*np.pi)
#m.visual.vertex_colors = trimesh.visual.interpolate(c, color_map='viridis')
#m.export("mean.obj")

# meshes.append(m)
# trimesh.Scene(meshes).show()

# c=trimesh.curvature.discrete_mean_curvature_measure(m, m.vertices, 500.0)
# c=np.clip(c,np.percentile(c,0), np.percentile(c,90))
# m.visual.vertex_colors = trimesh.visual.interpolate(c, color_map='viridis')
# m.export("/groups/cosem/cosem/ackermand/meshesForAubrey/1_mc_low_res_theirs.obj")



