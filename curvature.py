import trimesh
import numpy as np
import networkx as nx
import time
import math
import copy
#Discrete Differential-Geometry Operators for Triangulated 2-Manifolds
#https://arxiv.org/pdf/1706.02413.pdf
#https://computergraphics.stackexchange.com/questions/1718/what-is-the-simplest-way-to-compute-principal-curvature-for-a-mesh-triangle
import csv
def write_vertex_attributes(path,c):
	with open(path, mode='w') as file_to_write:
	    writer = csv.writer(file_to_write, delimiter=',')
	    for v,curv in enumerate(c):
	    	writer.writerow([v,curv])

def write_mesh(m,c,path):
    m.visual.vertex_colors = trimesh.visual.interpolate(c, color_map='viridis')
    m.export(path)

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

        current_ring = one_rings[vertex_id]
        n_ring_vertex_ids = current_ring
        for ring in range(n-1):
            #can speed up so not rechecking
            next_ring = []
            for i in current_ring:
                next_ring.extend(one_rings[i])
            current_ring = np.unique(next_ring)
            n_ring_vertex_ids.extend(current_ring)
        avg[vertex_id] = np.mean(c[np.unique(n_ring_vertex_ids)])

    return np.array(avg)


def my_discrete_mean_curvature_measure_cleaned(mesh):
    g = nx.from_edgelist(mesh.edges_unique)
    one_rings = [list(g[i].keys()) for i in range(len(m.vertices))]
   
    face_angles = m.face_angles_sparse
    cotangents = { f"{vertex},{face}":1/np.tan(angle) for vertex,face,angle in zip(face_angles.row, face_angles.col, face_angles.data)}
    end = time.time()

    fa = mesh.face_adjacency
    fae = mesh.face_adjacency_edges
    edge_measure = {f"{fae[i][0]},{fae[i][1]}":(mesh.vertices[fae[i][1]] - mesh.vertices[fae[i][0]])*(cotangents[f"{v[0]},{fa[i][0]}"]+cotangents[f"{v[1]},{fa[i][1]}"]) for i,v in enumerate(mesh.face_adjacency_unshared) }
  
    mean_curv = [0]*len(mesh.vertices)

    for vertex_id,face_ids in enumerate(mesh.vertex_faces):
        face_ids = face_ids[face_ids!=-1] #faces associated with vertex_id
        one_ring = one_rings[vertex_id]
        delta_s = 0
        for one_ring_vertex_id in one_ring:
            if f"{vertex_id},{one_ring_vertex_id}" in edge_measure:
                delta_s += edge_measure[f"{vertex_id},{one_ring_vertex_id}"]
            elif f"{one_ring_vertex_id},{vertex_id}"  in edge_measure:
                delta_s -= edge_measure[f"{one_ring_vertex_id},{vertex_id}"]
 
        delta_s *= 1/(2*sum(mesh.area_faces[face_ids])/3) #use 1/3 of the areas
        mean_curv[vertex_id] = 0.5*np.linalg.norm(delta_s)
       
    return np.array(mean_curv)

def my_discrete_mean_curvature_measure(mesh):
    """Calculate discrete mean curvature of mesh using one-ring neighborhood."""

    # one-rings (immediate neighbors of) each vertex
    g = nx.from_edgelist(mesh.edges_unique)
    one_rings = [list(g[i].keys()) for i in range(len(mesh.vertices))]
   
    # cotangents of angles and store in dictionary based on corresponding vertex and face
    face_angles = mesh.face_angles_sparse
    cotangents = { f"{vertex},{face}":1/np.tan(angle) for vertex,face,angle in zip(face_angles.row, face_angles.col, face_angles.data)}

    # discrete Laplace-Beltrami contribution of the shared edge of adjacent faces:
    #        /*\
    #       / * \
    #      /  *  \
    #    vi___*___vj
    #
    # store results in dictionary with vertex ids as keys     
    fa = mesh.face_adjacency
    fae = mesh.face_adjacency_edges
    edge_measure = {f"{fae[i][0]},{fae[i][1]}":(mesh.vertices[fae[i][1]] - mesh.vertices[fae[i][0]])*(cotangents[f"{v[0]},{fa[i][0]}"]+cotangents[f"{v[1]},{fa[i][1]}"]) for i,v in enumerate(mesh.face_adjacency_unshared) }
  
    # calculate mean curvature using one-ring
    mean_curv = [0]*len(mesh.vertices)
    for vertex_id,face_ids in enumerate(mesh.vertex_faces):
        face_ids = face_ids[face_ids!=-1] #faces associated with vertex_id
        one_ring = one_rings[vertex_id]
        delta_s = 0

        for one_ring_vertex_id in one_ring:
            if f"{vertex_id},{one_ring_vertex_id}" in edge_measure:
                delta_s += edge_measure[f"{vertex_id},{one_ring_vertex_id}"]
            elif f"{one_ring_vertex_id},{vertex_id}"  in edge_measure:
                delta_s -= edge_measure[f"{one_ring_vertex_id},{vertex_id}"]
        
        delta_s *= 1/(2*sum(mesh.area_faces[face_ids])/3) #use 1/3 of the areas
        mean_curv[vertex_id] = 0.5*np.linalg.norm(delta_s)
       
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

    g = nx.from_edgelist(mesh.edges_unique)
    #nearest = mesh.kdtree.query_ball_point(points, radius)
    one_ring = [list(g[i].keys()) for i in range(len(mesh.vertices))]
    #print(one_ring)

    #gauss_curv = [mesh.vertex_defects[vertices].sum() for vertices in one_ring]
    print(type(mesh.vertex_defects))
    #[print(mesh.vertex_faces[vertex]) for vertex in range(len(mesh.vertices))]
    
    #FACTOR OF 3 since using 1/3 area?
    gauss_curv = [3*mesh.vertex_defects[vertex]/mesh.area_faces[ mesh.vertex_faces[vertex][mesh.vertex_faces[vertex]!=-1] ].sum() for vertex in range(len(mesh.vertices))]

    #one_ring = [list(g[i].keys()) for i in range(len(m.vertices))]
    #gauss_curv = [gauss_curv[vertices] for vertices in one_ring]
    return np.asarray(gauss_curv)

def split_mesh_in_two(mesh, vertex_curvature, threshold):
    above_threshold_mesh = copy.deepcopy(mesh)
    below_threshold_mesh = copy.deepcopy(mesh)
    avg_face_curvature = [0]*len(mesh.faces)
    for face_id,vertex_ids in enumerate(mesh.faces):
        avg_face_curvature[face_id] = np.mean(vertex_curvature[vertex_ids])

    above_threshold_mesh.update_faces(avg_face_curvature>threshold)
    below_threshold_mesh.update_faces(avg_face_curvature<=threshold)
    return [above_threshold_mesh, below_threshold_mesh]


radius = 5
m= trimesh.creation.icosphere(radius=radius)

print(my_discrete_gaussian_curvature_measure(m)/(1/(radius*radius)))
print(my_discrete_mean_curvature_measure(m)/(0.5*(1/radius + 1/radius)))
c=my_discrete_mean_curvature_measure(m)
split_mesh_in_two(m, c, np.mean(c))
average_over_n_ring(m,c,3)
print("done with test")

ids_dictionary = {'410_roi1':[1,2],
'493_roi2':[1,3],
'494_roi2':[1,2],
'494_roi5':[2,1]
}
threshold = 0
for roi,ids in ids_dictionary.items():
    for object_id in ids:

        m = trimesh.load(f'/groups/cosem/cosem/ackermand/paperResultsWithFullPaths/collected/jrc_ctl-id8_a01.n5/{roi}/meshRescaleLevel2Smoothing2/{object_id}.obj')#'multiresolutionMeshes/test/mito_obj_meshes_s2/345809856042.obj')'multiresolutionMeshes/test/mito_obj_meshes_s2/345809856042.obj')#'/groups/cosem/cosem/ackermand/paperResultsWithFullPaths/collected/jrc_ctl-id8_a01.n5/mesh/1.obj')#'multiresolutionMeshes/test/mito_obj_meshes_s2/345809856042.obj')
        m.remove_duplicate_faces()
        m.remove_unreferenced_vertices()
        print(m.is_watertight)

        c = my_discrete_gaussian_curvature_measure(m)
        write_vertex_attributes(f'/groups/cosem/cosem/ackermand/meshesForAubrey/{roi}/lowres/{object_id}_gaussianCurvature.csv',c)
        
        c=np.clip(c,np.percentile(c,0), np.percentile(c,90))
        write_mesh(m,c,f"/groups/cosem/cosem/ackermand/meshesForAubrey/{roi}/lowres/{object_id}_gaussianCurvature.ply")
       
        avg = average_over_one_ring(m,c)
        write_mesh(m,avg,f"/groups/cosem/cosem/ackermand/meshesForAubrey/{roi}/lowres/{object_id}_gaussianCurvature_avg.ply")
       
        c = my_discrete_mean_curvature_measure(m)
        write_vertex_attributes(f'/groups/cosem/cosem/ackermand/meshesForAubrey/{roi}/lowres/{object_id}_meanCurvature.csv',c)

        c=np.clip(c,np.percentile(c,0), np.percentile(c,90))
        c+=0.000001
        write_mesh(m,c,f"/groups/cosem/cosem/ackermand/meshesForAubrey/{roi}/lowres/{object_id}_meanCurvature.ply")

        avg = average_over_one_ring(m,c)
        avg = average_over_one_ring(m,avg)
        avg = average_over_one_ring(m,avg)
        avg+=0.000001
        avg=np.clip(avg,np.percentile(avg,30), np.percentile(avg,90))

        write_mesh(m,avg,f"/groups/cosem/cosem/ackermand/meshesForAubrey/{roi}/lowres/{object_id}_meanCurvature_avg3.ply")
        threshold = 0.5*(np.amax(avg)+np.amin(avg))
        m_split = split_mesh_in_two(m,avg,threshold)
        m_split[0].export(f"/groups/cosem/cosem/ackermand/meshesForAubrey/{roi}/lowres/{object_id}_meanCurvature_avg3_highCurvature.ply")
        m_split[1].export(f"/groups/cosem/cosem/ackermand/meshesForAubrey/{roi}/lowres/{object_id}_meanCurvature_avg3_lowCurvature.ply")

        avg=avg>threshold
        write_mesh(m,avg,f"/groups/cosem/cosem/ackermand/meshesForAubrey/{roi}/lowres/{object_id}_meanCurvature_avg3_thresholded.ply")
       
        