import trimesh
import numpy as np
import networkx as nx

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
        delta_s = 0;

        for one_ring_vertex_id in one_ring:
            if f"{vertex_id},{one_ring_vertex_id}" in edge_measure:
                delta_s += edge_measure[f"{vertex_id},{one_ring_vertex_id}"]
            elif f"{one_ring_vertex_id},{vertex_id}"  in edge_measure:
                delta_s -= edge_measure[f"{one_ring_vertex_id},{vertex_id}"]
        
        delta_s *= 1/(2*sum(mesh.area_faces[face_ids])/3) #use 1/3 of the areas
        mean_curv[vertex_id] = 0.5*np.linalg.norm(delta_s)
       
    return np.array(mean_curv)


def main():
    # test on sphere of radius 5
    radius = 5
    m = trimesh.creation.icosphere(radius=radius)
    print(my_discrete_mean_curvature_measure(m)/(0.5*(1/radius + 1/radius)))

if __name__ == "__main__":
    main()