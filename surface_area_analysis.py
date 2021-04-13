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

def get_area_close_to_synapse(curvature_mesh, contacting_mesh, contact_distance,path):

    synapse_area = 0
    nonsynapse_area = 0

    contact_distance_squared = contact_distance*contact_distance

    #for all points in mesh of interest, calculate how far it is from contacting mesh   

    distance,_ = contacting_mesh.nearest.vertex(curvature_mesh.vertices)
    distance_thresholded = distance <= contact_distance

    face_areas = curvature_mesh.area_faces
    face_in_synapse = [False]*len(curvature_mesh.faces)
    face_in_nonsynapse = [False]*len(curvature_mesh.faces)
    for idx,face in enumerate(curvature_mesh.faces):
        if len([vertex_id for vertex_id in face if distance_thresholded[vertex_id]])>=2: #then at least two of the face vertices are within 1 um of the contacting mesh
            synapse_area += face_areas[idx]
            face_in_synapse[idx] = True
        else:
            nonsynapse_area += face_areas[idx]
            face_in_nonsynapse[idx] = True

    synapse_mesh = copy.deepcopy(curvature_mesh)
    nonsynapse_mesh = copy.deepcopy(curvature_mesh)

    synapse_mesh.update_faces(face_in_synapse)
    nonsynapse_mesh.update_faces(face_in_nonsynapse)

    synapse_mesh.export(f"{path}_synapse.ply")
    nonsynapse_mesh.export(f"{path}_nonsynapse.ply")

    return nonsynapse_area,synapse_area

ids_dictionary = {'410_roi1':[1,2],
'493_roi2':[1,3],
'494_roi2':[1,2],
'494_roi5':[2,1]
}

for key,value in ids_dictionary.items():
    high_curvature_mesh = trimesh.load(f"/groups/cosem/cosem/ackermand/meshesForAubrey/{key}/lowres/{value[0]}_meanCurvature_avg3_highCurvature.ply")
    low_curvature_mesh = trimesh.load(f"/groups/cosem/cosem/ackermand/meshesForAubrey/{key}/lowres/{value[0]}_meanCurvature_avg3_lowCurvature.ply")
    contacting_cell = trimesh.load(f"/groups/cosem/cosem/ackermand/meshesForAubrey/{key}/lowres/{value[1]}_meanCurvature_avg3.ply")
    nonsynapse_area_high, synapse_area_high = get_area_close_to_synapse(high_curvature_mesh, contacting_cell, 1000, f"/groups/cosem/cosem/ackermand/meshesForAubrey/{key}/lowres/{value[0]}_meanCurvature_avg3_highCurvature")
    nonsynapse_area_low, synapse_area_low = get_area_close_to_synapse(low_curvature_mesh, contacting_cell, 1000, f"/groups/cosem/cosem/ackermand/meshesForAubrey/{key}/lowres/{value[0]}_meanCurvature_avg3_lowCurvature")
    #temp_mesh = trimesh.load("/groups/cosem/cosem/ackermand/meshesForAubrey/410_roi1/lowres/1_meanCurvature_avg3_thresholded.ply")
    #print(f"{high_curvature_mesh.area + low_curvature_mesh.area},{temp_mesh.area},{nonsynapse_area_high+nonsynapse_area_low+synapse_area_high+synapse_area_low}")
    synapse_area = synapse_area_high+synapse_area_low
    synapse_protrusion_percent = synapse_area_high/synapse_area
    nonsynapse_area = nonsynapse_area_high+nonsynapse_area_low
    nonsynapse_protrusion_percent = nonsynapse_area_high/nonsynapse_area
    print(f"{key}")
    print(f"Nonsynapse Surface Area (nm^2) (Total, Protrusion, Protrusion Percent):\n {nonsynapse_area},{nonsynapse_area_high},{nonsynapse_protrusion_percent*100}")
    print(f"Synapse Surface Area (nm^2) (Total, Protrusion, Protrusion Percent):\n {synapse_area},{synapse_area_high},{synapse_protrusion_percent*100}")
    print("\n\n")