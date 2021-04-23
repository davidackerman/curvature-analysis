# Curvature Analysis of Interacting Cells
`curvature.py` was used to calculate the mean curvature along a the mesh of a cell surface. For a given cellular mesh, the discrete mean curvature is calculate at each vertex. Several rounds of averaging are performed in which a vertices curvature is set to the average of itself and its nearest neighbors. This helps smooth the measurements over the cell surface. The mesh is then separated into separate meshes based on a curvature cutoff.

`surface_area_analysis.py` can be used to determine how much area of a mesh is high curvature and how much is low curvature, and how much of the high and low curvature is within a certain distance of a neighboring cell.
