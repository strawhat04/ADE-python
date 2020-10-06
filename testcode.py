import meshplex
import numpy as np
import pygmsh
from MeshTopology import MeshTopoogy

def meshFaceID(mesh):
	#perfect
	k=0
	cells=mesh.cells[1].data
	cellFACE = []
	#print("#1", len(cells), len(cells[0,:]))
	for i in range(len(cells)-1):
		for j in range(len(cells[0,:])):
			#print(cells[i,:], cells[i,j-1], cells[i,j], cellFACE)
			if [ cells[i,j-1], cells[i,j]] and [cells[i,j], cells[i,j-1]] not in cellFACE:
				#print("in if")
				cellFACE.append([cells[i,j-1], cells[i,j]])

	return np.array(cellFACE)


# with pygmsh.geo.Geometry() as geom:
#     geom.add_circle([0.0, 0.0], 1.0, mesh_size=0.7)
#     mesh1 = geom.generate_mesh()

with pygmsh.geo.Geometry() as geom:
    geom.add_polygon(
        [
            [0.0, 0.0],
            [0, 1],
            [1.0, 1],
            [1, .0],
        ],
        mesh_size=0.025,
    )
    mesh1= geom.generate_mesh()
print("start")
print(mesh1.points[:,0:2])

print("\n---------------\n")
print(mesh1.cells)

# c=np.sort(mesh1.cells[1].data, axis=1)
# c=c[c[:,0].argsort(kind='mergesort')]
# print(c)

cellFACE=meshFaceID(mesh1)
#print(mesh1.cells)
# print(mesh1.cells)
# print("\n--------sorting cell:1-------------\n")

mesh = MeshTopoogy(mesh1,cellFACE)
mesh.get_neighbourCELL()
mesh.get_meshDATA()
print("Neighbour: ", mesh.neighbourCELL)
print("Mesh Volume: ",mesh.volume)
print(len(mesh1.cells[1].data))
mesh3=meshplex.MeshTri(mesh1.points, mesh1.cells[1].data)
#print(mesh.cells)

mesh3.show()