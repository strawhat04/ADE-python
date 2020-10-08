import numpy
import math

def tr(a):
	return np.transpose(a)

#advection flux coeff


def F(rho, u,dcv):
	return rho*u*dcv
#diffusion flux coeff
def gDiff(MeshPoints,facepoints, centroid):
	[x1,y1]=MeshPoints[facepoints[0]]
	[x2,y2]=MeshPoints[facepoints[1]]
	[xg1 ,yg2]=centroid[0]
	[xg1 , yg2]=centroid[1]
	return math.hypot(x1-x2,y1-y2)/math.hypot(x1-x2,y1-y2)
#Power law advection scheme
def funA(F,D):
	if D==0:
		return 0
	else:
		return max(0,(1-0.1*abs(F/(D)))**5 )

class linMatrix:

	@staticmethod
	def get_linMatrix(pymesh, custom_mesh, cellface):

		MeshCells=pymesh.cells[1].data
		MeshPoints=pymesh.points[:,0:2]
		neighbourID=custom_mesh.neighCellID
		commonFace=custom_mesh.commonFaceID
		cellCentroid = custom_mesh.centroid

		totalMeshCells=len(MeshCells)
		fluxMat=numpy.zeros((totalMeshCells, totalMeshCells))
		bMat=numpy.zeros((totalMeshCells,1))
		
		gamma=1
		for eno in range(totalMeshCells):

			index=0
			for cell in neighbourID[eno]:

				#Boundary Cells need to be treated differently
				if None in neighbourID[eno]:
					#BOUNDARY CELL
					pass

				else:
					print("indexing", eno, cell)
					print("commonface",commonFace[eno,index])
					print(cellface[commonFace[eno,index]])

					fluxMat[eno,cell]=gamma*gDiff(MeshPoints,cellface[commonFace[eno,index]],[cellCentroid[eno],cellCentroid[cell]])
					print("fluxmat:",fluxMat[eno,cell])

