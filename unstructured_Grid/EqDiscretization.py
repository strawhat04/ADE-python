import numpy
import math


#advection flux coeff
def F(rho, v, MeshPoints,facepoints,centroid):
	return -rho*(v @ fnormal(MeshPoints,facepoints,centroid))

def fnormal(MeshPoints, facepoints, centroid):
	[x1,y1]=MeshPoints[facepoints[0]]
	[x2,y2]=MeshPoints[facepoints[1]]

	[xg1 ,yg1]=centroid[0]
	[xg2 , yg2]=centroid[1]

	if numpy.array([-(y2-y1),x2-x1]) @ numpy.array([xg2-xg1,yg2-yg1]) >0:
		return numpy.array([-(y2-y1),x2-x1])
	else:
		return numpy.array([(y2-y1),-(x2-x1)])

def fArea(MeshPoints,facepoints):
	[x1,y1]=MeshPoints[facepoints[0]]
	[x2,y2]=MeshPoints[facepoints[1]]

	return math.hypot(x1-x2,y1-y2) 

#diffusion flux coeff
def gDiff(MeshPoints,facepoints, centroid):
	
	[xg1 ,yg1]=centroid[0]
	[xg2 , yg2]=centroid[1]
	return fArea(MeshPoints,facepoints)/math.hypot(xg1-xg2,yg1-yg2)

#Power law advection scheme
def funA(F,D):
	if D==0:
		return 0
	else:
		return max(0,(1-0.1*abs(F/(D)))**5 )

class linMatrix:

	@staticmethod
	def get_linMatrix(pymesh, custom_mesh, cellFaceID):

		MeshCells=pymesh.cells[1].data
		MeshPoints=pymesh.points[:,0:2]
		neighbourID=custom_mesh.neighCellID
		commonFace=custom_mesh.commonFaceID
		cellCentroid = custom_mesh.centroid

		totalMeshCells=len(MeshCells)
		fluxMat=numpy.zeros((totalMeshCells, totalMeshCells))
		bMat=numpy.zeros((totalMeshCells,1))
		bflux=.5
		gamma=1
		rho=1
		v=1
		for eno in range(totalMeshCells):
			index=0
			for cell in neighbourID[eno]:
				if cell == None:
					bMat[eno]=bflux*fArea(MeshPoints,cellFaceID[commonFace[eno,index]])
					index+=1
				else:
					D=gamma*gDiff(MeshPoints,cellFaceID[commonFace[eno,index]],[cellCentroid[eno],cellCentroid[cell]])
					A=F(rho,numpy.array([1,0]),MeshPoints,cellFaceID[commonFace[eno,index]],[cellCentroid[eno], cellCentroid[cell]])
					fluxMat[eno,cell]=D*funA(A,D) + max(A,0) 
					index+=1
					if eno==6:
						print(cell,D,A,funA(A,D))
		#print(neighbourID)	

		#print("flux Matrix: ",fluxMat)
		#print("B matrix:", bMat)

		return fluxMat, bMat