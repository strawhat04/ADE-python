import numpy
import math

#this function returns normal of face of elements with the sign 
def fnormal(MeshCoordi, facenodes, centroid):
	[x1,y1]=MeshCoordi[facenodes[0]]
	[x2,y2]=MeshCoordi[facenodes[1]]

	[xg1 ,yg1]=centroid[0]
	[xg2 , yg2]=centroid[1]

	#since face normal have 2 directions, this condition checks the outward vector
	#thi condition check whether cos() betwwen [-(y2-y1),x2-x1] and vector from current mess to neighbouring cell
	#if its value is positive means acute angle hence outward angle 
	if numpy.array([-(y2-y1),x2-x1]) @ numpy.array([xg2-xg1,yg2-yg1]) >0:
		return numpy.array([-(y2-y1),x2-x1])
	else:
		return numpy.array([(y2-y1),-(x2-x1)])

#this function compute area of face by taking face points 
def fArea(MeshCoordi,facenodes):
	
	[x1,y1]=MeshCoordi[facenodes[0]]
	[x2,y2]=MeshCoordi[facenodes[1]]

	return math.hypot(x1-x2,y1-y2) 

#advection flux coeff
def F(rho, v, MeshCoordi,facenodes,centroid):
	return -rho*(v @ fnormal(MeshCoordi,facenodes,centroid))

#diffusion flux coeff 
def gDiff(MeshCoordi,facenodes, centroid):
	#gdiff = area of face/ distanec between centroid of face  [ this is for non-skew elements ]
	[xg1 ,yg1]=centroid[0]
	[xg2 , yg2]=centroid[1]
	#fArea compute area of face and function hypot finds distance between two points
	return fArea(MeshCoordi,facenodes)/math.hypot(xg1-xg2,yg1-yg2)

#Power law advection scheme
def funA(F,D):
	if D==0:
		return 0
	else:
		return max(0,(1-0.1*abs(F/(D)))**5 )

#this function will only work for pygmsh, for different module you need to change algo to find the mesh element face location 
def get_Bcondition(MeshCoordi,c_cell,c_ni_faceNodes,cellCentroid):
	phif=numpy.array([2,5,2,2])
	bc='d'
	#cordinates of boundary face of current cell
	[bp1,bp2]=c_ni_faceNodes
	#check which face this face belongs
	
	for bindx in range(len(bfaceArray)):
		if (bp1 in bfaceArray[bindx,:]) and (bp2 in bfaceArray[bindx,:]):
			break

	for i in range(bindx+1):
		if bfaceArray[i,0] in [0,1,2,3]:
			findx=bfaceArray[i,0]
	if c_cell==12:
		print("in12 cell")
		print(bindx)
	if bc=='d':
		Db=gamma*gDiff(MeshCoordi,c_ni_faceNodes,[cellCentroid[c_cell], (MeshCoordi[bp1]+MeshCoordi[bp2])/2])
		Ab=F(rho,numpy.array([ux,uy]),MeshCoordi,c_ni_faceNodes,[cellCentroid[c_cell], (MeshCoordi[bp1]+MeshCoordi[bp2])/2])
		fmat=Db
		bmat=(Db*funA(Ab,Db) + max(Ab,0))*phif[findx] 

# fluxMat[c_cell,c_cell]=gamma*gDiff(MeshCoordi,cellFaceID[commonFace[c_cell,n_index]],[cellCentroid[c_cell], (MeshCoordi[bp1]+MeshCoordi[bp2])/2])
					# bMat[c_cell]=fluxMat[c_cell,c_cell]*bphi

	if bc=='n':
		fmat=0
		bmat=phif[findx]*fArea(MeshCoordi, c_ni_facenodes) 

	return fmat, bmat
		
class linMatrix:

	@staticmethod
	def get_linMatrix(pymesh, custom_mesh):

		MeshCells=pymesh.cells[1].data
		MeshCoordi=pymesh.points[:,0:2]
		neighbourID=custom_mesh.neighCellID
		commonFace=custom_mesh.commonFaceID
		cellCentroid = custom_mesh.centroid
		cellFaceID=pymesh.cellFaceID
		global bfaceArray
		bfaceArray=pymesh.cells[0].data
		print("bface",len(bfaceArray))
		global fface
		fface=len(pymesh.cells[2].data)

		totalMeshCells=len(MeshCells)
		fluxMat=numpy.zeros((totalMeshCells, totalMeshCells))
		bMat=numpy.zeros((totalMeshCells,1))
		bphi=5
		global gamma, rho, ux, uy
		gamma=1
		rho=1
		ux,uy=[0,0]
		#loop over all control volume elements, to discretize the eqaution
		#c_cell = current cell index of mesh
		for c_cell in range(totalMeshCells):
			#this n_index variable to track the neighbour elements and common face
			#NOTE the squential order of neighbour cell array is same as order of common face array 
			#like for ei [ep, eq, el] are neighbouring elements and [fp,fq,fl] are face share by ei cell, so fp face is shared by ei and ep and fq face is share by ei and eq , so on
			n_index=0
			#loop over each neighbouring elements
			
			for n_cell in neighbourID[c_cell]: #n_cell neighbouring cell index
				#c_ni_faceNodes gives cordinates index of common face between current cell and 1st neighbouring cell
				# like if c_ni_facenodes is [6,8] this means 6 and 8th node 
				c_ni_faceNodes=cellFaceID[commonFace[c_cell,n_index]]
				#this will handle the boundary cell elements
				if n_cell == None:
					fluxMat[c_cell,c_cell], bMat[c_cell]=get_Bcondition(MeshCoordi,c_cell,c_ni_faceNodes,cellCentroid)
				
				
				#non boundary elements
				else:
					
					#D is diffusion flux contribution
					D=gamma*gDiff(MeshCoordi,c_ni_faceNodes,[cellCentroid[c_cell],cellCentroid[n_cell]])
					#A is advection flux contribution
					A=F(rho,numpy.array([ux,uy]),MeshCoordi,c_ni_faceNodes,[cellCentroid[c_cell], cellCentroid[n_cell]])
					#General scheme by Patankar 
					fluxMat[c_cell,n_cell]=-(D*funA(A,D) + max(A,0)) 
					
				n_index+=1

			fluxMat[c_cell,c_cell]+= -( numpy.sum(fluxMat[c_cell,:]) - fluxMat[c_cell,c_cell] )

		return fluxMat, bMat