#Mesh Class to coumpute different attribute of Mesh Topology like storing neighbouring elements, faces of each mesh cells
#

import numpy 

class MeshTopoogy:

	#initiate instances 
	def __init__(self, mesh):
		print("in class")
		self.MeshCells=mesh.cells[1].data
		self.Bnodes=mesh.cells[0].data
		self.MeshNodeCoordi=mesh.points[:,0:2]
		self.cellFaceID=mesh.cellFaceID

		self.BoundaryFace = []
		self.neighCellID=[]
		self.commonFaceID=[]

		self.TotalMeshElement=len(self.MeshCells)
		#these are buffere list to store temporary elements later
		self.buff=[]
		self.buff2=[]

		self.centroid=numpy.zeros((self.TotalMeshElement, 2))
		self.volume=numpy.zeros((self.TotalMeshElement))

		self.xcordi=numpy.zeros((3))
		self.ycordi=numpy.zeros((3))


	def get_neighbourCELL(self):
		# now to store neighnbouring element, we'll loop over all mesh elemeents first
		for i in range(self.TotalMeshElement):
			#loop over the each face of the  ith element
			for j in range(len(self.MeshCells[i,:])):
				#searching the face points of the ith mesh elements [Pi1,Pi2] with every elements face points
				for k in range(self.TotalMeshElement):
					# we are checking the index(ID) of mesh element which have both p1 point and p2 point present 
					if (self.MeshCells[i,j] in self.MeshCells[k,:]) and (self.MeshCells[i,j-1]  in self.MeshCells[k,:]):
						# avoiding the storing ID of itself and 
						if k!=i : #and k not in self.buff
							self.buff.append(k)
							break
					
					#if those [p1,p2] is not shared by any other element this means it's a boundary element
					if k==self.TotalMeshElement-1:
						self.buff.append(None)

				#same algorithm as bove to store common faceID
				for l in range(len(self.cellFaceID)):
					if (self.MeshCells[i,j] in self.cellFaceID[l,:]) and (self.MeshCells[i,j-1]  in self.cellFaceID[l,:]):
						self.buff2.append(l)
						break

				# for m in range(len(self.Bnodes)):
				# 	if (self.MeshCells[i,j] in self.cellFaceID[m,:]) and (self.MeshCells[i,j-1]  in self.cellFaceID[m,:]):
				# 		if m not in self.BoundaryFace:
				# 			self.BoundaryFace.append(m)

			#https://ask.sagemath.org/question/25998/why-does-append-overwriteclobber-every-existing-element-of-a-list-with-the-one-that-was-just-appended/
			
			self.neighCellID.append(self.buff[:])
			self.commonFaceID.append(self.buff2[:])

			self.buff.clear()
			self.buff2.clear()

			#This part will compute the mesh data
			for m in range(len(self.MeshCells[i,:])):
				#store ith  mesh coordinates in xcordi and  y cordi array
				[self.xcordi[m], self.ycordi[m]] = self.MeshNodeCoordi[self.MeshCells[i,m]]
				
				self.centroid[i]=[numpy.sum(self.xcordi)/3, numpy.sum(self.ycordi)/3]
				#in 2D its equal t area of mesh, by vector product for finding area for triangle
				self.volume[i]=0.5*abs(self.xcordi[0]*(self.ycordi[1]-self.ycordi[2])+self.xcordi[1]*(self.ycordi[2]-self.ycordi[0])+self.xcordi[2]*(self.ycordi[0]-self.ycordi[1]))


		self.neighCellID=numpy.array(self.neighCellID)
		self.commonFaceID=numpy.array(self.commonFaceID)

			
			


