import numpy 

class MeshTopoogy:

	def __init__(self, mesh, cellFACE):
		print("in class")
		self.MeshElements=mesh.cells[1].data
		self.Boundarypoints=mesh.cells[0].data
		self.MeshPointsCoordinates=mesh.points[:,0:2]
		self.cellFACE=cellFACE

		self.BoundaryFace = []
		self.neighbourCellID=[]
		self.commonFaceID=[]

		self.TotalMeshElement=len(self.MeshElements)
		self.neighbourCELL=[]
		self.connectingFACE=[]
		self.buff=[]
		self.centroid=numpy.zeros((self.TotalMeshElement, 2))
		self.volume=numpy.zeros((self.TotalMeshElement))

		self.xcordi=numpy.zeros((3))
		self.ycordi=numpy.zeros((3))


	test=1

	def get_neighbourCELL(self):
		for i in range(self.TotalMeshElement):
			for j in range(len(self.MeshElements[i,:])):
				#print(i,j,[self.MeshElements[i,j-1], self.MeshElements[i,j] ])
				for k in range(self.TotalMeshElement):
					if (self.MeshElements[i,j] in self.MeshElements[k,:]) and (self.MeshElements[i,j-1]  in self.MeshElements[k,:]):
						if k!=i and k not in self.buff:
							#print(i,j,k,self.buff, [self.MeshElements[i,j-1], self.MeshElements[i,j] ], self.MeshElements[k,:])
							self.buff.append(k)
							print("*******************************")
							break


			#https://ask.sagemath.org/question/25998/why-does-append-overwriteclobber-every-existing-element-of-a-list-with-the-one-that-was-just-appended/
			self.neighbourCELL.append(self.buff[:])
		#	print("neighbour",self.neighbourCELL)
			self.buff.clear()
		#	print("ckear", self.buff)
		#print(self.neighbourCELL)
		
	
	def get_meshDATA(self):
		for i in range(self.TotalMeshElement):
			for j in [0,1,2]:
				#print(self.MeshPointsCoordinates[self.MeshElements[i,j]])
				[self.xcordi[j], self.ycordi[j]] = self.MeshPointsCoordinates[self.MeshElements[i,j]]
			self.centroid[i]=[numpy.sum(self.xcordi)/3, numpy.sum(self.ycordi)/3]
			self.volume[i]=0.5*abs(self.xcordi[0]*(self.ycordi[1]-self.ycordi[2])+self.xcordi[1]*(self.ycordi[2]-self.ycordi[0])+self.xcordi[2]*(self.ycordi[0]-self.ycordi[1]))



