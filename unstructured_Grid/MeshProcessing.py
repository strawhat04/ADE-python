import numpy 

class MeshTopoogy:

	def __init__(self, mesh, cellFace):
		print("in class")
		self.MeshElements=mesh.cells[1].data
		self.Boundarypoints=mesh.cells[0].data
		self.MeshPointsCoordinates=mesh.points[:,0:2]
		self.cellFace=cellFace

		self.BoundaryFace = []
		self.neighCellID=[]
		self.commonFaceID=[]

		self.TotalMeshElement=len(self.MeshElements)
		self.buff=[]
		self.buff2=[]
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
							break
					if k==self.TotalMeshElement-1:
						self.buff.append(None)

				for l in range(len(self.cellFace)):
					if (self.MeshElements[i,j] in self.cellFace[l,:]) and (self.MeshElements[i,j-1]  in self.cellFace[l,:]):
						self.buff2.append(l)
						break

				for m in range(len(self.Boundarypoints)):
					if (self.MeshElements[i,j] in self.cellFace[m,:]) and (self.MeshElements[i,j-1]  in self.cellFace[m,:]):
						if m not in self.BoundaryFace:
							self.BoundaryFace.append(m)

			#https://ask.sagemath.org/question/25998/why-does-append-overwriteclobber-every-existing-element-of-a-list-with-the-one-that-was-just-appended/
			self.neighCellID.append(self.buff[:])
			self.commonFaceID.append(self.buff2[:])
			#print(self.commonFaceID)
		#	print("neighbour",self.neighbourCELL)
			self.buff.clear()
			self.buff2.clear()
		self.neighCellID=numpy.array(self.neighCellID)
		self.commonFaceID=numpy.array(self.commonFaceID)
		#print(self.neighbourCELL)
		#	print("ckear", self.buff)
		#print(self.neighbourCELL)
		#self.commonFaceID=numpy.array(self.commonFaceID)
	def get_meshDATA(self):
		for i in range(self.TotalMeshElement):
			for j in [0,1,2]:
				#print(self.MeshPointsCoordinates[self.MeshElements[i,j]])
				[self.xcordi[j], self.ycordi[j]] = self.MeshPointsCoordinates[self.MeshElements[i,j]]
			self.centroid[i]=[numpy.sum(self.xcordi)/3, numpy.sum(self.ycordi)/3]
			self.volume[i]=0.5*abs(self.xcordi[0]*(self.ycordi[1]-self.ycordi[2])+self.xcordi[1]*(self.ycordi[2]-self.ycordi[0])+self.xcordi[2]*(self.ycordi[0]-self.ycordi[1]))



