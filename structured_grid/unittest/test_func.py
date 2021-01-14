import numpy 
import unittest
import dill 
import sys
import os

cwd = os.getcwd()
sys.path.append(cwd+'/..')

import StructuralDiscre
from Solver import solverCLS


'''

put dimensions=[2,2,2], mesh_grid=[20,18,15] in main.py
runtime=1, dt=0.1, flow_u=1, flow_v=1, flow_w=-1, G0=1 and r0=20 in paramater.py module
iphi=-1, bc,x0,y0,z0,xn,yn,zn=['d',5,4,1,6,2,9] in Solver.py

befor running test module

'''

class Testfunction(unittest.TestCase):


	def test_LinFlux(self):
		global testmesh
		with open('mesh20-18-15.pickle', 'rb') as input:
			testmesh=dill.load(input)
		
		tap=numpy.load('ap.npy')
		taw=numpy.load('aw.npy')
		tae=numpy.load('ae.npy')
		tad=numpy.load('ad.npy')
		tau=numpy.load('au.npy')
		tab=numpy.load('ab.npy')
		taf=numpy.load('af.npy')
		tap0=numpy.load('ap0.npy')

		ap,aw,ae,ad,au,ab,af,ap0,b=StructuralDiscre.LinFlux(testmesh)
		# print(ap-tap)
		numpy.testing.assert_array_equal(ap,tap)
		numpy.testing.assert_array_equal(aw,taw)
		numpy.testing.assert_array_equal(ae,tae)
		numpy.testing.assert_array_equal(ad,tad)
		numpy.testing.assert_array_equal(au,tau)
		numpy.testing.assert_array_equal(ab,tab)
		numpy.testing.assert_array_equal(af,taf)

	def test_solveByVector(self):

		tphi=numpy.load('phi.npy')

		solver=solverCLS(testmesh)
		solver.set_conditions(testmesh)
		phi=solver.solveByVector(testmesh)
		print(phi[5,:,:,:]-tphi[5,:,:,:])
		numpy.testing.assert_array_almost_equal(phi[:,2:-2,2:-2,2:-2],tphi[:,2:-2,2:-2,2:-2], decimal=3)

		
	

