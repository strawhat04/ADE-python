
#import module 
from Solver import solverCLS
from Meshing import struct_mesh
from parameters import runtime, dt
#Import python default module
import numpy as np
from vtk.util.numpy_support import numpy_to_vtk
import vtk
import time
from os import getcwd ,rename
from os.path import join
import os
import shutil
from glob import glob  
# from paraview.simple import *
    
if __name__=='__main__':
      

    dimensions=[1,1,2]
    mesh_elem=[10,15,20]
    print(mesh_elem)


    meshg=struct_mesh(dimensions,mesh_elem)
    meshg.generate_mesh()
  
    # solver=solverCLS(meshg)
    # solver.set_conditions(meshg)
    # tp=time.perf_counter()
    # phi=solver.solveByVector(meshg)
    # print("fast tdma taken time:",time.perf_counter()-tp)
    # print(phi.shape)
    # phi = np.swapaxes(phi,1,3)
    # print(phi.shape)

    # # saving the data 
    # out_path = "out_Paraview/"
    # filen='test_grid.npz'
       
    # np.savez('test_grid.npz',phi)

    # src=join(getcwd(),filen)
    # dest=join(getcwd()+"/out_Paraview/",filen)
    # print(src,getcwd()+"/out_Paraview/")
    # shutil.move(filen,"/out_Paraview/"+f)
    # os.replace(filen,dest)
    
    data_st=np.load('test_grid.npz')
    phi=  data_st['arr_0']
    del data_st



    # # """ VTK IMAGES"""
    colors = vtk.vtkNamedColors()
    grid = vtk.vtkStructuredGrid()
    grid.SetDimensions(phi.shape[1:])
    field = vtk.vtkDoubleArray()
    field.SetNumberOfComponents(1)
    field.SetNumberOfTuples(phi.shape[1] * phi.shape[2] * phi.shape[3])
  
    dims = grid.GetDimensions()
    meshg.gridDims()

    points = vtk.vtkPoints()
    points.Allocate(dims[0] * dims[1] * dims[2])
    cordi = [0.0] * 3
    for i in range(len(phi)):
        for z in range(dims[2]):
            cordi[2]=meshg.zgrid[z]
            zid = z * dims[0]*dims[1]
            for y in range(dims[1]):
                cordi[1]=meshg.ygrid[y]
                yid = y * dims[0]
                for x in range(dims[0]):
                    cordi[0]=meshg.xgrid[0]
                    ids = x + yid+zid
                    points.InsertPoint(ids, cordi)
                    field.SetValue(ids,phi[i,x,y,z])
                    
        grid.SetPoints(points)
        grid.GetPointData().AddArray(field)
        ff="0"*(len(str(len(phi)))-len(str(i)))+str(i)
        write = vtk.vtkStructuredGridWriter()
        filen= "Sgrid_{}.vts".format(ff)
        write.SetFileName(filen)
        write.SetInputData(grid)
        write.Write()
    print("done")

    filelist = sorted(glob('*.vts'))
    print(filelist)
    # reader = vtk.vtkXMLImageDataReader()

    # writer = vtk.vtkXMLImageDataWriter()
    # writer.SetFileName('outfile.vti')
    # writer.SetNumberOfTimeSteps(len(filelist))
    # writer.SetTimeStepRange(0,len(filelist)-1)
    # writer.SetInputConnection(reader.GetOutputPort())
    # writer.Start()
    # for file,i in zip(filelist,range(len(filelist))):

    #     print(file)
    #     reader.SetFileName(file)
    #     reader.Modified()
    #     writer.WriteNextTime(i*dt)
    # writer.Stop()