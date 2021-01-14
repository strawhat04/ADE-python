
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
    
if __name__=='__main__':
      
 
    dimensions=[1,1,1]
    mesh_elem=[10,10,10]
    meshg=struct_mesh(dimensions,mesh_elem)
    meshg.generate_mesh()
  
    solver=solverCLS(meshg)
    solver.set_conditions(meshg)
    tp=time.perf_counter()
    phi=solver.solveByVector(meshg)
    print("fast tdma taken time:",time.perf_counter()-tp)
    print(phi.shape)
    phi = np.swapaxes(phi,1,3)
    print(phi.shape)

    # saving the data 
    # out_path = "out_Paraview/"
    filen='test_grid.npz'
    # for ff in os.listdir():
    #     print(ff)
    
    np.savez('test_grid.npz',phi)
    # src=join(getcwd(),filen)
    # dest=join(getcwd()+"/out_Paraview/",filen)
    # print(src,getcwd()+"/out_Paraview/")
    # shutil.move(filen,getcwd()+"/out_Paraview/"+filen)
    # # os.replace(src,dest)
    
    # data_st=np.load('gr20_Vel_MID_Sc.npz')
    # phi=  data_st['arr_0']
    # del data_st
    
    

    # """ VTK IMAGES"""
    colors = vtk.vtkNamedColors()
    img= vtk.vtkImageData()
    img.SetDimensions(phi.shape[1]-1, phi.shape[2]-1, phi.shape[3]-1)
    img.SetSpacing(1,1,1)
    img.AllocateScalars(vtk.VTK_DOUBLE, 1)
    img.SetOrigin(0.0, 0.0, 0.0)
    print(img.GetNumberOfPoints(),img.GetNumberOfCells())
    # print(img.GetCellData())
    # print(img.GetData())
    dims = img.GetDimensions()
    print(dims)
    print(phi.shape)
    
    for i in range(len(phi)):
        for z in range(dims[2]):
            for y in range(dims[1]):
                for x in range(dims[0]):
                    img.SetScalarComponentFromDouble(x, y, z, 0,phi[i,x,y,z])
        writer = vtk.vtkXMLImageDataWriter()
        filen= "Test_{}.vti".format(i)
        writer.SetFileName(filen)
        writer.SetInputData(img)
        writer.Write()
    print("done")
        