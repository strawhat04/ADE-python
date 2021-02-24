
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
       
    dimensions=[5,3,2]
    mesh_grid=[21,26,31]
    mesh_elem=[mesh_grid[0]-1, mesh_grid[1]-1, mesh_grid[2]-1]
    meshg=struct_mesh(dimensions,mesh_grid)
    meshg.generate_mesh()
    solver=solverCLS(meshg)
    solver.set_conditions(meshg)
    tp=time.perf_counter()
    phi=solver.solveByVector(meshg)
    print("xgrid", meshg.cvxgrid)
    print("ygrid", meshg.cvygrid)
    print("zgrid", meshg.cvzgrid)
    print("fast tdma taken time:",time.perf_counter()-tp)
    print(phi.shape)
    phi = np.swapaxes(phi,1,3)
    
    pwd=getcwd()
    print(pwd)
    out_path =pwd+"/out_Paraview"

    
    if not os.path.exists(out_path):
        os.makedirs(out_path)

    os.chdir(out_path)
    print(getcwd())

    # with open("sim_index", "r+") as file:
    #   sim=file.read()
    #   if sim==None:
    #       file.write(0)
    #   else:
       #        file.write(sim+1)


    if os.path.isfile(out_path+"/sim_index"):
        with open("sim_index", "r+") as file:
            sim=file.read()
            sim=str(int(float(sim)))
            file.seek(0)
            file.truncate()
            file.write(sim)

    else:
        sim='0'
        with open("sim_index", "w") as file:
            file.write(sim)
        

    out_case=out_path+"/sim" + sim
    
    if not os.path.exists(out_case):
        os.makedirs(out_case)

    os.chdir(out_case)
    print(getcwd())
    #shutil.rmtree(out_path)
    #saving the data 
    # filen='test_grid.npz'
    np.savez('ScalarMatrixValue',phi)

    # src=join(getcwd(),filen)
    # dest=join(getcwd()+"/out_Paraview/",filen)
    # print(src,getcwd()+"/out_Paraview/")
    # shutil.move(filen,getcwd()+"/out_Paraview/"+filen)
    # # os.replace(src,dest)
    
    # data_st=np.load('gr20_Vel_MID_Sc.npz')
    # phi=  data_st['arr_0']
    # del data_st
    
    

    """ VTK IMAGES"""
    colors = vtk.vtkNamedColors()
    img= vtk.vtkImageData()
    img.SetDimensions(phi.shape[1], phi.shape[2], phi.shape[3])
    img.SetSpacing(np.divide(dimensions,mesh_elem))
    # img.SetSpacing([1,1,1])00000
    img.AllocateScalars(vtk.VTK_DOUBLE, 1)
    img.SetOrigin(0.0, 0.0, 0.0)
    print(img.GetNumberOfPoints(),img.GetNumberOfCells())
    # print(img.GetCellData())
    # print(img.GetData())
    dims = img.GetDimensions()
    print(dims)
    print(phi.shape)
    

    for i in range((phi.shape[0])):
        for z in range(dims[2]):
            for y in range(dims[1]):
                for x in range(dims[0]):
                    img.SetScalarComponentFromDouble(x, y, z, 0,phi[i,x,y,z])
        writer = vtk.vtkXMLImageDataWriter()
        filen= "Test_{}.vti".format(i*dt)
        writer.SetFileName(filen)
        writer.SetInputData(img)
        writer.Write()
    
    print("done")
    os.chdir(pwd)
        