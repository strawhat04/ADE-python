from Solver import Solver
from Meshing import struct_mesh
import numpy as np
# from vtk.util.numpy_support import numpy_to_vtk
# import vtk
import time
from mayavi import mlab

    
if __name__=='__main__':
    # tp=time.perf_counter()
    # p1=solver.solve2(meshg)
    # print("slow tdma taken time:",time.perf_counter()-tp)
    # tp=time.perf_counter()
    # p2=solver.solve(meshg)
    # print("fast tdma taken time:",time.perf_counter()-tp)
    # print(np.linalg.norm(p1-p2))
    dimensions=[1,1,1]
    mesh_elem=[20,20,20]
    # arr=np.arange(64).reshape((4,4,4))
    # print(arr[0,0,:])

    print("for shape",mesh_elem)
    meshg=struct_mesh(dimensions,mesh_elem)
    meshg.generate_mesh()
    # solver=Solver(meshg)
    # solver.set_conditions(meshg)
    # tp=time.perf_counter()
    # phi=solver.solve(meshg)
    # print("fast tdma taken time:",time.perf_counter()-tp)

    data_st=np.load('phi_10x10x10.npz')
    phi=  data_st['arr_0']
    del data_st
    print(phi.shape)

    phi = np.swapaxes(phi,1,3)
    xx,yy,zz= np.meshgrid(meshg.xgrid,meshg.ygrid,meshg.zgrid,indexing='ij')
    print(xx.shape)
    
    mlab.clf()
    pt = mlab.figure(size=(800,700))
    cphi=phi[0,1:-1,1:-1,1:-1]
    pt=mlab.contour3d(xx,yy,zz,cphi,contours=1000,colormap='YlOrBr',opacity=0.15,vmax=np.max(phi),vmin=0)
    mlab.colorbar(orientation="vertical")
#     pt.scene.movie_maker.record = True
    @mlab.animate(delay=10,ui=True)
    def anim():
        fig=mlab.gcf()
        for t in range(phi.shape[0]):
            pt.mlab_source.scalars=phi[t,1:-1,1:-1,1:-1]
            mlab.title("3D ADE")
#             mlab.outline()
            mlab.axes()
            mlab.savefig(filename ='ani_%02d.png'%t)
            yield          
    anim()
    mlab.view(distance=5)
    mlab.show()

    import glob
    import os
    files=glob.glob("*.png")
    list.sort(files, key=lambda x: int(x.split('_')[1].split('.png')[0]))
    print(len(files))
    with open('images.txt', 'w') as f:
        for item in files:
            f.write("%s\n" % item)
        f.close()
    os.system('convert @images.txt {}.gif'.format("tttt_YlOrBr")) # On windows convert is 'magick'
    [os.remove(ff) for ff in files]
    os.remove('images.txt')

"""
#     np.savez('phi_10x10x10.npz', phi)
   
    data_st=np.load('phi_10x10x10.npz')
    phi=  data_st['arr_0']
    del data_st
    phi=phi[:,:-1,:-1,:-1]
    dims= phi.shape[1:]
    print(phi.shape)
    print(np.max(phi))
    print(np.min(phi))


    colors = vtk.vtkNamedColors()
    grid = vtk.vtkStructuredPoints()
    grid.SetDimensions(dims)
    grid.SetOrigin(-0.5, -0.5, -0.5)
    sp=1/200;
    grid.SetSpacing(sp, sp, sp)
    # print(sgrid.GetPoins())
    # print(sgrid.GetNumberOfCells())
    data=vtk.vtkDoubleArray()
    data.SetNumberOfComponents(1)
    data.SetNumberOfTuples((dims[0]-1)*(dims[1]-1) *(dims[2]-1))

    # cells=vtk.vtkCellData()
    for i in range(dims[0]):
        ids=i*(dims[1]*dims[2])
        for j in range(dims[1]):
            ids+=j*dims[2]
            for k in range(dims[2]):
                data.InsertTuple1(ids+k,phi[4,i,j,k])
                # points.InsertNextPoint(i,j,k)
                # data.InsertNextValue(phi[4,i,j,k])
                
    # grid.SetPoints(points)
    data.SetName("phi_val");
    grid.GetPointData().SetScalars(data)
    colorLookupTable = vtk.vtkLookupTable()
    colorLookupTable.SetTableRange(0, np.max(phi))
    colorLookupTable.Build()    
    contour = vtk.vtkContourFilter()
    contour.SetInputData(grid)
    contour.SetValue(0, 5)

    renderer = vtk.vtkRenderer()
    renWin = vtk.vtkRenderWindow()
    renWin.AddRenderer(renderer)
    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renWin)

    volMapper = vtk.vtkPolyDataMapper()
    volMapper.SetInputConnection(contour.GetOutputPort())
    volMapper.ScalarVisibilityOff()
    volActor = vtk.vtkActor()
    volActor.SetMapper(volMapper)
    volActor.GetProperty().EdgeVisibilityOn()
    volActor.GetProperty().SetColor(colors.GetColor3d("Red"))
    renderer.AddActor(volActor)
    renderer.SetBackground(colors.GetColor3d("SlateGray"))
    renWin.SetSize(512, 512)
    # Interact with the data.
    # writer = vtk.vtkStructuredPointsWriter()
    # writer.SetInputData(grid)
    # writer.SetFileName('data.vtk')
    # writer.Write()

    renWin.Render()
    iren.Start()

    # print(sgrid.GetPointData())
    # print(sgrid.GetNumberOfCells())

    # print(sgrid.GetNumberOfPoints())  
"""

    
    
#     ani = ArtistAnimation(files, interval=50, blit=True,repeat_delay=1000)
#     ani.save('dynamic_images.gif')
    
#     fig = plt.figure()
#     ax = plt.axes(projection='3d')
#     ax.contourf(phi[2,1:-1,1:-1,1:-1], 20 , cmap = "RdGy")
#     plt.colorbar()
#     ax.set_xlabel('x')
#     ax.set_ylabel('y')
#     ax.set_zlabel('z')
#     plt.show()
