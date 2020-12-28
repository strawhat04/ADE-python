from Solver import Solver
from Meshing import struct_mesh
import numpy as np
from vtk.util.numpy_support import numpy_to_vtk
import vtk
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

    # arr=np.arange(64).reshape((4,4,4))
    # print(arr[0,0,:])

    # print("for shape",mesh_elem)
   
    # data_st=np.load('grid21.npz')
    # phi=  data_st['arr_0']
    # del data_st
    # print(phi.shape)

    # import glob
    # import os
    # files=glob.glob("*.png")
    # list.sort(files
    #           , key=lambda x: int(x.split('_')[1].split('.png')[0]))
    # # print(len(files))
    # with open('images.txt', 'w') as f:
    #     for item in files:
    #         f.write("%s\n" % item)
    #     f.close()
    # os.system('convert @images.txt {}.gif'.format("vtk test")) 
    
    # for i in range(len(files)):
    #     if(i%5!=0):
    #         os.remove(files[i])
    # os.remove('images.txt')
    dimensions=[1,1,1]
    mesh_elem=[20,20,20]
    meshg=struct_mesh(dimensions,mesh_elem)
    meshg.generate_mesh()
    # print(meshg.xgrid)
    # print(meshg.cvxgrid)
    # solver=Solver(meshg)
    # solver.set_conditions(meshg)
    # tp=time.perf_counter()
    # phi=solver.solve2(meshg)
    # print("fast tdma taken time:",time.perf_counter()-tp)
    # np.savez('grid21.npz',phi)

    data_st=np.load('gr20_Vel_MID_Sc.npz')
    phi=  data_st['arr_0']
    del data_st
    print(phi.shape)
    phi = np.swapaxes(phi,1,3)
    print(phi.shape)


    xx,yy,zz= np.meshgrid(meshg.xgrid,meshg.ygrid,meshg.zgrid,indexing='ij')
    print(xx.shape)

    mlab.clf()
    pt = mlab.figure(size=(800,700), bgcolor = (1,1,1), fgcolor = (0.5, 0.5, 0.5))
    cphi=phi[0,:-2,:-2,0:-2]
    # print(cphi.shape)
    pt=mlab.contour3d(xx,yy,zz,cphi,contours=1000,colormap='YlOrBr',opacity=0.15,vmax=np.max(phi),vmin=0)
    mlab.colorbar(orientation="vertical",title="range")
   
    print(mlab.view())
    # mlab.view(180,200,2)
    mlab.view(distance=4)
    mlab.title("test ",line_width=1,size=1)
    mlab.outline()
    mlab.axes()


    # mlab.text(0,0.3,text="V=0.5 in Y direction")
    @mlab.animate(delay=10,ui=True)
    def anim():
        fig=mlab.gcf()
        for t in range(phi.shape[0]):
            pt.mlab_source.scalars=phi[t,:-2,:-2,:-2]
            x=t*0.1
            mlab.ylabel("at t=%.3fsec, V=0.5 m/s"%x)
            # mlab.savefig(filename ='ani_%02d.png'%t)
            # mlab.screenshot(figure=pt, mode='rgba', antialiased=True)
            yield          
    anim()
    mlab.show()
    print("Done")

    # """ VTK IMAGES"""
    # colors = vtk.vtkNamedColors()
    # img= vtk.vtkImageData()
    # img.SetDimensions(21, 21, 21)
    # img.SetSpacing(1,1,1)
    # img.AllocateScalars(vtk.VTK_DOUBLE, 1)
    # img.SetOrigin(0.0, 0.0, 0.0)
    # print(img.GetNumberOfPoints(),img.GetNumberOfCells())
    # # print(img.GetCellData())
    # # print(img.GetData())
    # dims = img.GetDimensions()
    # print(dims)
    # for i in range(len(phi)):
    #     for z in range(dims[2]):
    #         for y in range(dims[1]):
    #             for x in range(dims[0]):
    #                 img.SetScalarComponentFromDouble(x, y, z, 0,phi[i,x,y,z])
    #     writer = vtk.vtkXMLImageDataWriter()
    #     writer.SetFileName("Test_{}.vti".format(i))
    #     writer.SetInputData(img)
    #     writer.Write()
    #     print("done")
        
    
    

    # writer = vtk.vtkXMLImageDataWriter()
    # writer.SetFileName("Test.vti")
    # writer.SetInputData(img)
    # writer.Write()

    # # Read the file (to test that it was written correctly)
    # reader = vtk.vtkXMLImageDataReader()
    # reader.SetFileName("Test.v")
    # reader.Update()

    # # Convert the image to a polydata
    # imageDataGeometryFilter = vtk.vtkImageDataGeometryFilter()
    # imageDataGeometryFilter.SetInputConnection(reader.GetOutputPort())
    # imageDataGeometryFilter.Update()

    # mapper = vtk.vtkPolyDataMapper()
    # mapper.SetInputConnection(imageDataGeometryFilter.GetOutputPort())

    # actor = vtk.vtkActor()
    # actor.SetMapper(mapper)
    # actor.GetProperty().SetPointSize(3)

    # # Setup rendering
    # renderer = vtk.vtkRenderer()
    # renderer.AddActor(actor)
    # renderer.SetBackground(colors.GetColor3d('White'))
    # renderer.ResetCamera()

    # renderWindow = vtk.vtkRenderWindow()
    # renderWindow.AddRenderer(renderer)

    # renderWindowInteractor = vtk.vtkRenderWindowInteractor()

    # renderWindowInteractor.SetRenderWindow(renderWindow)
    # renderWindowInteractor.Initialize()
    # renderWindowInteractor.Start()


    

    """
#     np.savez('phi_10x10x10.npz', phi)

    data_st=np.load('Vel2x_grid20.npz')
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
