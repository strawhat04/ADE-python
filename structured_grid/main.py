from Solver import Solver
from Meshing import struct_mesh
import numpy as np
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
    mesh_elem=[10,10,10]
    # arr=np.arange(64).reshape((4,4,4))
    # print(arr[0,0,:])

#     print("for shape",mesh_elem)


    meshg=struct_mesh(dimensions,mesh_elem)
    meshg.generate_mesh()
    solver=Solver(meshg)
    solver.set_conditions(meshg)
    tp=time.perf_counter()
    phi=solver.solve(meshg)
    print("fast tdma taken time:",time.perf_counter()-tp)

#     np.savez('phi_10x10x10.npz', phi)
    
#     data=np.load('phi_10x10x10.npz')
#     phi=  data['arr_0']


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
    list.sort(files
              , key=lambda x: int(x.split('_')[1].split('.png')[0]))
    print(len(files))
    with open('images.txt', 'w') as f:
        for item in files:
            f.write("%s\n" % item)
        f.close()
    os.system('convert @images.txt {}.gif'.format("10x10x10_X_axix_YlOrBr")) # On windows convert is 'magick'
    [os.remove(ff) for ff in files]
    os.remove('images.txt')
    
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



    
