
import numpy as np
from StructuralDiscre import LinFlux
from Meshing import struct_mesh
from parameters import runtime, dt

print("runttime",runtime)
#Class to solver
class Solver(struct_mesh):
    runtime=runtime
    dt=dt
    def  __init__(self,struct_mesh):
        self.dim=np.array([struct_mesh.nz,struct_mesh.ny,struct_mesh.nx])
        self.ap,self.aw,self.ae,self.ad,self.au,self.ab,self.af,self.ap0,self.b=LinFlux(struct_mesh)
        # self.ap=ap
        # self.aw=aw
        # self.ae=ae
        # self.ad=ad
        # self.au=au
        # self.ab=ab
        # self.af=af
        # self.ap0=ap0
        # self.b=b



    # def fun(self,x,y,z):
    #     return x*self.ap.shape[1]*self.ap.shape[2]+y*self.ap.shape[2]+z

    def set_conditions(self,mesh):
        iphi=0
        self.phi=iphi*np.ones((int(self.runtime/self.dt)+1,mesh.nz+2,mesh.ny+2,mesh.nx+2)) #include 2 ghost nodes
        print(self.phi.shape)


        #Put initital condition here
        #self.phi[0,int(mesh.nz/2-mesh.nz/5):int(mesh.nz/2+mesh.nz/5),int(mesh.ny/2-mesh.ny/5):int(mesh.ny/2+mesh.ny/5),int(mesh.nx/2-mesh.nx/5):int(mesh.nx/2+mesh.nx/5)]=10
        
        #BOUNDARY CONDITION

        #array storing face area of each control volume element in face normal to x axis
        Sx=np.dot((mesh.cvygrid[1:]-mesh.cvygrid[:-1]).reshape(-1,1), (mesh.cvzgrid[1:]-mesh.cvzgrid[:-1]).reshape(1,-1))
        #array storing face area of each control volume element in face normal to y axis
        Sy=np.dot((mesh.cvxgrid[1:]-mesh.cvxgrid[:-1]).reshape(-1,1), (mesh.cvzgrid[1:]-mesh.cvzgrid[:-1]).reshape(1,-1))
        #array storing face area of each control volume element in face normal to z axis
        Sz=np.dot((mesh.cvygrid[1:]-mesh.cvygrid[:-1]).reshape(-1,1), (mesh.cvxgrid[1:]-mesh.cvxgrid[:-1]).reshape(1,-1))

        #print("Imput Boundary Conditions as: D/N x0 y0 z0 xn yn zn")
        bc,x0,y0,z0,xn,yn,zn=['d',4,9,0,2,7,0]
        #Computing flux coefficient for dirichlet condition
        if bc=='D' or bc=='d':
            self.phi[:,:,:,1]=x0
            self.phi[:,:,1,:]=y0
            self.phi[:,1,:,:]=z0
            self.phi[:,:,:,mesh.nx]=xn
            self.phi[:,:,mesh.ny,:]=yn
            self.phi[:,mesh.nz,:,:]=zn

        #Computing flux coefficient for Neumann condition
        if bc=='N' or bc=='n':
            self.aw[:,:,0]=self.ae[:,:,nx]=0
            self.b[:,:,0]+=float(x0)*Sx
            self.b[:,:,nx]+=float(xn)*Sx
        
            self.ad[:,0,:]=self.au[:,ny,:]=0
            self.b[:,0,:]+=float(y0)*Sy
            self.b[0,ny,:]+=float(yn)*Sy
        
            self.ab[0,:,:]=self.af[nx,:,:]=0
            self.b[0,:,:]+=float(z0)*Sz
            self.b[nx,:,:]+=float(zn)*Sz
        print(self.phi[0,1:3,1:-1,1:-1])


    #this function solve linear equation array by array using line by line method
    #which is a combination of Gauss-Seidel and TDMA, we are implementing TDMA in x direction and Gauss-seidel in rest
    #Refer this solver function rather the vector solver for ease of understanding
    def solve(self,mesh):
        #since this is a iterative solver 
        #rdue solve the residual value of phi at each iteration means 
        rdue=np.ones((self.dim[0]-2,self.dim[1]-2,self.dim[2]-2),)
        
        #these 3 term are TDMA term 
        d=np.zeros((mesh.nz,mesh.ny,mesh.nx),)
        R=np.zeros(mesh.nx+2)
        Q=np.zeros(mesh.nx+2)
        print(" running slow TDMA Solver")
        # print(self.phi[0,:-2,:-2,:-2])
        print(self.phi[1,1,1,1])

        #special TDMA boundary contion 
        Q[1]=self.phi[0,0,0,1]
        Q[mesh.nx]=self.phi[0,0,0,mesh.nx]
        for t in range(1,self.phi.shape[0]):
            # print("calculating for t=",t*dt," sec")
            itr=1
            rdue.fill(1)
            error_init=0
            #we simply guessing the value of n(th) itera equal to (n-1)th  
            self.phi[t,:,:,:]=self.phi[t-1,:,:,:]
            
            while True:
                for k in range(1, self.dim[0]):
                    for j in range(1, self.dim[1]):
                        for i in range(1, self.dim[2]-1):
                            #implementing Gauss-Seidel in coefficient of y and z direction and merging it with the constant term
                            d[k,j,i]=self.b[k,j,i]+self.ap0[k,j,i]*self.phi[t-1,k+1,j+1,i+1]+self.au[k,j,i]*self.phi[t,k+1,j+1+1,i+1]
                            d[k,j,i]+=self.ad[k,j,i]*self.phi[t,k+1,j-1+1,i+1] + self.af[k,j,i]*self.phi[t,k+1+1,j+1,i+1] +self.ab[k,j,i]*self.phi[t,k+1-1,j+1,i+1]

                            R[i+1]=self.ae[k,j,i]/(self.ap[k,j,i]-self.aw[k,j,i]*R[i])
                            Q[i+1]=(d[k,j,i]+self.aw[k,j,i]*Q[i])/(self.ap[k,j,i]-self.aw[k,j,i]*R[i])
                        for i in range(self.dim[2]-1,0,-1):
                            self.phi[t,k+1,j+1,i+1]=R[i+1]*self.phi[t,k+1,j+1,i+2]+Q[i+1]
					#print("afd", Q[i+1], R[i+1])
                
                rdue=self.phi[t,2:-2,2:-2,2:-2]-np.divide((np.multiply(self.aw[1:-1,1:-1,1:-1],self.phi[t,2:-2,2:-2,1:-3])+np.multiply(self.ae[1:-1,1:-1,1:-1],self.phi[t,2:-2,2:-2,3:-1])+np.multiply(self.au[1:-1,1:-1,1:-1],self.phi[t,2:-2,3:-1,2:-2])+np.multiply(self.ad[1:-1,1:-1,1:-1], self.phi[t,2:-2,1:-3,2:-2])+np.multiply(self.af[1:-1,1:-1,1:-1],self.phi[t,3:-1,2:-2,2:-2])+np.multiply(self.ab[1:-1,1:-1,1:-1], self.phi[t,1:-3,2:-2,2:-2])+self.b[1:-1,1:-1,1:-1]+ np.multiply(self.ap0[1:-1,1:-1,1:-1],self.phi[t-1,2:-2,2:-2,2:-2])),self.ap[1:-1,1:-1,1:-1])
                


                #to to evaluate the no. of iteration for a particular time step the iterative solver should run
                #we are computing the 2nd order norm of the residual array
                #the lesser the value of err the more accurate and itera the solver will take
                err=np.linalg.norm(rdue)
                if err>1e-9:
                    # print("At itr=",itr," Error:",err)
                    itr+=1
                    if(np.absolute(err-error_init)<1e-12):
                        break
                    error_init=err
                    if itr>1000:
                        break
                else:
                    # print("done for time t:",t*self.dt," sec, Total taken itr:",itr)
                    break
        return self.phi

    #this function solve vectorized array which is much faster than array by array computation to understand this function refer solver3
    def solveByVector(self,mesh):
        rdue=np.ones((self.dim[0]-2,self.dim[1]-2,self.dim[2]-2),)
        d=np.zeros((mesh.nz,mesh.ny,mesh.nx),)
        R=np.zeros((mesh.nz,mesh.ny,mesh.nx+2),)
        Q=np.zeros((mesh.nz,mesh.ny,mesh.nx+2),)
        print("vectorized  TDMA solver")
        print(self.phi.shape)
        
        #these Bc need to implement for TDMA specifically
        Q[:,:,0]=R[:,:,0]=0
        
        #This Boundary condtion need to be implemented due to nature of TDMA Solver
        Q[:,:,1]=self.phi[0,0,0,1]
        Q[:,:,mesh.nx]=self.phi[0,0,0,mesh.nx]

        for t in range(1,self.phi.shape[0]):
            # print("calculating for t=",t*dt," sec")
            itr=1
            rdue.fill(1)
            error_init=0
            self.phi[t,:,:,:]=self.phi[t-1,:,:,:]
            while True:
                for i in range(1, self.dim[2]-1):
                    #                               self.b[k,j,i]+self.ap0[k,j,i]*self.phi[t-1,k+1,j+1,i+1]
                    d[:,1:-1,i]=self.b[:,1:-1,i]+np.multiply(self.ap0[:,1:-1,i],self.phi[t-1,1:self.dim[0]+1,2:self.dim[1],i+1])
                    #                               self.au[k,j,i]*self.phi[t,k+1,j+1+1,i+1]
                    d[:,1:-1,i]+=np.multiply(self.au[:,1:-1,i],self.phi[t,1:self.dim[0]+1,3:self.dim[1]+1,i+1])

                    #                               self.ad[k,j,i]*self.phi[t,k+1,j,i+1]
                    d[:,1:-1,i]+=np.multiply(self.ad[:,1:-1,i],self.phi[t,1:self.dim[0]+1,1:self.dim[1]-1,i+1])
                                                    # self.af[k,j,i]*self.phi[t,k+1+1,j+1,i+1]
                    d[:,1:-1,i]+=np.multiply(self.af[:,1:-1,i],self.phi[t,2:self.dim[0]+2,2:self.dim[1],i+1])

                                             # self.ab[k,j,i]*self.phi[t,k+1-1,j+1,i+1]
                    d[:,1:-1,i]+=np.multiply(self.ab[:,1:-1,i],self.phi[t,:self.dim[0],2:self.dim[1],i+1])
                    R[:,1:-1,i+1]=self.ae[:,1:-1,i]/(self.ap[:,1:-1,i]-self.aw[:,1:-1,i]*R[:,1:-1,i])
                    Q[:,1:-1,i+1]=(d[:,1:-1,i]+self.aw[:,1:-1,i]*Q[:,1:-1,i])/(self.ap[:,1:-1,i]-self.aw[:,1:-1,i]*R[:,1:-1,i])

                
                for i in range(self.dim[2]-1,-1,-1):
                    # self.phi[t,k+1,j+1,i+1]=R[i+1]*self.phi[t,k+1,j+1,i+2]+Q[i+1]
                    self.phi[t,1:self.dim[0]+1,2:self.dim[1],i+1]=np.multiply(R[:,1:-1,i+1],self.phi[t,1:self.dim[0]+1,2:self.dim[1],i+2]) +Q[:,1:-1,i+1]
                
                self.phi[t,:,1,:]=self.phi[t-1,:,1,:]

                rdue=self.phi[t,2:-2,2:-2,2:-2]-np.divide((np.multiply(self.aw[1:-1,1:-1,1:-1],self.phi[t,2:-2,2:-2,1:-3])+np.multiply(self.ae[1:-1,1:-1,1:-1],self.phi[t,2:-2,2:-2,3:-1])+np.multiply(self.au[1:-1,1:-1,1:-1],self.phi[t,2:-2,3:-1,2:-2])+np.multiply(self.ad[1:-1,1:-1,1:-1], self.phi[t,2:-2,1:-3,2:-2])+np.multiply(self.af[1:-1,1:-1,1:-1],self.phi[t,3:-1,2:-2,2:-2])+np.multiply(self.ab[1:-1,1:-1,1:-1], self.phi[t,1:-3,2:-2,2:-2])+self.b[1:-1,1:-1,1:-1]+ np.multiply(self.ap0[1:-1,1:-1,1:-1],self.phi[t-1,2:-2,2:-2,2:-2])),self.ap[1:-1,1:-1,1:-1])

                err=np.linalg.norm(rdue)
                if err>1e-8:
                    # print("At itr=",itr," Error:",err)
                    itr+=1
                    if(np.absolute(err-error_init)<1e-12):
                        break
                    error_init=err
                    if itr>1000:
                        break
                else:
                    print("done for time t:",t*self.dt," sec, Total taken itr:",itr)
                    break
        return self.phi

#     def solve1(self):
#         print("Using Gauss seedlel method")
#         runtime=1
#         dt=0.1
#         self.phi=np.zeros((int(runtime/dt)+1,self.ap.shape[0],self.ap.shape[1],self.ap.shape[2]))
#         self.phi[:,0,:,:]=10
# #         print(self.phi[1,0,:,:])
#         for t in range(1,self.phi.shape[0]):
#             xi=self.phi[t,:,:,:].flatten()
#             xf=np.copy(xi)
#             itr=1
# #             print("calculating for t=",t*dt," sec")
#             while True:
#                 for x in range(1,self.phi.shape[1]):
#                     for y in range(1,self.phi.shape[2]):
#                         for z in range(1,self.phi.shape[3]):
#                             i=self.fun(x,y,z)
#                             xf[i]=self.b[x,y,z]+self.ap0[x,y,z]*self.phi[t-1,x,y,z]
#                             if(x>0):
#                                 xf[i]+=self.aw[x,y,z]*xf[self.fun(x-1,y,z)]
#                             if(y>0):
#                                 xf[i]+=self.ad[x,y,z]*xf[self.fun(x,y-1,z)]
#                             if(z>0):
#                                 xf[i]+=self.ab[x,y,z]*xf[self.fun(x,y,z-1)]
#                             if x < self.phi.shape[1]-1:
#                                  xf[i]+=self.ae[x,y,z]*xi[self.fun(x+1,y,z)]
#                             if y < self.phi.shape[2]-1:
#                                  xf[i]+=self.au[x,y,z]*xi[self.fun(x,y+1,z)]
#                             if z < self.phi.shape[3]-1:
#                                  xf[i]+=self.af[x,y,z]*xi[self.fun(x,y,z+1)]
#                             xf[i]/=self.ap[x,y,z]
#                 if np.max(np.absolute(xf-xi))>1e-6:
#                     xi=np.copy(xf)
#                     itr+=1
#                     if itr>100:
#                         break
#                 else:
#                     print("for time t:",t," Itr:",itr)
#                     break
#             self.phi[t,:,:,:]=np.resize(xf,(self.phi.shape[1],self.phi.shape[2],self.phi.shape[3]))
#         return self.phi
